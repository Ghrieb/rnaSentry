#' Fetch and process GSE20685 breast cancer cohort
#'
#' Downloads GSE20685 (Li *et al.*, 2010) from GEO via
#' \code{GEOquery::getGEO()}, performs probe-to-gene collapse (largest mean
#' expression per symbol), subsets to the top 3000 most variable genes, and
#' returns a \code{SummarizedExperiment} with \code{logcounts} assay and
#' clinical metadata (\code{time}, \code{event}, \code{age}, \code{subtype}).
#' Results are cached via \code{BiocFileCache} for subsequent calls.
#'
#' @param cache Logical. If \code{TRUE} (default), cache the processed
#'   \code{SummarizedExperiment} via \code{BiocFileCache}. If \code{FALSE},
#'   always re-download and re-process.
#'
#' @return A \code{SummarizedExperiment} or \code{NULL} if the download or
#'   processing fails (e.g. network unavailable, \code{GEOquery} not
#'   installed).
#'
#' @details
#' This helper packages the processing pipeline from
#' \code{inst/scripts/repro_gse20685.R} for use in vignettes. The full
#' series matrix is downloaded from GEO, probes are collapsed to gene symbols
#' by largest mean expression, and the top 3000 most variable genes are
#' retained. Clinical metadata includes overall-survival follow-up time
#' (years), event indicator (1 = death), age at diagnosis, and breast-cancer
#' subtype.
#'
#' @examples
#' \donttest{
#' se <- fetch_gse20685()
#' se
#' }
#'
#' @importFrom Biobase exprs fData pData
#' @importFrom SummarizedExperiment assay rowData colData
#' @export
fetch_gse20685 <- function(cache = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    message("GEOquery not installed; returning NULL.")
    return(NULL)
  }

  bfc <- NULL
  rname <- "rnaSentry_gse20685"
  if (cache && requireNamespace("BiocFileCache", quietly = TRUE)) {
    bfc <- BiocFileCache::BiocFileCache()
    cached <- BiocFileCache::bfcquery(bfc, rname, "rname", exact = TRUE)
    if (nrow(cached) > 0L) {
      message("Loading cached GSE20685 ...")
      cached_se <- tryCatch(
        readRDS(BiocFileCache::bfcrpath(bfc, rname)),
        error = function(e) {
          message("Cached GSE20685 unreadable: ", conditionMessage(e))
          NULL
        }
      )
      if (!is.null(cached_se)) return(cached_se)
    }
  }

  message("Downloading GSE20685 from GEO ...")
  # NOTE: newer GEOquery versions return a (Ranged)SummarizedExperiment
  # instead of an ExpressionSet for some series. Handle both, and return
  # NULL (vignette synthetic fallback) on any failure so a GEO format
  # change can never hard-fail the vignette build.
  se <- tryCatch(
    {
      # Bound the download time (Bioconductor Appendix C: web queries must
      # fail quickly on nightly builders); never shorten a user-configured
      # longer timeout.
      geo <- withr::with_options(
        list(timeout = max(300, getOption("timeout"))),
        GEOquery::getGEO("GSE20685", GSEMatrix = TRUE, AnnotGPL = TRUE)
      )
      if (is.null(geo) || length(geo) == 0L) stop("empty GEO result")

      eset <- geo[[1]]

      if (methods::is(eset, "ExpressionSet")) {
        expr <- Biobase::exprs(eset)
        fd <- Biobase::fData(eset)
        p <- Biobase::pData(eset)
      } else if (methods::is(eset, "SummarizedExperiment")) {
        expr <- as.matrix(SummarizedExperiment::assay(eset, 1))
        fd <- as.data.frame(SummarizedExperiment::rowData(eset))
        p <- as.data.frame(SummarizedExperiment::colData(eset))
      } else {
        stop(
          "unsupported GEO object class: ",
          paste(class(eset), collapse = ", ")
        )
      }

      # --- probe -> gene collapse (largest-mean probe per symbol) ---
      # Column is "Gene symbol" (lowercase s) on GPL570 via ExpressionSet,
      # but S4Vectors::DataFrame() sanitizes names ("Gene.symbol") on the
      # SummarizedExperiment path. Match normalization-insensitively so
      # both object types (and future annotation tweaks) resolve.
      norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
      sym_hit <- which(norm_nm(colnames(fd)) == "genesymbol")
      if (length(sym_hit) == 0L) {
        stop("gene-symbol column not found in feature data")
      }
      sym <- as.character(fd[[sym_hit[1]]])

  has_sym <- !is.na(sym) & sym != "" & sym != "---"
  expr <- expr[has_sym, , drop = FALSE]
  sym <- sym[has_sym]

  means <- rowMeans(expr)
  o <- order(sym, -means)
  sym_o <- sym[o]
  expr_o <- expr[o, , drop = FALSE]
  dups <- duplicated(sym_o)
  expr_gene <- expr_o[!dups, , drop = FALSE]
  rownames(expr_gene) <- sym_o[!dups]

  # --- clinical metadata ---
  # `p` already holds the sample metadata from the class branch above.
  # Names may be sanitized on the SE path
  # ("follow_up_duration (years):ch1" -> "follow_up_duration..years..ch1"),
  # so match normalization-insensitively.
  norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
  find_col <- function(key) {
    hit <- which(norm_nm(colnames(p)) == key)
    if (length(hit) == 0L) return(NA_character_)
    colnames(p)[hit[1]]
  }
  time_col <- find_col("followupdurationyearsch1")
  event_col <- find_col("eventdeathch1")
  age_col <- find_col("ageatdiagnosisch1")
  subtype_col <- find_col("subtypech1")
  if (any(is.na(c(time_col, event_col, age_col, subtype_col)))) {
    stop("expected clinical columns not found in GEO sample metadata")
  }
  time_years <- as.numeric(as.character(p[[time_col]]))
  event_death <- as.integer(as.character(p[[event_col]]))
  age <- as.numeric(as.character(p[[age_col]]))
  subtype <- as.character(p[[subtype_col]])

  keep <- !is.na(time_years) & !is.na(event_death) & time_years >= 0 &
    event_death %in% c(0L, 1L)
  expr_gene <- expr_gene[, keep, drop = FALSE]
  time_years <- time_years[keep]
  event_death <- event_death[keep]
  age <- age[keep]
  subtype <- subtype[keep]

  # --- subset to top variance genes ---
  n_top <- min(3000L, nrow(expr_gene))
  v <- apply(expr_gene, 1, stats::var, na.rm = TRUE)
  v[!is.finite(v)] <- 0
  ord <- order(v, decreasing = TRUE)
  genes_top <- rownames(expr_gene)[ord[seq_len(n_top)]]
  expr_sub <- expr_gene[genes_top, , drop = FALSE]

  # --- build SummarizedExperiment ---
  coldata <- S4Vectors::DataFrame(
    time = time_years,
    event = event_death,
    age = age,
    subtype = factor(subtype),
    row.names = colnames(expr_sub)
  )
      SummarizedExperiment::SummarizedExperiment(
        assays = list(logcounts = expr_sub),
        colData = coldata
      )
    },
    error = function(e) {
      message("fetch_gse20685 failed: ", conditionMessage(e), "; returning NULL.")
      NULL
    }
  )
  if (is.null(se)) return(NULL)

  # --- cache if available ---
  if (!is.null(bfc)) {
    tryCatch(
      BiocFileCache::bfcadd(bfc, rname, rfile = function() {
        tf <- tempfile(fileext = ".rds")
        saveRDS(se, tf, compress = "xz")
        tf
      }),
      error = function(e) {
        message("BiocFileCache write failed: ", conditionMessage(e))
      }
    )
  }

  se
}

#' Synthetic BRCA fallback cohort (internal)
#'
#' Offline fallback for the breast-cancer case-study vignette when the live
#' GSE20685 download via \code{fetch_gse20685()} fails. Same structure as the
#' real output (3000 genes x 327 samples, \code{logcounts} assay, \code{time} /
#' \code{event} / \code{age} / \code{subtype} metadata) with a fixed seed so
#' the fallback is reproducible. Not exported.
#'
#' @return A \code{SummarizedExperiment}.
#' @noRd
.brca_fallback_cohort <- function() {
  # Scoped seed (no global RNG side effects; a bare seed call in R/ code is a
  # BiocCheck WARNING).
  withr::with_seed(20685, {
    n_genes <- 3000
    n_samples <- 327
    genes <- paste0("GENE", seq_len(n_genes))
    samples <- paste0("Sample", seq_len(n_samples))
    mat <- matrix(rnorm(n_genes * n_samples, mean = 6, sd = 1.5),
                  nrow = n_genes, ncol = n_samples,
                  dimnames = list(genes, samples))
    risk0 <- colMeans(mat[seq_len(20), , drop = FALSE])
    event_time <- rexp(n_samples, rate = 0.08 * exp(0.6 * scale(risk0)[, 1]))
    censor_time <- rexp(n_samples, rate = 0.04)
    time <- pmin(event_time, censor_time)
    event <- as.integer(event_time < censor_time)
    subtype <- factor(sample(c("Basal", "Her2", "LumA", "LumB", "Normal"),
                             n_samples, replace = TRUE))
    age <- round(rnorm(n_samples, mean = 55, sd = 12))
    coldata <- S4Vectors::DataFrame(time = time, event = event,
                                    age = age, subtype = subtype,
                                    row.names = samples)
    SummarizedExperiment::SummarizedExperiment(
      assays = list(logcounts = mat), colData = coldata
    )
  })
}
