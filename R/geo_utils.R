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
#' @export
fetch_gse20685 <- function(cache = TRUE) {
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    message("GEOquery not installed; returning NULL.")
    return(NULL)
  }
  if (!requireNamespace("Biobase", quietly = TRUE)) {
    message("Biobase not installed; returning NULL.")
    return(NULL)
  }

  bfc <- NULL
  rname <- "rnaSentry_gse20685"
  if (cache && requireNamespace("BiocFileCache", quietly = TRUE)) {
    bfc <- BiocFileCache::BiocFileCache()
    cached <- BiocFileCache::bfcquery(bfc, rname, "rname", exact = TRUE)
    if (nrow(cached) > 0L) {
      message("Loading cached GSE20685 ...")
      return(readRDS(BiocFileCache::bfcrpath(bfc, rname)))
    }
  }

  message("Downloading GSE20685 from GEO ...")
  geo <- tryCatch(
    GEOquery::getGEO("GSE20685", GSEMatrix = TRUE, AnnotGPL = TRUE),
    error = function(e) {
      message("GEOquery download failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(geo)) return(NULL)

  eset <- geo[[1]]

  # --- probe -> gene collapse (largest-mean probe per symbol) ---
  expr <- Biobase::exprs(eset)
  fd <- Biobase::fData(eset)
  sym <- as.character(fd[["Gene symbol"]])

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
  p <- Biobase::pData(eset)
  time_years <- as.numeric(as.character(p[["follow_up_duration (years):ch1"]]))
  event_death <- as.integer(as.character(p[["event_death:ch1"]]))
  age <- as.numeric(as.character(p[["age at diagnosis:ch1"]]))
  subtype <- as.character(p[["subtype:ch1"]])

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
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = expr_sub),
    colData = coldata
  )

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
