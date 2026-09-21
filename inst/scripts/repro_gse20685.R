# Reproduce GSE20685 subset: a deterministic processing of the GSE20685
# breast-cancer cohort via GEOquery (output NOT committed to inst/extdata
# per Bioconductor data guidelines; vignettes use synthetic offline fallback
# with the same structure — run this script locally to regenerate).
#
# This is the shipped copy of validation/repro_gse20685.R (kept in sync for
# review). Previously the output was inst/extdata/gse20685_case_study.rds
# for the case-study-brca vignette; now use tempdir()/cache (do not commit).
#
# Source and licensing:
#   GSE20685 (Li et al., 2010, "A five-gene molecular grade index and
#   HOXB13:IL17BR are complementary prognostic factors in early stage breast
#   cancer") is an Affymetrix GPL570 (HG-U133 Plus 2.0) series deposited in
#   the NCBI Gene Expression Omnibus (https://www.ncbi.nlm.nih.gov/geo/query/
#   acc.cgi?acc=GSE20685). It is redistributed here as a processed gene-level
#   subset under NCBI GEO's public-domain terms of use.
#
# GSE20685 is an Affymetrix GPL570 series with overall-survival follow-up. The
# full probe-level series matrix is ~60 MB and is NOT bundled; only the
# top-variance gene subset below is shipped so the vignette builds without
# network access. Full-cohort face-validity statistics live in
# validation/face_validity_review.md (CV C = 0.783, random-gene control 0.582,
# log-rank p = 3.8e-17).
#
# Usage (dev only; GEOquery is not a package dependency):
#   Rscript inst/scripts/repro_gse20685.R   # run from the package root
#
# If the raw eset was already downloaded (e.g. cached), point
# GSE20685_ESET_RDS at it to skip the download.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(SummarizedExperiment)
})

cache <- Sys.getenv("GSE20685_ESET_RDS",
                    file.path(tempdir(), "gse20685_eset.rds"))

if (file.exists(cache)) {
  cat(sprintf("[cache] reading eset from %s\n", cache))
  geo <- readRDS(cache)
  if (is.list(geo) && !is(geo, "ExpressionSet") && !is.null(geo$eset)) {
    eset <- geo$eset
  } else {
    eset <- geo
  }
} else {
  cat("[download] fetching GSE20685 series matrix from GEO ...\n")
  geo <- GEOquery::getGEO("GSE20685", GSEMatrix = TRUE, AnnotGPL = TRUE)
  eset <- geo[[1]]
}

# ---- probe -> gene collapse (largest-mean probe per symbol), as in Gate 5 --
# NOTE: newer GEOquery versions return a (Ranged)SummarizedExperiment
# instead of an ExpressionSet. Handle both (mirrors R/geo_utils.R).
if (is(eset, "ExpressionSet")) {
  expr <- exprs(eset)
  fd <- fData(eset)
  p <- pData(eset)
} else if (is(eset, "SummarizedExperiment")) {
  expr <- as.matrix(SummarizedExperiment::assay(eset, 1))
  fd <- as.data.frame(SummarizedExperiment::rowData(eset))
  p <- as.data.frame(SummarizedExperiment::colData(eset))
} else {
  stop("unsupported GEO object class: ", paste(class(eset), collapse = ", "))
}
# Normalization-insensitive matching: DataFrame() sanitizes names on the
# SE path ("Gene symbol" -> "Gene.symbol"). Mirrors R/geo_utils.R.
norm_nm <- function(x) gsub("[^a-z0-9]", "", tolower(x))
sym_hit <- which(norm_nm(colnames(fd)) == "genesymbol")
if (length(sym_hit) == 0L) stop("gene-symbol column not found in feature data")
sym <- as.character(fd[[sym_hit[1]]])
cat(sprintf("probes: %d, symbols mapped: %d (%.1f%%)\n",
            nrow(expr), sum(!is.na(sym) & sym != "" & sym != "---"),
            100 * sum(!is.na(sym) & sym != "" & sym != "---") / nrow(expr)))

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
cat(sprintf("samples with usable OS metadata: %d / %d (%d deaths)\n",
            sum(keep), ncol(eset), sum(event_death[keep] == 1)))
expr <- expr[, keep]
age <- age[keep]
subtype <- subtype[keep]

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
cat(sprintf("gene-symbol matrix: %d genes x %d samples\n",
            nrow(expr_gene), ncol(expr_gene)))

# ---- subset to the most variable genes (bounded size for inst/extdata) -----
n_top <- 3000
v <- apply(expr_gene, 1, stats::var)
v[!is.finite(v)] <- 0
ord <- order(v, decreasing = TRUE)
genes_top <- rownames(expr_gene)[ord[seq_len(n_top)]]
expr_sub <- expr_gene[genes_top, , drop = FALSE]
cat(sprintf("bundled subset: %d genes x %d samples (top %d by variance)\n",
            nrow(expr_sub), ncol(expr_sub), n_top))

coldata <- S4Vectors::DataFrame(
  time = time_years[keep],
  event = event_death[keep],
  age = age,
  subtype = factor(subtype),
  row.names = colnames(expr_sub)
)
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(logcounts = expr_sub), colData = coldata
)

out <- file.path(tempdir(), "gse20685_case_study.rds")
# Previously: "inst/extdata/gse20685_case_study.rds" (not committed; see NEWS 0.99.4)
# Cache via BiocFileCache for reuse: BiocFileCache::BiocFileCache()$add(...)
saveRDS(se, out, compress = "xz")
cat(sprintf("wrote %s (%.1f MB) — not committed to inst/extdata; copy manually if needed\n",
            out, file.info(out)$size / 1e6))
cat("colData columns:", paste(colnames(colData(se)), collapse = ", "), "\n")
cat(sprintf("events = %d, time range = %.1f-%.1f years\n",
            sum(se$event), min(se$time), max(se$time)))
