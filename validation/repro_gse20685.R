# Build inst/extdata/gse20685_case_study.rds: a deterministic, offline-safe
# subset of the GSE20685 breast-cancer cohort for the case-study vignette.
#
# GSE20685 (Li et al., 2010) is an Affymetrix GPL570 series with overall-
# survival follow-up. The full probe-level series matrix is ~60 MB and is NOT
# bundled; only the top-variance gene subset below is shipped so the vignette
# builds without network access. Full-cohort face-validity statistics live in
# validation/face_validity_review.md (CV C = 0.783, random-gene control 0.582,
# log-rank p = 3.8e-17).
#
# Usage (dev only; GEOquery is not a package dependency):
#   Rscript validation/repro_gse20685.R
#
# If the raw eset was already downloaded (e.g. cached), point
# GSE20685_ESET_RDS at it to skip the download.

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(SummarizedExperiment)
})

cache <- Sys.getenv("GSE20685_ESET_RDS",
                    "C:/Users/Hani/AppData/Local/Temp/opencode/geo_data/series.rds")

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
expr <- exprs(eset)
fd <- fData(eset)
sym <- as.character(fd[["Gene Symbol"]])
cat(sprintf("probes: %d, symbols mapped: %d (%.1f%%)\n",
            nrow(expr), sum(!is.na(sym) & sym != "" & sym != "---"),
            100 * sum(!is.na(sym) & sym != "" & sym != "---") / nrow(expr)))

p <- pData(eset)
time_years <- as.numeric(as.character(p[["follow_up_duration (years):ch1"]]))
event_death <- as.integer(as.character(p[["event_death:ch1"]]))
age <- as.numeric(as.character(p[["age at diagnosis:ch1"]]))
subtype <- as.character(p[["subtype:ch1"]])
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

out <- "inst/extdata/gse20685_case_study.rds"
dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
# xz keeps the bundled file under BiocCheck's 5MB data-file guideline
saveRDS(se, out, compress = "xz")
cat(sprintf("wrote %s (%.1f MB)\n", out,
            file.info(out)$size / 1e6))
cat("colData columns:", paste(colnames(colData(se)), collapse = ", "), "\n")
cat(sprintf("events = %d, time range = %.1f-%.1f years\n",
            sum(se$event), min(se$time), max(se$time)))
