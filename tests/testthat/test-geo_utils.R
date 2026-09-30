# Mock GEO objects mirroring GSE20685 structure at small scale: 60 probes
# mapping to 30 gene symbols (2 probes each), 40 samples with survival and
# clinical metadata. The SummarizedExperiment mock round-trips feature and
# sample metadata through S4Vectors::DataFrame, reproducing the name
# sanitization seen on builders with newer GEOquery.
mock_geo_objects <- function() {
  set.seed(1)
  n_probes <- 60L
  n_samples <- 40L
  mat <- matrix(
    rnorm(n_probes * n_samples, mean = 6, sd = 1.5),
    nrow = n_probes, ncol = n_samples,
    dimnames = list(
      paste0("probe", seq_len(n_probes)),
      paste0("S", seq_len(n_samples))
    )
  )
  fd <- data.frame(
    "Gene symbol" = rep(paste0("G", seq_len(30)), each = 2),
    check.names = FALSE,
    row.names = rownames(mat)
  )
  pd <- data.frame(
    "follow_up_duration (years):ch1" = runif(n_samples, 1, 10),
    "event_death:ch1" = sample(0:1, n_samples, replace = TRUE),
    "age at diagnosis:ch1" = round(rnorm(n_samples, 55, 10)),
    "subtype:ch1" = sample(c("Basal", "LumA"), n_samples, replace = TRUE),
    check.names = FALSE,
    row.names = colnames(mat)
  )
  eset <- Biobase::ExpressionSet(
    assayData = mat,
    featureData = Biobase::AnnotatedDataFrame(fd),
    phenoData = Biobase::AnnotatedDataFrame(pd)
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat),
    rowData = S4Vectors::DataFrame(fd),
    colData = S4Vectors::DataFrame(pd)
  )
  ranged <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat),
    colData = S4Vectors::DataFrame(pd),
    rowRanges = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(
        start = seq_len(n_probes) * 1000L,
        width = 500L
      ),
      # mcols go through DataFrame sanitization ("Gene symbol" ->
      # "Gene.symbol"), exactly as on builders with newer GEOquery
      `Gene symbol` = rep(paste0("G", seq_len(30)), each = 2)
    )
  )
  list(eset = eset, se = se, ranged = ranged)
}

test_that("fetch_gse20685 handles ExpressionSet input", {
  skip_if_not_installed("GEOquery")
  skip_if_not_installed("Biobase")
  m <- mock_geo_objects()
  out <- testthat::with_mocked_bindings(
    getGEO = function(...) list(m$eset),
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_s4_class(out, "SummarizedExperiment")
  expect_identical(SummarizedExperiment::assayNames(out), "exprs")
  # 60 probes collapse to 30 unique gene symbols
  expect_equal(nrow(out), 30L)
  expect_equal(ncol(out), 40L)
  expect_true(all(c("time", "event", "age", "subtype") %in%
    colnames(SummarizedExperiment::colData(out))))
  expect_true(all(!is.na(SummarizedExperiment::colData(out)$time)))
  expect_true(all(SummarizedExperiment::colData(out)$event %in% c(0L, 1L)))
})

test_that("fetch_gse20685 agrees on ExpressionSet and SummarizedExperiment input", {
  skip_if_not_installed("GEOquery")
  skip_if_not_installed("Biobase")
  m <- mock_geo_objects()
  out_es <- testthat::with_mocked_bindings(
    getGEO = function(...) list(m$eset),
    .package = "GEOquery",
    fetch_gse20685()
  )
  # DataFrame round-trip sanitizes names ("Gene symbol" -> "Gene.symbol",
  # "follow_up_duration (years):ch1" -> "follow_up_duration..years..ch1"),
  # exactly as on builders with newer GEOquery.
  out_se <- testthat::with_mocked_bindings(
    getGEO = function(...) list(m$se),
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_true(!is.null(out_es) && !is.null(out_se))
  expect_equal(
    SummarizedExperiment::assay(out_es, "exprs"),
    SummarizedExperiment::assay(out_se, "exprs")
  )
  expect_equal(
    as.data.frame(SummarizedExperiment::colData(out_es)),
    as.data.frame(SummarizedExperiment::colData(out_se))
  )
})

test_that("fetch_gse20685 handles RangedSummarizedExperiment input", {
  skip_if_not_installed("GEOquery")
  skip_if_not_installed("Biobase")
  skip_if_not_installed("GenomicRanges")
  m <- mock_geo_objects()
  out_se <- testthat::with_mocked_bindings(
    getGEO = function(...) list(m$se),
    .package = "GEOquery",
    fetch_gse20685()
  )
  out_ranged <- testthat::with_mocked_bindings(
    getGEO = function(...) list(m$ranged),
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_true(!is.null(out_ranged))
  expect_equal(
    SummarizedExperiment::assay(out_se, "exprs"),
    SummarizedExperiment::assay(out_ranged, "exprs")
  )
})

test_that("fetch_gse20685 returns NULL (not error) on bad input or failed download", {
  skip_if_not_installed("GEOquery")
  out_bad <- testthat::with_mocked_bindings(
    getGEO = function(...) list(data.frame(x = 1)),
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_null(out_bad)
  out_dl <- testthat::with_mocked_bindings(
    getGEO = function(...) stop("network down"),
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_null(out_dl)
})

test_that("fetch_gse20685 retries transient failures then succeeds", {
  skip_if_not_installed("GEOquery")
  skip_if_not_installed("Biobase")
  m <- mock_geo_objects()
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  flaky_getGEO <- function(...) {
    calls$n <- calls$n + 1L
    if (calls$n < 3L) stop("transient FTP timeout")
    list(m$eset)
  }
  out <- testthat::with_mocked_bindings(
    getGEO = flaky_getGEO,
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_true(!is.null(out))
  expect_equal(calls$n, 3L)
  expect_equal(nrow(out), 30L)
})

test_that("fetch_gse20685 gives up after 3 attempts and returns NULL", {
  skip_if_not_installed("GEOquery")
  calls <- new.env(parent = emptyenv())
  calls$n <- 0L
  dead_getGEO <- function(...) {
    calls$n <- calls$n + 1L
    stop("network down")
  }
  out <- testthat::with_mocked_bindings(
    getGEO = dead_getGEO,
    .package = "GEOquery",
    fetch_gse20685()
  )
  expect_null(out)
  expect_equal(calls$n, 3L)
})

test_that(".brca_fallback_cohort matches the fetch contract and is deterministic", {
  fb1 <- rnaSentry:::.brca_fallback_cohort()
  fb2 <- rnaSentry:::.brca_fallback_cohort()
  expect_s4_class(fb1, "SummarizedExperiment")
  expect_equal(nrow(fb1), 3000L)
  expect_equal(ncol(fb1), 327L)
  expect_identical(SummarizedExperiment::assayNames(fb1), "exprs")
  expect_true(all(c("time", "event", "age", "subtype") %in%
    colnames(SummarizedExperiment::colData(fb1))))
  expect_equal(
    SummarizedExperiment::assay(fb1, "exprs"),
    SummarizedExperiment::assay(fb2, "exprs")
  )
  # the helper must not disturb the caller's RNG stream
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(rnaSentry:::.brca_fallback_cohort())
  expect_identical(before, runif(1))
})

test_that("exprs assay is treated as log-ready across consumers", {
  # Same log-scale matrix under two names: build_signature and km_curve must
  # behave identically. Without "exprs" in the dispatcher, the microarray
  # matrix would be double-logged via log2(mat + 1).
  set.seed(11)
  n_genes <- 60L
  n_samples <- 50L
  mat <- matrix(rnorm(n_genes * n_samples, mean = 6, sd = 1.5),
                nrow = n_genes, ncol = n_samples,
                dimnames = list(paste0("g", seq_len(n_genes)),
                                paste0("S", seq_len(n_samples))))
  sig_expr <- colMeans(mat[seq_len(5), , drop = FALSE])
  risk <- scale(sig_expr)[, 1] * 0.6
  event_time <- rexp(n_samples, rate = 0.05 * exp(risk))
  censor_time <- rexp(n_samples, rate = 0.03)
  cd <- S4Vectors::DataFrame(time = pmin(event_time, censor_time),
                             event = as.integer(event_time < censor_time),
                             row.names = colnames(mat))
  se_log <- SummarizedExperiment::SummarizedExperiment(
    assays = list(logcounts = mat), colData = cd)
  se_ex <- SummarizedExperiment::SummarizedExperiment(
    assays = list(exprs = mat), colData = cd)
  sig_log <- build_signature(se_log, "time", "event", top_n = 5,
                             repeats = 1, folds = 3, seed = 11)
  sig_ex <- build_signature(se_ex, "time", "event", top_n = 5,
                            repeats = 1, folds = 3, seed = 11)
  expect_identical(sig_ex$genes, sig_log$genes)
  expect_equal(sig_ex$cv_summary, sig_log$cv_summary)
  km_log <- km_curve(sig_log, se_log)
  km_ex <- km_curve(sig_ex, se_ex)
  expect_equal(km_ex$cutpoint, km_log$cutpoint)
  expect_equal(km_ex$log_rank_p, km_log$log_rank_p)
})

