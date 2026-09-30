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
  expect_identical(SummarizedExperiment::assayNames(out), "logcounts")
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
    SummarizedExperiment::assay(out_es, "logcounts"),
    SummarizedExperiment::assay(out_se, "logcounts")
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
    SummarizedExperiment::assay(out_se, "logcounts"),
    SummarizedExperiment::assay(out_ranged, "logcounts")
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

test_that(".brca_fallback_cohort matches the fetch contract and is deterministic", {
  fb1 <- rnaSentry:::.brca_fallback_cohort()
  fb2 <- rnaSentry:::.brca_fallback_cohort()
  expect_s4_class(fb1, "SummarizedExperiment")
  expect_equal(nrow(fb1), 3000L)
  expect_equal(ncol(fb1), 327L)
  expect_identical(SummarizedExperiment::assayNames(fb1), "logcounts")
  expect_true(all(c("time", "event", "age", "subtype") %in%
    colnames(SummarizedExperiment::colData(fb1))))
  expect_equal(
    SummarizedExperiment::assay(fb1, "logcounts"),
    SummarizedExperiment::assay(fb2, "logcounts")
  )
  # the helper must not disturb the caller's RNG stream
  set.seed(99)
  before <- runif(1)
  set.seed(99)
  invisible(rnaSentry:::.brca_fallback_cohort())
  expect_identical(before, runif(1))
})

