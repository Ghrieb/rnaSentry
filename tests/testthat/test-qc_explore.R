library(SummarizedExperiment)

make_test_se <- function() {
  set.seed(42)
  counts <- matrix(rpois(60, lambda = 200), nrow = 6, ncol = 10,
                    dimnames = list(paste0("gene", 1:6), paste0("S", 1:10)))
  coldata <- S4Vectors::DataFrame(condition = rep(c("A", "B"), each = 5),
                        row.names = paste0("S", 1:10))
  SummarizedExperiment(assays = list(counts = counts), colData = coldata)
}

test_that("qc_explore rejects non-SummarizedExperiment input", {
  expect_error(qc_explore(matrix(1:4, 2, 2)), "SummarizedExperiment")
})

test_that("qc_explore returns expected structure on clean data", {
  se <- make_test_se()
  qc <- qc_explore(se)
  expect_s3_class(qc, "rnaSentry_qc")
  expect_equal(qc$n_samples, 10)
  expect_equal(qc$n_genes, 6)
  expect_length(qc$duplicate_samples, 0)
  expect_false(qc$non_integer_counts)
})

test_that("qc_explore flags duplicate sample IDs", {
  se <- make_test_se()
  colnames(se)[2] <- colnames(se)[1]
  qc <- qc_explore(se)
  expect_true("duplicate_samples" %in% qc$flags$check)
})

test_that("qc_explore flags library-size outliers", {
  se <- make_test_se()
  assay(se)[, 1] <- assay(se)[, 1] * 100L
  qc <- qc_explore(se, mad_threshold = 2)
  expect_true(colnames(se)[1] %in% qc$library_size_outliers)
})

test_that("qc_explore flags non-integer counts", {
  se <- make_test_se()
  assay(se)[1, 1] <- 3.5
  qc <- qc_explore(se)
  expect_true(qc$non_integer_counts)
})

test_that("print.rnaSentry_qc runs without error", {
  se <- make_test_se()
  qc <- qc_explore(se)
  expect_output(print(qc), "rnaSentry QC audit")
})
