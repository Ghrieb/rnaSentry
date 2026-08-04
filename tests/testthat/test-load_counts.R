test_that("load_counts constructs a valid SummarizedExperiment", {
  counts <- matrix(1:12, nrow = 3,
                   dimnames = list(c("TP53", "BRCA1", "MYC"),
                                   c("S1", "S2", "S3", "S4")))
  coldata <- data.frame(
    time = c(10, 12, 8, 15),
    event = c(1, 0, 1, 0),
    row.names = colnames(counts))
  se <- load_counts(counts, coldata)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(assayNames(se), "counts")
  expect_equal(nrow(se), 3)
  expect_equal(ncol(se), 4)
  expect_true(all(colnames(se) == colnames(counts)))
})

test_that("load_counts works without colData", {
  counts <- matrix(1:6, nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2", "S3")))
  se <- load_counts(counts)
  expect_s4_class(se, "SummarizedExperiment")
  expect_equal(ncol(se), 3)
})

test_that("load_counts rejects non-matrix input", {
  expect_error(load_counts(data.frame(x = 1:3)), "'counts' must be a matrix")
})

test_that("load_counts rejects non-numeric matrix", {
  m <- matrix(letters[1:6], nrow = 2,
              dimnames = list(c("A", "B"), c("S1", "S2", "S3")))
  expect_error(load_counts(m), "'counts' must be a numeric matrix")
})

test_that("load_counts rejects empty matrix", {
  m <- matrix(numeric(0), nrow = 0, ncol = 0)
  expect_error(load_counts(m), "no rows or no columns")
})

test_that("load_counts rejects missing dimnames", {
  m <- matrix(1:6, nrow = 2)
  expect_error(load_counts(m), "rownames.*column names")
})

test_that("load_counts detects duplicate sample IDs", {
  counts <- matrix(1:6, nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S1", "S2")))
  coldata <- S4Vectors::DataFrame(x = 1:3, row.names = c("S1", "S1", "S2"))
  expect_warning(se <- load_counts(counts, coldata), "Duplicate")
  expect_equal(ncol(se), 2)
  expect_true(any(S4Vectors::metadata(se)$flags$check == "duplicate_samples"))
})

test_that("load_counts warns on non-integer counts", {
  counts <- matrix(c(1.5, 2, 3, 4), nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_warning(se <- load_counts(counts), "non-integer")
  expect_true(any(S4Vectors::metadata(se)$flags$check == "non_integer_counts"))
})

test_that("load_counts warns on all-zero rows", {
  counts <- matrix(c(0, 5, 0, 3), nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  expect_warning(se <- load_counts(counts), "all-zero")
  expect_true(any(S4Vectors::metadata(se)$flags$check == "empty_rows"))
})

test_that("load_counts rejects colData with missing sample IDs", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  coldata <- data.frame(x = 1:2, row.names = c("S1", "S3"))
  expect_error(load_counts(counts, coldata), "missing row names")
})

test_that("load_counts detects Ensembl IDs", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("ENSG000001", "ENSG000002"),
                                   c("S1", "S2")))
  expect_warning(se <- load_counts(counts), "Ensembl")
  expect_true(any(S4Vectors::metadata(se)$flags$check == "ensembl_ids_detected"))
})

test_that("load_counts uses custom assay_name", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  se <- load_counts(counts, assay_name = "expression")
  expect_equal(assayNames(se), "expression")
})

test_that("load_counts returns flags as a data.frame", {
  counts <- matrix(1:4, nrow = 2,
                   dimnames = list(c("A", "B"), c("S1", "S2")))
  se <- load_counts(counts)
  expect_true(is.data.frame(S4Vectors::metadata(se)$flags))
})
