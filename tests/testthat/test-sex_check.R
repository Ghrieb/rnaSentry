library(SummarizedExperiment)

make_sex_se <- function(mislabel_last = TRUE) {
  genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "GAPDH")
  samples <- paste0("S", 1:6)
  expr <- matrix(0, nrow = length(genes), ncol = length(samples),
                  dimnames = list(genes, samples))
  expr["XIST", 1:3] <- c(950, 900, 875)
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 1:3] <- 5
  expr["XIST", 4:6] <- 5
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 4:6] <- 400
  expr["GAPDH", ] <- 1000

  reported <- c("F", "F", "F", "M", "M", if (mislabel_last) "F" else "M")
  coldata <- DataFrame(sex = reported, row.names = samples)
  SummarizedExperiment(assays = list(counts = expr), colData = coldata)
}

test_that("sex_check rejects non-SummarizedExperiment input", {
  expect_error(sex_check(matrix(1, 2, 2)), "SummarizedExperiment")
})

test_that("sex_check errors when sex_col is missing", {
  se <- make_sex_se()
  expect_error(sex_check(se, sex_col = "nope"), "no column")
})

test_that("sex_check errors when no marker genes are present", {
  se <- make_sex_se()
  rownames(se)[1] <- "NOT_XIST"
  expect_error(sex_check(se), "XIST/Y-gene markers")
})

test_that("sex_check correctly flags a mislabeled sample", {
  se <- make_sex_se(mislabel_last = TRUE)
  res <- sex_check(se)
  expect_equal(res$status[res$sample_id == "S6"], "MISMATCH")
  expect_true(all(res$status[res$sample_id %in% paste0("S", 1:5)] == "OK"))
})

test_that("sex_check reports OK when metadata matches expression", {
  se <- make_sex_se(mislabel_last = FALSE)
  res <- sex_check(se)
  expect_true(all(res$status == "OK"))
})

test_that("sex_check reports AMBIGUOUS when neither signal is dominant", {
  genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D")
  samples <- paste0("S", 1:6)
  expr <- matrix(0, nrow = 4, ncol = 6, dimnames = list(genes, samples))
  expr["XIST", ] <- c(1000, 900, 600, 5, 5, 5)
  y_vals <- matrix(rep(c(5, 5, 400, 400, 400, 400), each = 3),
                   nrow = 3, byrow = TRUE)
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), ] <- y_vals
  coldata <- DataFrame(sex = c("F", "F", "F", "M", "M", "M"),
                        row.names = samples)
  se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
  res <- sex_check(se)
  # S3 has comparable XIST (600) and Y-gene (400) expression: rank score
  # within one unit of zero -> neither signal clearly dominant.
  expect_equal(res$status[res$sample_id == "S3"], "AMBIGUOUS")
  expect_true(all(res$status[res$sample_id %in% c("S1", "S2")] == "OK"))
})
