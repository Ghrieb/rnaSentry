library(SummarizedExperiment)

test_that("pca_audit returns expected structure on clean data", {
  se <- make_pca_se()
  res <- pca_audit(se)
  expect_s3_class(res, "rnaSentry_pca")
  expect_s3_class(res$pca, "prcomp")
  expect_equal(length(res$percent_variance), ncol(res$pca$x))
  expect_true(abs(sum(res$percent_variance) - 100) < 1e-6)
  expect_true(all(names(res$percent_variance) ==
                    paste0("PC", seq_along(res$percent_variance))))
  expect_equal(dim(res$scores), c(30, 5))
  expect_equal(nrow(res$col_data), 30)
  expect_equal(res$n_genes_filtered, 0)
  expect_equal(res$filtered_reason, "")
  expect_null(res$batch_tests)
  expect_equal(res$assay_used, "log2(counts+1)")
})

test_that("pca_audit prefers logcounts then vst assays", {
  se <- make_pca_se()
  SummarizedExperiment::assay(se, "vst") <-
    as.matrix(SummarizedExperiment::assay(se, 1)) + 0.5
  expect_equal(pca_audit(se)$assay_used, "vst")
  SummarizedExperiment::assay(se, "logcounts") <-
    as.matrix(SummarizedExperiment::assay(se, 1)) + 1
  expect_equal(pca_audit(se)$assay_used, "logcounts")
})

test_that("pca_audit flags PCs associated with a separating batch", {
  se <- make_pca_se()
  res <- pca_audit(se, batch_col = "batch")
  expect_true("batch_associated_pc" %in% res$flags$check)
  expect_true(any(res$batch_tests$flagged))
  expect_true(all(res$batch_tests$effect_size >= 0))
  expect_equal(unique(res$batch_tests$effect_size_type), "eta_squared")
})

test_that("pca_audit tests a numeric batch via linear regression", {
  se <- make_pca_se()
  se$continuous_batch <- rep(c(0, 1), each = ncol(se) / 2) +
    stats::rnorm(ncol(se), 0, 0.05)
  res <- pca_audit(se, batch_col = "continuous_batch")
  expect_true(any(res$batch_tests$test == "lm"))
  expect_true(any(res$batch_tests$effect_size_type == "r_squared"))
})

test_that("pca_audit reports a single-level batch without error", {
  se <- make_pca_se()
  se$one_level <- "x"
  res <- pca_audit(se, batch_col = "one_level")
  expect_true("batch_single_level" %in% res$flags$check)
  expect_null(res$batch_tests)
})

test_that("pca_audit does not flag a null batch assignment", {
  set.seed(99)
  se <- make_pca_se()
  se$random_batch <- factor(sample(c("r1", "r2"), ncol(se), replace = TRUE))
  res <- pca_audit(se, batch_col = "random_batch", batch_alpha = 0.01)
  expect_false(any(res$batch_tests$flagged, na.rm = TRUE))
})

test_that("pca_audit rejects invalid input", {
  expect_error(pca_audit(matrix(1:4, 2, 2)), "SummarizedExperiment")
  se <- make_pca_se()
  expect_error(pca_audit(se, batch_col = "nope"), "no column")
  expect_error(pca_audit(se, top_n_pcs = 0), "positive integer")
  expect_error(pca_audit(se, top_n_pcs = "x"), "positive integer")
  expect_error(pca_audit(se, batch_alpha = 1), "\\(0, 1\\)")
})

test_that("pca_audit drops constant genes and reports the filter", {
  se <- make_pca_se()
  assay(se)[1, ] <- 100
  res <- pca_audit(se)
  expect_equal(res$n_genes_filtered, 1)
  expect_true("pca_gene_filter" %in% res$flags$check)
  expect_true(grepl("zero variance", res$filtered_reason))
  expect_equal(res$n_genes_analyzed, 29)
})

test_that("pca_audit stops when too few variable genes remain", {
  se <- make_pca_se(n_genes = 2, n_samples = 10)
  expect_error(pca_audit(se), "Too few variable genes")
})

test_that("pca_audit reports a BH-adjusted p-value alongside raw p", {
  se <- make_pca_se()
  res <- pca_audit(se, batch_col = "batch")
  expect_true("adj_p" %in% colnames(res$batch_tests))
  expect_equal(res$batch_tests$adj_p,
               stats::p.adjust(res$batch_tests$p_value, method = "BH"))
})

test_that("print.rnaSentry_pca runs without error", {
  se <- make_pca_se()
  res <- pca_audit(se, batch_col = "batch")
  expect_output(print(res), "rnaSentry PCA audit")
})

test_that("plot_pca_audit returns a ggplot and validates input", {
  skip_if_not_installed("ggplot2")
  se <- make_pca_se()
  res <- pca_audit(se, batch_col = "batch")
  p <- plot_pca_audit(res, color_by = "batch")
  expect_s3_class(p, "ggplot")
  expect_error(plot_pca_audit(list()), "rnaSentry_pca")
  expect_error(plot_pca_audit(res, color_by = "nope"), "no column")
  expect_error(plot_pca_audit(res, pc_x = 99), "pc_x")
})
