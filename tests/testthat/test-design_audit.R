library(SummarizedExperiment)

test_that("design_audit rejects non-SummarizedExperiment input", {
  expect_error(design_audit(matrix(1, 4, 4), design_vars = "x"),
               "SummarizedExperiment")
})

test_that("design_audit validates design_vars", {
  se <- make_pca_se()
  expect_error(design_audit(se, design_vars = character(0)),
               "non-empty")
  expect_error(design_audit(se, design_vars = 1L), "character vector")
  expect_error(design_audit(se, design_vars = "nope"), "no column")
  expect_error(design_audit(se, design_vars = c("batch", "batch")),
               "duplicates")
})

test_that("design_audit validates outcome_col", {
  se <- make_pca_se()
  expect_error(design_audit(se, design_vars = "batch", outcome_col = "nope"),
               "no column")
  expect_error(design_audit(se, design_vars = c("batch", "condition"),
                            outcome_col = "condition"),
               "must not also be listed")
})

test_that("design_audit validates thresholds", {
  se <- make_pca_se()
  expect_error(design_audit(se, "batch", surrogate_pc = 0), "positive integer")
  expect_error(design_audit(se, "batch", surrogate_pc = 1.5), "positive integer")
  expect_error(design_audit(se, "batch", alpha = 1), "single number in")
  expect_error(design_audit(se, "batch", interaction_alpha = 0),
               "single number in")
  expect_error(design_audit(se, "batch", redundant_effect_size = 2),
               "single number in")
})

test_that("design_audit flags a batch-driven surrogate and builds a formula", {
  se <- make_pca_se()
  res <- design_audit(se, design_vars = c("batch", "age"),
                      outcome_col = "condition")
  expect_s3_class(res, "rnaSentry_design")
  expect_true(all(c("formula", "formula_text", "rationale",
                    "confounder_table", "interaction_tested",
                    "interaction_table", "pairwise_table", "flags") %in%
                    names(res)))
  expect_s3_class(res$formula, "formula")
  expect_true(grepl("condition", res$formula_text))
  expect_true(grepl("batch", res$formula_text))
  expect_true(is.logical(res$interaction_tested))
  expect_true(!is.null(res$interaction_table))
  expect_true(!is.null(res$pairwise_table))
  expect_equal(res$surrogate_pc, 1)
  expect_true(all(c("variable", "type", "n", "test", "statistic", "df",
                    "p_value", "effect_size", "effect_size_type",
                    "flagged") %in% colnames(res$confounder_table)))
  flagged <- res$confounder_table$variable[res$confounder_table$flagged]
  expect_true("batch" %in% flagged)
  expect_true(length(res$rationale) >= 4)
})

test_that("design_audit treats a numeric variable as a regression test", {
  se <- make_pca_se()
  res <- design_audit(se, design_vars = "age")
  expect_equal(res$confounder_table$test, "lm")
  expect_equal(res$confounder_table$effect_size_type, "r_squared")
  expect_false(grepl("condition", res$formula_text))
})

test_that("design_audit reports a single-level variable as untested", {
  se <- make_pca_se()
  cd <- SummarizedExperiment::colData(se)
  cd$constant <- factor(rep("X", ncol(se)))
  SummarizedExperiment::colData(se) <- cd
  res <- design_audit(se, design_vars = c("batch", "constant"))
  row <- res$confounder_table[res$confounder_table$variable == "constant", ]
  expect_equal(row$test, "none")
  expect_false(isTRUE(row$flagged))
})

test_that("design_audit flags a fully-missing variable", {
  se <- make_pca_se()
  cd <- SummarizedExperiment::colData(se)
  cd$gone <- rep(NA_real_, ncol(se))
  SummarizedExperiment::colData(se) <- cd
  res <- design_audit(se, design_vars = "gone")
  expect_equal(res$confounder_table$test, "none")
  expect_true(any(res$flags$check == "variable_no_data"))
})

test_that("design_audit drops redundant numeric variables from the formula", {
  se <- make_pca_se()
  cd <- SummarizedExperiment::colData(se)
  cd$age2 <- cd$age
  SummarizedExperiment::colData(se) <- cd
  res <- design_audit(se, design_vars = c("batch", "age", "age2"))
  expect_true("age2" %in% res$dropped_vars)
  expect_false(grepl("age2", res$formula_text))
  expect_true(grepl("age", res$formula_text))
  expect_true(any(res$flags$check == "redundant_variable"))
})

test_that("design_audit clamps surrogate_pc to the available number of PCs", {
  se <- make_pca_se()
  res <- design_audit(se, design_vars = "batch", surrogate_pc = 100)
  expect_true(any(res$flags$check == "surrogate_pc_clamped"))
  expect_true(res$surrogate_pc < 100)
})

test_that("design_audit skips the interaction scan with a single PC", {
  counts <- matrix(stats::rpois(10, lambda = 100), nrow = 5, ncol = 2,
                   dimnames = list(paste0("g", 1:5), c("A", "B")))
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(batch = factor(c("B1", "B2")),
                                    row.names = c("A", "B")))
  res <- design_audit(se, design_vars = "batch")
  expect_false(res$interaction_tested)
  expect_null(res$interaction_table)
})

test_that("print.rnaSentry_design is informative", {
  se <- make_pca_se()
  res <- design_audit(se, design_vars = "batch", outcome_col = "condition")
  expect_output(print(res), "Recommended design formula")
})
