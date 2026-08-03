library(SummarizedExperiment)

test_that("cox_model validates its inputs", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_error(cox_model(list(genes = "a"), se), "rnaSentry_signature")
  expect_error(cox_model(sig, matrix(1, 4, 4)), "SummarizedExperiment")
  bad <- sig
  bad$genes <- c(bad$genes, "not_a_gene")
  expect_error(cox_model(bad, se), "not_a_gene")
  sig$time_col <- "nope"
  expect_error(cox_model(sig, se), "no column 'nope'")
})

test_that("cox_model validates confounders", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_error(cox_model(sig, se, confounders = 1), "character vector")
  expect_error(cox_model(sig, se, confounders = c("age", "age")), "duplicates")
  sig$design_terms <- "age"
  expect_error(cox_model(sig, se, confounders = "age"), "duplicate")
  sig$design_terms <- character(0)
  expect_error(cox_model(sig, se, confounders = sig$genes[1]), "signature genes")
  expect_error(cox_model(sig, se, confounders = "nope"), "no covariate")
  SummarizedExperiment::colData(se)$constant <- rep("x", ncol(se))
  expect_error(cox_model(sig, se, confounders = "constant"), "no variation")
})

test_that("cox_model returns the documented structure", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cm <- cox_model(sig, se, confounders = "age")
  expect_s3_class(cm, "rnaSentry_cox_model")
  expect_true(all(sig$genes %in% cm$terms))
  expect_true("age" %in% cm$terms)
  expect_true(all(c("term", "coefficient", "se", "HR", "HR_low",
                    "HR_high", "p") %in% colnames(cm$coef_table)))
  expect_true(all(c("HR", "HR_low", "HR_high", "p", "adj_p") %in%
                    colnames(cm$gene_summary)))
  expect_equal(nrow(cm$gene_summary), length(sig$genes))
  expect_true(is.finite(cm$concordance) && cm$concordance > 0 &&
                cm$concordance < 1)
  expect_true(is.logical(cm$ph_violated))
  expect_true("GLOBAL" %in% cm$zph_summary$term)
  expect_equal(cm$n, ncol(se))
  expect_equal(cm$events, sum(SummarizedExperiment::colData(se)$event))
  expect_false(cm$sig_locked)
  expect_true(nzchar(cm$created))
})

test_that("cox_model gene coefficients match the gene-only signature fit", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cm <- cox_model(sig, se)
  b <- stats::setNames(cm$coef_table$coefficient, cm$coef_table$term)
  expect_equal(unname(b[sig$genes]), unname(sig$coefficients[sig$genes]),
               tolerance = 1e-8)
})

test_that("cox_model includes design terms recorded in the signature", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7, design_terms = "batch")
  cm <- cox_model(sig, se)
  expect_true(any(grepl("^batch", cm$coef_table$term)))
})

test_that("cox_model flags missing covariate values and estimates by complete cases", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7, design_terms = "age")
  SummarizedExperiment::colData(se)$age[1] <- NA_real_
  cm <- cox_model(sig, se)
  expect_true(any(cm$flags$check == "missing_covariates"))
  expect_equal(cm$n, ncol(se) - 1)
})

test_that("cox_model records the lock state of the signature", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_false(cox_model(sig, se)$sig_locked)
  expect_true(cox_model(lock_signature(sig), se)$sig_locked)
})

test_that("print and plot are informative and return the object", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  cm <- cox_model(sig, se)
  expect_output(print(cm), "Concordance")
  expect_output(print(cm), "Signature genes")
  grDevices::pdf(NULL)
  res <- plot(cm)
  grDevices::dev.off()
  expect_identical(res, cm)
})
