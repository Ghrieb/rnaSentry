library(SummarizedExperiment)

test_that("run_rnaSentry runs the full pipeline and renders a report", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  se <- make_survival_se()
  run <- run_rnaSentry(se, "time", "event", top_n = 10, repeats = 2,
                       folds = 3, seed = 7,
                       design_vars = c("condition", "age"),
                       design_terms = "condition",
                       report_file = "run_report.html",
                       report_dir = tempdir())
  expect_s3_class(run, "rnaSentry_run")
  expect_true(all(c("pca_audit", "build_signature", "km_curve", "cox_model",
                    "survival_parametric") %in% names(run$stages)))
  expect_true("design_audit" %in% names(run$stages))
  expect_true(isTRUE(run$stages$km_curve$sig_locked))
  expect_true(isTRUE(run$stages$cox_model$sig_locked))
  expect_true(isTRUE(run$stages$build_signature$locked))
  expect_true(file.exists(run$report))
  expect_gt(file.info(run$report)$size, 0)
  unlink(run$report)
})

test_that("run_rnaSentry can skip the report and the design audit", {
  se <- make_survival_se()
  run <- run_rnaSentry(se, "time", "event", top_n = 10, repeats = 2,
                       folds = 3, seed = 7, render_report = FALSE)
  expect_null(run$report)
  expect_false("design_audit" %in% names(run$stages))
  expect_output(print(run), "pipeline run")
})

test_that("run_rnaSentry validates inputs", {
  expect_error(run_rnaSentry(list(), "time", "event"),
               "SummarizedExperiment")
  expect_error(run_rnaSentry(make_survival_se(), "time", "event",
                             report_dir = "nope"),
               "existing directory")
})
