library(SummarizedExperiment)

test_that("generate_report validates its inputs", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  km <- km_curve(sig, se)
  expect_error(generate_report(list(km)), "named list")
  expect_error(generate_report(list(km = km), output_dir = "nope"),
               "existing directory")
})

test_that("generate_report renders a pipeline report", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  stages <- list(
    km_curve = km_curve(sig, se),
    cox_model = cox_model(sig, se),
    survival_parametric = survival_parametric(sig, se)
  )
  out <- file.path(tempdir(), "rnaSentry_test_report.html")
  if (file.exists(out)) unlink(out)
  invisible(generate_report(stages, output_file = basename(out),
                            output_dir = tempdir()))
  expect_true(file.exists(out))
  expect_gt(file.info(out)$size, 0)
  unlink(out)
})

test_that("generate_report renders the pipeline assumptions section", {
  skip_if_not_installed("rmarkdown")
  skip_if_not_installed("knitr")
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  stages <- list(
    km_curve = km_curve(sig, se),
    cox_model = cox_model(sig, se)
  )
  out <- file.path(tempdir(), "rnaSentry_assumptions_report.html")
  if (file.exists(out)) unlink(out)
  invisible(generate_report(stages, output_file = basename(out),
                            output_dir = tempdir()))
  html <- paste(readLines(out, warn = FALSE), collapse = "\n")
  expect_match(html, "Pipeline assumptions and scope")
  expect_match(html, "right-censored")
  expect_match(html, "proportional-hazards")
  unlink(out)
})
