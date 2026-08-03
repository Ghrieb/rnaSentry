library(SummarizedExperiment)

test_that("km_curve validates its inputs", {
  expect_error(km_curve(list(genes = "a"), make_survival_se()),
               "rnaSentry_signature")
  sig <- make_signature_for_testing()
  expect_error(km_curve(sig, matrix(1, 4, 4)), "SummarizedExperiment")
})

test_that("km_curve requires the signature genes in the data", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  bad <- sig
  bad$genes <- c(bad$genes, "not_a_gene")
  expect_error(km_curve(bad, se), "not_a_gene")
})

test_that("km_curve reports missing survival columns", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  sig$time_col <- "nope"
  expect_error(km_curve(sig, se), "no column 'nope'")
})

test_that("km_curve errors on constant risk scores", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  sig$coefficients <- stats::setNames(rep(0, length(sig$genes)), sig$genes)
  expect_error(km_curve(sig, se), "constant")
})

test_that("km_curve returns the documented structure on a built signature", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  km <- km_curve(sig, se)
  expect_s3_class(km, "rnaSentry_km")
  expect_equal(length(km$score), ncol(se))
  expect_setequal(levels(km$groups), c("low", "high"))
  expect_true(all(table(km$groups) >= 1))
  expect_true(is.finite(km$log_rank_p) && km$log_rank_p >= 0 &&
                km$log_rank_p <= 1)
  expect_true(all(c("low", "high") %in% names(km$median_survival)))
  expect_false(km$sig_locked)
  expect_true(nzchar(km$created))
})

test_that("km_curve separates high-risk from low-risk survival", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 8)
  km <- km_curve(sig, se)
  ml <- km$median_survival[["low"]]
  mh <- km$median_survival[["high"]]
  if (is.finite(ml) && is.finite(mh)) {
    expect_gt(ml, mh)
  }
})

test_that("km_curve records the lock state of the signature", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_false(km_curve(sig, se)$sig_locked)
  expect_true(km_curve(lock_signature(sig), se)$sig_locked)
})

test_that("print and plot are informative and return the object", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  km <- km_curve(sig, se)
  expect_output(print(km), "Log-rank")
  expect_output(print(km), "Median survival")
  grDevices::pdf(NULL)
  res <- plot(km)
  grDevices::dev.off()
  expect_identical(res, km)
})
