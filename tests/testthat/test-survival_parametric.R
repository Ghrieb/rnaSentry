library(SummarizedExperiment)

test_that("survival_parametric validates its inputs", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_error(survival_parametric(list(genes = "a"), se),
               "rnaSentry_signature")
  expect_error(survival_parametric(sig, matrix(1, 4, 4)), "SummarizedExperiment")
  bad <- sig
  bad$genes <- c(bad$genes, "not_a_gene")
  expect_error(survival_parametric(bad, se), "not_a_gene")
  sig$coefficients <- stats::setNames(rep(0, length(sig$genes)), sig$genes)
  expect_error(survival_parametric(sig, se), "constant")
})

test_that("survival_parametric returns the documented structure", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  sp <- survival_parametric(sig, se)
  expect_s3_class(sp, "rnaSentry_parametric")
  expect_setequal(names(sp$fits),
                  c("weibull", "exponential", "lognormal", "loglogistic"))
  expect_true(all(c("dist", "loglik", "npar", "AIC", "delta_AIC",
                    "weight", "best") %in% colnames(sp$table)))
  expect_true(sp$best %in% sp$table$dist)
  expect_true(all(sp$table$best == (sp$table$AIC == min(sp$table$AIC))))
  expect_true(all(sp$table$delta_AIC >= 0))
  expect_true(abs(sum(sp$table$weight) - 1) < 1e-6)
  expect_true(all(sp$table$AIC == sort(sp$table$AIC)))
  expect_s3_class(sp$km_fit, "survfit")
  expect_setequal(names(sp$curves), names(sp$fits))
  invisible(lapply(sp$curves, function(cr) {
    expect_true(all(c("time", "survival") %in% colnames(cr)))
    expect_true(all(cr$survival >= 0 & cr$survival <= 1))
  }))
  expect_equal(length(sp$score), ncol(se))
  expect_false(sp$sig_locked)
  expect_true(nzchar(sp$created))
})

test_that("survival_parametric score matches the km_curve risk score", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  sp <- survival_parametric(sig, se)
  km <- km_curve(sig, se)
  expect_equal(unname(sp$score), unname(km$score))
})

test_that("exponential data is best explained by the exponential model", {
  set.seed(31)
  n <- 80
  x <- stats::rnorm(n)
  t <- stats::rexp(n, rate = 0.05 * exp(0.5 * x))
  c <- stats::rexp(n, rate = 0.02)
  time <- pmin(t, c)
  event <- as.integer(t < c)
  d <- data.frame(time = time, event = event, score = x)
  fit_exp <- survival::survreg(survival::Surv(time, event) ~ score, data = d,
                               dist = "exponential")
  fit_wei <- survival::survreg(survival::Surv(time, event) ~ score, data = d,
                               dist = "weibull")
  expect_lte(stats::AIC(fit_exp), stats::AIC(fit_wei) + 2)
})

test_that("survival_parametric records the lock state of the signature", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_false(survival_parametric(sig, se)$sig_locked)
  expect_true(survival_parametric(lock_signature(sig), se)$sig_locked)
})

test_that("print and plot are informative and return the object", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  sp <- survival_parametric(sig, se)
  expect_output(print(sp), "Best model")
  expect_output(print(sp), "AIC comparison")
  grDevices::pdf(NULL)
  res <- plot(sp)
  grDevices::dev.off()
  expect_identical(res, sp)
})
