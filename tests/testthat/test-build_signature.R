library(SummarizedExperiment)

test_that("build_signature rejects non-SummarizedExperiment input", {
  expect_error(build_signature(matrix(1, 4, 4), time_col = "t",
                               event_col = "e"), "SummarizedExperiment")
})

test_that("build_signature validates time and event columns", {
  se <- make_survival_se()
  expect_error(build_signature(se, time_col = "nope", event_col = "event"),
               "no column")
  expect_error(build_signature(se, time_col = "time", event_col = "nope"),
               "no column")

  cd <- SummarizedExperiment::colData(se)
  cd$bad_time <- cd$time
  cd$bad_time[1] <- NA_real_
  SummarizedExperiment::colData(se) <- cd
  expect_error(build_signature(se, time_col = "bad_time", event_col = "event"),
               "non-negative")
  cd$neg <- cd$time
  cd$neg[1] <- -1
  SummarizedExperiment::colData(se) <- cd
  expect_error(build_signature(se, time_col = "neg", event_col = "event"),
               "non-negative")

  cd$bad_event <- rep(0L, ncol(se))
  SummarizedExperiment::colData(se) <- cd
  expect_error(build_signature(se, time_col = "time", event_col = "bad_event"),
               "at least two events")

  cd$two_level <- ifelse(cd$event == 1, 2L, 0L)
  SummarizedExperiment::colData(se) <- cd
  expect_error(build_signature(se, time_col = "time", event_col = "two_level"),
               "0/1")
})

test_that("build_signature validates its options", {
  se <- make_survival_se()
  expect_error(build_signature(se, "time", "event", outcome_col = c("a", "b")),
               "single character")
  expect_error(build_signature(se, "time", "event", method = "bogus"),
               "should be one of")
  expect_error(build_signature(se, "time", "event", top_n = 0),
               "positive integer")
  expect_error(build_signature(se, "time", "event", p_threshold = 1),
               "single number in")
  expect_error(build_signature(se, "time", "event", repeats = 0),
               "positive integer")
  expect_error(build_signature(se, "time", "event", folds = 1.5),
               "positive integer")
  expect_error(build_signature(se, "time", "event", design_terms = "nope"),
               "no design term")
})

test_that("build_signature returns the documented signature shape", {
  se <- make_survival_se()
  sig <- build_signature(se, time_col = "time", event_col = "event",
                         top_n = 20, repeats = 2, folds = 3, seed = 1)
  expect_s3_class(sig, "rnaSentry_signature")
  expect_true(all(c("genes", "outcome", "time_col", "event_col",
                    "design_terms", "cox_stats", "coefficients", "cv_results",
                    "cv_summary", "selection", "locked", "created",
                    "lock_time", "flags") %in% names(sig)))
  expect_equal(sig$outcome, "overall_survival")
  expect_lte(length(sig$genes), 20)
  expect_gte(length(sig$genes), 1)
  expect_true(all(names(sig$coefficients) == sig$genes))
  expect_true(all(c("gene", "HR", "p", "adj_p") %in%
                    colnames(sig$cox_stats)))
  expect_equal(sig$selection$method, "top_n")
  expect_equal(nrow(sig$cv_results), 2 * 3)
  expect_true(all(sig$cv_results$c_index >= 0 &
                    sig$cv_results$c_index <= 1))
  expect_true(is.finite(sig$cv_summary[["mean"]]))
  expect_false(sig$locked)
  expect_null(sig$lock_time)
  expect_true(nzchar(sig$created))
})

test_that("build_signature supports adjusted-p-value selection", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", method = "p_value",
                         p_threshold = 0.5, repeats = 2, folds = 3, seed = 2)
  expect_equal(sig$selection$method, "p_value")
  expect_gte(length(sig$genes), 1)
})

test_that("build_signature flags fewer genes than requested", {
  set.seed(11)
  counts <- matrix(stats::rpois(4 * 15, lambda = 300), nrow = 4, ncol = 15,
                   dimnames = list(paste0("g", 1:4), paste0("S", 1:15)))
  evt <- stats::rexp(15, rate = 0.03)
  cens <- stats::rexp(15, rate = 0.02)
  time <- pmin(evt, cens)
  event <- as.integer(evt < cens)
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(time = time, event = event,
                                    row.names = colnames(counts)))
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 1,
                         folds = 2, seed = 3)
  expect_true(any(sig$flags$check == "fewer_genes_than_requested"))
  expect_lte(length(sig$genes), 4)
})

test_that("build_signature drops collinear genes with a flag", {
  set.seed(12)
  counts <- matrix(stats::rpois(2 * 20, lambda = 300), nrow = 2, ncol = 20,
                   dimnames = list(c("g1", "g2"), paste0("S", 1:20)))
  counts[1, ] <- counts[2, ]
  evt <- stats::rexp(20, rate = 0.03)
  cens <- stats::rexp(20, rate = 0.02)
  time <- pmin(evt, cens)
  event <- as.integer(evt < cens)
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(time = time, event = event,
                                    row.names = colnames(counts)))
  sig <- build_signature(se, "time", "event", top_n = 2, repeats = 1,
                         folds = 2, seed = 4)
  expect_lte(length(sig$genes), 1)
  expect_true(any(sig$flags$check == "coefficient_unstable"))
})

test_that("build_signature is reproducible under a fixed seed", {
  se <- make_survival_se()
  a <- build_signature(se, "time", "event", top_n = 5, repeats = 2, folds = 3,
                       seed = 42)
  b <- build_signature(se, "time", "event", top_n = 5, repeats = 2, folds = 3,
                       seed = 42)
  expect_identical(a$cv_results, b$cv_results)
  expect_identical(a$coefficients, b$coefficients)
})

test_that("print.rnaSentry_signature is informative", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                         folds = 2, seed = 5)
  expect_output(print(sig), "signature")
})
