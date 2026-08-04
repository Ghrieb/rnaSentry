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
  expect_error(build_signature(se, "time", "event",
                               min_events_per_parameter = 0),
               "positive number")
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

test_that("build_signature validates adjust_for_design and design terms", {
  se <- make_survival_se()
  expect_error(build_signature(se, "time", "event", design_terms = "batch",
                               adjust_for_design = "yes"),
               "single TRUE or FALSE")
  SummarizedExperiment::colData(se)$constant <- rep("x", ncol(se))
  expect_error(build_signature(se, "time", "event", design_terms = "constant"),
               "no variation")
})

test_that("design-adjusted screening excludes batch-driven genes", {
  set.seed(77)
  n_g <- 120
  n_s <- 60
  counts <- matrix(stats::rpois(n_g * n_s, lambda = 300), nrow = n_g,
                   ncol = n_s,
                   dimnames = list(paste0("gene", seq_len(n_g)),
                                   paste0("S", seq_len(n_s))))
  batch <- factor(rep(c("B1", "B2"), each = n_s / 2))
  counts[seq_len(30), batch == "B2"] <- counts[seq_len(30), batch == "B2"] * 3L
  rate <- ifelse(batch == "B2", 0.08, 0.02)
  evt <- stats::rexp(n_s, rate = rate)
  cens <- stats::rexp(n_s, rate = 0.02)
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      time = pmin(evt, cens), event = as.integer(evt < cens),
      batch = batch, row.names = colnames(counts)))
  block <- paste0("gene", 1:30)

  unadj <- build_signature(se, "time", "event", top_n = 3, repeats = 1,
                           folds = 2, seed = 1, design_terms = "batch",
                           adjust_for_design = FALSE)
  expect_true(any(unadj$flags$check == "unadjusted_screening"))
  expect_equal(unadj$screening_terms, character(0))
  p_block <- unadj$cox_stats$p[unadj$cox_stats$gene %in% block]
  p_rest <- unadj$cox_stats$p[!(unadj$cox_stats$gene %in% block)]
  expect_lt(stats::median(p_block, na.rm = TRUE),
            stats::median(p_rest, na.rm = TRUE))
  expect_lt(stats::median(p_block, na.rm = TRUE), 0.01)

  adj <- build_signature(se, "time", "event", top_n = 3, repeats = 1,
                         folds = 2, seed = 1, design_terms = "batch")
  expect_true(any(adj$flags$check == "adjusted_screening"))
  expect_equal(adj$screening_terms, "batch")
  adj_p_block <- adj$cox_stats$p[adj$cox_stats$gene %in% block]
  adj_p_rest <- adj$cox_stats$p[!(adj$cox_stats$gene %in% block)]
  expect_gt(stats::median(adj_p_block, na.rm = TRUE), 0.1)
  expect_gt(stats::median(adj_p_rest, na.rm = TRUE), 0.1)
})

test_that("print.rnaSentry_signature is informative", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                         folds = 2, seed = 5)
  expect_output(print(sig), "signature")
})

test_that("build_signature flags an already-log-scaled count assay", {
  set.seed(21)
  n_g <- 40
  n_s <- 30
  counts <- matrix(stats::rpois(n_g * n_s, lambda = 400), nrow = n_g,
                   ncol = n_s,
                   dimnames = list(paste0("g", seq_len(n_g)),
                                   paste0("S", seq_len(n_s))))
  prelogged <- log2(counts + 1)
  event_time <- stats::rexp(n_s, rate = 0.03)
  censor_time <- stats::rexp(n_s, rate = 0.02)
  se <- SummarizedExperiment(
    assays = list(counts = prelogged),
    colData = S4Vectors::DataFrame(
      time = pmin(event_time, censor_time),
      event = as.integer(event_time < censor_time),
      row.names = colnames(counts)))
  # min_events_per_parameter = 1 isolates the scale warning from the EPV
  # warning so the test targets exactly the guardrail under test.
  expect_warning(
    sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                           folds = 2, seed = 1,
                           min_events_per_parameter = 1),
    "already log-transformed")
  expect_true(any(sig$flags$check == "possibly_log_scaled"))
  expect_equal(sig$flags$severity[sig$flags$check == "possibly_log_scaled"],
               "warning")
})

test_that("build_signature does not flag raw integer counts", {
  se <- make_survival_se()
  expect_warning(
    sig <- build_signature(se, "time", "event", top_n = 3, repeats = 1,
                           folds = 2, seed = 6,
                           min_events_per_parameter = 1),
    NA)
  expect_false(any(sig$flags$check == "possibly_log_scaled"))
})

test_that("build_signature flags few events per parameter", {
  set.seed(22)
  n_g <- 60
  n_s <- 40
  counts <- matrix(stats::rpois(n_g * n_s, lambda = 300), nrow = n_g,
                   ncol = n_s,
                   dimnames = list(paste0("g", seq_len(n_g)),
                                   paste0("S", seq_len(n_s))))
  event_time <- stats::rexp(n_s, rate = 0.01)
  censor_time <- stats::rexp(n_s, rate = 0.04)
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(
      time = pmin(event_time, censor_time),
      event = as.integer(event_time < censor_time),
      row.names = colnames(counts)))
  n_events <- sum(SummarizedExperiment::colData(se)$event)
  expect_lt(n_events, 40)
  expect_warning(
    sig <- build_signature(se, "time", "event", top_n = 8, repeats = 1,
                           folds = 2, seed = 2),
    "events per parameter")
  expect_true(any(sig$flags$check == "events_per_parameter"))
  expect_equal(sig$flags$severity[sig$flags$check == "events_per_parameter"],
               "warning")
})

test_that("build_signature does not flag adequate event counts", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 3, repeats = 1,
                         folds = 2, seed = 6)
  expect_false(any(sig$flags$check == "events_per_parameter"))
  expect_false(any(sig$flags$check == "possibly_log_scaled"))
})
