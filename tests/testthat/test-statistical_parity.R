library(SummarizedExperiment)

# Step 2: statistical parity tests. Each statistic that rnaSentry computes
# itself is compared against the reference implementation from the survival /
# stats / base R toolchain on the same data. These pin the statistical
# contract, not just the output shape.

test_that("km_curve log-rank p matches survival::survdiff on the same split", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  km <- km_curve(sig, se)
  d <- data.frame(time = SummarizedExperiment::colData(se)$time,
                  event = SummarizedExperiment::colData(se)$event,
                  group = km$groups)
  dd <- survival::survdiff(survival::Surv(time, event) ~ group, data = d)
  ref_p <- stats::pchisq(dd$chisq, df = max(1, length(dd$n) - 1),
                         lower.tail = FALSE)
  expect_equal(km$log_rank_p, unname(ref_p))
})

test_that("validate_external log-rank p matches survival::survdiff", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  ext <- make_cohort_se(sig, seed = 202)
  ext_cutoff <- stats::median(km_curve(sig, ext)$score)
  val <- validate_external(sig, ext, cutpoint = ext_cutoff)
  d <- data.frame(time = SummarizedExperiment::colData(ext)$time,
                  event = SummarizedExperiment::colData(ext)$event,
                  group = val$groups)
  dd <- survival::survdiff(survival::Surv(time, event) ~ group, data = d)
  ref_p <- stats::pchisq(dd$chisq, df = max(1, length(dd$n) - 1),
                         lower.tail = FALSE)
  expect_equal(val$log_rank_p, unname(ref_p))
})

test_that("design_audit Cramer's V matches the textbook formula", {
  se <- make_pca_se()
  se$region <- factor(rep(c("R1", "R2"), length.out = ncol(se)))
  res <- design_audit(se, design_vars = c("batch", "region"))
  pr <- res$pairwise_table[
    ((res$pairwise_table$var1 == "batch" &
        res$pairwise_table$var2 == "region") |
       (res$pairwise_table$var1 == "region" &
          res$pairwise_table$var2 == "batch")), , drop = FALSE]
  expect_equal(nrow(pr), 1)
  x1 <- factor(SummarizedExperiment::colData(se)$batch)
  x2 <- factor(SummarizedExperiment::colData(se)$region)
  tab <- table(x1, x2)
  chi <- unname(stats::chisq.test(tab)$statistic)
  n <- sum(tab)
  k <- min(dim(tab)) - 1
  ref_v <- sqrt((chi / n) / k)
  expect_equal(pr$effect_size, ref_v, tolerance = 1e-8)
  expect_equal(pr$effect_size_type, "cramers_v")
})

test_that("cox_model hazard ratios and CIs match an independent survival::coxph", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cm <- cox_model(sig, se, confounders = "age")

  mat <- log2(as.matrix(SummarizedExperiment::assay(se, 1)) + 1)
  expr <- as.data.frame(t(mat[sig$genes, , drop = FALSE]))
  d <- cbind(data.frame(time = SummarizedExperiment::colData(se)$time,
                        event = SummarizedExperiment::colData(se)$event),
             expr)
  d$age <- SummarizedExperiment::colData(se)$age
  fit_ref <- survival::coxph(survival::Surv(time, event) ~ ., data = d)
  sm_ref <- summary(fit_ref)$coefficients
  ci_ref <- exp(stats::confint(fit_ref))

  expect_equal(cm$coef_table$coefficient,
               unname(sm_ref[cm$coef_table$term, "coef"]), tolerance = 1e-8)
  expect_equal(cm$coef_table$HR,
               unname(sm_ref[cm$coef_table$term, "exp(coef)"]), tolerance = 1e-8)
  expect_equal(cm$coef_table$HR_low, unname(ci_ref[cm$coef_table$term, 1]),
               tolerance = 1e-8)
  expect_equal(cm$coef_table$HR_high, unname(ci_ref[cm$coef_table$term, 2]),
               tolerance = 1e-8)
  expect_equal(cm$coef_table$p,
               unname(sm_ref[cm$coef_table$term, "Pr(>|z|)"]), tolerance = 1e-8)
})

test_that("cox_model Schoenfeld table matches survival::cox.zph", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cm <- cox_model(sig, se)
  zph_ref <- survival::cox.zph(cm$fit)
  zt <- zph_ref$table
  pcol <- intersect(c("p", "Pr(>|Chi|)", "Pr(>Chisq)"), colnames(zt))[1]
  rcol <- intersect(c("rho", "rho[1]"), colnames(zt))[1]
  rho_ref <- if (is.na(rcol)) rep(NA_real_, nrow(zt)) else unname(zt[, rcol])
  term_ref <- sub("^`(.*)`$", "\\1", rownames(zt))
  expect_equal(cm$zph_summary$p, unname(zt[, pcol]), tolerance = 1e-8)
  expect_equal(cm$zph_summary$rho, rho_ref, tolerance = 1e-8)
  expect_equal(cm$zph_summary$term, term_ref)
})

test_that("cox_model concordance is the risk-score direction (independently recomputed)", {
  # cox_model() takes its concordance from summary(fit)$concordance, i.e. the
  # coxph method's own C. Verify that value independently against
  # survival::concordance() on the fitted linear predictor with reverse = TRUE
  # (higher linear predictor = shorter survival), so a future direction flip in
  # the model construction cannot pass silently.
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cm <- cox_model(sig, se, confounders = "age")
  lp <- as.numeric(stats::predict(cm$fit, type = "lp"))
  ref <- survival::concordance(
    survival::Surv(SummarizedExperiment::colData(se)$time,
                   SummarizedExperiment::colData(se)$event) ~ lp,
    reverse = TRUE
  )
  expect_equal(cm$concordance, as.numeric(ref$concordance[1]),
               tolerance = 1e-6)
  # Direction agreement across stages: both cox_model's in-sample concordance
  # and build_signature's held-out CV mean must beat 0.5 on the same
  # signal-bearing fixture (a sign flip would send both well below 0.5).
  expect_gt(cm$concordance, 0.5)
  expect_gt(mean(sig$cv_results$c_index, na.rm = TRUE), 0.5)
})

test_that("survival_parametric AIC/loglik/npar match stats::AIC on the same fits", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  sp <- survival_parametric(sig, se)
  invisible(lapply(names(sp$fits), function(dist) {
    expect_equal(sp$table$AIC[sp$table$dist == dist],
                 stats::AIC(sp$fits[[dist]]))
    expect_equal(sp$table$loglik[sp$table$dist == dist],
                 as.numeric(stats::logLik(sp$fits[[dist]])))
    expect_equal(sp$table$npar[sp$table$dist == dist],
                 length(sp$fits[[dist]]$coefficients) + 1)
  }))
  expect_equal(sp$best, sp$table$dist[which.min(sp$table$AIC)])
})

test_that("build_signature fold concordance matches concordance on the same folds", {
  # Re-runs the documented CV protocol (event-stratified fold seeding,
  # per-fold gene coefficients, held-out concordance) with survival::concordance
  # on the exact fold assignments the package drew, and checks every fold.
  se <- make_survival_se()
  seed <- 7
  repeats <- 1
  folds <- 3
  sig <- build_signature(se, "time", "event", top_n = 10,
                         repeats = repeats, folds = folds, seed = seed)
  time_vec <- as.numeric(SummarizedExperiment::colData(se)$time)
  event_vec <- as.integer(SummarizedExperiment::colData(se)$event)
  mat <- log2(as.matrix(SummarizedExperiment::assay(se, 1)) + 1)

  set.seed(seed)
  n <- length(time_vec)
  fold_ids <- integer(n)
  strata <- lapply(list(which(event_vec == 0L), which(event_vec == 1L)),
                   function(grp) {
                     if (length(grp) == 0) return(NULL)
                     grp_shuffled <- grp[sample.int(length(grp))]
                     data.frame(idx = grp_shuffled,
                                f = rep(seq_len(folds), length.out = length(grp)))
                   })
  strata <- do.call(rbind, strata)
  if (!is.null(strata) && nrow(strata) > 0) fold_ids[strata$idx] <- strata$f

  ref <- vapply(seq_len(folds), function(f) {
    train <- which(fold_ids != f)
    test <- which(fold_ids == f)
    if (length(train) < 5 || length(test) < 2 ||
        sum(event_vec[train]) < 2 || sum(event_vec[test]) < 1) {
      return(NA_real_)
    }
    d_tr <- data.frame(time = time_vec[train], event = event_vec[train],
                       base::t(as.matrix(mat[sig$genes, train, drop = FALSE])))
    fit_tr <- survival::coxph(survival::Surv(time, event) ~ ., data = d_tr)
    b <- stats::coef(fit_tr)
    score_test <- as.vector(base::t(as.matrix(mat[names(b), test, drop = FALSE])) %*% b)
    if (isTRUE(stats::sd(score_test) == 0)) {
      return(NA_real_)
    }
    conc <- tryCatch(
      suppressWarnings(
        survival::concordance(survival::Surv(time_vec[test], event_vec[test]) ~
                                score_test, reverse = TRUE)
      ),
      error = function(e) NULL
    )
    val <- if (is.null(conc)) NA_real_ else as.numeric(conc$concordance[1])
    if (!is.finite(val)) NA_real_ else val
  }, numeric(1))

  expect_equal(sig$cv_results$c_index, ref)
  expect_equal(sig$cv_results$fold, seq_len(folds))
  # Direction guard: the fixture plants higher expression -> shorter survival,
  # so a correctly-oriented risk score must beat 0.5 on held-out folds.
  # (A sign-flipped concordance convention would push this well below 0.5.)
  expect_gt(mean(sig$cv_results$c_index, na.rm = TRUE), 0.5)
})

test_that("survival::concordance reverse convention satisfies C + C_rev = 1", {
  # Guards the assumption behind the risk-score convention used in
  # build_signature() and validate_external(): reversing the association of a
  # numeric predictor must give a concordance exactly complementary to 1.
  set.seed(9)
  t <- round(stats::rexp(40, rate = 0.05), 1)
  e <- stats::rbinom(40, 1, 0.7)
  x <- stats::rnorm(40)
  y <- survival::Surv(t, e)
  c1 <- as.numeric(survival::concordance(y ~ x)$concordance[1])
  c2 <- as.numeric(survival::concordance(y ~ x, reverse = TRUE)$concordance[1])
  expect_equal(c1 + c2, 1, tolerance = 1e-9)
})
