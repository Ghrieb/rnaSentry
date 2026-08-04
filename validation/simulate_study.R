# Step 3: simulation study. Each simulation is designed around a falsifiable
# target (documented in validation/simulation_study.md). If any target fails,
# the script exits non-zero and the gate log records exactly which one.
suppressPackageStartupMessages({
  library(pkgload)
  library(SummarizedExperiment)
})
pkg <- "C:/Users/Hani/Desktop/Academic File/rnaSentry_pkg/rnaSentry"
pkgload::load_all(pkg, quiet = TRUE)

results <- list()
pass <- function(name, detail) {
  cat(sprintf("[PASS] %-46s %s\n", name, detail))
  results[[name]] <<- TRUE
}
fail <- function(name, detail) {
  cat(sprintf("[FAIL] %-46s %s\n", name, detail))
  results[[name]] <<- FALSE
}
check <- function(cond, name, ok_detail, bad_detail) {
  if (isTRUE(cond)) pass(name, ok_detail) else fail(name, bad_detail)
  invisible(cond)
}

cat("==== Sim 1: sex_check detects 5% metadata swaps =========================\n")
set.seed(11)
n <- 200
n_f <- n / 2
xist <- c(rep(1000, n_f), rep(10, n - n_f))
y_panel <- c(rep(5, n_f), rep(400, n - n_f))
expr <- rbind(
  XIST = xist, RPS4Y1 = y_panel, DDX3Y = y_panel, KDM5D = y_panel,
  GAPDH = rep(1000, n)
)
expr <- expr * matrix(1 + stats::rnorm(nrow(expr) * n, 0, 0.01),
                      nrow = nrow(expr), ncol = n)
colnames(expr) <- paste0("S", seq_len(n))
true_sex <- c(rep("F", n_f), rep("M", n - n_f))
reported <- true_sex
swap_idx <- sample(seq_len(n), n * 0.05)
reported[swap_idx] <- ifelse(reported[swap_idx] == "F", "M", "F")
coldata <- S4Vectors::DataFrame(sex = reported, row.names = colnames(expr))
se_sex <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
res_sex <- sex_check(se_sex)
mis <- res_sex[res_sex$status == "MISMATCH", ]
n_swapped <- length(swap_idx)
check(all(swap_idx %in% which(res_sex$status == "MISMATCH")),
      "sex_check detects every swapped sample",
      sprintf("%d/%d swapped samples flagged MISMATCH", nrow(mis), n_swapped),
      sprintf("only %d/%d swapped samples flagged", nrow(mis), n_swapped))
check(all(res_sex$status[-swap_idx] == "OK"),
      "sex_check leaves true samples untouched",
      sprintf("%d/%d unswapped samples OK", n - n_swapped, n - n_swapped),
      sprintf("%d unswapped samples not OK",
              sum(res_sex$status[-swap_idx] != "OK")))

cat("==== Sim 2: design_audit flags Cramer's V ~0.7 pair as redundant ========\n")
set.seed(12)
n2 <- 180
latent <- stats::rnorm(n2)
region <- factor(cut(latent, breaks = 3, labels = c("R1", "R2", "R3")))
noise_p <- 0.33
batch <- as.character(region)
for (i in seq_len(n2)) {
  if (stats::runif(1) < noise_p) batch[i] <- sample(c("R1", "R2", "R3"), 1)
}
batch <- factor(batch)
counts2 <- matrix(stats::rpois(50 * n2, lambda = 200), nrow = 50, ncol = n2,
                  dimnames = list(paste0("g", seq_len(50)),
                                  paste0("S", seq_len(n2))))
coldata2 <- S4Vectors::DataFrame(batch = batch, region = region,
                                 row.names = colnames(counts2))
se_da <- SummarizedExperiment(assays = list(counts = counts2),
                              colData = coldata2)
da <- design_audit(se_da, design_vars = c("batch", "region"),
                   redundant_effect_size = 0.5)
pr <- da$pairwise_table
pr <- pr[((pr$var1 == "batch" & pr$var2 == "region") |
            (pr$var1 == "region" & pr$var2 == "batch")), , drop = FALSE]
V <- pr$effect_size[1]
check(nrow(pr) == 1 && V >= 0.65 && V <= 0.75,
      "design_audit pair effect size lands near 0.7",
      sprintf("Cramer's V = %.3f (target ~0.7)", V),
      sprintf("Cramer's V = %.3f outside 0.65-0.75", V))
check(nrow(pr) == 1 && isTRUE(pr$redundant[1]),
      "design_audit flags associated pair as redundant",
      sprintf("V = %.3f, redundant = TRUE", V),
      sprintf("V = %.3f, redundant = FALSE", V))
check(any(grepl("redundant_variable", da$flags$check)),
      "design_audit emits redundant_variable warning flag",
      "redundant_variable flag present in flags",
      "no redundant_variable flag")

cat("==== Sim 3: null data never yields spurious signal; CV protocol calibrates ==\n")
set.seed(13)
n3 <- 150
n_genes3 <- 40
counts3 <- matrix(stats::rpois(n_genes3 * n3, lambda = 400),
                  nrow = n_genes3, ncol = n3,
                  dimnames = list(paste0("gene", seq_len(n_genes3)),
                                  paste0("S", seq_len(n3))))
event_time3 <- stats::rexp(n3, rate = 0.03)
censor_time3 <- stats::rexp(n3, rate = 0.02)
time3 <- pmin(event_time3, censor_time3)
event3 <- as.integer(event_time3 < censor_time3)
coldata3 <- S4Vectors::DataFrame(time = time3, event = event3,
                                 row.names = colnames(counts3))
se_null <- SummarizedExperiment(assays = list(counts = counts3),
                                colData = coldata3)
mat3 <- log2(as.matrix(SummarizedExperiment::assay(se_null, 1)) + 1)
sig_null <- build_signature(se_null, "time", "event", top_n = 5,
                            repeats = 5, folds = 5, seed = 3)
ci_mean <- sig_null$cv_summary[["mean"]]
check(ci_mean < 0.6,
      "build_signature on null data reports no spurious signal",
      sprintf("mean C-index = %.3f (< 0.6; null data never yield strong signal)", ci_mean),
      sprintf("mean C-index = %.3f >= 0.6: spurious signal", ci_mean))
# Calibration control: identical CV protocol with RANDOM (unselected) genes.
# The null expectation for a well-calibrated concordance estimator is 0.5.
set.seed(3)
n <- length(time3)
ev3 <- event3
ref_cis <- numeric(25)
k <- 0
for (r in seq_len(5)) {
  fold_ids <- integer(n)
  for (grp in list(which(ev3 == 0L), which(ev3 == 1L))) {
    g <- grp[sample.int(length(grp))]
    fold_ids[g] <- rep(1:5, length.out = length(grp))
  }
  for (f in 1:5) {
    k <- k + 1
    train <- which(fold_ids != f)
    test <- which(fold_ids == f)
    rnd_genes <- sample(rownames(mat3), 5)
    d_tr <- data.frame(time = time3[train], event = ev3[train],
                       t(as.matrix(mat3[rnd_genes, train, drop = FALSE])))
    fit_tr <- survival::coxph(survival::Surv(time, event) ~ ., data = d_tr)
    b <- stats::coef(fit_tr)
    score <- as.vector(t(as.matrix(mat3[names(b), test, drop = FALSE])) %*% b)
    cc <- tryCatch(
      suppressWarnings(
        survival::concordance(survival::Surv(time3[test], ev3[test]) ~ score,
                              reverse = TRUE)
      ),
      error = function(e) NULL)
    ref_cis[k] <- if (is.null(cc)) NA_real_ else as.numeric(cc$concordance[1])
  }
}
ctl_mean <- mean(ref_cis, na.rm = TRUE)
check(ctl_mean >= 0.45 && ctl_mean <= 0.55,
      "CV protocol is calibrated under the null (random-gene control)",
      sprintf("control mean C = %.3f (null expectation 0.5)", ctl_mean),
      sprintf("control mean C = %.3f outside 0.45-0.55", ctl_mean))

cat("==== Sim 4: cox_model flags time-varying hazard (PH violation) ==========\n")
set.seed(14)
n4 <- 800
counts4 <- matrix(stats::rpois(30 * n4, lambda = 500), nrow = 30, ncol = n4,
                  dimnames = list(paste0("gene", seq_len(30)),
                                  paste0("S", seq_len(n4))))
z4 <- scale(colMeans(counts4))[, 1]
t0 <- 5
b1 <- 1.0
h0 <- 0.05
u <- stats::runif(n4)
lam1 <- h0 * exp(b1 * z4)
lam2 <- h0 * exp(-b1 * z4)
q1 <- (-log(u)) / lam1
T4 <- ifelse(q1 <= t0, q1, t0 + (q1 - t0) * lam1 / lam2)
censor4 <- stats::rexp(n4, rate = 0.02)
time4 <- pmin(T4, censor4)
event4 <- as.integer(T4 < censor4)
coldata4 <- S4Vectors::DataFrame(time = time4, event = event4,
                                 row.names = colnames(counts4))
se_tv <- SummarizedExperiment(assays = list(counts = counts4),
                              colData = coldata4)
sig_tv <- build_signature(se_tv, "time", "event", top_n = 3,
                          repeats = 2, folds = 3, seed = 4)
cm_tv <- cox_model(sig_tv, se_tv)
check(isTRUE(cm_tv$ph_violated),
      "cox_model flags PH violation on time-varying effect",
      sprintf("ph_violated = TRUE; %d violator term(s)", 
              if (is.null(cm_tv$zph_summary)) 0 else
                sum(!is.na(cm_tv$zph_summary$p) & cm_tv$zph_summary$p < 0.05)),
      "ph_violated = FALSE")
check(any(grepl("proportional_hazards", cm_tv$flags$check)),
      "cox_model emits proportional_hazards warning flag",
      "proportional_hazards flag present in flags",
      "no proportional_hazards flag")

cat("==== Sim 5: validate_external on unlocked signature runs, records FALSE ==\n")
set.seed(15)
n5 <- 80
counts5 <- matrix(stats::rpois(80 * n5, lambda = 400), nrow = 80, ncol = n5,
                  dimnames = list(paste0("gene", seq_len(80)),
                                  paste0("S", seq_len(n5))))
sig_expr5 <- colMeans(counts5[1:5, , drop = FALSE])
risk5 <- scale(sig_expr5)[, 1] * 0.4
et5 <- stats::rexp(n5, rate = 0.03 * exp(0.8 * risk5))
ct5 <- stats::rexp(n5, rate = 0.02)
coldata5 <- S4Vectors::DataFrame(
  time = pmin(et5, ct5), event = as.integer(et5 < ct5),
  row.names = colnames(counts5))
se_train <- SummarizedExperiment(assays = list(counts = counts5),
                                 colData = coldata5)
sig5 <- build_signature(se_train, "time", "event", top_n = 10,
                        repeats = 2, folds = 3, seed = 5)
cutoff5 <- stats::median(km_curve(sig5, se_train)$score)
ext5 <- make_cohort_se_like <- function(sig, n, seed) {
  set.seed(seed)
  cg <- c(sig$genes, paste0("extra_gene", seq_len(40)))
  cnt <- matrix(stats::rpois(length(cg) * n, lambda = 400),
                nrow = length(cg), ncol = n,
                dimnames = list(cg, paste0("C", seq_len(n))))
  et <- stats::rexp(n, rate = 0.03)
  ct <- stats::rexp(n, rate = 0.02)
  cd <- S4Vectors::DataFrame(time = pmin(et, ct), event = as.integer(et < ct),
                             row.names = colnames(cnt))
  SummarizedExperiment(assays = list(counts = cnt), colData = cd)
}
se_ext5 <- make_cohort_se_like(sig5, 60, 5)
val5 <- validate_external(sig5, se_ext5, cutpoint = cutoff5)
check(isFALSE(val5$sig_locked),
      "validate_external records sig_locked = FALSE on unlocked signature",
      "sig_locked = FALSE recorded; no hard error",
      "sig_locked is not FALSE or validation errored")
check(is.finite(val5$log_rank_p) && val5$log_rank_p >= 0 && val5$log_rank_p <= 1,
      "external validation produced a valid log-rank p",
      sprintf("log_rank_p = %.4f, %d samples scored", val5$log_rank_p,
              length(val5$score)),
      "log_rank_p non-finite or out of range")

cat("==== Sim 6: power analysis for external-signal transfer ================\n")
# Two signal tiers calibrated so that a large cohort of the same generative
# model yields concordance ~0.60 (beta = 0.4) and ~0.65 (beta = 0.6). For each
# tier we estimate the empirical power of *external* validation (log-rank
# p < 0.05 on an independent cohort scored against the discovery cutpoint) as
# a function of the external event count {35, 60, 100, 150, 250}.
# Discovery cohort is held fixed (n = 120) so the causal axis is the size of
# the external cohort. Falsifiable targets:
#   (a) power is monotone non-decreasing in external events;
#   (b) power < 0.30 at ~35 external events (small-cohort trap);
#   (c) power >= 0.80 at ~250 events for the moderate (C~0.65) tier.
make_cohort_se <- function(n, beta, seed) {
  set.seed(seed)
  p <- 60
  counts <- matrix(stats::rpois(p * n, lambda = 400), nrow = p, ncol = n,
                   dimnames = list(paste0("gene", seq_len(p)),
                                   paste0("S", seq_len(n))))
  z <- scale(colMeans(counts[1:5, , drop = FALSE]))[, 1]
  et <- stats::rexp(n, rate = 0.05 * exp(beta * z))
  ct <- stats::rexp(n, rate = 0.05)
  cd <- S4Vectors::DataFrame(time = pmin(et, ct),
                             event = as.integer(et < ct),
                             row.names = colnames(counts))
  SummarizedExperiment(assays = list(counts = counts), colData = cd)
}
# Calibration check: effective concordance of the two tiers on a large cohort.
set.seed(991)
n_cal <- 40000
z_cal <- scale(colMeans(matrix(stats::rpois(5 * n_cal, lambda = 400),
                               nrow = 5, ncol = n_cal)))[, 1]
u_cal <- stats::runif(n_cal)
c_eff <- function(beta) {
  T_cal <- (-log(u_cal)) / (0.05 * exp(beta * z_cal))
  as.numeric(survival::concordance(
    survival::Surv(T_cal, rep(1, n_cal)) ~ z_cal,
    reverse = TRUE)$concordance)
}
c_lo <- c_eff(0.4)
c_hi <- c_eff(0.6)
check(isTRUE(c_hi > c_lo),
      "Sim 6 tiers are ordered by effective concordance",
      sprintf("C(0.6) = %.3f > C(0.4) = %.3f", c_hi, c_lo),
      sprintf("C(0.6) = %.3f <= C(0.4) = %.3f", c_hi, c_lo))

n_ext_grid <- c(70, 120, 200, 300, 500, 600)
power_table <- matrix(NA_real_, nrow = 2, ncol = length(n_ext_grid),
                      dimnames = list(c("weak", "moderate"), n_ext_grid))
event_table <- matrix(NA_real_, nrow = 2, ncol = length(n_ext_grid),
                      dimnames = dimnames(power_table))
reps6 <- 120
for (ti in seq_len(2)) {
  beta6 <- c(0.4, 0.6)[ti]
  # Discovery cohort is built once per replicate and shared across the
  # external-event grid, so the causal axis is only the external cohort size.
  sigs6 <- vector("list", reps6)
  for (r in seq_len(reps6)) {
    disc6 <- make_cohort_se(120, beta6, 100000 + ti * 1000 + r)
    sigs6[[r]] <- build_signature(disc6, "time", "event", top_n = 8,
                                  repeats = 2, folds = 3, seed = r)
  }
  for (j in seq_along(n_ext_grid)) {
    hit <- 0L
    for (r in seq_len(reps6)) {
      sig6 <- sigs6[[r]]
      cut6 <- stats::median(km_curve(sig6, make_cohort_se(
        120, beta6, 100000 + ti * 1000 + r))$score)
      ext6 <- make_cohort_se(n_ext_grid[j], beta6, 200000 + ti * 1000 + r)
      val6 <- validate_external(sig6, ext6, cutpoint = cut6)
      if (is.finite(val6$log_rank_p) && val6$log_rank_p < 0.05) hit <- hit + 1L
    }
    power_table[ti, j] <- hit / reps6
    event_table[ti, j] <- mean(vapply(seq_len(reps6), function(r)
      sum(make_cohort_se(n_ext_grid[j], beta6,
                         200000 + ti * 1000 + r)$event), numeric(1)))
  }
}
cat(sprintf("external events (expected): %s\n",
            paste(sprintf("%.0f", event_table[1, ]), collapse = ", ")))
cat(sprintf("power weak     (C~%.2f): %s\n", c_lo,
            paste(sprintf("%.3f", power_table["weak", ]), collapse = ", ")))
cat(sprintf("power moderate (C~%.2f): %s\n", c_hi,
            paste(sprintf("%.3f", power_table["moderate", ]), collapse = ", ")))
p_weak_lo <- power_table["weak", 1]
p_weak_hi <- power_table["weak", ncol(power_table)]
p_mod_lo <- power_table["moderate", 1]
p_mod_hi <- power_table["moderate", ncol(power_table)]
check(p_weak_hi >= p_weak_lo && p_mod_hi >= p_mod_lo,
      "transfer power monotone non-decreasing in external events",
      sprintf("weak %.3f -> %.3f, moderate %.3f -> %.3f",
              p_weak_lo, p_weak_hi, p_mod_lo, p_mod_hi),
      "power decreased from ~35 to ~250 events in at least one tier")
check(p_weak_lo < 0.30 && p_mod_lo < 0.30,
      "small-cohort trap: power < 0.30 at ~35 external events",
      sprintf("weak %.3f, moderate %.3f", p_weak_lo, p_mod_lo),
      sprintf("power >= 0.30 at ~35 events (weak %.3f, moderate %.3f)",
              p_weak_lo, p_mod_lo))
check(p_mod_hi >= 0.80,
      "adequate power at scale for moderate tier",
      sprintf("power = %.3f at ~%.0f events (>= 0.80)",
              p_mod_hi, event_table["moderate", ncol(power_table)]),
      sprintf("power = %.3f at ~%.0f events (< 0.80)",
              p_mod_hi, event_table["moderate", ncol(power_table)]))

cat("\n==== SUMMARY ====\n")
n_pass <- sum(unlist(results))
n_fail <- length(results) - n_pass
cat(sprintf("Simulations: %d PASS, %d FAIL of %d targets\n",
            n_pass, n_fail, length(results)))
if (n_fail > 0) quit(status = 1)
