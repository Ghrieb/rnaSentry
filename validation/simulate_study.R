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

cat("\n==== SUMMARY ====\n")
n_pass <- sum(unlist(results))
n_fail <- length(results) - n_pass
cat(sprintf("Simulations: %d PASS, %d FAIL of %d targets\n",
            n_pass, n_fail, length(results)))
if (n_fail > 0) quit(status = 1)
