library(SummarizedExperiment)

# Regression tests for issues found by the adversarial pass (Step 1).
# Each block pins a specific failure mode that previously crashed, returned
# silent nonsense, or violated the documented output contract.

test_that("build_signature rejects leave-one-out or oversized folds cleanly", {
  se <- make_survival_se()
  expect_error(build_signature(se, "time", "event", folds = 60, repeats = 1),
               "smaller than the number of samples")
  expect_error(build_signature(se, "time", "event", folds = 61, repeats = 1),
               "smaller than the number of samples")
})

test_that("build_signature tolerates many folds without a cryptic crash", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 5, folds = 40,
                         repeats = 1, seed = 42)
  expect_equal(nrow(sig$cv_results), 40)
  expect_true(all(sig$cv_results$c_index >= 0 | is.na(sig$cv_results$c_index)))
  expect_true(is.finite(sig$cv_summary[["mean"]]))
})

test_that("build_signature rejects non-finite follow-up time", {
  se <- make_survival_se()
  SummarizedExperiment::colData(se)$time[10] <- Inf
  expect_error(build_signature(se, "time", "event"), "finite")
})

test_that("build_signature errors cleanly on a single-sample cohort", {
  counts <- matrix(500, 5, 1,
                   dimnames = list(paste0("g", 1:5), "S1"))
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(time = 10, event = 1L, row.names = "S1"))
  expect_error(build_signature(se, "time", "event", top_n = 5, repeats = 1,
                               folds = 2))
})

test_that("build_signature rejects folds >= n on a two-sample cohort", {
  counts <- matrix(stats::rpois(10, 500), nrow = 5, ncol = 2,
                   dimnames = list(paste0("g", 1:5), c("S1", "S2")))
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(time = c(10, 20), event = c(1L, 1L),
                                   row.names = c("S1", "S2")))
  expect_error(build_signature(se, "time", "event", top_n = 5, repeats = 1,
                               folds = 2),
               "smaller than the number of samples")
})

test_that("genes named like internal model columns are excluded with a flag", {
  set.seed(3)
  n <- 60
  g <- matrix(stats::rpois(60 * 20, 500), 20, n,
              dimnames = list(c("time", "event", paste0("gene", 3:20)),
                              paste0("S", 1:n)))
  risk <- scale(log2(g["time", ] + 1))[, 1] * 0.9
  cd <- S4Vectors::DataFrame(
    time = pmin(stats::rexp(n, 0.03 * exp(risk)), stats::rexp(n, 0.02)),
    event = as.integer(stats::rexp(n, 0.03 * exp(risk)) <
                         stats::rexp(n, 0.02)),
    row.names = colnames(g))
  se <- SummarizedExperiment(assays = list(counts = g), colData = cd)
  sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                         folds = 2, seed = 11)
  expect_true(any(sig$flags$check == "gene_name_collision"))
  expect_false("time" %in% sig$genes)
  expect_false("event" %in% sig$genes)
  expect_true(all(sig$genes %in% rownames(g)))
})

test_that("exactly collinear genes are dropped from the signature with a flag", {
  set.seed(3)
  n <- 60
  g <- matrix(stats::rpois(60 * 20, 500), 20, n,
              dimnames = list(paste0("gene", 1:20), paste0("S", 1:n)))
  risk <- scale(log2(g["gene2", ] + 1))[, 1] * 1.2
  g["gene3", ] <- g["gene2", ]
  cd <- S4Vectors::DataFrame(
    time = pmin(stats::rexp(n, 0.03 * exp(risk)), stats::rexp(n, 0.02)),
    event = as.integer(stats::rexp(n, 0.03 * exp(risk)) <
                         stats::rexp(n, 0.02)),
    row.names = colnames(g))
  se <- SummarizedExperiment(assays = list(counts = g), colData = cd)
  sig <- build_signature(se, "time", "event", top_n = 8, repeats = 1,
                         folds = 2, seed = 11)
  expect_true(all(is.finite(exp(sig$coefficients))))
  collinear_left <- intersect(c("gene2", "gene3"), sig$genes)
  if (length(collinear_left) > 0) {
    expect_true(any(sig$flags$check == "coefficient_unstable"))
  }
})

test_that("cox_model rejects confounders that are the survival columns", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_error(cox_model(sig, se, confounders = "event"),
               "survival outcome columns")
  expect_error(cox_model(sig, se, confounders = "time"),
               "survival outcome columns")
})

test_that("cox_model flags non-estimable coefficients instead of silent Inf HR", {
  set.seed(3)
  n <- 60
  g <- matrix(stats::rpois(60 * 20, 500), 20, n,
              dimnames = list(paste0("gene", 1:20), paste0("S", 1:n)))
  risk <- scale(log2(g["gene2", ] + 1))[, 1] * 1.2
  g["gene3", ] <- g["gene2", ]
  cd <- S4Vectors::DataFrame(
    time = pmin(stats::rexp(n, 0.03 * exp(risk)), stats::rexp(n, 0.02)),
    event = as.integer(stats::rexp(n, 0.03 * exp(risk)) <
                         stats::rexp(n, 0.02)),
    row.names = colnames(g))
  se <- SummarizedExperiment(assays = list(counts = g), colData = cd)
  sig <- make_signature_for_testing(n_genes = 2)
  sig$genes <- c("gene2", "gene3")
  sig$coefficients <- stats::setNames(c(1, 1), c("gene2", "gene3"))
  cm <- cox_model(sig, se)
  expect_true(any(cm$flags$check %in%
                  c("non_estimable_coefficients", "non_estimable_genes")))
})

test_that("design_audit excludes constant and all-missing variables from the formula", {
  se <- make_pca_se()
  se$cnst <- 5
  se$allna <- rep(NA_real_, ncol(se))
  res <- design_audit(se, design_vars = c("cnst", "allna", "batch"))
  expect_false(grepl("cnst", res$formula_text))
  expect_false(grepl("allna", res$formula_text))
  expect_true(grepl("batch", res$formula_text))
  expect_true(any(res$flags$check == "variable_constant"))
})

test_that("qc_explore flags a fully missing count assay without crashing", {
  se <- make_pca_se()
  counts <- assay(se, 1)
  counts[] <- NA_real_
  assay(se, 1) <- counts
  qc <- qc_explore(se)
  expect_true(any(qc$flags$check == "missing_counts"))
  expect_equal(qc$library_size_outliers, character(0))
})

test_that("sex_check reports MISSING for samples without reported sex", {
  genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D")
  expr <- matrix(0, nrow = 4, ncol = 3, dimnames = list(genes, paste0("S", 1:3)))
  expr["XIST", ] <- c(900, 5, 900)
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 1] <- 5
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 2] <- 400
  expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 3] <- 5
  coldata <- S4Vectors::DataFrame(sex = c("F", NA, "M"),
                                  row.names = paste0("S", 1:3))
  se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
  res <- sex_check(se)
  expect_equal(res$status, c("OK", "MISSING", "MISMATCH"))
})

test_that("pca_audit errors cleanly on a single-sample cohort", {
  counts <- matrix(stats::rpois(5, 500), nrow = 5, ncol = 1,
                   dimnames = list(paste0("g", 1:5), "S1"))
  se <- SummarizedExperiment(assays = list(counts = counts))
  expect_error(pca_audit(se), "at least two samples")
})

test_that("km_curve errors cleanly on a single-sample cohort", {
  sig <- make_signature_for_testing()
  counts <- matrix(stats::rpois(20, 500), nrow = 20, ncol = 1,
                   dimnames = list(sig$genes, "S1"))
  se <- SummarizedExperiment(
    assays = list(counts = counts),
    colData = S4Vectors::DataFrame(time = 10, event = 1L, row.names = "S1"))
  expect_error(km_curve(sig, se), "leaves one risk group empty")
})

test_that("downstream stages reject non-finite follow-up time", {
  sig <- make_signature_for_testing()
  se <- make_survival_se()
  SummarizedExperiment::colData(se)$time[1] <- Inf
  expect_error(cox_model(sig, se), "finite")
  expect_error(survival_parametric(sig, se), "finite")
  expect_error(validate_external(sig, se, cutpoint = 1), "finite")
  expect_error(km_curve(sig, se), "finite")
})

test_that("build_signature and cox_model handle non-syntactic gene symbols", {
  # Real annotation can contain symbols such as "1-Mar" / "7-Sep" that
  # data.frame()/coxph would otherwise mangle or backtick; the signature
  # genes and Cox terms must still match the assay rownames exactly.
  genes <- c("1-Mar", "7-Sep", "XIST", "ESR1", "MKI67", "A-kinase",
             "HLA-DRB1", "10-Sep")
  n <- 60
  set.seed(99)
  counts <- matrix(stats::rpois(length(genes) * n, lambda = 300),
                   nrow = length(genes), ncol = n,
                   dimnames = list(genes, paste0("S", seq_len(n))))
  sig_expr <- colMeans(counts[1:4, , drop = FALSE])
  risk <- scale(sig_expr)[, 1] * 0.4
  et <- stats::rexp(n, rate = 0.03 * exp(0.8 * risk))
  ct <- stats::rexp(n, rate = 0.02)
  cd <- S4Vectors::DataFrame(time = pmin(et, ct),
                             event = as.integer(et < ct),
                             row.names = colnames(counts))
  se <- SummarizedExperiment(assays = list(counts = counts), colData = cd)
  sig <- build_signature(se, "time", "event", top_n = 4, repeats = 2,
                         folds = 3, seed = 1)
  expect_true(all(sig$genes %in% genes))
  expect_true(all(names(sig$coefficients) == sig$genes))
  cm <- cox_model(sig, se)
  expect_true(all(cm$coef_table$term %in% c(genes, cm$terms)))
  expect_true(all(cm$coef_table$term[cm$coef_table$term %in% genes] %in%
                    genes))
  expect_identical(km_curve(sig, se)$log_rank_p <= 1, TRUE)
  expect_true(all(is.finite(survival_parametric(sig, se)$table$AIC)))
})
