library(SummarizedExperiment)

test_that("validate_external validates its inputs", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  expect_error(validate_external(list(genes = "a"), se, cutpoint = 0),
               "rnaSentry_signature")
  expect_error(validate_external(sig, matrix(1, 4, 4), cutpoint = 0),
               "SummarizedExperiment")
  expect_error(validate_external(sig, se), "never computed")
  expect_error(validate_external(sig, se, cutpoint = "high"),
               "single finite number")
  expect_error(validate_external(sig, se, cutpoint = 0,
                                 min_gene_overlap_frac = 0),
               "single number in")
  expect_error(validate_external(sig, se, cutpoint = 0,
                                 min_gene_overlap_frac = 1.5),
               "single number in")
  expect_error(validate_external(sig, se, cutpoint = 0, drop_missing = 1),
               "single TRUE or FALSE")
})

test_that("validate_external never recomputes cutpoint from external data", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cutoff <- km_curve(sig, se)$cutpoint
  ext <- make_survival_se(seed = 777)
  expect_error(validate_external(sig, ext), "never computed")
  val <- validate_external(sig, ext, cutpoint = cutoff)
  expect_equal(val$cutpoint, cutoff)
  expect_equal(val$cutpoint_type, "custom")
  expect_identical(unname(val$groups == "high"),
                   unname(val$score >= cutoff))
  ext_median <- stats::median(val$score)
  if (abs(ext_median - cutoff) > 1e-6) {
    expect_false(identical(unname(val$groups == "high"),
                           unname(val$score >= ext_median)))
  }
})

test_that("validate_external returns the documented structure", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cutoff <- km_curve(sig, se)$cutpoint
  ext <- make_cohort_se(sig, seed = 202)
  ext_cutoff <- stats::median(km_curve(sig, ext)$score)
  val <- validate_external(sig, ext, cutpoint = ext_cutoff)
  expect_s3_class(val, "rnaSentry_external")
  expect_setequal(val$genes_used, sig$genes)
  expect_length(val$genes_missing, 0)
  expect_equal(length(val$score), ncol(ext))
  expect_setequal(levels(val$groups), c("low", "high"))
  expect_true(all(table(val$groups) >= 1))
  expect_true(is.finite(val$log_rank_p) && val$log_rank_p >= 0 &&
                val$log_rank_p <= 1)
  expect_true(all(c("low", "high") %in% names(val$median_survival)))
  expect_false(val$sig_locked)
  expect_true(nzchar(val$created))
  expect_true(is.finite(val$concordance) || is.na(val$concordance))
  conc <- survival::concordance(
    survival::Surv(SummarizedExperiment::colData(ext)$time,
                   SummarizedExperiment::colData(ext)$event) ~ val$score
  )
  expect_equal(val$concordance, as.numeric(conc$concordance[1]))
})

test_that("validate_external errors strictly on missing genes by default", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cutoff <- km_curve(sig, se)$cutpoint
  ext <- make_cohort_se(sig, seed = 202)
  ext_cutoff <- stats::median(km_curve(sig, ext)$score)
  ext <- ext[rownames(ext) != sig$genes[1], ]
  expect_error(validate_external(sig, ext, cutpoint = ext_cutoff),
               sig$genes[1])
})

test_that("validate_external drops missing genes only when opted in", {
  se <- make_survival_se()
  sig <- build_signature(se, "time", "event", top_n = 10, repeats = 2,
                         folds = 3, seed = 7)
  cutoff <- km_curve(sig, se)$cutpoint
  ext <- make_cohort_se(sig, seed = 202)
  ext <- ext[rownames(ext) != sig$genes[1], ]
  present <- setdiff(sig$genes, sig$genes[1])
  mat19 <- log2(as.matrix(SummarizedExperiment::assay(ext, "counts")) + 1)
  score19 <- as.vector(sig$coefficients[present] %*% mat19[present, ])
  cutoff19 <- stats::median(score19)
  expect_error(validate_external(sig, ext, cutpoint = cutoff19,
                                 drop_missing = TRUE),
               "min_gene_overlap_frac")
  val <- validate_external(sig, ext, cutpoint = cutoff19,
                           drop_missing = TRUE, min_gene_overlap_frac = 0.9)
  expect_true(sig$genes[1] %in% val$genes_missing)
  expect_false(sig$genes[1] %in% val$genes_used)
  expect_true(any(val$flags$check == "missing_genes_dropped"))
  expect_match(val$flags$detail[val$flags$check == "missing_genes_dropped"],
               sig$genes[1])
})

test_that("validate_external records the lock state of the signature", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  cutoff <- stats::median(km_curve(sig, se)$score)
  expect_false(validate_external(sig, se, cutpoint = cutoff)$sig_locked)
  expect_true(validate_external(lock_signature(sig), se,
                                cutpoint = cutoff)$sig_locked)
})

test_that("print and plot are informative and return the object", {
  se <- make_survival_se()
  sig <- make_signature_for_testing()
  cutoff <- stats::median(km_curve(sig, se)$score)
  val <- validate_external(sig, se, cutpoint = cutoff)
  expect_output(print(val), "Log-rank")
  expect_output(print(val), "Concordance")
  grDevices::pdf(NULL)
  res <- plot(val)
  grDevices::dev.off()
  expect_identical(res, val)
})
