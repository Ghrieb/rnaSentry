library(SummarizedExperiment)

# Shared synthetic fixtures for the rnaSentry test suite. Constructors build
# small SummarizedExperiment objects inline, mirroring the style of the
# existing test files, so pipeline tests stay fast and network-free.

#' Build a synthetic survival SummarizedExperiment
#'
#' @param n_genes Number of genes (rows).
#' @param n_samples Number of samples (columns).
#' @param n_signal_genes Number of leading genes driving the survival signal.
#' @param seed RNG seed for reproducibility.
#'
#' @return A SummarizedExperiment with a raw-count assay and colData columns
#'   \code{time}, \code{event} (0/1), \code{condition} (factor A/B),
#'   \code{batch} (factor B1/B2/B3) and \code{age} (numeric).
make_survival_se <- function(n_genes = 120, n_samples = 60,
                             n_signal_genes = 5, seed = 101) {
  set.seed(seed)
  counts <- matrix(stats::rpois(n_genes * n_samples, lambda = 500),
                   nrow = n_genes, ncol = n_samples,
                   dimnames = list(paste0("gene", seq_len(n_genes)),
                                   paste0("S", seq_len(n_samples))))
  # engineer a survival signal through the leading genes so that Cox/CV
  # stages have something to recover and every CV fold sees events
  sig_expr <- colMeans(counts[seq_len(n_signal_genes), , drop = FALSE])
  risk <- scale(sig_expr)[, 1] * 0.4
  event_time <- stats::rexp(n_samples, rate = 0.03 * exp(0.8 * risk))
  censor_time <- stats::rexp(n_samples, rate = 0.02)
  time <- pmin(event_time, censor_time)
  event <- as.integer(event_time < censor_time)

  coldata <- S4Vectors::DataFrame(
    time = time,
    event = event,
    condition = factor(rep(c("A", "B"), length.out = n_samples)),
    batch = factor(rep(c("B1", "B2", "B3"), length.out = n_samples)),
    age = stats::runif(n_samples, 40, 80),
    row.names = colnames(counts)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = coldata
  )
}

#' Small synthetic survival SummarizedExperiment for guardrail tests
#'
#' @return A SummarizedExperiment with 30 genes and 30 samples.
make_small_survival_se <- function() {
  make_survival_se(n_genes = 30, n_samples = 30, n_signal_genes = 3, seed = 303)
}

#' Build a signature object in a controlled state
#'
#' Constructs an \code{rnaSentry_signature} list literal (same shape
#' \code{build_signature()} returns) so that later-stage tests can exercise
#' lock/guardrail logic without running the full build.
#'
#' @param n_genes Number of genes in the signature.
#' @param locked Logical; whether the signature is locked.
#' @param time_col,event_col colData column names recorded in the signature.
#'
#' @return An object of class \code{"rnaSentry_signature"}.
make_signature_for_testing <- function(n_genes = 20, locked = FALSE,
                                       time_col = "time", event_col = "event") {
  genes <- paste0("gene", seq_len(n_genes))
  sig <- list(
    genes = genes,
    outcome = "overall_survival",
    time_col = time_col,
    event_col = event_col,
    design_terms = character(0),
    cox_stats = data.frame(
      gene = genes,
      HR = stats::runif(n_genes, 0.8, 1.2),
      p = stats::runif(n_genes),
      adj_p = stats::runif(n_genes),
      stringsAsFactors = FALSE
    ),
    coefficients = stats::setNames(stats::runif(n_genes, -0.3, 0.3), genes),
    cv_results = data.frame(repeat_id = 1, fold = seq_len(5),
                            c_index = stats::runif(5, 0.55, 0.8)),
    cv_summary = c(mean = 0.65, sd = 0.05),
    selection = list(method = "top_n", n_genes = n_genes),
    locked = locked,
    created = "2026-01-01 00:00:00 UTC",
    lock_time = if (locked) "2026-01-01 00:00:00 UTC" else NULL
  )
  class(sig) <- "rnaSentry_signature"
  sig
}

#' Build a synthetic external-validation cohort
#'
#' The cohort shares all of \code{sig$genes} with the signature (so the
#' happy path re-scores fully) plus unrelated genes.
#'
#' @param sig An \code{rnaSentry_signature} object.
#' @param n_samples Number of cohort samples.
#' @param seed RNG seed.
#'
#' @return A SummarizedExperiment with rownames including \code{sig$genes}.
make_cohort_se <- function(sig, n_samples = 60, seed = 202) {
  set.seed(seed)
  cohort_genes <- c(sig$genes, paste0("cohort_gene", seq_len(40)))
  counts <- matrix(stats::rpois(length(cohort_genes) * n_samples, lambda = 400),
                   nrow = length(cohort_genes), ncol = n_samples,
                   dimnames = list(cohort_genes, paste0("C", seq_len(n_samples))))
  event_time <- stats::rexp(n_samples, rate = 0.03)
  censor_time <- stats::rexp(n_samples, rate = 0.02)
  time <- pmin(event_time, censor_time)
  event <- as.integer(event_time < censor_time)
  coldata <- S4Vectors::DataFrame(
    time = time,
    event = event,
    condition = factor(rep(c("A", "B"), length.out = n_samples)),
    row.names = colnames(counts)
  )
  SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts), colData = coldata
  )
}
