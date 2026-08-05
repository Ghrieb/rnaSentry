# Run the complete rnaSentry pipeline on a single cohort.

#' Run the complete rnaSentry pipeline
#'
#' \code{run_rnaSentry()} executes the full rnaSentry pipeline on a single
#' cohort and collects every stage result into one object:
#' \enumerate{
#'   \item \code{\link{design_audit}} over \code{design_vars} (when any are
#'     supplied).
#'   \item \code{\link{pca_audit}} (with \code{batch_col} when supplied).
#'   \item \code{\link{build_signature}}, then \code{\link{lock_signature}} so
#'     the survival-modeling stages operate on the published signature.
#'   \item \code{\link{km_curve}}, \code{\link{cox_model}} and
#'     \code{\link{survival_parametric}} on the locked signature.
#'   \item \code{\link{generate_report}} assembling all stages into an HTML
#'     report (when \code{render_report = TRUE}).
#' }
#'
#' The result's \code{stages} list can be passed directly to
#' \code{\link{generate_report}}, or any stage can be re-run individually
#' (for example \code{\link{validate_external}} against an independent
#' cohort, using the cutpoint recorded by the \code{km_curve} stage).
#'
#' @section Reproducibility lock:
#' The run locks the discovered signature (see \code{\link{lock_signature}}),
#' recording a session-level fingerprint. While any locked signature exists,
#' \code{run_rnaSentry()} refuses to run again, preventing silent re-selection
#' of signature genes after survival analysis. To run the pipeline a second
#' time in the same session (for example on another cohort), release the lock
#' first with \code{lock_signature(sig, lock = FALSE)} on the locked
#' signature. Individual stages remain re-runnable while the lock is set.
#'
#' @param se A \code{SummarizedExperiment} with a count or normalized
#'   expression assay and survival metadata.
#' @param time_col Character. Column of \code{colData(se)} with follow-up
#'   time.
#' @param event_col Character. Column of \code{colData(se)} with the event
#'   indicator (0/1 or logical).
#' @param outcome_col Character. Display label for the outcome. Defaults to
#'   \code{"overall_survival"}.
#' @param design_terms Character vector of \code{colData(se)} columns to
#'   adjust the signature's univariate screening by. Defaults to
#'   \code{character(0)}.
#' @param design_vars Character vector of \code{colData(se)} columns scanned
#'   by the \code{\link{design_audit}} confounder scan. Defaults to
#'   \code{design_terms}.
#' @param batch_col Optional character. \code{colData(se)} column tested by
#'   the \code{\link{pca_audit}} batch scan. Defaults to \code{NULL}.
#' @param method,top_n,p_threshold,repeats,folds,seed,adjust_for_design,min_events_per_parameter,BPPARAM Passed to \code{\link{build_signature}}.
#' @param report_file Character. Output file name for the HTML report.
#' @param report_dir Character. Directory (which must exist) to write the
#'   report into.
#' @param render_report Logical. When \code{TRUE} (default) the report is
#'   rendered as the final step.
#'
#' @section Assumptions and limitations:
#' The pipeline targets bulk RNA-seq with standard right-censored survival
#' and linear Cox risk scores (see the Assumptions sections of
#' \code{\link{build_signature}} and \code{\link{cox_model}}). It detects and
#' reports batch and design confounders but does not correct for them, and it
#' makes no differential-expression calls. On small cohorts an
#' \code{events_per_parameter} warning from \code{\link{build_signature}} is
#' expected and should be interpreted as a call for caution, not ignored.
#'
#' @return An object of class \code{"rnaSentry_run"} (a list) with elements:
#'   \describe{
#'     \item{stages}{Named list of the stage results, with entries
#'       \code{design_audit} (only when \code{design_vars} was supplied),
#'       \code{pca_audit}, \code{build_signature}, \code{km_curve},
#'       \code{cox_model} and \code{survival_parametric}.}
#'     \item{report}{Path of the rendered HTML report, or \code{NULL} when
#'       \code{render_report = FALSE}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(8)
#' counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
#'                   dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
#' sig_expr <- colMeans(counts[1:5, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(20, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  batch = factor(rep(c("B1", "B2"), 10)),
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' run <- run_rnaSentry(se, "time", "event", design_vars = "batch",
#'                       top_n = 5, repeats = 1, folds = 2, seed = 1,
#'                       report_dir = tempdir())
#' run
#'
#' @export
run_rnaSentry <- function(se, time_col, event_col,
                          outcome_col = "overall_survival",
                          design_terms = character(0),
                          design_vars = design_terms,
                          batch_col = NULL,
                          method = c("top_n", "p_value"), top_n = 20,
                          p_threshold = 0.05, repeats = 5, folds = 5,
                          seed = NULL, adjust_for_design = TRUE,
                          min_events_per_parameter = 5,
                          BPPARAM = NULL,
                          report_file = "rnaSentry_report.html",
                          report_dir = ".", render_report = TRUE) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  if (!is.character(report_dir) || length(report_dir) != 1 ||
      !dir.exists(report_dir)) {
    stop("'report_dir' must name an existing directory.", call. = FALSE)
  }
  if (!is.logical(render_report) || length(render_report) != 1 ||
      is.na(render_report)) {
    stop("'render_report' must be a single TRUE or FALSE.", call. = FALSE)
  }

  stages <- list()

  # Lock enforcement: if a locked signature already exists in this session,
  # refuse to re-run the pipeline.  The user must unlock first or start fresh.
  locked_env <- get(".rnaSentry_locked_sigs", envir = asNamespace("rnaSentry"))
  if (length(ls(locked_env)) > 0) {
    stop("A locked signature already exists in this session. ",
         "Call lock_signature(sig, lock = FALSE) on the locked signature ",
         "before re-running the pipeline, or start a fresh R session.",
         call. = FALSE)
  }

  if (length(design_vars) > 0) {
    stages$design_audit <- design_audit(se, design_vars = design_vars)
  }

  stages$pca_audit <- pca_audit(se, batch_col = batch_col)

  sig <- build_signature(se, time_col, event_col, outcome_col = outcome_col,
                         design_terms = design_terms, method = method,
                         top_n = top_n, p_threshold = p_threshold,
                         repeats = repeats, folds = folds, seed = seed,
                         adjust_for_design = adjust_for_design,
                         min_events_per_parameter = min_events_per_parameter,
                         BPPARAM = BPPARAM)
  sig_locked <- lock_signature(sig)
  stages$build_signature <- sig_locked

  stages$km_curve <- km_curve(sig_locked, se)
  stages$cox_model <- cox_model(sig_locked, se)
  stages$survival_parametric <- survival_parametric(sig_locked, se)

  report <- NULL
  if (render_report) {
    report <- generate_report(stages, output_file = report_file,
                              output_dir = report_dir)
  }

  result <- list(stages = stages, report = report)
  class(result) <- "rnaSentry_run"
  result
}

#' @export
print.rnaSentry_run <- function(x, ...) {
  cat("rnaSentry pipeline run.\n")
  for (nm in names(x$stages)) {
    cat(sprintf("  - %s\n", nm))
  }
  if (!is.null(x$report)) {
    cat(sprintf("Report: %s\n", x$report))
  }
  invisible(x)
}
