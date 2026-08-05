# Render the rnaSentry analysis report.

#' Generate the rnaSentry analysis report
#'
#' \code{generate_report()} assembles the results of any collection of
#' rnaSentry stages into a single self-contained HTML report. It renders the
#' plain \code{report_template.Rmd} shipped with the package (via
#' \code{\link[rmarkdown]{render}}), which prints a pipeline summary, a unified
#' audit-trail table collecting every flag from every supplied stage (all
#' rnaSentry stage results carry a \code{flags} data.frame with
#' \code{check}/\code{severity}/\code{detail}/\code{stage} columns), and one
#' section per stage using its \code{print} method.
#'
#' @param stages A named list of rnaSentry stage results (for example
#'   elements of the \code{stages} list returned by \code{\link{run_rnaSentry}}
#'   or a hand-assembled set of audit, signature, survival-modeling and
#'   validation outputs). Names are used as report section titles.
#' @param output_file Character. File name (without directory) of the rendered
#'   report. Defaults to \code{"rnaSentry_report.html"}.
#' @param output_dir Character. Directory to write the report into. Must
#'   exist. Defaults to the working directory.
#' @param quiet Logical. Passed to \code{rmarkdown::render()}; suppresses the
#'   rendering console output when \code{TRUE} (default).
#'
#' @return The path of the rendered report, invisibly.
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(9)
#' counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
#'                   dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
#' sig_expr <- colMeans(counts[1:5, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(20, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
#'                         folds = 2, seed = 1)
#' stages <- list(
#'   km_curve = km_curve(sig, se),
#'   cox_model = cox_model(sig, se)
#' )
#' generate_report(stages, output_dir = tempdir())
#'
#' @export
generate_report <- function(stages, output_file = "rnaSentry_report.html",
                            output_dir = ".", quiet = TRUE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("generate_report() requires the 'rmarkdown' package (in Suggests).",
         call. = FALSE)
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("generate_report() requires the 'knitr' package (in Suggests).",
         call. = FALSE)
  }
  if (!is.list(stages) || length(stages) == 0 || is.null(names(stages)) ||
      any(names(stages) == "")) {
    stop("'stages' must be a non-empty named list of stage results.",
         call. = FALSE)
  }
  if (!is.character(output_file) || length(output_file) != 1 ||
      is.na(output_file) || !nzchar(output_file)) {
    stop("'output_file' must be a single non-empty character string.",
         call. = FALSE)
  }
  if (!is.character(output_dir) || length(output_dir) != 1 ||
      !dir.exists(output_dir)) {
    stop("'output_dir' must name an existing directory.", call. = FALSE)
  }

  template <- system.file("report_template.Rmd", package = "rnaSentry")
  if (!nzchar(template)) {
    stop("The report template is not installed; reinstall rnaSentry.",
         call. = FALSE)
  }

  rendered <- rmarkdown::render(
    template,
    output_file = output_file,
    output_dir = output_dir,
    params = list(stages = stages),
    quiet = quiet
  )
  invisible(rendered)
}
