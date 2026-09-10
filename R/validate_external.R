# External validation of a signature against an independent cohort.

#' Validate a signature on an external cohort
#'
#' \code{validate_external()} applies a locked-down signature to an
#' independent cohort: it scores every external sample with the signature's
#' joint Cox coefficients (the gene-only model, so the score is portable to
#' cohorts that do not record the discovery cohort's design covariates), splits
#' the cohort into \code{"low"} and \code{"high"} risk groups at the supplied
#' \emph{discovery} cutpoint, and reports the external discriminant validity:
#' a log-rank test between the groups, per-group median survival and event
#' counts, and the concordance index of the continuous risk score against the
#' external survival outcome (Harrell's C in the risk-score convention: a
#' higher risk score is concordant with an earlier event).
#'
#' The stage \emph{never} recomputes a within-cohort cutpoint. The
#' \code{cutpoint} must be passed explicitly (for example the discovery-cohort
#' median recorded in the \code{cutpoint} element returned by
#' \code{\link{km_curve}}); passing no cutpoint is an error, so external risk
#' groups are always defined by the same threshold that defined the discovery
#' groups.
#'
#' Missing genes are handled strictly by default: \code{drop_missing = FALSE}
#' and \code{min_gene_overlap_frac = 1} (the defaults) require every signature
#' gene to be present, and error otherwise, naming the missing genes. Set
#' \code{drop_missing = TRUE} (with \code{min_gene_overlap_frac < 1} as
#' needed) to explicitly allow scoring on a subset of genes; the dropped genes
#' are named in a flag.
#'
#' The stage is read-only: it works on both locked and unlocked signatures and
#' records the signature's lock state in its output.
#'
#' @param sig An object of class \code{"rnaSentry_signature"} as returned by
#'   \code{\link{build_signature}}.
#' @param external_se A \code{SummarizedExperiment} of an independent cohort,
#'   with an expression assay and survival metadata.
#' @param time_col,event_col Names of the survival columns in
#'   \code{colData(external_se)}. When \code{NULL} (default) the column names
#'   recorded in \code{sig} are used.
#' @param cutpoint The risk-score threshold defining the \code{"high"} group
#'   (samples with score \code{>= cutpoint}). Must be a single finite number;
#'   it is never derived from the external cohort.
#' @param min_gene_overlap_frac Numeric in (0, 1]. The minimum fraction of
#'   signature genes that must be present in the external cohort when
#'   \code{drop_missing = TRUE}. Defaults to \code{1}.
#' @param drop_missing Logical. When \code{TRUE}, signature genes missing from
#'   the external cohort are dropped (and flagged by name), provided the
#'   remaining overlap is at least \code{min_gene_overlap_frac}. When
#'   \code{FALSE} (default) any missing gene is an error.
#'
#' @section Assumptions and limitations:
#' External risk groups are defined exclusively by the supplied discovery
#' cutpoint; the stage never re-estimates a threshold from the external
#' cohort. The cohort is assumed to be bulk RNA-seq with standard
#' right-censored survival and the same assay convention as the discovery
#' data (raw counts or a correctly-named \code{logcounts} assay; pre-scaled
#' data in a \code{counts}-named slot is flagged as
#' \code{possibly_log_scaled}). Competing risks and other non-standard
#' survival structures are not modeled. Unlike the discovery stage
#' (\code{\link{build_signature}}, which requires at least two events), a
#' one-event minimum applies here; the resulting instability is reported
#' through the \code{sparse_events} guardrail.
#'
#' @return An object of class \code{"rnaSentry_external"} (a list) with
#'   elements:
#'   \describe{
#'     \item{genes_used}{The signature genes present in the external cohort.}
#'     \item{genes_missing}{The signature genes absent from the external
#'       cohort (empty when none).}
#'     \item{score}{Per-sample external risk score.}
#'     \item{groups}{Factor with levels \code{"low"} and \code{"high"}
#'       from the discovery cutpoint.}
#'     \item{cutpoint, cutpoint_type}{The applied cutpoint and \code{"custom"}.}
#'     \item{km_fit}{The \code{survfit} object with one curve per group.}
#'     \item{log_rank_p}{p-value from the log-rank test between groups.}
#'     \item{median_survival}{Named vector of median survival per group.}
#'     \item{events_low, events_high}{Event counts per group.}
#'     \item{concordance}{Concordance index of the continuous risk score
#'       against the external survival outcome.}
#'     \item{sig_locked}{Whether the input signature was locked.}
#'     \item{created}{Provenance timestamp.}
#'     \item{flags}{Data.frame of issues raised (dropped genes, sparse
#'       events, non-finite concordance), with columns \code{check},
#'       \code{severity}, \code{detail} and \code{stage}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(11)
#' counts <- matrix(rpois(800, lambda = 500), nrow = 40, ncol = 20,
#'                  dimnames = list(paste0("gene", 1:40), paste0("S", 1:20)))
#' sig_expr <- colMeans(counts[1:5, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(20, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#' sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
#'                        folds = 2, seed = 1)
#' cutoff <- km_curve(sig, se)$cutpoint
#' ext <- se
#' ext <- ext[, sample(ncol(ext))]
#' val <- validate_external(sig, ext, cutpoint = cutoff)
#' val
#'
#' @export
validate_external <- function(sig, external_se, time_col = NULL,
                              event_col = NULL, cutpoint,
                              min_gene_overlap_frac = 1,
                              drop_missing = FALSE) {
  if (!inherits(sig, "rnaSentry_signature")) {
    stop("'sig' must be an rnaSentry_signature object.", call. = FALSE)
  }
  if (!methods::is(external_se, "SummarizedExperiment")) {
    stop("'external_se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  genes <- sig$genes
  coef <- sig$coefficients
  if (!is.character(genes) || length(genes) == 0 ||
      is.null(coef) || length(coef) == 0) {
    stop("'sig' has no genes or coefficients to score.", call. = FALSE)
  }
  if (missing(cutpoint) || is.null(cutpoint) || !is.numeric(cutpoint) ||
      length(cutpoint) != 1 || !is.finite(cutpoint)) {
    stop("'cutpoint' must be a single finite number, taken from the ",
         "discovery cohort (for example km_curve(sig, discovery)$cutpoint); ",
         "a cutpoint is never computed from the external cohort.",
         call. = FALSE)
  }
  if (!is.numeric(min_gene_overlap_frac) || length(min_gene_overlap_frac) != 1 ||
      is.na(min_gene_overlap_frac) || min_gene_overlap_frac <= 0 ||
      min_gene_overlap_frac > 1) {
    stop("'min_gene_overlap_frac' must be a single number in (0, 1].",
         call. = FALSE)
  }
  if (!is.logical(drop_missing) || length(drop_missing) != 1 ||
      is.na(drop_missing)) {
    stop("'drop_missing' must be a single TRUE or FALSE.", call. = FALSE)
  }

  cd <- SummarizedExperiment::colData(external_se)
  t_col <- if (is.null(time_col)) sig$time_col else time_col
  e_col <- if (is.null(event_col)) sig$event_col else event_col
  missing_ext <- setdiff(c(t_col, e_col), colnames(cd))
  if (length(missing_ext) > 0) {
    stop(sprintf("colData(external_se) has no column '%s'.", missing_ext[1]),
         call. = FALSE)
  }
  invalid_ext <- vapply(c(t_col, e_col), function(nm) {
    length(nm) != 1 || !is.character(nm) || is.na(nm)
  }, logical(1))
  if (any(invalid_ext)) {
    stop(sprintf("colData(external_se) has no column '%s'.",
                 c(t_col, e_col)[which(invalid_ext)[1]]), call. = FALSE)
  }
  time_vec <- as.numeric(cd[[t_col]])
  if (anyNA(time_vec) || any(time_vec < 0) || any(!is.finite(time_vec))) {
    stop(sprintf("'%s' must be a numeric, non-negative, finite, non-missing column.",
                 t_col), call. = FALSE)
  }
  event_raw <- cd[[e_col]]
  event_vec <- if (is.logical(event_raw)) {
    as.integer(event_raw)
  } else if (is.numeric(event_raw)) {
    as.integer(event_raw)
  } else {
    stop(sprintf("'%s' must be a 0/1 or logical event indicator.", e_col),
         call. = FALSE)
  }
  if (anyNA(event_vec) || !all(event_vec %in% c(0L, 1L))) {
    stop(sprintf("'%s' must contain only 0/1 values.", e_col), call. = FALSE)
  }
  if (sum(event_vec) < 1) {
    stop(sprintf("No events in '%s'; survival analysis requires at least one event.",
                 e_col), call. = FALSE)
  }

  flags <- .new_flags("validate_external")

  mat_info <- .get_analysis_matrix(external_se)
  mat <- mat_info$mat
  if (isTRUE(mat_info$log_scaled_possible)) {
    warning(mat_info$log_scaled_msg, call. = FALSE)
    flags <- .add_flag(flags, "possibly_log_scaled", "warning",
                       mat_info$log_scaled_msg)
  }
  present <- genes[genes %in% rownames(mat)]
  missing_genes <- setdiff(genes, rownames(mat))

  if (length(missing_genes) > 0) {
    if (!drop_missing) {
      stop(sprintf("The external cohort is missing signature gene(s): %s. To validate on the present subset, set drop_missing = TRUE (with an appropriate min_gene_overlap_frac).",
                   paste(missing_genes, collapse = ", ")), call. = FALSE)
    }
    overlap_frac <- length(present) / length(genes)
    if (overlap_frac < min_gene_overlap_frac) {
      stop(sprintf("Only %.3f of signature genes are present in the external cohort, below min_gene_overlap_frac = %.3f.",
                   overlap_frac, min_gene_overlap_frac), call. = FALSE)
    }
    flags <- .add_flag(flags, "missing_genes_dropped", "warning",
                       sprintf("Signature gene(s) missing from the external cohort and dropped: %s.",
                               paste(missing_genes, collapse = ", ")))
  }

  b <- coef[present]
  if (anyNA(b) || !all(is.finite(b))) {
    stop("The signature coefficients are not all finite.", call. = FALSE)
  }
  score <- as.vector(b %*% mat[present, , drop = FALSE])
  names(score) <- colnames(mat)
  if (isTRUE(stats::sd(score) == 0)) {
    stop("The external risk score is constant across samples; cannot define risk groups.",
         call. = FALSE)
  }
  if (sum(score >= cutpoint) == 0 || sum(score < cutpoint) == 0) {
    stop("'cutpoint' leaves one risk group empty in the external cohort.",
         call. = FALSE)
  }
  groups <- factor(ifelse(score >= cutpoint, "high", "low"),
                   levels = c("low", "high"))

  d <- data.frame(time = time_vec, event = event_vec, group = groups,
                  score = score)
  fit <- survival::survfit(survival::Surv(time, event) ~ group, data = d)
  dd <- survival::survdiff(survival::Surv(time, event) ~ group, data = d)
  log_rank_p <- stats::pchisq(dd$chisq, df = max(1, length(dd$n) - 1),
                              lower.tail = FALSE)

  sm <- summary(fit)$table
  strata_names <- sub("^.*=", "", rownames(sm))
  med_surv <- as.numeric(sm[, "median"])
  names(med_surv) <- strata_names
  median_survival <- c(
    low = unname(med_surv[["low"]]),
    high = unname(med_surv[["high"]])
  )

  events_low <- sum(d$event[d$group == "low"])
  events_high <- sum(d$event[d$group == "high"])
  if (events_low < 3 || events_high < 3) {
    flags <- .add_flag(flags, "sparse_events", "warning",
                       sprintf("Risk group event counts are low (low: %d, high: %d); median survival estimates may be unstable.",
                               events_low, events_high))
  }

  # Risk score: larger = higher hazard = shorter survival, so the concordance
  # must use reverse = TRUE (survival::concordance's default means "larger x =>
  # longer survival"). Regression guard: test-statistical_parity.R
  # ("survival::concordance reverse convention satisfies C + C_rev = 1").
  conc <- tryCatch(
    withCallingHandlers(
      survival::concordance(survival::Surv(time, event) ~ score,
                            data = d, reverse = TRUE),
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) NULL
  )
  concordance <- if (is.null(conc)) NA_real_ else
    as.numeric(conc$concordance[1])
  if (!is.finite(concordance)) {
    flags <- .add_flag(flags, "concordance_na", "warning",
                       "External concordance could not be estimated (survival::concordance returned a non-finite value); interpret the external validation cautiously.")
  }

  result <- list(
    genes_used = present,
    genes_missing = missing_genes,
    score = score,
    groups = groups,
    cutpoint = cutpoint,
    cutpoint_type = "custom",
    km_fit = fit,
    log_rank_p = log_rank_p,
    median_survival = median_survival,
    events_low = events_low,
    events_high = events_high,
    concordance = concordance,
    sig_locked = isTRUE(sig$locked),
    created = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     " UTC"),
    flags = flags
  )
  class(result) <- "rnaSentry_external"
  result
}

#' @export
print.rnaSentry_external <- function(x, ...) {
  cat(sprintf("rnaSentry external validation: %d/%d signature gene(s) used, %d samples.\n",
              length(x$genes_used),
              length(x$genes_used) + length(x$genes_missing),
              length(x$score)))
  if (length(x$genes_missing) > 0) {
    cat(sprintf("Missing genes: %s.\n",
                paste(x$genes_missing, collapse = ", ")))
  }
  cat(sprintf("Risk groups: low n = %d (%d events), high n = %d (%d events).\n",
              sum(x$groups == "low"), x$events_low,
              sum(x$groups == "high"), x$events_high))
  cat(sprintf("Cutpoint %.3g (custom, from discovery).\n", x$cutpoint))
  cat(sprintf("Log-rank p = %s. Concordance = %s.\n",
              format.pval(x$log_rank_p),
              if (is.na(x$concordance)) "NA" else
                sprintf("%.3f", x$concordance)))
  cat(if (x$sig_locked) "Signature locked.\n" else "Signature not locked.\n")
  if (nrow(x$flags) > 0) {
    cat(sprintf("%d issue(s) flagged:\n", nrow(x$flags)))
    cat(sprintf("  [%s] %s: %s\n", x$flags$severity, x$flags$check, x$flags$detail),
        sep = "")
  }
  invisible(x)
}

#' @export
plot.rnaSentry_external <- function(x, ...) {
  cols <- c(low = "steelblue", high = "firebrick")
  graphics::plot(x$km_fit, col = cols[c("low", "high")], lty = 1, lwd = 2,
       xlab = "Time", ylab = "Survival probability", ...)
  graphics::legend("topright", legend = c("low risk", "high risk"),
         col = cols[c("low", "high")], lty = 1, lwd = 2, bty = "n")
  graphics::mtext(sprintf("log-rank p = %s",
                format.pval(x$log_rank_p)),
        side = 3, line = 0.5, cex = 0.9)
  invisible(x)
}
