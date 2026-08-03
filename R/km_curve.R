#' Kaplan-Meier survival analysis of a signature's risk groups
#'
#' Scores each sample of a \code{SummarizedExperiment} with the signature's
#' joint Cox coefficients, splits the cohort into \code{"low"} and
#' \code{"high"} risk groups at the median risk score, and fits a Kaplan-Meier
#' survival curve for each group with a log-rank test. The analysis is
#' read-only: it works on both locked and unlocked signatures and records the
#' signature's lock state in its output.
#'
#' @param sig An object of class \code{"rnaSentry_signature"} as returned by
#'   \code{\link{build_signature}}.
#' @param se A \code{SummarizedExperiment} containing the signature genes in
#'   its rows and the survival metadata columns recorded in \code{sig}
#'   (\code{sig$time_col} and \code{sig$event_col}). Expression is taken from
#'   the same analysis assay used at signature-build time (see
#'   \code{\link{build_signature}}).
#'
#' @return An object of class \code{"rnaSentry_km"} (a list) with elements:
#'   \describe{
#'     \item{genes}{The signature genes scored.}
#'     \item{score}{Per-sample risk score (linear predictor).}
#'     \item{groups}{Factor with levels \code{"low"} and \code{"high"}
#'       from the median split.}
#'     \item{fit}{The \code{survfit} object with one curve per group.}
#'     \item{log_rank_p}{p-value from the log-rank test.}
#'     \item{median_survival}{Named vector of median survival per group.}
#'     \item{events_low, events_high}{Event counts per group.}
#'     \item{sig_locked}{Whether the input signature was locked.}
#'     \item{created}{Provenance timestamp.}
#'     \item{flags}{Data.frame of issues raised (e.g. sparse events per
#'       group).}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(5)
#' counts <- matrix(rpois(120, lambda = 500), nrow = 12, ncol = 10,
#'                  dimnames = list(paste0("gene", 1:12), paste0("S", 1:10)))
#' sig_expr <- colMeans(counts[1:3, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(10, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(10, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#' sig <- build_signature(se, time_col = "time", event_col = "event",
#'                        top_n = 3, repeats = 1, folds = 2, seed = 1)
#' km <- km_curve(sig, se)
#' km
#'
#' @export
km_curve <- function(sig, se) {
  if (!inherits(sig, "rnaSentry_signature")) {
    stop("'sig' must be an rnaSentry_signature object.", call. = FALSE)
  }
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  genes <- sig$genes
  coef <- sig$coefficients
  if (!is.character(genes) || length(genes) == 0 ||
      is.null(coef) || length(coef) == 0) {
    stop("'sig' has no genes or coefficients to score.", call. = FALSE)
  }

  cd <- SummarizedExperiment::colData(se)
  for (nm in c(sig$time_col, sig$event_col)) {
    if (length(nm) != 1 || !is.character(nm) || !nm %in% colnames(cd)) {
      stop(sprintf("colData(se) has no column '%s'.", nm), call. = FALSE)
    }
  }
  time_vec <- as.numeric(cd[[sig$time_col]])
  if (anyNA(time_vec) || any(time_vec < 0)) {
    stop(sprintf("'%s' must be a numeric, non-negative, non-missing column.",
                 sig$time_col), call. = FALSE)
  }
  event_raw <- cd[[sig$event_col]]
  event_vec <- if (is.logical(event_raw)) {
    as.integer(event_raw)
  } else if (is.numeric(event_raw)) {
    as.integer(event_raw)
  } else {
    stop(sprintf("'%s' must be a 0/1 or logical event indicator.",
                 sig$event_col), call. = FALSE)
  }
  if (anyNA(event_vec) || !all(event_vec %in% c(0L, 1L))) {
    stop(sprintf("'%s' must contain only 0/1 values.", sig$event_col),
         call. = FALSE)
  }
  if (sum(event_vec) < 1) {
    stop(sprintf("No events in '%s'; survival analysis requires at least one event.",
                 sig$event_col), call. = FALSE)
  }

  mat <- .get_analysis_matrix(se)$mat
  missing_genes <- setdiff(genes, rownames(mat))
  if (length(missing_genes) > 0) {
    stop(sprintf("The following signature gene(s) are missing from the SummarizedExperiment: %s.",
                 paste(missing_genes, collapse = ", ")), call. = FALSE)
  }
  b <- coef[genes]
  if (anyNA(b) || !all(is.finite(b))) {
    stop("The signature coefficients are not all finite.", call. = FALSE)
  }

  score <- as.vector(b %*% mat[genes, , drop = FALSE])
  names(score) <- colnames(mat)
  if (stats::sd(score) == 0) {
    stop("The risk score is constant across samples; cannot define risk groups.",
         call. = FALSE)
  }
  med <- stats::median(score)
  groups <- factor(ifelse(score >= med, "high", "low"),
                   levels = c("low", "high"))

  flags <- data.frame(check = character(0), severity = character(0),
                       detail = character(0), stringsAsFactors = FALSE)
  add_flag <- function(flags, check, severity, detail) {
    rbind(flags, data.frame(check = check, severity = severity,
                             detail = detail, stringsAsFactors = FALSE))
  }

  d <- data.frame(time = time_vec, event = event_vec, group = groups)
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
    flags <- add_flag(flags, "sparse_events", "warning",
                      sprintf("Risk group event counts are low (low: %d, high: %d); median survival estimates may be unstable.",
                              events_low, events_high))
  }

  result <- list(
    genes = genes,
    score = score,
    groups = groups,
    fit = fit,
    log_rank_p = log_rank_p,
    median_survival = median_survival,
    events_low = events_low,
    events_high = events_high,
    sig_locked = isTRUE(sig$locked),
    created = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     " UTC"),
    flags = flags
  )
  class(result) <- "rnaSentry_km"
  result
}

#' @export
print.rnaSentry_km <- function(x, ...) {
  cat(sprintf("rnaSentry KM analysis: %d-gene signature, %d samples.\n",
              length(x$genes), length(x$score)))
  cat(sprintf("Risk groups: low n = %d (%d events), high n = %d (%d events).\n",
              sum(x$groups == "low"), x$events_low,
              sum(x$groups == "high"), x$events_high))
  cat(sprintf("Median survival: low %.3g, high %s.\n",
              x$median_survival[["low"]],
              if (is.na(x$median_survival[["high"]])) "NA" else
                sprintf("%.3g", x$median_survival[["high"]])))
  cat(sprintf("Log-rank p = %s.\n", format.pval(x$log_rank_p)))
  cat(if (x$sig_locked) "Signature locked.\n" else "Signature not locked.\n")
  if (nrow(x$flags) > 0) {
    cat(sprintf("%d issue(s) flagged:\n", nrow(x$flags)))
    for (i in seq_len(nrow(x$flags))) {
      cat(sprintf("  [%s] %s: %s\n", x$flags$severity[i],
                  x$flags$check[i], x$flags$detail[i]))
    }
  }
  invisible(x)
}

#' @export
plot.rnaSentry_km <- function(x, ...) {
  cols <- c(low = "steelblue", high = "firebrick")
  plot(x$fit, col = cols[c("low", "high")], lty = 1, lwd = 2,
       xlab = "Time", ylab = "Survival probability", ...)
  graphics::legend("topright", legend = c("low risk", "high risk"),
         col = cols[c("low", "high")], lty = 1, lwd = 2, bty = "n")
  graphics::mtext(sprintf("log-rank p = %s",
                format.pval(x$log_rank_p)),
        side = 3, line = 0.5, cex = 0.9)
  invisible(x)
}
