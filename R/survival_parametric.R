# Parametric survival modeling of a signature's risk score.

#' Parametric survival models for a signature's risk score
#'
#' \code{survival_parametric()} scores each sample of a
#' \code{SummarizedExperiment} with the signature's joint Cox coefficients and
#' fits a set of parametric accelerated-failure-time models of the event time
#' on that score, using \code{\link[survival]{survreg}}. It returns an AIC
#' comparison across the fitted distributions, the fitted survival curves
#' (evaluated at the mean risk score) for overlay on the Kaplan-Meier
#' estimate, and a \code{plot()} method that draws the KM curve together with
#' the parametric curves.
#'
#' The stage is read-only: it works on both locked and unlocked signatures and
#' records the signature's lock state in its output.
#'
#' @param sig An object of class \code{"rnaSentry_signature"} as returned by
#'   \code{\link{build_signature}}.
#' @param se A \code{SummarizedExperiment} containing the signature genes in
#'   its rows and the survival metadata columns recorded in \code{sig}
#'   (\code{sig$time_col} and \code{sig$event_col}). Expression is taken from
#'   the same analysis assay used at signature-build time (see
#'   \code{\link{build_signature}}).
#' @param dists Character vector of parametric distributions to fit, any
#'   subset of \code{"weibull"}, \code{"exponential"}, \code{"lognormal"} and
#'   \code{"loglogistic"}. Defaults to all four.
#'
#' @return An object of class \code{"rnaSentry_parametric"} (a list) with
#'   elements:
#'   \describe{
#'     \item{genes}{The signature genes scored.}
#'     \item{score}{Per-sample risk score (linear predictor).}
#'     \item{fits}{Named list of \code{survreg} fits, one per distribution.}
#'     \item{table}{Data.frame comparing distributions on \code{loglik},
#'       \code{npar}, \code{AIC}, \code{delta_AIC}, \code{weight} and
#'       \code{best}, ordered by increasing AIC.}
#'     \item{best}{Name of the lowest-AIC distribution.}
#'     \item{km_fit}{The \code{survfit} object of the whole cohort.}
#'     \item{curves}{Named list of fitted survival curves (data.frame with
#'       \code{time} and \code{survival}), evaluated at the mean risk score.}
#'     \item{sig_locked}{Whether the input signature was locked.}
#'     \item{created}{Provenance timestamp.}
#'     \item{flags}{Data.frame of issues raised (e.g. a model that failed to
#'       fit, or near-tied AIC values), with columns \code{check},
#'       \code{severity}, \code{detail} and \code{stage}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(9)
#' counts <- matrix(rpois(450, lambda = 500), nrow = 30, ncol = 15,
#'                  dimnames = list(paste0("gene", 1:30), paste0("S", 1:15)))
#' sig_expr <- colMeans(counts[1:4, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(15, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(15, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#' sig <- build_signature(se, "time", "event", top_n = 4, repeats = 1,
#'                        folds = 2, seed = 1)
#' sp <- survival_parametric(sig, se)
#' sp
#' plot(sp)
#'
#' @export
survival_parametric <- function(sig, se,
                                dists = c("weibull", "exponential",
                                          "lognormal", "loglogistic")) {
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
  dists <- match.arg(dists, several.ok = TRUE,
                     choices = c("weibull", "exponential",
                                 "lognormal", "loglogistic"))
  cd <- SummarizedExperiment::colData(se)
  for (nm in c(sig$time_col, sig$event_col)) {
    if (length(nm) != 1 || !is.character(nm) || !nm %in% colnames(cd)) {
      stop(sprintf("colData(se) has no column '%s'.", nm), call. = FALSE)
    }
  }
  time_vec <- as.numeric(cd[[sig$time_col]])
  if (anyNA(time_vec) || any(time_vec < 0) || any(!is.finite(time_vec))) {
    stop(sprintf("'%s' must be a numeric, non-negative, finite, non-missing column.",
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
    stop(sprintf("No events in '%s'; survival modeling requires at least one event.",
                 sig$event_col), call. = FALSE)
  }

  flags <- .new_flags("survival_parametric")

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
  if (isTRUE(stats::sd(score) == 0)) {
    stop("The risk score is constant across samples; cannot fit a parametric model.",
         call. = FALSE)
  }

  d <- data.frame(time = time_vec, event = event_vec, score = score)
  res <- lapply(dists, function(dist) {
    fit <- tryCatch(
      withCallingHandlers(
        survival::survreg(survival::Surv(time, event) ~ score, data = d,
                          dist = dist),
        warning = function(w) invokeRestart("muffleWarning")
      ),
      error = function(e) NULL
    )
    list(fit = fit, dist = dist)
  })
  fits <- setNames(lapply(res, `[[`, "fit"), dists)
  keep <- !vapply(fits, is.null, logical(1))
  failed <- dists[!keep]
  fits <- fits[keep]
  for (dist in failed) {
    flags <- .add_flag(flags, "model_fit_failed", "warning",
                       sprintf("The %s model failed to fit and was excluded from the comparison.",
                               dist))
  }
  if (length(fits) == 0) {
    stop("None of the requested parametric models could be fitted.",
         call. = FALSE)
  }

  loglik_vals <- vapply(fits, function(f) as.numeric(stats::logLik(f)),
                        numeric(1))
  npar_vals <- vapply(fits, function(f) length(f$coefficients) + 1,
                      numeric(1))
  aic_vals <- vapply(fits, stats::AIC, numeric(1))
  aic_min <- min(aic_vals)
  delta <- aic_vals - aic_min
  weight <- exp(-0.5 * delta) / sum(exp(-0.5 * delta))
  best <- names(aic_vals)[which.min(aic_vals)]
  table <- data.frame(
    dist = names(fits),
    loglik = loglik_vals,
    npar = npar_vals,
    AIC = aic_vals,
    delta_AIC = unname(delta),
    weight = unname(weight),
    best = names(fits) == best,
    stringsAsFactors = FALSE
  )
  table <- table[order(table$AIC), ]
  if (nrow(table) >= 2 && diff(sort(aic_vals)[seq_len(2)]) < 2) {
    flags <- .add_flag(flags, "aic_models_close", "info",
                       "The two best-fitting models differ in AIC by less than 2; their evidence is nearly equivalent.")
  }

  km_fit <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
  ps <- seq(0.02, 0.98, length.out = 100)
  svals <- 1 - ps
  curves <- lapply(fits, function(f) {
    qt <- stats::predict(f, newdata = data.frame(score = mean(score)),
                         type = "quantile", p = ps)
    data.frame(time = as.numeric(qt), survival = svals)
  })

  result <- list(
    genes = genes,
    score = score,
    fits = fits,
    table = table,
    best = best,
    km_fit = km_fit,
    curves = curves,
    sig_locked = isTRUE(sig$locked),
    created = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     " UTC"),
    flags = flags
  )
  class(result) <- "rnaSentry_parametric"
  result
}

#' @export
print.rnaSentry_parametric <- function(x, ...) {
  cat(sprintf("rnaSentry parametric models: %d-gene signature, %d samples.\n",
              length(x$genes), length(x$score)))
  cat("AIC comparison:\n")
  for (i in seq_len(nrow(x$table))) {
    star <- if (x$table$best[i]) " *" else ""
    cat(sprintf("  %-12s AIC %8.2f  delta %6.2f  weight %5.3f%s\n",
                x$table$dist[i], x$table$AIC[i], x$table$delta_AIC[i],
                x$table$weight[i], star))
  }
  cat(sprintf("Best model: %s.\n", x$best))
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
plot.rnaSentry_parametric <- function(x, ...) {
  cols <- c(weibull = "firebrick", exponential = "darkgreen",
            lognormal = "darkorange", loglogistic = "purple")
  km_times <- x$km_fit$time
  fitted_max <- max(vapply(x$curves,
                           function(cr) max(cr$time, na.rm = TRUE), numeric(1)))
  xlim <- c(0, max(c(max(km_times), fitted_max)))
  plot(x$km_fit, lwd = 2, xlim = xlim, xlab = "Time",
       ylab = "Survival probability", ...)
  lty_map <- c(weibull = 1, exponential = 2, lognormal = 3, loglogistic = 4)
  for (dist in names(x$curves)) {
    graphics::lines(x$curves[[dist]]$time, x$curves[[dist]]$survival,
                    col = cols[[dist]], lty = lty_map[[dist]], lwd = 2)
  }
  graphics::legend("topright",
                   legend = c("KM", names(x$curves)),
                   col = c("black", cols[names(x$curves)]),
                   lty = c(1, lty_map[names(x$curves)]), lwd = 2, bty = "n")
  invisible(x)
}
