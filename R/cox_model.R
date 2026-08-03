# Adjusted Cox proportional-hazards modeling of a signature.

#' Adjusted Cox proportional-hazards model for a signature
#'
#' \code{cox_model()} fits a single multivariate Cox model of the survival
#' outcome on the signature genes together with any design terms recorded in
#' the signature at build time (\code{sig$design_terms}) and any extra
#' \code{confounders} supplied here. It is the companion to the gene-only
#' scoring model in \code{\link{build_signature}}: the signature's risk score
#' deliberately ignores the design covariates so it stays portable, while this
#' stage supplies inference \emph{adjusted} for them on a cohort that records
#' them.
#'
#' For every fitted term the model reports the hazard ratio with a Wald
#' confidence interval and p-value. For the signature genes it also reports a
#' Benjamini-Hochberg adjusted p-value across the genes, so that gene-level
#' multiplicity is controlled. Proportional-hazards diagnostics
#' (\code{\link[survival]{cox.zph}}) are run on every term and the GLOBAL test;
#' a violation flags the model.
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
#' @param confounders Character vector of additional \code{colData(se)} columns
#'   to adjust for, beyond \code{sig$design_terms}. Defaults to
#'   \code{character(0)}.
#'
#' @return An object of class \code{"rnaSentry_cox_model"} (a list) with
#'   elements:
#'   \describe{
#'     \item{genes}{The signature genes in the model.}
#'     \item{terms}{All fitted terms: the genes, then \code{sig$design_terms},
#'       then \code{confounders}.}
#'     \item{formula}{The model formula.}
#'     \item{fit}{The \code{coxph} fit object.}
#'     \item{coef_table}{Data.frame with one row per fitted term: \code{term},
#'       \code{coefficient}, \code{se}, \code{HR}, \code{HR_low},
#'       \code{HR_high} and \code{p}.}
#'     \item{gene_summary}{Data.frame of the signature-gene rows of
#'       \code{coef_table} with an added BH-adjusted \code{adj_p}.}
#'     \item{zph_summary}{Data.frame from \code{cox.zph} with one row per term
#'       plus \code{GLOBAL}; \code{NULL} when the test could not be run.}
#'     \item{concordance}{Harrell's concordance index of the fitted model.}
#'     \item{ph_violated}{Logical; \code{TRUE} when any term (including
#'       \code{GLOBAL}) violates proportional hazards at alpha 0.05.}
#'     \item{n, events}{Number of samples and events used by the fit.}
#'     \item{sig_locked}{Whether the input signature was locked.}
#'     \item{created}{Provenance timestamp.}
#'     \item{flags}{Data.frame of issues raised (missing covariate values,
#'       non-estimable genes, proportional-hazards violations).}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(7)
#' counts <- matrix(rpois(450, lambda = 500), nrow = 30, ncol = 15,
#'                  dimnames = list(paste0("gene", 1:30), paste0("S", 1:15)))
#' sig_expr <- colMeans(counts[1:4, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(15, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(15, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  batch = factor(rep(c("B1", "B2"),
#'                                                     length.out = 15)),
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#' sig <- build_signature(se, "time", "event", top_n = 4, repeats = 1,
#'                        folds = 2, seed = 1, design_terms = "batch")
#' cm <- cox_model(sig, se)
#' cm
#' plot(cm)
#'
#' @export
cox_model <- function(sig, se, confounders = character(0)) {
  if (!inherits(sig, "rnaSentry_signature")) {
    stop("'sig' must be an rnaSentry_signature object.", call. = FALSE)
  }
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  genes <- sig$genes
  if (!is.character(genes) || length(genes) == 0) {
    stop("'sig' has no genes to model.", call. = FALSE)
  }
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
  if (!is.character(confounders) || anyNA(confounders)) {
    stop("'confounders' must be a character vector.", call. = FALSE)
  }
  if (anyDuplicated(confounders)) {
    stop("'confounders' must not contain duplicates.", call. = FALSE)
  }
  cov_terms <- c(sig$design_terms, confounders)
  surv_cols <- c(sig$time_col, sig$event_col)
  collide_cov <- intersect(cov_terms, surv_cols)
  if (length(collide_cov) > 0) {
    stop(sprintf("Covariate(s) %s are the survival outcome columns; confounders must be separate from the time/event columns.",
                 paste(collide_cov, collapse = ", ")), call. = FALSE)
  }
  dup_cov <- cov_terms[duplicated(cov_terms)]
  if (length(dup_cov) > 0) {
    stop(sprintf("Confounder(s) duplicate a design term: %s.",
                 paste(dup_cov, collapse = ", ")), call. = FALSE)
  }
  overlap_genes <- intersect(confounders, genes)
  if (length(overlap_genes) > 0) {
    stop(sprintf("'confounders' must not contain signature genes: %s.",
                 paste(overlap_genes, collapse = ", ")), call. = FALSE)
  }
  missing_cov <- cov_terms[!cov_terms %in% colnames(cd)]
  if (length(missing_cov) > 0) {
    stop(sprintf("colData(se) has no covariate column(s): %s.",
                 paste(missing_cov, collapse = ", ")), call. = FALSE)
  }

  flags <- .new_flags("cox_model")

  mat <- .get_analysis_matrix(se)$mat
  missing_genes <- setdiff(genes, rownames(mat))
  if (length(missing_genes) > 0) {
    stop(sprintf("The following signature gene(s) are missing from the SummarizedExperiment: %s.",
                 paste(missing_genes, collapse = ", ")), call. = FALSE)
  }

  cd_df <- as.data.frame(cd)
  n_missing_cov <- 0L
  for (t in cov_terms) {
    x <- cd_df[[t]]
    u <- unique(x[!is.na(x)])
    if (length(u) < 2) {
      stop(sprintf("Covariate '%s' has no variation across samples.", t),
           call. = FALSE)
    }
    if (anyNA(x)) {
      n_missing_cov <- n_missing_cov + sum(is.na(x))
    }
  }
  if (n_missing_cov > 0) {
    flags <- .add_flag(flags, "missing_covariates", "info",
                       sprintf("Covariate values are missing for %d sample(s); those samples are excluded from the fit.",
                               n_missing_cov))
  }

  expr <- as.data.frame(t(mat[genes, , drop = FALSE]), check.names = FALSE)
  d <- cbind(data.frame(time = time_vec, event = event_vec), expr)
  for (t in cov_terms) d[[t]] <- cd_df[[t]]

  fit <- tryCatch(
    suppressWarnings(
      survival::coxph(survival::Surv(time, event) ~ ., data = d)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    stop("The adjusted Cox model failed to fit; check for insufficient events or collinear terms.",
         call. = FALSE)
  }

  sm <- summary(fit)$coefficients
  coef_table <- data.frame(
    term = .strip_backticks(rownames(sm)),
    coefficient = unname(sm[, "coef"]),
    se = unname(sm[, "se(coef)"]),
    HR = unname(sm[, "exp(coef)"]),
    p = unname(sm[, "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
  ci <- stats::confint(fit)
  coef_table$HR_low <- exp(ci[, 1])
  coef_table$HR_high <- exp(ci[, 2])

  nonfinite_est <- !is.finite(coef_table$coefficient) |
    !is.finite(coef_table$se) | !is.finite(coef_table$HR) |
    !is.finite(coef_table$HR_low) | !is.finite(coef_table$HR_high)
  if (any(nonfinite_est)) {
    flags <- .add_flag(flags, "non_estimable_coefficients", "warning",
                       sprintf("Term(s) %s have non-finite coefficient or hazard-ratio confidence intervals (possible rank deficiency or exact collinearity); interpret with caution.",
                               paste(coef_table$term[nonfinite_est],
                                     collapse = ", ")))
  }

  gene_rows <- coef_table$term %in% genes
  gene_summary <- coef_table[gene_rows, , drop = FALSE]
  gene_summary$adj_p <- stats::p.adjust(gene_summary$p, method = "BH")
  missing_est <- setdiff(genes, coef_table$term)
  if (length(missing_est) > 0) {
    flags <- .add_flag(flags, "non_estimable_genes", "warning",
                       sprintf("Signature gene(s) %s had no estimable coefficient in the joint model and are absent from the coefficient table.",
                               paste(missing_est, collapse = ", ")))
  }

  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  if (is.null(zph)) {
    zph_summary <- NULL
    ph_violated <- FALSE
    flags <- .add_flag(flags, "ph_test_failed", "warning",
                       "The proportional-hazards diagnostic (cox.zph) could not be computed.")
  } else {
    zt <- zph$table
    pcol <- intersect(c("p", "Pr(>|Chi|)", "Pr(>Chisq)"), colnames(zt))[1]
    rcol <- intersect(c("rho", "rho[1]"), colnames(zt))[1]
    zph_summary <- data.frame(
      term = .strip_backticks(rownames(zt)),
      rho = if (is.na(rcol)) NA_real_ else unname(zt[, rcol]),
      p = unname(zt[, pcol]),
      stringsAsFactors = FALSE
    )
    ph_violated <- any(!is.na(zph_summary$p) & zph_summary$p < 0.05)
    if (ph_violated) {
      violators <- zph_summary$term[!is.na(zph_summary$p) &
                                    zph_summary$p < 0.05]
      flags <- .add_flag(flags, "proportional_hazards", "warning",
                         sprintf("Proportional-hazards assumption may be violated (term(s): %s).",
                                 paste(violators, collapse = ", ")))
    }
  }

  concordance <- unname(summary(fit)$concordance[["C"]])

  result <- list(
    genes = genes,
    terms = c(genes, cov_terms),
    formula = stats::reformulate(.backquote_names(c(genes, cov_terms)),
                                 response = "survival::Surv(time, event)"),
    fit = fit,
    coef_table = coef_table,
    gene_summary = gene_summary,
    zph_summary = zph_summary,
    concordance = concordance,
    ph_violated = ph_violated,
    n = fit$n,
    events = fit$nevent,
    sig_locked = isTRUE(sig$locked),
    created = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     " UTC"),
    flags = flags
  )
  class(result) <- "rnaSentry_cox_model"
  result
}

#' @export
print.rnaSentry_cox_model <- function(x, ...) {
  cat(sprintf("rnaSentry Cox model: %d gene(s) + %d covariate term(s), %d events.\n",
              length(x$genes), max(0, length(x$terms) - length(x$genes)),
              x$events))
  cat(sprintf("Concordance: %.3f. Proportional hazards: %s.\n",
              x$concordance, if (x$ph_violated) "violated (see flags)" else "OK"))
  tab <- x$gene_summary
  if (nrow(tab) > 0) {
    cat("Signature genes (HR [95% CI], BH adj. p):\n")
    for (i in seq_len(nrow(tab))) {
      cat(sprintf("  %-14s %.2f [%.2f, %.2f]  p = %s  adj.p = %s\n",
                  tab$term[i], tab$HR[i], tab$HR_low[i], tab$HR_high[i],
                  format.pval(tab$p[i]), format.pval(tab$adj_p[i])))
    }
  }
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
plot.rnaSentry_cox_model <- function(x, ...) {
  tab <- x$gene_summary
  if (nrow(tab) == 0) {
    stop("No estimable signature genes to plot.", call. = FALSE)
  }
  tab <- tab[order(tab$p), ]
  n <- nrow(tab)
  lo <- log(tab$HR_low)
  hi <- log(tab$HR_high)
  xlim <- c(min(c(lo, hi, 0)), max(c(lo, hi, 0)))
  xr <- range(xlim)
  if (xr[1] == xr[2]) xr <- xr + c(-1, 1)
  graphics::plot(log(tab$HR), seq_len(n), xlim = xr, ylim = c(0.5, n + 0.5),
                 yaxt = "n", ylab = "", xlab = "log hazard ratio",
                 pch = 18, cex = 1.1, ...)
  graphics::segments(lo, seq_len(n), hi, seq_len(n))
  graphics::abline(v = 0, lty = 2, col = "grey50")
  graphics::axis(2, at = seq_len(n), labels = tab$term, las = 2, cex.axis = 0.8)
  invisible(x)
}
