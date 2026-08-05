# Construct a prognostic signature from univariate Cox screening plus
# cross-validated stability reporting.

#' Build a prognostic gene signature from a survival cohort
#'
#' \code{build_signature()} is the core signature-construction stage of the
#' rnaSentry pipeline. From a \code{SummarizedExperiment} with per-sample
#' survival metadata it:
#' \enumerate{
#'   \item \emph{Screens} every variable gene with a univariate Cox
#'     proportional-hazards model against the event indicator, recording the
#'     hazard ratio, Wald p-value and a Benjamini-Hochberg adjusted p-value.
#'   \item \emph{Selects} the signature genes by the most significant
#'     univariate p-values (\code{method = "top_n"}) or by an adjusted
#'     p-value threshold (\code{method = "p_value"}).
#'   \item \emph{Estimates joint coefficients} by refitting a single
#'     multivariate Cox model on the selected genes. Genes whose coefficient
#'     is not finite (rank deficiency) are dropped and the model refit; the
#'     signature is never silently degraded, the drop is reported as a flag.
#'   \item \emph{Reports stability} with repeated, event-stratified k-fold
#'     cross-validation: within each fold the coefficients are re-estimated on
#'     the training samples and the concordance index of the resulting risk
#'     score is evaluated on the held-out samples. The mean and standard
#'     deviation of the per-fold concordance are returned as
#'     \code{cv_summary}. Concordance is Harrell's C in the risk-score
#'     convention: a higher risk score is concordant with an earlier event.
#' }
#'
#' \strong{Reading the cross-validated concordance (important).} The univariate
#' screening that selects the signature runs on the full cohort \emph{before}
#' the cross-validation split, so the CV fold is \emph{not} fully untouched by
#' the selection step. On large candidate panels the resulting concordance is
#' therefore \emph{optimistic} (screening-internal): even data with permuted
#' survival yields a non-trivial null concordance (e.g. \code{~0.8} on a
#' ~21k-gene panel, driven by winner's-curse selection), far above the 0.5 a
#' fully out-of-sample estimate would give. Treat \code{cv_summary} as a
#' screening-internal stability metric and interpret it \emph{relative} to a
#' matched null — a permutation-null or random-gene control — rather than
#' against 0.5 in absolute terms. The only fully independent out-of-sample
#' estimate is \code{\link{validate_external}()} on a separate cohort.
#'
#' The return value has class \code{"rnaSentry_signature"} and is the input
#' expected by the survival modeling and validation stages
#' (\code{km_curve()}, \code{cox_model()}, \code{survival_parametric()},
#' \code{validate_external()}) and by \code{lock_signature()}. The signature
#' is returned \emph{unlocked}; downstream stages that must not be rerun after
#' publication can be protected with \code{lock_signature()}.
#'
#' @param se A \code{SummarizedExperiment} with a count or normalized
#'   expression assay and survival metadata.
#' @param time_col Character. Column of \code{colData(se)} with follow-up
#'   time (numeric, non-negative).
#' @param event_col Character. Column of \code{colData(se)} with the event
#'   indicator (0/1 or logical).
#' @param outcome_col Character. A display label for the outcome, e.g.
#'   \code{"overall_survival"}. This is recorded in the signature and used in
#'   report titles; it does \emph{not} need to name a column of
#'   \code{colData(se)}.
#' @param method Character. Selection rule, one of \code{"top_n"} (default;
#'   take the \code{top_n} genes with the most significant univariate p-value)
#'   or \code{"p_value"} (take every gene with adjusted p-value below
#'   \code{p_threshold}).
#' @param top_n Integer. Number of genes to keep when \code{method =
#'   "top_n"}. Defaults to \code{20}.
#' @param p_threshold Numeric in (0, 1). Adjusted p-value cutoff when
#'   \code{method = "p_value"}. Defaults to \code{0.05}.
#' @param repeats Integer. Number of cross-validation repeats. Defaults to
#'   \code{5}.
#' @param folds Integer. Number of folds per repeat. Defaults to \code{5}.
#' @param seed Optional integer. Seeds the cross-validation fold shuffling for
#'   reproducible results. Seeding is scoped with \code{withr::with_seed()}, so
#'   the caller's global random-number generator state is left unchanged.
#' @param design_terms Character vector of \code{colData(se)} columns to
#'   adjust the univariate screening by. When non-empty and
#'   \code{adjust_for_design = TRUE}, each gene is tested with
#'   \code{Surv ~ gene + design_terms}, so that genes whose survival
#'   association is driven entirely by a confounder (for example batch) are
#'   not selected. The final signature weights (the joint Cox coefficients)
#'   are deliberately estimated from a gene-only model so that the risk score
#'   remains portable to cohorts that do not record the design covariates;
#'   adjusted inference on the chosen signature is provided by
#'   \code{cox_model()}. Defaults to \code{character(0)}.
#' @param adjust_for_design Logical. When \code{TRUE} (default) and
#'   \code{design_terms} is non-empty, univariate screening is adjusted for
#'   the design terms as described above. Setting it to \code{FALSE} runs
#'   unadjusted screening and raises an \code{unadjusted_screening} warning,
#'   because the pipeline's confounder guardrails are then bypassed.
#' @param min_events_per_parameter Numeric. Minimum acceptable number of
#'   events per signature gene. When the cohort has fewer events per final
#'   signature gene, the signature is flagged with an
#'   \code{events_per_parameter} warning, because cross-validated concordance
#'   and hazard ratios are unstable at very low event counts. The default of
#'   \code{5} sits below the classic 10 events-per-variable rule of thumb
#'   (Peduzzi et al., \emph{J Clin Epidemiol} 1996) but at the lower end of
#'   the 5-9 events-per-variable range that Vittinghoff and McCulloch
#'   (\emph{Am J Epidemiol} 2007) argued can be adequate in some settings;
#'   it is a defensible middle choice, not a guarantee of estimability. Raise
#'   it to \code{10} for the stricter convention, or lower it to suppress the
#'   warning on small but well-behaved cohorts.
#' @param BPPARAM Optional \code{BiocParallelParam} from the \pkg{BiocParallel}
#'   package (for example \code{BiocParallel::SnowParam(2)}). When \code{NULL}
#'   (the default) the repeated cross-validation is evaluated serially. When a
#'   non-\code{NULL} backend is supplied, the per-fold evaluations are
#'   distributed with \code{BiocParallel::bplapply()}; \pkg{BiocParallel} is
#'   only a suggested package and is loaded on demand. Parallelism is
#'   strictly opt-in and is bit-identical to the serial path: all randomness
#'   is confined to the fold assignment (seeded with \code{seed}), and the
#'   per-fold Cox fits and concordance scores use no random numbers, so the
#'   returned signature, cross-validation table, and flag ledger are
#'   unchanged regardless of backend.
#'
#' @section Assumptions and limitations:
#' The signature is a linear Cox risk score (larger score = earlier event).
#' The pipeline assumes standard right-censored survival from bulk RNA-seq
#' with raw integer counts or a correctly-named \code{logcounts} assay
#' (pre-scaled data in a \code{counts}-named slot is flagged as
#' \code{possibly_log_scaled}). It does not model competing risks,
#' time-varying covariates, left truncation, single-cell data, or
#' non-survival end points, and it never corrects batch or design effects.
#' In underpowered cohorts (fewer than \code{min_events_per_parameter}
#' events per signature gene) the returned concordance and hazard ratios are
#' unstable and an \code{events_per_parameter} warning is raised; treat such
#' results as hypothesis-generating rather than confirmatory. The discovery
#' stage additionally refuses cohorts with fewer than two events; the
#' downstream survival stages (\code{\link{km_curve}},
#' \code{\link{cox_model}}, \code{\link{survival_parametric}},
#' \code{\link{validate_external}}) use a one-event minimum and instead
#' surface instability through their \code{sparse_events} guardrail.
#'
#' @return An object of class \code{"rnaSentry_signature"} (a list) with
#'   elements:
#'   \describe{
#'     \item{genes}{Character vector of the final signature genes.}
#'     \item{outcome}{The \code{outcome_col} display label.}
#'     \item{time_col, event_col, design_terms}{The metadata columns used.}
#'     \item{screening_terms}{Character vector of the design terms actually
#'       used to adjust univariate screening (empty when screening was
#'       unadjusted).}
#'     \item{cox_stats}{Data.frame with one row per screened gene (columns
#'       \code{gene}, \code{HR}, \code{p}, \code{adj_p}).}
#'     \item{coefficients}{Named numeric vector of joint Cox coefficients for
#'       the signature genes.}
#'     \item{cv_results}{Data.frame with columns \code{repeat_id},
#'       \code{fold} and \code{c_index} (concordance on held-out samples).}
#'     \item{cv_summary}{Named numeric vector \code{c(mean, sd)} of the
#'       per-fold concordance.}
#'     \item{selection}{List with \code{method} and \code{n_genes}.}
#'     \item{locked}{Always \code{FALSE}; set by \code{lock_signature()}.}
#'     \item{created, lock_time}{Provenance timestamps.}
#'     \item{assay_used, n_screened, n_filtered}{Analysis details.}
#'     \item{flags}{Data.frame summarizing every issue raised, with columns
#'       \code{check}, \code{severity}, \code{detail} and \code{stage}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(5)
#' counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
#'                   dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
#' sig_expr <- colMeans(counts[1:5, , drop = FALSE])
#' risk <- scale(sig_expr)[, 1] * 0.4
#' event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
#' censor_time <- rexp(20, rate = 0.02)
#' time <- pmin(event_time, censor_time)
#' event <- as.integer(event_time < censor_time)
#' coldata <- S4Vectors::DataFrame(time = time, event = event,
#'                                  condition = factor(rep(c("A", "B"), 10)),
#'                                  row.names = colnames(counts))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' sig <- build_signature(se, time_col = "time", event_col = "event",
#'                        top_n = 5, repeats = 1, folds = 2, seed = 1)
#' sig
#'
#' @export
build_signature <- function(se, time_col, event_col,
                            outcome_col = "overall_survival",
                            method = c("top_n", "p_value"), top_n = 20,
                             p_threshold = 0.05, repeats = 5, folds = 5,
                             seed = NULL, design_terms = character(0),
                             adjust_for_design = TRUE,
                             min_events_per_parameter = 5,
                             BPPARAM = NULL) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  cd <- SummarizedExperiment::colData(se)
  for (nm in c(time_col, event_col)) {
    if (length(nm) != 1 || !is.character(nm) || !nm %in% colnames(cd)) {
      stop(sprintf("colData(se) has no column '%s'.", nm), call. = FALSE)
    }
  }
  if (length(outcome_col) != 1 || !is.character(outcome_col) ||
      is.na(outcome_col)) {
    stop("'outcome_col' must be a single character string.", call. = FALSE)
  }
  method <- match.arg(method)
  if (!is.numeric(top_n) || length(top_n) != 1 || is.na(top_n) ||
      top_n < 1 || top_n != round(top_n)) {
    stop("'top_n' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(p_threshold) || length(p_threshold) != 1 ||
      is.na(p_threshold) || p_threshold <= 0 || p_threshold >= 1) {
    stop("'p_threshold' must be a single number in (0, 1).", call. = FALSE)
  }
  for (nm in c("repeats", "folds")) {
    val <- get(nm)
    if (!is.numeric(val) || length(val) != 1 || is.na(val) || val < 1 ||
        val != round(val)) {
      stop(sprintf("'%s' must be a single positive integer.", nm),
           call. = FALSE)
    }
  }
  if (!is.null(seed) &&
      (!is.numeric(seed) || length(seed) != 1 || is.na(seed))) {
    stop("'seed' must be NULL or a single number.", call. = FALSE)
  }
  if (!is.character(design_terms) || anyNA(design_terms)) {
    stop("'design_terms' must be a character vector.", call. = FALSE)
  }
  missing_terms <- design_terms[!design_terms %in% colnames(cd)]
  if (length(missing_terms) > 0) {
    stop(sprintf("colData(se) has no design term column(s): %s.",
                 paste(missing_terms, collapse = ", ")), call. = FALSE)
  }
  if (!is.logical(adjust_for_design) || length(adjust_for_design) != 1 ||
      is.na(adjust_for_design)) {
    stop("'adjust_for_design' must be a single TRUE or FALSE.", call. = FALSE)
  }
  if (anyDuplicated(design_terms)) {
    stop("'design_terms' must not contain duplicates.", call. = FALSE)
  }
  if (!is.numeric(min_events_per_parameter) ||
      length(min_events_per_parameter) != 1 ||
      is.na(min_events_per_parameter) || min_events_per_parameter <= 0) {
    stop("'min_events_per_parameter' must be a single positive number.",
         call. = FALSE)
  }
  if (!is.null(BPPARAM) &&
      requireNamespace("BiocParallel", quietly = TRUE) &&
      !methods::is(BPPARAM, "BiocParallelParam")) {
    stop("'BPPARAM' must be NULL or a BiocParallelParam object (for example BiocParallel::SnowParam(2)).",
         call. = FALSE)
  }

  time_vec <- as.numeric(cd[[time_col]])
  if (anyNA(time_vec) || any(time_vec < 0) || any(!is.finite(time_vec))) {
    stop(sprintf("'%s' must be a numeric, non-negative, finite, non-missing column.",
                 time_col), call. = FALSE)
  }
  if (isTRUE(stats::sd(time_vec) == 0)) {
    stop(sprintf("'%s' has no variation in follow-up time.", time_col),
         call. = FALSE)
  }
  event_raw <- cd[[event_col]]
  if (is.logical(event_raw)) {
    event_vec <- as.integer(event_raw)
  } else if (is.numeric(event_raw)) {
    event_vec <- as.integer(event_raw)
  } else {
    stop(sprintf("'%s' must be a 0/1 or logical event indicator.", event_col),
         call. = FALSE)
  }
  if (anyNA(event_vec) || !all(event_vec %in% c(0L, 1L))) {
    stop(sprintf("'%s' must contain only 0/1 values.", event_col),
         call. = FALSE)
  }
  if (sum(event_vec) < 2) {
    stop(sprintf("Fewer than two events in '%s'; survival modeling requires at least two events.",
                 event_col), call. = FALSE)
  }

  flags <- .new_flags("build_signature")

  use_adjust <- length(design_terms) > 0 && adjust_for_design
  if (length(design_terms) > 0) {
    cd_df <- as.data.frame(cd)
    for (t in design_terms) {
      x <- cd_df[[t]]
      u <- unique(x[!is.na(x)])
      if (length(u) < 2) {
        stop(sprintf("Design term '%s' has no variation across samples; drop it or use a variable design.",
                     t), call. = FALSE)
      }
      if (anyNA(x)) {
        flags <- .add_flag(flags, "design_terms_missing", "info",
                           sprintf("Design term '%s' has missing values; samples with missing values are excluded from screening.",
                                   t))
      }
    }
    if (use_adjust) {
      flags <- .add_flag(flags, "adjusted_screening", "info",
                         sprintf("Univariate screening adjusted for design term(s): %s.",
                                 paste(design_terms, collapse = ", ")))
    } else {
      flags <- .add_flag(flags, "unadjusted_screening", "warning",
                         paste0("design_terms were supplied but 'adjust_for_design' is FALSE; ",
                                "univariate screening is unadjusted and confounded genes may be selected."))
    }
  }

  mat_info <- .get_analysis_matrix(se)
  mat <- mat_info$mat
  if (isTRUE(mat_info$log_scaled_possible)) {
    warning(mat_info$log_scaled_msg, call. = FALSE)
    flags <- .add_flag(flags, "possibly_log_scaled", "warning",
                       mat_info$log_scaled_msg)
  }
  if (anyDuplicated(rownames(mat))) {
    stop("Duplicate gene names in the analysis matrix; use unique feature identifiers.",
         call. = FALSE)
  }
  filtered <- .drop_nonvariable(mat)
  mat <- filtered$mat
  n_filtered <- filtered$n_nonfinite + filtered$n_constant
  if (n_filtered > 0) {
    flags <- .add_flag(flags, "gene_filter", "warning",
                       sprintf("Removed %d gene(s) with missing values or zero variance before screening.",
                               n_filtered))
  }
  reserved <- intersect(rownames(mat), c("time", "event"))
  if (length(reserved) > 0) {
    flags <- .add_flag(flags, "gene_name_collision", "warning",
                       sprintf("Gene(s) %s share a name with the internal model columns 'time'/'event' and were excluded from screening and the signature; rename them to keep them in the model.",
                               paste(reserved, collapse = ", ")))
    mat <- mat[setdiff(rownames(mat), reserved), , drop = FALSE]
  }
  if (nrow(mat) < 1) {
    stop("No variable genes remain after filtering; cannot build a signature.",
         call. = FALSE)
  }

  # ---- univariate Cox screening -------------------------------------------
  d0 <- data.frame(time = time_vec, event = event_vec)
  if (use_adjust) {
    for (t in design_terms) d0[[t]] <- as.data.frame(cd)[[t]]
  }
  screen <- lapply(rownames(mat), function(g) {
    d0$g <- mat[g, ]
    form <- if (use_adjust) {
      stats::as.formula(paste("survival::Surv(time, event) ~ g +",
                              paste(.backquote_names(design_terms),
                                    collapse = " + ")))
    } else {
      stats::as.formula("survival::Surv(time, event) ~ g")
    }
    fit <- tryCatch(
      suppressWarnings(
        survival::coxph(form, data = d0)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(data.frame(gene = g, HR = NA_real_, p = NA_real_,
                        adj_p = NA_real_, stringsAsFactors = FALSE))
    }
    sm <- summary(fit)$coefficients
    data.frame(gene = g, HR = unname(sm[, "exp(coef)"][1]),
               p = unname(sm[, "Pr(>|z|)"][1]), adj_p = NA_real_,
               stringsAsFactors = FALSE)
  })
  cox_stats <- do.call(rbind, screen)
  cox_stats$adj_p <- stats::p.adjust(cox_stats$p, method = "BH")
  n_screened <- nrow(cox_stats)

  # ---- selection ----------------------------------------------------------
  if (method == "top_n") {
    sig_p <- cox_stats$p[is.finite(cox_stats$p)]
    if (length(sig_p) == 0) {
      stop("No gene produced a finite univariate Cox p-value; cannot select a signature.",
           call. = FALSE)
    }
    n_take <- min(top_n, length(sig_p))
    if (length(sig_p) < top_n) {
      flags <- .add_flag(flags, "fewer_genes_than_requested", "info",
                         sprintf("Only %d gene(s) had a finite p-value; using all of them instead of top_n = %d.",
                                 length(sig_p), top_n))
    }
    sel <- order(cox_stats$p, na.last = NA)[seq_len(n_take)]
    sel_method <- "top_n"
  } else {
    sel <- which(is.finite(cox_stats$adj_p) & cox_stats$adj_p < p_threshold)
    if (length(sel) == 0) {
      stop(sprintf("No gene passes the adjusted p-value threshold %.3g; relax 'p_threshold' or use method = 'top_n'.",
                   p_threshold), call. = FALSE)
    }
    sel <- sel[order(cox_stats$p[sel])]
    sel_method <- "p_value"
  }
  sel_genes <- rownames(mat)[sel]

  # ---- joint coefficients --------------------------------------------------
  coef_vec <- NULL
  for (iter in seq_len(5)) {
    d <- data.frame(time = time_vec, event = event_vec,
                    t(as.matrix(mat[sel_genes, , drop = FALSE])),
                    check.names = FALSE)
    fit <- tryCatch(
      suppressWarnings(
        survival::coxph(survival::Surv(time, event) ~ ., data = d)
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) break
    b <- stats::coef(fit)
    names(b) <- .strip_backticks(names(b))
    if (length(b) == 0) break
    bad <- !is.finite(b) | !is.finite(exp(b))
    if (any(bad)) {
      dropped <- names(b)[bad]
      flags <- .add_flag(flags, "coefficient_unstable", "warning",
                         sprintf("Gene(s) %s had a non-finite joint Cox coefficient and were dropped from the signature.",
                                 paste(dropped, collapse = ", ")))
      sel_genes <- setdiff(sel_genes, dropped)
      if (length(sel_genes) == 0) break
      next
    }
    coef_vec <- b
    break
  }
  if (is.null(coef_vec) || length(coef_vec) == 0) {
    stop("The joint Cox model for the selected genes did not yield finite ",
         "coefficients; reduce 'top_n' or check for collinear genes.",
         call. = FALSE)
  }
  genes <- names(coef_vec)

  n_events <- sum(event_vec)
  events_per_parameter <- n_events / length(genes)
  if (is.finite(events_per_parameter) &&
      events_per_parameter < min_events_per_parameter) {
    msg <- sprintf(paste0("Only %d event(s) for a %d-gene signature (%.1f events ",
                          "per parameter), below 'min_events_per_parameter' = %d. ",
                          "Cross-validated concordance and hazard ratios are unstable ",
                          "at this event count; consider fewer genes or a larger cohort."),
                   n_events, length(genes), events_per_parameter,
                   min_events_per_parameter)
    warning(msg, call. = FALSE)
    flags <- .add_flag(flags, "events_per_parameter", "warning", msg)
  }

  # ---- repeated stratified cross-validation --------------------------------
  n <- length(time_vec)
  if (folds >= n) {
    stop(sprintf("'folds' (%d) must be smaller than the number of samples (%d); leave-one-out cross-validation cannot estimate a concordance index on single-sample test folds.",
                 folds, n), call. = FALSE)
  }
  flag_fold <- function(flags, check, detail) {
    .add_flag(flags, check, "warning", detail)
  }
  # Fold assignment is the only source of randomness. Drawing every repeat's
  # partitions inside a single seeded scope (when `seed` is supplied) keeps the
  # original semantics - one seed, the RNG advancing across repeats so each
  # repeat gets a distinct partition - while restoring the caller's global RNG
  # state afterwards, which is the point of withr::with_seed().
  make_all_folds <- function() {
    lapply(seq_len(repeats), function(r) {
      fold_ids <- integer(n)
      for (grp in list(which(event_vec == 0L), which(event_vec == 1L))) {
        if (length(grp) == 0) next
        grp_shuffled <- grp[sample.int(length(grp))]
        fold_ids[grp_shuffled] <- rep(seq_len(folds), length.out = length(grp))
      }
      fold_ids
    })
  }
  all_fold_ids <- if (is.null(seed)) {
    make_all_folds()
  } else {
    withr::with_seed(seed, make_all_folds())
  }
  # Per-fold evaluation is deterministic: given a fold assignment it fits a
  # Cox model on the training rows and scores the held-out rows, and it uses
  # no random numbers. That is what makes the opt-in parallel path
  # (BPPARAM) bit-identical to the serial path - the RNG call sequence is
  # confined to make_all_folds() and is therefore unchanged regardless of
  # backend. Each element returns its cv row plus any guardrail flags so the
  # ledger merges in grid order, exactly reproducing the nested-loop order.
  eval_fold <- function(r, f) {
    fold_ids <- all_fold_ids[[r]]
    train <- which(fold_ids != f)
    test <- which(fold_ids == f)
    ci <- NA_real_
    fl <- list()
    if (length(train) < 5 || length(test) < 2 ||
        sum(event_vec[train]) < 2 || sum(event_vec[test]) < 1) {
      fl <- list(c("cv_fold_skipped",
                   sprintf("Fold %d of repeat %d skipped: too few samples or events.",
                           f, r)))
    } else {
      d_tr <- data.frame(time = time_vec[train], event = event_vec[train],
                         t(as.matrix(mat[genes, train, drop = FALSE])),
                         check.names = FALSE)
      fit_tr <- tryCatch(
        suppressWarnings(
          survival::coxph(survival::Surv(time, event) ~ ., data = d_tr)
        ),
        error = function(e) NULL
      )
      if (is.null(fit_tr)) {
        fl <- list(c("cv_fold_failed",
                     sprintf("Cox model failed to fit on fold %d of repeat %d.",
                             f, r)))
      } else {
        b <- stats::coef(fit_tr)
        names(b) <- .strip_backticks(names(b))
        if (length(b) == 0 || any(!is.finite(b))) {
          fl <- list(c("cv_fold_failed",
                       sprintf("Non-finite coefficients on fold %d of repeat %d.",
                               f, r)))
        } else {
          gn <- names(b)
          score_test <- as.vector(t(as.matrix(mat[gn, test, drop = FALSE])) %*% b)
          if (isTRUE(stats::sd(score_test) == 0)) {
            fl <- list(c("cv_fold_failed",
                         sprintf("Constant risk score on fold %d of repeat %d.",
                                 f, r)))
          } else {
            # Risk score: larger = higher hazard = shorter survival, so the
            # concordance must use reverse = TRUE (survival::concordance's
            # default means "larger x => longer survival"). Regression guard:
            # test-statistical_parity.R ("survival::concordance reverse
            # convention satisfies C + C_rev = 1") fails if this is reverted.
            conc <- tryCatch(
              suppressWarnings(
                survival::concordance(
                  survival::Surv(time_vec[test], event_vec[test]) ~ score_test,
                  reverse = TRUE
                )
              ),
              error = function(e) NULL
            )
            ci <- if (is.null(conc)) NA_real_ else
              as.numeric(conc$concordance[1])
            if (!is.finite(ci)) {
              ci <- NA_real_
              fl <- list(c("cv_fold_failed",
                           sprintf("Concordance undefined on fold %d of repeat %d.",
                                   f, r)))
            }
          }
        }
      }
    }
    list(cv_row = data.frame(repeat_id = r, fold = f, c_index = ci),
         flags = fl)
  }

  grid <- expand.grid(r = seq_len(repeats), f = seq_len(folds))
  run_fold <- function(i) eval_fold(grid$r[i], grid$f[i])
  fold_out <- if (is.null(BPPARAM)) {
    lapply(seq_len(nrow(grid)), run_fold)
  } else {
    if (!requireNamespace("BiocParallel", quietly = TRUE)) {
      stop("BiocParallel is not installed; set BPPARAM = NULL to run serially.",
           call. = FALSE)
    }
    BiocParallel::bplapply(seq_len(nrow(grid)), run_fold, BPPARAM = BPPARAM)
  }
  cv_rows <- lapply(fold_out, function(elt) elt$cv_row)
  cv_results <- do.call(rbind, cv_rows)
  rownames(cv_results) <- NULL
  for (elt in fold_out) {
    for (flag in elt$flags) {
      flags <- flag_fold(flags, flag[1], flag[2])
    }
  }
  ci_vals <- cv_results$c_index
  if (all(!is.finite(ci_vals))) {
    stop("Cross-validation could not produce a single concordance value; the model may be degenerate.",
         call. = FALSE)
  }
  n_finite <- sum(is.finite(ci_vals))
  cv_summary <- c(
    mean = mean(ci_vals, na.rm = TRUE),
    sd = if (n_finite < 2) NA_real_ else stats::sd(ci_vals, na.rm = TRUE)
  )

  result <- list(
    genes = genes,
    outcome = outcome_col,
    time_col = time_col,
    event_col = event_col,
    design_terms = design_terms,
    screening_terms = if (use_adjust) design_terms else character(0),
    cox_stats = cox_stats,
    coefficients = coef_vec[genes],
    cv_results = cv_results,
    cv_summary = cv_summary,
    selection = list(method = sel_method, n_genes = length(genes)),
    locked = FALSE,
    created = paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     " UTC"),
    lock_time = NULL,
    assay_used = mat_info$assay,
    n_screened = n_screened,
    n_filtered = n_filtered,
    flags = flags
  )
  class(result) <- "rnaSentry_signature"
  result
}

#' @export
print.rnaSentry_signature <- function(x, ...) {
  cat(sprintf("rnaSentry signature: %d gene(s), outcome '%s'.\n",
              length(x$genes), x$outcome))
  cat(sprintf("Selection: %s (%d gene(s)). CV concordance: mean %.3f, sd %s.\n",
              x$selection$method, x$selection$n_genes,
              x$cv_summary[["mean"]],
              if (is.na(x$cv_summary[["sd"]])) "NA" else
                sprintf("%.3f", x$cv_summary[["sd"]])))
  cat(if (isTRUE(x$locked)) "Locked.\n" else "Not locked.\n")
  if (!is.null(x$flags) && is.data.frame(x$flags) && nrow(x$flags) > 0) {
    cat(sprintf("%d issue(s) flagged:\n", nrow(x$flags)))
    for (i in seq_len(nrow(x$flags))) {
      cat(sprintf("  [%s] %s: %s\n", x$flags$severity[i],
                  x$flags$check[i], x$flags$detail[i]))
    }
  }
  invisible(x)
}
