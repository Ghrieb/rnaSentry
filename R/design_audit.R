# Audit candidate design covariates before differential-expression modeling.

# Classification helper: numeric variables with more than two distinct values
# are treated as continuous; everything else (factor, character, logical,
# binary numeric) is treated as categorical.
.design_var_type <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("empty")
  if (is.numeric(x) && length(unique(x)) > 2) "numeric" else "categorical"
}

#' Audit candidate design covariates before differential-expression modeling
#'
#' Given a \code{SummarizedExperiment}, \code{design_audit()} examines a set of
#' candidate design variables (technical or biological covariates such as
#' batch, sequencing lane or tissue) and recommends a design formula for
#' downstream differential-expression analysis. It does this without fitting
#' per-gene models: the leading principal component of the expression data
#' serves as a surrogate for global expression structure, and every candidate
#' variable is tested for association with that surrogate. A variable that
#' tracks the surrogate significantly is flagged as a potential confounder,
#' because omitting it from the design risks attributing expression structure
#' to the variable of interest.
#'
#' Three passes are performed:
#' \enumerate{
#'   \item \emph{Confounder scan}: each candidate variable is tested against
#'     the surrogate expression gradient (\code{surrogate_pc}). Categorical
#'     variables are tested with one-way ANOVA (or Kruskal-Wallis when the
#'     scores are not approximately normal) and numeric variables with linear
#'     regression; effect sizes are reported as eta-squared,
#'     epsilon-squared or R-squared.
#'   \item \emph{Redundancy scan}: every pair of variables (including
#'     \code{outcome_col} when supplied) is cross-tested with Cramer's V
#'     (categorical-categorical), Pearson correlation (numeric-numeric) or a
#'     one-way test (mixed), and pairs whose association is significant and
#'     whose effect size reaches \code{redundant_effect_size} are declared
#'     redundant.
#'   \item \emph{Interaction scan}: a likelihood-ratio test compares
#'     \code{other_pc ~ surrogate + v} with \code{other_pc ~ surrogate * v}
#'     for each candidate variable, where \code{other_pc} is the first PC that
#'     is not the surrogate. A significant interaction means the variable's
#'     influence on expression structure changes along the surrogate gradient;
#'     this is reported as interaction-type confounding rather than silently
#'     folded into the recommended formula.
#' }
#'
#' The recommended \code{formula} is \code{~ outcome_col + design_variables},
#' with variables declared redundant against a later-listed variable dropped.
#' Every statistical decision is returned alongside a plain-language
#' \code{rationale}; nothing is dropped or included silently.
#'
#' @param se A \code{SummarizedExperiment} with a count or normalized
#'   expression assay.
#' @param design_vars Character vector naming candidate design variables, all
#'   of which must be columns of \code{colData(se)}.
#' @param outcome_col Optional character. Name of the biological variable of
#'   interest (for example condition), also a column of \code{colData(se)}.
#'   When supplied it is placed first in the recommended formula, is never
#'   dropped as redundant, and is included in the pairwise redundancy scan.
#'   Defaults to \code{NULL}.
#' @param surrogate_pc Integer. Which principal component to use as the
#'   surrogate for global expression structure. Defaults to \code{1} (the
#'   leading PC); pass a different integer to override.
#' @param alpha Numeric in (0, 1). Significance threshold for the confounder
#'   and pairwise association tests. Defaults to \code{0.05}.
#' @param interaction_alpha Numeric in (0, 1). Significance threshold for the
#'   interaction likelihood-ratio test. Defaults to \code{alpha}.
#' @param redundant_effect_size Numeric in (0, 1). Minimum effect size
#'   (Cramer's V, absolute Pearson r or eta-squared) for a significant
#'   pairwise association to be declared redundant. Defaults to \code{0.8}.
#'
#' @return An object of class \code{"rnaSentry_design"} (a list) with
#'   elements:
#'   \describe{
#'     \item{formula}{The recommended design formula as a
#'       \code{stats::formula} object.}
#'     \item{formula_text}{Character version of \code{formula}.}
#'     \item{rationale}{Character vector of plain-language statements
#'       justifying every automated decision.}
#'     \item{confounder_table}{Data.frame with one row per candidate variable
#'       (columns \code{variable}, \code{type}, \code{n}, \code{test},
#'       \code{statistic}, \code{df}, \code{p_value}, \code{adj_p},
#'       \code{effect_size}, \code{effect_size_type}, \code{flagged}). The
#'       \code{p_value} column holds the raw association p-value; \code{adj_p}
#'       is the Benjamini-Hochberg adjusted p-value across all tested
#'       candidate variables, and \code{flagged} is based on \code{adj_p}.}
#'     \item{interaction_tested}{Logical; whether the interaction LRT could be
#'       performed (it requires at least two PCs and enough samples).}
#'     \item{interaction_table}{Data.frame with one row per candidate variable
#'       (columns \code{variable}, \code{test}, \code{statistic}, \code{df},
#'       \code{p_value}, \code{flagged}), or \code{NULL} when not tested.}
#'     \item{pairwise_table}{Data.frame with one row per variable pair
#'       (columns \code{var1}, \code{var2}, \code{type_pair}, \code{test},
#'       \code{statistic}, \code{df}, \code{p_value}, \code{effect_size},
#'       \code{effect_size_type}, \code{redundant}), or \code{NULL} when fewer
#'       than two variables are available.}
#'     \item{dropped_vars}{Character vector of design variables removed from
#'       the recommended formula as redundant.}
#'     \item{surrogate_pc}{The surrogate PC actually used (after clamping to
#'       the available number of PCs).}
#'     \item{alpha, interaction_alpha, redundant_effect_size}{The thresholds
#'       used.}
#'     \item{flags}{Data.frame summarizing every issue raised, with columns
#'       \code{check}, \code{severity}, \code{detail} and \code{stage}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(4)
#' counts <- matrix(rpois(360, lambda = 200), nrow = 30, ncol = 12,
#'                   dimnames = list(paste0("gene", 1:30), paste0("S", 1:12)))
#' batch <- factor(rep(c("B1", "B2"), each = 6))
#' condition <- factor(rep(c("A", "B"), times = 6))
#' counts[1:10, batch == "B2"] <- counts[1:10, batch == "B2"] * 3L
#' coldata <- S4Vectors::DataFrame(batch = batch, condition = condition,
#'                                  age = runif(12, 40, 80),
#'                                  row.names = paste0("S", 1:12))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' design_audit(se, design_vars = c("batch", "age"), outcome_col = "condition")
#'
#' @export
design_audit <- function(se, design_vars, outcome_col = NULL,
                         surrogate_pc = 1, alpha = 0.05,
                         interaction_alpha = alpha,
                         redundant_effect_size = 0.8) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  if (!is.character(design_vars) || length(design_vars) == 0 ||
      anyNA(design_vars)) {
    stop("'design_vars' must be a non-empty character vector.", call. = FALSE)
  }
  if (anyDuplicated(design_vars)) {
    stop("'design_vars' must not contain duplicates.", call. = FALSE)
  }
  cd <- SummarizedExperiment::colData(se)
  missing_vars <- design_vars[!design_vars %in% colnames(cd)]
  if (length(missing_vars) > 0) {
    stop(sprintf("colData(se) has no column(s): %s.",
                 paste(missing_vars, collapse = ", ")), call. = FALSE)
  }
  if (!is.null(outcome_col)) {
    if (length(outcome_col) != 1 || !is.character(outcome_col) ||
        anyNA(outcome_col)) {
      stop("'outcome_col' must be NULL or a single character string.",
           call. = FALSE)
    }
    if (!outcome_col %in% colnames(cd)) {
      stop(sprintf("colData(se) has no column '%s'.", outcome_col),
           call. = FALSE)
    }
    if (outcome_col %in% design_vars) {
      stop("'outcome_col' must not also be listed in 'design_vars'.",
           call. = FALSE)
    }
  }
  if (!is.numeric(surrogate_pc) || length(surrogate_pc) != 1 ||
      is.na(surrogate_pc) || surrogate_pc < 1 ||
      surrogate_pc != round(surrogate_pc)) {
    stop("'surrogate_pc' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(interaction_alpha) || length(interaction_alpha) != 1 ||
      is.na(interaction_alpha) || interaction_alpha <= 0 ||
      interaction_alpha >= 1) {
    stop("'interaction_alpha' must be a single number in (0, 1).",
         call. = FALSE)
  }
  if (!is.numeric(redundant_effect_size) ||
      length(redundant_effect_size) != 1 || is.na(redundant_effect_size) ||
      redundant_effect_size <= 0 || redundant_effect_size >= 1) {
    stop("'redundant_effect_size' must be a single number in (0, 1).",
         call. = FALSE)
  }

  flags <- .new_flags("design_audit")

  pa <- pca_audit(se, top_n_pcs = max(surrogate_pc + 1, 2))
  n_pcs <- ncol(pa$scores)
  s_idx <- min(surrogate_pc, n_pcs)
  if (s_idx != surrogate_pc) {
    flags <- .add_flag(flags, "surrogate_pc_clamped", "info",
                       sprintf("surrogate_pc = %d exceeds the %d available PC(s); using PC%d.",
                               surrogate_pc, n_pcs, s_idx))
  }
  surrogate <- pa$scores[[s_idx]]
  other_idx <- which(seq_len(n_pcs) != s_idx)[1]
  other_pc <- if (!is.na(other_idx)) pa$scores[[other_idx]] else NULL

  # ---- confounder scan ------------------------------------------------------
  types <- vapply(design_vars, function(v) .design_var_type(cd[[v]]),
                  character(1))
  conf_rows <- vector("list", length(design_vars))
  for (i in seq_along(design_vars)) {
    v <- design_vars[i]
    x <- cd[[v]]
    ok <- !is.na(x)
    base <- data.frame(variable = v, type = types[i], n = sum(ok),
                       stringsAsFactors = FALSE)
    if (sum(ok) < 2 || types[i] == "empty") {
      flags <- .add_flag(flags, "variable_no_data", "warning",
                         sprintf("Design variable '%s' has no usable values; not tested.",
                                 v))
      conf_rows[[i]] <- cbind(base, data.frame(
        test = "none", statistic = NA_real_, df = NA_character_,
        p_value = NA_real_, effect_size = NA_real_,
        effect_size_type = NA_character_, flagged = FALSE,
        stringsAsFactors = FALSE))
      next
    }
    s <- surrogate[ok]
    row <- if (types[i] == "numeric") {
      .regression_pc_test(s, as.numeric(x[ok]), alpha = alpha)
    } else {
      .oneway_pc_test(s, x[ok], alpha = alpha)
    }
    conf_rows[[i]] <- cbind(base, row)
  }
  confounder_table <- do.call(rbind, conf_rows)
  rownames(confounder_table) <- NULL
  confounder_table$adj_p <- stats::p.adjust(confounder_table$p_value,
                                            method = "BH")
  confounder_table$flagged <- !is.na(confounder_table$adj_p) &
    confounder_table$adj_p < alpha

  # ---- redundancy scan ------------------------------------------------------
  all_vars <- unique(c(outcome_col, design_vars))
  pair_rows <- list()
  if (length(all_vars) >= 2) {
    combos <- utils::combn(all_vars, 2)
    pair_rows <- lapply(seq_len(ncol(combos)), function(j) {
      v1 <- combos[1, j]
      v2 <- combos[2, j]
      x1 <- cd[[v1]]
      x2 <- cd[[v2]]
      ok <- !is.na(x1) & !is.na(x2)
      base <- data.frame(var1 = v1, var2 = v2, stringsAsFactors = FALSE)
      if (sum(ok) < 3) {
        return(cbind(base, data.frame(
          type_pair = NA_character_, test = "none", statistic = NA_real_,
          df = NA_character_, p_value = NA_real_, effect_size = NA_real_,
          effect_size_type = NA_character_, redundant = FALSE,
          stringsAsFactors = FALSE)))
      }
      t1 <- .design_var_type(x1)
      t2 <- .design_var_type(x2)
      if (t1 == "categorical" && t2 == "categorical") {
        tab <- table(factor(x1[ok]), factor(x2[ok]))
        chi <- tryCatch(stats::chisq.test(tab), error = function(e) NULL)
        p <- if (is.null(chi)) NA_real_ else unname(chi$p.value)
        V <- .cramers_v(tab)
        cbind(base, data.frame(
          type_pair = "categorical x categorical", test = "chisq.test",
          statistic = if (is.null(chi)) NA_real_ else unname(chi$statistic),
          df = if (is.null(chi)) NA_character_ else paste(chi$parameter, collapse = ", "),
          p_value = p, effect_size = V, effect_size_type = "cramers_v",
          redundant = !is.na(p) && p < alpha && !is.na(V) &&
            V >= redundant_effect_size, stringsAsFactors = FALSE))
      } else if (t1 == "numeric" && t2 == "numeric") {
        ct <- tryCatch(
          stats::cor.test(as.numeric(x1[ok]), as.numeric(x2[ok]),
                          method = "pearson"),
          error = function(e) NULL
        )
        r <- if (is.null(ct)) NA_real_ else unname(ct$estimate)
        p <- if (is.null(ct)) NA_real_ else unname(ct$p.value)
        cbind(base, data.frame(
          type_pair = "numeric x numeric", test = "cor.test",
          statistic = if (is.null(ct)) NA_real_ else unname(ct$statistic),
          df = as.character(sum(ok) - 2), p_value = p, effect_size = r,
          effect_size_type = "pearson_r",
          redundant = !is.na(p) && p < alpha && !is.na(r) &&
            abs(r) >= redundant_effect_size, stringsAsFactors = FALSE))
      } else {
        num <- if (t1 == "numeric") as.numeric(x1[ok]) else as.numeric(x2[ok])
        cat <- if (t1 == "categorical") x1[ok] else x2[ok]
        row <- .oneway_pc_test(num, cat, alpha = alpha)
        cbind(base, data.frame(
          type_pair = "categorical x numeric", test = row$test,
          statistic = row$statistic, df = row$df, p_value = row$p_value,
          effect_size = row$effect_size, effect_size_type = row$effect_size_type,
          redundant = !is.na(row$p_value) && row$p_value < alpha &&
            !is.na(row$effect_size) &&
            row$effect_size >= redundant_effect_size, stringsAsFactors = FALSE))
      }
    })
  }
  pairwise_table <- if (length(pair_rows) > 0) {
    out <- do.call(rbind, pair_rows)
    rownames(out) <- NULL
    out
  } else {
    NULL
  }

  # ---- interaction scan -----------------------------------------------------
  interaction_tested <- FALSE
  interaction_table <- NULL
  if (!is.null(other_pc) && length(other_pc) >= 6) {
    interaction_tested <- TRUE
    rows <- lapply(design_vars, function(v) {
      x <- cd[[v]]
      ok <- !is.na(x)
      none <- data.frame(
        variable = v, test = "none", statistic = NA_real_,
        df = NA_character_, p_value = NA_real_, flagged = FALSE,
        stringsAsFactors = FALSE)
      if (sum(ok) < 6) return(none)
      s <- surrogate[ok]
      t <- other_pc[ok]
      g <- x[ok]
      if (.design_var_type(g) == "categorical") {
        gf <- as.factor(g)
        if (nlevels(gf) < 2 || nlevels(gf) > sum(ok) - 3) return(none)
        red <- stats::lm(t ~ s + gf)
        full <- stats::lm(t ~ s * gf)
      } else {
        red <- stats::lm(t ~ s + g)
        full <- stats::lm(t ~ s * g)
      }
      lrt <- tryCatch(stats::anova(red, full), error = function(e) NULL)
      if (is.null(lrt) || is.na(lrt$F[2])) return(none)
      data.frame(
        variable = v, test = "lrt",
        statistic = unname(lrt$F[2]),
        df = sprintf("%d, %d", lrt$Df[2], lrt$Res.Df[2]),
        p_value = unname(lrt$`Pr(>F)`[2]),
        flagged = !is.na(lrt$`Pr(>F)`[2]) &&
          lrt$`Pr(>F)`[2] < interaction_alpha,
        stringsAsFactors = FALSE)
    })
    interaction_table <- do.call(rbind, rows)
    rownames(interaction_table) <- NULL
  }

  # ---- recommended formula --------------------------------------------------
  drop_vars <- character(0)
  if (!is.null(pairwise_table)) {
    for (j in seq_len(nrow(pairwise_table))) {
      pr <- pairwise_table[j, ]
      if (!isTRUE(pr$redundant)) next
      pos1 <- match(pr$var1, design_vars)
      pos2 <- match(pr$var2, design_vars)
      if (!is.na(pos1) && !is.na(pos2)) {
        drop_later <- if (pos1 > pos2) pr$var1 else pr$var2
        drop_vars <- unique(c(drop_vars, drop_later))
      }
    }
  }
  if (length(drop_vars) > 0) {
    for (dv in drop_vars) {
      pr <- pairwise_table[pairwise_table$var1 == dv |
                             pairwise_table$var2 == dv, , drop = FALSE]
      pr <- pr[isTRUE(pr$redundant), , drop = FALSE]
      partner <- if (nrow(pr) > 0 && pr$var1[1] == dv) pr$var2[1] else pr$var1[1]
      flags <- .add_flag(flags, "redundant_variable", "warning",
                         sprintf("Design variable '%s' is redundant with '%s' (effect size %.2f, p = %.3g) and was dropped from the formula.",
                                 dv, partner, pr$effect_size[1], pr$p_value[1]))
    }
  }
  constant_vars <- design_vars[vapply(design_vars, function(v) {
    x <- cd[[v]]
    length(unique(x[!is.na(x)])) < 2
  }, logical(1))]
  if (length(constant_vars) > 0) {
    flags <- .add_flag(flags, "variable_constant", "warning",
                       sprintf("Design variable(s) %s are constant or have no usable values; excluded from the recommended formula.",
                               paste(constant_vars, collapse = ", ")))
  }
  terms <- c(outcome_col,
             design_vars[!design_vars %in% c(drop_vars, constant_vars)])
  formula_text <- paste(c("~", if (length(terms) == 0) "1" else
    paste(terms, collapse = " + ")), collapse = " ")
  formula_obj <- stats::as.formula(formula_text)

  # ---- rationale ------------------------------------------------------------
  rationale <- character(0)
  rationale <- c(rationale, sprintf(
    "Surrogate expression gradient: PC%d (%.1f%% of variance).",
    s_idx, pa$percent_variance[s_idx]))
  flagged <- confounder_table$flagged[!is.na(confounder_table$flagged)]
  if (any(flagged)) {
    bad <- confounder_table$variable[confounder_table$flagged &
                                       !is.na(confounder_table$flagged)]
    rationale <- c(rationale, sprintf(
      "%d design variable(s) associate with the surrogate gradient at BH-adjusted alpha = %.2f and are flagged as potential confounders: %s.",
      length(bad), alpha, paste(bad, collapse = ", ")))
  } else {
    rationale <- c(rationale, sprintf(
      "No design variable associates with the surrogate gradient at BH-adjusted alpha = %.2f.",
      alpha))
  }
  if (length(drop_vars) > 0) {
    rationale <- c(rationale, sprintf(
      "Variable(s) %s dropped as redundant with a later-listed variable.",
      paste(sprintf("'%s'", drop_vars), collapse = ", ")))
  }
  if (length(constant_vars) > 0) {
    rationale <- c(rationale, sprintf(
      "Variable(s) %s excluded from the recommended formula (constant or no usable values).",
      paste(sprintf("'%s'", constant_vars), collapse = ", ")))
  }
  if (interaction_tested) {
    int_flagged <- interaction_table$flagged[!is.na(interaction_table$flagged)]
    if (any(int_flagged)) {
      int_bad <- interaction_table$variable[interaction_table$flagged &
                                              !is.na(interaction_table$flagged)]
      rationale <- c(rationale, sprintf(
        "Interaction-type confounding detected for variable(s) %s (interaction LRT p < %.2f): consider modeling the variable-condition interaction explicitly.",
        paste(sprintf("'%s'", int_bad), collapse = ", "), interaction_alpha))
    } else {
      rationale <- c(rationale, sprintf(
        "No significant surrogate-variable interaction (LRT p >= %.2f).",
        interaction_alpha))
    }
  } else {
    rationale <- c(rationale,
                   "Interaction scan skipped: fewer than two PCs (or too few samples) available.")
  }
  rationale <- c(rationale, sprintf("Recommended design formula: %s",
                                    formula_text))

  result <- list(
    formula = formula_obj,
    formula_text = formula_text,
    rationale = rationale,
    confounder_table = confounder_table,
    interaction_tested = interaction_tested,
    interaction_table = interaction_table,
    pairwise_table = pairwise_table,
    dropped_vars = drop_vars,
    surrogate_pc = s_idx,
    alpha = alpha,
    interaction_alpha = interaction_alpha,
    redundant_effect_size = redundant_effect_size,
    flags = flags
  )
  class(result) <- "rnaSentry_design"
  result
}

#' @export
print.rnaSentry_design <- function(x, ...) {
  cat(sprintf("rnaSentry design audit: %d candidate variable(s), surrogate PC%d.\n",
              nrow(x$confounder_table), x$surrogate_pc))
  cat(sprintf("Recommended design formula: %s\n", x$formula_text))
  n_conf <- sum(x$confounder_table$flagged, na.rm = TRUE)
  cat(sprintf("%d variable(s) flagged as potential confounders.\n", n_conf))
  if (length(x$dropped_vars) > 0) {
    cat(sprintf("Dropped as redundant: %s.\n",
                paste(x$dropped_vars, collapse = ", ")))
  }
  if (x$interaction_tested) {
    n_int <- sum(x$interaction_table$flagged, na.rm = TRUE)
    cat(sprintf("Interaction LRT run; %d variable(s) with significant surrogate interaction.\n",
                n_int))
  } else {
    cat("Interaction LRT not run (insufficient data).\n")
  }
  if (nrow(x$flags) == 0) {
    cat("No issues flagged.\n")
  } else {
    cat(sprintf("%d issue(s) flagged:\n", nrow(x$flags)))
    for (i in seq_len(nrow(x$flags))) {
      cat(sprintf("  [%s] %s: %s\n", x$flags$severity[i],
                  x$flags$check[i], x$flags$detail[i]))
    }
  }
  invisible(x)
}
