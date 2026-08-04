# Internal helpers shared across the rnaSentry pipeline. None are exported;
# they exist to keep per-stage statistics (normality checks, effect sizes,
# assay selection) in one place so every stage reports the same measures.

# Session-level environment for lock enforcement. When lock_signature() locks
# a signature, it stores a gene-set fingerprint here; build_signature() checks
# this environment and refuses to proceed if a locked signature already exists,
# preventing silent re-selection after survival analysis.
.rnaSentry_locked_sigs <- new.env(parent = emptyenv())

# Audit-flag ledger. Every stage reports issues through the same schema
# (check, severity, detail, stage) so downstream consumers such as
# generate_report() can assemble a uniform audit trail. The stage is carried
# as an attribute of the empty ledger and materialized into a column on the
# first append.
.new_flags <- function(stage) {
  structure(
    data.frame(check = character(0), severity = character(0),
               detail = character(0), stringsAsFactors = FALSE),
    stage = stage
  )
}

# Cox models quote non-syntactic predictor names (for example the gene symbol
# "1-Mar") with backticks; restore the original identifier so that coefficient
# names always match the analysis-matrix rownames.
.strip_backticks <- function(x) {
  sub("^`(.*)`$", "\\1", x)
}

# Backquote non-syntactic identifiers (e.g. "RP11-28F1.2") so that a formula
# built with paste/reformulate parses. Syntactic names are returned unchanged.
.backquote_names <- function(x) {
  vapply(x, function(nm) {
    if (!nzchar(nm)) return(nm)
    if (identical(make.names(nm), nm)) nm else paste0("`", nm, "`")
  }, character(1), USE.NAMES = FALSE)
}


.add_flag <- function(flags, check, severity, detail, stage = NULL) {
  st <- stage
  if (is.null(st)) st <- attr(flags, "stage")
  if (is.null(st)) st <- flags$stage[1]
  rbind(flags, data.frame(check = check, severity = severity,
                          detail = detail, stage = st,
                          stringsAsFactors = FALSE))
}

# Detect expression matrices that look already log-transformed but are not
# stored in a preferred assay slot. The pipeline log2-transforms the first
# assay when no "logcounts"/"vst" assay exists, so pre-scaled data sitting in
# a "counts"-named slot would be double-logged. Two signals must agree before
# we warn, so low-depth or filtered-but-correct integer count data is not
# false-flagged: (1) the values are non-integer, and (2) the maximum is below
# log2_scale_max (a ceiling far too small for raw bulk read counts).
# Returns a list with "possible" (logical) and a ready-to-report "message".
.log_detect_log_scaled <- function(mat) {
  log2_scale_max <- 40
  non_integer <- any(abs(mat - round(mat)) > 1e-8, na.rm = TRUE)
  possible <- isTRUE(non_integer) &&
    isTRUE(max(mat, na.rm = TRUE) < log2_scale_max)
  message <- if (possible) {
    sprintf(paste0("The first assay looks already log-transformed (non-integer ",
                   "values with a maximum below %d). rnaSentry will log2-transform ",
                   "it again, double-logging expression. Rename the assay to ",
                   "\"logcounts\" or supply raw integer counts."), log2_scale_max)
  } else {
    NA_character_
  }
  list(possible = possible, message = message)
}

# Select the analysis assay for PCA/confounder work. Prefers an existing
# "logcounts" assay, then "vst", otherwise computes log2(counts + 1) on the
# first assay. Returns a list with the matrix, the assay name used, and a
# flag/message from .log_detect_log_scaled() so stages can surface a
# double-log warning through the flag ledger.
.get_analysis_matrix <- function(se) {
  anames <- SummarizedExperiment::assayNames(se)
  preferred <- c("logcounts", "vst")
  chosen <- preferred[preferred %in% anames]
  if (length(chosen) > 0) {
    mat <- as.matrix(SummarizedExperiment::assay(se, chosen[1]))
    return(list(mat = mat, assay = chosen[1],
                log_scaled_possible = FALSE,
                log_scaled_msg = NA_character_))
  }
  mat <- as.matrix(SummarizedExperiment::assay(se, 1))
  if (!is.numeric(mat)) {
    stop("The first assay is not numeric; provide count or normalized expression data.",
         call. = FALSE)
  }
  detect <- .log_detect_log_scaled(mat)
  list(mat = log2(mat + 1), assay = "log2(counts+1)",
       log_scaled_possible = detect$possible,
       log_scaled_msg = detect$message)
}

# Drop rows with zero variance or any non-finite value. Returns a list with
# the filtered matrix and a per-reason count of dropped rows.
.drop_nonvariable <- function(mat) {
  finite_rows <- apply(is.finite(mat), 1, all)
  n_nonfinite <- sum(!finite_rows)

  keep_rsd <- if (ncol(mat) > 1) {
    col_scale <- apply(mat, 1, function(r) {
      if (any(!is.finite(r))) return(NA_real_)
      stats::sd(r)
    })
    is.finite(col_scale) & col_scale > 0
  } else {
    rep(TRUE, nrow(mat))
  }
  n_constant <- sum(!keep_rsd)

  keep <- finite_rows & keep_rsd
  if (!all(keep)) {
    mat <- mat[keep, , drop = FALSE]
  }
  list(mat = mat, n_nonfinite = n_nonfinite, n_constant = n_constant)
}

# Approximate-normality check via Shapiro-Wilk on residuals (if groups are
# given) or on the vector itself. FALSE when the test is unreliable (too few
# or too many observations) or when all values are identical.
.is_approx_normal <- function(x, groups = NULL, normality_p = 0.05) {
  if (is.null(groups)) {
    resid <- x
  } else {
    resid <- x - stats::ave(x, groups, FUN = function(v) mean(v, na.rm = TRUE))
  }
  resid <- resid[is.finite(resid)]
  n <- length(resid)
  if (n < 3 || n > 5000 || stats::sd(resid) == 0) {
    return(FALSE)
  }
  p <- tryCatch(
    stats::shapiro.test(resid)$p.value,
    error = function(e) NA_real_
  )
  !is.na(p) && p > normality_p
}

# Cramer's V for a contingency table. Returns NA for degenerate tables.
.cramers_v <- function(tab) {
  if (any(dim(tab) < 2)) return(NA_real_)
  chi <- tryCatch(
    stats::chisq.test(tab)$statistic,
    error = function(e) NA_real_
  )
  n <- sum(tab)
  k <- min(dim(tab)) - 1
  if (is.na(chi) || n == 0 || k == 0) return(NA_real_)
  sqrt(as.numeric((chi / n) / k))
}

# Eta-squared from the ANOVA table of a single-factor model.
.eta_squared_anova <- function(model) {
  tab <- stats::anova(model)
  total <- sum(tab$Sum)
  if (!is.finite(total) || total == 0) return(NA_real_)
  tab$Sum[1] / total
}

# Epsilon-squared for a Kruskal-Wallis test (H statistic, k groups, n obs).
# Clamped to >= 0.
.epsilon_squared_kw <- function(h, k, n) {
  if (!is.finite(h) || n <= k) return(NA_real_)
  max(0, (h - k + 1) / (n - k))
}

# Linear-regression association test of a numeric PC score against a numeric
# predictor. Reports the F statistic, degrees of freedom, p-value and
# R-squared effect size. Returns a one-row data.frame.
.regression_pc_test <- function(pc, x, alpha = 0.05) {
  if (stats::sd(x, na.rm = TRUE) == 0) {
    return(data.frame(
      test = "none", statistic = NA_real_, df = NA_character_,
      p_value = NA_real_, effect_size = NA_real_,
      effect_size_type = NA_character_, flagged = NA,
      stringsAsFactors = FALSE
    ))
  }
  model <- stats::lm(pc ~ x)
  f <- summary(model)$fstatistic
  p <- stats::pf(f[1], f[2], f[3], lower.tail = FALSE)
  data.frame(
    test = "lm", statistic = unname(f[1]),
    df = sprintf("%d, %d", f[2], f[3]),
    p_value = unname(p),
    effect_size = unname(summary(model)$r.squared),
    effect_size_type = "r_squared",
    flagged = !is.na(p) && p < alpha,
    stringsAsFactors = FALSE
  )
}

# One-way association test of a numeric PC score against a categorical
# group. ANOVA when the PC scores are approximately normal within groups,
# Kruskal-Wallis otherwise, with eta-squared or epsilon-squared effect size.
# Returns a one-row data.frame.
.oneway_pc_test <- function(pc, group, alpha = 0.05) {
  group <- droplevels(as.factor(group))
  levels_n <- length(unique(group))
  if (levels_n < 2) {
    return(data.frame(
      test = "none", statistic = NA_real_, df = NA_character_,
      p_value = NA_real_, effect_size = NA_real_,
      effect_size_type = NA_character_, flagged = NA,
      stringsAsFactors = FALSE
    ))
  }
  if (levels_n == 2 || .is_approx_normal(pc, groups = group)) {
    model <- stats::aov(pc ~ group)
    f <- summary(model)[[1]]
    if (is.null(f$`F value`) || length(f$Df) < 2 ||
        !is.finite(f$`F value`[1])) {
      return(data.frame(
        test = "aov", statistic = NA_real_, df = NA_character_,
        p_value = NA_real_, effect_size = NA_real_,
        effect_size_type = "eta_squared", flagged = FALSE,
        stringsAsFactors = FALSE
      ))
    }
    df <- sprintf("%d, %d", f$Df[1], f$Df[2])
    data.frame(
      test = "aov", statistic = unname(f$`F value`[1]), df = df,
      p_value = unname(f$`Pr(>F)`[1]),
      effect_size = .eta_squared_anova(model),
      effect_size_type = "eta_squared",
      flagged = !is.na(f$`Pr(>F)`[1]) && f$`Pr(>F)`[1] < alpha,
      stringsAsFactors = FALSE
    )
  } else {
    kw <- stats::kruskal.test(pc, group)
    h <- unname(kw$statistic)
    k <- length(unique(group))
    data.frame(
      test = "kruskal.test", statistic = h, df = as.character(k - 1),
      p_value = unname(kw$p.value),
      effect_size = .epsilon_squared_kw(h, k, length(pc)),
      effect_size_type = "epsilon_squared",
      flagged = !is.na(kw$p.value) && kw$p.value < alpha,
      stringsAsFactors = FALSE
    )
  }
}
