utils::globalVariables(".data")

#' PCA audit of expression data with optional batch association testing
#'
#' Runs principal component analysis on a variance-stabilized or
#' log-transformed version of the count assay and reports how much variance
#' each principal component (PC) explains. When \code{batch_col} is supplied,
#' each PC score is tested against the batch variable and any PC where batch
#' explains a significant share of variance is flagged. This feeds
#' \code{design_audit()} by surfacing which batch/technical variables
#' are associated with global expression structure.
#'
#' The analysis matrix is chosen by preference from an existing assay named
#' \code{"logcounts"} or \code{"vst"}; if neither exists, \code{log2(counts +
#' 1)} is computed from the first assay. Genes with zero variance or missing
#' values are removed because they carry no information for PCA and would
#' make \code{stats::prcomp} fail; the number removed and the reason are
#' reported in the result (nothing is dropped silently).
#'
#' @param se A \code{SummarizedExperiment} with a count or normalized
#'   expression assay.
#' @param batch_col Optional character. Name of a \code{colData(se)} column
#'   (categorical or numeric) to test against each PC score. Defaults to
#'   \code{NULL} (no batch testing).
#' @param top_n_pcs Integer. Number of leading PCs to include in
#'   \code{scores} and to test against \code{batch_col}. Defaults to \code{5}.
#' @param scale Logical. Whether to scale each gene to unit variance before
#'   PCA (\code{scale. = TRUE}). Defaults to \code{TRUE} so that genes are
#'   comparable on the log-expression scale.
#' @param batch_alpha Numeric in (0, 1). Significance threshold used to flag a
#'   PC whose score associates with \code{batch_col}. Defaults to \code{0.05}.
#'
#' @details
#' Unlike \code{\link{design_audit}}, no multiple-testing correction is
#' applied to the per-PC batch tests: \code{batch_alpha} is applied to each
#' test at its nominal level. This is a deliberate choice, because at most
#' \code{top_n_pcs} (default 5) PCs are tested against a single batch
#' variable, so the correction would be a near no-op; the raw p-values and a
#' Benjamini-Hochberg adjusted p-value (\code{adj_p}) are both reported in
#' \code{batch_tests} so users can apply their own threshold. Candidate
#' confounder screening in \code{\link{design_audit}}, which tests many
#' variables at once, is BH-corrected instead.
#'
#' @return An object of class \code{"rnaSentry_pca"} (a list) with elements:
#'   \describe{
#'     \item{pca}{The \code{stats::prcomp} object.}
#'     \item{percent_variance}{Named numeric vector of the percent of variance
#'       explained by every PC.}
#'     \item{scores}{A data.frame of the first \code{top_n_pcs} PC scores,
#'       with sample IDs as row names.}
#'     \item{col_data}{The \code{colData(se)} columns as a data.frame, kept so
#'       \code{\link{plot_pca_audit}} can color points by a metadata column.}
#'     \item{assay_used}{Name of the assay used for the analysis.}
#'     \item{n_genes_analyzed}{Number of genes retained after filtering.}
#'     \item{n_genes_filtered}{Number of genes removed before PCA.}
#'     \item{filtered_reason}{Character describing why genes were removed
#'       (empty string if none).}
#'     \item{batch_col}{The batch column name used, or \code{NULL}.}
#'     \item{batch_tests}{A data.frame with one row per tested PC (columns
#'       \code{pc}, \code{test}, \code{statistic}, \code{df}, \code{p_value},
#'       \code{adj_p}, \code{effect_size}, \code{effect_size_type},
#'       \code{flagged}), or \code{NULL} when \code{batch_col} is not
#'       supplied.}
#'     \item{flags}{A data.frame summarizing every issue raised, with columns
#'       \code{check}, \code{severity}, and \code{detail}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(1)
#' counts <- matrix(rpois(240, lambda = 200), nrow = 20, ncol = 12,
#'                   dimnames = list(paste0("gene", 1:20), paste0("S", 1:12)))
#' coldata <- DataFrame(condition = rep(c("A", "B"), each = 6),
#'                       batch = factor(rep(c("B1", "B2"), each = 6)),
#'                       row.names = paste0("S", 1:12))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' pca_audit(se, batch_col = "batch")
#'
#' @export
pca_audit <- function(se, batch_col = NULL, top_n_pcs = 5, scale = TRUE,
                      batch_alpha = 0.05) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  if (!is.numeric(top_n_pcs) || length(top_n_pcs) != 1 ||
      is.na(top_n_pcs) || top_n_pcs < 1 || top_n_pcs != round(top_n_pcs)) {
    stop("'top_n_pcs' must be a single positive integer.", call. = FALSE)
  }
  if (!is.numeric(batch_alpha) || length(batch_alpha) != 1 ||
      is.na(batch_alpha) || batch_alpha <= 0 || batch_alpha >= 1) {
    stop("'batch_alpha' must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(batch_col)) {
    if (length(batch_col) != 1 || !is.character(batch_col) ||
        !batch_col %in% colnames(SummarizedExperiment::colData(se))) {
      stop(sprintf("colData(se) has no column '%s'.", batch_col), call. = FALSE)
    }
  }

  mat_info <- .get_analysis_matrix(se)
  mat <- mat_info$mat

  flags <- .new_flags("pca_audit")

  filtered <- .drop_nonvariable(mat)
  mat <- filtered$mat
  n_genes_filtered <- filtered$n_nonfinite + filtered$n_constant
  reasons <- character(0)
  if (filtered$n_nonfinite > 0) {
    reasons <- c(reasons, sprintf("%d with missing/non-finite values",
                                  filtered$n_nonfinite))
  }
  if (filtered$n_constant > 0) {
    reasons <- c(reasons, sprintf("%d with zero variance",
                                  filtered$n_constant))
  }
  filtered_reason <- paste(reasons, collapse = "; ")
  if (n_genes_filtered > 0) {
    flags <- .add_flag(flags, "pca_gene_filter", "warning",
                       sprintf("Removed %d gene(s) before PCA: %s.",
                               n_genes_filtered, filtered_reason))
  }
  if (nrow(mat) < 3) {
    stop(sprintf("Too few variable genes after filtering (n = %d); cannot run PCA.",
                 nrow(mat)), call. = FALSE)
  }
  if (ncol(mat) < 2) {
    stop("Too few samples for PCA; at least two samples are required.",
         call. = FALSE)
  }

  pca <- stats::prcomp(t(mat), center = TRUE, scale. = scale)
  percent_variance <- pca$sdev^2 / sum(pca$sdev^2) * 100
  names(percent_variance) <- paste0("PC", seq_along(percent_variance))

  n_pcs <- min(top_n_pcs, ncol(pca$x))
  scores <- as.data.frame(pca$x[, seq_len(n_pcs), drop = FALSE])
  rownames(scores) <- colnames(se)

  batch_tests <- NULL
  if (!is.null(batch_col)) {
    bvar <- SummarizedExperiment::colData(se)[[batch_col]]
    if (is.numeric(bvar)) {
      bvar_num <- as.numeric(bvar)
      if (stats::sd(bvar_num, na.rm = TRUE) == 0 ||
          length(unique(bvar_num)) < 2) {
        flags <- .add_flag(flags, "batch_single_level", "warning",
                           paste0("Batch column '", batch_col,
                                  "' has no variation; no batch test performed."))
        bvar_num <- NULL
      }
      if (!is.null(bvar_num)) {
        batch_tests <- do.call(rbind, lapply(seq_len(n_pcs), function(i) {
          row <- .regression_pc_test(scores[, i], bvar_num, alpha = batch_alpha)
          cbind(pc = names(scores)[i], row, stringsAsFactors = FALSE)
        }))
      }
    } else {
      bfact <- as.factor(bvar)
      if (length(unique(bfact)) < 2) {
        flags <- .add_flag(flags, "batch_single_level", "warning",
                           paste0("Batch column '", batch_col,
                                  "' has a single level; no batch test performed."))
        bfact <- NULL
      }
      if (!is.null(bfact)) {
        batch_tests <- do.call(rbind, lapply(seq_len(n_pcs), function(i) {
          row <- .oneway_pc_test(scores[, i], bfact, alpha = batch_alpha)
          cbind(pc = names(scores)[i], row, stringsAsFactors = FALSE)
        }))
      }
    }
    if (!is.null(batch_tests)) {
      rownames(batch_tests) <- NULL
      batch_tests$adj_p <- stats::p.adjust(batch_tests$p_value, method = "BH")
      n_flagged <- sum(batch_tests$flagged, na.rm = TRUE)
      if (n_flagged > 0) {
        flags <- .add_flag(flags, "batch_associated_pc", "warning",
                           paste0("Batch '", batch_col,
                                  "' significantly associates with ",
                                  n_flagged, " PC(s): ",
                                  paste(batch_tests$pc[batch_tests$flagged],
                                        collapse = ", "), "."))
      }
    }
  }

  result <- list(
    pca = pca,
    percent_variance = percent_variance,
    scores = scores,
    col_data = as.data.frame(SummarizedExperiment::colData(se)),
    assay_used = mat_info$assay,
    n_genes_analyzed = nrow(mat),
    n_genes_filtered = n_genes_filtered,
    filtered_reason = filtered_reason,
    batch_col = batch_col,
    batch_tests = batch_tests,
    flags = flags
  )
  class(result) <- "rnaSentry_pca"
  result
}

#' @export
print.rnaSentry_pca <- function(x, ...) {
  cat(sprintf("rnaSentry PCA audit: %d samples, %d genes analyzed (assay: %s)\n",
              nrow(x$scores), x$n_genes_analyzed, x$assay_used))
  cat(sprintf("PC1 explains %.1f%% of variance.\n", x$percent_variance[1]))
  if (!is.null(x$batch_tests) && any(x$batch_tests$flagged, na.rm = TRUE)) {
    flagged_pcs <- paste(x$batch_tests$pc[x$batch_tests$flagged], collapse = ", ")
    cat(sprintf("Batch '%s' flagged as associated with PC(s): %s.\n",
                x$batch_col, flagged_pcs))
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

#' Plot PCA scores from a pca_audit result
#'
#' Produces a \code{ggplot2} scatter plot of two principal component score
#' axes, optionally colored by a sample-metadata column. This is a purely
#' visual convenience: the statistical results returned by
#' \code{\link{pca_audit}} do not depend on this function being called.
#'
#' @param pca_result An object of class \code{"rnaSentry_pca"} as returned by
#'   \code{\link{pca_audit}}.
#' @param pc_x,pc_y Integers selecting the two PC axes to plot. Defaults to
#'   \code{1} and \code{2}.
#' @param color_by Optional character. Name of a sample-metadata column in
#'   the original \code{colData(se)} to color points by.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(2)
#' counts <- matrix(rpois(240, lambda = 200), nrow = 20, ncol = 12,
#'                   dimnames = list(paste0("gene", 1:20), paste0("S", 1:12)))
#' coldata <- DataFrame(batch = factor(rep(c("B1", "B2"), each = 6)),
#'                       row.names = paste0("S", 1:12))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#' res <- pca_audit(se, batch_col = "batch")
#' plot_pca_audit(res)
#'
#' @export
plot_pca_audit <- function(pca_result, pc_x = 1, pc_y = 2, color_by = NULL) {
  if (!methods::is(pca_result, "rnaSentry_pca")) {
    stop("'pca_result' must be an rnaSentry_pca object from pca_audit().",
         call. = FALSE)
  }
  n_pcs <- ncol(pca_result$scores)
  if (!is.numeric(pc_x) || length(pc_x) != 1 || pc_x < 1 ||
      pc_x != round(pc_x) || pc_x > n_pcs) {
    stop(sprintf("'pc_x' must be an integer between 1 and %d.", n_pcs), call. = FALSE)
  }
  if (!is.numeric(pc_y) || length(pc_y) != 1 || pc_y < 1 ||
      pc_y != round(pc_y) || pc_y > n_pcs) {
    stop(sprintf("'pc_y' must be an integer between 1 and %d.", n_pcs), call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("plot_pca_audit() requires the 'ggplot2' package (in Suggests). ",
         "Install it to enable plotting.", call. = FALSE)
  }

  plot_data <- pca_result$scores
  x_var <- names(plot_data)[pc_x]
  y_var <- names(plot_data)[pc_y]
  if (!is.null(color_by)) {
    if (!color_by %in% colnames(pca_result$col_data)) {
      stop(sprintf("colData(se) has no column '%s'.", color_by), call. = FALSE)
    }
    plot_data[[color_by]] <- pca_result$col_data[[color_by]]
  }

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(.data[[x_var]], .data[[y_var]]))
  if (!is.null(color_by)) {
    p <- p + ggplot2::aes(color = .data[[color_by]])
  }
  p + ggplot2::geom_point() +
    ggplot2::labs(x = sprintf("%s (%.1f%%)", x_var,
                              pca_result$percent_variance[pc_x]),
                  y = sprintf("%s (%.1f%%)", y_var,
                              pca_result$percent_variance[pc_y]))
}
