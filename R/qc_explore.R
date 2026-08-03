#' Quality-control audit of a bulk RNA-seq SummarizedExperiment
#'
#' Runs a battery of pre-analysis checks on a raw-counts
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}}: duplicate
#' sample identifiers, missing values in sample metadata, non-integer
#' count values, and library-size outliers detected via a robust
#' (median absolute deviation-based) rule. Every check is reported with
#' the statistic that triggered it, not just a pass/fail flag, so the
#' result can be inspected or logged rather than trusted blindly.
#'
#' @param se A \code{SummarizedExperiment} with a raw-count assay as its
#'   first assay (\code{assay(se, 1)}). Counts should be non-negative
#'   integers (or values coercible to integers without loss).
#' @param mad_threshold Numeric scalar. Samples whose log10 library size
#'   deviates from the median by more than \code{mad_threshold} robust
#'   median absolute deviations are flagged as library-size outliers.
#'   Defaults to \code{3}.
#'
#' @return An object of class \code{"rnaSentry_qc"} (a list) with elements:
#'   \describe{
#'     \item{n_samples}{Number of samples (columns) in \code{se}.}
#'     \item{n_genes}{Number of genes (rows) in \code{se}.}
#'     \item{duplicate_samples}{Character vector of duplicated sample IDs,
#'       or \code{character(0)} if none.}
#'     \item{non_integer_counts}{Logical; \code{TRUE} if the count assay
#'       contains non-integer values.}
#'     \item{missing_summary}{A data.frame reporting the proportion of
#'       missing values per \code{colData(se)} column.}
#'     \item{library_size_outliers}{Character vector of sample IDs flagged
#'       as library-size outliers.}
#'     \item{flags}{A data.frame summarizing every issue raised, with
#'       columns \code{check}, \code{severity}, and \code{detail}.}
#'   }
#'
#' @examples
#' library(SummarizedExperiment)
#' set.seed(1)
#' counts <- matrix(rpois(60, lambda = 200), nrow = 6, ncol = 10,
#'                   dimnames = list(paste0("gene", 1:6), paste0("S", 1:10)))
#' counts[, 1] <- counts[, 1] * 50L  # inject one library-size outlier
#' coldata <- DataFrame(condition = rep(c("A", "B"), each = 5),
#'                       row.names = paste0("S", 1:10))
#' se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
#'
#' qc <- qc_explore(se)
#' qc$flags
#'
#' @export
qc_explore <- function(se, mad_threshold = 3) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }

  counts <- SummarizedExperiment::assay(se, 1)
  sample_ids <- colnames(se)
  cdata <- as.data.frame(SummarizedExperiment::colData(se))

  flags <- .new_flags("qc_explore")

  # 1. duplicate sample IDs
  dupes <- sample_ids[duplicated(sample_ids)]
  if (length(dupes) > 0) {
    flags <- .add_flag(flags, "duplicate_samples", "critical",
                       paste("Duplicated sample IDs:", paste(unique(dupes), collapse = ", ")))
  }

  # 2. non-integer counts
  non_integer <- any(abs(counts - round(counts)) > 1e-8, na.rm = TRUE)
  if (non_integer) {
    flags <- .add_flag(flags, "non_integer_counts", "warning",
                       "Count assay contains non-integer values; verify this is raw count data.")
  }

  # 3. missing values in colData
  missing_summary <- data.frame(
    column = colnames(cdata),
    pct_missing = vapply(cdata, function(x) mean(is.na(x)) * 100, numeric(1)),
    row.names = NULL, stringsAsFactors = FALSE
  )
  bad_cols <- missing_summary[missing_summary$pct_missing > 0, ]
  if (nrow(bad_cols) > 0) {
    for (i in seq_len(nrow(bad_cols))) {
      flags <- .add_flag(flags, "missing_metadata", "warning",
                         sprintf("Column '%s' has %.1f%% missing values.",
                                 bad_cols$column[i], bad_cols$pct_missing[i]))
    }
  }

  # 4. library-size outliers (robust, MAD-based, on log10 scale)
  lib_sizes <- colSums(counts, na.rm = TRUE)
  log_lib <- log10(lib_sizes + 1)
  med <- stats::median(log_lib)
  mad_val <- stats::mad(log_lib)
  outlier_idx <- if (mad_val == 0) logical(length(log_lib)) else
    abs(log_lib - med) / mad_val > mad_threshold
  outliers <- sample_ids[outlier_idx]
  if (length(outliers) > 0) {
    flags <- .add_flag(flags, "library_size_outlier", "warning",
                       paste("Samples with library size >", mad_threshold,
                             "MADs from the median:", paste(outliers, collapse = ", ")))
  }

  result <- list(
    n_samples = ncol(se),
    n_genes = nrow(se),
    duplicate_samples = unique(dupes),
    non_integer_counts = non_integer,
    missing_summary = missing_summary,
    library_size_outliers = outliers,
    flags = flags
  )
  class(result) <- "rnaSentry_qc"
  result
}

#' @export
print.rnaSentry_qc <- function(x, ...) {
  cat(sprintf("rnaSentry QC audit: %d samples x %d genes\n", x$n_samples, x$n_genes))
  if (nrow(x$flags) == 0) {
    cat("No issues flagged.\n")
  } else {
    cat(sprintf("%d issue(s) flagged:\n", nrow(x$flags)))
    for (i in seq_len(nrow(x$flags))) {
      cat(sprintf("  [%s] %s: %s\n", x$flags$severity[i], x$flags$check[i], x$flags$detail[i]))
    }
  }
  invisible(x)
}
