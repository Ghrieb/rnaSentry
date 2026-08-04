#' Validate and construct a SummarizedExperiment from count data
#'
#' Thin validation wrapper around
#' \code{SummarizedExperiment::SummarizedExperiment()} that catches the most
#' common input mistakes before they propagate into the pipeline: duplicate or
#' missing sample IDs, non-integer counts, empty rows or columns, and
#' non-syntactic gene identifiers.
#'
#' @param counts Integer or numeric matrix of gene expression values (genes in
#'   rows, samples in columns). Row names must be gene identifiers; column names
#'   must be sample identifiers.
#' @param colData A \code{data.frame} or \code{S4Vectors::DataFrame} of
#'   sample-level metadata. Row names must match the column names of
#'   \code{counts}. When omitted, an empty \code{DataFrame} is created.
#' @param assay_name Character. Name for the assay slot. Defaults to
#'   \code{"counts"}.
#' @param check_gene_ids Logical. When \code{TRUE} (default) the function
#'   warns if gene identifiers look like Ensembl accessions (starting with
#'   \code{"ENS"}) or contain non-syntactic characters, because downstream
#'   stages rely on gene-symbol rownames.
#'
#' @return A \code{SummarizedExperiment} with one assay named
#'   \code{assay_name}. Any issues found are recorded in a flag ledger
#'   accessible via \code{S4Vectors::metadata(se)$flags} and surfaced via
#'   \code{warning()}.
#'
#' @section Assumptions and limitations:
#' This function validates structure; it does not normalise, filter, or
#' transform the count matrix. Gene identifiers are checked but not converted.
#' Sample-metadata columns used later by the pipeline (e.g. \code{time_col},
#' \code{event_col}, \code{batch_col}) are not validated here — that happens
#' in the stage that first uses them.
#'
#' @examples
#' counts <- matrix(1:12, nrow = 3,
#'                  dimnames = list(c("TP53", "BRCA1", "MYC"),
#'                                  c("S1", "S2", "S3", "S4")))
#' coldata <- data.frame(
#'   time = c(10, 12, 8, 15),
#'   event = c(1, 0, 1, 0),
#'   row.names = colnames(counts))
#' se <- load_counts(counts, coldata)
#' se
#'
#' @export
load_counts <- function(counts, colData, assay_name = "counts",
                        check_gene_ids = TRUE) {
  flags <- .new_flags("load_counts")

  # --- counts matrix --------------------------------------------------------
  if (!is.matrix(counts)) {
    stop("'counts' must be a matrix.", call. = FALSE)
  }
  if (!is.numeric(counts)) {
    stop("'counts' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("'counts' has no rows or no columns.", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("'counts' must have both rownames (gene IDs) and column names ",
         "(sample IDs).", call. = FALSE)
  }

  # Duplicate sample IDs
  dup_samples <- duplicated(colnames(counts))
  if (any(dup_samples)) {
    dups <- unique(colnames(counts)[dup_samples])
    flags <- .add_flag(flags, "duplicate_samples", "critical",
                       sprintf("Duplicate sample IDs: %s.",
                               paste(dups, collapse = ", ")))
    warning(sprintf("Duplicate sample IDs found: %s. Removing duplicates.",
                    paste(dups, collapse = ", ")), call. = FALSE)
    counts <- counts[, !dup_samples, drop = FALSE]
  }

  # Non-integer counts
  if (any(abs(counts - round(counts)) > 1e-8, na.rm = TRUE)) {
    flags <- .add_flag(flags, "non_integer_counts", "warning",
                       "Count matrix contains non-integer values.")
    warning("Count matrix contains non-integer values. ",
            "Downstream stages expect integer counts.", call. = FALSE)
  }

  # Empty rows (all zero or all NA)
  row_sums <- rowSums(counts, na.rm = TRUE)
  empty_rows <- is.na(row_sums) | row_sums == 0
  if (any(empty_rows)) {
    n_empty <- sum(empty_rows)
    flags <- .add_flag(flags, "empty_rows", "warning",
                       sprintf("%d gene(s) have all-zero or all-NA counts.",
                               n_empty))
    warning(sprintf("%d gene(s) have all-zero or all-NA counts and will be ",
                    n_empty), "retained but may be filtered later.",
            call. = FALSE)
  }

  # --- colData --------------------------------------------------------------
  if (missing(colData) || is.null(colData)) {
    colData <- data.frame(row.names = colnames(counts))
  }
  if (!methods::is(colData, "DataFrame") && !is.data.frame(colData)) {
    stop("'colData' must be a data.frame or DataFrame.", call. = FALSE)
  }
  if (is.null(rownames(colData))) {
    stop("'colData' must have row names matching the column names of ",
         "'counts'.", call. = FALSE)
  }
  missing_ids <- setdiff(colnames(counts), rownames(colData))
  if (length(missing_ids) > 0) {
    stop(sprintf("colData is missing row names for sample(s): %s.",
                 paste(missing_ids[1:min(5, length(missing_ids))],
                       collapse = ", ")),
         call. = FALSE)
  }
  # Reorder colData to match counts columns
  colData <- colData[colnames(counts), , drop = FALSE]
  if (!methods::is(colData, "DataFrame")) {
    colData <- S4Vectors::DataFrame(colData)
  }

  # --- gene ID format check -------------------------------------------------
  if (check_gene_ids) {
    ens_pattern <- "^ENS[A-Z]*[0-9]"
    ens_genes <- grepl(ens_pattern, rownames(counts))
    if (mean(ens_genes, na.rm = TRUE) > 0.5) {
      flags <- .add_flag(flags, "ensembl_ids_detected", "info",
                         paste0("More than half of gene IDs look like ",
                                "Ensembl accessions; downstream stages ",
                                "expect gene symbols."))
      warning("More than half of gene identifiers look like Ensembl accessions ",
              "(e.g. 'ENSG00000...'). The pipeline expects gene symbols. ",
              "Consider mapping IDs before proceeding.", call. = FALSE)
    }
    syntactic <- make.names(rownames(counts)) == rownames(counts)
    if (mean(syntactic, na.rm = TRUE) < 0.9) {
      n_nonsyntactic <- sum(!syntactic)
      flags <- .add_flag(flags, "non_syntactic_ids", "info",
                         sprintf("%d gene ID(s) contain non-syntactic characters.",
                                 n_nonsyntactic))
      warning(sprintf("%d of %d gene identifiers contain non-syntactic ",
                      n_nonsyntactic, nrow(counts)),
              "characters. They will work but may appear mangled in plots.",
              call. = FALSE)
    }
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(as.matrix(counts)),
    colData = colData
  )
  SummarizedExperiment::assayNames(se) <- assay_name

  S4Vectors::metadata(se)$flags <- flags
  se
}
