#' Verify reported sex against inferred biological sex from expression
#'
#' Compares the sex recorded in sample metadata against sex inferred from
#' expression of \emph{XIST} and a panel of Y-chromosome marker genes.
#' Designed to catch sample swaps and metadata entry errors before they
#' propagate into downstream differential expression or survival analysis.
#'
#' Inference rule: a sample is called \code{"F"} when XIST expression is
#' high relative to the Y-gene panel, \code{"M"} when the reverse holds,
#' and \code{"ambiguous"} when neither signal is clearly dominant (for
#' example, low-quality or highly degraded samples). The decision uses
#' the sign of each sample's XIST-minus-Y-panel rank score: a score at
#' least one rank unit above zero is called \code{"F"}, at least one rank
#' unit below zero \code{"M"}, and a score within one rank unit of zero
#' \code{"ambiguous"}. Because the rule is based on ranks rather than
#' hard-coded expression cutoffs, it is robust across different
#' normalizations and sequencing depths.
#'
#' @param se A \code{SummarizedExperiment} with gene symbols as
#'   \code{rownames(se)} and a numeric (counts or normalized) assay as its
#'   first assay.
#' @param sex_col Character. Name of the \code{colData(se)} column holding
#'   reported sex, coded as \code{"M"} / \code{"F"}. Defaults to
#'   \code{"sex"}.
#' @param xist_gene Character. Gene symbol for XIST. Defaults to
#'   \code{"XIST"}.
#' @param y_genes Character vector of Y-chromosome marker gene symbols.
#'   Defaults to \code{c("RPS4Y1", "DDX3Y", "KDM5D")}.
#'
#' @return A data.frame with one row per sample and columns
#'   \code{sample_id}, \code{reported_sex}, \code{inferred_sex}, and
#'   \code{status} (one of \code{"OK"}, \code{"MISMATCH"},
#'   \code{"AMBIGUOUS"}, or \code{"MISSING"} when reported sex is absent).
#'
#' @examples
#' library(SummarizedExperiment)
#' genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "GAPDH")
#' samples <- paste0("S", 1:6)
#' expr <- matrix(0, nrow = length(genes), ncol = length(samples),
#'                 dimnames = list(genes, samples))
#' # samples 1-3 look female (high XIST, low Y genes)
#' expr["XIST", 1:3] <- c(950, 900, 875)
#' expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 1:3] <- 5
#' # samples 4-6 look male, but S6's metadata will be mislabeled below
#' expr["XIST", 4:6] <- 5
#' expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 4:6] <- 400
#' expr["GAPDH", ] <- 1000
#'
#' coldata <- DataFrame(sex = c("F", "F", "F", "M", "M", "F"),
#'                       row.names = samples)  # S6 mislabeled as F
#' se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)
#'
#' sex_check(se)
#'
#' @export
sex_check <- function(se, sex_col = "sex", xist_gene = "XIST",
                       y_genes = c("RPS4Y1", "DDX3Y", "KDM5D")) {
  if (!methods::is(se, "SummarizedExperiment")) {
    stop("'se' must be a SummarizedExperiment object.", call. = FALSE)
  }
  if (!is.character(sex_col) || length(sex_col) != 1 ||
      is.na(sex_col) || !nzchar(sex_col)) {
    stop("'sex_col' must be a single non-empty character string.", call. = FALSE)
  }
  if (!sex_col %in% colnames(SummarizedExperiment::colData(se))) {
    stop(sprintf("colData(se) has no column '%s'.", sex_col), call. = FALSE)
  }

  present_y <- intersect(y_genes, rownames(se))
  if (length(present_y) == 0 || !(xist_gene %in% rownames(se))) {
    stop("None of the required XIST/Y-gene markers were found in rownames(se). ",
         "Check that gene symbols match your reference annotation.", call. = FALSE)
  }
  if (length(present_y) < length(y_genes)) {
    warning(sprintf("Only %d of %d requested Y-gene markers found; proceeding with those.",
                     length(present_y), length(y_genes)))
  }

  expr <- SummarizedExperiment::assay(se, 1)
  xist_expr <- as.numeric(expr[xist_gene, ])
  y_expr <- if (length(present_y) == 1) as.numeric(expr[present_y, ]) else
    colMeans(expr[present_y, , drop = FALSE])

  # rank-based score: robust to scale/normalization differences
  score <- rank(xist_expr) - rank(y_expr)

  # A score within one rank unit of zero means XIST and the Y panel have
  # comparable ranks: neither signal is clearly dominant.
  inferred <- ifelse(score >= 1, "F", ifelse(score <= -1, "M", "ambiguous"))

  reported <- as.character(SummarizedExperiment::colData(se)[[sex_col]])
  status <- ifelse(inferred == "ambiguous", "AMBIGUOUS",
                   ifelse(is.na(reported), "MISSING",
                          ifelse(inferred == reported, "OK", "MISMATCH")))

  data.frame(
    sample_id = colnames(se),
    reported_sex = reported,
    inferred_sex = inferred,
    status = status,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
