#' rnaSentry: Guarded and Auditable Discovery of Prognostic RNA-Seq Signatures
#'
#' Orchestrates the full prognostic signature pipeline from bulk RNA-seq
#' \code{SummarizedExperiment} objects: sample-identity and quality-control
#' auditing, confounder and design-formula detection, signature construction
#' with cross-validation, survival modelling (Kaplan-Meier, Cox, parametric)
#' and external validation. Every automated decision is reported with its
#' supporting statistic and a plain-language rationale in a flag ledger and a
#' standalone HTML report via \code{generate_report()}.
#'
#' @section Main functions:
#' \itemize{
#'   \item \code{\link{load_counts}}: intake audit.
#'   \item \code{\link{qc_explore}}, \code{\link{sex_check}}: quality checks.
#'   \item \code{\link{design_audit}}, \code{\link{pca_audit}}: confounder detection.
#'   \item \code{\link{build_signature}}, \code{\link{km_curve}}, \code{\link{cox_model}}, \code{\link{survival_parametric}}, \code{\link{validate_external}}: signature and survival modelling.
#'   \item \code{\link{run_rnaSentry}}: end-to-end pipeline; \code{\link{lock_signature}}: single-use lock.
#' }
#'
#' @keywords internal
"_PACKAGE"
