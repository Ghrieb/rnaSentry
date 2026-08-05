# Package index

## End-to-end orchestration

Run the full guarded discovery pipeline from a SummarizedExperiment.

- [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  : Run the complete rnaSentry pipeline

## Signature construction and locking

Build a stable Cox-driven gene signature with explicit cross-validation
stability reporting, then lock it for export.

- [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  : Build a prognostic gene signature from a survival cohort
- [`lock_signature()`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md)
  : Lock or unlock an rnaSentry signature

## Auditing stages (gated, not silent)

Sample identity, quality control, confounding and batch design
detection, principal-component audit, and the documented XIST-based sex
check. Each stage reports its supporting statistic and a plain-language
rationale.

- [`qc_explore()`](https://ghrieb.github.io/rnaSentry/reference/qc_explore.md)
  : Quality-control audit of a bulk RNA-seq SummarizedExperiment
- [`sex_check()`](https://ghrieb.github.io/rnaSentry/reference/sex_check.md)
  : Verify reported sex against inferred biological sex from expression
- [`pca_audit()`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
  : PCA audit of expression data with optional batch association testing
- [`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
  : Audit candidate design covariates before differential-expression
  modeling
- [`plot_pca_audit()`](https://ghrieb.github.io/rnaSentry/reference/plot_pca_audit.md)
  : Plot PCA scores from a pca_audit result

## Survival modeling

Kaplan-Meier curves, log-rank tests, Cox proportional hazards, and
parametric alternatives.

- [`km_curve()`](https://ghrieb.github.io/rnaSentry/reference/km_curve.md)
  : Kaplan-Meier survival analysis of a signature's risk groups
- [`cox_model()`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md)
  : Adjusted Cox proportional-hazards model for a signature
- [`survival_parametric()`](https://ghrieb.github.io/rnaSentry/reference/survival_parametric.md)
  : Parametric survival models for a signature's risk score

## External validation and reporting

Validate a locked signature against public cohorts and render the
audited report.

- [`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
  : Validate a signature on an external cohort
- [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
  : Validate and construct a SummarizedExperiment from count data
- [`generate_report()`](https://ghrieb.github.io/rnaSentry/reference/generate_report.md)
  : Generate the rnaSentry analysis report
