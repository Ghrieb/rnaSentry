# Run the complete rnaSentry pipeline

`run_rnaSentry()` executes the full rnaSentry pipeline on a single
cohort and collects every stage result into one object:

1.  [`design_audit`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
    over `design_vars` (when any are supplied).

2.  [`pca_audit`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
    (with `batch_col` when supplied).

3.  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md),
    then
    [`lock_signature`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md)
    so the survival-modeling stages operate on the published signature.

4.  [`km_curve`](https://ghrieb.github.io/rnaSentry/reference/km_curve.md),
    [`cox_model`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md)
    and
    [`survival_parametric`](https://ghrieb.github.io/rnaSentry/reference/survival_parametric.md)
    on the locked signature.

5.  [`generate_report`](https://ghrieb.github.io/rnaSentry/reference/generate_report.md)
    assembling all stages into an HTML report (when
    `render_report = TRUE`).

## Usage

``` r
run_rnaSentry(
  se,
  time_col,
  event_col,
  outcome_col = "overall_survival",
  design_terms = character(0),
  design_vars = design_terms,
  batch_col = NULL,
  method = c("top_n", "p_value"),
  top_n = 20,
  p_threshold = 0.05,
  repeats = 5,
  folds = 5,
  seed = NULL,
  adjust_for_design = TRUE,
  min_events_per_parameter = 5,
  BPPARAM = NULL,
  report_file = "rnaSentry_report.html",
  report_dir = ".",
  render_report = TRUE
)
```

## Arguments

- se:

  A `SummarizedExperiment` with a count or normalized expression assay
  and survival metadata.

- time_col:

  Character. Column of `colData(se)` with follow-up time.

- event_col:

  Character. Column of `colData(se)` with the event indicator (0/1 or
  logical).

- outcome_col:

  Character. Display label for the outcome. Defaults to
  `"overall_survival"`.

- design_terms:

  Character vector of `colData(se)` columns to adjust the signature's
  univariate screening by. Defaults to `character(0)`.

- design_vars:

  Character vector of `colData(se)` columns scanned by the
  [`design_audit`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
  confounder scan. Defaults to `design_terms`.

- batch_col:

  Optional character. `colData(se)` column tested by the
  [`pca_audit`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
  batch scan. Defaults to `NULL`.

- method, top_n, p_threshold, repeats, folds, seed, adjust_for_design,
  min_events_per_parameter, BPPARAM:

  Passed to
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md).

- report_file:

  Character. Output file name for the HTML report.

- report_dir:

  Character. Directory (which must exist) to write the report into.

- render_report:

  Logical. When `TRUE` (default) the report is rendered as the final
  step.

## Value

An object of class `"rnaSentry_run"` (a list) with elements:

- stages:

  Named list of the stage results, with entries `design_audit` (only
  when `design_vars` was supplied), `pca_audit`, `build_signature`,
  `km_curve`, `cox_model` and `survival_parametric`.

- report:

  Path of the rendered HTML report, or `NULL` when
  `render_report = FALSE`.

## Details

The result's `stages` list can be passed directly to
[`generate_report`](https://ghrieb.github.io/rnaSentry/reference/generate_report.md),
or any stage can be re-run individually (for example
[`validate_external`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
against an independent cohort, using the cutpoint recorded by the
`km_curve` stage).

## Reproducibility lock

The run locks the discovered signature (see
[`lock_signature`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md)),
recording a session-level fingerprint. While any locked signature
exists, `run_rnaSentry()` refuses to run again, preventing silent
re-selection of signature genes after survival analysis. To run the
pipeline a second time in the same session (for example on another
cohort), release the lock first with `lock_signature(sig, lock = FALSE)`
on the locked signature. Individual stages remain re-runnable while the
lock is set.

## Assumptions and limitations

The pipeline targets bulk RNA-seq with standard right-censored survival
and linear Cox risk scores (see the Assumptions sections of
[`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
and
[`cox_model`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md)).
It detects and reports batch and design confounders but does not correct
for them, and it makes no differential-expression calls. On small
cohorts an `events_per_parameter` warning from
[`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
is expected and should be interpreted as a call for caution, not
ignored.

## Examples

``` r
library(SummarizedExperiment)
set.seed(8)
counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
                  dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
sig_expr <- colMeans(counts[1:5, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(20, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 batch = factor(rep(c("B1", "B2"), 10)),
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

run <- run_rnaSentry(se, "time", "event", design_vars = "batch",
                      top_n = 5, repeats = 1, folds = 2, seed = 1,
                      report_dir = tempdir())
#> Warning: Only 13 event(s) for a 5-gene signature (2.6 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
run
#> rnaSentry pipeline run.
#>   - design_audit
#>   - pca_audit
#>   - build_signature
#>   - km_curve
#>   - cox_model
#>   - survival_parametric
#> Report: /tmp/Rtmp4OQvut/rnaSentry_report.html
```
