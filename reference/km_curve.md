# Kaplan-Meier survival analysis of a signature's risk groups

Scores each sample of a `SummarizedExperiment` with the signature's
joint Cox coefficients and splits the cohort into `"low"` and `"high"`
risk groups. By default the split is at the median risk score of the
supplied cohort; pass `cutpoint` to re-apply a cutpoint established
elsewhere (for example the discovery-cohort median recorded in a
previous `km_curve()` call), so risk groups are comparable across
cohorts. Fits a Kaplan-Meier survival curve for each group with a
log-rank test. The analysis is read-only: it works on both locked and
unlocked signatures and records the signature's lock state in its
output.

## Usage

``` r
km_curve(sig, se, cutpoint = NULL)
```

## Arguments

- sig:

  An object of class `"rnaSentry_signature"` as returned by
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md).

- se:

  A `SummarizedExperiment` containing the signature genes in its rows
  and the survival metadata columns recorded in `sig` (`sig$time_col`
  and `sig$event_col`). Expression is taken from the same analysis assay
  used at signature-build time (see
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)).

- cutpoint:

  Numeric risk-score threshold: samples with score `>= cutpoint` are
  `"high"` and the rest `"low"`. When `NULL` (default) the split is at
  the median risk score of the supplied cohort.

## Value

An object of class `"rnaSentry_km"` (a list) with elements:

- genes:

  The signature genes scored.

- score:

  Per-sample risk score (linear predictor).

- groups:

  Factor with levels `"low"` and `"high"` from the risk-score split.

- cutpoint:

  The risk-score cutpoint applied.

- cutpoint_type:

  `"median"` when the cohort median was used, `"custom"` when `cutpoint`
  was supplied.

- fit:

  The `survfit` object with one curve per group.

- log_rank_p:

  p-value from the log-rank test.

- median_survival:

  Named vector of median survival per group.

- events_low, events_high:

  Event counts per group.

- sig_locked:

  Whether the input signature was locked.

- created:

  Provenance timestamp.

- flags:

  Data.frame of issues raised (e.g. sparse events per group), with
  columns `check`, `severity`, `detail` and `stage`.

## Examples

``` r
library(SummarizedExperiment)
set.seed(5)
counts <- matrix(rpois(120, lambda = 500), nrow = 12, ncol = 10,
                 dimnames = list(paste0("gene", 1:12), paste0("S", 1:10)))
sig_expr <- colMeans(counts[1:3, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(10, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(10, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
sig <- build_signature(se, time_col = "time", event_col = "event",
                       top_n = 3, repeats = 1, folds = 2, seed = 1)
#> Warning: Only 7 event(s) for a 3-gene signature (2.3 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
km <- km_curve(sig, se)
km
#> rnaSentry KM analysis: 3-gene signature, 10 samples.
#> Risk groups: low n = 5 (2 events), high n = 5 (5 events).
#> Cutpoint 264 (median).
#> Median survival: low NA, high 6.89.
#> Log-rank p = 0.20596.
#> Signature not locked.
#> 1 issue(s) flagged:
#>   [warning] sparse_events: Risk group event counts are low (low: 2, high: 5); median survival estimates may be unstable.
```
