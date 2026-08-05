# Validate a signature on an external cohort

`validate_external()` applies a locked-down signature to an independent
cohort: it scores every external sample with the signature's joint Cox
coefficients (the gene-only model, so the score is portable to cohorts
that do not record the discovery cohort's design covariates), splits the
cohort into `"low"` and `"high"` risk groups at the supplied *discovery*
cutpoint, and reports the external discriminant validity: a log-rank
test between the groups, per-group median survival and event counts, and
the concordance index of the continuous risk score against the external
survival outcome (Harrell's C in the risk-score convention: a higher
risk score is concordant with an earlier event).

## Usage

``` r
validate_external(
  sig,
  external_se,
  time_col = NULL,
  event_col = NULL,
  cutpoint,
  min_gene_overlap_frac = 1,
  drop_missing = FALSE
)
```

## Arguments

- sig:

  An object of class `"rnaSentry_signature"` as returned by
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md).

- external_se:

  A `SummarizedExperiment` of an independent cohort, with an expression
  assay and survival metadata.

- time_col, event_col:

  Names of the survival columns in `colData(external_se)`. When `NULL`
  (default) the column names recorded in `sig` are used.

- cutpoint:

  The risk-score threshold defining the `"high"` group (samples with
  score `>= cutpoint`). Must be a single finite number; it is never
  derived from the external cohort.

- min_gene_overlap_frac:

  Numeric in (0, 1\]. The minimum fraction of signature genes that must
  be present in the external cohort when `drop_missing = TRUE`. Defaults
  to `1`.

- drop_missing:

  Logical. When `TRUE`, signature genes missing from the external cohort
  are dropped (and flagged by name), provided the remaining overlap is
  at least `min_gene_overlap_frac`. When `FALSE` (default) any missing
  gene is an error.

## Value

An object of class `"rnaSentry_external"` (a list) with elements:

- genes_used:

  The signature genes present in the external cohort.

- genes_missing:

  The signature genes absent from the external cohort (empty when none).

- score:

  Per-sample external risk score.

- groups:

  Factor with levels `"low"` and `"high"` from the discovery cutpoint.

- cutpoint, cutpoint_type:

  The applied cutpoint and `"custom"`.

- km_fit:

  The `survfit` object with one curve per group.

- log_rank_p:

  p-value from the log-rank test between groups.

- median_survival:

  Named vector of median survival per group.

- events_low, events_high:

  Event counts per group.

- concordance:

  Concordance index of the continuous risk score against the external
  survival outcome.

- sig_locked:

  Whether the input signature was locked.

- created:

  Provenance timestamp.

- flags:

  Data.frame of issues raised (dropped genes, sparse events, non-finite
  concordance), with columns `check`, `severity`, `detail` and `stage`.

## Details

The stage *never* recomputes a within-cohort cutpoint. The `cutpoint`
must be passed explicitly (for example the discovery-cohort median
recorded in the `cutpoint` element returned by
[`km_curve`](https://ghrieb.github.io/rnaSentry/reference/km_curve.md));
passing no cutpoint is an error, so external risk groups are always
defined by the same threshold that defined the discovery groups.

Missing genes are handled strictly by default: `drop_missing = FALSE`
and `min_gene_overlap_frac = 1` (the defaults) require every signature
gene to be present, and error otherwise, naming the missing genes. Set
`drop_missing = TRUE` (with `min_gene_overlap_frac < 1` as needed) to
explicitly allow scoring on a subset of genes; the dropped genes are
named in a flag.

The stage is read-only: it works on both locked and unlocked signatures
and records the signature's lock state in its output.

## Assumptions and limitations

External risk groups are defined exclusively by the supplied discovery
cutpoint; the stage never re-estimates a threshold from the external
cohort. The cohort is assumed to be bulk RNA-seq with standard
right-censored survival and the same assay convention as the discovery
data (raw counts or a correctly-named `logcounts` assay; pre-scaled data
in a `counts`-named slot is flagged as `possibly_log_scaled`). Competing
risks and other non-standard survival structures are not modeled. Unlike
the discovery stage
([`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md),
which requires at least two events), a one-event minimum applies here;
the resulting instability is reported through the `sparse_events`
guardrail.

## Examples

``` r
library(SummarizedExperiment)
set.seed(11)
counts <- matrix(rpois(800, lambda = 500), nrow = 40, ncol = 20,
                 dimnames = list(paste0("gene", 1:40), paste0("S", 1:20)))
sig_expr <- colMeans(counts[1:5, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(20, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                       folds = 2, seed = 1)
#> Warning: Only 12 event(s) for a 5-gene signature (2.4 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
cutoff <- km_curve(sig, se)$cutpoint
ext <- se
ext <- ext[, sample(ncol(ext))]
val <- validate_external(sig, ext, cutpoint = cutoff)
val
#> rnaSentry external validation: 5/5 signature gene(s) used, 20 samples.
#> Risk groups: low n = 10 (5 events), high n = 10 (7 events).
#> Cutpoint -94.3 (custom, from discovery).
#> Log-rank p = 0.009113. Concordance = 0.791.
#> Signature not locked.
```
