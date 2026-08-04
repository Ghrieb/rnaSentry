# rnaSentry

Guarded and auditable discovery of prognostic RNA-seq signatures.

`rnaSentry` orchestrates the statistical pipeline that goes from a raw bulk
RNA-seq count matrix to a validated prognostic gene signature. Every
automated decision is reported with the statistic that justified it, and the
pipeline gates rather than silently degrades when its guardrails are not
met. The package is designed to make a biomarker discovery workflow
transparent, reproducible, and reviewable end to end.

## Installation

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rnaSentry")
```

## Quick start

```r
library(rnaSentry)
library(SummarizedExperiment)

# counts: a gene x sample integer matrix, gene symbols as rownames.
# coldata: a data.frame of sample metadata (one row per sample).
# load_counts() validates the intake and records an audit flag ledger.
se <- load_counts(counts, coldata)

# Run the discovery pipeline (time/event are survival columns in colData).
run <- run_rnaSentry(se, time_col = "time", event_col = "event",
                     design_vars = c("batch", "age"))

# Validate the signature against an independent cohort.
val <- validate_external(run$stages$build_signature, external_se,
                         cutpoint = run$stages$km_curve$cutpoint)
```

See `vignette("rnaSentry")` for a worked walkthrough and
`?run_rnaSentry`, `?build_signature`, `?validate_external` for details.

## Reproducibility guardrail

`run_rnaSentry()` is single-use per session: after it runs, the signature is
locked and a re-run is refused until the lock is released with
`lock_signature(sig, lock = FALSE)`. This prevents silent re-selection of
signature genes after survival analysis, which would otherwise inflate
reported concordance. Individual stages (for example `validate_external()`)
remain re-runnable on the locked signature.

## Reading the CV concordance

The concordance reported by `build_signature()` is a **screening-internal**
metric: genes are selected on the full cohort *before* the CV split, so the
CV fold is not fully untouched by selection. On a large candidate panel the
value is *optimistic* — even survival-permuted data yields a high null
concordance (≈0.8 on ~21k genes) rather than 0.5. Always interpret it
**relative to a matched null** (a permutation-null or random-gene control),
and treat `validate_external()` on an independent cohort as the only fully
out-of-sample estimate.

## Input requirements

- **Bulk RNA-seq** expression, with **gene symbols** as rownames.
- Either **raw integer counts** (log2-transformed internally) or an assay
  already named `logcounts`/`vst`.
- **Standard right-censored survival metadata**: a numeric follow-up time
  column and a 0/1 event indicator.
- Use `load_counts()` to validate the count matrix and metadata before the
  pipeline runs; it records any intake issues (duplicate samples,
  non-integer counts, all-zero gene rows, non-syntactic gene symbols) in the
  object's flag ledger.

## When to use rnaSentry

- You have a bulk RNA-seq cohort with survival follow-up and want to
  construct and validate a linear prognostic risk score.
- You want every selection, modeling, and validation decision recorded with
  its statistic so the analysis is auditable and reviewable.
- You want confounder and design-formula issues surfaced **before** you
  invest in differential expression or signature modeling.

## When NOT to use rnaSentry

- **Single-cell or spatial data.** The pipeline is designed for bulk
  matrices; there is no pseudobulk or cell-level modeling.
- **Non-survival end points.** The modeling stages assume standard,
  right-censored survival. Continuous, binary, or competing-risks outcomes
  are out of scope.
- **Data already log-transformed but stored as `counts`.** The pipeline
  detects this and warns (`possibly_log_scaled`), but you should store such
  data under the name `logcounts` instead.
- **Cohorts too small to support the signature size.** With fewer than
  roughly 5 events per signature gene, cross-validated concordance and
  hazard ratios are unstable; `build_signature()` warns
  (`events_per_parameter`). Consider fewer genes or a larger cohort, or a
  regularized approach.
- **Competing risks, time-varying covariates, or left truncation.** None of
  these are modeled.
- **Batch-effect correction or differential-expression analysis.** The
  pipeline *detects and reports* batch/confounder association; it does not
  correct or remove effects, and it makes no differential-expression calls.

## Limitations and caveats

- **Small-cohort trap.** The single most common misuse is underpowered
  discovery. `build_signature()` warns when events per signature gene fall
  below `min_events_per_parameter` (default 5, mid-range of the 5-9
  events-per-variable guidance; see `?build_signature` for the sources).
  Treat concordance and hazard ratios from such cohorts as indicative, not
  confirmatory.
- **Linear Cox risk scores only.** The signature is a weighted linear
  combination of expression; nonlinear or machine-learning score forms are
  not produced.
- **No correction step.** Confounders such as batch are flagged, and
  screening can be adjusted for design terms, but the pipeline never
  corrects or removes effects itself.
- **CPU-only and single-threaded.** Screening runs are sequential; memory,
  not parallelism, is the usual binding constraint on large cohorts.

## Documentation

- Vignette: `browseVignettes("rnaSentry")`
- Function reference: `help(package = "rnaSentry")`
- Validation dossier (sign-convention fix, parity tests, audit table,
  simulation, face validity, maintainer testing guide): see the
  `validation/` directory of the source package.

## Reporting issues

Please file bugs and feature requests at
<https://github.com/Ghrieb/rnaSentry/issues>.
