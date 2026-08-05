# rnaSentry

Guarded and auditable discovery of prognostic RNA-seq signatures.

`rnaSentry` orchestrates the statistical pipeline that goes from a raw
bulk RNA-seq count matrix to a validated prognostic gene signature.
Every automated decision is reported with the statistic that justified
it, and the pipeline gates rather than silently degrades when its
guardrails are not met. The package is designed to make a biomarker
discovery workflow transparent, reproducible, and reviewable end to end.

## Status

This package and repository are the working artifacts of an in-progress
manuscript and Bioconductor submission. Results are **illustrative and
provisional** — not peer-reviewed — and may change before a preprint is
published. Treat all case-study numbers as demonstrations, not clinical
claims.

## Installation

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rnaSentry")
```

## Quick start

``` r

library(rnaSentry)
library(SummarizedExperiment)

# counts: a gene x sample integer matrix, gene symbols as rownames.
# coldata: a data.frame of sample metadata (one row per sample).
# load_counts() validates the intake and records an audit flag ledger.
se <- load_counts(counts, coldata)

# Run the discovery pipeline (time/event are survival columns in colData).
run <- run_rnaSentry(se, time_col = "time", event_col = "event",
                     design_vars = c("batch", "age"))

# Optional: parallelize the repeated cross-validation folds.
# run <- run_rnaSentry(se, time_col = "time", event_col = "event",
#                      design_vars = c("batch", "age"),
#                      BPPARAM = BiocParallel::SnowParam(2))


# Validate the signature against an independent cohort.
val <- validate_external(run$stages$build_signature, external_se,
                         cutpoint = run$stages$km_curve$cutpoint)
```

See
[`vignette("getting-started")`](https://ghrieb.github.io/rnaSentry/articles/getting-started.md)
to bridge your own CSV files into the pipeline,
[`vignette("rnaSentry")`](https://ghrieb.github.io/rnaSentry/articles/rnaSentry.md)
for a worked walkthrough, and
[`?run_rnaSentry`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md),
[`?build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md),
[`?validate_external`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
for details.

## Getting started with your own data (CSV -\> SummarizedExperiment)

Most users start with two plain files: a **counts matrix** (genes in
rows, samples in columns) and a **clinical table** (one row per sample).
rnaSentry accepts a `SummarizedExperiment`, so the only bridge you need
is the small block below. It reads both files, keeps exactly the samples
present in both (mismatched sample IDs are the \#1 pitfall), and
validates the result with
[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md),
which records any intake issues in a flag ledger on the object.

``` r

library(rnaSentry)
library(SummarizedExperiment)

# 1. Read the two CSV files.
counts   <- as.matrix(read.csv("counts.csv",   row.names = 1))
clinical <- read.csv("clinical.csv", row.names = 1)

# 2. Keep samples present in BOTH files (the #1 pitfall: mismatched IDs).
common   <- intersect(colnames(counts), rownames(clinical))
counts   <- counts[,  common, drop = FALSE]
clinical <- clinical[common, ,   drop = FALSE]
# check: dim(counts)[2] == nrow(clinical)

# 3. Build + validate the SummarizedExperiment (issues -> flag ledger).
se <- load_counts(counts, clinical)
S4Vectors::metadata(se)$flags

# 4. Audit candidate confounders before any modeling.
da <- design_audit(se, design_vars = c("batch", "stage"))
da$formula_text

# 5. Run the full pipeline (needs 'time' + 'event' columns in clinical.csv).
# run <- run_rnaSentry(se, time_col = "time", event_col = "event",
#                      design_vars = c("batch", "stage"))
```

Notes on your input files:

- **Counts go in a matrix, not a data.frame.**
  [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
  expects a matrix, so wrap the
  [`read.csv()`](https://rdrr.io/r/utils/read.table.html) result in
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html). Column names must
  be sample IDs; row names should be gene symbols (if your file uses
  Ensembl accessions, map them to symbols first – the pipeline warns
  otherwise).
- **Sample IDs must match between the two files.**
  [`intersect()`](https://generics.r-lib.org/reference/setops.html)
  keeps only samples found in both; samples missing from either file are
  silently dropped. Check the dimensions after intersecting – if you
  lose hundreds of samples you may have a versioning or delimiter
  mismatch.
- **Raw counts are log2-transformed internally.**
  [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
  stores your matrix under the assay name `counts`, and the pipeline
  computes `log2(counts + 1)` when it runs. If your file is *already
  normalized* (e.g. CPM, TPM, or VST), name that assay `logcounts`
  instead – either pass `load_counts(..., assay_name = "logcounts")` or
  rename after construction.
- **Survival columns.** The discovery pipeline needs a numeric follow-up
  `time` column and a `0/1` `event` indicator in the clinical table.
  [`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
  and the other auditing stages do not need them.

## Reproducibility guardrail

[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
is single-use per session: after it runs, the signature is locked and a
re-run is refused until the lock is released with
`lock_signature(sig, lock = FALSE)`. This prevents silent re-selection
of signature genes after survival analysis, which would otherwise
inflate reported concordance. Individual stages (for example
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md))
remain re-runnable on the locked signature.

## Reading the CV concordance

The concordance reported by
[`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
is a **screening-internal** metric: genes are selected on the full
cohort *before* the CV split, so the CV fold is not fully untouched by
selection. On a large candidate panel the value is *optimistic* — even
survival-permuted data yields a high null concordance (≈0.8 on ~21k
genes) rather than 0.5. Always interpret it **relative to a matched
null** (a permutation-null or random-gene control), and treat
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
on an independent cohort as the only fully out-of-sample estimate.

## Input requirements

- **Bulk RNA-seq** expression, with **gene symbols** as rownames.
- Either **raw integer counts** (log2-transformed internally) or an
  assay already named `logcounts`/`vst`.
- **Standard right-censored survival metadata**: a numeric follow-up
  time column and a 0/1 event indicator.
- Use
  [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
  to validate the count matrix and metadata before the pipeline runs; it
  records any intake issues (duplicate samples, non-integer counts,
  all-zero gene rows, non-syntactic gene symbols) in the object’s flag
  ledger.

## When to use rnaSentry

- You have a bulk RNA-seq cohort with survival follow-up and want to
  construct and validate a linear prognostic risk score.
- You want every selection, modeling, and validation decision recorded
  with its statistic so the analysis is auditable and reviewable.
- You want confounder and design-formula issues surfaced **before** you
  invest in differential expression or signature modeling.

## When NOT to use rnaSentry

- **Single-cell or spatial data.** The pipeline is designed for bulk
  matrices; there is no pseudobulk or cell-level modeling.
- **Non-survival end points.** The modeling stages assume standard,
  right-censored survival. Continuous, binary, or competing-risks
  outcomes are out of scope.
- **Data already log-transformed but stored as `counts`.** The pipeline
  detects this and warns (`possibly_log_scaled`), but you should store
  such data under the name `logcounts` instead.
- **Cohorts too small to support the signature size.** With fewer than
  roughly 5 events per signature gene, cross-validated concordance and
  hazard ratios are unstable;
  [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  warns (`events_per_parameter`). Consider fewer genes or a larger
  cohort, or a regularized approach.
- **Competing risks, time-varying covariates, or left truncation.** None
  of these are modeled.
- **Batch-effect correction or differential-expression analysis.** The
  pipeline *detects and reports* batch/confounder association; it does
  not correct or remove effects, and it makes no differential-expression
  calls.

## Limitations and caveats

- **Small-cohort trap.** The single most common misuse is underpowered
  discovery.
  [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  warns when events per signature gene fall below
  `min_events_per_parameter` (default 5, mid-range of the 5-9
  events-per-variable guidance; see
  [`?build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  for the sources). Treat concordance and hazard ratios from such
  cohorts as indicative, not confirmatory.
- **Linear Cox risk scores only.** The signature is a weighted linear
  combination of expression; nonlinear or machine-learning score forms
  are not produced.
- **No correction step.** Confounders such as batch are flagged, and
  screening can be adjusted for design terms, but the pipeline never
  corrects or removes effects itself.
- **Serial by default; parallelism opt-in.** Screening runs are
  sequential unless you pass a `BiocParallelParam`
  (e.g. `BiocParallel::SnowParam(2)`) to
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)/[`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  via `BPPARAM`. Parallel fold evaluation is bit-identical to serial
  (same seeded folds), and memory is usually the binding constraint on
  large cohorts regardless.

## Case studies

Every case study is written as a three-tier narrative: a **naive
analysis** (what a standard, unguarded pipeline would report),
**rnaSentry standing guard** (the same data through the pipeline, with
the flag ledger), and the **counterfactual** (the error prevented and
the true signal recovered).
[`vignette("case-study-impact")`](https://ghrieb.github.io/rnaSentry/articles/case-study-impact.md)
is the landing page with the master impact table. Three studies are
bundled so the pipeline’s behaviour can be inspected without pulling
data from the network; the fourth (LUAD) is reproduced from live GEO and
documented here.

| Failure mode | Naive headline | Guarded headline |
|----|----|----|
| Wrong direction (GSE20685) | inverted-convention read (C ~ 0.2) looks like a failure | held-out C = 0.783 vs 0.582 random-gene control |
| Confounded design (synthetic cohort) | six batch-driven genes are “prognostic”; batch + region are “independent” | batch/region redundant (Cramer’s V = 0.82); design reduced to `~ batch` |
| Underpowered discovery (30 samples, 23 events) | “mean CV C = 0.74 - decent” | fold range 0.50-1.00, sd 0.14, 4.6 events/parameter |
| False transfer (GSE31210 -\> GSE50081) | “CV C = 0.86, validated” | external C = 0.540, p = 0.266: transfer not demonstrated |

### 1. Breast cancer signature discovery (GSE20685)

A 3000-gene x 327-sample subset of GSE20685 (Li *et al.*, 2010;
Affymetrix GPL570, 83 deaths) ships as
`inst/extdata/gse20685_case_study.rds` and is walked through in
[`vignette("case-study-brca")`](https://ghrieb.github.io/rnaSentry/articles/case-study-brca.md).
On the bundled subset the pipeline finds a 20-gene signature with
screening-internal CV C = 0.80 (sd 0.03) and log-rank p = 3.0e-15 at the
discovery cutpoint; the design audit flags `subtype` (eta-squared =
0.78) and recommends adjusting for it. The narrative contrasts this with
the naive default-convention concordance read (C ~ 0.2, i.e. “below
chance”), which the `reverse = TRUE` convention and the `C + C_rev = 1`
parity test prevent – the exact mechanism behind the historical
concordance-inversion bug documented in
`validation/face_validity_review.md`. The signature is provisional: at
83 events for 20 genes the `events_per_parameter` guardrail fires (4.2
\< 5), so this demonstrates the pipeline, not a validated biomarker.
There is no independent breast-cohort validation here; external
validation is shown on an independent cohort in the LUAD narrative (case
study 4).

### 2. Confounder stress test

[`vignette("case-study-confounder-audit")`](https://ghrieb.github.io/rnaSentry/articles/case-study-confounder-audit.md)
plants a batch effect and a redundant, near-collinear design variable
(Cramer’s V = 0.82), plus a survival signal that runs through the batch.
The naive tier shows the six batch-driven genes ranking as the most
“prognostic” genes in an unadjusted screen; the guarded tier shows
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
flagging both variables as associated with the expression surrogate and
recommending the minimal design `~ batch`, after which the
batch-artifact genes drop out of screening.

### 3. The small-cohort guardrail

[`vignette("case-study-small-cohort")`](https://ghrieb.github.io/rnaSentry/articles/case-study-small-cohort.md)
shows why the `events_per_parameter` guardrail exists: on 30 samples /
23 events the naive reading quotes the mean CV C = 0.74, while the
fold-level table exposes a range of **0.50-1.00** (sd 0.14). The
guardrail fires (`4.6 < 5`) rather than quietly reporting unstable
numbers.

### 4. LUAD honest negative (GSE31210 -\> GSE50081)

A full live-GEO reproduction of a discovery-to-external cross-cohort
attempt: discovery in GSE31210 (LUAD, 35 deaths), external validation in
GSE50081 (n = 128, 52 deaths). Mechanics pass end to end, but external C
= 0.540 with log-rank p = 0.266 — transfer was **not demonstrated**,
reported as an honest negative rather than a tool failure. The companion
finding is the **screening-internal CV optimism**: on this ~21k-gene
panel even survival-permuted data yields a null CV C ~ 0.8 (not 0.5), so
CV concordance must be read relative to a matched null and only
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
is fully out-of-sample (see “Reading the CV concordance”). The power
analysis (`validation/simulate_study.R`, Sim 6) quantifies why a
35-event external cohort is structurally underpowered for weak real
signatures: transfer power \< 0.30 at ~35 events, and \>= 0.80 only near
~300 events for a C ~ 0.65 signature.

## Positioning and prior art

rnaSentry’s niche is **signature discovery with an enforced audit trail
and honest failure modes**: quality control, sample-identity checks,
confounder detection, and explicit gating are first-class product
features rather than side effects. Two Bioconductor neighbours are
closest:

- **asuri** (Bioc 3.23, `10.18129/B9.bioc.asuri`) — de novo
  survival-marker discovery from a `SummarizedExperiment` via
  subsampling glmnet + univariate Cox. The closest methodological
  analogue; rnaSentry differs by design in its sparse, directly
  inspectable linear score (univariate screen + joint Cox + repeated
  stratified CV), its concordance/external-validation focus, and its
  flag-ledger audit layer.
- **signifinder** (Bioc 3.16, `10.18129/B9.bioc.signifinder`) — applies
  60+ published cancer signatures as single-sample scores (GSVA;
  single-cell and spatial support). Complementary rather than competing:
  it answers “what do known signatures say about this sample?” where
  rnaSentry answers “can a new survival signature be learned from this
  cohort?”. It can be used downstream to benchmark a discovered
  signature.

Two prior-art anchors for the single-sample scoring approach itself:
**SurvMarker** (Gammune & Gu, 2025, DOI 10.64898/2025.12.31.697184) uses
PCA-based weighted scoring, and the **mRNAsi** stemness index (Malta *et
al.*, *Cell* 2018;173:338-354.e15) uses OCLR-based scoring. Full
profiles and a comparison table live in `validation/related_tools.md`.

## Contributing

Contributions are welcome — see `CONTRIBUTING.md` for the reporting
process, the developer workflow, and the gate suite that every change
must pass.

## Documentation

- Vignette: `browseVignettes("rnaSentry")` —
  [`vignette("getting-started")`](https://ghrieb.github.io/rnaSentry/articles/getting-started.md)
  bridges two CSV files into a `SummarizedExperiment`;
  [`vignette("rnaSentry")`](https://ghrieb.github.io/rnaSentry/articles/rnaSentry.md)
  is the main walkthrough;
  [`vignette("case-study-impact")`](https://ghrieb.github.io/rnaSentry/articles/case-study-impact.md)
  is the naive-vs-guarded landing page;
  [`vignette("case-study-brca")`](https://ghrieb.github.io/rnaSentry/articles/case-study-brca.md),
  [`vignette("case-study-confounder-audit")`](https://ghrieb.github.io/rnaSentry/articles/case-study-confounder-audit.md),
  and
  [`vignette("case-study-small-cohort")`](https://ghrieb.github.io/rnaSentry/articles/case-study-small-cohort.md)
  are the bundled case studies.
- Function reference:
  [`help(package = "rnaSentry")`](https://ghrieb.github.io/rnaSentry/reference)
- Validation dossier (sign-convention fix, parity tests, audit table,
  simulation incl. the power analysis, face validity, maintainer testing
  guide): see the `validation/` directory of the source package.

## Reporting issues

Please file bugs and feature requests at
<https://github.com/Ghrieb/rnaSentry/issues>.
