# Build a prognostic gene signature from a survival cohort

`build_signature()` is the core signature-construction stage of the
rnaSentry pipeline. From a `SummarizedExperiment` with per-sample
survival metadata it:

1.  *Screens* every variable gene with a univariate Cox
    proportional-hazards model against the event indicator, recording
    the hazard ratio, Wald p-value and a Benjamini-Hochberg adjusted
    p-value.

2.  *Selects* the signature genes by the most significant univariate
    p-values (`method = "top_n"`) or by an adjusted p-value threshold
    (`method = "p_value"`).

3.  *Estimates joint coefficients* by refitting a single multivariate
    Cox model on the selected genes. Genes whose coefficient is not
    finite (rank deficiency) are dropped and the model refit; the
    signature is never silently degraded, the drop is reported as a
    flag.

4.  *Reports stability* with repeated, event-stratified k-fold
    cross-validation: within each fold the coefficients are re-estimated
    on the training samples and the concordance index of the resulting
    risk score is evaluated on the held-out samples. The mean and
    standard deviation of the per-fold concordance are returned as
    `cv_summary`. Concordance is Harrell's C in the risk-score
    convention: a higher risk score is concordant with an earlier event.

## Usage

``` r
build_signature(
  se,
  time_col,
  event_col,
  outcome_col = "overall_survival",
  method = c("top_n", "p_value"),
  top_n = 20,
  p_threshold = 0.05,
  repeats = 5,
  folds = 5,
  seed = NULL,
  design_terms = character(0),
  adjust_for_design = TRUE,
  min_events_per_parameter = 5,
  BPPARAM = NULL
)
```

## Arguments

- se:

  A `SummarizedExperiment` with a count or normalized expression assay
  and survival metadata.

- time_col:

  Character. Column of `colData(se)` with follow-up time (numeric,
  non-negative).

- event_col:

  Character. Column of `colData(se)` with the event indicator (0/1 or
  logical).

- outcome_col:

  Character. A display label for the outcome, e.g. `"overall_survival"`.
  This is recorded in the signature and used in report titles; it does
  *not* need to name a column of `colData(se)`.

- method:

  Character. Selection rule, one of `"top_n"` (default; take the `top_n`
  genes with the most significant univariate p-value) or `"p_value"`
  (take every gene with adjusted p-value below `p_threshold`).

- top_n:

  Integer. Number of genes to keep when `method = "top_n"`. Defaults to
  `20`.

- p_threshold:

  Numeric in (0, 1). Adjusted p-value cutoff when `method = "p_value"`.
  Defaults to `0.05`.

- repeats:

  Integer. Number of cross-validation repeats. Defaults to `5`.

- folds:

  Integer. Number of folds per repeat. Defaults to `5`.

- seed:

  Optional integer. Seeds the cross-validation fold shuffling for
  reproducible results. Seeding is scoped with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
  so the caller's global random-number generator state is left
  unchanged.

- design_terms:

  Character vector of `colData(se)` columns to adjust the univariate
  screening by. When non-empty and `adjust_for_design = TRUE`, each gene
  is tested with `Surv ~ gene + design_terms`, so that genes whose
  survival association is driven entirely by a confounder (for example
  batch) are not selected. The final signature weights (the joint Cox
  coefficients) are deliberately estimated from a gene-only model so
  that the risk score remains portable to cohorts that do not record the
  design covariates; adjusted inference on the chosen signature is
  provided by
  [`cox_model()`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md).
  Defaults to `character(0)`.

- adjust_for_design:

  Logical. When `TRUE` (default) and `design_terms` is non-empty,
  univariate screening is adjusted for the design terms as described
  above. Setting it to `FALSE` runs unadjusted screening and raises an
  `unadjusted_screening` warning, because the pipeline's confounder
  guardrails are then bypassed.

- min_events_per_parameter:

  Numeric. Minimum acceptable number of events per signature gene. When
  the cohort has fewer events per final signature gene, the signature is
  flagged with an `events_per_parameter` warning, because
  cross-validated concordance and hazard ratios are unstable at very low
  event counts. The default of `5` sits below the classic 10
  events-per-variable rule of thumb (Peduzzi et al., *J Clin
  Epidemiol* 1996) but at the lower end of the 5-9 events-per-variable
  range that Vittinghoff and McCulloch (*Am J Epidemiol* 2007) argued
  can be adequate in some settings; it is a defensible middle choice,
  not a guarantee of estimability. Raise it to `10` for the stricter
  convention, or lower it to suppress the warning on small but
  well-behaved cohorts.

- BPPARAM:

  Optional `BiocParallelParam` from the BiocParallel package (for
  example `BiocParallel::SnowParam(2)`). When `NULL` (the default) the
  repeated cross-validation is evaluated serially. When a non-`NULL`
  backend is supplied, the per-fold evaluations are distributed with
  [`BiocParallel::bplapply()`](https://rdrr.io/pkg/BiocParallel/man/bplapply.html);
  BiocParallel is only a suggested package and is loaded on demand.
  Parallelism is strictly opt-in and is bit-identical to the serial
  path: all randomness is confined to the fold assignment (seeded with
  `seed`), and the per-fold Cox fits and concordance scores use no
  random numbers, so the returned signature, cross-validation table, and
  flag ledger are unchanged regardless of backend.

## Value

An object of class `"rnaSentry_signature"` (a list) with elements:

- genes:

  Character vector of the final signature genes.

- outcome:

  The `outcome_col` display label.

- time_col, event_col, design_terms:

  The metadata columns used.

- screening_terms:

  Character vector of the design terms actually used to adjust
  univariate screening (empty when screening was unadjusted).

- cox_stats:

  Data.frame with one row per screened gene (columns `gene`, `HR`, `p`,
  `adj_p`).

- coefficients:

  Named numeric vector of joint Cox coefficients for the signature
  genes.

- cv_results:

  Data.frame with columns `repeat_id`, `fold` and `c_index` (concordance
  on held-out samples).

- cv_summary:

  Named numeric vector `c(mean, sd)` of the per-fold concordance.

- selection:

  List with `method` and `n_genes`.

- locked:

  Always `FALSE`; set by
  [`lock_signature()`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md).

- created, lock_time:

  Provenance timestamps.

- assay_used, n_screened, n_filtered:

  Analysis details.

- flags:

  Data.frame summarizing every issue raised, with columns `check`,
  `severity`, `detail` and `stage`.

## Details

**Reading the cross-validated concordance (important).** The univariate
screening that selects the signature runs on the full cohort *before*
the cross-validation split, so the CV fold is *not* fully untouched by
the selection step. On large candidate panels the resulting concordance
is therefore *optimistic* (screening-internal): even data with permuted
survival yields a non-trivial null concordance (e.g. `~0.8` on a
~21k-gene panel, driven by winner's-curse selection), far above the 0.5
a fully out-of-sample estimate would give. Treat `cv_summary` as a
screening-internal stability metric and interpret it *relative* to a
matched null — a permutation-null or random-gene control — rather than
against 0.5 in absolute terms. The only fully independent out-of-sample
estimate is
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
on a separate cohort.

The return value has class `"rnaSentry_signature"` and is the input
expected by the survival modeling and validation stages
([`km_curve()`](https://ghrieb.github.io/rnaSentry/reference/km_curve.md),
[`cox_model()`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md),
[`survival_parametric()`](https://ghrieb.github.io/rnaSentry/reference/survival_parametric.md),
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md))
and by
[`lock_signature()`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md).
The signature is returned *unlocked*; downstream stages that must not be
rerun after publication can be protected with
[`lock_signature()`](https://ghrieb.github.io/rnaSentry/reference/lock_signature.md).

## Assumptions and limitations

The signature is a linear Cox risk score (larger score = earlier event).
The pipeline assumes standard right-censored survival from bulk RNA-seq
with raw integer counts or a correctly-named `logcounts` assay
(pre-scaled data in a `counts`-named slot is flagged as
`possibly_log_scaled`). It does not model competing risks, time-varying
covariates, left truncation, single-cell data, or non-survival end
points, and it never corrects batch or design effects. In underpowered
cohorts (fewer than `min_events_per_parameter` events per signature
gene) the returned concordance and hazard ratios are unstable and an
`events_per_parameter` warning is raised; treat such results as
hypothesis-generating rather than confirmatory. The discovery stage
additionally refuses cohorts with fewer than two events; the downstream
survival stages
([`km_curve`](https://ghrieb.github.io/rnaSentry/reference/km_curve.md),
[`cox_model`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md),
[`survival_parametric`](https://ghrieb.github.io/rnaSentry/reference/survival_parametric.md),
[`validate_external`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md))
use a one-event minimum and instead surface instability through their
`sparse_events` guardrail.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
set.seed(5)
counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
                  dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
sig_expr <- colMeans(counts[1:5, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(20, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 condition = factor(rep(c("A", "B"), 10)),
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

sig <- build_signature(se, time_col = "time", event_col = "event",
                       top_n = 5, repeats = 1, folds = 2, seed = 1)
#> Warning: Only 13 event(s) for a 5-gene signature (2.6 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
sig
#> rnaSentry signature: 5 gene(s), outcome 'overall_survival'.
#> Selection: top_n (5 gene(s)). CV concordance: mean 0.412, sd 0.099.
#> Not locked.
#> 1 issue(s) flagged:
#>   [warning] events_per_parameter: Only 13 event(s) for a 5-gene signature (2.6 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
```
