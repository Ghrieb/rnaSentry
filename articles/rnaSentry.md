# Introduction to rnaSentry

## Introduction

`rnaSentry` automates the statistical pipeline used to go from a raw
bulk RNA-seq count matrix to a validated prognostic gene signature:
quality control and sample-identity auditing, confounder detection prior
to differential expression, signature construction with
cross-validation, survival modeling, and external validation. Every
automated decision is reported with the statistic that justified it.

This vignette demonstrates the pipeline’s first stage: quality control
and biological-sex verification.

Three bundled case studies show the pipeline on real and adversarial
data and are meant to be read alongside this walkthrough:
[`vignette("case-study-brca")`](https://ghrieb.github.io/rnaSentry/articles/case-study-brca.md)
(the full pipeline on a breast-cancer cohort),
[`vignette("case-study-confounder-audit")`](https://ghrieb.github.io/rnaSentry/articles/case-study-confounder-audit.md)
(a planted batch confound), and
[`vignette("case-study-small-cohort")`](https://ghrieb.github.io/rnaSentry/articles/case-study-small-cohort.md)
(why small cohorts fail silently).

## Installation

``` r

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("rnaSentry")
```

## Design decisions

Two deliberate choices in the pipeline are worth stating explicitly,
because they are easy to mistake for oversights.

First,
[`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
screens genes *adjusted* for the design terms supplied via
`design_terms` (by default), so that genes whose survival association is
driven entirely by a confounder such as batch are not selected. The
final signature weights, however, are estimated from a gene-only model.
This asymmetry is intentional: the risk score must remain portable to
external cohorts that do not record the discovery cohort’s design
covariates. Adjusted inference on the chosen signature is provided
separately by
[`cox_model()`](https://ghrieb.github.io/rnaSentry/reference/cox_model.md).

Second,
[`pca_audit()`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
reports both raw and Benjamini-Hochberg adjusted p-values for its per-PC
batch tests but flags on the raw p-value, whereas
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
flags candidate confounders on BH-adjusted p-values. This, too, is
intentional:
[`pca_audit()`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
tests at most a handful of PCs against a single batch variable, where a
multiple-testing correction would be a near no-op, while
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
scans many candidate confounders at once and needs the correction to
control spurious flags.

## Data intake

[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
validates a raw count matrix and its sample metadata before the pipeline
runs. It returns a `SummarizedExperiment` and records any intake issues
in a flag ledger in the object’s metadata, so problems such as duplicate
sample identifiers or non-integer counts are surfaced up front and
carried into the rendered report’s audit trail.

``` r

library(rnaSentry)
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: 'MatrixGenerics'
#> The following objects are masked from 'package:matrixStats':
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
#> Attaching package: 'generics'
#> The following objects are masked from 'package:base':
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: 'BiocGenerics'
#> The following objects are masked from 'package:stats':
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from 'package:base':
#> 
#>     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
#>     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
#>     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
#>     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
#>     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: 'S4Vectors'
#> The following object is masked from 'package:utils':
#> 
#>     findMatches
#> The following objects are masked from 'package:base':
#> 
#>     expand.grid, I, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: 'Biobase'
#> The following object is masked from 'package:MatrixGenerics':
#> 
#>     rowMedians
#> The following objects are masked from 'package:matrixStats':
#> 
#>     anyMissing, rowMedians

cnt <- matrix(rpois(600, lambda = 200), nrow = 60, ncol = 10,
              dimnames = list(paste0("gene", 1:60), paste0("S", 1:10)))
cnt <- cnt * 1.0
cnt[1, ] <- 0L                 # gene1 has no counts anywhere (empty row)
cnt[2, 1] <- cnt[2, 1] + 0.5   # a non-integer value
colnames(cnt)[10] <- "S1"      # a duplicate sample identifier
cd <- data.frame(condition = rep(c("tumor", "normal"), each = 5),
                 row.names = paste0("S", 1:10))

se_in <- load_counts(cnt, cd)
#> Warning: Duplicate sample IDs found: S1. Removing duplicates.
#> Warning: Count matrix contains non-integer values. Downstream stages expect
#> integer counts.
#> Warning: 1 gene(s) have all-zero or all-NA counts and will be retained but may
#> be filtered later.
S4Vectors::metadata(se_in)$flags
#>                check severity                                    detail
#> 1  duplicate_samples critical                 Duplicate sample IDs: S1.
#> 2 non_integer_counts  warning Count matrix contains non-integer values.
#> 3         empty_rows  warning 1 gene(s) have all-zero or all-NA counts.
#>         stage
#> 1 load_counts
#> 2 load_counts
#> 3 load_counts
```

Each row of the ledger is a `(check, severity, detail, stage)` record;
the duplicate sample is dropped (critical), while the non-integer and
empty-row issues are flagged as warnings and carried forward.

## Quality control

[`qc_explore()`](https://ghrieb.github.io/rnaSentry/reference/qc_explore.md)
takes a `SummarizedExperiment` of raw counts and audits it for duplicate
sample identifiers, missing metadata, non-integer count values, and
library-size outliers.

``` r

set.seed(1)
counts <- matrix(rpois(600, lambda = 200), nrow = 60, ncol = 10,
                  dimnames = list(paste0("gene", 1:60), paste0("S", 1:10)))
counts[, 1] <- counts[, 1] * 50L
coldata <- DataFrame(condition = rep(c("tumor", "normal"), each = 5),
                      row.names = paste0("S", 1:10))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

qc <- qc_explore(se)
qc
#> rnaSentry QC audit: 10 samples x 60 genes
#> 1 issue(s) flagged:
#>   [warning] library_size_outlier: Samples with library size > 3 MADs from the median: S1
```

## Biological sex verification

[`sex_check()`](https://ghrieb.github.io/rnaSentry/reference/sex_check.md)
compares reported sex in sample metadata against sex inferred from XIST
and Y-chromosome marker gene expression, to catch sample swaps or
metadata errors before they propagate downstream.

``` r

genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D")
expr <- matrix(c(950, 5, 5, 5,   # sample looks female
                  5, 400, 380, 410),  # sample looks male
                nrow = 4, dimnames = list(genes, c("S1", "S2")))
coldata2 <- DataFrame(sex = c("F", "F"), row.names = c("S1", "S2")) # S2 mislabeled
se2 <- SummarizedExperiment(assays = list(counts = expr), colData = coldata2)

sex_check(se2)
#>   sample_id reported_sex inferred_sex   status
#> 1        S1            F            F       OK
#> 2        S2            F            M MISMATCH
```

## The survival pipeline in one run

[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
orchestrates the full pipeline on a single cohort: confounder scan, PCA
audit, signature construction with cross-validation, Kaplan-Meier risk
groups, an adjusted Cox model, and parametric survival models. Here we
build a synthetic cohort with an engineered survival signal, split it
into discovery and validation halves, and run the pipeline on the
discovery half.

``` r

set.seed(1)
n_genes <- 120
n_samples <- 120
counts <- matrix(stats::rpois(n_genes * n_samples, lambda = 500),
                 nrow = n_genes, ncol = n_samples,
                 dimnames = list(paste0("gene", 1:n_genes),
                                 paste0("S", 1:n_samples)))
sig_expr <- colMeans(counts[1:5, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- stats::rexp(n_samples, rate = 0.03 * exp(0.8 * risk))
censor_time <- stats::rexp(n_samples, rate = 0.02)
coldata <- DataFrame(time = pmin(event_time, censor_time),
                     event = as.integer(event_time < censor_time),
                     condition = factor(rep(c("A", "B"), length.out = n_samples)),
                     batch = factor(rep(c("B1", "B2", "B3"), length.out = n_samples)),
                     row.names = colnames(counts))
se_all <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

discovery <- se_all[, 1:72]
external <- se_all[, 73:120]

run <- run_rnaSentry(discovery, time_col = "time", event_col = "event",
                     top_n = 5, repeats = 2, folds = 3, seed = 7,
                     design_vars = c("condition", "batch"),
                     design_terms = "condition")
run
#> rnaSentry pipeline run.
#>   - design_audit
#>   - pca_audit
#>   - build_signature
#>   - km_curve
#>   - cox_model
#>   - survival_parametric
#> Report: /home/runner/work/rnaSentry/rnaSentry/vignettes/rnaSentry_report.html
```

The run locked the discovered signature:
[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
is single-use per session, and a second call in the same session is
refused until the lock is released. This is the reproducibility
guardrail against silently re-selecting signature genes after survival
analysis.

``` r

# Release the lock so the pipeline can be run again in this session.
run$stages$build_signature <-
  lock_signature(run$stages$build_signature, lock = FALSE)
```

The result’s `stages` list holds each stage’s full result, including the
`flags` data.frame that records every issue raised:

``` r

run$stages$build_signature$flags
#>                check severity
#> 1 adjusted_screening     info
#> 2   signature_locked     info
#> 3 signature_unlocked     info
#>                                                         detail           stage
#> 1 Univariate screening adjusted for design term(s): condition. build_signature
#> 2                Signature locked against downstream mutation.  lock_signature
#> 3          Signature unlocked; downstream stages may be rerun.  lock_signature
```

### External validation

The discovery cutpoint is recorded by the `km_curve` stage. Re-scoring
the held-out samples against that same cutpoint, without re-estimating
anything from the external data, gives an honest portability check:

``` r

val <- validate_external(run$stages$build_signature, external,
                         cutpoint = run$stages$km_curve$cutpoint)
val
#> rnaSentry external validation: 5/5 signature gene(s) used, 48 samples.
#> Risk groups: low n = 23 (12 events), high n = 25 (11 events).
#> Cutpoint 55.5 (custom, from discovery).
#> Log-rank p = 0.8303. Concordance = 0.485.
#> Signature not locked.
```

A real discovery-to-external transfer on independent GEO cohorts — where
transfer was *not* demonstrated and the reasons (low event count,
CV-optimism) are quantified — is documented in the README’s case study 4
(LUAD, GSE31210 -\> GSE50081).

## Assumptions and limitations

`rnaSentry` is deliberately narrow. Read this before trusting any
output.

- **Bulk, right-censored survival only.** The pipeline expects bulk
  RNA-seq and standard right-censored (time, event) end points.
  Single-cell data, continuous or binary outcomes, competing risks,
  time-varying covariates, and left truncation are out of scope.
- **Linear Cox risk scores.** A signature is a weighted linear
  combination of gene expression. Higher scores mean higher hazard
  (shorter survival); all concordance estimates use this convention.
- **No batch correction, no differential expression.** Confounders are
  *detected and reported*, and univariate screening can be *adjusted*
  for design terms, but the pipeline never corrects or removes effects
  and makes no DE calls.
- **The small-cohort trap.** With fewer than ~5 events per signature
  gene, cross-validated concordance and hazard ratios are unstable.
  [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  warns (`events_per_parameter`); a flagged signature should be treated
  as hypothesis-generating, not confirmatory.
- **The CV concordance is screening-internal.** Gene selection runs on
  the full cohort before the CV split, so
  [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)’s
  concordance is optimistic: on a large panel even survival-permuted
  data yields a high null (≈0.8 on ~21k genes) rather than 0.5. Read it
  relative to a matched null (permutation or random-gene control), and
  rely on
  [`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
  on an independent cohort as the fully out-of-sample estimate.
- **Right input scale.** Raw counts are log2-transformed internally. If
  you store already-logged data under a name like `counts`, it will be
  double-logged; the pipeline warns (`possibly_log_scaled`). Store
  pre-scaled data under the name `logcounts`.
- **One discovery per session.**
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  locks the discovered signature and refuses to re-run until
  `lock_signature(sig, lock = FALSE)` releases the lock. This prevents
  silent gene-set re-selection after survival analysis. To analyze a
  second cohort in the same session, unlock the first signature first.

### Avoiding common mistakes

Two misuses are easy to hit and easy to avoid.

First, wrong input scale. If your expression matrix is already
log-transformed but stored under a `counts`-named assay, every
downstream stage sees a double-logged matrix. The pipeline catches this
and warns:

``` r

prelogged <- log2(as.matrix(assay(se_all, "counts")) + 1)
se_wrong <- se_all
assay(se_wrong, "counts") <- prelogged
build_signature(se_wrong, "time", "event", top_n = 5, repeats = 1,
                folds = 2, seed = 7, min_events_per_parameter = 1)
#> Warning: The first assay looks already log-transformed (non-integer values with
#> a maximum below 40). rnaSentry will log2-transform it again, double-logging
#> expression. Rename the assay to "logcounts" or supply raw integer counts.
#> rnaSentry signature: 5 gene(s), outcome 'overall_survival'.
#> Selection: top_n (5 gene(s)). CV concordance: mean 0.626, sd 0.014.
#> Not locked.
#> 1 issue(s) flagged:
#>   [warning] possibly_log_scaled: The first assay looks already log-transformed (non-integer values with a maximum below 40). rnaSentry will log2-transform it again, double-logging expression. Rename the assay to "logcounts" or supply raw integer counts.
```

The fix is to name such an assay `logcounts`, which the pipeline uses
as-is without a second transform.

Second, underpowered discovery. On a small cohort,
[`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
warns that the event count per signature gene is low:

``` r

small <- se_all[, 1:40]
build_signature(small, "time", "event", top_n = 10, repeats = 1,
                folds = 2, seed = 7)
#> Warning: Only 22 event(s) for a 10-gene signature (2.2 events per parameter),
#> below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard
#> ratios are unstable at this event count; consider fewer genes or a larger
#> cohort.
#> rnaSentry signature: 10 gene(s), outcome 'overall_survival'.
#> Selection: top_n (10 gene(s)). CV concordance: mean 0.647, sd 0.073.
#> Not locked.
#> 1 issue(s) flagged:
#>   [warning] events_per_parameter: Only 22 event(s) for a 10-gene signature (2.2 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
```

Both warnings are also recorded in each stage’s `flags` data.frame, so
they propagate into the rendered report’s audit trail.

## Session info

``` r

sessionInfo()
#> R version 4.6.1 (2026-06-24)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats4    stats     graphics  grDevices utils     datasets  methods  
#> [8] base     
#> 
#> other attached packages:
#>  [1] SummarizedExperiment_1.42.0 Biobase_2.72.0             
#>  [3] GenomicRanges_1.64.0        Seqinfo_1.2.0              
#>  [5] IRanges_2.46.0              S4Vectors_0.50.1           
#>  [7] BiocGenerics_0.58.1         generics_0.1.4             
#>  [9] MatrixGenerics_1.24.0       matrixStats_1.5.0          
#> [11] rnaSentry_0.99.0            BiocStyle_2.40.0           
#> 
#> loaded via a namespace (and not attached):
#>  [1] Matrix_1.7-5        jsonlite_2.0.0      compiler_4.6.1     
#>  [4] BiocManager_1.30.27 jquerylib_0.1.4     splines_4.6.1      
#>  [7] systemfonts_1.3.2   textshaping_1.0.5   yaml_2.3.12        
#> [10] fastmap_1.2.0       lattice_0.22-9      XVector_0.52.0     
#> [13] R6_2.6.1            S4Arrays_1.12.0     knitr_1.51         
#> [16] htmlwidgets_1.6.4   DelayedArray_0.38.2 bookdown_0.47      
#> [19] desc_1.4.3          bslib_0.12.0        rlang_1.3.0        
#> [22] cachem_1.1.0        xfun_0.60           fs_2.1.0           
#> [25] sass_0.4.10         otel_0.2.0          SparseArray_1.12.2 
#> [28] cli_3.6.6           withr_3.0.3         pkgdown_2.2.1.9000 
#> [31] grid_4.6.1          digest_0.6.39       lifecycle_1.0.5    
#> [34] evaluate_1.0.5      ragg_1.5.2          survival_3.8-6     
#> [37] abind_1.4-8         rmarkdown_2.31      tools_4.6.1        
#> [40] htmltools_0.5.9
```
