# Naive analysis vs rnaSentry: what the guardrails prevent

## The three-tier structure

The case studies are written as three acts. Tier 1 runs the analysis the
way a standard, unguarded pipeline would, and prints the headline a
researcher would have believed. Tier 2 runs the *same data* through
`rnaSentry` and shows which guardrail fired, with the statistic that
justified it. Tier 3 states the counterfactual: the error the guardrails
prevented, and the true signal the guarded analysis still recovered.

Four failure modes – the ways a prognostic-signature study most often
goes wrong – are covered by the case studies:

| \# | Failure mode | Where it is demonstrated |
|----|----|----|
| 1 | **Wrong direction**: an inverted concordance convention turns a strong signature into a “failure” | [`vignette("case-study-brca")`](https://ghrieb.github.io/rnaSentry/articles/case-study-brca.md) |
| 2 | **Confounded design**: batch-driven genes masquerade as prognostic; collinear terms are double-counted | [`vignette("case-study-confounder-audit")`](https://ghrieb.github.io/rnaSentry/articles/case-study-confounder-audit.md) |
| 3 | **Underpowered discovery**: a mean C-index hides coin-flip folds | [`vignette("case-study-small-cohort")`](https://ghrieb.github.io/rnaSentry/articles/case-study-small-cohort.md) |
| 4 | **False transfer**: screening-internal CV optimism is mistaken for cross-cohort validity | README / main vignette (LUAD honest negative, below) |

## The master impact table

| Failure mode | Naive headline | Guarded headline | Guardrails that fired | What a guarded manuscript can claim |
|----|----|----|----|----|
| Wrong direction (GSE20685) | “CV C = 0.217, signature fails to generalize (winner’s curse)” | held-out C = **0.783**, direction correct, far above the 0.582 random-gene control | `reverse = TRUE` convention; `C + C_rev = 1` parity tests; direction-sensitive regression tests | the signature is not a 0.22 collapse but a strong, directionally correct held-out signal |
| Confounded design (synthetic batch cohort) | “six batch-driven genes are prognostic; batch and region are two independent risk factors” | `batch`/`region` flagged redundant (Cramer’s V = 0.82); `batch` surrogate eta-squared 0.97; minimal design `~ batch`; adjusted screening drops the batch-artifact genes | `redundant_variable`, surrogate-association flags, `adjusted_screening` | the discovered signature is not a batch echo |
| Underpowered discovery (30 samples, 23 events) | “cross-validated C = 0.74 – decent” | fold range **0.50-1.00**, sd 0.14, 4.6 events/parameter \< 5 | `events_per_parameter`, fold-level reporting | the concordance is an average over coin-flip folds; hypothesis-generating only |
| False transfer (GSE31210 -\> GSE50081, LUAD) | “CV C = 0.86, the signature is validated” | CV C = 0.862 vs permutation null 0.836 (delta +0.026); external C = **0.540**, log-rank p = **0.266** – transfer *not* demonstrated | permutation-null calibration; honest-negative external validation; power analysis | no false claim of cross-cohort validity; the negative is a power finding (35 events), not a tool failure |

## Why the first failure mode is worth dwelling on

A sign-convention bug is the one error class that flips a paper’s
conclusion without any red flag in the data.
[`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html)
defaults to *larger score means longer survival*; a Cox risk score means
the opposite. Calling the default on a risk score returns
`1 - Harrell's C`. The identity below is what the parity test suite now
asserts on every fold:

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
library(survival)

se <- readRDS(system.file("extdata", "gse20685_case_study.rds",
                          package = "rnaSentry"))
score <- colMeans(assays(se)$logcounts[1:5, , drop = FALSE])
c_default <- as.numeric(concordance(Surv(se$time, se$event) ~
                                    score)$concordance[1])
c_reverse <- as.numeric(concordance(Surv(se$time, se$event) ~ score,
                                    reverse = TRUE)$concordance[1])
c(c_default = c_default, c_reverse = c_reverse, sum = c_default + c_reverse)
#> c_default c_reverse       sum 
#> 0.4876938 0.5123062 1.0000000
```

The first real-data run of this package reported 0.217 and the review
briefly rationalized it as winner’s curse; `1 - 0.217 = 0.783` was the
correct, biologically expected value all along (full diagnosis in
`validation/face_validity_review.md`). The same class of error has
appeared in published survival literature; rnaSentry pins the convention
in code, in parity tests, and in the rendered report.

## The LUAD honest negative (GSE31210 -\> GSE50081)

The fourth failure mode is demonstrated on independent GEO cohorts in
the README and main vignette; the numbers are reproduced by the dev-only
script `validation/repro_live_geo.R`. Because the data are downloaded
live, the walkthrough is narrative here rather than a build-time
vignette:

- **Discovery** GSE31210 (LUAD, 35 deaths), **external validation**
  GSE50081 (n = 128, 52 deaths).
- Screening-internal CV C = 0.862 against a permutation null of 0.836
  (delta +0.026) – even survival-permuted data yields a high null CV
  (about 0.8 on ~21k genes), so the CV must be calibrated against a
  matched null.
- External validation: **C = 0.540, log-rank p = 0.266** – transfer was
  *not* demonstrated, and the tool reports an honest negative rather
  than a “validated signature.”
- The power analysis (`validation/simulate_study.R`, Sim 6) quantifies
  why: transfer power \< 0.30 at ~35 discovery events, and \>= 0.80 only
  near ~300 events for a C ~ 0.65 signature.

## Reading this package’s claims

- **Every number above is pinned.** The case-study vignettes build
  offline and their outputs match the values recorded in the validation
  dossier (`validation/paper_qa_summary.md`); the stale-number sweep
  (`validation/grep_stale_numbers.ps1`) fails the build if a superseded
  number leaks back into shipped docs outside the allowed correction
  narrative.
- **The guardrails are falsifiable.** Each flag in the report’s audit
  trail is a `(check, severity, detail, stage)` record backed by a
  statistic, and each has a regression test asserting the condition that
  fires it.
- **This is a demonstration, not a clinical claim.** All signatures are
  provisional until a preprint is published.

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
#>  [1] survival_3.8-6              SummarizedExperiment_1.42.0
#>  [3] Biobase_2.72.0              GenomicRanges_1.64.0       
#>  [5] Seqinfo_1.2.0               IRanges_2.46.0             
#>  [7] S4Vectors_0.50.1            BiocGenerics_0.58.1        
#>  [9] generics_0.1.4              MatrixGenerics_1.24.0      
#> [11] matrixStats_1.5.0           rnaSentry_0.99.0           
#> [13] BiocStyle_2.40.0           
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
#> [28] cli_3.6.6           pkgdown_2.2.1.9000  grid_4.6.1         
#> [31] digest_0.6.39       lifecycle_1.0.5     evaluate_1.0.5     
#> [34] ragg_1.5.2          abind_1.4-8         rmarkdown_2.31     
#> [37] tools_4.6.1         htmltools_0.5.9
```
