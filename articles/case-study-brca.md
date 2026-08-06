# Case study: a breast cancer prognostic signature (GSE20685)

## Status note

This case study is an illustrative walkthrough of `rnaSentry` on a
public dataset, provided to show how the pipeline behaves on real data.
It is *not* peer-reviewed and the signature reported here is
provisional. Treat every number below as a demonstration, not a clinical
claim.

The narrative follows the *three-tier* structure used across the case
studies: a **naive analysis** first, then **rnaSentry standing guard**,
then the **counterfactual** – what a researcher would have concluded
without the guardrails, and what the guarded analysis actually recovers.

## Data

We use a bundled subset of GSE20685 (Li *et al.*, 2010), an Affymetrix
GPL570 (U133 Plus 2.0) microarray series of 327 primary breast tumors
with overall-survival follow-up. The series matrix was downloaded from
GEO and processed by the repository’s validation script
`validation/repro_gse20685.R`:

- probe-level log2 intensities from the MAS5-summarized series matrix,
- probe-to-gene collapse by largest mean expression,
- the top 3000 most variable genes retained for the bundled subset,
- a `SummarizedExperiment` with the subset under the `logcounts` assay
  name, plus `time` (years), `event` (1 = death), `age`, and `subtype`
  in the sample metadata.

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
se
#> class: SummarizedExperiment 
#> dim: 3000 327 
#> metadata(0):
#> assays(1): logcounts
#> rownames(3000): CRISP3 S100A7 ... LINC00858 DSG4
#> rowData names(0):
#> colnames(327): GSM519117 GSM519118 ... GSM519442 GSM519443
#> colData names(4): time event age subtype
table(se$subtype)
#> 
#>   type I  type II type III  type IV   type V  type VI 
#>       37       34       41       81       41       93
sum(se$event)
#> [1] 83
```

## Tier 1: the naive analysis

A standard, unguarded analysis of this cohort would proceed along these
lines: pick a few highly variable genes, screen them against survival
with univariate Cox regressions, average the survivors into a risk
score, and quote the concordance that R prints by default. Let us run
exactly that.

``` r

lg <- assays(se)$logcounts
vars <- apply(lg, 1, stats::var, na.rm = TRUE)
cand <- names(sort(vars, decreasing = TRUE))[seq_len(200)]
stat <- do.call(rbind, lapply(cand, function(g) {
    s <- summary(coxph(Surv(se$time, se$event) ~ lg[g, ]))
    c(beta = s$coefficients[1, 1], p = s$coefficients[1, 5])
}))
rownames(stat) <- cand
top <- stat[order(stat[, "p"]), , drop = FALSE][seq_len(5), , drop = FALSE]
naive_score <- colMeans(sign(top[, "beta"]) * lg[rownames(top), , drop = FALSE])

c_default <- as.numeric(concordance(Surv(se$time, se$event) ~
                                    naive_score)$concordance[1])
c_reverse <- as.numeric(concordance(Surv(se$time, se$event) ~ naive_score,
                                    reverse = TRUE)$concordance[1])
c(c_default = c_default, c_reverse = c_reverse, sum = c_default + c_reverse)
#> c_default c_reverse       sum 
#> 0.3156447 0.6843553 1.0000000
```

[`survival::concordance()`](https://rdrr.io/pkg/survival/man/concordance.html)
defaults to the convention *larger score means longer survival*
(`reverse = FALSE`). A Cox risk score means the opposite – larger score
means *shorter* survival – so the default read is inverted:
`C_default = 1 - Harrell's C`. The naive analyst sees `C_default < 0.5`
and concludes “below chance, poor generalization, the signature is a
failure.”

This is not hypothetical. The first real-data face-validity run of this
very cohort reported CV concordance **0.217** and a naive reading would
rationalize it as extreme selection-induced winner’s curse; the correct
value is `1 - 0.217 = 0.783`, a strong held-out result in the
biologically expected direction. The full diagnosis, including the
`C + C_rev = 1` invariant you can see holding above, is recorded in the
repository’s validation dossier. The parity test suite now asserts that
invariant on every fold so the convention cannot silently flip again.

Even read in the correct direction, a naive analyst stops at the cross-
validated mean and treats it as an out-of-sample estimate. As the LUAD
honest-negative case study documents, a screening-internal CV on a large
candidate panel is *optimistic*: even survival-permuted data yields a
null CV around 0.8 on ~21k genes, not 0.5. Any concordance below must be
read relative to a matched null, and the only fully out-of-sample number
is
[`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md)
on an independent cohort.

## Tier 2: rnaSentry standing guard

Now the same cohort through the guarded pipeline: confounder scan, PCA
audit, signature construction with repeated cross-validation, and
Kaplan-Meier risk groups. Because the cohort metadata records `age` and
`subtype`, both are supplied as design variables so the audit can check
them against the expression surrogate.

``` r

run <- run_rnaSentry(se, time_col = "time", event_col = "event",
                     outcome_col = "overall_survival",
                     design_vars = c("age", "subtype"),
                     top_n = 20, repeats = 3, folds = 3, seed = 42,
                     render_report = FALSE)
#> Warning: Only 83 event(s) for a 20-gene signature (4.2 events per parameter),
#> below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard
#> ratios are unstable at this event count; consider fewer genes or a larger
#> cohort.
sig <- run$stages$build_signature
head(sig$genes)
#> [1] "C20orf96"    "RP11-28F1.2" "SPRR3"       "EPHX2"       "CD1E"       
#> [6] "CASP14"
sig$cv_summary
#>       mean         sd 
#> 0.79736356 0.02655456
sig$flags
#>                  check severity
#> 1 events_per_parameter  warning
#> 2     signature_locked     info
#>                                                                                                                                                                                                                                detail
#> 1 Only 83 event(s) for a 20-gene signature (4.2 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
#> 2                                                                                                                                                                                       Signature locked against downstream mutation.
#>             stage
#> 1 build_signature
#> 2  lock_signature
```

The run locks the discovered signature, so a second discovery in the
same session is refused until the lock is released (see the main
vignette). Because all vignettes build in one session, we release it
here so the remaining vignettes can run their own pipeline calls:

``` r

run$stages$build_signature <-
    lock_signature(run$stages$build_signature, lock = FALSE)
```

The guardrails that the naive analysis missed are now explicit in the
flag ledger:

- **Direction.** Every concordance in the pipeline is computed with
  `reverse = TRUE`; the sign convention cannot silently flip (a parity
  test asserts `C + C_rev = 1` on every fold).
- **Events per parameter.** The `events_per_parameter` warning is
  expected and matches the documented limitation: 83 deaths for a
  20-gene signature is 4.2 events per parameter, below the guardrail
  of 5. The signature should be read as hypothesis-generating, and the
  cross-validated concordance and hazard ratios as unstable at this
  event count.
- **Screening-internal CV.** The cross-validated concordance (mean C =
  0.80) is a *screening-internal* estimate: gene selection happens on
  the full cohort before the CV split, so it is optimistic. Read it
  relative to a matched null, as discussed in the main vignette and in
  the README.
- **Design.** The design audit below decides which confounders the
  screening adjusts for, instead of an ad-hoc pick.

### Risk groups

``` r

km <- run$stages$km_curve
km$log_rank_p
#> [1] 2.99146e-15
```

The discovered signature separates the cohort into risk groups whose
survival curves differ dramatically at the discovery cutpoint.

### Design audit

``` r

da <- run$stages$design_audit
da$confounder_table
#>   variable        type   n test   statistic     df       p_value effect_size
#> 1      age     numeric 327   lm   0.7155663 1, 325  3.982245e-01 0.002196906
#> 2  subtype categorical 327  aov 225.5113213 5, 321 9.806085e-103 0.778400099
#>   effect_size_type flagged         adj_p
#> 1        r_squared   FALSE  3.982245e-01
#> 2      eta_squared    TRUE 1.961217e-102
da$formula_text
#> [1] "~ age + subtype"
```

`subtype` is flagged as associated with the expression surrogate (very
large eta-squared) and the recommended design formula therefore keeps
it, so that univariate gene screening is adjusted for subtype. `age` is
not flagged. This is exactly the kind of design decision `rnaSentry` is
meant to surface before any differential-expression or signature work.

## Tier 3: the counterfactual

|  | Naive analysis (Tier 1) | rnaSentry (Tier 2) |
|----|----|----|
| Score direction | default convention reads the score as *protective*; `C < 0.5` looks like failure | `reverse = TRUE`; `C_rev = 1 - C_default`, direction correct |
| Headline concordance | `C_default` quoted as gospel | `C_rev` quoted with its convention, and flagged screening-internal |
| Out-of-sample claim | none / implicit | only [`validate_external()`](https://ghrieb.github.io/rnaSentry/reference/validate_external.md) counts; matched-null calibration required |
| Events per parameter | not reported | 4.2 \< 5, warned in the ledger |
| Design | ad-hoc | `subtype` flagged and adjusted for |

**What a guarded manuscript can claim:** on the full-cohort
face-validity run of this dataset rnaSentry reports a directionally
correct held-out concordance of 0.783 against a 0.582 random-gene
control (and a 0.490 simulated null) – genuine survival signal, not a
0.22 collapse – while the naive pipeline either reports a “failed”
signature (inverted convention) or an uncalibrated 0.80
(screening-internal optimism). The events-per-parameter limitation is
stated, not hidden.

## Limitations

- **No independent breast cohort is validated here.** External
  validation is demonstrated on an independent cohort in the LUAD
  honest-negative narrative in the README and in the main vignette; this
  walkthrough deliberately stops at in-cohort analysis so the numbers
  stay interpretable.
- **Screening-internal CV.** The concordance above is optimistic by
  design; see the main vignette’s “Assumptions and limitations”.
- **Low events per parameter.** The guardrail warning above is a real
  signal; this signature is a demonstration, not a published result.
- GSE20685 is an *Affymetrix* series; its `logcounts` are already log2
  intensities, so no second transform is applied (the `logcounts` assay
  name is used).

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
#> [28] cli_3.6.6           withr_3.0.3         pkgdown_2.2.1.9000 
#> [31] grid_4.6.1          digest_0.6.39       lifecycle_1.0.5    
#> [34] evaluate_1.0.5      ragg_1.5.2          abind_1.4-8        
#> [37] rmarkdown_2.31      tools_4.6.1         htmltools_0.5.9
```
