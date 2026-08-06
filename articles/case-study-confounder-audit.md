# Case study: the confounder audit on a batch-confounded cohort

## Status note

This is a synthetic-data walkthrough used to verify
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
behavior on a deliberately confounded cohort. It is a vignette, not a
study. The narrative follows the *three-tier* structure used across the
case studies: a **naive analysis** first, then **rnaSentry standing
guard**, then the **counterfactual**.

## Setup

We build a cohort of 180 samples with two design variables that are
almost interchangeable:

- `region` is a 3-level factor derived by cutting a latent normal
  variable;
- `batch` is a *noisy copy* of `region` (80% of samples keep their
  `region` label; the rest are shuffled).

We also plant a real batch effect in expression: six genes are
multiplied by a batch-dependent factor. Finally, survival depends on the
latent variable (so `region`, and through its overlap `batch`, are
genuine survival correlates) – a realistic setting in which an
unadjusted screen can mistake batch-driven genes for prognostic ones.

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

set.seed(12)
n <- 180
latent <- rnorm(n)
region <- factor(cut(latent, breaks = 3, labels = c("R1", "R2", "R3")))
batch <- as.character(region)
for (i in seq_len(n))
    if (runif(1) < 0.2) batch[i] <- sample(c("R1", "R2", "R3"), 1)
batch <- factor(batch)

counts <- matrix(rpois(50 * n, lambda = 200), nrow = 50, ncol = n,
                 dimnames = list(paste0("g", 1:50), paste0("S", 1:n)))
for (b in 1:3) {
    idx <- which(batch == levels(batch)[b])
    counts[1:6, idx] <- counts[1:6, idx] * (1 + 0.6 * b)
}

event_time <- rexp(n, rate = 0.05 * exp(0.8 * latent))
censor_time <- rexp(n, rate = 0.03)
cd <- DataFrame(batch = batch, region = region,
                time = pmin(event_time, censor_time),
                event = as.integer(event_time < censor_time),
                row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = cd)
se
#> class: SummarizedExperiment 
#> dim: 50 180 
#> metadata(0):
#> assays(1): counts
#> rownames(50): g1 g2 ... g49 g50
#> rowData names(0):
#> colnames(180): S1 S2 ... S179 S180
#> colData names(4): batch region time event
```

## Tier 1: the naive analysis

A naive analyst runs univariate Cox regressions for every gene against
survival and takes the top hits – without checking whether the candidate
genes are driven by batch, and without auditing the design.

``` r

lg <- log2(counts + 1)
stat <- do.call(rbind, lapply(rownames(lg), function(g) {
    s <- summary(coxph(Surv(cd$time, cd$event) ~ lg[g, ]))
    c(beta = s$coefficients[1, 1], p = s$coefficients[1, 5])
}))
rownames(stat) <- rownames(lg)
naive_top <- head(stat[order(stat[, "p"]), , drop = FALSE], 6)
data.frame(gene = rownames(naive_top),
           planted_batch_effect = rownames(naive_top) %in% paste0("g", 1:6),
           naive_top)
#>    gene planted_batch_effect     beta            p
#> g3   g3                 TRUE 2.155854 2.099216e-08
#> g5   g5                 TRUE 2.044039 6.186780e-08
#> g6   g6                 TRUE 1.896791 9.984742e-08
#> g4   g4                 TRUE 1.858827 7.106999e-07
#> g1   g1                 TRUE 1.778014 7.113127e-07
#> g2   g2                 TRUE 1.826484 9.162488e-07
```

The six genes whose expression we multiplied by the batch factor rank as
the most “prognostic” genes in the cohort. An unguarded screen would put
exactly these batch artifacts into a signature. A naive model that then
includes both of the near-identical design variables compounds the
problem by double-counting collinear terms:

``` r

summary(coxph(Surv(time, event) ~ batch + region, data = as.data.frame(cd)))
#> Call:
#> coxph(formula = Surv(time, event) ~ batch + region, data = as.data.frame(cd))
#> 
#>   n= 180, number of events= 103 
#> 
#>            coef exp(coef) se(coef)     z Pr(>|z|)    
#> batchR2  0.2114    1.2354   0.5220 0.405  0.68551    
#> batchR3  0.5957    1.8144   0.4873 1.222  0.22154    
#> regionR2 0.9230    2.5168   0.5293 1.744  0.08118 .  
#> regionR3 1.6504    5.2090   0.4966 3.323  0.00089 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#>          exp(coef) exp(-coef) lower .95 upper .95
#> batchR2      1.235     0.8095    0.4441     3.437
#> batchR3      1.814     0.5512    0.6981     4.716
#> regionR2     2.517     0.3973    0.8919     7.102
#> regionR3     5.209     0.1920    1.9679    13.788
#> 
#> Concordance= 0.695  (se = 0.023 )
#> Likelihood ratio test= 45.62  on 4 df,   p=3e-09
#> Wald test            = 42.64  on 4 df,   p=1e-08
#> Score (logrank) test = 50.55  on 4 df,   p=3e-10
```

`batch` and `region` carry almost the same information; fitting both
inflates standard errors and produces unstable coefficients that look
like two “independent” effects when there is really one.

## Tier 2: rnaSentry standing guard

### Design audit

``` r

da <- design_audit(se, design_vars = c("batch", "region"),
                   redundant_effect_size = 0.5)
```

#### Pairwise redundancy

``` r

da$pairwise_table
#>    var1   var2                 type_pair       test statistic df      p_value
#> 1 batch region categorical x categorical chisq.test  241.7757  4 3.846267e-51
#>   effect_size effect_size_type redundant
#> 1   0.8195115        cramers_v      TRUE
```

`batch` and `region` are almost interchangeable (Cramer’s V = 0.82), and
the pair is flagged `redundant`. Including both in a model would
double-count nearly the same information.

#### Confounder scan

``` r

da$confounder_table
#>   variable        type   n         test  statistic     df       p_value
#> 1    batch categorical 180          aov 2577.97331 2, 177 1.285386e-131
#> 2   region categorical 180 kruskal.test   98.67523      2  3.740642e-22
#>   effect_size effect_size_type flagged         adj_p
#> 1   0.9668101      eta_squared    TRUE 2.570773e-131
#> 2   0.5461877  epsilon_squared    TRUE  3.740642e-22
```

Both variables are flagged as associated with the expression surrogate:
`batch` with eta-squared 0.97, `region` with epsilon-squared 0.55 (it
leaks the same signal through its overlap with `batch`).

#### What the surrogate gradient looks like

The surrogate the confounder scan tests against is the leading principal
component of the expression data.
[`pca_audit()`](https://ghrieb.github.io/rnaSentry/reference/pca_audit.md)
performs the same scan against a single batch variable, one PC at a
time:

``` r

pa <- pca_audit(se, batch_col = "batch")
pa$batch_tests
#>    pc test    statistic     df       p_value effect_size effect_size_type
#> 1 PC1  aov 2577.9733053 2, 177 1.285386e-131 0.966810093      eta_squared
#> 2 PC2  aov    1.4053241 2, 177  2.480108e-01 0.015631155      eta_squared
#> 3 PC3  aov    0.1040554 2, 177  9.012304e-01 0.001174387      eta_squared
#> 4 PC4  aov    0.2356731 2, 177  7.902864e-01 0.002655900      eta_squared
#> 5 PC5  aov    0.2489179 2, 177  7.799165e-01 0.002804743      eta_squared
#>   flagged         adj_p
#> 1    TRUE 6.426932e-131
#> 2   FALSE  6.200270e-01
#> 3   FALSE  9.012304e-01
#> 4   FALSE  9.012304e-01
#> 5   FALSE  9.012304e-01
```

`batch` is flagged as associated with the leading PCs. Plotting the
scores makes the structure visible: the three batches separate along
PC1, which is exactly the surrogate gradient
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
measures:

``` r

plot_pca_audit(pa, color_by = "batch")
```

![PC1-PC2 scatter of the batch-confounded cohort, colored by batch. The
three batch clusters separate along PC1 (the surrogate gradient); this
is the structure the confounder scan
detects.](case-study-confounder-audit_files/figure-html/pca-batch-plot-1.png)

PC1-PC2 scatter of the batch-confounded cohort, colored by batch. The
three batch clusters separate along PC1 (the surrogate gradient); this
is the structure the confounder scan detects.

Coloring the same scores by `region` gives a nearly identical picture,
because `batch` is a noisy copy of `region`:

``` r

plot_pca_audit(pa, color_by = "region")
```

![The same PC1-PC2 scatter colored by region. The near-identical
structure is why the redundancy scan reports Cramer's V = 0.82 and the
recommended formula drops
\`region\`.](case-study-confounder-audit_files/figure-html/pca-region-plot-1.png)

The same PC1-PC2 scatter colored by region. The near-identical structure
is why the redundancy scan reports Cramer’s V = 0.82 and the recommended
formula drops `region`.

The two candidate variables track the same expression structure, so the
redundancy scan collapses them: keeping only `batch` in the design
captures the shared confound without double-counting collinear terms.

#### Recovered design

``` r

da$formula_text
#> [1] "~ batch"
```

The audit therefore recommends the minimal formula `~ batch`, dropping
the redundant `region` term. This is the behavior verified by the test
suite (`confounder-detection` and `design-audit` blocks): a redundant
design term is reported, not silently kept.

### Batch-adjusted screening

Armed with the recommended design, the guarded discovery run adjusts its
univariate screening for `batch`. The batch-driven genes are no longer
the top hits:

``` r

adj <- do.call(rbind, lapply(rownames(lg), function(g) {
    s <- summary(coxph(Surv(cd$time, cd$event) ~ lg[g, ] + batch,
                       data = as.data.frame(cd)))
    c(beta = s$coefficients[1, 1], p = s$coefficients[1, 5])
}))
rownames(adj) <- rownames(lg)
guarded_top <- head(adj[order(adj[, "p"]), , drop = FALSE], 6)
data.frame(gene = rownames(guarded_top),
           planted_batch_effect = rownames(guarded_top) %in% paste0("g", 1:6),
           guarded_top)
#>     gene planted_batch_effect      beta          p
#> g35  g35                FALSE -2.421041 0.01122042
#> g15  g15                FALSE  2.543141 0.01764476
#> g31  g31                FALSE  2.107756 0.02446470
#> g47  g47                FALSE -2.160102 0.06155177
#> g13  g13                FALSE  1.600798 0.09347755
#> g23  g23                FALSE -1.665415 0.11027676
```

Once `batch` is in the model, the planted genes lose their survival
association. The pipeline records this decision in the flag ledger:

``` r

sig_guarded <- build_signature(se, "time", "event", top_n = 6,
                               repeats = 3, folds = 3,
                               design_terms = "batch", seed = 12)
sig_guarded$genes
#> [1] "g35" "g15" "g31" "g47" "g13" "g23"
sig_guarded$flags
#>                check severity
#> 1 adjusted_screening     info
#>                                                     detail           stage
#> 1 Univariate screening adjusted for design term(s): batch. build_signature
```

## Tier 3: the counterfactual

|  | Naive analysis (Tier 1) | rnaSentry (Tier 2) |
|----|----|----|
| Screening | unadjusted; the 6 batch-driven genes rank as most “prognostic” | adjusted for `batch`; batch-driven genes no longer top hits |
| Design | `batch + region` fitted together; collinear terms double-counted | `redundant` flag on the pair (Cramer’s V = 0.82); minimal `~ batch` recommended |
| Surrogate association | never checked | `batch` eta-squared 0.97, `region` epsilon-squared 0.55 reported |
| Ledger | none | `redundant_variable`, `adjusted_screening` recorded |

**What a guarded analysis can claim:** without
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md),
this cohort would produce a six-gene “prognostic” signature that is
entirely a batch artifact, plus a model reporting two collinear
“independent” risk factors. With it, the batch confounder is surfaced
with its effect size, the redundant term is dropped, and screening is
adjusted so the discovered signature is not a batch echo.

## Reading

The lesson is that the recommended formula is only as good as the design
variables you supply.
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
will not invent a missing confounder; it audits what it is given. Supply
every variable that could plausibly confound expression (batch, plate,
sex, subtype), then trust the recommended formula over an ad-hoc pick.

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
#>  [1] sass_0.4.10         SparseArray_1.12.2  lattice_0.22-9     
#>  [4] magrittr_2.0.5      digest_0.6.39       RColorBrewer_1.1-3 
#>  [7] evaluate_1.0.5      grid_4.6.1          bookdown_0.47      
#> [10] fastmap_1.2.0       jsonlite_2.0.0      Matrix_1.7-5       
#> [13] BiocManager_1.30.27 scales_1.4.0        textshaping_1.0.5  
#> [16] jquerylib_0.1.4     abind_1.4-8         cli_3.6.6          
#> [19] rlang_1.3.0         XVector_0.52.0      splines_4.6.1      
#> [22] withr_3.0.3         cachem_1.1.0        DelayedArray_0.38.2
#> [25] yaml_2.3.12         otel_0.2.0          S4Arrays_1.12.0    
#> [28] tools_4.6.1         dplyr_1.2.1         ggplot2_4.0.3      
#> [31] vctrs_0.7.3         R6_2.6.1            lifecycle_1.0.5    
#> [34] fs_2.1.0            htmlwidgets_1.6.4   ragg_1.5.2         
#> [37] pkgconfig_2.0.3     desc_1.4.3          pillar_1.11.1      
#> [40] pkgdown_2.2.1.9000  bslib_0.12.0        gtable_0.3.6       
#> [43] glue_1.8.1          systemfonts_1.3.2   tidyselect_1.2.1   
#> [46] tibble_3.3.1        xfun_0.60           knitr_1.51         
#> [49] farver_2.1.2        htmltools_0.5.9     labeling_0.4.3     
#> [52] rmarkdown_2.31      compiler_4.6.1      S7_0.2.2
```
