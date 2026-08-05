# Getting started: from two CSVs to a SummarizedExperiment

## Motivation

rnaSentry operates on a `SummarizedExperiment`: the expression matrix in
the assay slot and the sample metadata in `colData`. Most wet-lab and
clinical users, however, start with two plain text files:

- a **counts matrix** (genes in rows, samples in columns), and
- a **clinical table** (one row per sample, with follow-up time and
  event plus any batch, stage, or treatment covariates).

This vignette shows the exact bridge between those two files and the
pipeline: read the files, keep only the samples present in both (the
single most common pitfall), validate the intake with
[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md),
audit the design with
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md),
and run the full pipeline.

The example below writes two small synthetic files to a temporary
directory as stand-ins for *your* `counts.csv` and `clinical.csv`;
replace the file paths and column names with your own and the block
works unchanged.

## Build two example CSV files

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

set.seed(42)
n_genes <- 300
n_samples <- 30
genes <- paste0("GENE", seq_len(n_genes))
samples <- paste0("S", seq_len(n_samples))

# Raw integer counts (bulk RNA-seq style), genes x samples.
cnt <- matrix(rpois(n_genes * n_samples, lambda = 200),
              nrow = n_genes, ncol = n_samples,
              dimnames = list(genes, samples))

# Clinical metadata, one row per sample.
clin <- data.frame(
    time  = round(rexp(n_samples, rate = 0.08), 1),
    event = rbinom(n_samples, size = 1, prob = 0.5),
    batch = factor(rep(c("B1", "B2", "B3"), length.out = n_samples)),
    stage = factor(sample(c("I", "II", "III"), n_samples, replace = TRUE)),
    row.names = samples
)

# --- Deliberately introduce the two classic file problems -------------------
# A sample present in the counts file but not the clinical file,
cnt <- cbind(cnt, S31 = cnt[, "S2"])
# and a clinical row with no matching counts column.
clin <- rbind(clin, S32 = clin["S1", ])
# Plus two intake issues load_counts() should catch: a non-integer count...
cnt["GENE1", "S1"] <- cnt["GENE1", "S1"] + 0.5
# ...and a gene with no counts anywhere.
cnt["GENE300", ] <- 0L

td <- tempfile("csvdemo")
dir.create(td)
write.csv(cnt, file.path(td, "counts.csv"))
write.csv(clin, file.path(td, "clinical.csv"))
list.files(td)
#> [1] "clinical.csv" "counts.csv"
```

Every reader sees the same two files a real user exports from their
pipeline: `counts.csv` (genes in rows, samples in columns) and
`clinical.csv` (one row per sample).

## Step 1: read both files

``` r

counts   <- as.matrix(read.csv(file.path(td, "counts.csv"),   row.names = 1))
clinical <- read.csv(file.path(td, "clinical.csv"), row.names = 1)

dim(counts)
#> [1] 300  31
dim(clinical)
#> [1] 31  4
```

Note the [`as.matrix()`](https://rdrr.io/r/base/matrix.html):
[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
expects a matrix, and
[`read.csv()`](https://rdrr.io/r/utils/read.table.html) returns a
`data.frame`.

## Step 2: intersect the sample IDs

The most common pitfall is sample IDs that do not match between the two
files: one file may contain samples the other lacks (a failed RNA
extraction, a sample re-labeled, or a version mismatch). Passing them to
[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
un-checked produces a hard error or, worse, silently mis-aligned
columns. The fix is a single
[`intersect()`](https://generics.r-lib.org/reference/setops.html) that
keeps only samples found in **both** files:

``` r

n_counts_before <- ncol(counts)
n_clin_before   <- nrow(clinical)

common <- intersect(colnames(counts), rownames(clinical))
counts   <- counts[,  common, drop = FALSE]
clinical <- clinical[common, ,   drop = FALSE]

data.frame(
    "counts columns" = c(before = n_counts_before, after = ncol(counts)),
    "clinical rows"  = c(before = n_clin_before,   after = nrow(clinical)),
    row.names = c("before intersect", "after intersect")
)
#>                  counts.columns clinical.rows
#> before intersect             31            31
#> after intersect              30            30
```

This example planted one sample in each file that the other lacks, so
the 30 truly shared samples survive and `S31`/`S32` are dropped. Always
check the dimensions after intersecting: if you lose hundreds of
samples, suspect a version or delimiter mismatch between the two files.

## Step 3: build and validate the SummarizedExperiment

[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
constructs the object and records every intake issue it finds in a flag
ledger stored in `metadata(se)$flags`:

``` r

se <- load_counts(counts, clinical)
#> Warning: Count matrix contains non-integer values. Downstream stages expect
#> integer counts.
#> Warning: 1 gene(s) have all-zero or all-NA counts and will be retained but may
#> be filtered later.
se
#> class: SummarizedExperiment 
#> dim: 300 30 
#> metadata(1): flags
#> assays(1): counts
#> rownames(300): GENE1 GENE2 ... GENE299 GENE300
#> rowData names(0):
#> colnames(30): S1 S2 ... S29 S30
#> colData names(4): time event batch stage
```

``` r

S4Vectors::metadata(se)$flags
#>                check severity                                    detail
#> 1 non_integer_counts  warning Count matrix contains non-integer values.
#> 2         empty_rows  warning 1 gene(s) have all-zero or all-NA counts.
#>         stage
#> 1 load_counts
#> 2 load_counts
```

The ledger flags the two problems we planted (the non-integer count and
the all-zero gene) and carries them into the rendered report’s audit
trail. Duplicate sample IDs, missing metadata, and Ensembl-style gene
IDs are caught the same way. The pipeline still runs on a flagged object
– the ledger’s job is to make sure you saw the issue first.

## Step 4: audit the design before modeling

[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
tests every candidate covariate against the leading expression principal
component and flags potential confounders, before any per-gene modeling
is invested in:

``` r

da <- design_audit(se, design_vars = c("batch", "stage"))
#> Warning in stats::chisq.test(tab): Chi-squared approximation may be incorrect
#> Warning in stats::chisq.test(tab): Chi-squared approximation may be incorrect
da
#> rnaSentry design audit: 2 candidate variable(s), surrogate PC1.
#> Recommended design formula: ~ batch + stage
#> 0 variable(s) flagged as potential confounders.
#> Interaction LRT run; 0 variable(s) with significant surrogate interaction.
#> No issues flagged.
```

## Step 5: run the full pipeline

[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
orchestrates design audit, PCA audit, signature construction with
repeated cross-validation, and the survival-modeling stages on the same
object. It needs a numeric follow-up `time` column and a `0/1` `event`
indicator in `colData`:

``` r

run <- run_rnaSentry(se, time_col = "time", event_col = "event",
                     top_n = 5, repeats = 1, folds = 2, seed = 42,
                     design_vars = c("batch", "stage"),
                     render_report = FALSE)
#> Warning in stats::chisq.test(tab): Chi-squared approximation may be incorrect
#> Warning in stats::chisq.test(tab): Chi-squared approximation may be incorrect
#> Warning: Only 14 event(s) for a 5-gene signature (2.8 events per parameter),
#> below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard
#> ratios are unstable at this event count; consider fewer genes or a larger
#> cohort.
run
#> rnaSentry pipeline run.
#>   - design_audit
#>   - pca_audit
#>   - build_signature
#>   - km_curve
#>   - cox_model
#>   - survival_parametric
```

On 30 samples the `events_per_parameter` guardrail warns that
cross-validated concordance is unstable – exactly as it should: this
cohort is a demonstration, not a discovery set. See
[`vignette("case-study-small-cohort")`](https://ghrieb.github.io/rnaSentry/articles/case-study-small-cohort.md)
for what that warning means and
[`?generate_report`](https://ghrieb.github.io/rnaSentry/reference/generate_report.md)
for how to render the audited HTML report of a run.

The run locked the discovered signature:
[`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
is single-use per session, so a second call in the same session is
refused until the lock is released. This is the reproducibility
guardrail against silently re-selecting signature genes after survival
analysis.

``` r

run$stages$build_signature <-
    lock_signature(run$stages$build_signature, lock = FALSE)
```

## Assay naming: `counts` vs `logcounts`

[`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
stores your matrix under the assay name `counts`, and the pipeline
computes `log2(counts + 1)` when it runs. If your input file is *already
normalized* (CPM, TPM, VST, or a microarray matrix), store it under the
name `logcounts` instead so it is used as-is and not double-logged:

``` r

se_norm <- se
SummarizedExperiment::assayNames(se_norm) <- "logcounts"
```

Either way the bridge is identical; only the assay name tells the
pipeline which scale your data is in.

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
