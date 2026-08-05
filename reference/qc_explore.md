# Quality-control audit of a bulk RNA-seq SummarizedExperiment

Runs a battery of pre-analysis checks on a raw-counts
[`SummarizedExperiment`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html):
duplicate sample identifiers, missing values in sample metadata,
non-integer count values, and library-size outliers detected via a
robust (median absolute deviation-based) rule. Every check is reported
with the statistic that triggered it, not just a pass/fail flag, so the
result can be inspected or logged rather than trusted blindly.

## Usage

``` r
qc_explore(se, mad_threshold = 3)
```

## Arguments

- se:

  A `SummarizedExperiment` with a raw-count assay as its first assay
  (`assay(se, 1)`). Counts should be non-negative integers (or values
  coercible to integers without loss).

- mad_threshold:

  Numeric scalar. Samples whose log10 library size deviates from the
  median by more than `mad_threshold` robust median absolute deviations
  are flagged as library-size outliers. Defaults to `3`.

## Value

An object of class `"rnaSentry_qc"` (a list) with elements:

- n_samples:

  Number of samples (columns) in `se`.

- n_genes:

  Number of genes (rows) in `se`.

- duplicate_samples:

  Character vector of duplicated sample IDs, or `character(0)` if none.

- non_integer_counts:

  Logical; `TRUE` if the count assay contains non-integer values.

- missing_summary:

  A data.frame reporting the proportion of missing values per
  `colData(se)` column.

- library_size_outliers:

  Character vector of sample IDs flagged as library-size outliers.

- flags:

  A data.frame summarizing every issue raised, with columns `check`,
  `severity`, `detail` and `stage`.

## Examples

``` r
library(SummarizedExperiment)
set.seed(1)
counts <- matrix(rpois(60, lambda = 200), nrow = 6, ncol = 10,
                  dimnames = list(paste0("gene", 1:6), paste0("S", 1:10)))
counts[, 1] <- counts[, 1] * 50L  # inject one library-size outlier
coldata <- DataFrame(condition = rep(c("A", "B"), each = 5),
                      row.names = paste0("S", 1:10))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

qc <- qc_explore(se)
qc$flags
#>                  check severity
#> 1 library_size_outlier  warning
#>                                                   detail      stage
#> 1 Samples with library size > 3 MADs from the median: S1 qc_explore
```
