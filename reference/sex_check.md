# Verify reported sex against inferred biological sex from expression

Compares the sex recorded in sample metadata against sex inferred from
expression of *XIST* and a panel of Y-chromosome marker genes. Designed
to catch sample swaps and metadata entry errors before they propagate
into downstream differential expression or survival analysis.

## Usage

``` r
sex_check(
  se,
  sex_col = "sex",
  xist_gene = "XIST",
  y_genes = c("RPS4Y1", "DDX3Y", "KDM5D")
)
```

## Arguments

- se:

  A `SummarizedExperiment` with gene symbols as `rownames(se)` and a
  numeric (counts or normalized) assay as its first assay.

- sex_col:

  Character. Name of the `colData(se)` column holding reported sex,
  coded as `"M"` / `"F"`. Defaults to `"sex"`.

- xist_gene:

  Character. Gene symbol for XIST. Defaults to `"XIST"`.

- y_genes:

  Character vector of Y-chromosome marker gene symbols. Defaults to
  `c("RPS4Y1", "DDX3Y", "KDM5D")`.

## Value

A data.frame with one row per sample and columns `sample_id`,
`reported_sex`, `inferred_sex`, and `status` (one of `"OK"`,
`"MISMATCH"`, `"AMBIGUOUS"`, or `"MISSING"` when reported sex is
absent).

## Details

Inference rule: a sample is called `"F"` when XIST expression is high
relative to the Y-gene panel, `"M"` when the reverse holds, and
`"ambiguous"` when neither signal is clearly dominant (for example,
low-quality or highly degraded samples). The decision uses the sign of
each sample's XIST-minus-Y-panel rank score: a score at least one rank
unit above zero is called `"F"`, at least one rank unit below zero
`"M"`, and a score within one rank unit of zero `"ambiguous"`. Because
the rule is based on ranks rather than hard-coded expression cutoffs, it
is robust across different normalizations and sequencing depths.

## Examples

``` r
library(SummarizedExperiment)
genes <- c("XIST", "RPS4Y1", "DDX3Y", "KDM5D", "GAPDH")
samples <- paste0("S", 1:6)
expr <- matrix(0, nrow = length(genes), ncol = length(samples),
                dimnames = list(genes, samples))
# samples 1-3 look female (high XIST, low Y genes)
expr["XIST", 1:3] <- c(950, 900, 875)
expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 1:3] <- 5
# samples 4-6 look male, but S6's metadata will be mislabeled below
expr["XIST", 4:6] <- 5
expr[c("RPS4Y1", "DDX3Y", "KDM5D"), 4:6] <- 400
expr["GAPDH", ] <- 1000

coldata <- DataFrame(sex = c("F", "F", "F", "M", "M", "F"),
                      row.names = samples)  # S6 mislabeled as F
se <- SummarizedExperiment(assays = list(counts = expr), colData = coldata)

sex_check(se)
#>   sample_id reported_sex inferred_sex   status
#> 1        S1            F            F       OK
#> 2        S2            F            F       OK
#> 3        S3            F            F       OK
#> 4        S4            M            M       OK
#> 5        S5            M            M       OK
#> 6        S6            F            M MISMATCH
```
