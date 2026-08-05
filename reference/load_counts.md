# Validate and construct a SummarizedExperiment from count data

Thin validation wrapper around
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
that catches the most common input mistakes before they propagate into
the pipeline: duplicate or missing sample IDs, non-integer counts, empty
rows or columns, and non-syntactic gene identifiers.

## Usage

``` r
load_counts(counts, colData, assay_name = "counts", check_gene_ids = TRUE)
```

## Arguments

- counts:

  Integer or numeric matrix of gene expression values (genes in rows,
  samples in columns). Row names must be gene identifiers; column names
  must be sample identifiers.

- colData:

  A `data.frame` or
  [`S4Vectors::DataFrame`](https://rdrr.io/pkg/S4Vectors/man/DataFrame-class.html)
  of sample-level metadata. Row names must match the column names of
  `counts`. When omitted, an empty `DataFrame` is created.

- assay_name:

  Character. Name for the assay slot. Defaults to `"counts"`.

- check_gene_ids:

  Logical. When `TRUE` (default) the function warns if gene identifiers
  look like Ensembl accessions (starting with `"ENS"`) or contain
  non-syntactic characters, because downstream stages rely on
  gene-symbol rownames.

## Value

A `SummarizedExperiment` with one assay named `assay_name`. Any issues
found are recorded in a flag ledger accessible via
`S4Vectors::metadata(se)$flags` and surfaced via
[`warning()`](https://rdrr.io/r/base/warning.html).

## Assumptions and limitations

This function validates structure; it does not normalise, filter, or
transform the count matrix. Gene identifiers are checked but not
converted. Sample-metadata columns used later by the pipeline (e.g.
`time_col`, `event_col`, `batch_col`) are not validated here — that
happens in the stage that first uses them.

## Examples

``` r
counts <- matrix(1:12, nrow = 3,
                 dimnames = list(c("TP53", "BRCA1", "MYC"),
                                 c("S1", "S2", "S3", "S4")))
coldata <- data.frame(
  time = c(10, 12, 8, 15),
  event = c(1, 0, 1, 0),
  row.names = colnames(counts))
se <- load_counts(counts, coldata)
se
#> class: SummarizedExperiment 
#> dim: 3 4 
#> metadata(1): flags
#> assays(1): counts
#> rownames(3): TP53 BRCA1 MYC
#> rowData names(0):
#> colnames(4): S1 S2 S3 S4
#> colData names(2): time event
```
