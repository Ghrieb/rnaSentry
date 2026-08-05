# PCA audit of expression data with optional batch association testing

Runs principal component analysis on a variance-stabilized or
log-transformed version of the count assay and reports how much variance
each principal component (PC) explains. When `batch_col` is supplied,
each PC score is tested against the batch variable and any PC where
batch explains a significant share of variance is flagged. This feeds
[`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md)
by surfacing which batch/technical variables are associated with global
expression structure.

## Usage

``` r
pca_audit(
  se,
  batch_col = NULL,
  top_n_pcs = 5,
  scale = TRUE,
  batch_alpha = 0.05
)
```

## Arguments

- se:

  A `SummarizedExperiment` with a count or normalized expression assay.

- batch_col:

  Optional character. Name of a `colData(se)` column (categorical or
  numeric) to test against each PC score. Defaults to `NULL` (no batch
  testing).

- top_n_pcs:

  Integer. Number of leading PCs to include in `scores` and to test
  against `batch_col`. Defaults to `5`.

- scale:

  Logical. Whether to scale each gene to unit variance before PCA
  (`scale. = TRUE`). Defaults to `TRUE` so that genes are comparable on
  the log-expression scale.

- batch_alpha:

  Numeric in (0, 1). Significance threshold used to flag a PC whose
  score associates with `batch_col`. Defaults to `0.05`.

## Value

An object of class `"rnaSentry_pca"` (a list) with elements:

- pca:

  The [`stats::prcomp`](https://rdrr.io/r/stats/prcomp.html) object.

- percent_variance:

  Named numeric vector of the percent of variance explained by every PC.

- scores:

  A data.frame of the first `top_n_pcs` PC scores, with sample IDs as
  row names.

- col_data:

  The `colData(se)` columns as a data.frame, kept so
  [`plot_pca_audit`](https://ghrieb.github.io/rnaSentry/reference/plot_pca_audit.md)
  can color points by a metadata column.

- assay_used:

  Name of the assay used for the analysis.

- n_genes_analyzed:

  Number of genes retained after filtering.

- n_genes_filtered:

  Number of genes removed before PCA.

- filtered_reason:

  Character describing why genes were removed (empty string if none).

- batch_col:

  The batch column name used, or `NULL`.

- batch_tests:

  A data.frame with one row per tested PC (columns `pc`, `test`,
  `statistic`, `df`, `p_value`, `adj_p`, `effect_size`,
  `effect_size_type`, `flagged`), or `NULL` when `batch_col` is not
  supplied.

- flags:

  A data.frame summarizing every issue raised, with columns `check`,
  `severity`, `detail` and `stage`.

## Details

The analysis matrix is chosen by preference from an existing assay named
`"logcounts"` or `"vst"`; if neither exists, `log2(counts + 1)` is
computed from the first assay. Genes with zero variance or missing
values are removed because they carry no information for PCA and would
make [`stats::prcomp`](https://rdrr.io/r/stats/prcomp.html) fail; the
number removed and the reason are reported in the result (nothing is
dropped silently).

Unlike
[`design_audit`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md),
no multiple-testing correction is applied to the per-PC batch tests:
`batch_alpha` is applied to each test at its nominal level. This is a
deliberate choice, because at most `top_n_pcs` (default 5) PCs are
tested against a single batch variable, so the correction would be a
near no-op; the raw p-values and a Benjamini-Hochberg adjusted p-value
(`adj_p`) are both reported in `batch_tests` so users can apply their
own threshold. Candidate confounder screening in
[`design_audit`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md),
which tests many variables at once, is BH-corrected instead.

## Examples

``` r
library(SummarizedExperiment)
set.seed(1)
counts <- matrix(rpois(240, lambda = 200), nrow = 20, ncol = 12,
                  dimnames = list(paste0("gene", 1:20), paste0("S", 1:12)))
coldata <- DataFrame(condition = rep(c("A", "B"), each = 6),
                      batch = factor(rep(c("B1", "B2"), each = 6)),
                      row.names = paste0("S", 1:12))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

pca_audit(se, batch_col = "batch")
#> rnaSentry PCA audit: 12 samples, 20 genes analyzed (assay: log2(counts+1))
#> PC1 explains 21.3% of variance.
#> No issues flagged.
```
