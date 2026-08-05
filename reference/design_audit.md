# Audit candidate design covariates before differential-expression modeling

Given a `SummarizedExperiment`, `design_audit()` examines a set of
candidate design variables (technical or biological covariates such as
batch, sequencing lane or tissue) and recommends a design formula for
downstream differential-expression analysis. It does this without
fitting per-gene models: the leading principal component of the
expression data serves as a surrogate for global expression structure,
and every candidate variable is tested for association with that
surrogate. A variable that tracks the surrogate significantly is flagged
as a potential confounder, because omitting it from the design risks
attributing expression structure to the variable of interest.

## Usage

``` r
design_audit(
  se,
  design_vars,
  outcome_col = NULL,
  surrogate_pc = 1,
  alpha = 0.05,
  interaction_alpha = alpha,
  redundant_effect_size = 0.8
)
```

## Arguments

- se:

  A `SummarizedExperiment` with a count or normalized expression assay.

- design_vars:

  Character vector naming candidate design variables, all of which must
  be columns of `colData(se)`.

- outcome_col:

  Optional character. Name of the biological variable of interest (for
  example condition), also a column of `colData(se)`. When supplied it
  is placed first in the recommended formula, is never dropped as
  redundant, and is included in the pairwise redundancy scan. Defaults
  to `NULL`.

- surrogate_pc:

  Integer. Which principal component to use as the surrogate for global
  expression structure. Defaults to `1` (the leading PC); pass a
  different integer to override.

- alpha:

  Numeric in (0, 1). Significance threshold for the confounder and
  pairwise association tests. Defaults to `0.05`.

- interaction_alpha:

  Numeric in (0, 1). Significance threshold for the interaction
  likelihood-ratio test. Defaults to `alpha`.

- redundant_effect_size:

  Numeric in (0, 1). Minimum effect size (Cramer's V, absolute Pearson r
  or eta-squared) for a significant pairwise association to be declared
  redundant. Defaults to `0.8`.

## Value

An object of class `"rnaSentry_design"` (a list) with elements:

- formula:

  The recommended design formula as a
  [`stats::formula`](https://rdrr.io/r/stats/formula.html) object.

- formula_text:

  Character version of `formula`.

- rationale:

  Character vector of plain-language statements justifying every
  automated decision.

- confounder_table:

  Data.frame with one row per candidate variable (columns `variable`,
  `type`, `n`, `test`, `statistic`, `df`, `p_value`, `adj_p`,
  `effect_size`, `effect_size_type`, `flagged`). The `p_value` column
  holds the raw association p-value; `adj_p` is the Benjamini-Hochberg
  adjusted p-value across all tested candidate variables, and `flagged`
  is based on `adj_p`.

- interaction_tested:

  Logical; whether the interaction LRT could be performed (it requires
  at least two PCs and enough samples).

- interaction_table:

  Data.frame with one row per candidate variable (columns `variable`,
  `test`, `statistic`, `df`, `p_value`, `flagged`), or `NULL` when not
  tested.

- pairwise_table:

  Data.frame with one row per variable pair (columns `var1`, `var2`,
  `type_pair`, `test`, `statistic`, `df`, `p_value`, `effect_size`,
  `effect_size_type`, `redundant`), or `NULL` when fewer than two
  variables are available.

- dropped_vars:

  Character vector of design variables removed from the recommended
  formula as redundant.

- surrogate_pc:

  The surrogate PC actually used (after clamping to the available number
  of PCs).

- alpha, interaction_alpha, redundant_effect_size:

  The thresholds used.

- flags:

  Data.frame summarizing every issue raised, with columns `check`,
  `severity`, `detail` and `stage`.

## Details

Three passes are performed:

1.  *Confounder scan*: each candidate variable is tested against the
    surrogate expression gradient (`surrogate_pc`). Categorical
    variables are tested with one-way ANOVA (or Kruskal-Wallis when the
    scores are not approximately normal) and numeric variables with
    linear regression; effect sizes are reported as eta-squared,
    epsilon-squared or R-squared.

2.  *Redundancy scan*: every pair of variables (including `outcome_col`
    when supplied) is cross-tested with Cramer's V
    (categorical-categorical), Pearson correlation (numeric-numeric) or
    a one-way test (mixed), and pairs whose association is significant
    and whose effect size reaches `redundant_effect_size` are declared
    redundant.

3.  *Interaction scan*: a likelihood-ratio test compares
    `other_pc ~ surrogate + v` with `other_pc ~ surrogate * v` for each
    candidate variable, where `other_pc` is the first PC that is not the
    surrogate. A significant interaction means the variable's influence
    on expression structure changes along the surrogate gradient; this
    is reported as interaction-type confounding rather than silently
    folded into the recommended formula.

The recommended `formula` is `~ outcome_col + design_variables`, with
variables declared redundant against a later-listed variable dropped.
Every statistical decision is returned alongside a plain-language
`rationale`; nothing is dropped or included silently.

## Examples

``` r
library(SummarizedExperiment)
set.seed(4)
counts <- matrix(rpois(360, lambda = 200), nrow = 30, ncol = 12,
                  dimnames = list(paste0("gene", 1:30), paste0("S", 1:12)))
batch <- factor(rep(c("B1", "B2"), each = 6))
condition <- factor(rep(c("A", "B"), times = 6))
counts[1:10, batch == "B2"] <- counts[1:10, batch == "B2"] * 3L
coldata <- S4Vectors::DataFrame(batch = batch, condition = condition,
                                 age = runif(12, 40, 80),
                                 row.names = paste0("S", 1:12))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

design_audit(se, design_vars = c("batch", "age"), outcome_col = "condition")
#> Warning: Chi-squared approximation may be incorrect
#> Warning: Chi-squared approximation may be incorrect
#> rnaSentry design audit: 2 candidate variable(s), surrogate PC1.
#> Recommended design formula: ~ condition + batch + age
#> 1 variable(s) flagged as potential confounders.
#> Interaction LRT run; 0 variable(s) with significant surrogate interaction.
#> No issues flagged.
```
