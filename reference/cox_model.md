# Adjusted Cox proportional-hazards model for a signature

`cox_model()` fits a single multivariate Cox model of the survival
outcome on the signature genes together with any design terms recorded
in the signature at build time (`sig$design_terms`) and any extra
`confounders` supplied here. It is the companion to the gene-only
scoring model in
[`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md):
the signature's risk score deliberately ignores the design covariates so
it stays portable, while this stage supplies inference *adjusted* for
them on a cohort that records them.

## Usage

``` r
cox_model(sig, se, confounders = character(0))
```

## Arguments

- sig:

  An object of class `"rnaSentry_signature"` as returned by
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md).

- se:

  A `SummarizedExperiment` containing the signature genes in its rows
  and the survival metadata columns recorded in `sig` (`sig$time_col`
  and `sig$event_col`). Expression is taken from the same analysis assay
  used at signature-build time (see
  [`build_signature`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)).

- confounders:

  Character vector of additional `colData(se)` columns to adjust for,
  beyond `sig$design_terms`. Defaults to `character(0)`.

## Value

An object of class `"rnaSentry_cox_model"` (a list) with elements:

- genes:

  The signature genes in the model.

- terms:

  All fitted terms: the genes, then `sig$design_terms`, then
  `confounders`.

- formula:

  The model formula.

- fit:

  The `coxph` fit object.

- coef_table:

  Data.frame with one row per fitted term: `term`, `coefficient`, `se`,
  `HR`, `HR_low`, `HR_high` and `p`.

- gene_summary:

  Data.frame of the signature-gene rows of `coef_table` with an added
  BH-adjusted `adj_p`.

- zph_summary:

  Data.frame from `cox.zph` with one row per term plus `GLOBAL`; `NULL`
  when the test could not be run.

- concordance:

  Harrell's concordance index of the fitted model.

- ph_violated:

  Logical; `TRUE` when any term (including `GLOBAL`) violates
  proportional hazards at alpha 0.05.

- n, events:

  Number of samples and events used by the fit.

- sig_locked:

  Whether the input signature was locked.

- created:

  Provenance timestamp.

- flags:

  Data.frame of issues raised (missing covariate values, non-estimable
  genes, proportional-hazards violations), with columns `check`,
  `severity`, `detail` and `stage`.

## Details

For every fitted term the model reports the hazard ratio with a Wald
confidence interval and p-value. For the signature genes it also reports
a Benjamini-Hochberg adjusted p-value across the genes, so that
gene-level multiplicity is controlled. Proportional-hazards diagnostics
([`cox.zph`](https://rdrr.io/pkg/survival/man/cox.zph.html)) are run on
every term and the GLOBAL test; a violation flags the model.

The stage is read-only: it works on both locked and unlocked signatures
and records the signature's lock state in its output.

## Assumptions and limitations

The model assumes proportional hazards (tested and flagged via
`cox.zph`) and standard right-censored survival. It does not model
competing risks, time-varying covariates, or left truncation. Like the
rest of the pipeline it makes no batch-effect correction; adjustment for
recorded design covariates is supplied through `sig$design_terms` and
`confounders`.

## Examples

``` r
library(SummarizedExperiment)
set.seed(7)
counts <- matrix(rpois(450, lambda = 500), nrow = 30, ncol = 15,
                 dimnames = list(paste0("gene", 1:30), paste0("S", 1:15)))
sig_expr <- colMeans(counts[1:4, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(15, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(15, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 batch = factor(rep(c("B1", "B2"),
                                                    length.out = 15)),
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)
sig <- build_signature(se, "time", "event", top_n = 4, repeats = 1,
                       folds = 2, seed = 1, design_terms = "batch")
#> Warning: Only 10 event(s) for a 4-gene signature (2.5 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
cm <- cox_model(sig, se)
cm
#> rnaSentry Cox model: 4 gene(s) + 1 covariate term(s), 10 events.
#> Concordance: 0.852. Proportional hazards: violated (see flags).
#> Signature genes (HR [95% CI], BH adj. p):
#>   gene14         35.43 [0.00, 27783335525367.37]  p = 0.79849  adj.p = 0.79849
#>   gene7          0.00 [0.00, 6.34]  p = 0.070145  adj.p = 0.28058
#>   gene3          20828683.23 [0.00, 439369724642329664.00]  p = 0.16471  adj.p = 0.32943
#>   gene22         1768.75 [0.00, 20237523690717.25]  p = 0.52685  adj.p = 0.70246
#> Signature not locked.
#> 1 issue(s) flagged:
#>   [warning] proportional_hazards: Proportional-hazards assumption may be violated (term(s): batch).
plot(cm)

```
