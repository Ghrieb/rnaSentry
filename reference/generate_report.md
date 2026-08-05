# Generate the rnaSentry analysis report

`generate_report()` assembles the results of any collection of rnaSentry
stages into a single self-contained HTML report. It renders the plain
`report_template.Rmd` shipped with the package (via
[`render`](https://pkgs.rstudio.com/rmarkdown/reference/render.html)),
which prints a pipeline summary, a unified audit-trail table collecting
every flag from every supplied stage (all rnaSentry stage results carry
a `flags` data.frame with `check`/`severity`/`detail`/`stage` columns),
and one section per stage using its `print` method.

## Usage

``` r
generate_report(
  stages,
  output_file = "rnaSentry_report.html",
  output_dir = ".",
  quiet = TRUE
)
```

## Arguments

- stages:

  A named list of rnaSentry stage results (for example elements of the
  `stages` list returned by
  [`run_rnaSentry`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  or a hand-assembled set of audit, signature, survival-modeling and
  validation outputs). Names are used as report section titles.

- output_file:

  Character. File name (without directory) of the rendered report.
  Defaults to `"rnaSentry_report.html"`.

- output_dir:

  Character. Directory to write the report into. Must exist. Defaults to
  the working directory.

- quiet:

  Logical. Passed to
  [`rmarkdown::render()`](https://pkgs.rstudio.com/rmarkdown/reference/render.html);
  suppresses the rendering console output when `TRUE` (default).

## Value

The path of the rendered report, invisibly.

## Examples

``` r
library(SummarizedExperiment)
set.seed(9)
counts <- matrix(rpois(400, lambda = 500), nrow = 20, ncol = 20,
                  dimnames = list(paste0("gene", 1:20), paste0("S", 1:20)))
sig_expr <- colMeans(counts[1:5, , drop = FALSE])
risk <- scale(sig_expr)[, 1] * 0.4
event_time <- rexp(20, rate = 0.03 * exp(0.8 * risk))
censor_time <- rexp(20, rate = 0.02)
time <- pmin(event_time, censor_time)
event <- as.integer(event_time < censor_time)
coldata <- S4Vectors::DataFrame(time = time, event = event,
                                 row.names = colnames(counts))
se <- SummarizedExperiment(assays = list(counts = counts), colData = coldata)

sig <- build_signature(se, "time", "event", top_n = 5, repeats = 1,
                        folds = 2, seed = 1)
#> Warning: Only 11 event(s) for a 5-gene signature (2.2 events per parameter), below 'min_events_per_parameter' = 5. Cross-validated concordance and hazard ratios are unstable at this event count; consider fewer genes or a larger cohort.
stages <- list(
  km_curve = km_curve(sig, se),
  cox_model = cox_model(sig, se)
)
generate_report(stages, output_dir = tempdir())
```
