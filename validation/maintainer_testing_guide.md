# rnaSentry — maintainer testing guide

How to reproduce every validation gate for the rnaSentry Bioconductor
submission, in order. Each gate has a script, an expected result, and the
committed log that records the passing run.

Environment used by the gate: Windows 11, R 4.5.2 at
`C:\Program Files\R\R-4.5.2`, pandoc via RStudio at
`C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools` (point
`RSTUDIO_PANDOC` at that **directory**, not the exe). The R session does not
see that variable from a plain `Rscript -e`; it is honored by `R CMD build` /
`R CMD check` when exported in the same shell. There is **no Rtools** on this
machine, so `devtools::check()` aborts — use `R CMD build` + `R CMD check` on
the tarball instead. PowerShell mangles `$` in `Rscript -e` strings — put
multi-statement R into a temp script under
`C:\Users\Hani\AppData\Local\Temp\opencode\`.

## Gate 0 — hygiene and submission basics

```powershell
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
R CMD build rnaSentry
R CMD check rnaSentry_0.99.0.tar.gz --as-cran --output=.\chk
Rscript -e "BiocCheck::BiocCheck('.')"
Rscript -e "BiocCheck::BiocCheckGitClone()"   # run on a tree with the
                                              # rnaSentry.BiocCheck stamp removed
```

Expected: `R CMD check --as-cran` → 0 ERROR, 1 WARNING (`qpdf`), 1 NOTE
(`tidy`), both external tools; the previous run without `--no-build-vignettes`
also verifies the compiled vignette (`inst/doc`) and `browseVignettes`.
BiocCheck: 1 environmental ERROR (support-site email 404 — register the
maintainer email on https://support.bioconductor.org), 1 justified WARNING
(`set.seed` in the documented `seed` argument), 12 advisory NOTES. Every item
is explained in `logs/notes_documented.md`. Logs: `logs/00.2_check.txt`,
`00.3_bioccheck.txt`, `00.3_bioccheck_gitclone.txt`, `00.5_bioccheck_tarball.txt`,
`00.4_as_cran.txt`, `00.6_as_cran.txt`, `05_as_cran.txt`, `06_as_cran.txt`,
`07_as_cran.txt` (post sign-convention fix; 06 re-run after the
`NEWS`/man-page updates, 07 re-run after the cox_model concordance parity test —
both 0 ERROR / 1 WARNING / 1 NOTE).

## Gate 1 — full test suite

```powershell
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
Rscript <temp>/run_tests_parity.R   # devtools::test() with load_all
```

Expected: **120 files / 393 passed / 0 failed / 0 error / 0 warning**.
The suite must be run via `devtools::test()` (loading environment differs from
a plain `test_dir`). Log: `logs/00.1_test.txt`.

## Gate 2 — adversarial review

The `test-audit.md` document lists every check that failed or passed only after
deliberate fixes (seeds, fold stratification, edge cases like `n_pcs == 1`,
non-syntactic gene symbols). After any change to `R/*.R`, re-run Gate 1 and
`R CMD check --as-cran` before proceeding.

## Gate 3 — statistical parity

`tests/testthat/test-statistical_parity.R` re-implements every statistic in the
package with an independent reference implementation and asserts equality:

- `km_curve` / `validate_external` log-rank ↔ `survival::survdiff`
- `design_audit` Cramér's V ↔ `sqrt(chisq/(n*(min(r,c)-1)))`
- `cox_model` HR / CI / p ↔ independent `survival::coxph` + `stats::confint`
- `cox_model` concordance ↔ `survival::concordance(Surv ~ lp, reverse = TRUE)`
  on the fitted linear predictor (independently recomputes the direction that
  `summary(fit)$concordance` reports)
- Schoenfeld rho / p ↔ `survival::cox.zph`
- `survival_parametric` AIC / loglik / npar ↔ `stats::AIC` / `logLik`
- `build_signature` fold concordance ↔ `survival::concordance` on replayed
  event-stratified seeded folds, all with `reverse = TRUE` (the risk-score
  convention). Direction guards: on a signal-bearing fixture the CV mean must
  beat 0.5, and `validate_external` concordance must beat 0.5 on a cohort
  whose survival is engineered to follow the signature score. A dedicated
  invariant asserts `concordance(Surv ~ x) + concordance(Surv ~ x,
  reverse = TRUE) = 1`. These tests were added after a reviewer caught a
  sign-convention bug in `build_signature()` and `validate_external()` that
  had reported `1 − Harrell's C` (see `face_validity_review.md` and
  `test_audit.md`).

## Gate 4 — simulation study (11 falsifiable targets)

```r
source("validation/simulate_study.R")   # exits non-zero on any failure
```

Covers planted sex-label swaps, planted confounder redundancy, the null-data
CV-concordance expectation (with a random-gene calibration control), a
planted time-varying-hazard scenario, and the external-validation lock check.
Expected: **11 targets / 11 PASS**. Results and interpretation:
`simulation_study.md`. Log: `logs/03_simulation.txt`.

## Gate 5 — real-data face validity (GSE20685)

The repro script downloads GSE20685 via GEOquery, maps GPL570 probes to symbols,
collapses to a gene × sample matrix, runs the full pipeline, and extracts the
four-point evidence:

1. high-risk median OS finite (11.4 y) vs low-risk median not reached; log-rank
   p = 3.8e-17;
2. MKI67 HR = 1.21, ESR1 HR = 0.86 (p = 0.0017) — literature directions;
3. `design_audit` flags subtype (eta² = 0.76) — biologically plausible;
4. ≥ 1 limitation warning in the report (`proportional_hazards`).

Also documents the out-of-sample result: held-out CV concordance of the
selected signature is **0.783** (in-sample 0.840, random-gene control 0.582 on
real data / 0.490 on simulated nulls). An earlier reading of 0.217 was the
symptom of the sign-convention bug described under Gate 3, now fixed. Known
limitation: the discovery cohort has 83 events for a 20-gene signature (~4.2
events/parameter, below the ~10 rule of thumb) — see the "Known limitation"
note in `face_validity_review.md`. Log: `logs/04_face_validity.txt`, report:
`logs/face_validity_report.html`.

## Gate 6 — cross-platform (not yet run on this machine)

Needs a Docker-capable machine (`validation/docker_check.sh` runs the devel
container: build + BiocCheck + `--as-cran` + tests), or rhub/win-builder
(needs a GitHub token / an email address). Status and instructions:
`cross_platform.md`. When results arrive, append the per-platform
`Status:` lines to that file and to the submission notes.

## Full regression loop after any code change

1. `R CMD build` → `R CMD check --as-cran` on the tarball (0 ERROR; the only
   WARNING/NOTE are the `qpdf`/`tidy` external tools).
2. Gate 1 suite (120/393).
3. If statistics/tests changed: Gate 3 parity suite, Gate 4 simulation.
4. If feature/annotation handling changed: Gate 5 GSE20685 repro.
5. Commit logs with the change at a logical checkpoint.
