# rnaSentry — maintainer testing guide

How to reproduce every validation gate for the rnaSentry Bioconductor
submission, in order. Each gate has a script, an expected result, and the
committed log that records the passing run.

Environment used by the gate: Windows 11, R 4.5.2 at
`C:\Program Files\R\R-4.5.2`, pandoc via RStudio at
`C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools` (point
`RSTUDIO_PANDOC` at that **directory**, not the exe). The R session does not
see that variable from a plain `Rscript -e`; it is honored by `R CMD build` /
`R CMD check` when exported in the same shell. As of 2026-08-05 the qpdf CLI
12.3.2 (MSVC64) is installed locally for the PDF size-reduction check; export
`R_QPDF` to its `qpdf.exe` and set `_R_CHECK_DOC_SIZES_=true` (otherwise
`--as-cran` emits the unconditional `'qpdf' is needed` WARNING). There is
**no Rtools** on this machine, so `devtools::check()` aborts — use `R CMD build` + `R CMD check` on
the tarball instead. PowerShell mangles `$` in `Rscript -e` strings — put
multi-statement R into a temp script under
`C:\Users\Hani\AppData\Local\Temp\opencode\`.

## Gate 0 — hygiene and submission basics

```powershell
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
$env:R_QPDF = "$env:LOCALAPPDATA\Programs\qpdf\qpdf-12.3.2-msvc64\bin\qpdf.exe"
$env:_R_CHECK_DOC_SIZES_ = "true"
R CMD build rnaSentry
$env:_R_CHECK_CRAN_INCOMING_ = "false"   # no reliable CRAN connectivity on this machine
R CMD check rnaSentry_0.99.0.tar.gz --as-cran --output=.\chk
Rscript -e "BiocCheck::BiocCheck('.')"
Rscript -e "BiocCheck::BiocCheckGitClone()"   # run on a tree with the
                                              # rnaSentry.BiocCheck stamp removed
```

Expected: `R CMD check --as-cran` → 0 ERROR, 0 WARNING, 1 NOTE
(`tidy`; the HTML-validator binary is deliberately not installed — see
`logs/notes_documented.md`). The `qpdf` WARNING of prior gates cleared after
the 2026-08-05 qpdf CLI install (`R_QPDF` + `_R_CHECK_DOC_SIZES_=true`). The
previous run without `--no-build-vignettes`
also verifies the compiled vignette (`inst/doc`) and `browseVignettes`.
BiocCheck: 1 environmental ERROR (support-site email 404 — register the
maintainer email on https://support.bioconductor.org), 0 WARNINGs, 8 advisory
NOTEs at the `16_bioccheck.txt`/`16_bioccheck_tarball.txt` run (2026-08-05
Ultimate Pre-Flight; the `set.seed` WARNING and the R-version / `Avoid 1:` /
`\dontrun` NOTEs were resolved in the same round). The Phase-1 BiocParallel
round (2026-08-05) re-ran the whole gate unchanged: `17_as_cran.txt` is
0 ERROR / 0 WARNING / 1 NOTE (`tidy` only) and
`17_bioccheck.txt`/`17_bioccheck_tarball.txt` are 1 ERROR (support-site
email 404, environmental) / 0 WARNING / 8 advisory NOTES; the GitClone
re-run `17_bioccheck_gitclone.txt` is 0 ERROR / 1 WARNING (CITATION `doi`) /
0 NOTES. Every item
is explained in `logs/notes_documented.md`. Logs: `logs/00.2_check.txt`,
`00.3_bioccheck.txt`, `00.3_bioccheck_gitclone.txt`, `00.5_bioccheck_tarball.txt`,
`00.4_as_cran.txt`, `00.6_as_cran.txt`, `05_as_cran.txt`, `06_as_cran.txt`,
`07_as_cran.txt` (post sign-convention fix; 06 re-run after the
`NEWS`/man-page updates, 07 re-run after the cox_model concordance parity test —
  both 0 ERROR / 1 WARNING / 1 NOTE),   `08_as_cran.txt` (2026-08-04, after the
  guardrail + README/vignette/report-assumptions round — 0 ERROR / 1 WARNING /
  1 NOTE, both environmental), `09_as_cran.txt` (2026-08-04, after the
  `load_counts()` + enforceable-lock round — 0 ERROR / 1 WARNING / 2 NOTEs
  because of a one-off `unable to verify current time` note),
  `10_as_cran.txt` and `11_as_cran.txt` (2026-08-04 — back to 0 ERROR /
  1 WARNING / 1 NOTE, both environmental; log `11_live_geo.txt` is the Phase-3
  live GEO run, see Gate 5b), and `12_as_cran.txt` (2026-08-04, after the
  CV-optimism reframe doc round — clean baseline again), and
  `13_as_cran.txt` (2026-08-04, Phase-3 case-study/power round: three new
  vignettes + `inst/extdata/gse20685_case_study.rds` — 0 ERROR / 1 WARNING /
  1 NOTE, both environmental), and `14_as_cran.txt` (2026-08-04, Batch-A
  ship-blocking round: LICENSE holder, report-template `results='asis'`,
  S4Vectors → Imports, xz-recompressed extdata — 0 ERROR / 1 WARNING /
  1 NOTE, both environmental), and `15_as_cran.txt` (2026-08-04, Batch C
   design-hardening round: targeted lock unlock, `min_events_per_parameter`
   forwarding, `concordance_na` flag, flags `stage` column docs — 0 ERROR /
   1 WARNING / 1 NOTE, both environmental), and `16_as_cran.txt` (2026-08-05,
   Ultimate Pre-Flight: `withr::with_seed()` seeding refactor, `seq_len()`
   fix, self-contained `\donttest` examples, qpdf CLI installed — **0 ERROR /
   0 WARNING / 1 NOTE**, the sole NOTE being the deliberately-skipped `tidy`).
   Note: prior gates run `R CMD check --as-cran` with
   `_R_CHECK_CRAN_INCOMING_=false` (this machine has no reliable CRAN
   connectivity; the remote/incoming block is skipped, so no examples-timing
   or pandoc NOTEs appear — the documented baseline was 1 WARNING / 1 NOTE
   before the 2026-08-05 qpdf install and is now 0 WARNING / 1 NOTE).

## Gate 1 — full test suite

```powershell
$env:RSTUDIO_PANDOC = "C:\Program Files\RStudio\resources\app\bin\quarto\bin\tools"
Rscript <temp>/run_tests_parity.R   # devtools::test() with load_all
```

Expected: **147 blocks / 457 passed / 0 failed / 0 error / 32 warnings**.
The 32 warnings are the `events_per_parameter` guardrail firing on
deliberately small synthetic fixtures used by tests that exercise other
behavior; each guardrail has a dedicated test that asserts its own firing
(see Gate 3 and `test_audit.md`). The suite must be run via
`devtools::test()` (loading environment differs from a plain `test_dir`).
Logs: `logs/00.1_test.txt`; the suite grew to 147/457 with the 2026-08-05
BiocParallel opt-in round (4 new blocks in `test-build_signature.R`: serial
default == explicit `SerialParam`, parallel == serial bit-identical, RNG
state preserved, non-`BPPARAM` backend rejected) — log `logs/17_test.txt`.

## Gate 2 — adversarial review

The `test_audit.md` document lists every check that failed or passed only after
deliberate fixes (seeds, fold stratification, edge cases like `n_pcs == 1`,
non-syntactic gene symbols). After any change to `R/*.R`, re-run Gate 1 and
`R CMD check --as-cran` before proceeding.

Gate 2 also runs the **literal stale-number sweep** so that no superseded value
or phrase survives a refactor by memory alone:

```powershell
powershell -ExecutionPolicy Bypass -File validation/grep_stale_numbers.ps1
```

It scans `.R` (code and `@examples`), `.Rmd` (vignette and report template),
`.md`, and `.Rd` (man pages) for the superseded numbers (`0.217`, `0.424`,
`0.510`, `0.160`) and the old "poor generalization" / "winner's curse"
phrasing. Hits are permitted only under `validation/` (the intentional
before/after correction narrative) and the two narrative vignettes that
deliberately quote the superseded numbers as naive-tier counterfactuals
(`vignettes/case-study-brca.Rmd`, `vignettes/case-study-impact.Rmd`); any hit
elsewhere fails the gate (exit 1).
First run (2026-08-04): 15 hits, all inside `validation/`, 0 elsewhere.
After the 2026-08-04 guardrail/documentation round the sweep scanned 52 files
and found 19 hits, all still inside `validation/` (the intentional
before/after correction narrative), exit 0. After the `load_counts()` +
enforceable-lock round it scans 58 files (19 hits, same narrative), exit 0.
After the CV-optimism reframe round it scans 60 files (19 hits, same
narrative), exit 0. After the Phase-3 case-study/power round (2026-08-04) it
scans **157 files (20 hits, same narrative), exit 0**. After the Phase-1 BiocCheck-prep round (2026-08-05) it scans **67 files (20 hits, same
narrative), exit 0**. After the Phase-1 BiocParallel round (2026-08-05) it
scans **68 files (20 hits, same narrative), exit 0** — the one additional
file is the new `validation/make_logo.R` logo generator (no stale patterns).
After the three-tier case-study round (2026-08-05) it scans **69 files
(29 hits), exit 0** — the vignette allow-list above was added (8 naive-tier
counterfactual quotes in `case-study-brca.Rmd` and `case-study-impact.Rmd`,
plus 1 new quote in the four-vignette section of `paper_qa_summary.md`),
build-artifact directories (`*.Rcheck/`, `chk*/`) are now excluded from the
scan, and the README and the other two case-study vignettes stay clean.

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
- `sex_check` rank score ↔ the documented XIST-vs-Y rule re-derived from the
  assay (closes parity item 7).
- `build_signature` serial↔parallel equivalence (2026-08-05 BiocParallel
  opt-in): an explicit `BiocParallel::SerialParam()` run is bit-identical to
  the default, a `BiocParallel::SnowParam(2)` run is bit-identical to the
  serial run (same seeded folds; per-fold evaluation is deterministic Cox fit
  + `survival::concordance`), the caller's RNG state is left untouched, and a
  non-`BPPARAM` backend errors cleanly. These live in
  `test-build_signature.R` and are skipped automatically when BiocParallel is
  not installed.

Gate 3 also covers the two runtime guardrails added 2026-08-04: a fixture
with deliberately pre-logged, non-integer `counts` must fire
`possibly_log_scaled` (and a raw integer matrix must not), and a fixture
with few events relative to `top_n` must fire `events_per_parameter` (and an
adequately powered one must not). These live in `test-build_signature.R`,
`test-validate_external.R`, and `test-pca_audit.R`.

Gate 3 also covers the 2026-08-04 intake and lock features:
`test-load_counts.R` asserts every intake flag (duplicate samples, non-integer
counts, empty rows, Ensembl/non-syntactic IDs) fires only on its trigger and
that the flag ledger survives into the returned object; `test-lock_signature.R`
asserts that a locked signature makes a second `run_rnaSentry()` call refuse
to run and that `lock_signature(sig, lock = FALSE)` restores the ability to
re-run.

## Gate 4 — simulation study (15 falsifiable targets)

```r
source("validation/simulate_study.R")   # exits non-zero on any failure
```

Covers planted sex-label swaps, planted confounder redundancy, the null-data
CV-concordance expectation (with a random-gene calibration control), a
planted time-varying-hazard scenario, the external-validation lock check, and
the Phase-3 **external-transfer power analysis** (Sim 6: two tiers calibrated
to effective C = 0.608 / 0.654; power monotone in external events, < 0.30 at
~35 events, >= 0.80 at ~300 events for the moderate tier). Expected:
**15 targets / 15 PASS**. Results and interpretation:
`simulation_study.md`. Logs: `logs/03_simulation.txt` (Sims 1-5),
`logs/13_simulation.txt` (full suite); re-run unchanged 2026-08-05 after the
BiocParallel round — `logs/17_simulation.txt`.

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
note in `face_validity_review.md`. As of 2026-08-04 `build_signature()`
*catches this automatically*: 4.2 < `min_events_per_parameter` (default 5),
so a re-run of the repro flags `events_per_parameter` and the report audit
trail surfaces it. Log: `logs/04_face_validity.txt`, report:
`logs/face_validity_report.html`.

## Gate 5b — live GEO cross-cohort test (Phase 3)
`validation/repro_live_geo.R` downloads (or reads from the local cache under
`validation/cache/`) two independent GEO LUAD cohorts, GSE31210 (discovery,
226 tumors, 35 deaths) and GSE50081 (external validation, 128
adenocarcinomas, 52 deaths), then runs the full pipeline on the discovery
cohort and re-scores the external cohort against the discovery cutpoint.

Two cross-cohort pitfalls are solved deterministically in the script:

1. **Scale mismatch.** GSE31210's deposited series matrix is RMA-*linear*
   (median ≈ 99) despite the paper's "log2" wording; GSE50081 is log2 RMA.
   The script log2-transforms at the probe level when the median exceeds 20
   so both cohorts are on the same scale.
2. **Batch drift vs. the discovery-median cutpoint.** Raw-score transfer
   leaves the external cohort 3–6 discovery-SD below the cutpoint (empty risk
   group). Per-gene **z-scoring within each cohort** makes the cutpoint
   scale-free; probes collapse to the **mean over all probes per gene**
   (a per-cohort best-probe rule picked different probes per cohort and was
   rejected).

Expected output and the gate results (2026-08-04): discovery CV C = 0.862 vs
permutation-null 0.836 (calibrated PASS), external C = 0.540 / log-rank
p = 0.266 (**scientific transfer not demonstrated — honest negative**),
mechanical gates all PASS. The run also documents **selection leakage**: gene
selection happens on the full cohort before the CV split, so the CV
concordance is screening-internal and even survival-permuted data yields a
high null (≈0.8) on this ~21k-gene panel; the package's gates calibrate
against it and `validate_external()` is the only fully out-of-sample
estimate. Log: `validation/logs/11_live_geo.txt`, report:
`validation/logs/live_geo_discovery_report.html`. Result write-up:
`face_validity_review.md` (live GEO section),
`submission_success_criteria.md` (Section A gate).

## Gate 5c — case-study vignettes and bundled data (Phase 3)

`validation/repro_gse20685.R` builds `inst/extdata/gse20685_case_study.rds`
from the GEO series matrix (reads the cached `series.rds` when available, else
downloads via GEOquery; the raw file is wrapped in a list, so the script
unwraps with the same `eset[[1]]` pattern as `repro_live_geo.R`). Expected
output: 3000 genes x 327 samples, 83 events, colData `time/event/age/subtype`.

The three case-study vignettes (`vignettes/case-study-brca.Rmd`,
`case-study-confounder-audit.Rmd`, `case-study-small-cohort.Rmd`) must build
**offline** via `rmarkdown::render()` (or, as part of Gate 0, inside
`R CMD build`). Their inline numbers are pinned in
`validation/paper_qa_summary.md`:
BRCA subset CV C = 0.797 (sd 0.027), log-rank p = 2.99e-15, recommended
formula `~ age + subtype`; confounder case Cramer's V = 0.82 with recommended
`~ batch`; small-cohort case 23 events / 5 genes (4.6 < 5, guardrail fires),
mean CV C = 0.743 with per-fold spread 0.50-1.00.

## Gate 6 — cross-platform (skipped by decision 2026-08-05)

Needs a Docker-capable machine (`validation/docker_check.sh` runs the devel
container: build + BiocCheck + `--as-cran` + tests), or rhub/win-builder
(needs a GitHub token / an email address). Phase 4 (rhub) was explicitly
cancelled by the maintainer on 2026-08-05; Gate 6 is therefore NOT met and is
the only outstanding "ready to submit" criterion (see
`submission_success_criteria.md`, Section D). Status and instructions:
`cross_platform.md`. If a later run happens, append the per-platform
`Status:` lines to that file and to the submission notes.

## Full regression loop after any code change

1. `R CMD build` → `R CMD check --as-cran` on the tarball (0 ERROR / 0 WARNING;
   the sole NOTE is the deliberately-skipped `tidy`). Remember to export
   `RSTUDIO_PANDOC`, `R_QPDF`, `_R_CHECK_DOC_SIZES_=true`, and
   `_R_CHECK_CRAN_INCOMING_=false`.
2. Gate 1 suite (147/457; 32 expected guardrail warnings on small fixtures).
3. If statistics/tests changed: Gate 3 parity suite, Gate 4 simulation
   (15/15).
4. If feature/annotation handling changed: Gate 5 GSE20685 repro.
5. If vignettes or bundled data changed: Gate 5c offline vignette build.
6. Gate 2 stale-number sweep (`validation/grep_stale_numbers.ps1`).
7. If the DESCRIPTION URL, `_pkgdown.yml`, the logo, README, NEWS, or
   vignettes changed: re-run `pkgdown::build_site()` (with `RSTUDIO_PANDOC`
   set) and confirm the sitrep is clean (`URLs ok`, `Favicons ok`). The site
   lives in `docs/` (excluded from the tarball via `.Rbuildignore`); the logo
   is regenerated from `validation/make_logo.R` (base R only). **Publishing is
   automatic**: `.github/workflows/pkgdown.yaml` rebuilds and deploys to
   `gh-pages` on every push to `main` (2026-08-05) — never push to `gh-pages`
   by hand; the manual local `build_site()` is a preview only.
8. If any of the above touched `main`: confirm the pushed commit's
   `.github/workflows/pkgdown.yaml` run is green (GitHub Actions → pkgdown);
   it installs the package with Bioc release deps, builds the full site, and
    redeploys `gh-pages` (Gate 7).
9. Commit logs with the change at a logical checkpoint.
