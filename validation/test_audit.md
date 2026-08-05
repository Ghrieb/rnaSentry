# rnaSentry test-suite audit

Scope: the 15 `tests/testthat/` files (457 expectations across 147 test
blocks per `devtools::test()`) reviewed for what they actually verify,
against the documented guarantees in `man/`. This audit feeds (a) the Step-1
adversarial pass and (b) the Step-2 statistical-parity tests.

## Inventory (file -> expectations -> what is verified)

| File | expect_ calls | Verifies |
|---|---|---|
| test-build_signature.R | 71 | input type/col guards, option ranges, signature shape, `top_n`/`p_value` selection, fewer-genes flag, collinear-drop flag, fixed-seed reproducibility, design-adjusted screening exclusion, `min_events_per_parameter` validation, `possibly_log_scaled` flag on pre-logged assays, `events_per_parameter` flag on low-event cohorts, print, **serial↔parallel equivalence (BPPARAM opt-in, 2026-08-05): default == explicit `SerialParam`, `SnowParam` == serial bit-identical, caller RNG state preserved, non-`BPPARAM` backend rejected** |
| test-cox_model.R | 34 | input/confounder guards, table schema, HR<=CI bounds, gene coefficients == signature coefficients, design-term inclusion, NA covariate complete-cases, lock state, print/plot |
| test-design_audit.R | 44 | type/range guards, formula text, confounder-table schema + BH, lm-vs-factor test dispatch, single-level untested, fully-missing flag, redundant numeric drop, PC clamp, single-PC interaction skip, print |
| test-validate_external.R | 39 | input guards, cutpoint never recomputed from external data, structure, log-rank finite, **concordance == survival::concordance**, strict missing-gene default, `drop_missing` opt-in + overlap frac, `possibly_log_scaled` flag on a pre-logged external cohort, lock state, print/plot |
| test-km_curve.R | 30 | input guards, constant-score error, median/custom cutpoint, cutpoint range error, high-vs-low median order, lock state, print/plot |
| test-load_counts.R | 23 | matrix/colData/ID guards, duplicate-sample critical flag + dedup, non-integer warning flag, empty-row warning flag, Ensembl/non-syntactic info flags, flag-ledger structure, SE round-trip |
| test-survival_parametric.R | 26 | input guards, dist set, AIC table self-consistency (weights sum 1, delta>=0, sorted), km_fit, curves bounds, score == km score, indirect exp-vs-weibull AIC check, lock state, print/plot |
| test-pca_audit.R | 41 | schema, assay preference, batch-associated PC flag, lm for numeric batch, single-level batch, null batch not flagged, guards, constant-gene filter, too-few-genes stop, BH parity, `possibly_log_scaled` flag on a pre-logged assay, print, ggplot |
| test-sex_check.R | 13 | input/col guards, missing markers error, mismatch/OK/AMBIGUOUS scenarios, rank-score parity against the documented XIST-vs-Y rule |
| test-lock_signature.R | 31 | type/arg guards, lock/unlock state transitions, audit entries, stage schema, round-trip, print, session-lock enforcement (`run_rnaSentry()` refuses to re-run while a lock is set; unlock permits re-run), targeted unlock leaves other locked signatures protected |
| test-qc_explore.R | 10 | type guard, structure, duplicate-sample flag, library-size outlier, non-integer counts, print |
| test-run_rnaSentry.R | 13 | full pipeline + report render, skip-report, input guards, lock-refusal at the run level |
| test-generate_report.R | 7 | named-list guard, output-dir guard, report renders non-empty, report HTML contains the Pipeline assumptions section |
| test-adversarial_regressions.R | 36 | regression pins for REAL-BUG A/B/C and soft-warns below |
| test-statistical_parity.R | 24 | independent-reference replays for parity targets 1-6b below |

## Strengths

- **Guardrail depth.** Every public function has type/range/missing-input
  guards exercised (`test-build_signature.R:3-54`, `test-design_audit.R:8-36`,
  etc.), and every downstream stage is tested against the invalid-data family.
- **Scientific non-regressions are encoded as scenarios**, not shapes:
  - batch-driven genes excluded only when adjusted (`test-build_signature.R:146-185`);
  - cutpoint is never recomputed from external data (`test-validate_external.R:23-40`);
  - collinear genes dropped with `coefficient_unstable` flag (`test-build_signature.R:107-124`);
  - high-risk group median survival lower than low (`test-km_curve.R:73-83`).
- **Reproducibility**: fixed seed => identical `cv_results`/`coefficients`
  (`test-build_signature.R:126-134`).
- **Lock-state guardrails** run through km/cox/parametric/validate
  (`test-km_curve.R:85-90`, `test-cox_model.R:84-89`, `test-validate_external.R:104-111`),
  and `run_rnaSentry()` refuses to re-run while any session lock exists
  (`test-lock_signature.R` enforcement test).
- **One true parity test already exists**: `validate_external$concordance`
  against `survival::concordance` (`test-validate_external.R:62-67`).

## Reviewer-caught defect: concordance sign convention (fixed)

A real-data review of GSE20685 (CV C = 0.217, ~9 SD below chance) showed the
`survival::concordance(Surv ~ x)` default (`reverse = FALSE`) means "larger x ⇒
longer survival", so both risk-score call sites reported `1 − Harrell's C`:

- `R/build_signature.R` (held-out CV concordance) — fixed with `reverse = TRUE`.
- `R/validate_external.R` — same fix.
- `R/cox_model.R` was already correct (uses the coxph method's own concordance).

The existing parity tests could not catch this because they replayed the same
inverted call. The blind spots are now closed with:
- `reverse = TRUE` in both parity replays (`test-statistical_parity.R`,
  `test-validate_external.R`).
- Direction guards: CV mean > 0.5 on a signal-bearing fixture; `validate_external`
  concordance > 0.5 on a cohort whose survival is engineered to follow the
  signature score.
- A `C + C_rev = 1` invariant test in `test-statistical_parity.R`.

Impact on previously reported numbers (all `1 − C`): GSE20685 CV 0.217 → 0.783;
null-data simulation 0.424 → 0.576; control 0.510 → 0.490; `validate_external`
0.160 → 0.840. The "winner's-curse below 0.5" narrative in the early
face-validity review was an artifact of this bug and has been corrected.

## Statistical-parity coverage (Step-2 targets)

`tests/testthat/test-statistical_parity.R` re-implements each statistic the
package computes itself with an independent reference implementation. Status of
the Step-2 targets:

| # | Statistic | Independent reference | Status |
|---|-----------|-----------------------|--------|
| 1 | `km_curve`/`validate_external` log-rank p | `survival::survdiff` | closed |
| 2 | `design_audit` Cramér's V | `sqrt(chisq/(n*(min(r,c)-1)))` | closed |
| 3 | `cox_model` HR / CI / p | independent `survival::coxph` + `stats::confint` | closed |
| 4 | `survival_parametric` AIC / loglik / npar | `stats::AIC` / `logLik` | closed |
| 5 | `cox_model` Schoenfeld rho / p | `survival::cox.zph` | closed |
| 6 | `build_signature` CV C-index | `survival::concordance` on replayed event-stratified folds (`reverse = TRUE`) | closed |
| 6b | `cox_model` concordance direction | `survival::concordance(Surv ~ lp, reverse = TRUE)` on the fitted linear predictor | closed |
| 7 | `sex_check` rank score | documented XIST-vs-Y scoring rule | closed |

Item 7 is now closed: `test-sex_check.R` re-derives the documented rule
(score = `rank(XIST) - rank(Y panel)`, `>= 1` -> F, `<= -1` -> M, else
ambiguous) independently from the assay and asserts `sex_check()` reproduces
the inference on a fixture spanning all three branches.

## Parallel-equivalence coverage (2026-08-05 BiocParallel opt-in)

`build_signature()` / `run_rnaSentry()` now accept `BPPARAM` (Bioc-native,
strictly opt-in; `NULL` keeps the serial default). Because all fold
partitioning stays inside `withr::with_seed()` and per-fold evaluation is a
deterministic Cox fit + `survival::concordance` with no RNG calls, parallel
fold evaluation is bit-identical to serial. `test-build_signature.R` pins
this:

- serial default `==` explicit `BiocParallel::SerialParam()` run
  (`expect_identical` across `cv_results`/`cv_summary`/`genes`/`coefficients`/`flags`);
- `BiocParallel::SnowParam(2)` run `==` serial run (same five fields);
- caller's `.Random.seed` is untouched after a parallel run;
- a non-`BPPARAM` backend errors cleanly.

These tests `skip_if_not(requireNamespace("BiocParallel"))`, so the suite is
green without BiocParallel installed. `BiocParallel` is Suggests-only and all
calls are namespace-qualified behind a `requireNamespace()` guard.

## Guardrail coverage (2026-08-04 additions)

Two runtime guardrails were added so misuse is caught at run time, not just
described in prose. Each fires only on a specific, test-pinned condition.

| Check | Trigger | Where | Test that pins it |
|---|---|---|---|
| `possibly_log_scaled` | first assay is non-integer **and** max < 40 while no `logcounts`/`vst` assay exists (double-log risk) | `build_signature`, `validate_external`, `pca_audit` (warning + flag) | `test-build_signature.R` (pre-logged fires, raw integer does not), `test-validate_external.R`, `test-pca_audit.R` |
| `events_per_parameter` | events / final signature genes < `min_events_per_parameter` (default 5; source documented in `man/build_signature.Rd`) | `build_signature` (warning + flag) | `test-build_signature.R` (low-event fires, adequate does not) |

## Edge-case / coverage gaps

- Time values: `time = 0`, negative time (guarded for NA/negative at first
  element only), all-equal times, ties in event times.
- Expression: NA/NaN/Inf in the count matrix, all-zero rows, single-column SE,
  duplicate rownames, extreme counts (integer overflow in `*3L`-style scalings).
- Folds/repeats: `folds > ncol(se)`, `repeats` producing identical fold seeds,
  fold splits leaving a fold with <2 events.
- `cox_model()` with zero variance in a confounder handled; zero variance in a
  *signature gene* across the cohort not explicitly tested.
- `km_curve()`/`validate_external()` with a signature whose genes include the
  survival/design columns by name collision.
- `design_audit()` with an ordered-factor design var, and with factor levels
  containing spaces/special chars.
- `pca_audit()` with `n_samples < 3` (prcomp degeneracy) and with
  `top_n_pcs > n_samples - 1`.

## Verification status (gate logs, see logs/)

- tests: **147 blocks / 457 passed / 0 failed / 0 error / 32 warnings**
  (2026-08-04, after the `load_counts()` intake wrapper + enforceable-lock
  round and the Batch C lock-targeted-unlock + concordance-flag round;
  re-confirmed at `logs/17_test.txt` on 2026-08-05 after the Phase-1
  BiocCheck-prep + BiocParallel rounds — withr-scoped seeding, `seq_len()`,
  runnable examples, plus the 4 new serial↔parallel blocks that took the
  suite from 143/445 to 147/457). The
  32 warnings are the
  `events_per_parameter` guardrail firing on deliberately small synthetic
  fixtures in tests that exercise other behavior; each guardrail's own
  dedicated test asserts that firing. They are expected, honest noise, not
  regressions.
- R CMD check tarball: 0 ERROR (00.2 pre-fix; re-confirmed 0 ERROR after the
  sign-convention fix — 05/06/07)
- BiocCheck source + tarball: 1 ERROR (support-site email, environmental),
  12 NOTES at the 00.3/00.5 runs and 13 NOTES at the 2026-08-04 re-run
  (`logs/13_bioccheck.txt`, which adds the `Avoid 1:` note from
  `load_counts.R`; re-confirmed at `logs/14_bioccheck.txt` after the Batch-A
  round — also cleared the one-off "data files exceed 5MB" warning via the
  extdata xz recompress; re-confirmed at `logs/15_bioccheck.txt` after the
  Batch C round); at the Phase-1 BiocCheck-prep round (2026-08-05) the
  justified `set.seed` WARNING is gone — reproducibility is now scoped with
  `withr::with_seed()` (`logs/16_bioccheck.txt`: 1 ERROR (support-site email,
  environmental) / 0 WARNING / 8 advisory NOTES); GitClone 0 ERROR / 1 WARNING
  (CITATION doi) — all explained in `logs/notes_documented.md`;
  re-confirmed unchanged at the BiocParallel round
  (`logs/17_bioccheck.txt`, `logs/17_bioccheck_tarball.txt`,
  `logs/17_bioccheck_gitclone.txt`)
- `--as-cran`: 0 ERROR / 1 WARNING (qpdf, environmental) / 1 NOTE (tidy,
  environmental) (05/06/07, post sign-convention fix; 08 after the
  2026-08-04 guardrail + docs round; 09 after the `load_counts()` +
  enforceable-lock round — 1 WARNING / 2 NOTEs because of a one-off
  `unable to verify current time` note; 10/11 back to 1 WARNING / 1 NOTE;
  12 after the CV-optimism reframe doc round; 13 after the
  Phase-3 case-study/power round; 14 after the Batch-A ship-blocking round
  (LICENSE holder, report-template `results='asis'`, S4Vectors → Imports,
  xz-recompressed extdata); 15 after the Batch C round (targeted lock unlock,
  `min_events_per_parameter` forwarding, `concordance_na` flag, flags `stage`
  column docs) — all clean 0 ERROR / 1 WARNING / 1 NOTE; 16 after the
  Phase-1 BiocCheck-prep round — `qpdf` now installed, so 0 ERROR / 0
  WARNING / 1 NOTE, the residual HTML `tidy` NOTE environmental;
  re-confirmed at the BiocParallel round — `logs/17_as_cran.txt`,
  0 ERROR / 0 WARNING / 1 NOTE (`tidy` only))
- stale-number sweep (2026-08-04; re-run 2026-08-05): `validation/grep_stale_numbers.ps1` scans
  `.R`/`.Rmd`/`.md`/`.Rd` for `0.217`/`0.424`/`0.510`/`0.160` and the old
  "poor generalization" / "winner's curse" phrasing — 20 hits (68 files
  scanned at the 2026-08-05 BiocParallel round, +1 from the new
  `validation/make_logo.R`), all inside `validation/` (the intentional
  correction narrative), 0 in code, `@examples`, the vignette, the report
  template, or the man pages (PASS, exit 0); after the three-tier case-study
  round the sweep scans **69 files (29 hits)**, the extra hits being the
  naive-tier counterfactual quotes in `vignettes/case-study-brca.Rmd` and
  `vignettes/case-study-impact.Rmd` (now explicit allow-list entries
  alongside `validation/`) and one new quote in `paper_qa_summary.md`;
  `*.Rcheck/`/`chk*/` build-artifact directories are excluded from the scan
  (PASS, exit 0); re-run 2026-08-05 after `CONTRIBUTING.md` gained the
  site-deployment section — still **69 files / 29 hits**, all inside
  `validation/` or the two allow-listed narrative vignettes, exit 0
- CI / site round (2026-08-05): **no R code or test changes**, so Gate 1
  stands at 147 blocks / 457 expectations / 0 fail / 0 error / 32 expected
  guardrail warnings, and the as-cran / BiocCheck / simulation numbers are
  unchanged. New automation verified green: `.github/workflows/pkgdown.yaml`
  run #1 installed the package with Bioc release deps on ubuntu-latest,
  built the full pkgdown site (all 4 case-study vignettes + main vignette)
  with the same pinned numbers as the local runs, and deployed to `gh-pages`
  (`3893b77`; prior manual publish `d36a20f`). `.Rbuildignore` excludes
  `^\.github$`, so the tarball is unaffected (verified by re-building the
  tarball 2026-08-05 — no `.github` in the archive). Log:
  `logs/18_ci_site.txt`.
- Phase-3 case-study/power gate (2026-08-04, commit `b7e080d`): simulation
  suite **15/15 falsifiable targets PASS** (Sims 1-5 + Sim 6 external-transfer
  power analysis; `logs/13_simulation.txt`); three offline case-study
  vignettes + bundled `inst/extdata/gse20685_case_study.rds`
  (`repro_gse20685.R`, Gate 5c); `logs/13_as_cran.txt` clean baseline; the
  "TCGA-BRCA GSE20685" mislabel fixed to GSE20685 (Li et al., 2010;
  Affymetrix GPL570)

## Priorities

1. Step-2 parity tests for gaps 1-6 above (direct Bioconductor reviewer value).
2. Adversarial pass (Step 1): fold/event degeneracies, NA/Inf matrices,
   one-column data, duplicate rownames, extreme counts.
3. Simulation study (Step 3) re-uses the engineered-signal fixtures above and
   extends them to the falsifiable targets listed in simulation_study.md.

## Step 1: adversarial pass (2026-08-03)

A fresh-context agent probed ~90 degenerate inputs across all public
functions. 78+ probes were correctly rejected by existing guardrails or
returned flagged output. The pass surfaced 3 bug families + 4 soft-warns,
all fixed and pinned by tests in `tests/testthat/test-adversarial_regressions.R`
(+30 expectations; the suite stood at 393 passed / 0 failed at this
checkpoint; later rounds raised it to 414).

### Fixed: REAL-BUG A — `stats::sd(x) == 0` NA-crash
`sd()` is NA for length-1 vectors or vectors containing `Inf`, and `if(NA)`
threw `missing value where TRUE/FALSE needed`. Sites:
- `build_signature.R` time check (single-sample cohort, `Inf` time)
- `build_signature.R` fold risk-score check (leave-one-out / large folds)
- `km_curve.R` constant-score check (single-sample cohort)

Fixes: reject non-finite follow-up time in all five survival stages;
`isTRUE(sd(x) == 0)` at the remaining sites; skip single-sample test folds;
reject `folds >= n` with a clear message. Covered by
`test-adversarial_regressions.R` (LOO/folds>n, Inf time, n=1, n=2).

### Fixed: REAL-BUG B — gene names colliding with internal model columns
A gene literally named `time` or `event` was renamed to `time.1`/`event.1`
by `data.frame(check.names)`, corrupting the signature's gene identity and
crashing CV with `subscript out of bounds`. Fix: such genes are excluded
before screening with a `gene_name_collision` warning flag.

### Fixed: REAL-BUG C — exactly collinear genes silently degraded
The joint Cox fit flagged only non-finite coefficients, but exact
collinearity produced finite coefficients ~1.8e5 (exp -> Inf), so
`build_signature` returned HR-overflow garbage with zero flags and
`cox_model` returned HR = Inf/0 with zero flags. Fixes: the joint-fit
instability check now also treats `!is.finite(exp(b))` as unstable;
`cox_model` adds a `non_estimable_coefficients` warning when any
coefficient/SE/HR/CI is non-finite.

### Fixed: SOFT-WARNs
- `design_audit` recommended a formula containing constant/all-NA variables
  -> such variables are now excluded from the formula and flagged
  `variable_constant`.
- `qc_explore` silently passed an all-NA count assay -> new `missing_counts`
  critical flag; library-size outlier detection is skipped (not crashed)
  when library sizes are not finite.
- `cox_model` silently dropped a confounder equal to the survival columns ->
  now a clear error ("confounders must be separate from the time/event
  columns").
- `sex_check` returned NA status for missing reported sex -> new `MISSING`
  status (documented in the man page).

### Minor
- `pca_audit` on a single-sample cohort now errors cleanly ("at least two
  samples") instead of a raw `prcomp` error.
