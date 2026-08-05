# rnaSentry — submission success criteria

Definition of done for the Bioconductor submission. Two perspectives must
both hold: **statistical validity** (the analysis a reviewer can trust) and
**Bioconductor compliance** (the package a reviewer can install and check).
Every criterion below maps to an artifact and a gate; the status column is
updated as gates are re-run. Maintained alongside `maintainer_testing_guide.md`
(how to run each gate) and `paper_qa_summary.md` (the numbers for the paper).

---

## A. Statistical validity — what the tool must demonstrably do

| # | Criterion | How it is proven | Status |
|---|---|---|---|
| S1 | Every statistic the package computes matches an independent reference implementation | Statistical-parity suite (`tests/testthat/test-statistical_parity.R`): log-rank vs `survival::survdiff`; Cramér's V vs closed form; Cox HR/CI/p vs independent `coxph`; AIC/loglik vs `stats::AIC`; Schoenfeld vs `cox.zph`; held-out CV C vs replayed seeded folds; `C + C_rev = 1` invariant; `sex_check` score vs documented XIST-vs-Y rule (parity items 1–7 closed) | PASS (all parity tests green in Gate 1) |
| S2 | Guardrails fire only when their trigger is present | Dedicated fire/no-fire tests for `possibly_log_scaled`, `events_per_parameter`, `duplicate_samples`, `non_integer_counts`, `empty_rows`, `ensembl_ids_detected`, `non_syntactic_ids`, `signature_locked`, and the run-level lock refusal | PASS (Gate 1) |
| S3 | No winner's-curse inversion on null data | Simulation study: null-data CV concordance ≈ 0.5 (mean 0.576 full-pipeline on nulls, random-gene control 0.490 on simulated nulls), direction guarded by `reverse = TRUE` convention | PASS (Gate 4, 15/15 falsifiable targets incl. Sim 6) |
| S4 | Real data discriminates beyond chance | GSE20685: held-out CV C = 0.783 vs random-gene control 0.582 (6.3 SD above control, logged CV sd 0.032); log-rank p = 3.8e-17; MKI67 HR > 1, ESR1 HR < 1 (literature directions); `design_audit` flags subtype (eta² = 0.76) | PASS (Gate 5; extended live test below) |
| S5 | External validation is honest | `validate_external()` never recomputes the cutpoint from external data (test-pinned); discovery cutpoint applied unchanged; external C reported on an independent cohort | PASS (Gate 5 + live GEO pair, Phase 3) |
| S6 | Results are reproducible | Fixed seed ⇒ identical `cv_results`/`coefficients` (test-pinned); the enforceable lock prevents silent gene-set re-selection after survival analysis; `run_rnaSentry()` refuses a second run until `lock_signature(lock = FALSE)`; `BPPARAM` parallel fold evaluation is **bit-identical to serial** (test-pinned, 2026-08-05) | PASS (Gate 1) |
| S7 | Degenerate inputs are rejected or flagged, never silently wrong | Adversarial pass: ~90 probes; REAL-BUG A/B/C + 4 soft-warns found and fixed, pinned by regression tests; NAs/Inf/zero-variance/single-sample/fold degeneracies covered | PASS (Gate 2, 2026-08-03) |
| S8 | The audit trail is complete | Every automated decision carries `(check, severity, detail, stage)`; report renders the full flag ledger and the pipeline-assumptions section | PASS (Gate 1 `test-generate_report.R`) |
| S9 | Power/transfer behaviour is quantified and falsifiable | Sim 6: external-transfer power monotone in external events; < 0.30 at ~35 events (small-cohort trap); >= 0.80 at ~300 events for a C≈0.65 signature | PASS (Gate 4, 15/15 targets) |
| S10 | Case-study narrative runs on real bundled data, offline | `inst/extdata/gse20685_case_study.rds` + three offline-safe vignettes (BRCA pipeline, confounder audit, small-cohort guardrail) | PASS (Phase-3 vignettes build locally) |

**Phase-3 live acceptance gate (GSE31210 → GSE50081, LUAD)** — executed
2026-08-04, log `validation/logs/11_live_geo.txt`:
- `load_counts()` on real GEO data returns a flag ledger (recorded:
  `non_integer_counts`, `empty_rows`; no false-criticals). **PASS**
- Discovery `run_rnaSentry()` completes and renders a report (7-gene signature;
  CV C = 0.862). **PASS**
- Lock refuses a second run; `lock_signature(lock = FALSE)` releases it
  (0 session-lock entries after unlock). **PASS**
- `events_per_parameter` respected (35 events / 7 genes = 5.0, guardrail not
  tripped). **PASS**
- CV C above the matched **permutation-null** (0.862 vs 0.836, delta +0.026).
  **PASS (calibrated)** — note the CV is screening-internal (see below).
- External validation on GSE50081: **C = 0.540, log-rank p = 0.266** →
  scientific transfer **NOT demonstrated** on this low-power dataset (35
  discovery events). Reported as an **honest negative** (mechanics validated;
  no significant independent-cohort stratification), not a tool failure.
- **CV-optimism finding (selection leakage):** genes are selected on the full
  cohort before the CV split, so even survival-permuted data gives null CV C
  ≈ 0.8 (not 0.5) on this ~21k-gene panel. The package's calibration handles
  this, and the package docs (Rd/README/vignette) now state it explicitly.
  `validate_external()` is the only fully out-of-sample estimate.

**Phase-3 case-study + power gate** — executed 2026-08-04, logs
`validation/logs/13_simulation.txt`:
- Sim 6 (external-transfer power): two tiers calibrated to effective
  C = 0.608 / 0.654; transfer power rises with external events (weak
  0.067 → 0.417; moderate 0.275 → 0.842 over ~35 → ~300 events); monotone in
  events; **< 0.30 at ~35 events**, **≥ 0.80 at ~300 events** for the
  moderate tier. **PASS (4/4 new targets, 15/15 suite).**
- Three offline-safe case-study vignettes build locally:
  `case-study-brca`, `case-study-confounder-audit`, `case-study-small-cohort`
  (BRCA subset: CV C = 0.797 sd 0.027; log-rank p = 2.99e-15; `subtype`
  flagged, recommended formula `~ age + subtype`).
- README now carries the Status note, Case studies 1–4 (incl. the LUAD
  honest-negative), Positioning and prior art, and Contributing sections;
  `CONTRIBUTING.md` and `NEWS.md` added.

## B. Bioconductor rules — what the package must pass

| # | Rule | Current status |
|---|---|---|
| B1 | `R CMD build` produces a clean tarball | PASS — vignette compiles, `inst/doc` present |
| B2 | `R CMD check --as-cran` ≤ 1 WARNING / 1 NOTE, both environmental (`qpdf`, `tidy`) | PASS — 0 ERROR throughout logs 05–15; log 16 (`16_as_cran.txt`, 2026-08-05 Ultimate Pre-Flight) is **0 ERROR / 0 WARNING / 1 NOTE** — the `qpdf` WARNING cleared after installing the qpdf CLI 12.3.2 and running with `R_QPDF` + `_R_CHECK_DOC_SIZES_=true`; the sole remaining NOTE is `tidy` (HTML validation), deliberately skipped per maintainer decision; re-confirmed unchanged at `17_as_cran.txt` (2026-08-05 BiocParallel round) |
| B3 | BiocCheck: no package errors; only documented items | 1 environmental ERROR (support-site email 404 → fixed by registering `ghriebabdelkarimhani@gmail.com` on https://support.bioconductor.org); at `16_bioccheck.txt`/`16_bioccheck_tarball.txt` (2026-08-05 Ultimate Pre-Flight) **0 WARNING / 8 advisory NOTES** — the `set.seed` WARNING (→ `withr::with_seed()` refactor), the R-version NOTE (→ `R >= 4.5.0`), the `Avoid 1:` NOTE (→ `seq_len()`), and the two `\dontrun`/runnable-examples NOTES (→ self-contained `\donttest`) are all resolved; the remaining 8 are justified in `logs/notes_documented.md`; re-confirmed unchanged at the BiocParallel round (`17_bioccheck.txt`/`17_bioccheck_tarball.txt`); GitClone `16_bioccheck_gitclone.txt` is **0 ERROR / 1 WARNING (CITATION `doi`, kept by decision) / 0 NOTES**, re-confirmed at `17_bioccheck_gitclone.txt` |
| B4 | Cross-platform (devel: Linux + macOS + Windows, R-devel) | Gate 6 — **SKIPPED by maintainer decision 2026-08-05** (Phase 4 / rhub was cancelled; no Docker on this machine); local Windows/R 4.5.2 is green; **partial Linux evidence added 2026-08-05**: the CI pkgdown workflow (ubuntu-latest, R 4.6.1 / Bioc release) installs the package and compiles all six vignettes successfully (run #1 with the five original, run #3 `11d575a` adding `getting-started`) — not a full `R CMD check`/test-suite run; see `cross_platform.md` |
| B5 | Vignette builds and is informative | PASS — `inst/doc/rnaSentry.html` builds, plus three offline case-study vignettes (`case-study-brca`, `case-study-confounder-audit`, `case-study-small-cohort`); covers intake, QC, sex check, pipeline, external validation, lock, limitations, case studies |
| B6 | NEWS is complete and truthful | PASS — covers all features incl. `load_counts()` and the enforceable lock |
| B7 | DESCRIPTION fields complete (`Authors@R` with maintainer email, `biocViews`, `License`, `URL`, `BugReports`) | PASS — `biocViews`: Software, GeneExpression, RNASeq, DifferentialExpression, Survival, QualityControl, BatchEffect, Normalization |
| B8 | Version 0.99.x targets Bioc devel | PASS — 0.99.0 |
| B9 | GitHub: package-only default branch, maintainer SSH key, `gh` access | PENDING — user task (push + open Contributions issue) |
| B10 | Bioconductor Support Site account registered with the maintainer email | PENDING — user task |

## C. Evidence map

| Gate | What it runs | Expected | Artifact |
|---|---|---|---|
| Gate 1 | Full test suite | 147 blocks / 457 expectations / 0 fail / 0 error (32 expected guardrail warnings) | `logs/00.1_test.txt`; re-confirmed 2026-08-05 at `logs/17_test.txt` (post-`withr` refactor + BiocParallel opt-in round) |
| Gate 2 | Adversarial review + stale-number sweep | 0 stale numbers outside `validation/`; hits inside (correction narrative) | `test_audit.md`, `grep_stale_numbers.ps1`; sweep re-run 2026-08-05 (BiocParallel round): **68 files / 20 hits**, all inside `validation/`, exit 0; re-run 2026-08-05 (three-tier case-study round): **69 files / 29 hits**, all inside `validation/` or the two narrative-vignette allow-list entries (`case-study-brca.Rmd`, `case-study-impact.Rmd`), `*.Rcheck/`/`chk*/` build artifacts excluded, exit 0 |
| Gate 3 | Statistical parity | items 1–7 closed + serial↔parallel `BPPARAM` equivalence (2026-08-05) | `test-statistical_parity.R`, `test-build_signature.R` |
| Gate 4 | Simulation study | 15/15 falsifiable targets PASS (Sims 1-5 + Sim 6 power) | `simulation_study.md`, `logs/03_simulation.txt`, `logs/13_simulation.txt`, `logs/16_simulation.txt`, `logs/17_simulation.txt` (re-run 2026-08-05) |
| Gate 5 | Real-data face validity | discrimination beyond chance, correct biology, honest external check | `face_validity_review.md`, `logs/04_face_validity.txt` |
| Gate 5b | Live GEO pair (Phase 3) | see Section A acceptance gate | `repro_live_geo.R`, `logs/11_live_geo.txt` |
| Gate 5c | Case-study vignettes + bundled data (Phase 3) | offline build; numbers match documented values | `repro_gse20685.R`, `inst/extdata/gse20685_case_study.rds`, `vignettes/case-study-*.Rmd` |
| Gate 6 | Cross-platform | per-platform `Status:` lines | `cross_platform.md` |
| Gate 7 | Site + CI deploy | pkgdown build green + `gh-pages` updated on every push to `main` | `.github/workflows/pkgdown.yaml`, `logs/18_ci_site.txt` (2026-08-05: run #1 green, `gh-pages` `3893b77`) |

## D. Definition of "ready to submit"

All of S1–S10 pass, B1–B8 green, B9–B10 completed by the maintainer, Gate 6
records at least one non-Windows platform pass, and the Phase-3 live gate is
logged. As of 2026-08-05 all S1–S10 and B1–B8 are green on this Windows
machine; **Gate 6 (cross-platform) is the only outstanding "ready to submit"
criterion** — it was explicitly skipped by maintainer decision (Phase 4 /
rhub cancelled), so submission readiness depends on either a later
non-Windows check (rhub/win-builder/Docker) or an explicit decision to
accept that risk. At that point the GitHub push and the Contributions issue
(`Bioconductor/Contributions#...`, title `rnaSentry`) can be opened.

## Step-5 closure (local verification fully closed, 2026-08-05)

The Phase-1 BiocCheck-prep + BiocParallel go-live round re-ran **every** gate
on this machine and committed fresh logs (`17_*`), so the "verify before
Step 5 is closed" criterion is met in full:

- `R CMD check --as-cran` on a fresh tarball: 0 ERROR / 0 WARNING / 1 NOTE
  (`tidy` only) — `logs/17_as_cran.txt`.
- Full `devtools::test()`: 147 files / 457 expectations / 0 fail / 0 error /
  32 expected guardrail warnings / 0 skipped — `logs/17_test.txt`
  (serial↔parallel bit-identical parity pinned in `test-build_signature.R`).
- BiocCheck source + tarball: 1 ERROR (support-site email 404, environmental)
  / 0 WARNING / 8 advisory NOTES — `logs/17_bioccheck.txt`,
  `logs/17_bioccheck_tarball.txt`; GitClone 0 ERROR / 1 WARNING (CITATION
  `doi`) / 0 NOTES — `logs/17_bioccheck_gitclone.txt`.
- Simulation study: 15/15 falsifiable targets PASS — `logs/17_simulation.txt`.
- Stale-number sweep: 68 files scanned, 20 hits, all inside `validation/`,
  exit 0.

The remaining items before submission are **user actions only**, none of them
code gates: register the maintainer email on the Bioconductor Support Site
(clears the last BiocCheck ERROR), and open the Contributions issue. The
"enable GitHub Pages / verify the pushed site" action from earlier rounds is
**done (2026-08-05)**: GitHub Pages serves from the `gh-pages` branch and is
now updated automatically on every push to `main` (Gate 7).
