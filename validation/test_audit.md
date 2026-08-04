# rnaSentry test-suite audit

Scope: the 13 `tests/testthat/` files (322 expectations across 94 files per
`devtools::test()`) reviewed for what they actually verify, against the
documented guarantees in `man/`. This audit feeds (a) the Step-1 adversarial
pass and (b) the Step-2 statistical-parity tests.

## Inventory (file -> expectations -> what is verified)

| File | expect_ calls | Verifies |
|---|---|---|
| test-build_signature.R | 47 | input type/col guards, option ranges, signature shape, `top_n`/`p_value` selection, fewer-genes flag, collinear-drop flag, fixed-seed reproducibility, design-adjusted screening exclusion, print |
| test-cox_model.R | 34 | input/confounder guards, table schema, HR<=CI bounds, gene coefficients == signature coefficients, design-term inclusion, NA covariate complete-cases, lock state, print/plot |
| test-design_audit.R | 44 | type/range guards, formula text, confounder-table schema + BH, lm-vs-factor test dispatch, single-level untested, fully-missing flag, redundant numeric drop, PC clamp, single-PC interaction skip, print |
| test-validate_external.R | 35 | input guards, cutpoint never recomputed from external data, structure, log-rank finite, **concordance == survival::concordance**, strict missing-gene default, `drop_missing` opt-in + overlap frac, lock state, print/plot |
| test-km_curve.R | 30 | input guards, constant-score error, median/custom cutpoint, cutpoint range error, high-vs-low median order, lock state, print/plot |
| test-survival_parametric.R | 26 | input guards, dist set, AIC table self-consistency (weights sum 1, delta>=0, sorted), km_fit, curves bounds, score == km score, indirect exp-vs-weibull AIC check, lock state, print/plot |
| test-pca_audit.R | 39 | schema, assay preference, batch-associated PC flag, lm for numeric batch, single-level batch, null batch not flagged, guards, constant-gene filter, too-few-genes stop, BH parity, print, ggplot |
| test-sex_check.R | 11 | input/col guards, missing markers error, mismatch/OK/AMBIGUOUS scenarios |
| test-lock_signature.R | 23 | type/arg guards, lock/unlock state transitions, audit entries, stage schema, round-trip, print |
| test-qc_explore.R | 10 | type guard, structure, duplicate-sample flag, library-size outlier, non-integer counts, print |
| test-run_rnaSentry.R | 13 | full pipeline + report render, skip-report, input guards |
| test-generate_report.R | 4 | named-list guard, output-dir guard, report renders non-empty |

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
  (`test-km_curve.R:85-90`, `test-cox_model.R:84-89`, `test-validate_external.R:104-111`).
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

## Statistical-parity gaps (Step-2 targets)

The suite verifies *shape and self-consistency* but not agreement with
reference implementations for the statistics the package computes itself:

1. **Log-rank p** in `km_curve()` / `validate_external()`: only checked
   finite/in [0,1]. No parity vs `survival::survdiff`.
2. **Cramér's V** effect size in `design_audit()` confounder scan: only checked
   present/flagged. No parity vs `sqrt(chisq / (n*(min(r,c)-1)))`.
3. **Cox HR/CI** in `cox_model()`: coefficients compared to the signature's
   own internal Cox fit, but never to an independent `survival::coxph`
   (`exp(coef)`, `confint`).
4. **Parametric AIC** in `survival_parametric()`: table is self-consistent but
   not compared to `stats::AIC(survreg(...))` per distribution.
5. **Schoenfeld PH diagnostics** in `cox_model()` `$zph_summary`: shape only,
   no parity vs `survival::cox.zph`.
6. **CV C-index** in `build_signature()` `cv_results`: bounded, not compared to
   `survival::concordance` on the same held-out folds.
7. **sex_check rank score**: scenario tests only; no parity vs the documented
   XIST-vs-Y rank scoring rule.

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

- tests: 322 passed / 0 failed (00.1_test.txt, re-run on final code)
- R CMD check tarball: OK 0/0/0 (00.2)
- BiocCheck source + tarball: 1 ERROR (support-site email, environmental),
  1 WARNING (`set.seed`, justified), 12 NOTES; GitClone 0 ERROR / 1 WARNING
  (CITATION doi) (00.3, 00.5)
- `--as-cran`: 0 ERROR / 1 WARNING (qpdf, environmental) / 2 NOTES
  (pandoc README/NEWS; examples >5s = package-load overhead) (00.6)

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
(+30 expectations; suite now 352 passed / 0 failed).

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
