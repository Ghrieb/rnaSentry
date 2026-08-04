# rnaSentry real-data face-validity review (GSE20685)

Pipeline: `Rscript validation/simulate_study.R` is the simulation gate; this
review covers the real-data run. Reproduction script (data download, processing,
pipeline, evidence extraction) and full log are in
`validation/logs/04_face_validity.txt`; the rendered report is
`validation/logs/face_validity_report.html`.

## Dataset

- **GSE20685** (Li et al., 2010), Affymetrix GPL570, 327 primary breast-cancer
  samples, overall-survival follow-up.
- Preprocessing: probe-level MAS5 log2 intensities (already normalized; the
  assay is stored as `logcounts` so the pipeline does not log-transform again),
  probes mapped to gene symbols via the platform annotation, collapsed to the
  highest-mean probe per symbol → **23,520 genes × 327 samples, 83 deaths**.
- colData: `time` = follow-up duration in years (0.4–14.1), `event` = death
  indicator, `age`, `subtype` (paper's intrinsic subtypes I–VI).

## Pipeline

`run_rnaSentry(se, "time", "event", outcome_col = "overall_survival",
design_vars = c("age", "subtype"), top_n = 20, repeats = 3, folds = 3,
seed = 42, render_report = TRUE)` completed end-to-end and wrote
`face_validity_report.html`. The signature is reproducible: same 20 genes under
the same seed before and after the concordance fix below.

## Four-point face-validity checklist

| # | Check | Evidence | Verdict |
|---|-------|----------|---------|
| 1 | High-risk group has shorter OS | `km_curve`: median OS **high = 11.4 y** (finite) vs **low = not reached** (beyond the 14.1-y follow-up); log-rank **p = 3.8e-17** | PASS |
| 2 | MKI67 / ESR1 hazard-ratio direction | Univariate screening: **MKI67 HR = 1.21** (proliferation → worse outcome, p = 0.14), **ESR1 HR = 0.86** (ER signaling → better outcome, **p = 0.0017**) | PASS |
| 3 | Confounder / design flags plausible | `design_audit` flags **subtype** as a confounder of the expression surrogate (eta-squared = 0.76, BH-adjusted p < 1e-96) — biologically expected, since intrinsic subtype is a dominant driver of global expression. `age` is not flagged (eta² ≈ 0.005). | PASS |
| 4 | ≥ 1 limitation warning surfaced in the report | `cox_model` emits a `proportional_hazards` warning for terms BCL2, RP11-28F1.2, LEPROTL1, GLOBAL; the rendered report includes the flag ledger. | PASS |

All four face-validity criteria pass on real data.

## Out-of-sample concordance (sign-convention bug found and fixed)

An early run reported **CV concordance 0.217 (sd 0.032)** and the review
initially rationalized it as extreme selection-induced winner's curse. A
methodological review challenged that: a C-index ~9 SD below chance with a
tight sd is not noise, `1 − 0.217 = 0.783` is exactly the biologically expected
value, and the parity tests could not catch a sign error because they replayed
the same call.

Diagnosis (recorded in `validation/logs/04_face_validity.txt`):
`survival::concordance(Surv ~ x)` defaults to **"larger x ⇒ longer survival"**
(reverse = FALSE). A Cox risk score (larger = higher hazard = shorter survival)
must be passed with `reverse = TRUE`. Both risk-score call sites in the package
used the default, so they reported `1 − Harrell's C`:

| Quantity | Reported (bug) | Correct |
|----------|----------------|---------|
| CV concordance, GSE20685 | 0.217 | **0.783** (sd 0.032) |
| representative fold (r1,f1) | 0.190 | 0.810 |
| in-sample full-model concordance | — | 0.840 |
| random-gene control (real data) | 0.418 | 0.582 |

The identity `C(default) + C(reverse=TRUE) = 1` held exactly on every fold, and
the held-out CV (0.783) sits below the in-sample value (0.840) by a modest,
healthy amount — a strong and directionally correct out-of-sample result, not
a 0.22 collapse. The same bug affected `validate_external()` (its reported
0.160 was really 0.840); `cox_model()` was already correct (it uses the coxph
method's own concordance). Both call sites now pass `reverse = TRUE`, and new
direction-sensitive regression tests (CV mean > 0.5 on a signal-bearing
fixture; `validate_external` concordance > 0.5 on a cohort whose survival is
engineered to follow the signature score; the `C + C_rev = 1` invariant) would
have caught the original error.

## Interpretation of the corrected numbers

- The selected 20-gene signature has honest held-out discrimination of
  **0.783** on GSE20685.
- The random-gene control under the same protocol reads **0.582** on real data
  (random real genes carry weak signal via shared biology / abundance
  structure); on simulated null data the same control reads **0.490** (the
  calibration reference). The signature's 0.783 is far above both → genuine
  held-out risk signal.
- No winner's-curse inversion is present; the earlier "poor generalization"
  reading of 0.217 was an artifact of the inverted concordance convention.

## Additional observations

- **Tied survival times:** follow-up is recorded at 0.1-year precision, so
  concordance estimates on this cohort carry ties; this affects variance but
  not the direction of the statistic.
- **Annotation artifacts:** probe symbols containing ` /// ` (multi-gene
  mappings, e.g. `LOC101928198 /// MFAP3L`) were kept verbatim; the package
  handled them correctly (non-syntactic-name support exercised on real data),
  and a production preprocessing step should split them.
- **Real-data bugs found & fixed:** (1) non-syntactic gene symbols
  (`RP11-28F1.2`, `1-Mar`) broke data-frame construction and formula parsing —
  fixed with `check.names = FALSE` and `.strip_backticks()`/`.backquote_names()`
  helpers, pinned by a regression test; (2) the concordance sign-convention bug
  described above, pinned by direction-sensitive tests.

## Known limitation: events per parameter in the discovery cohort

GSE20685 has 83 deaths and the selected signature has 20 genes — about
**4.2 events per model parameter**, below the ~10 events-per-parameter rule of
thumb commonly cited for stable Cox-model estimation. The held-out CV
concordance (0.783) far exceeds the random-gene control (0.582), and the
risk-group split is highly significant (log-rank p = 3.8e-17), which together
support genuine held-out signal; but the parameterization is intentionally
small (`top_n = 20`), and the per-gene hazard ratios should be interpreted with
this event count in mind. A larger discovery cohort — or a regularized /
smaller signature — would tighten per-parameter stability. This is a
face-validity gate, not a claim of clinical utility.

## Verdict

The package processes a real breast-cancer cohort end-to-end, recovers the
expected survival-signature behavior (risk groups separate, ESR1 protective,
MKI67 adverse), flags plausible confounders, surfaces its own limitations in
the report, and now reports a directionally correct, strong held-out
concordance (0.783). Submission gate: **PASS** for real-data face validity.
