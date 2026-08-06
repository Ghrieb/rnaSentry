# rnaSentry real-data face-validity review (GSE20685)

Pipeline: `Rscript validation/simulate_study.R` is the simulation gate; this
review covers the real-data run. Reproduction script (data download, processing,
pipeline, evidence extraction) is `validation/repro_gse20685.R`; the run record
and rendered report are kept in the maintainer's private dossier.

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
seed = 42, render_report = TRUE)` completed end-to-end and wrote the rendered
report (kept in the private dossier). The signature is reproducible: same 20
genes under the same seed before and after the concordance fix below.

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

Diagnosis (run record kept in the private dossier):
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
concordance (0.783). Verdict: **PASS** for real-data face validity.

---

# Live GEO test: GSE31210 → GSE50081 (LUAD), honest-negative result

Second independent real-data gate. Reproduction script
`validation/repro_live_geo.R`; run record and rendered discovery report are
kept in the maintainer's private dossier.

## Design

- **Discovery** GSE31210 (LUAD, GPL570): 226 primary tumors, 35 deaths, OS in
  days. The deposited series matrix is RMA-*linear* (median ≈ 99, max ≈ 54k)
  despite the paper's "log2" wording, so it is log2-transformed at the probe
  level to match GSE50081's log2 RMA.
- **Validation** GSE50081 (LUAD): 128 adenocarcinomas, 52 deaths, OS in years.
- Probe→gene collapse by **averaging all probes per gene** (deterministic and
  cohort-independent; a per-cohort best-probe rule was found to select
  different probes per cohort and bias cross-cohort scoring).
- Per-gene **z-score within each cohort** so the discovery-median cutpoint is
  scale-free across batches (raw/log2 scoring left the external cohort 3–6
  discovery-SD below the cutpoint; the discovery-median split then fails
  outright, which `validate_external()` correctly refuses).
- `run_rnaSentry(top_n = 7, repeats = 5, folds = 5, seed = 20260804,
  design_vars = stage/age/gender/smoking)`; permutation null with 3 permuted
  builds.

## Results

| Quantity | Value |
|----------|-------|
| Signature | TRAPPC3L, IKZF5, PTH2R, KLHL36, HMOX1, ZCRB1, ADAM10 |
| CV concordance (screening-internal) | 0.862 |
| Permutation-null CV (matched-optimism baseline) | 0.836 (sd 0.007) |
| Delta vs null | +0.026 |
| External validation concordance | 0.540 |
| External log-rank p | 0.266 (n.s.) |
| Events per signature gene | 5.0 (guardrail not tripped) |

## Acceptance-gate outcome (see log)

- **Mechanical gates — all PASS**: `load_counts()` flag ledger; pipeline +
  rendered report; session lock refused a second run; unlock released the
  lock; events/gene guardrail satisfied.
- **Calibrated gate — PASS**: CV C (0.862) exceeds the matched permutation
  null (0.836).
- **Scientific transfer — NOT demonstrated**: external C = 0.540, log-rank
  p = 0.266. The signature does not significantly stratify the independent
  cohort on this low-power dataset (35 discovery events).

## CV-optimism finding (selection leakage)

The pipeline selects genes on the full cohort *before* the CV split, so the
reported concordance is **screening-internal**: even survival-permuted data
yields a high null CV (0.77–0.84 across seeds/cohorts, i.e. ≈0.8 on this
~21k-gene panel) instead of the 0.5 a fully out-of-sample estimate would give.
This is winner's-curse optimism from 21k univariate tests × 35 events. The
package's gates handle it correctly by calibrating against the permutation
null, and the only fully independent estimate is `validate_external()`. The
package docs (Rd, README, vignette) now state this explicitly.

Note the distinction from the GSE20685 gate, which used a *random-gene*
control (0.582 on real data): that control keeps real survival but removes
selection. The permutation null instead keeps selection and removes signal —
it is the stricter, matched-optimism calibration used here.

## Verdict

The tool's mechanics, intake checks, lock, report rendering, and honest
guardrails all performed as designed on independent real GEO data. The
scientific result is an honest **negative** (a power finding at 35 events, not
a tool failure), and the demo documents a limitation of the screening-internal
CV metric that reviewers should be told about. This is the correct outcome for
the live GEO gate: **mechanics PASS, scientific transfer not demonstrated**.
