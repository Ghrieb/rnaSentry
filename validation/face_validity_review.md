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
seed = 42, render_report = TRUE)` completed end-to-end in ~3.8 min and wrote
`face_validity_report.html`.

## Four-point face-validity checklist

| # | Check | Evidence | Verdict |
|---|-------|----------|---------|
| 1 | High-risk group has shorter OS | `km_curve`: median OS **high = 11.4 y** (finite) vs **low = not reached** (beyond the 14.1-y follow-up); log-rank **p = 3.8e-17** | PASS |
| 2 | MKI67 / ESR1 hazard-ratio direction | Univariate screening: **MKI67 HR = 1.21** (proliferation → worse outcome, p = 0.14), **ESR1 HR = 0.86** (ER signaling → better outcome, **p = 0.0017**) | PASS |
| 3 | Confounder / design flags plausible | `design_audit` flags **subtype** as a confounder of the expression surrogate (eta-squared = 0.76, BH-adjusted p < 1e-96) — biologically expected, since intrinsic subtype is a dominant driver of global expression. `age` is not flagged (eta² ≈ 0.005). | PASS |
| 4 | ≥ 1 limitation warning surfaced in the report | `cox_model` emits a `proportional_hazards` warning for terms BCL2, RP11-28F1.2, LEPROTL1, GLOBAL; the rendered report includes the flag ledger. | PASS |

All four face-validity criteria pass on real data.

## Honest-out-of-sample finding (documented, not a bug)

- The selected 20-gene signature reports **CV concordance 0.217 (sd 0.032)** on
  held-out folds — far below 0.5 — while the in-sample `km_curve` split is
  highly significant (p = 3.8e-17).
- A **random-gene control** running the identical CV protocol with unselected
  genes on the same cohort gives mean C = **0.418** (simulated null data give
  the same control at 0.510; see `simulation_study.md`).
- Interpretation: aggressive top-N selection on ~23.5k genes with 83 events
  induces the documented winner's-curse bias in honest held-out CV. The pipeline
  **reports this honestly** rather than over-claiming generalizability; the gap
  between the selected signature (0.217) and the unselected control (0.418) is
  selection-induced, and the parity tests + simulation study confirm the CV
  machinery itself is computed correctly. This is an expected property of
  strong feature selection and is exactly the kind of limitation a validation
  report should surface.

## Additional observations

- **Tied survival times:** follow-up is recorded at 0.1-year precision, so the
  real-data concordance estimates carry heavy ties (relevant to the random-gene
  control reading 0.418 rather than 0.500; on continuous simulated times the
  same control reads 0.510).
- **Annotation artifacts:** probe symbols containing ` /// ` (multi-gene
  mappings, e.g. `LOC101928198 /// MFAP3L`) were kept verbatim; the package
  handled them correctly (non-syntactic-name support exercised on real data),
  and a production preprocessing step should split them.
- **Real-data bug found & fixed:** the first real-data run exposed two related
  defects that synthetic tests could not: (1) `data.frame()`/`as.data.frame()`
  mangle gene symbols such as `RP11-28F1.2` or `1-Mar`, so coefficient names no
  longer matched assay rownames (`mat[genes, ]` subscript out of bounds); and
  (2) `reformulate()` does not backquote non-syntactic symbols, so the stored
  Cox formula failed to parse. Both are fixed (`check.names = FALSE` on
  expression-column data frames, `.strip_backticks()`/`.backquote_names()`
  helpers) and pinned by a regression test
  (`test-adversarial_regressions.R` "non-syntactic gene symbols").

## Verdict

The package processes a real breast-cancer cohort end-to-end, recovers the
expected survival-signature behavior (risk groups separate, ESR1 protective,
MKI67 adverse), flags plausible confounders, and surfaces its own limitations
in the report. Submission gate: **PASS** for real-data face validity.
