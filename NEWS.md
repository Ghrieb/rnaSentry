# rnaSentry 0.99.0 (development)

## New in this release

### Case studies and documentation

- New vignette `getting-started`: the exact CSV -> `SummarizedExperiment`
  bridge for your own data (read `counts.csv` + `clinical.csv`, intersect
  sample IDs, `load_counts()`, `design_audit()`, then `run_rnaSentry()`).
  README gains the same block under "Getting started with your own data".
- The pkgdown site (https://ghrieb.github.io/rnaSentry/) is rebuilt and
  deployed automatically on every push to `main` by
  `.github/workflows/pkgdown.yaml` (Bioc release in CI); `docs/` stays
  local-only and `.github/` is excluded from the source tarball. The
  case-study impact article is live on the site.
- New bundled case-study data: `inst/extdata/gse20685_case_study.rds`, a
  3000-gene x 327-sample subset of GSE20685 (Li et al., 2010; Affymetrix
  GPL570) with overall-survival metadata, built by the dev-only repro script
  `validation/repro_gse20685.R`.
- Three new offline-safe vignettes: `case-study-brca` (full pipeline on the
  bundled breast-cancer subset), `case-study-confounder-audit` (planted
  batch/region redundancy recovered as a minimal design), and
  `case-study-small-cohort` (the `events_per_parameter` guardrail and
  fold-level instability).
- The case-study vignettes were re-written as a three-tier naive-vs-guarded
  narrative (naive analysis -> rnaSentry standing guard -> counterfactual
  impact), and a landing page `case-study-impact` adds a master impact table
  across the four failure modes (wrong direction, confounded design,
  underpowered discovery, false transfer) plus the reproducible
  `C + C_rev = 1` invariant. The naive tiers are scripted inline (manual
  univariate screening, default-convention concordance, unadjusted design)
  and build fully offline.
- README now carries a status note, four case studies (including the LUAD
  discovery-to-external honest negative, GSE31210 -> GSE50081), a
  positioning/prior-art section (asuri, signifinder, SurvMarker, mRNAsi), and
  a Contributing pointer. Added `CONTRIBUTING.md` and `NEWS.md`.

### Engine reproducibility and check hygiene (2026-08-05)

- Reproducibility is now scoped with `withr::with_seed()` in
  `build_signature()` / `run_rnaSentry()`: the documented `seed` argument is
  still fully reproducible, but the caller's RNG state is left untouched
  (previously a `set.seed()` call leaked into the session).
- R dependency raised to `R (>= 4.5.0)` to track the Bioconductor 3.21 devel
  build.
- `build_signature()` and `run_rnaSentry()` examples are now self-contained
  and runnable (`\donttest`).

### Opt-in parallel cross-validation (2026-08-05)

- `build_signature()` and `run_rnaSentry()` gain a `BPPARAM` argument
  (Bioc-native, `BiocParallelParam`). Cross-validation stays **strictly
  serial by default** (`BPPARAM = NULL`); parallel evaluation is opt-in and
  only engaged when a `BiocParallelParam` object (e.g.
  `BiocParallel::SnowParam(2)`) is supplied.
- Reproducibility is preserved under parallelism: the documented `seed`
  continues to scope all fold partitioning inside `withr::with_seed()`, so
  per-fold evaluation is deterministic and a parallel run is **bit-identical**
  to the serial run (pinned by dedicated tests in `test-build_signature.R`).
- `BiocParallel` is a **Suggests-only** dependency: all calls are
  namespace-qualified behind a `requireNamespace()` guard, so the package
  works unchanged on installations without BiocParallel.

### Power analysis gate

- `validation/simulate_study.R` gains **Sim 6**, an external-transfer power
  analysis: two signal tiers calibrated to effective concordance 0.608 /
  0.654 show transfer power monotone in external event count, < 0.30 at ~35
  events (the small-cohort trap), and >= 0.80 at ~300 events for the moderate
  tier. The simulation suite is now 15 falsifiable targets, all passing.

### Earlier in the 0.99.0 cycle

- Guardrails: `possibly_log_scaled` (detects pre-logged data stored as
  `counts`), `events_per_parameter` (warns below ~5 events per signature
  gene), `load_counts()` intake audit (duplicate samples, non-integer counts,
  empty rows, Ensembl/non-syntactic IDs) with a flag ledger carried into the
  rendered report.
- Signature lock: `run_rnaSentry()` is single-use per session and refuses a
  second run until `lock_signature(sig, lock = FALSE)` releases the lock,
  preventing silent gene-set re-selection after survival analysis.
- External validation never recomputes the cutpoint from external data; the
  discovery cutpoint is applied unchanged.
- Validation dossier: statistical-parity suite (147 blocks / 457
  expectations), adversarial pass with regression tests, simulation study,
  real-data face validity, and a live GEO cross-cohort gate (GSE31210 ->
  GSE50081), including the documented CV-optimism (selection-leakage) finding.
- Sign-convention fix: concordance now consistently uses `reverse = TRUE`
  (higher score = higher hazard), pinned by parity tests and the
  `C + C_rev = 1` invariant.
