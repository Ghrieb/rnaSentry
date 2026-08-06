# Changelog

## rnaSentry 0.99.0 (development)

### New in this release

#### Case studies and documentation

- New vignette `getting-started`: the exact CSV -\>
  `SummarizedExperiment` bridge for your own data (read `counts.csv` +
  `clinical.csv`, intersect sample IDs,
  [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md),
  [`design_audit()`](https://ghrieb.github.io/rnaSentry/reference/design_audit.md),
  then
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)).
  README gains the same block under “Getting started with your own
  data”.
- The pkgdown site (<https://ghrieb.github.io/rnaSentry/>) is deployed
  automatically on every push to `main`; the case-study vignettes are
  live on the site.
- New bundled case-study data: `inst/extdata/gse20685_case_study.rds`, a
  3000-gene x 327-sample subset of GSE20685 (Li et al., 2010; Affymetrix
  GPL570) with overall-survival metadata, built via
  `validation/repro_gse20685.R` in the GitHub repository.
- Three new offline-safe vignettes: `case-study-brca` (full pipeline on
  the bundled breast-cancer subset), `case-study-confounder-audit`
  (planted batch/region redundancy recovered as a minimal design), and
  `case-study-small-cohort` (the `events_per_parameter` guardrail and
  fold-level instability).
- The four pipeline vignettes now carry figures, all generated from
  bundled offline data at build time (no network access needed):
  Kaplan-Meier risk groups, an adjusted Cox hazard-ratio forest, a
  parametric survival overlay, and a subtype-colored PCA audit scatter
  in `case-study-brca`; batch- and region-colored PC1-PC2 scatters in
  `case-study-confounder-audit`; an external-cohort Kaplan-Meier split
  in the main vignette; and a per-fold cross-validated concordance
  stripchart in `case-study-small-cohort`.
- The case-study vignettes were re-written as a three-tier
  naive-vs-guarded narrative (naive analysis -\> rnaSentry standing
  guard -\> counterfactual impact), and a landing page
  `case-study-impact` adds a master impact table across the four failure
  modes (wrong direction, confounded design, underpowered discovery,
  false transfer) plus the reproducible `C + C_rev = 1` invariant. The
  naive tiers are scripted inline (manual univariate screening,
  default-convention concordance, unadjusted design) and build fully
  offline.
- README now carries a status note, four case studies (including the
  LUAD discovery-to-external honest negative, GSE31210 -\> GSE50081), a
  positioning/prior-art section (asuri, signifinder, SurvMarker,
  mRNAsi), and a Contributing pointer. Added `CONTRIBUTING.md` and
  `NEWS.md`.
- README gains an “Upcoming case studies” section listing three planned
  case studies: a clinical-covariate trap on the bundled GSE20685
  subset, a batch catastrophe on the GSE31210 + GSE30219
  identical-GPL570 merge, and a cross-histology LUAD -\> LUSC transfer
  on GSE30219. They will ship as new vignettes once the analyses are
  finalized.

#### Engine reproducibility

- Reproducibility is now scoped with
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html)
  in
  [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  /
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md):
  the documented `seed` argument is still fully reproducible, but the
  caller’s RNG state is left untouched (previously a
  [`set.seed()`](https://rdrr.io/r/base/Random.html) call leaked into
  the session).
- R dependency raised to `R (>= 4.5.0)` to track the Bioconductor 3.21
  devel build.
- [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  and
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  examples are now self-contained and runnable (`\donttest`).

#### Opt-in parallel cross-validation

- [`build_signature()`](https://ghrieb.github.io/rnaSentry/reference/build_signature.md)
  and
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  gain a `BPPARAM` argument (Bioc-native, `BiocParallelParam`).
  Cross-validation stays **strictly serial by default**
  (`BPPARAM = NULL`); parallel evaluation is opt-in and only engaged
  when a `BiocParallelParam` object (e.g. `BiocParallel::SnowParam(2)`)
  is supplied.
- Reproducibility is preserved under parallelism: the documented `seed`
  continues to scope all fold partitioning inside
  [`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
  so per-fold evaluation is deterministic and a parallel run is
  **bit-identical** to the serial run (pinned by dedicated tests in
  `test-build_signature.R`).
- `BiocParallel` is a **Suggests-only** dependency: all calls are
  namespace-qualified behind a
  [`requireNamespace()`](https://rdrr.io/r/base/ns-load.html) guard, so
  the package works unchanged on installations without BiocParallel.

#### Power analysis

- The simulation study (`validation/simulate_study.R`) gains a power
  analysis of external-cohort transfer: two signal tiers calibrated to
  effective concordance 0.608 / 0.654 show transfer power monotone in
  external event count, \< 0.30 at ~35 events (the small-cohort trap),
  and \>= 0.80 at ~300 events for the moderate tier. All 15 simulation
  targets pass.

#### Earlier in the 0.99.0 cycle

- Guardrails: `possibly_log_scaled` (detects pre-logged data stored as
  `counts`), `events_per_parameter` (warns below ~5 events per signature
  gene),
  [`load_counts()`](https://ghrieb.github.io/rnaSentry/reference/load_counts.md)
  intake audit (duplicate samples, non-integer counts, empty rows,
  Ensembl/non-syntactic IDs) with a flag ledger carried into the
  rendered report.
- Signature lock:
  [`run_rnaSentry()`](https://ghrieb.github.io/rnaSentry/reference/run_rnaSentry.md)
  is single-use per session and refuses a second run until
  `lock_signature(sig, lock = FALSE)` releases the lock, preventing
  silent gene-set re-selection after survival analysis.
- External validation never recomputes the cutpoint from external data;
  the discovery cutpoint is applied unchanged.
- Validation: a statistical-parity suite (147 blocks / 457
  expectations), regression tests, the simulation study
  (`validation/simulate_study.R`), real-data face validity
  (`validation/repro_gse20685.R`), and a live GEO cross-cohort run
  (GSE31210 -\> GSE50081, `validation/repro_live_geo.R`), including the
  documented CV-optimism (selection-leakage) finding.
- Sign-convention fix: concordance now consistently uses
  `reverse = TRUE` (higher score = higher hazard), pinned by parity
  tests and the `C + C_rev = 1` invariant.
