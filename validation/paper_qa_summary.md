# rnaSentry — QA collection digest for the paper

Single source for the supplementary-methods numbers. Every figure below is
traceable to files under `validation/` and logs under `validation/logs/`.
Last updated 2026-08-04 (commit `b161f3a`).

## Method summary

- Pipeline stages: data intake and validation (`load_counts`), sample QC and
  identity audit (`qc_explore`, `sex_check`, `pca_audit`), confounder and
  design-formula detection (`design_audit`), signature construction
  (`build_signature`, `lock_signature`), survival modeling (`km_curve`,
  `cox_model`, `survival_parametric`), external validation
  (`validate_external`), and an HTML report with a full audit trail
  (`generate_report`).
- Signature construction: univariate Cox screening followed by a joint
  multivariate Cox fit producing a weighted linear score; stability assessed
  by repeated stratified cross-validation. Higher score = higher hazard.
- Concordance convention: `C` computed with `reverse = TRUE` so values above
  0.5 indicate correct discrimination (this convention is the one validated
  against biology; see `test_audit.md`, sign-convention section).
- Guardrail model: every automated statistical decision is reported with its
  supporting statistic and a plain-language rationale; validation stages are
  gated rather than silently degraded (flag ledger surfaced in the report).
- Signature locking: `lock_signature()` records a session-level fingerprint
  that prevents `run_rnaSentry()` from being called again until the
  signature is explicitly unlocked, guarding against silent gene-set
  re-selection after survival analysis.

## Statistical routing in design_audit()

The design audit performs three scans, each routing between tests based on
variable type and distributional properties. Every raw statistic is returned
in the result object; the print method shows aggregate counts, not routing
decisions.

**Confounder scan** (each candidate design variable vs. the surrogate PC):
- Numeric variable: `lm(PC ~ variable)`, F-test, R-squared effect size.
- Categorical variable, 2 levels: `aov(PC ~ group)`, F-test, eta-squared.
- Categorical variable, >2 levels, residuals approximately normal:
  `aov(PC ~ group)`, F-test, eta-squared. Normality tested via
  Shapiro-Wilk on group-centred residuals (`.is_approx_normal()`).
- Categorical variable, >2 levels, residuals not normal:
  `kruskal.test(PC, group)`, H-statistic, epsilon-squared.
- P-values are Benjamini-Hochberg adjusted across all tested variables.

**Redundancy scan** (pairwise variable associations):
- Categorical x categorical: `chisq.test()` + Cramer's V.
- Numeric x numeric: `cor.test(method = "pearson")`, Pearson r.
- Mixed (one of each): reuses the ANOVA/Kruskal-Wallis router above.
- A pair is flagged redundant when `p < alpha` AND the effect size exceeds
  `redundant_effect_size` (default 0.8).

**Interaction scan** (surrogate PC x design variable):
- Nested linear models: `lm(PC ~ v + s)` vs. `lm(PC ~ v * s)`.
- F-test via `anova(reduced, full)`.
- Uses a separate `interaction_alpha` threshold.
- Skipped when fewer than 2 PCs or fewer than 6 samples are available.

## Test suite (Gate 1)

- **142 test blocks / 441 expectations / 0 failures / 0 errors** (32
  warnings, all expected: the `events_per_parameter` guardrail firing on
  deliberately small synthetic fixtures).
- Layers:
  1. Unit tests per stage (`tests/testthat/test-*.R`).
  2. Statistical-parity tests: independent recomputation of published rules
     (e.g. `sex_check` XIST-vs-Y rank rule) pinned by `expect_identical`.
  3. Adversarial pass (2026-08-03): ~90 degenerate inputs probed; 78+ rejected
     or flagged; 3 bug families + 4 soft-warnings found and fixed, all pinned
     by regression tests.
  4. Guardrail tests (2026-08-04): `possibly_log_scaled` and
     `events_per_parameter` fire/no-fire cases; report renders assumptions
     section.
- Run command (Windows/R 4.5.2): `Rscript run_tests_parity.R` with
  `RSTUDIO_PANDOC` set. Log: `validation/logs/00.1_test.txt`.
- Evidence: `validation/test_audit.md`, `validation/maintainer_testing_guide.md`.

## Simulation study (Gate 3)

- **11 falsifiable targets, 11/11 PASS.**
- Engineered-signal fixtures with known survival association; null-data
  simulations for type-I-error behaviour. Concordance on nulls ~0.49-0.51
  (no winner's-curse inversion; see sign-convention note).
- Evidence: `validation/simulation_study.md`, `validation/logs/03_simulation.txt`,
  `validation/simulate_study.R`.

## Real-data face validity (Gate 4)

- Dataset: TCGA-BRCA GSE20685 (bulk RNA-seq, survival).
- 20-gene signature splits risk groups with log-rank p = 3.8e-17.
- Cross-validated concordance **0.783** vs random-gene control **0.582**
  (8.3 SD above control mean).
- Recovered expected biology: MKI67 high-risk / ESR1 high-risk hazard-ratio
  directions correct.
- Evidence: `validation/face_validity_review.md`,
  `validation/logs/04_face_validity.txt`, `validation/logs/face_validity_report.html`.

## R CMD check (Gate 5)

- `R CMD check --as-cran` on fresh tarball: **0 ERROR / 1 WARNING / 1-2 NOTEs**,
  all environmental (`qpdf` missing; HTML `tidy` missing; occasionally
  "unable to verify current time"), no package issues.
- Reproduced across six consecutive gates: `logs/05_as_cran.txt`,
  `06_as_cran.txt`, `07_as_cran.txt`, `08_as_cran.txt`, `09_as_cran.txt`,
  `10_as_cran.txt` (plus earlier `00.*`).
- BiocCheck equivalent run locally via `validation/docker_check.sh` when a
  Docker-capable machine is available (see `validation/cross_platform.md`).

## Guardrail catalog (flag ledger)

All flags carry `(check, severity, detail, stage)` and are rendered into the
report's audit trail. Severities: critical / warning / info.

| Check | Severity | Stage | Meaning |
|---|---|---|---|
| `duplicate_samples` | critical | load_counts / qc_explore | duplicated sample identifiers |
| `missing_counts` | critical | qc_explore | samples/genes with no counts |
| `non_integer_counts` | warning | load_counts / qc_explore | counts assay not integer |
| `empty_rows` | warning | load_counts | gene(s) have all-zero or all-NA counts |
| `ensembl_ids_detected` | info | load_counts | more than half of gene IDs look like Ensembl accessions |
| `non_syntactic_ids` | info | load_counts | gene IDs contain non-syntactic characters |
| `missing_metadata` | warning | qc_explore | empty required metadata columns |
| `library_size_outlier` | warning | qc_explore | extreme library sizes |
| `signature_locked` | info | lock_signature | signature locked against downstream mutation |
| `signature_unlocked` | info | lock_signature | signature unlocked; downstream stages may be rerun |
| `surrogate_pc_clamped` | info | design_audit | surrogate-variable PC capped |
| `variable_no_data` | warning | design_audit | design variable has no usable data |
| `redundant_variable` | warning | design_audit | collinear/duplicate design terms |
| `variable_constant` | warning | design_audit | design variable has no variation |
| `possibly_log_scaled` | warning | pca_audit / build_signature / validate_external | assay looks pre-log-transformed (non-integer, max < 40) and would be double-logged |
| `pca_gene_filter` | warning | pca_audit | genes dropped before PCA |
| `batch_single_level` | warning | pca_audit | batch column has a single level |
| `batch_associated_pc` | warning | pca_audit | top PC associates with batch |
| `design_terms_missing` | info | build_signature | requested design term absent |
| `adjusted_screening` | info | build_signature | screening adjusted for design |
| `unadjusted_screening` | warning | build_signature | design adjustment bypassed |
| `gene_filter` | warning | build_signature | genes filtered pre-model |
| `gene_name_collision` | warning | build_signature | duplicate/non-syntactic gene symbols |
| `fewer_genes_than_requested` | info | build_signature | requested top_n reduced |
| `coefficient_unstable` | warning | build_signature | coefficient near-non-estimable |
| `events_per_parameter` | warning | build_signature | n_events / n_genes < min_events_per_parameter (default 5) |
| `sparse_events` | warning | km_curve / validate_external | sparse risk-group events |
| `missing_covariates` | info | cox_model | covariates dropped |
| `non_estimable_coefficients` | warning | cox_model | coefficients non-estimable |
| `non_estimable_genes` | warning | cox_model | signature genes non-estimable in joint fit |
| `ph_test_failed` | warning | cox_model | proportional-hazards test error |
| `proportional_hazards` | warning | cox_model | PH assumption violated (cox.zph) |
| `model_fit_failed` | warning | survival_parametric | parametric fit did not converge |
| `aic_models_close` | info | survival_parametric | competing models within AIC |
| `missing_genes_dropped` | warning | validate_external | external cohort lacks signature genes |

## How to cite the QA evidence in the paper

- Software: cite the Bioconductor package page and DOI once accepted
  (`10.18129/B9.bioc.rnaSentry` pattern); until then cite the repo
  `https://github.com/Ghrieb/rnaSentry`.
- Supplementary methods: reference this file and `validation/test_audit.md`,
  `simulation_study.md`, `face_validity_review.md` with the commit hash.
- Methodological references for guardrails: Peduzzi et al., *J Clin Epidemiol*
  1996 (10 events-per-variable rule); Vittinghoff & McCulloch, *Am J Epidemiol*
  2007 (5-9 EPV adequate in some settings); Grambsch & Therneau for the
  scaled Schoenfeld residual PH test (cox.zph).
