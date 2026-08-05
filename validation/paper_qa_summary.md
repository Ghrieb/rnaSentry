# rnaSentry — QA collection digest for the paper

Single source for the supplementary-methods numbers. Every figure below is
traceable to files under `validation/` and logs under `validation/logs/`.
Last updated 2026-08-05 (getting-started CSV->SE bridge round; also CI auto-deploy + site-sync).

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

- **147 test blocks / 457 expectations / 0 failures / 0 errors** (32
  warnings, all expected: the `events_per_parameter` guardrail firing on
  deliberately small synthetic fixtures).
- Layers:
  1. Unit tests per stage (`tests/testthat/test-*.R`).
  2. Statistical-parity tests: independent recomputation of published rules
     (e.g. `sex_check` XIST-vs-Y rank rule) pinned by `expect_identical`;
     serial↔parallel equivalence for `build_signature` cross-validation
     (2026-08-05, `BPPARAM` opt-in — parallel runs are bit-identical to
     serial).
  3. Adversarial pass (2026-08-03): ~90 degenerate inputs probed; 78+ rejected
     or flagged; 3 bug families + 4 soft-warnings found and fixed, all pinned
     by regression tests.
  4. Guardrail tests (2026-08-04): `possibly_log_scaled` and
     `events_per_parameter` fire/no-fire cases; report renders assumptions
     section.
- Run command (Windows/R 4.5.2): `Rscript <temp>/run_tests_parity.R` (the
  dev-only parity runner, same path as in `maintainer_testing_guide.md`) with
   `RSTUDIO_PANDOC` set. Logs: `validation/logs/00.1_test.txt`,
   `validation/logs/17_test.txt` (2026-08-05, post-BiocParallel round).
- Evidence: `validation/test_audit.md`, `validation/maintainer_testing_guide.md`.

## Simulation study (Gate 3)

- **15 falsifiable targets, 15/15 PASS.** (Sims 1-5: 11 targets as before;
  Sim 6 adds the external-transfer power analysis.)
- Engineered-signal fixtures with known survival association; null-data
  simulations for type-I-error behaviour. Concordance on nulls ~0.49-0.51
  (no winner's-curse inversion; see sign-convention note).
- Sim 6 (power): two signal tiers calibrated to effective C = 0.608 / 0.654
  transfer power rises with external events (weak 0.067 -> 0.417; moderate
  0.275 -> 0.842 over ~35 -> ~300 events); power < 0.30 at ~35 events
  (small-cohort trap) and >= 0.80 at ~300 events for the moderate tier.
- Evidence: `validation/simulation_study.md`, `validation/logs/03_simulation.txt`,
  `validation/logs/13_simulation.txt`, `validation/simulate_study.R`.

## Case-study vignettes and bundled data (Phase 3)

- `inst/extdata/gse20685_case_study.rds`: 3000 x 327 subset (top-3000 most
  variable genes) of GSE20685 with `time/event/age/subtype` metadata
  (83 events), built by `validation/repro_gse20685.R` from the GEO series
  matrix (`logcounts` = MAS5 log2 intensities; probe->gene collapse by largest
  mean).
- Four vignettes, all offline-safe (no network during build), all written as
  a three-tier naive-vs-guarded narrative (naive analysis -> rnaSentry
  standing guard -> counterfactual impact):
  - `case-study-impact.Rmd` — landing page with the master impact table
    across the four failure modes (wrong direction, confounded design,
    underpowered discovery, false transfer) and the reproducible
    `C + C_rev = 1` invariant (0.488 + 0.512 on the bundled subset).
  - `case-study-brca.Rmd` (pipeline on GSE20685 subset: CV C = 0.797 sd
    0.027, log-rank p = 2.99e-15, `design_audit` flags `subtype` —
    eta-squared 0.778 on the 3000-gene subset — and recommends
    `~ age + subtype`). The naive tier reproduces the inverted-convention
    read (default-convention concordance 0.316 vs 0.684 with
    `reverse = TRUE`, sum = 1) — the mechanism behind the historical
    0.217/0.783 face-validity episode.
  - `case-study-confounder-audit.Rmd` (synthetic batch/region redundancy:
    Cramer's V = 0.82, recommended `~ batch`): naive screening ranks the six
    planted batch-driven genes as the top "prognostic" hits; batch-adjusted
    screening drops all six.
  - `case-study-small-cohort.Rmd` (n = 30 / 23 events: guardrail fires,
    fold-level C ranges 0.50-1.00 while the mean reads 0.743).
- Gate 4's dataset label corrected: GSE20685 is **not** TCGA-BRCA (it is
  Li et al., 2010, Affymetrix GPL570); the old "TCGA-BRCA GSE20685" wording
  was removed from this file.

## Site, CI deployment, and cross-environment build (2026-08-05)

- The live site https://ghrieb.github.io/rnaSentry/ is served from the
  `gh-pages` branch and is rebuilt + deployed **automatically** on every push
  to `main` by `.github/workflows/pkgdown.yaml` (`r-lib` actions: `setup-r`
  with `bioc-version: release`, `setup-r-dependencies` with `needs: website`
  + `r-lib/pkgdown` + `local::.`, then `build_site_github_pages()` and
  `deploy_to_branch()`). `docs/` is gitignored (local preview only) and
  `^\.github$` is excluded from the source tarball (`.Rbuildignore`).
- Verified end-to-end 2026-08-05: workflow run #1 green (run
  31026245981); all four case-study vignettes + the main vignette compiled in
  a clean Ubuntu container on **R 4.6.1 / Bioc release** (SummarizedExperiment
  1.42, Biobase 2.72, BiocStyle 2.40) and produced the **same pinned numbers**
  as the local R 4.5.2 / Bioc 3.22 runs (0.316/0.684/1.0 convention read, CV
  C = 0.797, `C + C_rev = 1` = 0.488 + 0.512) — a cross-environment
  reproducibility data point. The getting-started round added a sixth vignette
  (CSV -> `SummarizedExperiment` bridge); run #3 (commit `11d575a` ->
  `gh-pages` `f9646cd`) compiled **all six vignettes** on the same Linux/Bioc
  release stack.
- Historical drift fixed: before this round the live site had been stale at
  gh-pages `46d9d06` (docs/ is gitignored, so a push to `main` alone did not
  publish). Manual `deploy_to_branch()` published `d36a20f` (case-study
  impact article went live), then CI published `3893b77`; future pushes sync
  the site with no manual step.
- Round record: `validation/logs/18_ci_site.txt`.

## Real-data face validity (Gate 4)

- Dataset: **GSE20685** (Li et al., 2010; Affymetrix GPL570, 327 primary
  breast tumors, 83 deaths). Not TCGA-BRCA.
- 20-gene signature splits risk groups with log-rank p = 3.8e-17.
- Cross-validated concordance **0.783** vs random-gene control **0.582**
  (6.3 SD above control mean, using the logged CV sd 0.032).
- Recovered expected biology: MKI67 high-risk / ESR1 high-risk hazard-ratio
  directions correct.
- Evidence: `validation/face_validity_review.md`,
  `validation/logs/04_face_validity.txt`, `validation/logs/face_validity_report.html`.

## Live GEO cross-cohort test (Gate 5b, Phase 3)

- Discovery GSE31210 (LUAD, n = 226, 35 deaths) → external GSE50081 (LUAD,
  n = 128, 52 deaths). Series matrices cached under `validation/cache/`; repro
  script `validation/repro_live_geo.R`, log `validation/logs/11_live_geo.txt`,
  rendered discovery report `validation/logs/live_geo_discovery_report.html`.
- Cross-cohort scoring required per-gene **z-scoring within each cohort**
  (raw/log2 scales differ across batches and the discovery-median cutpoint
  then leaves the external cohort with an empty risk group); probes collapsed
  to the **mean over all probes per gene** (a per-cohort best-probe rule
  selected different probes per cohort and biased the transfer).
- Result: discovery CV C = **0.862** vs a matched **permutation-null** of
  **0.836** (delta +0.026); external validation C = **0.540**, log-rank
  **p = 0.266**. Mechanical gates all PASS; scientific transfer **not
  demonstrated** on this low-power discovery (35 events) — reported as an
  honest negative.
- CV-optimism finding: gene selection on the full cohort before the CV split
  makes the CV concordance **screening-internal** — even survival-permuted
  data yields null CV C ≈ 0.8 (not 0.5) on this ~21k-gene panel. The tool's
  gates calibrate against that null, and the only fully out-of-sample
  estimate is `validate_external()`. Package docs (Rd / README / vignette)
  state this explicitly.

## R CMD check (Gate 5)

- `R CMD check --as-cran` on fresh tarball: **0 ERROR / 0 WARNING / 1 NOTE**
  (HTML `tidy` missing, external tool) since the 2026-08-05 qpdf CLI install;
  earlier gates carried an additional environmental `qpdf` WARNING, now
  resolved. No package issues.
- Reproduced across six consecutive gates: `logs/05_as_cran.txt`,
  `06_as_cran.txt`, `07_as_cran.txt`, `08_as_cran.txt`, `09_as_cran.txt`,
  `10_as_cran.txt` (plus earlier `00.*`), and `12_as_cran.txt` (2026-08-04,
  post CV-reframe doc round — 0 ERROR / 1 WARNING / 1 NOTE, both
  environmental).
- `logs/13_as_cran.txt` (2026-08-04, Phase-3 case-study/power round with the
  three new vignettes and `inst/extdata/gse20685_case_study.rds`): **0 ERROR /
  1 WARNING / 1 NOTE**, both environmental (`qpdf` WARNING, `tidy` NOTE).
- `logs/14_as_cran.txt` (2026-08-04, Batch-A ship-blocking round — LICENSE
  holder, report-template `results='asis'`, S4Vectors → Imports,
  xz-recompressed extdata): **0 ERROR / 1 WARNING / 1 NOTE**, both
  environmental.
- `logs/15_as_cran.txt` (2026-08-04, Batch C design-hardening round — targeted
  lock unlock, `min_events_per_parameter` forwarding, `concordance_na` flag,
  flags `stage` column docs): **0 ERROR / 1 WARNING / 1 NOTE**, both
  environmental.
- `logs/16_as_cran.txt` (2026-08-05, Phase-1 BiocCheck-prep round — withr-
  scoped seeding, `seq_len()`, runnable examples): **0 ERROR / 0 WARNING /
  1 NOTE** (`qpdf` installed; the residual HTML `tidy` NOTE environmental).
- BiocCheck (`logs/13_bioccheck.txt`, re-confirmed at `logs/14_bioccheck.txt`,
  `logs/15_bioccheck.txt`, and `logs/16_bioccheck.txt`):
  1 environmental ERROR (support-site
  email 404 — register
  `ghriebabdelkarimhani@gmail.com` on https://support.bioconductor.org),
  0 WARNING (the justified `set.seed` warning was resolved at the Phase-1
  BiocCheck-prep round by scoping reproducibility with `withr::with_seed()`),
  8 advisory NOTES (down from 13) — all explained in
  `logs/notes_documented.md`.
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
| `concordance_na` | warning | validate_external | external concordance non-estimable (`survival::concordance` returned non-finite); interpret validation cautiously |

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
