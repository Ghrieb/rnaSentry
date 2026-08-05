# rnaSentry — Plain-Language Project Description

*A feature-by-feature overview of the `rnaSentry` R package (v0.99.0), written so
the whole codebase is understandable at a glance.*

---

## 1. What rnaSentry is

`rnaSentry` is an R package that automates the statistical pipeline for turning a
bulk RNA-seq count matrix into a **validated prognostic gene signature**: a small
set of genes whose weighted expression predicts survival (risk of an event such
as death or recurrence).

It takes you from raw data to a written report in one reproducible flow:

```
raw counts + sample metadata
   → QC and sample-identity audits
   → confounder / batch detection
   → signature construction with cross-validation
   → survival modeling (Kaplan-Meier, Cox, parametric)
   → external validation on an independent cohort
   → self-contained HTML report
```

Two words capture the package's philosophy: **auditable** and **gated**.

- **Auditable.** Every automated decision is recorded in a per-stage *flag
  ledger* — rows of `(check, severity, detail, stage)` — that is carried into
  the final HTML report's audit trail. Nothing that changes an analysis is done
  silently.
- **Gated.** Pipeline stages do not quietly degrade when something is off. When
  data quality is poor or an assumption is violated, the stage raises a warning
  flag (or, for truly fatal issues, refuses to proceed) instead of emitting a
  number that looks fine but is meaningless.

The package is deliberately narrow: **bulk** RNA-seq, **right-censored survival**
outcomes, **linear Cox risk scores**. Single-cell data, binary/continuous
endpoints, competing risks, and batch correction are explicitly out of scope
(see §7).

---

## 2. The pipeline in one run

`run_rnaSentry(se, time_col, event_col, ...)` is the one-call entry point. It
runs the stages in order on a single cohort:

```
design_audit → pca_audit → build_signature → (auto lock_signature)
   → km_curve → cox_model → survival_parametric → generate_report (HTML)
```

It returns a `run_rnaSentry` object whose `stages` list holds each stage's full
result, including the `flags` data frame. The pipeline is **single-use per R
session**: once a signature is discovered it is locked, and a second
`run_rnaSentry()` call is refused until the lock is released (see §4).

---

## 3. Feature-by-feature

Each exported function is described as **input → what it does → output →
guardrails** (the flags and gates that protect you).

### 3.1 `load_counts(counts, colData, assay_name = "counts", ...)`

- **Input:** a raw count matrix and its sample metadata.
- **What it does:** validates the structure and builds a `SummarizedExperiment`.
  It does *not* normalize — it is pure intake.
- **Output:** a `SummarizedExperiment` with a flag ledger stored in its metadata.
- **Guardrails:** records flags for duplicate sample IDs (the duplicate is
  dropped — critical), non-integer count values, empty gene rows, and
  non-syntactic or Ensembl-style gene IDs. Issues are surfaced up front and
  carried into the report's audit trail.

### 3.2 `qc_explore(se)`

- **Input:** a `SummarizedExperiment` of raw counts.
- **What it does:** basic sample-level quality control audit.
- **Guardrails:** checks duplicate sample identifiers, missing sample metadata,
  non-integer count values, and library-size outliers (median-absolute-deviation
  based, default threshold 3 MAD).

### 3.3 `pca_audit(se, batch_cols, ...)` and `plot_pca_audit()`

- **Input:** the expression object plus candidate batch/technical columns.
- **What it does:** runs PCA on the log-scale expression and tests each leading
  PC against each candidate batch variable to find technical structure hiding in
  the data.
- **Output:** a `rnaSentry_pca` object with per-PC batch test results; the
  companion `plot` method visualizes the PC loadings/clustering.
- **Guardrails:** prefers an assay named `logcounts` or `vst`; otherwise it
  log2-transforms counts internally (`log2(counts + 1)`). Genes with zero
  variance or missing values are dropped **with a report**, never silently.
  Batch tests report both raw and BH-adjusted p-values but flag on the **raw**
  p-value — intentional, because only a handful of PCs are tested against one
  variable at a time, so multiple-testing correction would be a near no-op.

### 3.4 `design_audit(se, design_vars, ...)`

- **Input:** the expression object and a set of candidate design/confounder
  columns (batch, tissue, platform, etc.).
- **What it does:** a confounder sweep. It builds a leading-PC surrogate for
  expression and tests each candidate confounder against it using the
  appropriate test per variable type (ANOVA / Kruskal-Wallis / linear
  regression) with effect sizes (eta, epsilon, R-squared).
- **Output:** a `rnaSentry_design` object ranking candidate confounders.
- **Guardrails:** in contrast to `pca_audit`, it flags candidate confounders on
  **BH-adjusted** p-values — intentional, because it scans many candidate
  variables at once and needs the correction to avoid spurious flags.

### 3.5 `sex_check(se, sex_col = "sex", xist_gene = "XIST", y_genes, ...)`

- **Input:** the expression object with a reported-sex column.
- **What it does:** verifies biological sex by comparing reported sex against
  sex inferred from **XIST** (female-marker) versus **Y-chromosome marker**
  expression (defaults `RPS4Y1`, `DDX3Y`, `KDM5D`), using a rank-based score.
- **Output:** per-sample calls of female / male / ambiguous, with the basis for
  each call.
- **Purpose:** catches sample swaps or metadata errors before they propagate
  into survival analysis.

### 3.6 `build_signature(se, time_col, event_col, top_n, repeats, folds, seed, ...)`

The core discovery stage.

- **Input:** the expression object and survival columns, plus selection and
  cross-validation parameters.
- **What it does, step by step:**
  1. Screens every variable gene with a univariate Cox model, recording hazard
     ratio, Wald p-value, and BH-adjusted p-value. Screening can be *adjusted*
     for design terms supplied via `design_terms` (default), so genes whose
     survival association is driven entirely by a confounder (e.g. batch) are
     not selected.
  2. Selects the signature genes by `top_n` or by a p-value cutoff.
  3. Refits a joint multivariate Cox model on the selected genes; if any
     coefficient is non-finite, it is dropped **with an explicit flag**.
  4. Repeats an **event-stratified k-fold cross-validation** and reports
     Harrell's C mean and SD (`reverse = TRUE`: higher risk score = earlier
     event).
- **Output:** an `rnaSentry_signature` object holding the gene list, per-gene
  Cox statistics, final coefficients, CV results, and the selection record.
- **Guardrails — the screening/weight asymmetry:** genes are screened adjusted
  for design terms, but the final signature **weights are estimated from a
  gene-only model**. This is intentional: the risk score must stay portable to
  external cohorts that don't record the discovery cohort's design covariates.
  Adjusted inference on the chosen genes is provided separately by `cox_model()`.
- **Guardrails — event count per parameter:** if the event count per signature
  gene falls below `min_events_per_parameter` (default 5), it raises the
  `events_per_parameter` warning. A flagged signature should be treated as
  hypothesis-generating, not confirmatory.
- **Known caveat:** CV concordance here is *screening-internal* (gene selection
  happens before the CV split), so it is optimistic — on a ~21k-gene panel even
  survival-permuted data can show a high null (~0.8). Read it relative to a
  matched null; the honest out-of-sample number is `validate_external()`.

### 3.7 `lock_signature(sig, lock = TRUE)`

- **Input:** a signature object; `lock = FALSE` releases.
- **What it does:** toggles the reproducibility lock on a discovered signature.
  `run_rnaSentry()` locks automatically after `build_signature()`; later stages
  record the lock state.
- **Output:** the signature object with `locked`/`lock_time` updated.
- **Guardrails:** while any locked signature exists, `run_rnaSentry()` refuses
  to run again — preventing silent gene-set re-selection after survival
  analysis. To analyze a second cohort in the same session, unlock the first.

### 3.8 `km_curve(sig, se, time_col, event_col, cutpoint = "median", ...)`

- **Input:** a signature and a cohort to score.
- **What it does:** scores each sample with the signature's joint Cox
  coefficients, splits samples into low/high risk groups (median or supplied
  cutpoint), and fits Kaplan-Meier survival curves with a log-rank test.
- **Output:** an `rnaSentry_km` object with the KM fit, per-group medians and
  event counts, the log-rank p-value, and the **cutpoint used**.
- **Guardrails:** the cutpoint chosen here is recorded on the object because it
  is the value `validate_external()` must reuse (it never re-derives one).
  Sparse events in a risk group raise a `sparse_events` warning.

### 3.9 `cox_model(sig, se, time_col, event_col, ...)`

- **Input:** a signature and the cohort.
- **What it does:** fits the adjusted multivariate Cox model — signature genes +
  the recorded design terms + any supplied confounders — and runs
  proportional-hazards diagnostics (`cox.zph`) per gene.
- **Output:** an `rnaSentry_cox_model` object (with `print`/`plot` methods)
  containing the fit, per-gene HRs with BH-adjusted p-values, and the PH check.
- **Guardrails:** a potential proportional-hazards violation raises the
  `proportional_hazards` flag, which the report tells you to interpret
  cautiously.

### 3.10 `survival_parametric(sig, se, time_col, event_col, dists, ...)`

- **Input:** a signature and the cohort.
- **What it does:** fits parametric survival models (`survreg`) — Weibull,
  exponential, lognormal, loglogistic — and compares them by AIC.
- **Output:** an `rnaSentry_parametric` object with per-distribution fits, the
  AIC comparison, and a plot method overlaying the fitted survival curves.

### 3.11 `validate_external(sig, external_se, cutpoint, ...)`

The honest portability check.

- **Input:** the discovery signature, an **independent** cohort, and the
  discovery cutpoint.
- **What it does:** re-scores the external samples against the signature and the
  *given* cutpoint, then reports log-rank test, per-group median survival, event
  counts, and concordance — **without re-estimating anything** from the external
  data.
- **Guardrails — strict by default:** it requires the cutpoint to be passed
  explicitly (it never recomputes one), and requires full gene overlap between
  signature and external cohort by default (`drop_missing = FALSE`,
  `min_gene_overlap_frac = 1`). If concordance is non-finite, it raises a
  `concordance_na` warning flag rather than reporting a misleading number.

### 3.12 `generate_report(run, output_file, ...)`

- **Input:** a `run` object from `run_rnaSentry()`.
- **What it does:** renders the shipped `report_template.Rmd` into a
  **self-contained HTML** report.
- **Output:** the report file, containing a pipeline summary, the assumptions
  and scope section, a **unified audit-trail table** built from every stage's
  flag ledger, and a section printing each stage's full result.

### 3.13 S3 methods

Each stage object has a tailored `print` method (human-readable summary) and,
where relevant, a `plot` method (`plot.rnaSentry_cox_model`,
`plot.rnaSentry_external`, `plot.rnaSentry_km`, `plot.rnaSentry_parametric`).
`plot_pca_audit` is exported as a standalone function for the PCA stage.

---

## 4. Cross-cutting architecture

Two design elements run through every stage and are worth understanding on their
own.

### The flag ledger (`utils.R`)

A shared internal API (`.new_flags` / `.add_flag`) appends rows to each stage's
`flags` data frame: `(check, severity, detail, stage)`. The stage is attached as
an attribute until the first append. This single mechanism is what makes every
automated decision auditable — and what `generate_report()` flattens into the
report's audit-trail table. Tests assert both the flag *names* and the
*conditions* that trigger them (e.g. `events_per_parameter`, `concordance_na`,
`sparse_events`, `proportional_hazards`, `possibly_log_scaled`).

### The discovery lock (`.rnaSentry_locked_sigs`)

A session-level environment records discovered signatures. `run_rnaSentry()`
locks automatically after discovery and refuses a second run while any locked
signature exists; `lock_signature(sig, lock = FALSE)` releases. This is the
reproducibility guardrail against silently re-selecting signature genes after
survival analysis — the classic way to overfit a prognostic paper.

### Known-good input handling

Stages prefer a `logcounts`/`vst` assay and otherwise log2-transform raw counts
internally. If you store pre-transformed data under a `counts` name it gets
double-logged — the pipeline detects and warns with `possibly_log_scaled`, and
the vignette documents the fix (rename the assay `logcounts`).

---

## 5. Quality gates

`rnaSentry` does not ship without evidence that it works.

### Automated tests
- **445 tests passing across 143 test files**, 0 failures, 0 errors, with 32
  expected `events_per_parameter` warnings (these are deliberate guardrail
  assertions). Shared fixtures (`helper-fixtures.R`) build synthetic survival
  cohorts, adversarial PCA cohorts, and in-memory signature objects so tests
  stay fast and network-free.

### Validation dossier (`validation/`)
- **Simulation study** (`simulate_study.R`): 15 falsifiable targets, 15 PASS —
  including the "small-cohort trap" (why low event counts should fail loudly).
- **Real-data reproduction**: `repro_gse20685.R` (breast cancer) and
  `repro_live_geo.R` reproducing a discovery→external transfer
  (GSE31210 → GSE50081, LUAD) where transfer was *not* demonstrated and the
  reasons (low event count, CV optimism) are quantified.
- **Cross-platform check** (`docker_check.sh`), a **face-validity review**, a
  **maintainer testing guide**, and submission-success criteria tied to
  Bioconductor.

### Build checks
- `R CMD check` on the built tarball (`rnaSentry_0.99.0.tar.gz`): **0 errors /
  1 warning (qpdf) / 1 note (tidy)**; `checkRd` clean (14/14).
- BiocCheck baseline: 1 error (404) / 1 warning / 13 notes.

---

## 6. What's in the repo

```
R/           14 source files (one per stage + shared utils.R)
tests/       testthat suite (15 test files + helper-fixtures.R)
man/         14 Rd help pages (generated from roxygen)
inst/        report_template.Rmd, CITATION, extdata/gse20685_case_study.rds
vignettes/   introduction + case studies (breast cancer, confounder audit,
             small cohort)
validation/  the dossier described above
```

Exported API (14 functions): `build_signature`, `cox_model`, `design_audit`,
`generate_report`, `km_curve`, `load_counts`, `lock_signature`, `pca_audit`,
`plot_pca_audit`, `qc_explore`, `run_rnaSentry`, `sex_check`,
`survival_parametric`, `validate_external`.

---

## 7. Scope, assumptions, and provisional status

**Assumptions (stated up front in the vignette):**
- Bulk RNA-seq, standard right-censored (time, event) endpoints. Not single-cell,
  not binary/continuous outcomes, not competing risks, not time-varying
  covariates, not left truncation.
- Signature = weighted linear combination of expression; higher score = higher
  hazard (shorter survival); all concordance uses this convention.
- Confounders are **detected and reported**; the pipeline never batch-corrects
  and makes no differential-expression calls.
- CV concordance is screening-internal (optimistic); rely on
  `validate_external()` for the fully out-of-sample estimate.

**Provisional status:** v0.99.0 is a pre-1.0 release with no Bioconductor
release yet. The naming (`rnaSentry`, "sentry" = guard/watch) matches the
design: the package guards the user against silent data or modeling failures at
every step.
