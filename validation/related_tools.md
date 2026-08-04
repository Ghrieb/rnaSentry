# rnaSentry — related tools and comparison notes

Living notes for the paper's "related work" section. Each tool gets a
short profile; the table summarises the comparison. New tools found during
the project should be added as a row here (with date and source) so the
paper reflects the full landscape.

Reviewed 2026-08-04.

## Comparison table

| | **rnaSentry** (0.99.0, submitting) | **asuri** (Bioc 3.23, v1.0.0) | **signifinder** (Bioc 3.19, v1.6.0) |
|---|---|---|---|
| Purpose | Discover a new prognostic signature from the user's cohort | Discover survival marker genes + risk score from the user's cohort | Apply 60+ published cancer signatures as single-sample scores |
| Selection method | Univariate Cox screen + joint multivariate Cox fit; repeated stratified CV | Subsampling glmnet (Lasso) + univariate Cox | GSVA-based scoring of literature signatures |
| Survival output | KM, log-rank, Cox PH, concordance, external validation | Risk score, KM stratification, ROC curves (ROCR) | Signature scores (+ survival association plots) |
| Data scope | Bulk RNA-seq | Bulk (SummarizedExperiment) | Bulk + single-cell + spatial |
| Guardrails / auditability | Central feature (flag ledger + HTML audit trail) | Statistical robustness via subsampling | None (scoring-focused) |
| Dependencies | 5 Imports (minimal) | ~13 Imports (glmnet, siggenes, survcomp, ROCR) | Heavy (~30, TxDb/annotation/GSVA/ComplexHeatmap) |
| License | MIT | LGPL-3 | AGPL-3 |
| In Bioconductor since | — | BioC 3.23 (< 6 months) | BioC 3.16 (~2 years) |
| Source | `https://github.com/Ghrieb/rnaSentry` | `https://github.com/jdelasrivas-lab/asuri` | `https://github.com/CaluraLab/signifinder` |

## asuri

- Full name: ASURI — Analysis of SUrvival and patients RIsk prediction based
  on gene signatures.
- Method: two main steps, subsampling glmnet and univariate Cox, to discover
  survival markers related to a clinical variable and to predict a risk score.
  Provides patient-stratification plots, gene-relevance plots, and a robust
  version of Kaplan-Meier curves.
- Maintainer: Alberto Berral-Gonzalez (`aberralgonzalez@gmail.com`), Instituto
  de Biología Molecular y Celular del Cáncer (IBMCC) / jdelasrivas-lab.
- Released: BioC 3.23 (R >= 4.5.0), version 1.0.0, < 6 months in Bioconductor.
- Overlap with rnaSentry: the closest methodological analogue — de novo
  survival-signature discovery from a `SummarizedExperiment` with survival
  outcome. Key differences: (i) Lasso-based selection with subsampling vs our
  univariate screen + joint Cox with repeated stratified CV; (ii) ROC/ROCR
  evaluation vs our concordance + external-validation focus; (iii) no
  guardrail/audit layer comparable to rnaSentry's flag ledger.
- How to position: acknowledge asuri as the methodological neighbour and state
  our design choice explicitly — interpretability of a sparse, directly
  inspectable linear score with stability reported by repeated CV, plus
  enforced auditability and honest failure modes (events-per-parameter and
  sparse-event warnings) rather than silent degradation.
- Cite: package page `https://bioconductor.org/packages/asuri/`, DOI
  `10.18129/B9.bioc.asuri`.

## signifinder

- Full name: signifinder — collection and implementation of public
  transcriptional cancer signatures.
- Method: computes single-sample scores for > 60 distinct signatures collected
  from the literature, relating to multiple tumors and cancer processes.
  Supports single-cell and spatial data (SpatialExperiment). Imports GSVA,
  ensembldb, TxDb annotation packages, ComplexHeatmap, survminer, etc.
- Maintainer: Stefania Pirrotta (`stefania.pirrotta@phd.unipd.it`), CaluraLab.
- Released: BioC 3.16 (R >= 4.3.0), version 1.6.0, ~2 years in Bioconductor.
- Overlap with rnaSentry: minimal — it answers "what do known signatures say
  about this sample?" whereas rnaSentry answers "can a new survival signature
  be learned from this cohort?". Complementary rather than competing; the
  overlap is limited to working on gene expression and some survival-association
  plots.
- How to position: no direct competition; can be cited as the existing
  infrastructure for applying published signatures, with rnaSentry filling the
  discovery + auditability niche.
- Cite: package page `https://bioconductor.org/packages/signifinder/`, DOI
  `10.18129/B9.bioc.signifinder`.

## Positioning narrative for the paper

- Distinct niche claimed: **signature discovery with an enforced audit trail
  and honest failure modes**. Neither asuri nor signifinder treats
  quality-control, confounder detection, and explicit gating/flags as a
  first-class product feature.
- Methodological honesty: our linear-score design (univariate screen + joint
  Cox + repeated stratified CV) is deliberately more interpretable than Lasso
  alternatives; stability is reported, and underpowered cohorts are flagged
  (`events_per_parameter`, `sparse_events`) rather than quietly returning
  unstable numbers.
- Scope boundaries: bulk RNA-seq survival discovery only; single-cell/spatial
  scoring and DE-calling are explicitly out of scope (see "Assumptions and
  limitations" in README, vignette, and man pages).
- Reuse opportunity: signifinder could be used downstream of rnaSentry to
  benchmark the discovered signature against known published signatures.

## Template for new entries

When a new related tool is found (date, URL):

- Name, maintainer, version, Bioc tenure, source URL.
- Method in 2-3 sentences.
- Overlap with rnaSentry and how to position against it.
- Citation (package page / DOI).
