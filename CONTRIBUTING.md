# Contributing to rnaSentry

Thanks for considering a contribution. rnaSentry is a Bioconductor package
with an unusual requirement: **every statistical decision is audited, and
every guardrail is falsifiable**. This file explains how to keep that property
intact while you change code, fix bugs, or add features.

## Scope

rnaSentry is deliberately narrow. It does bulk RNA-seq
discovery-to-external-validation of linear prognostic risk scores with
standard right-censored survival. Out of scope (and likely rejected):

- single-cell, spatial, or pseudobulk modeling;
- non-survival outcomes (continuous, binary, competing risks);
- batch-effect *correction* or differential-expression calling (detection
  and reporting only);
- nonlinear/machine-learning score forms.

If your contribution adds a new statistical stage, it must come with an
audit record and a falsifiable gate (see below).

## Reporting issues

- Open an issue at <https://github.com/Ghrieb/rnaSentry/issues>.
- Include a minimal reproducible example: the input `SummarizedExperiment`
  (or a way to build it), the exact function call with all arguments, the R
  and package version (`sessionInfo()`), and the output you expected vs. what
  you got.
- Note any flag ledger entries printed by the stage — the `(check, severity,
  detail)` records are usually the fastest way to localize a problem.

## Development workflow

1. Fork and branch; keep the default branch package-only (Bioconductor rule).
2. Follow the existing conventions: `R/*.R` one stage family per file,
   roxygen comments with `@return` documenting the result list, guardrail
   flags registered as `(check, severity, detail, stage)` in the flag ledger.
3. No new hard dependencies: the package intentionally imports only
   `SummarizedExperiment, S4Vectors, survival, withr, methods, stats, graphics`.
   If you must add
   an Imports/Suggests entry, justify it in the PR.
4. Write tests before or with the code (see below).
5. Do **not** add comments to shipped code unless they document a statistical
   decision; the audit trail is carried by flag details and the validation
   docs instead.
6. Vignettes and the report template must build **offline** (Bioconductor
   vignette builds disallow network access). Never fetch data at build time;
   bundle small `inst/extdata` subsets instead.

## The gate suite

Every change must pass, in order, the same gates the maintainer runs. Full
instructions and expected values are in `validation/maintainer_testing_guide.md`.

| Gate | What it runs | How |
|---|---|---|
| 0 | Hygiene + check | `R CMD build`, `R CMD check --as-cran` on the tarball (expect 0 ERROR / 0 WARNING / 1 NOTE - the residual NOTE is the HTML `tidy` external tool; the `qpdf` WARNING cleared 2026-08-05), `BiocCheck` |
| 1 | Full test suite | `Rscript <temp>/run_tests_parity.R` (devtools::test) — 147 blocks / 457 expectations / 0 fail / 0 error / 32 expected guardrail warnings |
| 2 | Adversarial review + stale-number sweep | `powershell -ExecutionPolicy Bypass -File validation/grep_stale_numbers.ps1` (0 hits outside the allow-list: `validation/*`, `vignettes/case-study-brca.Rmd`, `vignettes/case-study-impact.Rmd`) |
| 3 | Statistical parity | every statistic re-implemented independently and asserted equal (log-rank, Cramér's V, Cox HR/CI/p, Schoenfeld, AIC/loglik, fold CV C, `C + C_rev = 1`) |
| 4 | Simulation study | `Rscript validation/simulate_study.R` — **15 falsifiable targets / 15 PASS** (Sims 1-5 + Sim 6 power analysis) |
| 5 | Real-data face validity | `Rscript validation/repro_gse20685.R` (bundled subset) and the live GEO pair (Gate 5b) |
| 5c | Case-study vignettes | `vignettes/case-study-*.Rmd` must build offline with numbers matching the validation dossier's recorded values |

### Rules that make the gates meaningful

- **Never weaken a guardrail to make a test pass.** If a flag should fire,
  argue why it should not — with data — rather than silencing it.
- **Never change the sign convention.** Concordance uses
  `reverse = TRUE` (higher risk score = shorter survival); a parity test
  asserts `C + C_rev = 1`. This convention is load-bearing for every
  downstream number.
- **Keep the simulation falsifiable.** Each new Sim must state its target as
  a pass/fail predicate before it is written, and the script must exit
  non-zero on failure (it is CI-gated).
- **Reproducibility is enforced, not requested.** New code paths that select
  genes or fit models must be seed-deterministic; the lock
  (`lock_signature()`) must not be bypassable silently.

## Full regression loop after any code change

1. Gate 1 suite (147/457, 0 fail, 0 error).
2. `R CMD build` + `R CMD check --as-cran` on the tarball (0 ERROR / 0
   WARNING / 1 NOTE, `tidy` only).
3. If statistics or tests changed: Gate 3 parity suite + Gate 4 simulation
   (15/15).
4. If intake/annotation handling changed: Gate 5 repro on the bundled subset.
5. Gate 2 stale-number sweep.
6. Commit logs together with the change at a logical checkpoint (see
   `maintainer_testing_guide.md` for the log naming convention).

## Site deployment

The pkgdown site (https://ghrieb.github.io/rnaSentry/) is rebuilt and
deployed **automatically** by `.github/workflows/pkgdown.yaml` on every push
to `main`: it installs the package with its Bioconductor dependencies, runs
`pkgdown::build_site_github_pages()`, and pushes the result to the `gh-pages`
branch. You do **not** need to touch `docs/` or `gh-pages` manually — a normal
commit to `main` syncs the site. `docs/` stays gitignored (local preview
only); the workflow is the single source of truth for what is published.
