# rnaSentry — cross-platform check status

Target: Bioconductor 3.21 devel (R 4.5), Linux + macOS + Windows. The gate
machine is Windows 11 / R 4.5.2 only, so the devel/other-platform runs need a
container or a remote service. This document records what is done, what is
blocked, and exactly how to finish it.

## Done locally (Windows, R 4.5.2)

- `R CMD build` — tarball `rnaSentry_0.99.0.tar.gz` builds cleanly (LF/UTF-8,
  no stray dirs), and the compiled vignette `inst/doc/rnaSentry.html` is in the
  tarball (verified: `vignette(package="rnaSentry")` lists it, `browseVignettes`
  resolves it).
- `R CMD check --as-cran` on the fresh tarball → **0 ERROR, 0 WARNING,
  1 NOTE (HTML `tidy` missing, external tool)**. The former `qpdf` WARNING
  is gone — `qpdf` was installed 2026-08-05. Logs: `logs/16_as_cran.txt`
  (Ultimate Pre-Flight), re-confirmed at the BiocParallel round
  `logs/17_as_cran.txt`.
- Full `devtools::test()` → **147 blocks / 457 passed / 0 failed / 0 error**
  (32 expected warnings from the `events_per_parameter` guardrail on
  deliberately small fixtures; includes the `load_counts()` intake and
  enforceable-lock tests added 2026-08-04, the withr-scoped seeding refactor
  of 2026-08-05, and the 4 serial↔parallel `BPPARAM` blocks added in the
  same round). Logs: `logs/00.1_test.txt`, re-run `logs/17_test.txt`.

## Partial Linux evidence via CI (2026-08-05)

The pkgdown workflow (`.github/workflows/pkgdown.yaml`) runs on
`ubuntu-latest` with **R 4.6.1 / Bioc release** (SummarizedExperiment 1.42,
Biobase 2.72, BiocStyle 2.40). Run #1 was green end-to-end: the package
installed with all Bioconductor dependencies, and **all five vignettes
compiled offline** with the same pinned numbers as the local R 4.5.2 / Bioc
3.22 runs (0.316/0.684/1.0, CV C = 0.797, `C + C_rev = 1`). The
getting-started round added a sixth vignette that run #3 compiled the same
way (commit `11d575a` -> `gh-pages` `f9646cd`, 2026-08-05): **all six
vignettes now compile on Linux/Bioc release**. This is genuine
non-Windows evidence that the package installs and the vignettes build on a
Linux/Bioc release environment — but it is **not** a full `R CMD check` or a
test-suite run, which still require the container/service options below.

## Blocked (needs a Docker-capable machine or remote service)

The local machine has **no Docker** and no `R CMD BiocCheck` CLI launcher
(Windows), and R 4.5.2 ≠ the Bioc 3.21 devel R-devel build.

### Status note (2026-08-05)

Phase 4 of the submission checklist (cross-platform / devel-container check)
was **skipped by decision on 2026-08-05** — rationale and follow-up are
recorded in `maintainer_testing_guide.md` (Gate 6) and
`submission_success_criteria.md`. The options below
remain the documented finishing path.

### Option A — Bioconductor devel container (recommended)

`validation/docker_check.sh` mounts this repo into
`bioconductor/bioconductor_docker:devel` and runs, in that exact devel
environment:

1. `R CMD build` (source tarball)
2. `R CMD BiocCheck <tarball>` (Linux CLI, not the Windows-broken launcher)
3. `R CMD check --as-cran <tarball>`
4. the full test suite

All output goes to `validation/logs/docker/`. Run:

```bash
bash validation/docker_check.sh
```

on any Docker-capable machine. This resolves both blockers at once and is the
closest reproduction of the Bioc build machines.

### Option B — rhub (Windows + macOS + Linux flavors)

```r
rhub::rhub_check("rnaSentry_0.99.0.tar.gz")
```

Requires a GitHub account/token (used to read package sources) — see
`rhub::rs_connect()`.

### Option C — win-builder (Windows only)

Upload `rnaSentry_0.99.0.tar.gz` at
https://win-builder.r-project.org/ (devel + release). Requires an email address
to which the check logs are sent.

## Procedure when the results arrive

1. Copy the container/service log excerpts into this document (replace this
   status text with the actual `Status:` lines from each platform).
2. Confirm BiocCheck is clean except the documented items in
   `logs/notes_documented.md` (support-site email 404 → register
   `ghriebabdelkarimhani@gmail.com` on https://support.bioconductor.org; the
   `CITATION` DOI placeholder).
3. Re-run `bash validation/docker_check.sh` if anything fails and iterate.
4. Add a `PASS/FAIL` line per platform to the submission notes.
