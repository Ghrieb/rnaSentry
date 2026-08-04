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
- `R CMD check --as-cran` on the fresh tarball → **0 ERROR, 1 WARNING
  (`qpdf` missing, external tool), 1 NOTE (HTML `tidy` missing, external
  tool)**; both external-tool only. Log: `logs/00.4_as_cran.txt`.
- Full `devtools::test()` → **143 blocks / 445 passed / 0 failed / 0 error**
  (32 expected warnings from the `events_per_parameter` guardrail on
  deliberately small fixtures; includes the `load_counts()` intake and
  enforceable-lock tests added 2026-08-04). Log: `logs/00.1_test.txt`.

## Blocked (needs a Docker-capable machine or remote service)

The local machine has **no Docker** and no `R CMD BiocCheck` CLI launcher
(Windows), and R 4.5.2 ≠ the Bioc 3.21 devel R-devel build.

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
   `CITATION` DOI placeholder; the justified `set.seed` warning).
3. Re-run `bash validation/docker_check.sh` if anything fails and iterate.
4. Add a `PASS/FAIL` line per platform to the submission notes.
