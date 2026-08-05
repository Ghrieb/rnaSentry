# rnaSentry — BiocCheck NOTES / WARNINGS documentation

Date: 2026-08-05 | BiocCheck 1.44.2 (Bioc 3.21) | R 4.5.2 | Windows 11
(gate logs: `16_bioccheck_tarball.txt`, `16_bioccheck.txt`,
`16_bioccheck_gitclone.txt`, `16_as_cran.txt` — the Phase-3 "Ultimate
Pre-Flight" round, run after the `withr::with_seed()` seeding refactor,
the `seq_len()`/`\donttest` fixes, and the local install of the qpdf CLI)

This document records every BiocCheck ERROR / WARNING / NOTE that the gate
does not (or cannot) resolve, with a justification for each. This is the
`cran-comments.md`-equivalent for the Bioconductor submission.

## BiocCheck::BiocCheck(tarball) and BiocCheck::BiocCheck(".") — 1 ERROR / 0 WARNING / 8 NOTES

### 1 ERROR — support site registration (environmental, action required)
`Unable to find your email in the Support Site: HTTP 404 Not Found.`
BiocCheck queries `https://support.bioconductor.org` for the maintainer email
(`ghriebabdelkarimhani@gmail.com`). The 404 means that email address is not
registered on the Bioconductor Support Site.
**Action for maintainer:** register `ghriebabdelkarimhani@gmail.com` on the
Bioconductor Support Site before submission. Steps:
1. Open https://support.bioconductor.org and click **Sign Up** (top right).
2. Complete the form with the **same email address** used in `DESCRIPTION`
   `Authors@R` `cre` (`ghriebabdelkarimhani@gmail.com`); choose a username and
   password.
3. Confirm the activation email from the support site.
4. (Optional but recommended) complete the profile — name, ORCID, affiliation.
5. Re-run `BiocCheck::BiocCheck(".")`; this ERROR should clear, and the
   "maintainer subscribed to Bioc-Devel" note (which hits the same site)
   resolves too. Record the clean run and update this document.
This check is not a code defect.

### 0 WARNING — `set.seed` warning RESOLVED
The pre-Phase-3 gate carried `1 WARNING — Remove set.seed usage ...
R/build_signature.R (line 279)`. On 2026-08-05 `build_signature()` was
refactored to scope its seeding with `withr::with_seed()` (called once around
a `make_all_folds()` that draws every repeat's partitions in the same order
as before), so the caller's global RNG state is left unchanged and the
warning no longer fires. Gate 1 re-confirmed 143 blocks / 445 passed / 0
failed / 0 errors (32 expected guardrail warnings on deliberately small
fixtures); the seeded results are bit-identical to the pre-refactor runs.

### 8 NOTES — all advisory
- **Consider adding automatically suggested biocViews: SingleCell, KEGG**:
  single-cell/spatial modeling and KEGG pathway terms are not used by the
  package (both are explicitly out of scope); current `biocViews` cover the
  actual functionality. No action.
- **Consider adding maintainer's ORCID iD**: no ORCID was supplied by the
  maintainer (deliberately skipped 2026-08-05). To be added when available.
- **No 'fnd' role found in Authors@R**: the work is not grant-funded; the
  `fnd` role does not apply.
- **Avoid 'suppressWarnings'/'*Messages' if possible (7)**: locations
  flagged are `R/build_signature.R:330` (wraps a per-gene `coxph` fit during
  univariate screening; thousands of fits, convergence warnings expected and
  handled) and `R/validate_external.R:248` (wraps `survival::concordance`).
  Both are deliberate, narrowly-scoped suppression around a single robust
  statistic. Justified.
- **Function length > 50 lines (11 functions)**: long functions implement
  multi-step guarded pipelines (e.g., `build_signature`) with explicit
  per-step audit flags; splitting would scatter the audit logic. Accepted
  stylistic note.
- **Consider shorter lines (106 lines, 3% > 80 chars)** and
  **multiples-of-4 indents (1493 lines, 36%)**: cosmetic; indentation is
  2-space by house style. Accepted.
- **Cannot determine whether maintainer is subscribed to Bioc-Devel**:
  network-only check against the support site; see ERROR above.

Resolved relative to the pre-Phase-3 gate (`15_bioccheck.txt`, 13 NOTES):
`set.seed` WARNING (→ `withr::with_seed()`, above); R-version dependency
NOTE 4.4.0→4.5.0 (DESCRIPTION now `R (>= 4.5.0)` for Bioc 3.21);
`Avoid 1:` NOTE (`R/load_counts.R:112` → `seq_len()`); two `\dontrun`/`\donttest`
NOTEs and the "add runnable examples" NOTE (`generate_report` and
`run_rnaSentry` examples are now self-contained `\donttest` blocks that
construct a small `SummarizedExperiment` and run under
`R CMD check --as-cran --run-donttest`).

## BiocCheck::BiocCheckGitClone() — 0 ERROR / 1 WARNING / 0 NOTES

### 1 WARNING — CITATION `doi` argument missing or empty (justified)
The package is pre-publication; no DOI exists yet. The CITATION entry
carries the maintainer's real identity and repository URL; the `doi` field
will be filled when a preprint/paper DOI is assigned. The maintainer decided
on 2026-08-05 to keep `inst/CITATION` with `doi = ""` and this justified
warning. BiocCheck's companion hint ("only include a CITATION file if there
is a preprint or publication") is noted; we keep a correct, honest CITATION.
**Action for maintainer:** update `inst/CITATION` `doi` once a DOI exists.

### Git-system-file ERROR RESOLVED
A pre-Phase-3 GitClone run reported `System files found that should not be
Git tracked` (the stray `chk16/rnaSentry.Rcheck/` output tree inside the
package directory). The stray directory was deleted from the working tree on
2026-08-05; `.Rbuildignore` now also excludes `^chk[0-9]+$` so check-output
directories cannot leak into `R CMD build`.

### Procedure note
`BiocCheck(".")` writes a `rnaSentry.BiocCheck/` stamp directory into the
current working directory; `BiocCheckGitClone()` then reports it as an
untracked "system file". The gate therefore deletes the stamp between the
two checks and runs the GitClone check on the clean tree. The stamp is
git-ignored and excluded from the tarball via `.Rbuildignore`
(`^rnaSentry\.BiocCheck$`).

## R CMD check --as-cran (tarball, log `16_as_cran.txt`) — 0 ERROR / 0 WARNING / 1 NOTE

The qpdf WARNING from prior gates (`'qpdf' is needed for checks on size
reduction of PDFs`) is RESOLVED: on 2026-08-05 the qpdf CLI 12.3.2
(MSVC64) was installed and checks are run with
`R_QPDF=<...>\qpdf.exe` and `_R_CHECK_DOC_SIZES_=true`, which makes
`R CMD check` run its real PDF size-reduction check (`check_doc_size`)
instead of emitting the unconditional as-cran warning.

### 1 NOTE — environmental, deliberately retained
- **HTML version of manual: no command `tidy` found**: the HTML validator
  `tidy` is not installed on this Windows machine. The maintainer decided on
  2026-08-05 to skip installing it (the package ships no HTML; R CMD check
  only invokes tidy on its generated manual HTML). It is an optional
  external tool; Bioconductor's build machines have it. No package action.

### Prior env-dependent items now stable
- **Examples with CPU/elapsed > 5s**: the cost is `library(SummarizedExperiment)`
  in a fresh check session (loading GenomicRanges/IRanges/GenomeInfoDb/
  Biobase/S4Vectors/matrixStats/MatrixGenerics). Known, acceptable for
  Bioconductor packages whose examples construct `SummarizedExperiment`
  objects; the threshold is advisory and Bioc build machines are faster.
- The examples now run cleanly under `--run-donttest` (self-contained
  fixtures; report output is directed to `tempdir()` so no stray HTML is
  left in the check directory).
