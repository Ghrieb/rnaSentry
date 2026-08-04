# rnaSentry — BiocCheck NOTES / WARNINGS documentation

Date: 2026-08-03 | BiocCheck 1.44.2 (Bioc 3.21) | R 4.5.2 | Windows 11

This document records every BiocCheck ERROR / WARNING / NOTE that the gate
does not (or cannot) resolve, with a justification for each. This is the
`cran-comments.md`-equivalent for the Bioconductor submission.

## BiocCheck::BiocCheck(".")

### 1 ERROR — support site registration (environmental, action required)
`Unable to find your email in the Support Site: HTTP 404 Not Found.`
BiocCheck queries `https://support.bioconductor.org` for the maintainer email
(`ghriebabdelkarimhani@gmail.com`). The 404 means that email address is not
registered on the Bioconductor Support Site (this run also intermittently
timed out against that site).
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
   resolves too. Record the clean run as `logs/08_bioccheck.txt` and update
   the "Verification status" table in `test_audit.md`.
This check is not a code defect.

### 1 WARNING — `set.seed` usage (justified)
`Remove set.seed usage (found 1 times): R/build_signature.R (line 233)`
The call is `if (!is.null(seed)) set.seed(seed)` inside `build_signature()`:
a public, documented `seed` argument (default `NULL`) that makes the
repeated cross-validation reproducible. This is the standard, recommended
pattern for seed-parameterized functions and is exercised by the test suite
(identical results under identical seed). We intentionally keep it.

### 12 NOTES — all advisory
- **Update R version dependency 4.4.0 -> 4.5.0**: rnaSentry requires R >= 4.4.0.
  Bioc 3.21 builds on R 4.5; the lower bound is a superset and is left
  intentionally permissive. No action.
- **Consider adding automatically suggested biocViews: KEGG**: KEGG pathway
  terms are not used by the package; current `biocViews` cover the actual
  functionality. No action.
- **Consider adding maintainer's ORCID iD**: no ORCID was supplied by the
  maintainer. To be added when available.
- **No 'fnd' role found in Authors@R**: the work is not grant-funded; the
  `fnd` role does not apply.
- **Avoid 'suppressWarnings'/'*Messages' if possible (7)**: two locations
  flagged — `R/build_signature.R:277` wraps a per-gene `coxph` fit during
  univariate screening (thousands of fits; convergence warnings are expected
  and handled), and `R/validate_external.R:224` wraps `survival::concordance`.
  Both are deliberate, narrowly-scoped suppression around a single robust
  statistic. Justified.
- **Function length > 50 lines (10 functions)**: long functions implement
  multi-step guarded pipelines (e.g., `build_signature`) with explicit
  per-step audit flags; splitting would scatter the audit logic. Accepted
  stylistic note.
- **Consider adding runnable examples to exported man pages**: all 13 man
  pages contain `\examples`; `generate_report` and `run_rnaSentry` wrap
  theirs in `\dontrun{}` because a full run requires a prepared
  `SummarizedExperiment` and, for `generate_report`, report rendering. The
  same end-to-end paths are covered by `tests/testthat/` (incl. report
  rendering) and by `validation/`. Converting to `\donttest` without
  embedding fixture construction would break `R CMD check`.
- **dontrun/donttest usage (15% of man pages)**: see previous item; the two
  examples use `\dontrun` deliberately.
- **Consider shorter lines (2% > 80 chars)** and **multiples-of-4 indents
  (40%)**: cosmetic; indentation is 2-space by house style. Accepted.
- **Cannot determine whether maintainer is subscribed to Bioc-Devel**:
  network-only check against support site; see ERROR above.

## BiocCheck::BiocCheckGitClone()

### 1 WARNING — CITATION `doi` argument missing or empty (justified)
The package is pre-publication; no DOI exists yet. The CITATION entry
carries the maintainer's real identity and repository URL; the `doi` field
will be filled when a preprint/paper DOI is assigned. BiocCheck's companion
hint ("only include a CITATION file if there is a preprint or publication")
is noted; we keep a correct, honest CITATION.
**Action for maintainer:** update `inst/CITATION` `doi` once a DOI exists.

### Procedure note
`BiocCheck(".")` writes a `rnaSentry.BiocCheck/` stamp directory into the
package directory; `BiocCheckGitClone()` then reports it as an untracked
"system file". The gate therefore deletes the stamp between the two checks
and runs the GitClone check on the clean tree. The stamp is git-ignored.

## R CMD check --as-cran (tarball) — remaining items

### 1 WARNING — `qpdf` is needed for checks on size reduction of PDFs
`qpdf` is not installed on this Windows machine. It is an optional external
tool used only for PDF size-reduction checks; Bioconductor's build machines
have it. No action required from the package.

### 2 NOTEs — environmental
- **Top-level files: `README.md`/`NEWS` cannot be checked without
  `pandoc`**: the check's own pandoc discovery did not honor the
  `RSTUDIO_PANDOC` environment variable used elsewhere in the gate.
  Bioconductor build machines have pandoc installed; the vignette and this
  gate build it successfully with `RSTUDIO_PANDOC` set.
- **Examples with CPU/elapsed > 5s (build_signature, 7.4s user)**:
  investigated — the entire cost is `library(SummarizedExperiment)` in the
  fresh check session (measured 7.60s loading GenomicRanges, IRanges,
  GenomeInfoDb, Biobase, S4Vectors, matrixStats, MatrixGenerics). The
  example's own computation is 0.17s (build_signature call) + 0.10s (data
  construction). This is a known, acceptable NOTE for Bioconductor packages
  whose examples construct `SummarizedExperiment` objects; the threshold is
  advisory and Bioc build machines are faster than this Windows session.
  Example logic is kept minimal (`20x20` matrix, `repeats = 1, folds = 2`).
