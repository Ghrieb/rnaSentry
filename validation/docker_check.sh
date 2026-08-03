#!/usr/bin/env bash
#
# rnaSentry — Bioconductor devel container check.
#
# The machine used for the local gate has no Docker and no `R CMD BiocCheck`
# launcher (Windows). Run this on a Docker-capable machine to reproduce the
# Bioconductor devel build exactly:
#
#   1. Make sure the rnaSentry source is inside the container mount.
#      By default this script mounts this repository (the directory one
#      level above the package) read-only.
#   2. bash validation/docker_check.sh
#
# It runs, inside bioconductor/bioconductor_docker:devel (tracks R-devel +
# the current Bioc devel branch):
#     * R CMD build  -- builds the source tarball (excludes validation/)
#     * R CMD BiocCheck <tarball>          (CLI, Linux container)
#     * R CMD check  --as-cran <tarball>
#     * the full test suite
# All output is written to validation/logs/docker/.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PKG_DIR="$REPO_ROOT/rnaSentry"
IMAGE="${IMAGE:-bioconductor/bioconductor_docker:devel}"
LOGDIR="$PKG_DIR/validation/logs/docker"
mkdir -p "$LOGDIR"

echo "== rnaSentry Docker gate =="
echo "repo root : $REPO_ROOT"
echo "image     : $IMAGE"
echo "logs      : $LOGDIR"
echo

docker pull "$IMAGE"

docker run --rm -v "$REPO_ROOT:/mnt/pkg:ro" -w /mnt/pkg/rnaSentry "$IMAGE" bash -c '
  set -e
  export RSTUDIO_PANDOC="$(command -v pandoc || true)"
  echo "R: $(R --version | head -n1)"
  echo "PANDOC: ${RSTUDIO_PANDOC:-none}"

  echo
  echo "=== 1/4 R CMD build ==="
  cd /mnt/pkg
  R CMD build --no-manual rnaSentry

  echo
  echo "=== 2/4 R CMD BiocCheck (tarball) ==="
  R CMD BiocCheck rnaSentry_0.99.0.tar.gz || true

  echo
  echo "=== 3/4 R CMD check --as-cran (tarball) ==="
  R CMD check --no-manual --as-cran rnaSentry_0.99.0.tar.gz || true

  echo
  echo "=== 4/4 devtools::test() ==="
  Rscript -e "suppressPackageStartupMessages(library(devtools)); devtools::load_all(\".\", quiet=TRUE); print(devtools::test(reporter=\"summary\"))"
' 2>&1 | tee "$LOGDIR/container_run.log"

echo
echo "Done. Full container log: $LOGDIR/container_run.log"
echo "Copy these into validation/cross_platform.md and the submission notes."
