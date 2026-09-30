#!/bin/bash
# Build MkPrime_<version>.tar.gz from the committed HEAD with `RemoteSha:` in
# DESCRIPTION, so the installed package can say which commit it is
# (packageDescription("MkPrime")$RemoteSha). A plain tarball install leaves it NULL.
#
# Usage, from a git checkout:  bash data-raw/hamilton/build_mkp_tarball.sh [outdir]
# Builds HEAD, not the working tree: commit first.
set -euo pipefail

OUTDIR=$(cd "${1:-.}" && pwd)
SHA=$(git rev-parse HEAD)
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

git archive HEAD | tar -x -C "$TMP"
sed -i '/^RemoteSha:/d' "$TMP/DESCRIPTION"
[ -n "$(tail -c1 "$TMP/DESCRIPTION")" ] && echo >> "$TMP/DESCRIPTION"
echo "RemoteSha: $SHA" >> "$TMP/DESCRIPTION"

(cd "$OUTDIR" && R CMD build --no-build-vignettes --no-manual --no-resave-data "$TMP")
echo "Built from $SHA in $OUTDIR"
