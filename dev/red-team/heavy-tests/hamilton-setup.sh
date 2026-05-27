#!/bin/bash
# Run on Hamilton login node: set up isolated red-team scratch under
# /nobackup/$USER/mkp-study/red-team/, with main-HEAD worktree and isolated lib.
set -euo pipefail

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
LIB=${RT}/lib
DEPS_LIB=${PROJECT}/lib

mkdir -p "${RT}/heavy-tests" "${RT}/results" "${RT}/logs" "${LIB}"

cd "${PROJECT}/mkp"
git fetch origin main
if [ ! -d "${SRC}" ]; then
  git worktree add "${SRC}" origin/main
else
  cd "${SRC}"
  git fetch origin main
  git checkout origin/main
fi

cd "${SRC}"
echo "--- main HEAD ---"
git log --oneline -1
echo "--- packages in deps lib ---"
ls "${DEPS_LIB}" | wc -l

module load r/4.5.1
module load gcc/14.2 || true   # advisor: ensure C++17

# Build and install MkPrime from this worktree into red-team/lib,
# with deps lib chained via R_LIBS so we don't rebuild ape/Rcpp/TreeTools.
export R_LIBS_USER="${LIB}:${DEPS_LIB}"
export R_LIBS="${LIB}:${DEPS_LIB}"

rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="${LIB}" MkPrime_*.tar.gz
rm -f MkPrime_*.tar.gz
echo "--- installed MkPrime ---"
ls "${LIB}/MkPrime/DESCRIPTION" && head -3 "${LIB}/MkPrime/DESCRIPTION"
