#!/bin/bash
# setup_project.sh — one-time setup of the sim3-multirep-v3 project on Hamilton
#
# Run once before the first dispatch:
#   ssh hamilton8.dur.ac.uk "bash -s" < setup_project.sh
# Or upload + execute via ssh::ssh_exec_wait().
#
# Idempotent: safe to re-run after a git pull or after pulling new
# package changes from a different branch.

set -euo pipefail

PROJECT_DIR=/nobackup/$USER/mkp-sim3-multirep-v3
LIB=$PROJECT_DIR/lib
REPO=$PROJECT_DIR/MkPrime
BRANCH=${MKP_BRANCH:-worktree-ecology-aware}

mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/results" "$LIB"

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=$LIB

cd "$PROJECT_DIR"

if [ ! -d "$REPO/.git" ]; then
  echo "Cloning MkPrime ($BRANCH)..."
  git clone --depth 1 --branch "$BRANCH" \
    https://github.com/Mk-prime/r.git "$REPO"
else
  echo "Updating MkPrime ($BRANCH)..."
  cd "$REPO"
  git fetch origin "$BRANCH"
  git checkout "$BRANCH"
  git reset --hard "origin/$BRANCH"
  cd "$PROJECT_DIR"
fi

# Install/update package dependencies into project lib (Ncpus=4 for build speed).
echo "Installing dependencies into $LIB ..."
Rscript -e "\
  .libPaths(c('$LIB', .libPaths())); \
  needed <- c('Rcpp', 'TreeTools', 'TreeDist', 'TreeSearch', 'ape', \
              'phangorn', 'cli'); \
  to_install <- needed[!vapply(needed, requireNamespace, logical(1), \
                                quietly = TRUE)]; \
  if (length(to_install)) install.packages(to_install, \
                                           lib = '$LIB', \
                                           Ncpus = 4L, \
                                           repos = 'https://cloud.r-project.org')"

# Build MkPrime from source
echo "Building MkPrime from source..."
cd "$REPO"
rm -f src/*.o src/*.so src/*.dll
R CMD build --no-build-vignettes --no-manual --no-resave-data .
R CMD INSTALL --library="$LIB" MkPrime_*.tar.gz
echo "Install exit: $?"
rm -f MkPrime_*.tar.gz

echo "=== Setup complete ==="
ls "$LIB" | head
