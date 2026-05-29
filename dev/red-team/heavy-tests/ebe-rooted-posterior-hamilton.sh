#!/bin/bash
# SLURM submission for the EBE rooted-topology posterior heavy test (full scale).
#
# Verifies that the ecology-mode topology moves (NNI/SPR/TBR) sample the correct
# posterior over ROOTED trees, by comparing the chain's rooted-topology marginal
# to the exact prior-predictive evidence per topology.
#
# Pre-build note (feedback_pkgload_prebuild): the script calls
# devtools::load_all(), which compiles src/ on first use.  When running multiple
# array tasks against a shared ${SRC}, PRE-BUILD ONCE on the login node before
# sbatch so tasks do not race-clobber src/.  This single-task job is safe, but
# keep the pattern if you fan out.
#
# Resource ask: 1 node, 8 cores, ~10 h, < 2 GB.  The work is CPU-light per call
# (tiny trees) but high call-count; cores help if you parallelise the inner
# (move x replicate) loop (not done by default — single-threaded is fine within
# walltime).
#
#SBATCH --job-name=ebe-rootpost
#SBATCH --output=ebe-rootpost-%j.out
#SBATCH --error=ebe-rootpost-%j.err
#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=4G
#SBATCH --partition=shared

set -euo pipefail

# SRC = the worktree root containing the MkPrime package.  Set explicitly; do
# NOT rely on the submit CWD (feedback_slurm_inscript_path).
SRC="${SRC:-$HOME/MkPrime/ecology-aware}"
cd "${SRC}"

module load R 2>/dev/null || true

echo "[$(date)] host=$(hostname) SRC=${SRC}"
echo "[$(date)] R: $(which Rscript)"

# Invoke the IN-TREE script via ${SRC}/ (never an out-of-tree scratch copy).
Rscript "${SRC}/dev/red-team/heavy-tests/ebe-rooted-posterior.R" --full

echo "[$(date)] done. Verdict:"
cat "${SRC}/dev/red-team/heavy-tests/ebe-rooted-posterior-results/verdict.txt" || true
