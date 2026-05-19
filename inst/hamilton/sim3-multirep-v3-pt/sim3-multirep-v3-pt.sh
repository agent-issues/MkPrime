#!/bin/bash
#SBATCH --job-name=mkp-mr3-pt
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=06:00:00
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-sim3-multirep-v3-pt/logs/%x_%j.out
#SBATCH --error=/nobackup/%u/mkp-sim3-multirep-v3-pt/logs/%x_%j.err
#SBATCH --export=ALL

# Sim 3 multirep-v3 PT — single rep 05, blind + aware, 4 PT chains each.
#
# Tests whether PT overcomes the mixing problem in single-chain aware
# (P(AC) mean 0.012 across 8 reps, 94-167 unique topologies).
#
# PT is sequential in MkPrime so cpus-per-task=1 is correct.
# mem=1G: note this is 4x the chain state of single-chain; if OOM
# in first ~10 min raise to 4G and resubmit (checkpoint resumes).
#
# Submit: sbatch sim3-multirep-v3-pt.sh

set -euo pipefail

PROJECT_DIR=/nobackup/$USER/mkp-sim3-multirep-v3-pt
LIB=/nobackup/$USER/mkp-sim3-multirep-v3/lib
REPO=/nobackup/$USER/mkp-sim3-multirep-v3/MkPrime
mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/results/rep05"

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKP_LIB=$LIB
export MKP_REPO_ROOT=$REPO

echo "=== Sim 3 multirep-v3 PT: rep 05, nIter 100000, nChains 4 ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "TMPDIR:  $TMPDIR"
echo "Library: $LIB"
echo "Repo:    $REPO"

Rscript "$REPO/inst/hamilton/sim3-multirep-v3-pt/run_pt.R"

echo "=== Finished: $(date) ==="
du -hs "$TMPDIR" 2>/dev/null || true
