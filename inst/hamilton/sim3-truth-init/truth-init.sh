#!/bin/bash
#SBATCH --job-name=sim3ti
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --gres=tmp:2G
#SBATCH --output=/nobackup/%u/mkp-sim3-multirep-v3/logs/%x_%j.out
#SBATCH --error=/nobackup/%u/mkp-sim3-multirep-v3/logs/%x_%j.err
#SBATCH --export=ALL

# Sim 3 truth-init diagnostic: 10k iter, single chain, aware model.
#
# Verifies whether the chain holds (phi=4, pi0=0.75, theta=1, TL=13.5,
# truth tree). If it stays put, topology mixing is the only bottleneck.
# If it drifts in <1k iter, params need tightening first.

set -euo pipefail

PROJECT_DIR=/nobackup/$USER/mkp-sim3-multirep-v3
LIB=$PROJECT_DIR/lib
REPO=$PROJECT_DIR/MkPrime
mkdir -p "$PROJECT_DIR/logs"

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKP_LIB=$LIB
export MKP_REPO_ROOT=$REPO

N_ITER=${N_ITER:-10000}

echo "=== truth-init: nIter $N_ITER ==="
echo "Started: $(date) on $(hostname)"

Rscript "$REPO/inst/hamilton/sim3-truth-init/truth_init.R" "$N_ITER"

echo "=== Finished: $(date) ==="
