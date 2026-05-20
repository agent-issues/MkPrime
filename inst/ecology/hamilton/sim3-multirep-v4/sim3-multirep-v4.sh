#!/bin/bash
#SBATCH --job-name=sim3v4
#SBATCH --partition=shared
#SBATCH --array=1-8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=05:00:00
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-sim3-multirep-v4/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-sim3-multirep-v4/logs/%x_%A_%a.err
#SBATCH --export=ALL

# Sim 3 v4 multirep — 8-chain parallel tempering, array job.
#
# Reuses mkp-sim3-multirep-v3/{lib,MkPrime} (no reinstall needed).
# Outputs to mkp-sim3-multirep-v4/results/rep%02d/.
# nChains = 8 (PT chains sequential, ~2-3h per rep; 5h wall is safe).
#
# Pre-flight (login node):
#   mkdir -p /nobackup/$USER/mkp-sim3-multirep-v4/{logs,results}
#   cd /nobackup/$USER/mkp-sim3-multirep-v3/MkPrime
#   git pull
#   R CMD build --no-build-vignettes --no-manual --no-resave-data .
#   R CMD INSTALL --library=/nobackup/$USER/mkp-sim3-multirep-v3/lib MkPrime_*.tar.gz
#   rm -f MkPrime_*.tar.gz
#
# Submit:        sbatch sim3-multirep-v4.sh
# Resume single: sbatch --array=5 sim3-multirep-v4.sh

set -euo pipefail

V3_DIR=/nobackup/$USER/mkp-sim3-multirep-v3
V4_DIR=/nobackup/$USER/mkp-sim3-multirep-v4

LIB=$V3_DIR/lib
REPO=$V3_DIR/MkPrime

mkdir -p "$V4_DIR/logs" "$V4_DIR/results"

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKP_LIB=$LIB
export MKP_REPO_ROOT=$REPO

REP_ID=$SLURM_ARRAY_TASK_ID
N_ITER=${N_ITER:-100000}
SEED_BASE=${SEED_BASE:-20260601}

echo "=== Sim 3 v4 multirep: rep $REP_ID, nIter $N_ITER, nChains=8, seed base $SEED_BASE ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "TMPDIR:  $TMPDIR"
echo "Library: $LIB"

Rscript "$REPO/inst/ecology/hamilton/sim3-multirep-v4/run_rep.R" \
        "$REP_ID" "$N_ITER" "$SEED_BASE"

echo "=== Finished: $(date) ==="
du -hs "$TMPDIR" 2>/dev/null || true
