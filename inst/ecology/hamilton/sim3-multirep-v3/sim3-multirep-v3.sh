#!/bin/bash
#SBATCH --job-name=sim3v3
#SBATCH --partition=shared
#SBATCH --array=1-8
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-sim3-multirep-v3/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-sim3-multirep-v3/logs/%x_%A_%a.err
#SBATCH --export=ALL

# Sim 3 v3 multirep — array job, one rep per task
#
# Resource profile: single-core R, 8 GB mem, 2h wall (pilot ~50 min).
# Checkpoints written to per-rep dirs; rerun the same task to resume
# from where it left off.
#
# Submit:        sbatch sim3-multirep-v3.sh
# Submit larger: sbatch --array=1-20 sim3-multirep-v3.sh
# Resume single: sbatch --array=5 sim3-multirep-v3.sh

set -euo pipefail

PROJECT_DIR=/nobackup/$USER/mkp-sim3-multirep-v3
LIB=$PROJECT_DIR/lib
REPO=$PROJECT_DIR/MkPrime
mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/results"

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKP_LIB=$LIB
export MKP_REPO_ROOT=$REPO

REP_ID=$SLURM_ARRAY_TASK_ID
N_ITER=${N_ITER:-100000}
SEED_BASE=${SEED_BASE:-20260601}

echo "=== Sim 3 v3 multirep: rep $REP_ID, nIter $N_ITER, seed base $SEED_BASE ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "TMPDIR:  $TMPDIR"
echo "Library: $LIB"

Rscript "$REPO/inst/ecology/hamilton/sim3-multirep-v3/run_rep.R" \
        "$REP_ID" "$N_ITER" "$SEED_BASE"

echo "=== Finished: $(date) ==="
du -hs "$TMPDIR" 2>/dev/null || true
