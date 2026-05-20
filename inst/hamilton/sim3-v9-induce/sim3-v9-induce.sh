#!/bin/bash
#SBATCH --job-name=mkp-v9
#SBATCH --partition=shared
#SBATCH --array=1-3
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=04:00:00
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-sim3-v9-induce/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-sim3-v9-induce/logs/%x_%A_%a.err
#SBATCH --export=ALL

# Sim 3 v9-induce: 3-rep pilot array — does the v9c parsimony trap
# also fool BLIND ML MCMC, with AWARE ML rescuing truth?
#
# Each task: one rep, blind + aware, 4 PT chains each, 100k iter.
# PT is sequential -> cpus-per-task=1 is correct.
# mem 2G is conservative for 4 PT chains; bump to 4G if OOM and resubmit
# (resumes from checkpoint).
#
# Auto-expSteps: run_v9.R does NOT pass expSteps to MkPrimeModel(); the
# bc24d52 default resolves expSteps from 1.05 * parsimony at startup.
# Watch the .out log for "Tree length prior: parsimony score = N;
# expSteps = ..." to confirm.
#
# Pilot:  sbatch sim3-v9-induce.sh
# Scale:  sbatch --array=1-10 sim3-v9-induce.sh
# Resume: sbatch --array=2 sim3-v9-induce.sh  (single rep)

set -euo pipefail

PROJECT_DIR=/nobackup/$USER/mkp-sim3-v9-induce
LIB=/nobackup/$USER/mkp-sim3-multirep-v3/lib       # reuse v3 library
REPO=/nobackup/$USER/mkp-sim3-multirep-v3/MkPrime  # reuse v3 repo checkout
mkdir -p "$PROJECT_DIR/logs" "$PROJECT_DIR/results"

module load r/4.5.1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKP_LIB=$LIB
export MKP_REPO_ROOT=$REPO

REP_ID=$SLURM_ARRAY_TASK_ID
N_ITER=${N_ITER:-100000}
N_CHAINS=${N_CHAINS:-4}
SEED_BASE=${SEED_BASE:-20260520}

echo "=== Sim 3 v9-induce: rep $REP_ID, nIter $N_ITER, nChains $N_CHAINS, seed base $SEED_BASE ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "TMPDIR:  $TMPDIR"
echo "Library: $LIB"
echo "Repo:    $REPO"

Rscript "$REPO/inst/hamilton/sim3-v9-induce/run_v9.R" \
        "$REP_ID" "$N_ITER" "$N_CHAINS" "$SEED_BASE"

echo "=== Finished: $(date) ==="
du -hs "$TMPDIR" 2>/dev/null || true
