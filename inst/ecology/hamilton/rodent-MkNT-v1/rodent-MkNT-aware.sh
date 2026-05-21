#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v1/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v1/aware/rodent-MkNT-aware.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Rodent MkNT AWARE: validated MkNT likelihood + ecology-aware layer.
# nRuns=4 nChains=4 nCore=4 nIter=100000 — ~13h walltime (per v4 throughput);
# 24h is comfortable margin. Resumable.

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-100000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v1/aware/results

echo "=== Rodent MkNT AWARE: nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"
Rscript inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R \
        "$NEX_FILE" aware "$N_ITER"
echo "=== Finished: $(date) ==="
