#!/bin/bash
#SBATCH --job-name=mkp-rod-v4-aw
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-v4/aware/rodent-v4-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-v4/aware/rodent-v4-aware.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Rodent AWARE v4 PRODUCTION.
# nRuns=4 nChains=4 nCore=4 nIter=100000 — ~13h walltime per the v4 validation
# throughput (~7.5k iter/h cold-chain with nChains=4 PT). Parallel runs are
# resumable (per-run ckps under feat/parallel-checkpointing) so 24h walltime
# is a comfortable margin; if walltime is hit, resubmit the same script and
# ResumeMkPrime() synthesises the master checkpoint automatically.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-100000}

mkdir -p /nobackup/pjjg18/mkp-rodent-v4/aware/results

echo "=== Rodent AWARE v4 PRODUCTION: nIter=$N_ITER nRuns=4 nChains=4 nCore=4 ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/rodent-v4/run_rodent_aware_v4.R \
        "$NEX_FILE" "$N_ITER"

echo "=== Finished: $(date) ==="
