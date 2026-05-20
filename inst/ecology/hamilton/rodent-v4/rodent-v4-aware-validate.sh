#!/bin/bash
#SBATCH --job-name=mkp-rod-v4-val
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-v4/aware-validate/rodent-v4-aware-validate.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-v4/aware-validate/rodent-v4-aware-validate.err
#SBATCH --time=08:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1

# Rodent AWARE v4 VALIDATION — short serial PT run to confirm the streaming
# fix (commit 326bb13) holds end-to-end on rodent data. Resumable.
# Production parallel run follows (rodent-v4-aware.sh) only after this passes.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-20000}

mkdir -p /nobackup/pjjg18/mkp-rodent-v4/aware-validate/results

echo "=== Rodent AWARE v4 VALIDATION: nIter=$N_ITER nChains=4 PT serial ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/rodent-v4/run_rodent_aware_v4_validate.R \
        "$NEX_FILE" "$N_ITER"

echo "=== Finished: $(date) ==="
