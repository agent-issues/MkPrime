#!/bin/bash
#SBATCH --job-name=mkp-rod-aware
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.err
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent AWARE Mk' MCMC v2 — paired with rodent-blind-v2 (blind).
# Reuses the mkp-sim3-multirep-v3 lib and MkPrime repo on Hamilton.
# 1M-iter run (local pilot terminated at 127k still in warmup).
# Nexus reused from mkp-rodent-blind-v2/data/ — not re-uploaded.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-1000000}

mkdir -p /nobackup/pjjg18/mkp-rodent-aware-v2/results

echo "=== Rodent AWARE v2: nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/rodent-aware-v2/run_rodent_aware_v2.R \
        "$NEX_FILE" "$N_ITER"

echo "=== Finished: $(date) ==="
