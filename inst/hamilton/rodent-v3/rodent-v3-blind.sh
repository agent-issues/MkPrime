#!/bin/bash
#SBATCH --job-name=mkp-rod-v3-bl
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-v3/blind/rodent-v3-blind.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-v3/blind/rodent-v3-blind.err
#SBATCH --time=04:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent BLIND Mk' MCMC v3 — auto-expSteps from parsimony (commit bc24d52).
# Reuses the mkp-sim3-multirep-v3 lib and MkPrime repo on Hamilton.
# 1M-iter run; resume via checkpoint if ESS < 200 (extend walltime).

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-1000000}

mkdir -p /nobackup/pjjg18/mkp-rodent-v3/blind/results

echo "=== Rodent BLIND v3: nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/rodent-v3/run_rodent_blind_v3.R \
        "$NEX_FILE" "$N_ITER"

echo "=== Finished: $(date) ==="
