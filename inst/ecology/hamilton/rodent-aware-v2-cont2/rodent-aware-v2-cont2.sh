#!/bin/bash
#SBATCH --job-name=mkp-rod-awar2
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2-cont2.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2-cont2.err
#SBATCH --time=36:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent AWARE Mk' MCMC v2 second continuation — 1M reached (minESS=88),
# extending to 2.3M total (~1.3M additional) to push minESS >= 200.
# Patches checkpoint mcmc$nIter to 2300000L; thin/treeThin unchanged.
# Output appends to existing log/trees files.
# Requires checkpoint at /nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.ckp

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex

echo "=== Rodent AWARE v2 CONTINUATION 2 (1M->2.3M) ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/rodent-aware-v2-cont2/run_rodent_aware_v2_cont2.R \
        "$NEX_FILE"

echo "=== Finished: $(date) ==="
