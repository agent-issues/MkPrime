#!/bin/bash
#SBATCH --job-name=mkp-rod-blind2
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2-cont2.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2-cont2.err
#SBATCH --time=04:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent BLIND Mk' MCMC v2 second continuation — 1M reached (minESS=38),
# extending to 5M total (~4M additional) to push minESS >= 200.
# Patches checkpoint mcmc$nIter to 5000000L; thin/treeThin unchanged.
# Output appends to existing log/trees files.
# Requires checkpoint at /nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2.ckp

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex

echo "=== Rodent BLIND v2 CONTINUATION 2 (1M->5M) ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/rodent-blind-v2-cont2/run_rodent_blind_v2_cont2.R \
        "$NEX_FILE"

echo "=== Finished: $(date) ==="
