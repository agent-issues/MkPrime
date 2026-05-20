#!/bin/bash
#SBATCH --job-name=mkp-rod-aware-cont
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2-cont.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2-cont.err
#SBATCH --time=24:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent AWARE Mk' MCMC v2 continuation — extends 371k pilot to 1M total.
# The checkpoint already has mcmc$nIter = 1000000L; ResumeMkPrime runs
# the remaining ~629k iterations.  Output appends to existing log/trees.
# Requires checkpoint at /nobackup/pjjg18/mkp-rodent-aware-v2/rodent-aware-v2.ckp

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex

echo "=== Rodent AWARE v2 CONTINUATION (371k->1M) ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/rodent-aware-v2-cont/run_rodent_aware_v2_cont.R \
        "$NEX_FILE"

echo "=== Finished: $(date) ==="
