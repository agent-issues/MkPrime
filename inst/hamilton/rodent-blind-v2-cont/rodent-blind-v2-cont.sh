#!/bin/bash
#SBATCH --job-name=mkp-rod-blind-cont
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2-cont.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2-cont.err
#SBATCH --time=04:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent BLIND Mk' MCMC v2 continuation — extends 200k pilot to 1M total.
# Patches checkpoint mcmc$nIter to 1000000L, then ResumeMkPrime runs
# the remaining ~800k iterations.  Output appends to existing results dir.
# Requires checkpoint at /nobackup/pjjg18/mkp-rodent-blind-v2/rodent-blind-v2.ckp

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex

echo "=== Rodent BLIND v2 CONTINUATION (pilot→1M) ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"

cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/rodent-blind-v2-cont/run_rodent_blind_v2_cont.R \
        "$NEX_FILE"

echo "=== Finished: $(date) ==="
