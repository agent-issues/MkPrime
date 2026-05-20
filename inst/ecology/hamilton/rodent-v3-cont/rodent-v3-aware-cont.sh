#!/bin/bash
#SBATCH --job-name=mkp-rod-v3-aw-cont
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-v3/aware/rodent-v3-aware-cont.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-v3/aware/rodent-v3-aware-cont.err
#SBATCH --time=12:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

# Rodent AWARE Mk' MCMC v3 — CONTINUATION.
# Resumes from /nobackup/pjjg18/mkp-rodent-v3/aware/rodent-aware-v3.ckp.
# At previous timeout: iter ~346k of 1M, minESS=27. Needs ~7x more
# samples to reach minESS=200, so multiple back-to-back continuations
# may be necessary (queue cont2 via --dependency=afterany).
# Appends to the existing log/checkpoint in place.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex

echo "=== Rodent AWARE v3 CONTINUATION ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"
echo "Library: $MKP_LIB"
echo "Nexus:   $NEX_FILE"
echo "SLURM job: ${SLURM_JOB_ID:-?}"

cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/rodent-v3-cont/run_rodent_aware_v3_cont.R \
        "$NEX_FILE"

echo "=== Finished: $(date) ==="
