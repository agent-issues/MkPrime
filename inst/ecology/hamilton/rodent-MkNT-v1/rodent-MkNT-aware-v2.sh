#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw2
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v2/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v2/aware/rodent-MkNT-aware.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Rodent MkNT AWARE v2: smaller config after the v1 stall (job 17254914).
# Hamilton lib has now been rebuilt with the T-012/T-013 ecology CL cache
# patches (commit 76023d3); local rodent timing showed aware = 10.6 iter/s
# single-chain serial, 21x ratio over blind — the expected slowdown.
#
# Config: nIter=20000 (vs 100000 in v1), nChains=4 PT, nCore=4 parallel runs.
# Projected wall: ~2-3h per worker. 12h walltime gives plenty of margin.
# Output dir kept distinct from v1 so the failed v1 logs are untouched.

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-20000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v2/aware/results

echo "=== Rodent MkNT AWARE v2: nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"
# Override the v1 R script's hardcoded outRoot via MKP_OUT env var. The v1
# script reads outRoot from /nobackup/$USER/mkp-rodent-MkNT-v1/$mode; we
# point to v2 via a small inline wrapper that exec's the v1 R but with the
# outRoot path pre-created and a symlink so the v1 script writes to v2 paths.
# Simpler: copy the R, change one line.
TMPSCRIPT=$(mktemp --suffix=.R)
sed 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v2"|g' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"
Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
