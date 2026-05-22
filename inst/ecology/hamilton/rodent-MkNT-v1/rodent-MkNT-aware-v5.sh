#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw5
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v5/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v5/aware/rodent-MkNT-aware.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Rodent MkNT AWARE v5: same config as v4 but with the fixed PT temperature
# adapter (PT-RT-001).
#
# v4 diagnostic showed adapter narrowed the ladder to [1, 0.96, 0.92, 0.88]
# - heated chains were barely heated, so PT couldn't bridge modes. Two PT
# runs converged on DIFFERENT topologies, neither close to blind, both far
# from the Fabre reference.
#
# v5 fix (R/RunMkPrime.R::.AdaptTemperatures): clamp upper bound on heat
# tightened from 0.95 to 0.5, so the hot chain is always meaningfully
# heated. Also surfaces a warning if the clamp binds AND swap rates are
# low (recommends more chains). Local sanity check on 500 iter showed
# final betas [1, 0.61, 0.37, 0.22] and 1 round trip - real PT mixing.
#
# round_trip_count is now exposed in the result; this is the canonical PT
# diagnostic.

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKPRIME_ECO_RESYNC_EVERY=200

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-30000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v5/aware/results

echo "=== Rodent MkNT AWARE v5 (PT-RT-001 fix): nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"

TMPSCRIPT=$(mktemp --suffix=.R)
sed -e 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v5"|g' \
    -e 's|nRuns          = 4L,|nRuns          = 2L,|' \
    -e 's|nCore          = 4L,|nCore          = 2L,|' \
    -e 's|maxWarmup      = 30000L,|maxWarmup      = 20000L,|' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"

grep -q 'mkp-rodent-MkNT-v5' "$TMPSCRIPT" || { echo "outRoot patch FAILED"; exit 1; }
grep -q 'nRuns          = 2L'   "$TMPSCRIPT" || { echo "nRuns patch FAILED"; exit 1; }
grep -q 'nCore          = 2L'   "$TMPSCRIPT" || { echo "nCore patch FAILED"; exit 1; }
grep -q 'maxWarmup      = 20000L' "$TMPSCRIPT" || { echo "maxWarmup patch FAILED"; exit 1; }
echo "Patched script: nChains=4 PT, nRuns=2, nCore=2, maxWarmup=20000"
echo "(PT clamp fix is in the installed lib at $MKP_LIB)"

Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
