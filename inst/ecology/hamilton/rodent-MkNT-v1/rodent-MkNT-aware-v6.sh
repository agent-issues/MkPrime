#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw6
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v6/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v6/aware/rodent-MkNT-aware.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Rodent MkNT AWARE v6: nChains=8, heat<=0.25, nRuns=4, nCore=4.
#
# v5 diagnosis (PT-RT-001 fix applied):
#   - Betas correctly spread: [1.00, 0.794, 0.630, 0.500]
#   - But: cold-end swap rate 0.3%, round_trip_count=0 for BOTH runs at 30k iter
#   - Inter-run RF = 72 / 114 (63%): two runs, same logLik (~-5371), different topology
#   - Ladder warning: "heat pinned to 0.50; consider nChains >= 6"
#   - Conclusion: 4-chain ladder is too coarse; swap rate is MUCH too low even at heat=0.5
#
# v6 changes:
#   1. nChains=8 (finer ladder, better adjacent-pair swap rates)
#   2. heatMax=0.25 in .AdaptTemperatures() — hot chain beta=0.25 means a much
#      flatter likelihood; warm-start diagnostic (job 17272483) will confirm
#      whether the blind topology is accessible. Pending that result, 0.25 is
#      the recommended floor based on cold-end swap rate analysis.
#   3. nRuns=4, nCore=4 (4 independent PT chains for Rhat diagnostics)
#   4. nIter=60000 (aware nChains=8 estimated ~12h for 60k iter; fits 24h)
#   5. maxWarmup=25000 (tight upper bound to leave >=35k iter for sampling)
#   6. All 4 runs start from AdditionTree (default) — warm-start from blind
#      consensus handled separately if v6 still mode-traps.
#
# heatMax=0.25 betas with nChains=8:
#   [1.00, 0.89, 0.79, 0.71, 0.63, 0.56, 0.50, 0.25]
#   Wait -- the ladder builder is:
#     betas[i] = heat^(i / (nChains-1))  for i = 0..nChains-1
#   With heat=0.25, nChains=8:
#     betas = [1.00^0, 0.25^(1/7), 0.25^(2/7), ... 0.25^1]
#           = [1.000, 0.820, 0.672, 0.551, 0.451, 0.370, 0.303, 0.250]
#   Much better coverage of [0.25, 1.0] than the old [0.88, 1.0] range.
#
# Expected performance on Hamilton (scaling from v5):
#   v5: 4.83h for 30k iter at nChains=4  => ~0.58 s/iter/chain-unit
#   v6: nChains=8 => ~1.16 s/iter => 60k iter => ~19h (fits 24h wall)
#   With nCore=4 all runs parallel, wall stays ~19h.
#
# DO NOT SUBMIT until warm-start job 17272483 completes and confirms
# the blind topology is accessible (mode-trap) vs genuinely preferred
# (aware posterior different). If warm-start shows blind is NOT accessible,
# v6 needs a different design (initialise from diverse starts, consider
# replica-exchange with more runssss).

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKPRIME_ECO_RESYNC_EVERY=200
# PT-RT-001 v6: clamp hot chain beta at 0.25 (lower = hotter = more mode-hopping).
# The adapter will push heat toward this ceiling; starting at heat=0.2 (in
# run_rodent_MkNT.R) means the hot chain starts BELOW this ceiling so early
# warmup can escape basins before swap rates have adapted.
export MKPRIME_PT_HEAT_MAX=0.25

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-60000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v6/aware/results

echo "=== Rodent MkNT AWARE v6 (nChains=8, heat<=0.25): nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"

TMPSCRIPT=$(mktemp --suffix=.R)
sed -e 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v6"|g' \
    -e 's|nChains        = 4L,|nChains        = 8L,|' \
    -e 's|maxWarmup      = 30000L,|maxWarmup      = 25000L,|' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"

grep -q 'mkp-rodent-MkNT-v6' "$TMPSCRIPT" || { echo "outRoot patch FAILED"; exit 1; }
grep -q 'nChains        = 8L'   "$TMPSCRIPT" || { echo "nChains patch FAILED"; exit 1; }
grep -q 'maxWarmup      = 25000L' "$TMPSCRIPT" || { echo "maxWarmup patch FAILED"; exit 1; }
echo "Patched script: nChains=8 PT, nRuns=4, nCore=4, maxWarmup=25000, heatMax=0.25 (via env)"

Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
