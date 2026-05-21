#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw4
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v4/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v4/aware/rodent-MkNT-aware.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

# Rodent MkNT AWARE v4: long PT for a trustworthy paper result.
#
# Config rationale:
#  - nChains=4 PT  : PT mode-trap resilience (T-015 confirmed PT cost is
#                    genuinely linear in nChains; no fix; pay the 4x).
#  - nRuns=2 nCore=2: 2 independent PT chains in parallel callr workers,
#                    enough for cross-run Rhat. Smaller queue footprint =
#                    faster SLURM scheduling than nCore=4. (v3's nRuns=4 only
#                    got 1 run through warmup in 20k iter; this gives each
#                    PT chain a full 30k iter budget.)
#  - nIter=30000   : ~10-25k iter for sampling after warmup completes
#                    somewhere in [5000, 20000].
#  - maxWarmup=20000: tight enough that sampling phase always fires before
#                     nIter; loose enough that Geweke has a real shot.
#  - 24h walltime  : ~14h projected per worker (30k iter / 0.6 iter/s on
#                     Hamilton); 24h is safe margin. Per-run checkpoints
#                     support resume if it overruns (parallel-mode resumable
#                     since feat/parallel-checkpointing).
#
# T-014 env-gated drift cadence is the new default /200; explicit here for clarity.

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

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v4/aware/results

echo "=== Rodent MkNT AWARE v4 (nChains=4 PT, nRuns=2, nCore=2): nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"

# Inline-patch: swap outRoot to v4, nRuns/nCore 4 -> 2, maxWarmup 30000 -> 20000.
# nChains stays at 4 (PT). The v1 R script's run_rodent_MkNT.R is the source.
TMPSCRIPT=$(mktemp --suffix=.R)
sed -e 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v4"|g' \
    -e 's|nRuns          = 4L,|nRuns          = 2L,|' \
    -e 's|nCore          = 4L,|nCore          = 2L,|' \
    -e 's|maxWarmup      = 30000L,|maxWarmup      = 20000L,|' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"

grep -q 'mkp-rodent-MkNT-v4' "$TMPSCRIPT" || { echo "outRoot patch FAILED"; exit 1; }
grep -q 'nRuns          = 2L'   "$TMPSCRIPT" || { echo "nRuns patch FAILED"; exit 1; }
grep -q 'nCore          = 2L'   "$TMPSCRIPT" || { echo "nCore patch FAILED"; exit 1; }
grep -q 'maxWarmup      = 20000L' "$TMPSCRIPT" || { echo "maxWarmup patch FAILED"; exit 1; }
echo "Patched script: nChains=4 PT, nRuns=2, nCore=2, maxWarmup=20000, outRoot=mkp-rodent-MkNT-v4"

Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
