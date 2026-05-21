#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw3
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v3/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v3/aware/rodent-MkNT-aware.err
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4

# Rodent MkNT AWARE v3: nChains=1 (no PT) + nRuns=4 (Rhat across runs) +
# nCore=4 (parallel via callr). Workaround for the PT amplification
# diagnosed in T-014 (aware nChains=4 PT is ~4x slower on Hamilton than
# nChains=1; 30k iter warmup with PT projects ~12h, blowing walltime).
# nChains=1 nRuns=4 gives 4 independent posterior chains for Rhat from a
# single Hamilton submission. Loses PT mode-trap mitigation; on rodent with
# 4 independent AdditionTree starts the mode-trap risk should be tolerable.
#
# Projected wall: ~3-4h based on local 10.6 iter/s aware nChains=1 single-chain,
# scaled ~4x slower on Hamilton => ~2.5 iter/s per worker. 4 parallel workers
# at 20k iter each => ~2.2h. 12h walltime leaves wide margin.

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# T-014: explicit resync cadence; new default is /200 anyway but be explicit.
export MKPRIME_ECO_RESYNC_EVERY=200

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-20000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v3/aware/results

echo "=== Rodent MkNT AWARE v3 (nChains=1, nRuns=4, nCore=4): nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"

# Inline-patch the v1 R script: swap outRoot to v3 + nChains to 1.
TMPSCRIPT=$(mktemp --suffix=.R)
sed -e 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v3"|g' \
    -e 's|nChains        = 4L,|nChains        = 1L,|' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"
# Verify the substitutions landed (sed silently no-ops if pattern misses).
grep -q 'mkp-rodent-MkNT-v3' "$TMPSCRIPT" || { echo "outRoot patch FAILED"; exit 1; }
grep -q 'nChains        = 1L' "$TMPSCRIPT" || { echo "nChains patch FAILED"; exit 1; }
echo "Patched script: nChains=1, outRoot=mkp-rodent-MkNT-v3"

Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
