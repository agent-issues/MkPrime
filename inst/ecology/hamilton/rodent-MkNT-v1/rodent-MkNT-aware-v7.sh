#!/bin/bash
#SBATCH --job-name=rod-MkNT-aw7
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-MkNT-v7/aware/rodent-MkNT-aware.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-MkNT-v7/aware/rodent-MkNT-aware.err
#SBATCH --time=36:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4

# Rodent MkNT AWARE v7: nChains=16, heat<=0.05, nRuns=4, nCore=4.
#
# v6 diagnosis (60k iter, nChains=8, heatMax=0.25):
#   - Betas: [1.00, 0.82, 0.67, 0.55, 0.45, 0.37, 0.30, 0.25]
#   - round_trip_count=0 for ALL 4 runs
#   - swap_cold median=0 all runs
#   - 4 runs converged to 4 different topologies (RF 44-61% pairwise)
#   - RF vs blind: 100% for all runs (no shared splits)
#   - Run 3 found logLik ~46 units better than other runs (max -5290.9)
#     but its topology is still biologically nonsensical
#   - Conclusion: topology energy barriers too large for beta=0.25
#
# v7 changes:
#   1. nChains=16 (finer ladder over wider beta range)
#   2. heatMax=0.05 -- hot chain is 20x flatter than cold
#      With nChains=16, betas:
#        0.05^(i/15) for i=0..15
#        = [1.000, 0.857, 0.735, 0.630, 0.540, 0.463, 0.397, 0.340,
#           0.291, 0.250, 0.214, 0.183, 0.157, 0.135, 0.116, 0.050]
#      Adjacent Δbeta ≈ 0.13 at cold end, 0.07 at hot end.
#   3. nIter=40000 (wall budget: 16-chain cost ≈ 2x v6 per iter;
#      40k * 2x / 60k * 1x = 1.33x v6's 18.5h ≈ 24.6h; 36h wall is safe)
#   4. maxWarmup=20000 (half nIter, leaving >=20k for sampling)
#   5. All 4 runs start from AdditionTree (default) -- fresh starts
#      so any run that stumbles into run3's better mode is identifiable.
#
# Note: if round_trip_count is still 0 after v7, the barrier is
# >1/0.05 = 20 log-likelihood units per unit beta, and PT is not
# the right tool for this problem.

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKPRIME_ECO_RESYNC_EVERY=200
export MKPRIME_PT_HEAT_MAX=0.05

NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
N_ITER=${N_ITER:-40000}

mkdir -p /nobackup/pjjg18/mkp-rodent-MkNT-v7/aware/results

echo "=== Rodent MkNT AWARE v7 (nChains=16, heat<=0.05): nIter=$N_ITER ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

cd "$MKP_REPO_ROOT"

TMPSCRIPT=$(mktemp --suffix=.R)
sed -e 's|"mkp-rodent-MkNT-v1"|"mkp-rodent-MkNT-v7"|g' \
    -e 's|nChains        = 4L,|nChains        = 16L,|' \
    -e 's|maxWarmup      = 30000L,|maxWarmup      = 20000L,|' \
    inst/ecology/hamilton/rodent-MkNT-v1/run_rodent_MkNT.R > "$TMPSCRIPT"

grep -q 'mkp-rodent-MkNT-v7'   "$TMPSCRIPT" || { echo "outRoot patch FAILED"; exit 1; }
grep -q 'nChains        = 16L'  "$TMPSCRIPT" || { echo "nChains patch FAILED"; exit 1; }
grep -q 'maxWarmup      = 20000L' "$TMPSCRIPT" || { echo "maxWarmup patch FAILED"; exit 1; }
echo "Patched script: nChains=16 PT, nRuns=4, nCore=4, maxWarmup=20000, heatMax=0.05 (via env)"

Rscript "$TMPSCRIPT" "$NEX_FILE" aware "$N_ITER"
rm -f "$TMPSCRIPT"
echo "=== Finished: $(date) ==="
