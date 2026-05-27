#!/bin/bash
#SBATCH --job-name=rt-swap-db
#SBATCH --partition=shared
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/swap-db_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/swap-db_%j.err
#
# Red-team Lane D2(b): β=0 detailed-balance test for GibbsSubtreeSwap and
# WeightedSubtreeSwap at n=5, 1M iter. See tree-ess-and-swap.md.

set -euo pipefail

module load r/4.5.1
module load gcc/14.2 || true

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
export R_LIBS_USER="${RT}/lib:${PROJECT}/lib"
export R_LIBS="${RT}/lib:${PROJECT}/lib"

echo "[$(date)] subtree-swap-db full (n=5, 1M iter)"

cd "${SRC}"
Rscript "${RT}/heavy-tests/subtree-swap-db.R"

echo "[$(date)] complete"
