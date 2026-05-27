#!/bin/bash
#SBATCH --job-name=rt-swap-n6
#SBATCH --partition=shared
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/swap-db-n6_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/swap-db-n6_%j.err
#
# D2(b) escalation: n=6 (105 topologies) with the same 1M iter / thin=200
# budget gives ~48 samples / bucket — sufficient power for chi² to detect
# the math-prover predicted SWAP-001/SWAP-002 bias if it exists.
# n=7 (945 topos / ~5 per bucket) is underpowered at this iter count;
# skipped per advisor.

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

echo "[$(date)] subtree-swap-db n=6 (1M iter)"
cd "${SRC}"
Rscript "${RT}/heavy-tests/subtree-swap-db.R" --n 6
echo "[$(date)] complete"
