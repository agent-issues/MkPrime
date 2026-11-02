#!/bin/bash
#SBATCH --job-name=ss-longrich
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --array=1-16
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/ss-longrich_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/ss-longrich_%A_%a.err
#
# Reference log Z for Longrich2010 (#66): tasks 1-8 run 20 stones x 50 000
# moves, tasks 9-16 x 10 000, one descending and one ascending ladder each.
# At about 2 ms a move, a 50 000 task (2.5 million moves) takes 1.5-3 h.
#
#   BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)
#   sbatch --dependency=afterok:$BUILD dev/stepping-stone/submit-longrich-reference.sh

set -euo pipefail

module load r/4.5.1
module load gcc/14.2 || true

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
OUT=${RT}/stepping-stone
export R_LIBS_USER="${RT}/lib:${PROJECT}/lib"
export R_LIBS="${RT}/lib:${PROJECT}/lib"
mkdir -p "${OUT}"

cd "${SRC}"
if [ "${SLURM_ARRAY_TASK_ID}" -le 8 ]; then N_ITER=50000; else N_ITER=10000; fi
Rscript "${SRC}/dev/stepping-stone/longrich-reference.R" \
  "${SLURM_ARRAY_TASK_ID}" "${N_ITER}" "${OUT}"
