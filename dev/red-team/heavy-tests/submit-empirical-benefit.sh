#!/bin/bash
#SBATCH --job-name=eg-benefit
#SBATCH --partition=shared
#SBATCH --time=02:00:00
#SBATCH --array=1-200
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/eg_benefit_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/eg_benefit_%A_%a.err
#
# Empirical-vs-geometric k' prior study (#179): 4 configs x 50 datasets, each
# fitted under both priors (empirical-benefit/mcmc.R). About 13 min a task
# (two fits of 2 x 100k iterations, 30 tips x 100 characters): ~45 CPU-hours.
#
# Needs TreeDist in ${R_LIBS}. Build the shared tree first:
#
#   cd ${SRC}
#   BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)
#   sbatch --dependency=afterok:$BUILD \
#     dev/red-team/heavy-tests/submit-empirical-benefit.sh
#
# A task whose output exists is skipped, so the array can be resubmitted
# after a walltime kill. Summarize with
#   Rscript dev/red-team/heavy-tests/empirical-benefit/aggregate.R ${OUT}

set -euo pipefail

module load r/4.5.1
module load gcc/14.2 || true

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
OUT=${RT}/empirical-benefit
export R_LIBS_USER="${RT}/lib:${PROJECT}/lib"
export R_LIBS="${RT}/lib:${PROJECT}/lib"
export MKP_EB_OUT="${OUT}"
mkdir -p "${OUT}"

echo "[$(date)] eg-benefit task ${SLURM_ARRAY_TASK_ID} job ${SLURM_JOB_ID}"
cd "${SRC}"
Rscript "${SRC}/dev/red-team/heavy-tests/empirical-benefit/mcmc.R" \
  "${SLURM_ARRAY_TASK_ID}"
echo "[$(date)] done"
