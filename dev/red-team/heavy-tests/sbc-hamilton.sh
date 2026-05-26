#!/bin/bash
#SBATCH --job-name=sbc-mkprime
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/logs/sbc_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/logs/sbc_%A_%a.err
#
# Lane D1 — Simulation-Based Calibration heavy test, Hamilton submission.
#
# One array task per arm (6 arms). Each task runs 200 SBC simulations,
# fixed-topology MCMC of ~6000 iter post warmup (L = 100 thinned).
#
# Resource ask justification:
#   - Per-sim wall: ~30–60 s at N_TIP=8, N_CHAR=30, N_ITER=6000.
#   - 200 sims/arm × 45 s = ~2.5 h per arm at single thread.
#   - 8 h walltime budgeted (project memory: feedback_resumable_runs).
#   - 8 GB memory: MkPrime state for 8-tip × 30-char data is ~tens of MB;
#     8 GB is conservative for buffers + libraries.
#   - 4 GB scratch: per-sim checkpoint files.
#
# Usage:
#   sbatch /nobackup/$USER/mkp-study/sbc-hamilton.sh
#
# Or run a single arm interactively to sanity-check before sbatch:
#   srun --partition=test.q --time=00:05:00 --mem=4G --pty bash
#   module load r/4.5.1
#   cd /nobackup/$USER/mkp-study
#   ARM_NAME=Mkp_geometric Rscript sbc.R --quick --arm $ARM_NAME

set -euo pipefail

module load r/4.5.1

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT_ROOT=/nobackup/${USER}/mkp-study
LIB=${PROJECT_ROOT}/lib
.libPaths_addendum=".libPaths(c('${LIB}', .libPaths()))"
RESULTS=${PROJECT_ROOT}/sbc-results
mkdir -p "${RESULTS}" "${PROJECT_ROOT}/logs"

# Arm names — match ALL_ARMS in sbc.R
ARMS=(
  "MkNT_geometric"
  "Mkp_geometric"
  "Mkp_beta_geometric"
  "Mkp_empirical_geometric"
  "Mkp_logseries"
  "MkNT_logseries"
)

ARM_NAME="${ARMS[${SLURM_ARRAY_TASK_ID}]}"

echo "[$(date)] Starting SBC arm: ${ARM_NAME}"
echo "  job_id=${SLURM_JOB_ID} array_task=${SLURM_ARRAY_TASK_ID}"
echo "  TMPDIR=${TMPDIR}"
echo "  results=${RESULTS}/${ARM_NAME}"

cd "${PROJECT_ROOT}"

# Pass --full to run at production scale (N_SIM=200, N_ITER=6000)
Rscript sbc.R \
  --full \
  --arm "${ARM_NAME}" \
  --out "${RESULTS}" \
  --seed $((20260526 + SLURM_ARRAY_TASK_ID))

# Record per-task TMPDIR usage for capacity planning
du -hs "${TMPDIR}" > "${PROJECT_ROOT}/logs/sbc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_tmpdir.log" || true

echo "[$(date)] Arm ${ARM_NAME} complete"
