#!/bin/bash
#SBATCH --job-name=rt-sbc-mixed
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/sbc_mixed_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/sbc_mixed_%j.err
#
# Red-team A3 — SBC on mixed-partition (neo + trans) MkNT data.
# Validates partition-rate normalisation fix (main 2026-05-27):
#   neoScale = r/(1+r)*n_total/n_neo;  transScale = 1/(1+r)*n_total/n_trans
# Single job (not array): one arm, MkNT_mixed, 4 neo + 8 trans chars.
#
# IMPORTANT: pre-build the package once before submitting:
#   cd ${SRC} && Rscript -e 'devtools::load_all(".")'
# (avoids race in shared src/ directory if re-submitted alongside other jobs)

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

echo "[$(date)] SBC mixed arm: MkNT_mixed"
echo "  job_id=${SLURM_JOB_ID}"
echo "  TMPDIR=${TMPDIR}"
echo "  R_LIBS=${R_LIBS}"

cd "${SRC}"

Rscript "${SRC}/dev/red-team/heavy-tests/sbc-mixed.R" \
  --full \
  --out "${RT}/results/sbc-mixed" \
  --seed 20260528

du -hs "${TMPDIR}" > "${RT}/logs/sbc_mixed_${SLURM_JOB_ID}_tmpdir.log" || true

echo "[$(date)] MkNT_mixed complete"
