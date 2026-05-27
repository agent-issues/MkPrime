#!/bin/bash
#SBATCH --job-name=rt-sbc
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/sbc_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/sbc_%A_%a.err
#
# Red-team Lane D1 SBC heavy test — 6 arms in parallel array tasks.
# Authored 2026-05-26; v10 (2026-05-27): reverted Option γ — kSim=2 for all
# arms; kPrime_pooled and p both excluded as structural. See sbc.md +
# dev/red-team/findings.md::SBC-KPRIME-STRUCTURAL.
#
# IMPORTANT: pre-build the package once before submitting this script:
#   cd ${SRC} && Rscript -e 'devtools::load_all(".")'
# Otherwise 6 array tasks race-compile in the shared src/ directory and
# half die with "MkPrime.so: file too short" (v9 lost 3/6 to this race).

set -euo pipefail

module load r/4.5.1
module load gcc/14.2 || true

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PROJECT=/nobackup/${USER}/mkp-study
RT=${PROJECT}/red-team
SRC=${RT}/mkp-source
# Chained libpath: fresh MkPrime first, deps lib second.
export R_LIBS_USER="${RT}/lib:${PROJECT}/lib"
export R_LIBS="${RT}/lib:${PROJECT}/lib"

ARMS=(
  "MkNT_geometric"
  "Mkp_geometric"
  "Mkp_beta_geometric"
  "Mkp_empirical_geometric"
  "Mkp_logseries"
  "MkNT_logseries"
)

ARM_NAME="${ARMS[${SLURM_ARRAY_TASK_ID}]}"

echo "[$(date)] SBC arm: ${ARM_NAME}"
echo "  job_id=${SLURM_JOB_ID} array_task=${SLURM_ARRAY_TASK_ID}"
echo "  TMPDIR=${TMPDIR}"
echo "  R_LIBS=${R_LIBS}"

# Run from the worktree root so pkgload::load_all(".") works AND
# so that any package-relative file paths in sbc.R resolve correctly.
cd "${SRC}"

# Use the source-controlled sbc.R, NOT the stale out-of-tree
# ${RT}/heavy-tests/sbc.R (which is not updated when the source tree changes).
# All v6/v7 runs prior to 2026-05-27 silently used that stale copy.
Rscript "${SRC}/dev/red-team/heavy-tests/sbc.R" \
  --full \
  --arm "${ARM_NAME}" \
  --out "${RT}/results/sbc-v10" \
  --seed $((20260528 + SLURM_ARRAY_TASK_ID))

du -hs "${TMPDIR}" > "${RT}/logs/sbc_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_tmpdir.log" || true

echo "[$(date)] Arm ${ARM_NAME} complete"
