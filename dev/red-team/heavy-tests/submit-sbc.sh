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
# The 6 array tasks share one src/ tree, and pkgload::load_all() compiles into
# it in place, so an unbuilt tree means 6 concurrent compiles corrupting each
# other's objects -- v9 lost 3/6 tasks to "MkPrime.so: file too short" (#15).
#
# sbc.R now REFUSES to compile inside an array task, so a forgotten pre-build
# stops the job immediately instead of producing a partial array. Submit the
# build as a dependency and SLURM enforces the ordering:
#
#   BUILD=$(sbatch --parsable dev/red-team/heavy-tests/submit-build.sh)
#   sbatch --dependency=afterok:$BUILD dev/red-team/heavy-tests/submit-sbc.sh
#
# Building by hand on the login node still works:
#   cd ${SRC} && Rscript -e 'pkgload::load_all(getwd())'

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

# ---- Aggregate the per-arm verdicts -----------------------------------------
# Each task writes verdict-<arm>.txt; the run-level verdict.txt is produced
# only by the reduction below, so it can never be one task's output wearing the
# aggregate's name (#16). Submit it after the array:
#
#   sbatch --dependency=afterany:$ARRAY_JOB_ID --wrap #     "Rscript ${SRC}/dev/red-team/heavy-tests/aggregate-verdicts.R #        ${RT}/results/sbc-v10 ${ARMS[*]}"
