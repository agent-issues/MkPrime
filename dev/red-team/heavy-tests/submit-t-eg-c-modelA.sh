#!/bin/bash
#SBATCH --job-name=t-eg-c
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/t_eg_c_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/t_eg_c_%j.err
#
# T-EG-C — empirical_geometric Model A SBC.
# Single non-array task (one configuration: Model A, kprimeTruncK=30,
# N_SIM=200, N_ITER=12000, N_TIP=8, N_CHAR=30).
#
# IMPORTANT: pre-build the package once on the login node before submitting:
#
#   cd ${SRC} && Rscript -e 'devtools::load_all(".")'
#
# Otherwise pkgload::load_all races to recompile in the shared src/
# directory (feedback_pkgload_prebuild memory).  This is a single task so
# the race is not as catastrophic as a 6-task array, but it is still wasted
# work on the compute node.
#
# Submit:
#
#   cd ${SRC}
#   sbatch dev/red-team/heavy-tests/submit-t-eg-c-modelA.sh
#
# Output artefacts land in
# ${SRC}/dev/red-team/heavy-tests/kprime-viability/T-EG-C/:
#   - sims.rds, ranks.rds, rank-histograms.png, verdict.txt
# Pull back to local with `scp -r` or `rsync`.

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

echo "[$(date)] T-EG-C empirical_geometric Model A SBC"
echo "  job_id=${SLURM_JOB_ID}"
echo "  TMPDIR=${TMPDIR}"
echo "  R_LIBS=${R_LIBS}"

cd "${SRC}"

# Source-controlled harness path (feedback_slurm_inscript_path memory:
# always invoke ${SRC}/dev/red-team/..., never ${RT}/heavy-tests/...).
Rscript "${SRC}/dev/red-team/heavy-tests/kprime-viability/T-EG-C-modelA-sbc.R"

du -hs "${TMPDIR}" \
  > "${RT}/logs/t_eg_c_${SLURM_JOB_ID}_tmpdir.log" || true

echo "[$(date)] T-EG-C complete"
