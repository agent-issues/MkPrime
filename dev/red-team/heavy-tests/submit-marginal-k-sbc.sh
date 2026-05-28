#!/bin/bash
#SBATCH --job-name=rt-marg-k-sbc
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/marg-k-sbc_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/marg-k-sbc_%j.err
#
# PR-C / plan §7.3 — Hamilton SBC for the marginal-k geometric arm.
# Mirrors `submit-sbc.sh` but single-task (v1 = geometric arm only).
#
# IMPORTANT — PRE-BUILD ON THE LOGIN NODE BEFORE sbatch:
#
#     cd /nobackup/${USER}/mkp-study/red-team/mkp-source
#     module load r/4.5.1
#     module load gcc/14.2 || true
#     R_LIBS_USER=/nobackup/${USER}/mkp-study/red-team/lib \
#       Rscript -e 'devtools::load_all(".")'
#
# Per `feedback_pkgload_prebuild`: pkgload::load_all() inside the Rscript
# this script invokes will re-compile MkPrime if a stale src/MkPrime.so is
# missing — and for an array job this races between tasks and half die
# with "file too short" (sbc v9 lost 3/6 to exactly this). Single-task
# jobs avoid the race but a stale .so still wastes an 8h walltime window
# rebuilding. Pre-build once on the login node so the worktree's src/ has
# a fresh .so before sbatch.
#
# Per `feedback_slurm_inscript_path`: the Rscript path below is
# ${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R
# — NOT ${RT}/heavy-tests/... which was the stale-script bug that silently
# bit sbc v6/v7 (running an out-of-tree old copy of sbc.R).

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

mkdir -p "${RT}/logs"
mkdir -p "${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results-hamilton"

echo "[$(date)] marginal-k SBC (geometric arm)"
echo "  job_id=${SLURM_JOB_ID}"
echo "  TMPDIR=${TMPDIR}"
echo "  R_LIBS=${R_LIBS}"
echo "  SRC=${SRC}"

# Run from the worktree root so pkgload::load_all(".") works AND so the
# in-script OUT_DIR path resolves correctly (it's relative).
cd "${SRC}"

# Use the source-controlled driver — NOT ${RT}/heavy-tests/... (stale copy
# bug — feedback_slurm_inscript_path).
Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R"

# Mirror artefacts into the Hamilton-results dir for easier collection
# (the driver writes to sbc-results/; we copy to sbc-results-hamilton/ so
# the SRC tree stays clean if re-run interactively).
RESULTS_SRC="${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results"
RESULTS_DST="${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results-hamilton"
if [ -d "${RESULTS_SRC}" ]; then
  cp -r "${RESULTS_SRC}"/* "${RESULTS_DST}/" || true
fi

du -hs "${TMPDIR}" > "${RT}/logs/marg-k-sbc_${SLURM_JOB_ID}_tmpdir.log" || true

echo "[$(date)] marginal-k SBC complete"
