#!/bin/bash
#SBATCH --job-name=rt-mk-sbc-samp
#SBATCH --partition=shared
#SBATCH --array=0-39
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=tmp:2G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/mk-sbc-samp_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/mk-sbc-samp_%A_%a.err
#
# PR-C / plan §7.3 — Hamilton SBC for the SAMPLED_K geometric arm (Stage 2).
#
# Mirrors submit-marginal-k-sbc.sh but drives T-SBC-sampled-geometric.R
# (likelihoodMode = "sampled_k", kprimeTruncK = K). Validates the Stage-2
# claim that the truncated+renormalised sampled_k prior targets the same
# posterior as marginal_k (which already passed SBC at low p).
#
# TWO-BATCH design (pre-registered criterion, MARGINAL-K-CACHE-002-resume.md):
# the per-run all-3-AD>0.4 strict gate rejects a perfect sampler ~78% of the
# time (0.6^3), so the operative bar is two disjoint seed-batches pooled to
# N=400. Select the batch with --export=ALL,BATCH={1,2}:
#   BATCH=1 -> seedBase 20260528 -> sbc-results-sampled-b1/
#   BATCH=2 -> seedBase 20260901 -> sbc-results-sampled-b2/
# Disjoint OUTDIR per batch is REQUIRED: each array writes shards/sims-shard-
# NNN.rds, so a shared dir would collide when both batches run concurrently.
#
# IMPORTANT — PRE-BUILD ON THE LOGIN NODE BEFORE sbatch (feedback_pkgload_prebuild):
#   cd /nobackup/${USER}/mkp-study/red-team/mkp-source
#   git pull
#   module load r/4.5.1 gcc/14.2
#   R_LIBS_USER=/nobackup/${USER}/mkp-study/red-team/lib:/nobackup/${USER}/mkp-study/lib \
#     Rscript -e 'pkgload::load_all(getwd())'
# The stale-/missing-.so guard below refuses to recompile in the 40-way array.
#
# Submit (both batches + their aggregators run in parallel):
#   A1=$(sbatch --parsable --export=ALL,BATCH=1 submit-marginal-k-sbc-sampled.sh)
#   sbatch --dependency=afterany:${A1} --export=ALL,BATCH=1 submit-marginal-k-agg-sampled.sh
#   A2=$(sbatch --parsable --export=ALL,BATCH=2 submit-marginal-k-sbc-sampled.sh)
#   sbatch --dependency=afterany:${A2} --export=ALL,BATCH=2 submit-marginal-k-agg-sampled.sh

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

# ---- Batch -> seedBase + output dir ---------------------------------------
BATCH="${BATCH:-1}"
case "${BATCH}" in
  1) SEEDBASE=20260528; OUTSUB=sbc-results-sampled-b1 ;;
  2) SEEDBASE=20260901; OUTSUB=sbc-results-sampled-b2 ;;
  *) echo "Unknown BATCH=${BATCH} (expected 1 or 2)"; exit 1 ;;
esac
OUTDIR="${SRC}/dev/red-team/heavy-tests/marginal-k/${OUTSUB}"

mkdir -p "${RT}/logs"
mkdir -p "${OUTDIR}/shards"

# ---- Stale-/missing-.so guard: never recompile in array context -----------
SO="${SRC}/src/MkPrime.so"
if [ ! -f "${SO}" ]; then
  echo "[shard ${SLURM_ARRAY_TASK_ID:-?}] MkPrime.so MISSING — pre-build on the login node first; refusing to run."
  exit 1
fi
if find "${SRC}/src" \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) -newer "${SO}" | grep -q .; then
  echo "[shard ${SLURM_ARRAY_TASK_ID:-?}] STALE .so vs src — refusing to recompile in a 40-way array; pre-build on the login node first."
  exit 1
fi

echo "[$(date)] sampled_k SBC BATCH=${BATCH} seedBase=${SEEDBASE} shard ${SLURM_ARRAY_TASK_ID}"
echo "  array_job=${SLURM_ARRAY_JOB_ID} task=${SLURM_ARRAY_TASK_ID}  OUTDIR=${OUTDIR}"
echo "  TMPDIR=${TMPDIR}  R_LIBS=${R_LIBS}  SRC=${SRC}"

cd "${SRC}"

export MARGINAL_K_SBC_NSHARD=40
export MARGINAL_K_SBC_SHARD=${SLURM_ARRAY_TASK_ID}
export MARGINAL_K_SBC_SEEDBASE=${SEEDBASE}
export MARGINAL_K_SBC_OUTDIR="${OUTDIR}"
Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-sampled-geometric.R"

du -hs "${TMPDIR}" > "${RT}/logs/mk-sbc-samp_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_tmpdir.log" || true
echo "[$(date)] sampled_k SBC BATCH=${BATCH} shard ${SLURM_ARRAY_TASK_ID} complete"
