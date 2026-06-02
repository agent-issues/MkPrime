#!/bin/bash
#SBATCH --job-name=rt-mk-agg-samp
#SBATCH --partition=shared
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=tmp:2G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/mk-agg-samp_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/mk-agg-samp_%j.err
#
# Aggregation step for the sampled_k SBC array (submit-marginal-k-sbc-sampled.sh).
# Submit with an afterany dependency on the matching batch's array job and the
# SAME BATCH value so it reads that batch's shards/OUTDIR:
#
#   A1=$(sbatch --parsable --export=ALL,BATCH=1 submit-marginal-k-sbc-sampled.sh)
#   sbatch --dependency=afterany:${A1} --export=ALL,BATCH=1 submit-marginal-k-agg-sampled.sh
#
# afterany (not afterok): a single dead shard then costs its 5 sims, not the
# whole verdict — the harness' min_good filter tolerates missing shards. Runs
# the SAME harness in AGGREGATE mode: merges OUTDIR/shards/sims-shard-*.rds and
# runs the identical AD/verdict code, so the verdict is identical-by-construction.

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

# ---- Batch -> seedBase + output dir (MUST match the array job) ------------
BATCH="${BATCH:-1}"
case "${BATCH}" in
  1) SEEDBASE=20260528; OUTSUB=sbc-results-sampled-b1 ;;
  2) SEEDBASE=20260901; OUTSUB=sbc-results-sampled-b2 ;;
  *) echo "Unknown BATCH=${BATCH} (expected 1 or 2)"; exit 1 ;;
esac
OUTDIR="${SRC}/dev/red-team/heavy-tests/marginal-k/${OUTSUB}"

mkdir -p "${RT}/logs"

echo "[$(date)] sampled_k SBC AGGREGATE BATCH=${BATCH} (job ${SLURM_JOB_ID})  OUTDIR=${OUTDIR}"
cd "${SRC}"

export MARGINAL_K_SBC_NSHARD=40
export MARGINAL_K_SBC_AGGREGATE=1
export MARGINAL_K_SBC_SEEDBASE=${SEEDBASE}
export MARGINAL_K_SBC_OUTDIR="${OUTDIR}"
Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-sampled-geometric.R"

echo "[$(date)] sampled_k SBC aggregate BATCH=${BATCH} complete"
echo "--- verdict.txt (BATCH=${BATCH}) ---"
cat "${OUTDIR}/verdict.txt" || true
