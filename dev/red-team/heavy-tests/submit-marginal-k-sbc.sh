#!/bin/bash
#SBATCH --job-name=rt-marg-k-sbc
#SBATCH --partition=shared
#SBATCH --array=0-39
#SBATCH --time=00:15:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --gres=tmp:2G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/marg-k-sbc_%A_%a.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/marg-k-sbc_%A_%a.err
#
# PR-C / plan §7.3 — Hamilton SBC for the marginal-k geometric arm.
#
# ARRAY (40 shards x 5 sims = 200) — replaces the single 8 h task. Rationale:
# under fair-share depression (many of the user's own jobs running) a single
# long job sits at (Priority) and cannot be backfilled into idle CPUs that the
# scheduler reserves for higher-priority jobs. A ~3-min/task array backfills
# into those gaps. Each task runs a contiguous block of sims via the harness'
# MARGINAL_K_SBC_NSHARD / _SHARD env vars (seeds seedBase+i unchanged, so the
# sharded run reproduces the monolith's 200 sims EXACTLY). A separate
# aggregate job (submit-marginal-k-agg.sh, afterany dependency) merges the
# shards and computes the verdict.
#
# IMPORTANT — PRE-BUILD ON THE LOGIN NODE BEFORE sbatch (per
# feedback_pkgload_prebuild): the Rscript below runs pkgload::load_all(),
# which recompiles MkPrime if src/MkPrime.so is stale/missing. With 40 tasks
# that recompile races and clobbers the shared src/ (sbc v9 lost 3/6 exactly
# this way). Pre-build once so the .so is fresh:
#
#     cd /nobackup/${USER}/mkp-study/red-team/mkp-source
#     git pull
#     module load r/4.5.1 gcc/14.2
#     R_LIBS_USER=/nobackup/${USER}/mkp-study/red-team/lib:/nobackup/${USER}/mkp-study/lib \
#       Rscript -e 'pkgload::load_all(getwd())'   # NOT devtools (not installed)
#     mkdir -p dev/red-team/heavy-tests/marginal-k/sbc-results/shards
#
# The stale-.so guard below is belt-and-braces: if the pre-build did not take,
# every task refuses to run rather than 40 of them racing to recompile.
#
# Per feedback_slurm_inscript_path: the Rscript path is ${SRC}/dev/red-team/...
# NOT ${RT}/heavy-tests/... (the stale out-of-tree copy that bit sbc v6/v7).

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
# Submit-side shard dir (idempotent; harness also dir.create()s defensively).
mkdir -p "${SRC}/dev/red-team/heavy-tests/marginal-k/sbc-results/shards"

# ---- Stale-/missing-.so guard (advisor): never recompile in array context --
SO="${SRC}/src/MkPrime.so"
if [ ! -f "${SO}" ]; then
  echo "[shard ${SLURM_ARRAY_TASK_ID:-?}] MkPrime.so MISSING — pre-build on the login node first; refusing to run."
  exit 1
fi
if find "${SRC}/src" \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) -newer "${SO}" | grep -q .; then
  echo "[shard ${SLURM_ARRAY_TASK_ID:-?}] STALE .so vs src sources — refusing to recompile in a 40-way array; pre-build on the login node first."
  exit 1
fi

echo "[$(date)] marginal-k SBC shard ${SLURM_ARRAY_TASK_ID} (geometric arm)"
echo "  array_job=${SLURM_ARRAY_JOB_ID} task=${SLURM_ARRAY_TASK_ID}"
echo "  TMPDIR=${TMPDIR}  R_LIBS=${R_LIBS}  SRC=${SRC}"

cd "${SRC}"

export MARGINAL_K_SBC_NSHARD=40
export MARGINAL_K_SBC_SHARD=${SLURM_ARRAY_TASK_ID}
Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-SBC-marginal-geometric.R"

du -hs "${TMPDIR}" > "${RT}/logs/marg-k-sbc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}_tmpdir.log" || true
echo "[$(date)] marginal-k SBC shard ${SLURM_ARRAY_TASK_ID} complete"
