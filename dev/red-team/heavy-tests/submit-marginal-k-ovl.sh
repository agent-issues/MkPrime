#!/bin/bash
#SBATCH --job-name=rt-mk-ovl
#SBATCH --partition=shared
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --gres=tmp:4G
#SBATCH --output=/nobackup/%u/mkp-study/red-team/logs/mk-ovl_%j.out
#SBATCH --error=/nobackup/%u/mkp-study/red-team/logs/mk-ovl_%j.err
#
# PR-C / plan §7.2 (Rung 2) — Hamilton posterior-overlap test:
# sampled_k vs marginal_k geometric arm. By Rao-Blackwell the two modes target
# the SAME marginal posterior on (tree_length, rate_log_sd, p); this driver runs
# matched-seed chains across a 2x2 (nTip x nChar) grid x N_REP and KS-tests the
# marginals. Pass bar: KS p > 0.01 on every parameter in every cell.
#
# Single task (the driver loops the 16-cell grid serially, ~10-30 min total at
# 12k iter; it saveRDS()s after every cell so a timeout keeps completed cells).
# A FAIL here, given the RB proof + both-variant prior bit-check, indicates
# MIXING/convergence (try longer chains), not a target-posterior difference.
#
# PRE-BUILD ON THE LOGIN NODE first (feedback_pkgload_prebuild) — the driver
# pkgload::load_all()s the in-tree build; the guard below refuses to recompile.
#
# Submit (runs in parallel with the sampled_k SBC arrays):
#   sbatch submit-marginal-k-ovl.sh

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

mkdir -p "${RT}/logs"

# ---- Stale-/missing-.so guard: never recompile under a scheduled job -------
SO="${SRC}/src/MkPrime.so"
if [ ! -f "${SO}" ]; then
  echo "[ovl] MkPrime.so MISSING — pre-build on the login node first; refusing to run."
  exit 1
fi
if find "${SRC}/src" \( -name '*.cpp' -o -name '*.h' -o -name '*.hpp' \) -newer "${SO}" | grep -q .; then
  echo "[ovl] STALE .so vs src — refusing to recompile under sbatch; pre-build on the login node first."
  exit 1
fi

echo "[$(date)] marginal-k posterior-overlap (T-OVL) job ${SLURM_JOB_ID}"
echo "  TMPDIR=${TMPDIR}  R_LIBS=${R_LIBS}  SRC=${SRC}"

cd "${SRC}"

Rscript "${SRC}/dev/red-team/heavy-tests/marginal-k/T-OVL-sampled-vs-marginal.R"

echo "[$(date)] marginal-k T-OVL complete"
echo "--- T-OVL-verdict.txt ---"
cat "${SRC}/dev/red-team/heavy-tests/marginal-k/T-OVL-verdict.txt" || true
