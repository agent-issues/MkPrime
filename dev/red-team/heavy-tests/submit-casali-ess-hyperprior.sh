#!/bin/bash
# submit-casali-ess-hyperprior.sh
#
# SLURM array submitter for casali-ess-hyperprior-vs-gamma.R.
# Each array task runs ONE (cell, prior) MCMC at production scale.
#
# Layout on Hamilton (set up by hand before running this):
#   SRC=/nobackup/pjjg18/mkp-hyperprior-bench/repos/mkp     # worktree at bench/hyperprior-ess
#   LIB=/nobackup/pjjg18/mkp-hyperprior-bench/lib           # MkPrime + AutoPart + coda installed here
#   OUT=/nobackup/pjjg18/mkp-hyperprior-bench/results
#   PARTS=/nobackup/pjjg18/auto-part-skel-pp/repos/auto-part/dev/benchmarks/casali/cache/partitions
#   LOGS=/nobackup/pjjg18/mkp-hyperprior-bench/logs
#
# Index → (cell, prior) mapping is built from casali-ess-cells.tsv:
#   For N cells in the tsv, array indices are 0..(2N-1):
#     task_idx // 2 -> cell row     (0-based, skipping header)
#     task_idx %  2 -> 0 = hyperprior_pooled, 1 = gamma_independent
#
# Submit (after pre-building MkPrime on the LOGIN node first):
#   sbatch --array=0-13 dev/red-team/heavy-tests/submit-casali-ess-hyperprior.sh
#
#SBATCH --job-name=mkp-hp-ess
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=23:59:00
#SBATCH --output=/nobackup/pjjg18/mkp-hyperprior-bench/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-hyperprior-bench/logs/%x_%A_%a.err

set -euo pipefail

# --- Paths (override via env to point at a different scratch tree) ------
SRC="${SRC:-/nobackup/pjjg18/mkp-hyperprior-bench/repos/mkp}"
LIB="${LIB:-/nobackup/pjjg18/mkp-hyperprior-bench/lib}"
OUT="${OUT:-/nobackup/pjjg18/mkp-hyperprior-bench/results}"
PARTS="${PARTS:-/nobackup/pjjg18/auto-part-skel-pp/repos/auto-part/dev/benchmarks/casali/cache/partitions}"
CELLS="${CELLS:-${SRC}/dev/red-team/heavy-tests/casali-ess-cells.tsv}"

mkdir -p "$OUT"

module load r/4.5.1
module load gcc/14.2 || true

export R_LIBS_USER="$LIB"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

# --- Pick cell+prior from array index -----------------------------------
TIDX=${SLURM_ARRAY_TASK_ID:-0}
PRIOR_IDX=$(( TIDX % 2 ))
CELL_IDX=$(( TIDX / 2 ))

# Read non-header line $((CELL_IDX + 2)) — skip TSV header
LINE=$(awk -v n=$((CELL_IDX + 2)) 'NR==n' "$CELLS")
if [ -z "$LINE" ]; then
  echo "No TSV row at index $CELL_IDX (file: $CELLS)" >&2
  exit 2
fi
MATRIX=$(echo "$LINE" | cut -f1)
TREAT=$(echo  "$LINE" | cut -f2)

if [ "$PRIOR_IDX" -eq 0 ]; then
  PRIOR="hyperprior_pooled"
else
  PRIOR="gamma_independent"
fi

OUTFILE="${OUT}/${MATRIX}__${TREAT}__${PRIOR}.rds"

echo "== mkp-hp-ess =="
echo "  job          : ${SLURM_JOB_ID:-local} (array task ${TIDX})"
echo "  host         : $(hostname)"
echo "  date         : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  matrix       : ${MATRIX}"
echo "  treatment    : ${TREAT}"
echo "  prior        : ${PRIOR}"
echo "  out          : ${OUTFILE}"
echo "  R_LIBS_USER  : ${R_LIBS_USER}"
echo "  SRC          : ${SRC}"

cd "$SRC"

# Skip if already produced (idempotent re-submission)
if [ -f "$OUTFILE" ]; then
  echo "Already exists, skipping: $OUTFILE"
  exit 0
fi

Rscript "${SRC}/dev/red-team/heavy-tests/casali-ess-hyperprior-vs-gamma.R" \
  --matrix          "${MATRIX}" \
  --treatment       "${TREAT}" \
  --prior           "${PRIOR}" \
  --partitions-dir  "${PARTS}" \
  --out             "${OUTFILE}" \
  --autopart-lib    "${LIB}" \
  --cores           "${SLURM_CPUS_PER_TASK:-1}"
