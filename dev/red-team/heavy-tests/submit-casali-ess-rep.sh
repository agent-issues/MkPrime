#!/bin/bash
# submit-casali-ess-rep.sh
#
# Replication runs for Burns2011/T4 and CarranoSampson2008/T2a at two
# fresh seeds under each prior. Validates whether the class3 σ_c stuck
# pattern observed in the 17304194 run on Burns2011/T4 is a real
# pathology of the non-centred parameterisation on K=3 balanced data,
# or a single-seed artefact. CarranoSampson2008/T2a is the K=3 control
# that was "neutral" in the first batch.
#
# Driven by casali-ess-rep-tasks.tsv (8 rows). Submit:
#   sbatch --array=0-7 dev/red-team/heavy-tests/submit-casali-ess-rep.sh
#
#SBATCH --job-name=mkp-hp-ess-rep
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=23:59:00
#SBATCH --output=/nobackup/pjjg18/mkp-hyperprior-bench/logs/%x_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-hyperprior-bench/logs/%x_%A_%a.err

set -euo pipefail

SRC="${SRC:-/nobackup/pjjg18/mkp-hyperprior-bench/repos/mkp}"
LIB="${LIB:-/nobackup/pjjg18/mkp-hyperprior-bench/lib}"
OUT="${OUT:-/nobackup/pjjg18/mkp-hyperprior-bench/results}"
PARTS="${PARTS:-/nobackup/pjjg18/auto-part-skel-pp/repos/auto-part/dev/benchmarks/casali/cache/partitions}"
TASKS="${TASKS:-${SRC}/dev/red-team/heavy-tests/casali-ess-rep-tasks.tsv}"

mkdir -p "$OUT"
module load r/4.5.1
module load gcc/14.2 || true

export R_LIBS_USER="$LIB"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}

TIDX=${SLURM_ARRAY_TASK_ID:-0}
LINE=$(awk -v n=$((TIDX + 2)) 'NR==n' "$TASKS")
if [ -z "$LINE" ]; then
  echo "No TSV row at index $TIDX (file: $TASKS)" >&2
  exit 2
fi
MATRIX=$(echo "$LINE" | cut -f1)
TREAT=$(echo  "$LINE" | cut -f2)
PRIOR=$(echo  "$LINE" | cut -f3)
SEED=$(echo   "$LINE" | cut -f4)
TAG=$(echo    "$LINE" | cut -f5)

OUTFILE="${OUT}/${MATRIX}__${TREAT}__${PRIOR}__${TAG}.rds"

echo "== mkp-hp-ess-rep =="
echo "  job          : ${SLURM_JOB_ID:-local} (array task ${TIDX})"
echo "  host         : $(hostname)"
echo "  date         : $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  matrix       : ${MATRIX}  treatment: ${TREAT}  prior: ${PRIOR}"
echo "  seed         : ${SEED}    tag: ${TAG}"
echo "  out          : ${OUTFILE}"

cd "$SRC"

if [ -f "$OUTFILE" ]; then
  echo "Already exists, skipping: $OUTFILE"
  exit 0
fi

Rscript "${SRC}/dev/red-team/heavy-tests/casali-ess-hyperprior-vs-gamma.R" \
  --matrix          "${MATRIX}" \
  --treatment       "${TREAT}" \
  --prior           "${PRIOR}" \
  --seed            "${SEED}" \
  --partitions-dir  "${PARTS}" \
  --out             "${OUTFILE}" \
  --autopart-lib    "${LIB}" \
  --cores           "${SLURM_CPUS_PER_TASK:-1}"
