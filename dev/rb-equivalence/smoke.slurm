#!/bin/bash
# Smoke / per-cell driver: run MkPrime then RB on the same Hamilton node,
# write rds outputs to /nobackup/pjjg18/rb-equivalence/out/, copy back.
#
# Submit:
#   sbatch --export=ALL,PID=950,MODEL=by_nt_9v,TARGET_ESS=128,TARGET_RHAT=1.025,MAX_TIME=1800 \
#          dev/rb-equivalence/smoke.slurm
#
# Or via the dispatcher helper run_cell.ps1 (see plan).
#
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --gres=tmp:4G
#SBATCH --job-name=rbeq
#SBATCH --output=/nobackup/pjjg18/rb-equivalence/logs/%x_%A_%j.out
#SBATCH --error=/nobackup/pjjg18/rb-equivalence/logs/%x_%A_%j.err
#SBATCH -p shared
#SBATCH --export=ALL

set -e

if [ -z "${PID:-}" ] || [ -z "${MODEL:-}" ]; then
  echo "ERROR: PID and MODEL env vars required" >&2
  exit 1
fi
TARGET_ESS=${TARGET_ESS:-128}
TARGET_RHAT=${TARGET_RHAT:-1.02}    # intra-sampler convergence target
MAX_TIME=${MAX_TIME:-3600}

export REPO=/nobackup/pjjg18/rb-equivalence
export NEOTRANS=/nobackup/pjjg18/neotrans
OUT=$REPO/out
CELL=$REPO/cells/${PID}_${MODEL}
mkdir -p $OUT $CELL $REPO/logs

echo "[smoke] cell=${PID}_${MODEL}  target=(rhat<${TARGET_RHAT}, ess>${TARGET_ESS})  max_time=${MAX_TIME}s  host=$(hostname)"

module load r/4.5.1 gcc/14.2 boost/1.88.0
export R_LIBS_USER=$REPO/lib

# --- MkPrime ---------------------------------------------------------------

echo "[smoke] === MkPrime ==="
t0=$SECONDS
Rscript $REPO/run_mkprime.R $PID $MODEL \
  --rhat=$TARGET_RHAT --ess=$TARGET_ESS --max-time=$MAX_TIME \
  --out-dir=$OUT
echo "[smoke] MkPrime wall: $((SECONDS - t0))s"

# --- RB --------------------------------------------------------------------

echo "[smoke] === RB ==="
Rscript $REPO/render_rev.R $PID $MODEL \
  --ess=$TARGET_ESS --max-time=$MAX_TIME \
  --neo-path=$NEOTRANS/inst/projects/project${PID}.neo.nex \
  --trans-path=$NEOTRANS/inst/projects/project${PID}.trans.nex \
  --out-dir=$CELL --cellinfo-dir=$OUT

RB_BIN=${RB_BIN:-/nobackup/pjjg18/revbayes/projects/cmake/build-pr816/rb}
cd $CELL
t0=$SECONDS
$RB_BIN long_${MODEL}.Rev
rb_wall=$((SECONDS - t0))
echo "[smoke] RB wall: ${rb_wall}s"

Rscript $REPO/post_rb.R $PID $MODEL --wall=$rb_wall --out-dir=$OUT --rb-dir=$CELL \
  --cellinfo-dir=$OUT --rb-bin=$RB_BIN

echo "[smoke] Cell done: ${PID}_${MODEL}"
