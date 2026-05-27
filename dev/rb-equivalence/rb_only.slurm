#!/bin/bash
# RB-only cell driver: skips MkPrime, runs only the RB side.
# Use when MkPrime rds is already in $REPO/out/ and we want to validate RB
# independently or rerun RB only.
#
# Submit:
#   sbatch --export=ALL,PID=635,MODEL=by_nt_9v,TARGET_ESS=128,MAX_TIME=1800 \
#          dev/rb-equivalence/rb_only.slurm
#
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH --time=2:00:00
#SBATCH --gres=tmp:4G
#SBATCH --job-name=rbeqrb
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
MAX_TIME=${MAX_TIME:-3600}

export REPO=/nobackup/pjjg18/rb-equivalence
export NEOTRANS=/nobackup/pjjg18/neotrans
OUT=$REPO/out
CELL=$REPO/cells/${PID}_${MODEL}
mkdir -p $OUT $CELL $REPO/logs

echo "[rb-only] cell=${PID}_${MODEL} target=ess>${TARGET_ESS} max_time=${MAX_TIME}s host=$(hostname)"

module load r/4.5.1 gcc/14.2 boost/1.88.0
export R_LIBS_USER=$REPO/lib

# Render RB cell
Rscript $REPO/render_rev.R $PID $MODEL \
  --ess=$TARGET_ESS --max-time=$MAX_TIME \
  --neo-path=$NEOTRANS/inst/projects/project${PID}.neo.nex \
  --trans-path=$NEOTRANS/inst/projects/project${PID}.trans.nex \
  --out-dir=$CELL

cd $CELL
t0=$SECONDS
set +e
${RB_BIN:-/nobackup/pjjg18/revbayes/projects/cmake/build-pr816/rb} long_${MODEL}.Rev
rb_exit=$?
set -e
rb_wall=$((SECONDS - t0))
echo "[rb-only] RB exit=$rb_exit wall=${rb_wall}s"

if [ $rb_exit -ne 0 ]; then
  echo "[rb-only] RB FAILED -- skipping post_rb"
  exit $rb_exit
fi

Rscript $REPO/post_rb.R $PID $MODEL --wall=$rb_wall --out-dir=$OUT --rb-dir=$CELL
echo "[rb-only] cell done: ${PID}_${MODEL}"
