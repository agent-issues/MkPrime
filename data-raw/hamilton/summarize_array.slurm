#!/bin/bash
#SBATCH --job-name=mkp-summarize
#SBATCH --partition=shared
#SBATCH --time=00:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/summarize_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/summarize_%A_%a.err

# Usage:  sbatch --array=0-259 summarize_array.slurm mk
#         sbatch --array=0-25  summarize_array.slurm mkp_geo
#
# tree_idx = TASK_ID / 10 + 1,  rep_idx = TASK_ID % 10 + 1 (mk arm)
# For mkp_geo the user has only run rep_01, so we map TASK_ID -> tree only,
# rep_idx = 1.

ARM="${1:?usage: sbatch --array=... summarize_array.slurm ARMNAME}"

# This is the only script here with no in-file `#SBATCH --array`; the range has
# to come from the command line. Submitted without one, SLURM_ARRAY_TASK_ID is
# unset, the arithmetic below reads the empty string as 0, and the job
# summarises t01_r01 alone -- then reports success, indistinguishable in sacct
# from a run that covered all 260 cells (#101).
: "${SLURM_ARRAY_TASK_ID:?must be submitted as an array job, e.g. sbatch --array=0-259}"

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

case "$ARM" in
  mk|mkp|mkp_eg|mk_kp1|mk_kp2|mk_k9|mk_k15|mk_k24|mk_k40|mk_ktrue|mkp_highk|mkp_logs|mk_tlshrink)
    TREE=$(( SLURM_ARRAY_TASK_ID / 10 + 1 ))
    REP=$(( SLURM_ARRAY_TASK_ID % 10 + 1 ))
    ;;
  mkp_geo)
    TREE=$(( SLURM_ARRAY_TASK_ID + 1 ))
    REP=1
    ;;
  *)
    echo "Unknown arm: $ARM" >&2; exit 2 ;;
esac

RESULTS_ROOT=/nobackup/pjjg18/mkp-study/results
DATA_ROOT=/nobackup/pjjg18/mkprime-files/tree-inference
OUT_DIR=/nobackup/pjjg18/mkp-study/summary

mkdir -p "$OUT_DIR"

echo "[summarize] arm=$ARM tree=$TREE rep=$REP task=$SLURM_ARRAY_TASK_ID"

Rscript /nobackup/pjjg18/mkp-study/summarize_streamed.R \
  "$TREE" "$REP" "$ARM" "$RESULTS_ROOT" "$DATA_ROOT" "$OUT_DIR"
