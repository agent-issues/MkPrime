#!/bin/bash
#SBATCH --job-name=mkp-validate
#SBATCH --partition=shared
#SBATCH --time=8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-259
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/mkp_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/mkp_%A_%a.err

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Map linear index (0-259) to tree (1-26) and rep (1-10)
TREE=$(( SLURM_ARRAY_TASK_ID / 10 + 1 ))
REP=$(( SLURM_ARRAY_TASK_ID % 10 + 1 ))

DATA_ROOT=/nobackup/pjjg18/mkprime-files/tree-inference
OUT_DIR=/nobackup/pjjg18/mkp-study/results

Rscript /nobackup/pjjg18/mkp-study/run_one.R \
  "$TREE" "$REP" "$DATA_ROOT" "$OUT_DIR"
