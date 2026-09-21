#!/bin/bash
#SBATCH --job-name=mkp-geo
#SBATCH --partition=shared
#SBATCH --time=08:00:00
# run_one.R sets maxTime = 7.5h, below this 8h wall, so the R-level graceful
# stop fires and `saveRDS(partial, ...)` runs before SIGKILL. (This note used to
# say maxTime was 10h and pointed at dev/red-team/findings.md; HARNESS-001 was
# fixed long ago and findings.md is retired in favour of GitHub issues.)
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --array=0-25
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/mkp_geo_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/mkp_geo_%A_%a.err

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

# Pilot: 26 trees × 1 rep (task IDs 0–25, one rep per tree)
TREE=$(( SLURM_ARRAY_TASK_ID + 1 ))
REP=1

DATA_ROOT=/nobackup/pjjg18/mkprime-files/tree-inference
OUT_DIR=/nobackup/pjjg18/mkp-study/results

Rscript /nobackup/pjjg18/mkp-study/run_one.R \
  "$TREE" "$REP" "$DATA_ROOT" "$OUT_DIR" mkp_geo
