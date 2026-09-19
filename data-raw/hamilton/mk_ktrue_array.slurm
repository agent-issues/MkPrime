#!/bin/bash
#SBATCH --job-name=mkp-mk-ktrue
#SBATCH --partition=shared
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --array=0-259
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/mk_ktrue_%A_%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/mk_ktrue_%A_%a.err

# RECONSTRUCTED 2026-09-19 from its sibling array scripts (#100).
#
# This file was committed at 0 bytes in 714d3bc and has been 0 bytes at every
# commit since, so `sbatch` on it fails with "Batch script is empty" and the
# arm has never been runnable from the repository. The 260 tasks whose result
# 714d3bc reports -- `mk_ktrue 0.2624 (ORACLE LOSES to mk_k40 by 0.024,
# p=5.9e-45)` -- were submitted from a script that is not in version control.
#
# This reconstruction makes the arm runnable again. It does NOT recover the
# configuration that produced the recorded number: the original's --time,
# --mem and --array range are unknown, so the recorded mk_ktrue CID remains
# unauditable and any comparison against it should say so.

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

TREE=$(( SLURM_ARRAY_TASK_ID / 10 + 1 ))
REP=$(( SLURM_ARRAY_TASK_ID % 10 + 1 ))

DATA_ROOT=/nobackup/pjjg18/mkprime-files/tree-inference
OUT_DIR=/nobackup/pjjg18/mkp-study/results

Rscript /nobackup/pjjg18/mkp-study/run_one.R \
  "$TREE" "$REP" "$DATA_ROOT" "$OUT_DIR" mk_ktrue
