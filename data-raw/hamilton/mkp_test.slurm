#!/bin/bash
#SBATCH --job-name=mkp-test
#SBATCH --partition=shared
#SBATCH --time=0:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --output=/nobackup/pjjg18/mkp-study/logs/mkp_test_%j.out
#SBATCH --error=/nobackup/pjjg18/mkp-study/logs/mkp_test_%j.err

module load r/4.5.1 gcc/14.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

DATA_ROOT=/nobackup/pjjg18/mkprime-files/tree-inference
OUT_DIR=/nobackup/pjjg18/mkp-study/results

Rscript /nobackup/pjjg18/mkp-study/run_one.R \
  1 1 "$DATA_ROOT" "$OUT_DIR" --test
