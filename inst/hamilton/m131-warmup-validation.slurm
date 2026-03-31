#!/bin/bash
#SBATCH --job-name=m131-warmup
#SBATCH --partition=shared
#SBATCH --array=0-31
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/m131_%A_%a.out

# M-131: Warmup stabilisation validation study
# 8 datasets x 4 seeds = 32 tasks

DATASETS=(Longrich2010 Vinther2008 DeAssis2011 Wortley2006
          Sun2018 Wilson2003 Zhu2013 Dikow2009)
SEEDS=(2847 5193 7621 9034)

IDX=$SLURM_ARRAY_TASK_ID
DS_IDX=$((IDX / 4))
SEED_IDX=$((IDX % 4))

export DATASET=${DATASETS[$DS_IDX]}
export SEED=${SEEDS[$SEED_IDX]}
export OUTDIR=results

echo "Task $IDX: dataset=$DATASET seed=$SEED"

module load r/4.5.1
module load gcc/14.2

BASEDIR=/nobackup/$USER/m131
export R_LIBS_USER=$BASEDIR/lib

cd $BASEDIR
mkdir -p logs results

Rscript m131-warmup-validation.R

echo "Done: exit $?"
