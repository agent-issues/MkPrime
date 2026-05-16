#!/bin/bash
# sim3-v4b.sh — Sim 3 v4 balanced sweet-spot test (blind + aware, 1 rep).
# Target: blind P(eco12) high; aware recovers true clade A, C, AC.

#SBATCH --job-name=mkp-v4b
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v4b/v4b.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v4b/v4b.err
#SBATCH --time=03:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v4b/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-v4b/run_v4b.R
