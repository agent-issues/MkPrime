#!/bin/bash
# sim3-v4c.sh — Sim 3 v4 stress regime (blind + aware, 1 rep).
# Target: blind P(eco12) > 0 (active false clade support); aware recovers truth.

#SBATCH --job-name=mkp-v4c
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v4c/v4c.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v4c/v4c.err
#SBATCH --time=03:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v4c/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-v4c/run_v4c.R
