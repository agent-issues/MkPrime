#!/bin/bash
# sim3-v4a.sh — Sim 3 v4 conservative smoke test (blind + aware, 1 rep).
# Expect: both chains recover true topology; verifies v4 harness.

#SBATCH --job-name=mkp-v4a
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v4a/v4a.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v4a/v4a.err
#SBATCH --time=03:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v4a/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-v4a/run_v4a.R
