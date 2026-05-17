#!/bin/bash
# sim3-v4cross-c.sh — Sim 3 v4-cross STRESS regime (blind + aware, 1 rep).
# Uses SIM3V4_PARAMS$v4c: stemBrEco=0.15, stemBrClade=0.05 — max eco signal.
# Target: blind P(falseInner) or P(falseAB) > 0 (active false eco clade);
# aware recovers true (A,C),(B,D) topology.

#SBATCH --job-name=mkp-v4xc
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v4cross-c/v4xc.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v4cross-c/v4xc.err
#SBATCH --time=03:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v4cross-c/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-v4cross-c/run_v4cross_c.R
