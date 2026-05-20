#!/bin/bash
# sim3-v4cross-b-pt.sh — Sim 3 v4-cross balanced regime, PT rerun (nChains=4).
# PT chains are sequential in MkPrimeMCMC; cpus-per-task=1 is correct.
# Single-chain aware was 22 min → 4 sequential PT chains ≈ 90 min; 6h is safe.

#SBATCH --job-name=mkp-v4xb-pt
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v4cross-b-pt/v4xb-pt.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v4cross-b-pt/v4xb-pt.err
#SBATCH --time=06:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v4cross-b-pt/results
cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/sim3-v4cross-b-pt/run_v4cross_b_pt.R
