#!/bin/bash
# sim3-phi2.sh — softened-landscape (phi=2) blind+aware single-chain test.
#
# Posture 2 probe: does standard MCMC recover AC at half the ecology magnitude?

#SBATCH --job-name=mkp-phi2
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-phi2/phi2.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-phi2/phi2.err
#SBATCH --time=03:00:00
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-phi2/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-phi2/run_phi2.R
