#!/bin/bash
# sim3-phi2-fixedpriors.sh — phi=2 rerun with tightened priors (R5-1/R5-2/R5-3).
#
# Rerun of job 17185756 with explicit prior hyperparameters to break the pi0
# feedback loop identified in red-team audit rounds R5-1, R5-2, R5-3:
#   sigmaPhi=1.5, thetaAlpha=2, thetaBeta=2, rho0Alpha=360, rho0Beta=120.
# No package reinstall required.

#SBATCH --job-name=mkp-phi2-fp
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-phi2/phi2-fp.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-phi2/phi2-fp.err
#SBATCH --time=03:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-phi2/results-fixedpriors
cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/sim3-phi2/run_phi2_fixedpriors.R
