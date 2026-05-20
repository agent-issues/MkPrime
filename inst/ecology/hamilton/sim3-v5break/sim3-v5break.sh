#!/bin/bash
# sim3-v5break.sh -- Sim 3 v5-break (BREAK-BLIND) regime, PT 3-rep array.
#
# Design intent: actually break the ecology-blind chain at 16 tips so the
# ecology-aware chain has something to rescue. Across 7+ prior v4 / v4-cross
# regimes the correctly-scored result has been that blind essentially
# recovers truth; v5break shrinks ancestry signal (stemBrClade=0.03,
# 80 chars) and keeps the eco confound in the unsaturated band
# (stemBrEco=0.10, phi=6, pi0=0.45).
#
# 3-rep array hedges against single-rep stochasticity (which we have
# observed in v4cross-b-pt: rep 1 collapsed differently from rep 5).
#
# Walltime: PT 4-chain x 100k iter ~ 90 min per chain on v4cross-b-pt5;
# blind + aware sequential -> ~3-4 h. 4 h walltime is safe.

#SBATCH --job-name=mkp-v5break
#SBATCH --array=1-3
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-v5break/v5break-%a.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-v5break/v5break-%a.err
#SBATCH --time=04:00:00
#SBATCH --mem=512M
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-v5break/results/rep$(printf "%02d" "$SLURM_ARRAY_TASK_ID")
cd "$MKP_REPO_ROOT"

Rscript inst/ecology/hamilton/sim3-v5break/run_v5break.R
