#!/bin/bash
# sim3-pt-test.sh — single-rep PT (nChains=4) test job
#
# Re-uses the already-bootstrapped /nobackup/$USER/mkp-sim3-multirep-v3
# project (R library + MkPrime install). Drops a dedicated output dir
# under /nobackup/$USER/mkp-sim3-pt-test/results so the multirep
# checkpoints aren't touched.

#SBATCH --job-name=mkp-pt-test
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-pt-test/pt-test.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-pt-test/pt-test.err
#SBATCH --time=04:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

# Re-use the existing multirep-v3 install.
export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-pt-test/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-pt-test/run_pt_test.R
