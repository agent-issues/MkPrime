#!/bin/bash
# sim3-truth-start.sh — chain started at truth tree (rep 05 seed)
#
# Diagnostic: does the aware chain stay at AC, or drift away?
# Answers whether AC is genuinely the posterior mode (fix = eco proposals)
# or merely high-logPost at MAP nuisance (fix = something else).

#SBATCH --job-name=mkp-truth-start
#SBATCH --output=/nobackup/pjjg18/mkp-sim3-truth-start/truth-start.out
#SBATCH --error=/nobackup/pjjg18/mkp-sim3-truth-start/truth-start.err
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1

set -euo pipefail

module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export MKP_REPO_ROOT=/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime

mkdir -p /nobackup/pjjg18/mkp-sim3-truth-start/results
cd "$MKP_REPO_ROOT"

Rscript inst/hamilton/sim3-pt-test/run_truth_start.R
