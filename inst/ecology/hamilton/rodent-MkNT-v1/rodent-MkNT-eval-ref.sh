#!/bin/bash
#SBATCH --job-name=rod-eval-ref
#SBATCH --output=/nobackup/pjjg18/mkp-rodent-eval-ref/eval.out
#SBATCH --error=/nobackup/pjjg18/mkp-rodent-eval-ref/eval.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1

# Evaluate aware logLik at two reference topologies (fixed-topology MCMC):
#   A) Blind MR consensus -- biologically plausible, converged blind run
#   B) v6 run-3 best tree -- highest logLik found so far (-5290.9),
#      but biologically nonsensical
#
# Answers: "would a good tree be accepted if only it could be found?"
#   |diff| < 10    -> mode-trap; blind topology roughly equivalent in logLik
#   blind << run3  -> aware model genuinely prefers nonsensical topology
#   blind >> run3  -> run3 is itself a trap; true peak not yet found

set -euo pipefail
module load r/4.5.1
module load gcc/14.2

export MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKPRIME_ECO_RESYNC_EVERY=200

export NEX_FILE=/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex
export BLIND_DIR=/nobackup/pjjg18/mkp-rodent-MkNT-v1/blind
export V6_DIR=/nobackup/pjjg18/mkp-rodent-MkNT-v6/aware
export OUT_DIR=/nobackup/pjjg18/mkp-rodent-eval-ref

mkdir -p "$OUT_DIR"

echo "=== Rodent aware logLik evaluation at reference topologies ==="
echo "Started: $(date)"
echo "Host:    $(hostname)"

Rscript /nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime/inst/ecology/hamilton/rodent-MkNT-v1/eval_aware_at_reference.R

echo "=== Finished: $(date) ==="
