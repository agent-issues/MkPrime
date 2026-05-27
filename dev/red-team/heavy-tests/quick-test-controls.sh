#!/bin/bash
# Quick-mode SBC on the two MkNT control arms. Per advisor decision rule:
# if either control still fails at p < 0.05, harness has a third bug —
# stop, document, fall back to L6 analytic EG-001 proof. Do not loop.
set -e
module load r/4.5.1 >/dev/null 2>&1
module load gcc/14.2 >/dev/null 2>&1 || true
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS=/nobackup/pjjg18/mkp-study/red-team/lib:/nobackup/pjjg18/mkp-study/lib
cd /nobackup/pjjg18/mkp-study/red-team/mkp-source

for arm in MkNT_geometric MkNT_logseries; do
  echo "=== quick: $arm ==="
  time Rscript /nobackup/pjjg18/mkp-study/red-team/heavy-tests/sbc.R \
    --quick --arm $arm \
    --out /nobackup/pjjg18/mkp-study/red-team/results/sbc-quick2 \
    --seed 42
  echo
done
