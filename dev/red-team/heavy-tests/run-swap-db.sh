#!/bin/bash
# Run subtree-swap-db.R full-scale on the Hamilton login node.
# Per tree-ess-and-swap.md, n=5 1M iter is 5-10 min on a compute node;
# login node is plenty. The n=6/n=7 escalation needs the harness modified.
set -e
module load r/4.5.1 >/dev/null 2>&1
module load gcc/14.2 >/dev/null 2>&1 || true
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS=/nobackup/pjjg18/mkp-study/red-team/lib:/nobackup/pjjg18/mkp-study/lib
cd /nobackup/pjjg18/mkp-study/red-team/mkp-source
LOG=/nobackup/pjjg18/mkp-study/red-team/logs/swap-db_$(date +%Y%m%d_%H%M%S).log
{
  echo "=== subtree-swap-db full ==="
  date
  time Rscript /nobackup/pjjg18/mkp-study/red-team/heavy-tests/subtree-swap-db.R 2>&1
  echo "=== DONE ==="
  date
} > "$LOG" 2>&1
echo "Log: $LOG"
