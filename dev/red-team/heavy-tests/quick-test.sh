#!/bin/bash
set -e
module load r/4.5.1 >/dev/null 2>&1
module load gcc/14.2 >/dev/null 2>&1 || true
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export R_LIBS=/nobackup/pjjg18/mkp-study/red-team/lib:/nobackup/pjjg18/mkp-study/lib
cd /nobackup/pjjg18/mkp-study/red-team/mkp-source
echo "=== quick-mode SBC: MkNT_geometric ==="
time Rscript /nobackup/pjjg18/mkp-study/red-team/heavy-tests/sbc.R \
  --quick --arm MkNT_geometric \
  --out /nobackup/pjjg18/mkp-study/red-team/results/sbc-quicktest \
  --seed 1
