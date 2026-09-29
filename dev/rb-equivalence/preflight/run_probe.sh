#!/bin/bash
# Evaluate every probe_*.Rev in the current directory with RevBayes; seconds
# per cell, so a login node is fine.  Writes out_<cell>.txt.
module load gcc/14.2 boost/1.88.0 2>/dev/null
RB=${RB_BIN:-/nobackup/$USER/revbayes/projects/cmake/build-pr816/rb}
for f in probe_*.Rev; do
  cell=${f#probe_}; cell=${cell%.Rev}
  timeout 300 "$RB" "$f" > "out_$cell.txt" 2>&1
  echo "$cell exit=$? rows=$(grep -c '^ROW' "out_$cell.txt")"
done
