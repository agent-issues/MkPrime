#!/usr/bin/env bash
# Full grid: 8 datasets x {gibbsSpr on, off} x 3 seeds, 200k iterations each.
SP="${GSPR_DIR:?set GSPR_DIR to a scratch dir containing common.R, align.R and the numbered scripts}"
OUT="$SP/grid"
LOG="$SP/grid.log"
NITER=200000
PAR=2   # 2-core policy
mkdir -p "$OUT"
: > "$LOG"

jobs_list=()
for t in 1 2 3 4 5 6 7 8; do
  for arm in on off; do
    for s in 1 2 3; do
      jobs_list+=("$t $arm $s")
    done
  done
done

i=0
for j in "${jobs_list[@]}"; do
  set -- $j
  # outermost guard: above the 3600 s setTimeLimit in 02 (maxTime 1500 s x 2.1 overshoot)
  timeout 3900 Rscript "$SP/02-toggle-run-cell.R" "$1" "$2" "$3" "$NITER" "$OUT" \
    >> "$LOG" 2>&1 &
  i=$((i+1))
  if [ $((i % PAR)) -eq 0 ]; then wait; echo "--- batch $((i/PAR)) done ---" >> "$LOG"; fi
done
wait
echo "ALL DONE" >> "$LOG"
