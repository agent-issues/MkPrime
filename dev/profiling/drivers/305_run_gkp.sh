#!/usr/bin/env bash
# #305 A5-06: gibbs_kPrime at its auto-pin (weight ~ nTrans) vs pinned at one
# sweep per nMoves-iteration cycle, paired by seed, default model.
# Usage: bash 305_run_gkp.sh <outdir> <pkg> <lib> [warmup] [maxTime]
set -u
OUT="$1"; PKG="$2"; LIB="$3"; WARMUP="${4:-3000}"; MAXT="${5:-75}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"
for ds in Vinther2008 Longrich2010 Sansom2010; do
  for seed in 1 2; do
    base="$OUT/${ds}_s${seed}"
    timeout $((MAXT + 150)) Rscript "$HERE/305_move_shares.R" "$ds" adapted \
      "$seed" "$WARMUP" "$MAXT" "${base}_auto.rds" "$PKG" "$LIB" 2>&1 |
      tail -1
    pin=$(Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]);
      cat(signif(1 / length(x$finalWeights), 6))' "${base}_auto.rds")
    timeout $((MAXT + 150)) Rscript "$HERE/305_move_shares.R" "$ds" gkp \
      "$seed" "$WARMUP" "$MAXT" "${base}_cycle.rds" "$PKG" "$LIB" "$pin" \
      2>&1 | tail -1
  done
done
echo ALL DONE
