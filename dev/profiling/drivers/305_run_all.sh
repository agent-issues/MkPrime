#!/usr/bin/env bash
# #305 A5-07 batch: adapted vs capped mh_logit_p, paired by seed, on three
# small TreeSearch matrices, with and without the ascertainment term.
# Usage: bash 305_run_all.sh <outdir> <pkg> <lib> [warmup] [maxTime]
set -u
OUT="$1"; PKG="$2"; LIB="$3"; WARMUP="${4:-3000}"; MAXT="${5:-75}"
HERE="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$OUT"
for coding in variable none; do
  for ds in Vinther2008 Longrich2010 Sansom2010; do
    for seed in 1 2; do
      base="$OUT/${ds}_${coding}_s${seed}"
      timeout $((MAXT + 150)) Rscript "$HERE/305_move_shares.R" "$ds" adapted \
        "$seed" "$WARMUP" "$MAXT" "${base}_adapted.rds" "$PKG" "$LIB" NA \
        "$coding" 2>&1 | tail -1
      # Pin at the share the scheduler starts from (first batch's weight).
      pin=$(Rscript -e 'x <- readRDS(commandArgs(TRUE)[1]);
        cat(signif(x$batches[[1]]$weights[["mh_logit_p"]], 6))' \
        "${base}_adapted.rds")
      timeout $((MAXT + 150)) Rscript "$HERE/305_move_shares.R" "$ds" capped \
        "$seed" "$WARMUP" "$MAXT" "${base}_capped.rds" "$PKG" "$LIB" "$pin" \
        "$coding" 2>&1 | tail -1
    done
  done
done
echo ALL DONE
