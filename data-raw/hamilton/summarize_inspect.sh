#!/bin/bash
# Quick inspection helper for streamed-MCMC logs.
D=/nobackup/pjjg18/mkp-study/results/t01_r01
echo "===MK header==="
head -1 "$D/mk_run_1.log" | tr '\t' '\n' | nl
echo
echo "===MK ncols==="
head -1 "$D/mk_run_1.log" | awk -F'\t' '{print NF}'
echo
echo "===MKP_GEO header==="
head -1 "$D/mkp_geo_run_1.log" | tr '\t' '\n' | nl | head -20
echo "(...)"
head -1 "$D/mkp_geo_run_1.log" | tr '\t' '\n' | nl | tail -55 | head -10
echo
echo "===MKP_GEO ncols==="
head -1 "$D/mkp_geo_run_1.log" | awk -F'\t' '{print NF}'
echo
echo "===MK row 2 (comment)==="
sed -n '2p' "$D/mk_run_1.log" | head -c 200; echo
echo "===MK row 3 (first data)==="
sed -n '3p' "$D/mk_run_1.log" | tr '\t' '\n' | nl | head -10
echo
echo "===n lines (mk_run_1, trees)==="
wc -l "$D/mk_run_1.log" "$D/mk_trees.nwk" "$D/mkp_geo_run_1.log" "$D/mkp_geo_trees.nwk" 2>/dev/null
echo
echo "===n comment rows mk==="
head -5 "$D/mk_run_1.log" | nl | head -5 | sed 's/^/  /'
echo "===n comment rows geo==="
head -5 "$D/mkp_geo_run_1.log" | nl | head -5 | sed 's/^/  /'
echo "===dep check==="
module load r/4.5.1 gcc/14.2 >/dev/null 2>&1
Rscript -e '.libPaths(c("/nobackup/pjjg18/mkp-study/lib", .libPaths())); cat("data.table:", requireNamespace("data.table", quietly=TRUE), "\n"); cat("TreeDist:", requireNamespace("TreeDist", quietly=TRUE), "\n"); cat("ape:", requireNamespace("ape", quietly=TRUE), "\n")'
