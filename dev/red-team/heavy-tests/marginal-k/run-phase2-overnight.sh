#!/bin/bash
# Phase-2 overnight confirmations (sequential -> one .so lock at a time).
#  1. MOVES-ON OVERLAP: the four re-enabled weighted/block moves ON in BOTH modes;
#     confirms they do not bias the marginal_k posterior vs sampled_k (behavioral
#     gate on top of the deterministic gap-sweep + Test 3/4).
#  2. p-RECHECK: re-run the 16x48 row (incl. the d_sd=3.18 FAIL cell) at 200k iters
#     to push marginal-p-ESS >> floor; if d_sd resolves it was a low-ESS sd-MCSE
#     artifact, if it persists it is a real mh_logit_p spread bias to investigate.
cd /c/Users/pjjg18/GitHub/worktrees/mkp/marginal-k || exit 1
H=dev/red-team/heavy-tests/marginal-k
S=$H/T-OVL-sampled-vs-marginal.R

echo "[$(date)] === (1) MOVES-ON OVERLAP (weightedSpr+SubtreeSwap+BranchScale+blockGibbsBranch) ==="
MARGINAL_K_OVL_EXTRA=weightedSpr,weightedSubtreeSwap,weightedBranchScale,blockGibbsBranch \
  MARGINAL_K_OVL_NTIP=8,16 MARGINAL_K_OVL_NCHAR=24 MARGINAL_K_OVL_REP=2 \
  MARGINAL_K_OVL_ITER=60000 MARGINAL_K_OVL_WARM=8000 \
  Rscript "$S" > "$H/T-OVL-moveson.log" 2>&1
echo "[$(date)] moves-on exit=$?"

echo "[$(date)] === (2) p-RECHECK (gated, 16x48 row, 200k iters) ==="
MARGINAL_K_OVL_NTIP=16 MARGINAL_K_OVL_NCHAR=48 MARGINAL_K_OVL_REP=3 \
  MARGINAL_K_OVL_ITER=200000 MARGINAL_K_OVL_WARM=15000 \
  Rscript "$S" > "$H/T-OVL-precheck-200k.log" 2>&1
echo "[$(date)] p-recheck exit=$?"
echo "[$(date)] === OVERNIGHT DONE ==="
