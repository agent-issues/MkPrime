# Plan: M-112 Reassessment + M-113 End-to-End Benchmark

## Context

After completing M-111 (partial CL for Gibbs swap), the optimization
roadmap has three open tasks:
- M-112: Cache CLGroup across consecutive Gibbs invocations
- M-113: End-to-end mixing benchmark on production datasets
- M-114: Partial CL under Q-heterogeneity (low priority)

## Analysis: M-112 is not viable as originally described

The MCMC engine draws **one random move per chain per iteration** (lines
2428–2467 of mcmc.cpp). Gibbs SPR and swap are each ~13% of the move
weight for Sun2018. Between any two Gibbs calls, ~6–8 other moves
typically execute, and most of them modify state that the CLGroup
depends on:

| Move | Changes what? | Invalidates CLGroup? |
|------|--------------|---------------------|
| tree_length scale | treeLength → all absLen | Yes |
| beta_simplex | relBrLengths → all absLen | Yes |
| NNI/SPR/TBR | topology | Yes |
| kPrime int_walk | kPrime → group structure | Yes |
| rateLoss/rateNeo scale | rateLoss/rateNeo | Yes |
| rateLogSd scale | ACRV rates | Yes |
| gibbs_p | p only | **No** (p not in CLGroup) |

Estimated cache hit rate: **<5%**, making cross-invocation caching
effectively worthless.

**Alternative considered**: Combine Gibbs SPR + swap into one "compound
topology" move that shares the caching_downpass. Analysis shows this is
actually **worse** for throughput: the combined move (13.5ms) runs at
26% frequency (= 3.5ms expected/iter) vs separate moves at 13% each
(= 2.3ms expected/iter). Combining forces evaluation of both candidate
sets every time, whereas separate draws only evaluate one set per
iteration.

**Conclusion**: Close M-112 as WONTFIX. The per-invocation
caching_downpass (one full O(N×C) traversal) is the irreducible minimum
cost of the partial CL approach. No cross-invocation caching strategy
can help given the move schedule.

## Plan: M-113 + profiling

Since M-112 doesn't pan out, the next useful work is empirical:
understand where time actually goes and whether Gibbs moves improve
mixing efficiency (ESS/s).

### Step 1: Per-move timing breakdown (profiling pass)

Run a short MCMC on Sun2018 (54 tips) with Gibbs moves enabled.
Extract the per-move timing from `chain_time_ns` (already collected by
the C++ engine) to see:
- What % of wall time goes to each move type
- Mean cost per invocation for each move
- Whether any non-Gibbs move is unexpectedly expensive

This uses `RunMkPrime()` directly — no VTune needed.

### Step 2: ESS/s benchmark (M-113)

Run longer chains (≥50k iterations, ideally 100k) on Sun2018 with
three configurations:
1. **Baseline**: NNI + SPR only (no Gibbs)
2. **+Gibbs**: add gibbs_spr + gibbs_subtree_swap
3. **+Gibbs +Block**: add block_gibbs_branch

Measure:
- ESS for log_posterior, tree_length (via `coda::effectiveSize`)
- Wall time
- ESS/s for each configuration

Compare to decide whether Gibbs moves actually improve mixing
efficiency on this dataset. Previous M-091 results were debug build on
a 23-taxon dataset — this is the first production benchmark on a
realistic problem size.

### Step 3: Update coordination

- Close M-112 as WONTFIX in to-do.md and completed-tasks.md
- Record M-113 benchmark results
- Update S-PROF standing task with profiling findings
- File any new tasks if profiling reveals unexpected hotspots

## Execution

All work on the `main` branch (no feature branch needed for benchmarking).
