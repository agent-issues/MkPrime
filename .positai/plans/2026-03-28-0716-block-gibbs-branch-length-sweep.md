# Plan: Block Gibbs Branch-Length Sweep + Dimension-Adjusted Scheduler

**Date:** 2026-03-28
**Scope:** New move type + scheduler enhancement (Phase 7d / Phase 9d)
**Worktree:** `mkp-gibbs` (builds on M-087 infrastructure)

---

## Motivation

BetaSimplex and WeightedBranchLengthScale each update a single edge pair
per call. For trees with many edges, traversing the full branch-length
space requires many iterations of single-site updates. A block sweep
that updates all edge pairs in one MCMC "move" can make faster progress
through branch-length space at the cost of more likelihood evaluations
per call.

HMC was the original candidate (M-054) but is a poor fit during active
topology exploration: its gradient-based trajectory is invalidated by
subsequent topology changes, and the implementation cost (analytical
gradients, log-ratio simplex transform, leapfrog tuning / NUTS) is
substantial. A block Gibbs sweep reuses the existing bin-based
approximate-conditional machinery from `weighted_branch_scale_impl`
(M-087), avoids gradient computation, and is robust to topology changes
because each edge's conditional is recomputed from scratch.

The adaptive scheduler (M-092) currently scores moves as
`accept_rate / cost`. This systematically undervalues multi-dimensional
moves: a block sweep that moves 37 edge pairs counts as "1 acceptance"
just like a BetaSimplex that moves 1 pair. Adding a `dim` field to the
scheduler corrects this for any current or future block proposal.

---

## Design

### Part 1: Dimension-adjusted scheduler score

**Change to M-092's score formula:**

```
score[m] = accept_rate[m] × dim[m] / mean_cost_s[m]
```

where `dim[m]` is the number of parameters the move updates per call.
All existing moves have `dim = 1` except the new block sweep
(`dim = nEdge`).

**Implementation:**

- R side: add a `dim` field to each move spec returned by `.BuildMoves()`.
  All current moves get `dim = 1L`.
- R side: `.AdaptMoveWeights()` gains a `moveDim` parameter. The score
  computation on line 2040 becomes:
  ```r
  logScores <- log(pmax(acceptRate, 1e-12)) + log(moveDim[scoreableIdx])
               - log(meanCostS)
  ```
- C++ side: no change needed — `dim` is used only in the R-side
  scheduler, not in the C++ batch loop.

This is a backward-compatible change: all existing moves get `dim = 1`,
so `log(1) = 0` and the score is unchanged.

### Part 2: Block Gibbs branch-length sweep

**New move type 15: `block_gibbs_branch`**

A single call sweeps over all edge pairs in random-permutation order.
For each pair, it performs the same bin-based approximate-conditional
sampling as `weighted_branch_scale_impl` (M-087) and independently
accepts/rejects via MH. This is a **random-permutation-scan
MH-within-Gibbs** composition, which preserves the target distribution
because each component step leaves the target invariant.

#### Algorithm

```
block_gibbs_branch_sweep_impl(data, state, beta, nBins):
  nEdge = state->relBrLengths.size()
  bins  = get_branch_bins(nBins)       // reuse M-087 bin structure

  // Random permutation of edge indices
  perm = {0, 1, ..., nEdge-1}
  Fisher-Yates shuffle(perm)

  nAccepted = 0
  absLen = treeLength × relBrLengths   // working copy

  for idx in perm:
    // Pick partner: the "other" edge for this pair
    other = random edge ≠ idx
    oldRelA = relBrLengths[idx]
    oldRelB = relBrLengths[other]
    relTotal = oldRelA + oldRelB
    if relTotal ≤ 0: continue
    oldF = oldRelA / relTotal
    absTotal = absLen[idx] + absLen[other]

    // Evaluate LL at each bin midpoint (nBins evals)
    for b in 0..nBins-1:
      trialAbs = absLen (modified at [idx] and [other] only)
      trialAbs[idx]   = mids[b] × absTotal
      trialAbs[other] = (1 - mids[b]) × absTotal
      midLL[b] = compute_full_loglik_at(data, state, parent, child, trialAbs)

    // Weight, sample bin, draw fraction (identical to M-087)
    weights = exp(beta × (midLL - max(midLL)))
    chosenBin = categorical_sample(weights)
    newF ~ Beta(mids[chosenBin] × conc + 1, (1-mids[chosenBin]) × conc + 1)

    // Hastings ratio (same as M-087)
    oldBin = bin containing oldF
    logHR = log(w[oldBin]) + logBeta(oldF|oldBin)
          - log(w[chosenBin]) - logBeta(newF|chosenBin)

    // Evaluate LL + prior at proposed fraction
    trialAbs[idx]   = newF × absTotal
    trialAbs[other] = (1-newF) × absTotal
    proposedLL = compute_full_loglik_at(...)
    proposedPrior = cpp_log_prior(... with updated relBr ...)

    // MH accept/reject for this pair
    logAlpha = beta × (proposedLL - currentLL) + (proposedPrior - currentPrior) + logHR
    if log(U) < logAlpha:
      relBrLengths[idx]   = newF × relTotal
      relBrLengths[other] = (1-newF) × relTotal
      absLen[idx]   = trialAbs[idx]
      absLen[other] = trialAbs[other]
      state->logLik   = proposedLL
      state->logPrior = proposedPrior
      nAccepted++

  return nAccepted > 0
```

**Cost:** `nEdge × (nBins + 1)` full likelihood evaluations per sweep.
The `+1` is for the final proposed-state evaluation at each edge.
For a 20-tip tree (37 edges, 10 bins): ~407 evaluations.
For a 50-tip tree (97 edges, 10 bins): ~1067 evaluations.

#### Key differences from weighted_branch_scale_impl

| Aspect | WeightedBranchLengthScale | Block Gibbs Sweep |
|--------|---------------------------|-------------------|
| Edges per call | 1 pair | All pairs |
| Scan order | Random single edge | Random permutation |
| Accept/reject | Single MH for the pair | Independent MH per pair |
| State mutation | Deferred to do_move_impl | In-place within sweep |
| Cost | O(B) | O(nEdge × B) |
| dim | 1 | nEdge |

#### Integration into do_move_impl

The block sweep handles its own accept/reject loop internally (like
gibbs_spr_impl and weighted_spr_impl). The dispatch case returns early:

```cpp
case 15: { // block_gibbs_branch — M-054 reframed
  return block_gibbs_branch_sweep_impl(data, state, beta, data->nBranchBins);
}
```

No rollback needed in the outer function — the sweep manages its own
state mutations, accepting/rejecting each pair independently.

#### R-side wiring

- `.kMoveTypes`: add `block_gibbs_branch = 15L`
- `.BuildMoves()`: add entry when `mcmc$blockGibbsBranch == TRUE`,
  with `weight = max(1, nEdge / 4)`, `dim = nEdge`
- `MkPrimeMCMC()`: add `blockGibbsBranch` toggle (default `FALSE` —
  opt-in, pending M-091 mixing validation)
- `.AdaptTuning()`: add `block_gibbs_branch = NA_real_` (no tuning
  parameter to adapt — the move is tuning-free like other Gibbs moves)

#### Parallel tempering

The sweep passes `beta` (the chain's inverse temperature) into the bin
weights: `exp(beta × (midLL - maxLL))`. The MH test uses
`beta × (proposedLL - currentLL)` for the likelihood component. Prior
is unheated. This is consistent with all other tempered moves.

---

## Worktree decision

**No new worktree.** This work goes directly on `mkp-gibbs`
(`feature/gibbs-weighted-moves` branch) because:

1. It depends on M-087's bin infrastructure (`BranchBins`,
   `get_branch_bins`, `compute_full_loglik_at`)
2. It's a natural extension of Phase 9's move suite
3. The scheduler dim enhancement (Part 1) modifies `.AdaptMoveWeights()`
   which is already on this branch
4. The new move will be validated alongside the Phase 9 moves in M-091

---

## Task breakdown

### Task 1: Dimension-adjusted scheduler score

**Files:** `R/RunMkPrime.R` (`.BuildMoves()`, `.AdaptMoveWeights()`),
`R/MkPrimeMCMC.R`

1. Add `dim = 1L` field to every move spec in `.BuildMoves()`
2. Add `moveDim` parameter to `.AdaptMoveWeights()`; incorporate into
   log-score: `log(acceptRate) + log(dim) - log(cost)`
3. Extract `moveDim` alongside `moveWeights` in `.RunMkPrimeSingleRun()`
   and pass through to `.AdaptMoveWeights()` calls
4. Unit test: verify that a move with `dim = 10` and 10× the cost of a
   `dim = 1` move gets the same score (all else equal)

**Depends on:** M-092 (ASSIGNED B). If M-092 is not yet merged, this
change layers on top of it. Coordinate with Agent B: the `dim` field is
additive and doesn't conflict with existing scheduler work.

### Task 2: Block Gibbs branch-length sweep

**Files:** `src/mcmc.cpp` (new `block_gibbs_branch_sweep_impl` + case 15
dispatch), `R/RunMkPrime.R` (`.kMoveTypes`, `.BuildMoves()`),
`R/MkPrimeMCMC.R` (`blockGibbsBranch` parameter),
`tests/testthat/test-block-gibbs-branch.R` (new)

1. Implement `block_gibbs_branch_sweep_impl()` in `src/mcmc.cpp`,
   adjacent to `weighted_branch_scale_impl` (shared infrastructure)
2. Add `case 15` to `do_move_impl()` dispatch
3. Add R-side wiring: `.kMoveTypes`, `.BuildMoves()` (with
   `dim = nEdge`), `MkPrimeMCMC()` toggle, `.AdaptTuning()` entry
4. Tests:
   - **Stationary distribution:** Run block sweep only (all other moves
     disabled) on a 5-tip tree with fixed topology. Compare posterior
     branch-length distribution against BetaSimplex-only run. K-S test
     on marginal distributions of each `relBrLengths[i]`.
   - **Reversibility:** Verify sweep with nEdge=1 (degenerate case)
     matches WeightedBranchLengthScale behaviour
   - **Tempering:** Verify with β < 1 (heated chain): run does not
     crash, acceptance rates are reasonable
   - **Permutation:** Verify edge order is shuffled (statistical test:
     correlation between sweep index and edge index ≈ 0)

**Depends on:** Task 1 (dim field), M-087 (weighted branch scale —
already done on mkp-gibbs), M-090 (move wiring — in progress).

---

## Implementation order

```
Task 1 (dim field)  ──→  Task 2 (block sweep)
                              │
                              └──→ M-091 (validation — already planned)
```

M-091 already plans to compare ESS/s across move combinations. The
block sweep becomes another arm in that comparison:
(d) + blockGibbsBranch alongside the existing (a) baseline, (b) Gibbs,
(c) Weighted comparisons.

---

## Risks and mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| O(nEdge × B) evals too slow for large trees | Medium | Scheduler downweights via dim-adjusted score; user disables via toggle. Default OFF. |
| Single-site conditional ignores inter-branch correlations | Low | Accepted limitation. Random-permutation scan reduces autocorrelation. HMC upgrade path remains if this proves insufficient. |
| Sweep mutates state in-place; crash mid-sweep leaves inconsistent state | Medium | Each within-sweep MH step is self-contained. If any eval returns non-finite, skip that pair (don't corrupt state). |
| `dim` field complicates scheduler for minimal gain on existing moves | Low | All existing moves get dim=1; the formula reduces to the status quo. Only block moves benefit. |
| Interaction with partition cache (M-064) | Medium | The sweep modifies relBrLengths, which affects ALL partitions. Must invalidate the full partition cache after the sweep, or skip caching within the sweep. Simplest: use `compute_full_loglik_at` (no cache) inside the sweep, update `state->partLogLik` once at sweep end. |

---

## Files touched

### Modified
- `src/mcmc.cpp` — new `block_gibbs_branch_sweep_impl()`, case 15
  dispatch in `do_move_impl()`
- `R/RunMkPrime.R` — `.BuildMoves()` (dim field + block sweep entry),
  `.AdaptMoveWeights()` (moveDim parameter), `.kMoveTypes`,
  `.RunMkPrimeSingleRun()` (extract/pass moveDim)
- `R/MkPrimeMCMC.R` — `blockGibbsBranch` toggle parameter

### New
- `tests/testthat/test-block-gibbs-branch.R`

### Not touched
- `src/mcmc_state.h` — no struct changes (dim is R-only)
- `src/weighted_moves.cpp` — reuse, not modify
- Existing tests — must continue to pass unchanged
