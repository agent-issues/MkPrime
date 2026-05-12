# Plan: Gibbs & Weighted Tree Moves for Unrooted Trees

**Date:** 2026-03-27
**Scope:** New Phase 9 — advanced topology and branch-length proposals
**Branch:** `feature/gibbs-weighted-moves` (new worktree `mkp-gibbs`)

---

## Context

MkPrime currently has two topology moves (NNI, SPR) and one branch-length
redistribution move (BetaSimplex). All are standard MH proposals that pick
a single random candidate and accept/reject. This plan adds five new moves
from the RevBayes `mcmc_tree_moves` branch, adapted for MkPrime's unrooted
tree representation and flat-buffer likelihood engine. It also adds an
adaptive move scheduler to learn optimal move frequencies during warmup.

### Current moves

| Move | Type | Cost | Acceptance |
|------|------|------|------------|
| NNI | Topology | 1 likelihood eval | Symmetric (logHR = 0) |
| SPR | Topology | 1 likelihood eval | Asymmetric (edge-length Jacobian) |
| BetaSimplex | Branch fraction | 1 likelihood eval | Beta-ratio HR |

### New moves (from RevBayes `mcmc_tree_moves` branch)

| Move | Type | Cost per call | Key idea |
|------|------|---------------|----------|
| **GibbsSPR** | Topology | O(N) likelihood evals | Enumerate all SPR targets, sample ∝ likelihood |
| **GibbsSubtreeSwap** | Topology | O(N) likelihood evals | Enumerate all swap partners, sample ∝ likelihood |
| **WeightedBranchLengthScale** | Branch | O(B) likelihood evals | Discretize branch fraction, propose from approximate density |
| **WeightedSPR** | Topology+Branch | O(N×B) likelihood evals | Gibbs SPR + integrate over branch fractions (marginal lik) |
| **WeightedSubtreeSwap** | Topology+Branch | O(N×B) likelihood evals | Gibbs swap + integrate over branch fractions (marginal lik) |

Where N = number of candidate positions (~nEdge), B = number of discretization
breaks (default ~10).

---

## Architecture analysis

### Key blocker: likelihood evaluation inside proposals

The current `do_move_impl()` follows a strict pattern:
1. Propose (returns modified topology/parameters)
2. Evaluate likelihood ONCE
3. Accept/reject via MH

All five new moves break this pattern — they evaluate likelihood **multiple
times** within the proposal to compute sampling weights. This requires a
helper function that can evaluate the full likelihood for a given tree state
without going through the accept/reject machinery.

**Solution:** Extract a `compute_full_loglik()` static helper from the
likelihood-evaluation section of `do_move_impl()`. It takes
`(McmcData*, parent, child, edgeLen, kPrime, rateLoss, rateLogSd, rateNeo, ClWorkspace*)`
and returns `logLik`. This is essentially what lines 392–448 of `mcmc.cpp`
already do, just factored out.

### Subtree swap: new topology operation

MkPrime has no subtree swap operation. Need to implement the equivalent of
RevBayes's `swapNodes()`: given two non-nested, non-sibling nodes A and B,
detach A from its parent and B from its parent, then attach A to B's former
parent and B to A's former parent. This is simpler than SPR (no prune/suppress
/regraft) but requires careful handling of the unrooted tree's parent/child
vectors.

### MkPrime vs RevBayes tree representation

| Aspect | RevBayes | MkPrime |
|--------|----------|---------|
| Data structure | `TopologyNode*` linked objects | `IntegerVector parent, child` (1-indexed) |
| Branch lengths | Per-node `branchLength` double | `relBrLengths` simplex × `treeLength` scalar |
| Root | Virtual root (rooted representation of unrooted tree) | Node `nTip + 1` |
| Reorder | Implicit (DAG) | Explicit `postorder_order()` after topology change |

The "parent/sibling" decomposition in WeightedBranchLengthScale maps onto
MkPrime's relBrLengths: given edge `i` from parent `u` to child `v`,
the "parent branch" is the edge from `u` to `u`'s parent, and the
"sibling branch" is the other child of `u`. The total is
`treeLength * (relBr[parent_edge] + relBr[sibling_edge])`.

### Gibbs moves and parallel tempering

RevBayes's `AbstractGibbsMove` rejects heating (requires `beta == 1.0`).
The Gibbs moves here are technically **Metropolized Gibbs** — they compute
a Hastings ratio and go through normal MH acceptance. This means they
*can* be used with heated chains, though the weights computed inside the
proposal use the untempered likelihood. For tempered chains, the weights
should incorporate the temperature: `weight[i] = exp(beta * logLik[i])`.

**Decision:** Tempered Gibbs moves are correct as long as we pass `beta`
into the weight computation. The Hastings ratio formula
`log(backward_prob / forward_prob)` remains valid because the tempering
cancels in the MH ratio when both forward and backward use the same `beta`.

### Adaptive move scheduling

The new moves have wildly different cost/benefit profiles: a Gibbs SPR
costs ~35× more than a standard SPR but may not deliver 35× the ESS gain.
The optimal move mix is dataset-dependent (tree size, character count,
model complexity). Static weights will be suboptimal.

**Solution:** An adaptive scheduler that learns **ESS per wall-clock
second** for each move type during warmup and adjusts move frequencies
accordingly. See M-092 for full design.

---

## Task breakdown

### Phase 9a: Infrastructure (BLOCKER for all moves)

#### M-083: Extract `compute_full_loglik()` helper
**Priority:** P1 (blocks M-085–M-089)
**Files:** `src/mcmc.cpp`, `src/mcmc_state.h`

Extract the likelihood evaluation logic (mcmc.cpp lines 392–448) into a
standalone static function:

```cpp
static double compute_full_loglik(
    const McmcData* data, const McmcState* state,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& relBrLengths, double treeLength,
    ClWorkspace* ws);
```

This uses the state's scalar parameters (rateLoss, rateLogSd, rateNeo,
kPrime) and the given topology to compute the full log-likelihood. It
reconstructs absolute edge lengths internally (`treeLength * relBr[i]`).

Refactor `do_move_impl()` to call this helper for its default-path
likelihood evaluation. Existing behaviour must be unchanged (run existing
tests to verify).

Also add `compute_full_loglik_with_prior()` that returns logLik + logPrior,
since the Gibbs weights are proportional to the joint posterior, not just
the likelihood. (RevBayes computes `tree->getLnProbability()` +
`affected->getLnProbability()` which is likelihood + prior.)

**Tests:** Existing test suite must pass unchanged.

#### M-084: Subtree swap topology operation
**Priority:** P1 (blocks M-086, M-089)
**Files:** `src/tree_moves.cpp`

Implement `subtree_swap_impl()`:
```cpp
void subtree_swap_impl(
    IntegerVector& parent, IntegerVector& child,
    NumericVector& relBrLengths,
    int nodeA_edge, int nodeB_edge,
    int nTip);
```

Given two edge indices identifying nodes A and B (where A is not an
ancestor/descendant of B, and they are not siblings), swap A into B's
position and vice versa. Branch lengths stay with the edges (not the
nodes), so relBrLengths are unchanged — only parent/child adjacency
changes. Reorder to postorder after.

Also implement `get_valid_swap_partners()` that, given a node's edge
index, returns all valid partner edge indices (excluding ancestors,
descendants, siblings, root).

**Tests:**
- Unit test: 5-tip tree, specific swap, verify postorder-correct output
- Verify swap(A,B) followed by swap(A,B) restores original tree
- Verify invalid swaps (ancestor/descendant) are excluded

### Phase 9b: Gibbs moves

#### M-085: Gibbs SPR
**Priority:** P2 (depends on M-083)
**Files:** `src/gibbs_moves.cpp` (new), `src/mcmc.cpp`, `R/RunMkPrime.R`

Implement `gibbs_spr_impl()`:
1. Pick random non-root, non-root-child node
2. Enumerate all valid reattachment points (same logic as existing SPR,
   but exhaustive)
3. For each candidate:
   a. Perform SPR (prune + regraft)
   b. Reorder to postorder
   c. Call `compute_full_loglik_with_prior()`
   d. Store `weight[i] = exp(beta * logLik[i] + logPrior[i])`
   e. Undo the SPR (reverse regraft)
4. Sample candidate proportional to weights
5. Apply selected SPR
6. Compute Hastings ratio: `log(backward_prob / forward_prob)`

Returns: proposed parent/child/relBr vectors + logHastings.

Integration: add moveType 9 (`gibbs_spr`) to `do_move_impl()` dispatch.
The Gibbs move handles its own enumeration but returns the final state and
logHastings to the standard accept/reject machinery.

**Tests:**
- Correctness: on a 5-tip tree, verify all candidates are enumerated
- Hastings ratio: verify detailed balance by running forward/backward
- Acceptance rate: should be substantially higher than standard SPR

#### M-086: Gibbs SubtreeSwap
**Priority:** P2 (depends on M-083, M-084)
**Files:** `src/gibbs_moves.cpp`, `src/mcmc.cpp`, `R/RunMkPrime.R`

Implement `gibbs_subtree_swap_impl()`:
1. Pick random non-root node
2. Get all valid swap partners (from M-084)
3. For each candidate:
   a. Perform swap
   b. Call `compute_full_loglik_with_prior()`
   c. Store weight
   d. Undo swap
4. Sample proportional to weights
5. Apply selected swap
6. Compute Hastings ratio

Integration: add moveType 10 (`gibbs_subtree_swap`).

**Tests:** Same pattern as M-085.

### Phase 9c: Weighted moves

#### M-087: Weighted Branch Length Scale
**Priority:** P2 (depends on M-083)
**Files:** `src/weighted_moves.cpp` (new), `src/mcmc.cpp`, `R/RunMkPrime.R`

Implement `weighted_branch_length_scale_impl()`:
1. Pre-compute discretization points: Beta(0.25, 0.25) quantiles at
   `i/(num_breaks+1)` for `i = 1..num_breaks` (done once at MCMC init,
   stored in McmcData or a new WeightedMoveConfig struct)
2. Pick random internal node (not root, parent not root)
3. Identify parent edge and sibling edge; compute total branch mass
   `total = treeLength * (relBr[parent_edge] + relBr[sibling_edge])`
4. For each break point `tau_j`:
   a. Set `relBr[parent_edge] = tau_j * total / treeLength`
   b. Set `relBr[sibling_edge] = (1-tau_j) * total / treeLength`
   c. Evaluate `logLik[j]` via `compute_full_loglik()`
5. Extrapolate endpoints (j=0 and j=B+1) linearly
6. Convert to likelihoods: `lik[j] = exp(logLik[j] - max(logLik))`
7. Trapezoidal integration → normalize → CDF
8. Sample new fraction from piecewise-linear CDF (quadratic inversion)
9. Set proposed branch fractions
10. Compute Hastings ratio: `log(weight_old / weight_new)`

Integration: add moveType 11 (`weighted_branch_length_scale`). The
`num_breaks` parameter (default 10) is specified in `MkPrimeMCMC()`.

**MkPrime-specific adaptation:** RevBayes operates on absolute branch
lengths; MkPrime uses `relBrLengths` simplex + `treeLength`. The fraction
`tau = relBr[parent_edge] / (relBr[parent_edge] + relBr[sibling_edge])`
is the quantity being optimized. The total mass
`relBr[parent_edge] + relBr[sibling_edge]` stays constant.

**Tests:**
- Verify discretization points match Beta(0.25,0.25) quantiles
- Verify CDF sampling produces values in [0,1]
- Integration test: compare with BetaSimplex acceptance rates

#### M-088: Weighted SPR
**Priority:** P3 (depends on M-085, M-087)
**Files:** `src/weighted_moves.cpp`, `src/mcmc.cpp`, `R/RunMkPrime.R`

Implement `weighted_spr_impl()`:
1. Pick random node, enumerate reattachment points (like M-085)
2. For each candidate:
   a. Perform SPR
   b. Compute marginal likelihood by integrating over the branch fraction
      at the regraft point (same discretization as M-087):
      - For each break `tau_j`: set branch split, evaluate logLik
      - Trapezoidal integral → marginal `M_i`
   c. Store `weight[i] = M_i`
   d. Undo SPR
3. Sample candidate proportional to weights
4. Apply selected SPR
5. For the selected candidate, sample a branch fraction from its
   conditional density (same CDF sampling as M-087)
6. Compute Hastings ratio

**Cost:** O(N × B) likelihood evaluations per move. For a 20-tip tree
(~37 edges, ~35 candidates, 10 breaks) = ~350 evaluations. Expensive
but acceptance rate should be near 100%.

Integration: moveType 12 (`weighted_spr`).

**Tests:** Same correctness pattern. Validate marginal computation.

#### M-089: Weighted Subtree Swap
**Priority:** P3 (depends on M-086, M-087)
**Files:** `src/weighted_moves.cpp`, `src/mcmc.cpp`, `R/RunMkPrime.R`

Implement `weighted_subtree_swap_impl()`:
Same pattern as M-088 but using subtree swap instead of SPR. For each
candidate swap partner, compute the marginal likelihood by integrating
over the branch length of the swapped nodes.

Integration: moveType 13 (`weighted_subtree_swap`).

### Phase 9d: Integration, scheduling & validation

#### M-090: Wire new moves into move schedule
**Priority:** P2 (depends on M-085–M-089, but can be done incrementally)
**Files:** `R/RunMkPrime.R` (`.BuildMoves`, `.kMoveTypes`), `src/mcmc.cpp`,
`R/MkPrimeMCMC.R`

- Add new move types to `.kMoveTypes` vector (must match C++ exactly)
- Add entries to `.BuildMoves()` with initial weights (used as starting
  point for the adaptive scheduler):
  - Gibbs SPR: weight ~1
  - Gibbs SubtreeSwap: weight ~1
  - WeightedBranchLength: weight ~nEdge/5
  - Weighted SPR: weight ~0.5
  - Weighted SubtreeSwap: weight ~0.5
- Add `numBreaks` parameter to `MkPrimeMCMC()` (default: 10)
- Handle moveType dispatch in `run_mcmc_batch_cpp()` and `do_move_impl()`
- Update acceptance tracking for new move types

#### M-092: Adaptive move scheduler (cost-aware)
**Priority:** P2 (depends on M-090; can be developed in parallel with
individual moves once the scheduling interface is defined)
**Files:** `src/mcmc.cpp`, `src/mcmc_state.h`, `R/RunMkPrime.R`,
`R/MkPrimeMCMC.R`

##### Problem

The optimal move frequency mix depends on the dataset. A Gibbs SPR that
costs 35× more than a standard SPR is only worth calling if the ESS gain
per call exceeds 35× — and this varies with tree size, character count,
and model complexity. Static weights cannot account for this.

##### Design: periodic cost-aware weight re-estimation

The scheduler maximises **ESS per wall-clock second** across all move
types using a softmax score updated during warmup.

**State (per chain, per move type `m`):**

```cpp
struct MoveStats {
  int    n_accepted;        // cumulative accepts this window
  int    n_proposed;        // cumulative proposals this window
  double wall_time_ns;      // cumulative wall time this window (nanoseconds)
  double score;             // current softmax score (persists across windows)
};
```

Stored in a `std::vector<MoveStats>` of length `nMoveTypes` inside each
chain's state, or alongside the existing `accept_counts`/`propose_counts`
matrices in `run_mcmc_batch_cpp()`.

**Per-move timing:** Wrap each `do_move_impl()` call with a
`std::chrono::high_resolution_clock` before/after pair. Accumulate into
`wall_time_ns`. This adds ~10 ns overhead per call (negligible vs. a
single likelihood evaluation at ~1–100 μs).

**Update rule (every adaptation window, during warmup only):**

At each adaptation boundary (every 200 iterations, same cadence as the
existing proposal-width adaptation):

1. Compute per-move score:
   ```
   accept_rate[m] = n_accepted[m] / max(n_proposed[m], 1)
   mean_cost_s[m] = (wall_time_ns[m] / max(n_proposed[m], 1)) / 1e9
   score[m]       = accept_rate[m] / max(mean_cost_s[m], 1e-9)
   ```
   This is the **acceptance rate per second of wall time** — a proxy for
   ESS/s that doesn't require computing actual ESS (which is expensive
   and requires storing a trace window).

2. Update weights via softmax with temperature `tau`:
   ```
   w[m] = max(w_min, exp(score[m] / tau))
   normalise so sum(w) = 1
   ```
   - `tau` controls exploration vs exploitation. Start with `tau` high
     (more uniform) and anneal toward a smaller value as warmup
     progresses. E.g. `tau = tau_init * (1 - warmup_frac)`.
   - `w_min` is a floor ensuring every move gets proposed at least
     occasionally for continued learning. Default: `1 / (10 * nMoves)`.

3. Reset per-window counters: `n_accepted = n_proposed = wall_time_ns = 0`.

4. Write updated weights into the `moveWeights` vector used by the
   weighted-sampling loop in `run_mcmc_batch_cpp()`.

**At warmup end:** Freeze weights. The final adapted weights apply for
the entire post-warmup sampling phase. This preserves detailed balance
(the transition kernel is fixed).

**User overrides:**
- `MkPrimeMCMC(adaptMoveWeights = TRUE)` (default TRUE when Gibbs/Weighted
  moves are enabled, FALSE otherwise for backward compatibility).
- `MkPrimeMCMC(moveWeights = list(gibbs_spr = 2, spr = 1, ...))` pins
  specific weights and exempts them from adaptation. Only moves without
  a user-specified weight are adapted.

**Reporting:**
- During warmup progress display (`.PrintProgressTable()`), optionally
  show the current adapted move weights (one line per move, weight + accept%).
- `print.MkPosterior()` includes final adapted weights in the summary.

**Why not a formal bandit?**

- The arm set is small (~8–13 moves) and costs are directly measurable,
  so the exploration problem is mild.
- ESS attribution to individual moves is fundamentally noisy in MCMC —
  a good topology move unlocks mixing for all subsequent parameter moves.
  Acceptance-rate/cost is a cheaper, more stable proxy.
- EXP3 or Thompson sampling add machinery without clear benefit given the
  small arm count and the fact that we freeze at warmup end anyway.
- If the simple score proves inadequate (e.g. a move with high acceptance
  but zero ESS contribution dominates), the natural upgrade is to replace
  `accept_rate / cost` with a windowed ESS estimate from the trace
  buffer (but this is substantially more expensive to compute).

**Tests:**
- Unit: verify score computation and softmax normalisation
- Integration: run a short chain with `adaptMoveWeights = TRUE`; verify
  that expensive moves (Gibbs/Weighted) get downweighted relative to
  cheap moves unless their acceptance rate compensates
- Verify weights are frozen after warmup (sample from a known posterior;
  compare with fixed-weight run)

#### M-091: Validation and mixing comparison
**Priority:** P3 (depends on M-090, M-092)
**Files:** `tests/testthat/`, possibly `benchmark/`

- Test that Gibbs SPR produces same stationary distribution as standard
  SPR (run both, compare posterior summaries)
- Test that WeightedBranchLength produces same stationary distribution as
  BetaSimplex
- Benchmark ESS/second for old vs new moves on the hyoliths dataset
- Compare adaptive scheduler vs static weights on a moderate dataset
- Verify tempered chains (beta < 1) work correctly with Gibbs moves

---

## Implementation order (dependency graph)

```
M-083 (extract compute_full_loglik)  ─── BLOCKER for everything
  │
  ├── M-085 (Gibbs SPR)
  │     └── M-088 (Weighted SPR) ─── depends on M-087 too
  │
  ├── M-084 (subtree swap operation)
  │     └── M-086 (Gibbs SubtreeSwap)
  │           └── M-089 (Weighted SubtreeSwap) ─── depends on M-087 too
  │
  └── M-087 (Weighted BranchLength Scale)
        ├── M-088 (Weighted SPR)
        └── M-089 (Weighted SubtreeSwap)

M-090 (wiring) ─── incremental, after each move lands
M-092 (scheduler) ─── after M-090, or in parallel once interface is defined
M-091 (validation) ─── after M-090 + M-092
```

**Recommended build order:**
1. M-083 (infrastructure — unblocks everything)
2. M-084 (subtree swap — simple, unblocks M-086/M-089)
3. M-085 (Gibbs SPR — first Gibbs move, validates the approach)
4. M-086 (Gibbs SubtreeSwap — uses same pattern as M-085)
5. M-087 (Weighted BranchLength — first weighted move)
6. M-088 (Weighted SPR — combines M-085 + M-087)
7. M-089 (Weighted SubtreeSwap — combines M-086 + M-087)
8. M-090 (wiring — can be done incrementally after each move)
9. M-092 (adaptive scheduler)
10. M-091 (validation)

---

## Worktree setup

```bash
cd C:/Users/pjjg18/GitHub/mkp
git stash  # stash any uncommitted changes
git checkout -b feature/gibbs-weighted-moves
git worktree add ../mkp-gibbs feature/gibbs-weighted-moves
```

Update parent `AGENTS.md` worktree table and `mkp/coordination.md`.

---

## Risks and mitigations

| Risk | Impact | Mitigation |
|------|--------|------------|
| O(N) likelihood evals per Gibbs move too slow for large trees | High | Adaptive scheduler downweights expensive moves when ESS/s is low. User can disable via `moveWeights`. |
| Numerical underflow in weight computation | Medium | Use log-sum-exp: subtract max(logLik) before exponentiation (RevBayes does this). |
| Subtree swap may not improve mixing over SPR | Low | Adaptive scheduler will naturally reduce its weight if acceptance/cost is poor. |
| Tempered Gibbs weights incorrect | High | Pass `beta` into weight computation; verify with known posterior on simulated data. |
| ClWorkspace invalidation during multi-eval | Medium | Each candidate eval must reset ClWorkspace init flags; or use a dedicated scratch workspace. |
| Adaptive scheduler locks onto suboptimal mix early | Medium | Temperature annealing + `w_min` floor ensures continued exploration. Freeze only at warmup end. |
| `accept_rate / cost` is a poor proxy for ESS/s | Medium | Sufficient for small arm sets. Upgrade path: windowed ESS estimate (deferred). |

---

## Files touched

### New files
- `src/gibbs_moves.cpp` — Gibbs SPR and Gibbs SubtreeSwap implementations
- `src/weighted_moves.cpp` — WeightedBranchLength, WeightedSPR, WeightedSubtreeSwap
- `tests/testthat/test-gibbs-moves.R`
- `tests/testthat/test-weighted-moves.R`
- `tests/testthat/test-move-scheduler.R`

### Modified files
- `src/mcmc.cpp` — new moveType dispatch, `compute_full_loglik()` helper,
  per-move timing, adaptive weight update in batch loop
- `src/mcmc_state.h` — `MoveStats` struct, WeightedMoveConfig, new declarations
- `src/tree_moves.cpp` — subtree swap operation
- `R/RunMkPrime.R` — `.BuildMoves()`, `.kMoveTypes`, move dispatch,
  scheduler integration
- `R/MkPrimeMCMC.R` — new parameters (numBreaks, adaptMoveWeights, moveWeights)
- `R/MkPosterior.R` — report adapted move weights in `print()`/`summary()`
- `NAMESPACE` — (auto-generated by roxygen)
