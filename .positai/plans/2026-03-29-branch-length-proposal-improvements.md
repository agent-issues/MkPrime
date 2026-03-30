# Plan: Branch-Length Proposal Improvements (M-127)

## Problem

M-125 added a block Dirichlet branch-length proposal that improved
branch ESS by ~4.5×. But two limitations remain:

1. **No partial evaluation.** The Dirichlet move falls through to
   full likelihood evaluation, making it ~3–5× more expensive per
   step than beta_simplex (which uses partial CL caching). The
   adaptive scheduler doesn't know this — it sees one likelihood
   evaluation either way.

2. **Random edge selection.** Selecting K edges uniformly at random
   wastes proposal mass on uncorrelated edges. Branch lengths are
   most correlated *locally* (neighboring edges anti-correlate).

## Approach: Three-phase incremental improvement

All work stays on `main` — no worktree needed. The changes are
small, incremental, and all touch the same files as M-125. No
parallel agent work is active on these files.

---

### Phase A: K sweep (no code changes)

**Goal:** Find optimal K for current random Dirichlet.

Run Vinther2008 (fixed topology, 15k iter) with K ∈ {3, 5, 10, 15, nEdge}.
Compare ESS/s for `branch_lengths` and `log_likelihood`.

This informs the default K and tells us whether the Dirichlet is
cost-effective at all before investing in partial eval.

**Deliverable:** Table of K vs ESS/s; update default if warranted.

---

### Phase B: Partial evaluation for Dirichlet (`find_dirty_dirichlet`)

**Goal:** Make the Dirichlet proposal cheap enough to compete with
beta_simplex on a per-step basis.

#### B1. Pass modified edge indices out of `dirichlet_simplex_impl`

Currently the proposal shuffles indices internally and doesn't
expose which edges were modified. Add an output parameter:

```cpp
bool dirichlet_simplex_impl(NumericVector& x, int nCats, double alpha,
                            double& logHastings, NumericVector& snapshot,
                            std::vector<int>& modifiedEdges);  // NEW
```

The first `nCats` elements of the shuffled `indices` array are the
modified edge indices — just copy them to `modifiedEdges`.

#### B2. `find_dirty_dirichlet` in `node_cl_cache.h`

Generalization of `find_dirty_beta_simplex` from 2 paths to K paths:

1. For each modified edge index, get parent node from `parent[idx]`.
2. Walk each parent to root, marking nodes in a `bool` array.
3. Collect all marked nodes, sort by depth (deepest first = postorder).

This is the same pattern as `find_dirty_beta_simplex` but with K
starting points instead of 2. When K is small (3–5), the dirty set
is much smaller than the full tree. When K ≈ nEdge, it degrades to
full eval — add a heuristic: if `dirty.size() > 0.8 * nInternalNodes`,
fall back to full eval (avoids the overhead of save/restore for no
benefit).

#### B3. Wire into `do_move_impl` (mcmc.cpp)

Add a branch for `moveType == 23 && state->nodeCL.valid`:

```cpp
} else if (moveType == 23 && state->nodeCL.valid) {
    auto dirty = find_dirty_dirichlet(state->nodeCL.topo,
                                       evalParent, modifiedEdges);
    if ((int)dirty.size() > (int)(0.8 * nInternalNodes)) {
        // Fall back to full eval
        ...
    } else {
        newLogLik = partial_eval_dirty(..., dirty);
        usedPartialCL = true;
    }
}
```

On rejection with `usedPartialCL`, the existing `restore_dirty_cls`
handles rollback.

#### B4. Store `modifiedEdges` in McmcState

Need a scratch buffer in `McmcState` to hold the K modified indices
between proposal and evaluation. Add `std::vector<int> dirEdges`
alongside `brSnapshot`.

#### B5. Tests

- Unit test: `find_dirty_dirichlet` with K=2 matches
  `find_dirty_beta_simplex` output.
- Unit test: postorder invariant holds for K=3,5,10.
- Integration test: partial eval logLik matches full eval logLik
  for Dirichlet proposals (same pattern as existing beta_simplex
  partial eval tests).
- Regression: full test suite still passes.

**Deliverable:** Dirichlet proposal uses partial CL evaluation.
Re-run K sweep to measure speedup.

---

### Phase C: Localized Dirichlet (neighborhood selection)

**Goal:** Exploit local correlation structure by selecting connected
edges instead of random edges.

#### C1. `select_neighborhood` in proposals.cpp

Given a random starting node, collect K edges in its local
neighborhood using BFS:

1. Pick a random internal node `u`.
2. BFS outward from `u` along the tree (both toward tips and root),
   collecting edges until K edges are gathered.
3. Return the K edge indices.

This ensures all selected edges are topologically close, so:
- The proposal targets the most correlated branch lengths.
- The dirty set is compact (one connected region → one path to root),
  making partial eval maximally effective.

#### C2. New move type or parameter?

Two options:

**Option 1: New move type (e.g., `local_dirichlet`, case 24).**
Cleaner separation of concerns. The scheduler can independently
weight random vs. localized Dirichlet.

**Option 2: Parameter on existing move (e.g., `localized = TRUE`).**
Less code, but conflates two different proposals in one move type.

**Decision: Option 1.** The scheduler already handles multiple
branch-length proposals (beta_simplex, dirichlet_branch, Gibbs
branch). A third is architecturally consistent and lets adaptation
find the right mix.

#### C3. Integration

- Add case 24 in `do_move_impl`
- Add to `.kMoveTypes`, `.BuildMoves`, `.BuildScaleTuningMatrix`,
  `.AdaptTuning`
- Default: included alongside random Dirichlet, with
  `weight = max(1, nEdge / 4)` (same as current Dirichlet)
- K for local Dirichlet: `min(nEdge, 6)` (smaller than random
  Dirichlet because local edges are more correlated, so fewer are
  needed for effective mixing)

#### C4. Tests

- Unit test: `select_neighborhood` returns K connected edges.
- Unit test: dirty set from localized selection is smaller than from
  random selection (on average).
- Integration: MCMC with localized Dirichlet produces valid posterior
  (stationary distribution test or comparison with known result).
- Benchmark: ESS/s comparison with and without localized Dirichlet.

**Deliverable:** Localized Dirichlet move with partial eval. Final
ESS/s comparison across all branch-length proposals.

---

## What this does NOT include

- **Adaptive covariance (ALR-space Gaussian):** High implementation
  cost, uncertain payoff. Defer to Phase 5+.
- **HMC on branch lengths:** Major architectural change. Defer.
- **Full remainder-plus-selected redistribution:** The Jacobian
  penalty is intrinsic; this design is fundamentally flawed for
  n >> K.

## Worktree decision

**No worktree.** Rationale:
- Single-agent work, no parallel development on these files
- All changes are incremental on M-125 (same files)
- Phases A/B/C are sequential and small enough to commit directly
  to `main`

## Files modified

| File | Change |
|------|--------|
| `src/proposals.cpp` | `modifiedEdges` output param; `select_neighborhood` |
| `src/node_cl_cache.h` | `find_dirty_dirichlet` |
| `src/mcmc.cpp` | Case 23 partial eval; case 24 dispatch; `dirEdges` in state |
| `R/RunMkPrime.R` | `.kMoveTypes`, `.BuildMoves`, tuning, adaptation for local_dirichlet |
| `R/MkPrimeMCMC.R` | `localDirichlet` parameter (if separate from `dirichletBranch`) |
| `tests/testthat/test-m127-*.R` | Tests for all phases |

## Task ID

M-127: Branch-length proposal improvements (K sweep + partial eval + localized Dirichlet)
