# M-111: Partial CL reuse for Gibbs subtree swap

**Status:** Plan  
**Branch:** `feature/gibbs-weighted-moves` (merged up to main `c558b0b`)  
**Prerequisite:** M-105 (partial CL for Gibbs SPR) ✅, M-109 (clone elimination) ✅

## Problem

`gibbs_subtree_swap_impl` evaluates every candidate swap partner with a
full-tree Felsenstein pruning (`compute_full_loglik_at`): O(N × C) per
candidate, where N = nEdge ≈ 105 and C = total characters across all
partition-groups.

Benchmark (Sun2018, 54 tips, 180 chars):
- SPR-only (with M-105 partial CL): 670 iter/s
- Swap-only (full eval, M-109 clone-free): 367 iter/s

Swap is now the dominant Gibbs bottleneck.

## Goal

Port the partial CL pattern from `gibbs_spr_impl` to
`gibbs_subtree_swap_impl`, reducing per-candidate cost from O(N × C) to
O(depth × C).  Expected speedup: ~15× per candidate (depth ≈ 7 vs
N ≈ 105).

## Key difference from SPR

SPR detaches a subtree and reattaches it elsewhere.  One node (u) is
removed, and one path (regraft → root) is "dirty."  The ResidualCL
structure precomputes what the tree looks like *without* the subtree.

Subtree swap exchanges two subtrees (nodeA ↔ nodeB) by swapping their
parent assignments.  **No node is removed.**  Instead, *two* parent
positions have changed children, producing *two* dirty paths that merge
at their LCA.

Crucially, ResidualCL is not applicable here — there is no single
"residual tree" because both subtrees remain.

## Swap geometry (detailed)

Given nodeA (fixed for all candidates) and candidate nodeB:

```
pA = parent(nodeA),  slotA = childSlot(pA, nodeA),  lenA = edgeLen(pA→nodeA)
pB = parent(nodeB),  slotB = childSlot(pB, nodeB),  lenB = edgeLen(pB→nodeB)
```

After swap:
- At pA: child in slotA changes from nodeA to nodeB.
  Branch length stays lenA (lengths stay at the parent position).
  **New F(pA, slotA) = P(lenA) × I(nodeB)**

- At pB: child in slotB changes from nodeB to nodeA.
  Branch length stays lenB.
  **New F(pB, slotB) = P(lenB) × I(nodeA)**

CLs of nodeA and nodeB themselves are unchanged (subtrees below them
don't move).  Only CLs at and above pA/pB are affected.

### Dirty path structure

```
         root               pathA: pA → ... → LCA → ... → root
          |                 pathB: pB → ... → LCA
         ...
          |
         LCA ← paths merge here
        /   \
      ...    ...
      |        |
     pA       pB
     |         |
   [nodeB]   [nodeA]      (after swap)
```

Three relationship cases (pA ≠ pB always, since siblings excluded):

| Case | Condition | LCA |
|------|-----------|-----|
| Independent | neither pA nor pB is ancestor of the other | true LCA of pA, pB |
| pB ancestral to pA | pB on pathA | pB itself |
| pA ancestral to pB | pA on pathB | pA itself |

All three are handled uniformly by the same bottom-up algorithm (below).

## Algorithm: `evaluate_swap_candidate()`

### Setup (once per Gibbs step, shared across all candidates):

1. Build `TreeNav` from current topology
2. Build `CLGroup`s for each (partition, kStates) unit (same as SPR)
3. Run `caching_downpass()` for all groups (and pseudo-groups for
   ascertainment) — O(N × C)
4. Precompute for nodeA (fixed across candidates):
   - `pA`, `slotA`, `lenA`
   - `pathA` = [pA, parent(pA), ..., root]
   - `onPathA[node]` = bool flags for O(1) LCA detection

### Per-candidate evaluation:

For candidate nodeB:

```
1.  pB = parentNode[nodeB], slotB = childSlot(pB, nodeB), lenB = edgeLen[nodeB→pB]

2.  Find LCA: walk from pB upward; first node with onPathA[n]=true is LCA.
    Record pathB_below_LCA = [pB, ..., node_before_LCA].
    Record lcaIdxA = index of LCA in pathA.

3.  Build dirty list (bottom-up processing order):
      segment_A = pathA[0 .. lcaIdxA-1]      (pA side, below LCA)
      segment_B = pathB_below_LCA             (pB side, below LCA)
      segment_shared = pathA[lcaIdxA .. end]  (LCA to root)

4.  For each ACRV category:
      // --- Process segment A (pA → just below LCA) ---
      For each node n in segment_A (bottom-up):
        For each child slot s:
          if (n == pA && s == slotA):
            F_new = P(lenA × rate) × I(nodeB)       ← THE SWAP
          else if child_at(n, s) is previous dirty node (from segment_A):
            F_new = P(t) × newI[previous_dirty]       ← propagated
          else:
            F_new = cached grp.Fslot(cat, n, s)       ← unchanged
        newI[n] = product of all F_new

      // --- Process segment B (pB → just below LCA) ---
      (Same structure, but swap condition is n == pB && s == slotB,
       and the swap is: F_new = P(lenB × rate) × I(nodeA))

      // --- Process shared segment (LCA → root) ---
      For each node n in segment_shared (bottom-up):
        For each child slot s:
          child_n_s = child at slot s
          if child_n_s is last node of segment_A:
            F_new = P(t) × newI[last_seg_A]
          else if child_n_s is last node of segment_B:
            F_new = P(t) × newI[last_seg_B]
          else if child_n_s is previous node on shared segment:
            F_new = P(t) × newI[previous_shared]
          else:
            F_new = cached F
          // Also: if n == pA && s == slotA and LCA == pA (pA ancestral):
          //   apply the swap here too
          // Likewise for pB
          newI[n] = product of all F_new

      // --- Root likelihood ---
      siteLik[c] += sum_s( pi_s × newI[root][c,s] )

5.  logLik = sum_c log(siteLik[c] / nCat)
```

### Handling pA-is-LCA and pB-is-LCA cases

When LCA == pA (pA is ancestral to pB), segment_A is empty (lcaIdxA = 0).
At the LCA node, two slots change simultaneously:
- slotA → swap (nodeB replaces nodeA)
- the slot toward segment_B's last node → propagated from below

When LCA == pB (pB is ancestral to pA), segment_B is empty.
At the LCA node:
- slotB → swap (nodeA replaces nodeB)
- the slot toward segment_A's last node → propagated from below

The unified loop handles this by checking swap conditions AND dirty-child
conditions at every node.  Both conditions can fire at the same node
(when LCA == pA or LCA == pB).

### Ascertainment correction

Same pattern as SPR (M-105):
1. Create pseudo-character CLGroups (`create_const_pseudo_group`)
2. Run `caching_downpass` on them
3. For each candidate, evaluate with `evaluate_swap_candidate` on
   pseudo-groups to get P(constant site)
4. Adjust: `grpLL -= nChar × log(1 - constP)`

### Relabelling correction

Topology-independent constant — add once after all groups are summed
(same as SPR).

## Implementation plan

### Step 1: New function `evaluate_swap_candidate()`

In `gibbs_partial_cl.h`, add:

```cpp
static double evaluate_swap_candidate(
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int nodeA, int nodeB,
    int pA, int slotA, double lenA,
    const std::vector<int>& pathA,
    const std::vector<bool>& onPathA);
```

Implement the 5-step per-candidate algorithm above.

Internal buffers (contrib, curI, etc.) should be pre-allocated in the
caller and passed by pointer, or allocated once per candidate as
`std::vector`s (small: stride × 2-3 buffers).  Given that stride is
typically 500–2000 doubles and we evaluate ~50–100 candidates, stack
allocation via small vectors is fine.

### Step 2: New function `evaluate_swap_const_prob()`

Parallel to `evaluate_const_prob` for SPR.  Same algorithm as
`evaluate_swap_candidate` but returns P(constant site) instead of
log-likelihood.  Can likely be a mode flag or template parameter on the
same function to avoid code duplication.

### Step 3: Refactor `gibbs_subtree_swap_impl()`

Replace the per-candidate full-evaluation loop (lines 904–932) with:

```
// --- Shared setup ---
if (data->qHeterogeneity)
    return gibbs_subtree_swap_impl_full(data, state, beta);

TreeNav topo;
topo.build(state->parent, state->child, absLen, nTip);
rates = gibbs_acrv_rates(...);

// Build CLGroups (same as gibbs_spr_impl lines 514-509)
// caching_downpass() for all groups + pseudo-groups

// Precompute pathA, onPathA
int pA = topo.parentNode[nodeA];
int slotA = topo.childSlot(pA, nodeA);
double lenA = topo.edgeLen[topo.edgeToPar[nodeA]];
pathA = [pA → root];
onPathA[n] = true for n in pathA;

// --- Per-candidate loop ---
for each partner nodeB:
    totalLL = 0;
    for each group gi:
        grpLL = evaluate_swap_candidate(groups[gi], topo, rates,
                    nodeA, nodeB, pA, slotA, lenA, pathA, onPathA);
        if (coding != 0):
            constP = evaluate_swap_const_prob(pseudoGroups[gi], ...);
            grpLL -= nChar * log(1 - constP);
        totalLL += grpLL;
    candLL[pi] = totalLL;

// Add relabelling correction (same as SPR)
// Sampling, apply (unchanged from current M-109 code)
```

### Step 4: Keep current code as `gibbs_subtree_swap_impl_full()`

Rename the current full-evaluation version to
`gibbs_subtree_swap_impl_full()` (same pattern as SPR).  Use as
fallback for Q-heterogeneity and for validation.

### Step 5: Validation test

Add a test in `test-gibbs-spr.R` (or new file) that:
1. For several random trees and random nodeA choices:
   - Evaluates all candidates with both full and partial CL
   - Asserts log-likelihoods match within floating-point tolerance (1e-8)
2. Test all three model types: neomorphic, transformational, known-k
3. Test with and without ACRV
4. Test with ascertainment correction

### Step 6: Benchmark

Run controlled benchmarks:
- Swap-only Gibbs (pre/post M-111) on Sun2018 dataset
- Mixed Gibbs (SPR+swap) to measure combined throughput
- Compare with SPR-only to check relative cost

## Complexity analysis

| Phase | Cost | Notes |
|-------|------|-------|
| TreeNav + CLGroups | O(N × C) | Same as SPR setup |
| caching_downpass | O(N × C) per group | One-time, shared |
| pathA + onPathA | O(depth) | Precomputed once |
| Per-candidate: find LCA | O(depth) | Walk from pB up |
| Per-candidate: evaluate | O(depth × nCat × stride) per group | The key win |
| Sampling + apply | O(nCand) + O(N) | Unchanged |

**Total:** O(N × C) + O(nCand × depth × C)  
**Before:** O(nCand × N × C)  
**Speedup factor:** N / depth ≈ 105 / 7 ≈ 15×

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Ancestor/descendant path cases have subtle bugs | Validation test against full evaluation (Step 5) |
| Floating-point divergence from full eval | Use same JC/MkN transition code; tolerance test |
| CLGroup setup cost becomes significant at small nCand | nCand ≈ 50–100 for 54-tip tree; setup amortised |
| Code duplication with SPR's evaluate_candidate | Accept for now; refactor to shared infrastructure in a later task if both stabilise |

## Files to modify

| File | Changes |
|------|---------|
| `src/gibbs_partial_cl.h` | Add `evaluate_swap_candidate()`, `evaluate_swap_const_prob()` |
| `src/mcmc.cpp` | Refactor `gibbs_subtree_swap_impl()`, add `_full` fallback |
| `tests/testthat/test-gibbs-spr.R` | Add partial CL validation tests for swap |
| `to-do.md` (main worktree) | Update M-111 status |
| `completed-tasks.md` (main worktree) | Archive on completion |
