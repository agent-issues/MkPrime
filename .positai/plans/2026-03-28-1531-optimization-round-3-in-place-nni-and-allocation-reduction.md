# Optimization Round 3: In-Place NNI and Allocation Reduction

**Agent:** C  
**Date:** 2026-03-28  
**Status:** COMPLETE  
**Related:** M-106 (round 2 complete), M-108 (postorder cache), PROFILING.md  

## Context

After round 2 (1.80× cumulative speedup, 1895 iter/s on Sun2018 hyoliths),
remaining C++ time is dominated by Felsenstein pruning. The next tier of
optimizations targets structural overhead *around* the pruning: the
preorder reorder required by every topology proposal, and unnecessary
vector clones in the hot path.

## Key Insight: NNI Does Not Require Edge Reordering

The Felsenstein pruning traversal iterates edges in reverse order using
`initFlg` to track first-write. It requires only that the reverse edge
order is a **valid postorder** (children processed before parents). This
does NOT require canonical preorder — any topological sort works.

NNI only changes 2 parent assignments (`newParent[cRow] = u; newParent[wRow] = v`).
Verified (concrete examples + general argument) that this preserves the
valid-postorder property:
- cRow moves from under v to under u (u is v's parent → u's edge is
  earlier in the preorder → cRow is still processed before u in reversed order) ✓
- wRow moves from under u to under v (v is v → wRow was after v's subtree
  block → wRow is processed before v in reversed order) ✓
- Property holds across multiple consecutive NNIs ✓
- All other topology proposals (SPR, TBR, Gibbs) produce canonical preorder
  via `preorder_weighted_impl`, so they "reset" the ordering naturally ✓
- Downstream R code (tree sampling, `ape::write.tree`) handles any edge
  ordering ✓

## Optimizations

### 1. In-place NNI (OPP-6b) — Primary optimization

**Files:** `src/mcmc.cpp` (case 5 in `do_move_impl`)

**Current flow:**
1. Call `nni_proposal_impl(state->parent, state->child, ...)` → clones 2 vectors
2. Inside: find internal edge, find swap targets, clone, swap parents
3. Call `preorder_weighted_impl` → full O(nEdge) reorder + postorder_order
4. Pack result into `Rcpp::List` (4 named elements)
5. Return to `do_move_impl`, unpack List into proposedParent/proposedChild/proposedRelBr
6. Evaluate likelihood using proposed vectors
7. On rejection: proposed vectors discarded (no rollback needed)
8. On acceptance: `state = std::move(proposed)` (3 vector moves)

**New flow:**
1. Inline NNI logic in case 5: find internal edge, find swap targets
2. Save 2 original parent values: `savedP_cRow`, `savedP_wRow`
3. Modify `state->parent` in-place: 2 assignments
4. Set `nniInPlace = true` flag (new variable alongside `topologyChanged`)
5. Evaluate likelihood using `state->parent`/`state->child` directly
6. On rejection: restore 2 parent values (O(1))
7. On acceptance: nothing to do (state already modified); clear partLogLik

**Eliminates per NNI proposal:**
- 2× `clone()` (parent + child vectors, each ~106 ints for 54-tip tree)
- 1× `preorder_weighted_impl` call (includes `postorder_order` — VTune's #1 C++ hotspot at 3.2% total)
- 1× `Rcpp::List` construction + destruction
- 1× `absLen` vector construction
- 3× `std::move` on acceptance
- relBrLengths vector construction (`ordAbs / treeLength`)

**Does not touch:** `nni_proposal_impl` (kept for R-exported wrapper and test compatibility).

### 2. O(1) BetaSimplex rollback

**Files:** `src/mcmc.cpp` (cases 4 and 12 in `do_move_impl`)

**Current:** `oldRelBr = clone(state->relBrLengths)` — full vector clone.
**New:** Save 2 indices and original values:
```cpp
int bsIdx1, bsIdx2;
double bsOldVal1, bsOldVal2;
```
On rejection: `state->relBrLengths[bsIdx1] = bsOldVal1; state->relBrLengths[bsIdx2] = bsOldVal2;`

**Problem:** `beta_simplex_impl` internally picks the `other` index via
`unif_rand()`, so the caller doesn't know which 2 elements changed. Two options:
- (a) Add output parameters to `beta_simplex_impl` for the two indices, or
- (b) Wrap the call: save all of `x` before, compare after to find changes.

Option (a) is cleanest. Modify `beta_simplex_impl` signature to also accept
`int& outIdx1, int& outIdx2`.

For case 12 (weighted_branch_scale), the function modifies up to nBins+1
elements (it's a sweep, not a single pair), so clone rollback may still be
needed there. Leave case 12 alone if the sweep modifies >2 elements.

### 3. Ascertainment tip-init hoist

**Files:** `src/ascertainment.cpp` (all 4 functions)

The constant/singleton site prob functions loop over ACRV categories and
re-initialize tips every category (via `std::fill` + tip setup). Tips are
identical across categories.

**Fix:** Same pattern as the main pruning ACRV functions:
- Initialize tips once before the category loop
- In the category loop, only reset `cl_init` for internal nodes

For constant_site_prob_jc, nChar = kStates (small), so the gain is small.
For singleton_site_prob_jc/mkn, nChar = nTip or 2*nTip (larger).

### 4. TopologyProposal struct (deferred / if time permits)

Replace `Rcpp::List` returns from `spr_proposal_impl`, `tbr_proposal_impl`,
`swap_subtrees_impl` with a plain struct:
```cpp
struct TopologyProposal {
    IntegerVector parent, child;
    NumericVector relBrLengths;
    double logHastings;
};
```
This would eliminate 4× string-hashed named-element insertions per topology
proposal. Modest impact since these proposals are less frequent than NNI.

**Risk:** Requires updating all callers. Defer to a follow-up if round 3
is already large.

## Implementation Order

1. In-place NNI (highest value, most complex)
2. O(1) BetaSimplex rollback (easy, independent)
3. Ascertainment tip-init hoist (easy, independent)
4. TopologyProposal struct (time permitting)

## Testing Strategy

- Run full test suite (`test-mcmc-engine.R` — 335 tests including slow
  inference tests, plus `test-likelihood.R`, `test-proposals.R`, `test-priors.R`)
- In-place NNI correctness: the existing acceptance-rate tests and
  likelihood-validation tests will catch any postorder violations
- Benchmark: same Sun2018 workload (54 taxa, 225 chars, 20k iter) with
  `bench::mark` for before/after comparison

## Expected Impact

| Optimization | Est. savings | Confidence |
|---|---|---|
| In-place NNI | 3-8% of total | High — eliminates VTune's #1 C++ hotspot (postorder_order) + clone/List overhead |
| O(1) BetaSimplex rollback | <1% | High — eliminates one clone per BetaSimplex move |
| Ascertainment tip-init hoist | <1% | Medium — depends on nCat and nTip |
| TopologyProposal struct | <0.5% | Medium |

Combined: **~5-10%** additional speedup (~2000+ iter/s target).
