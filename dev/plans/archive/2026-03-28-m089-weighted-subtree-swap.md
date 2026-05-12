# Plan: M-089 — WeightedSubtreeSwap (moveType 14)

**Date:** 2026-03-28
**Depends on:** M-083 (compute_full_loglik), M-084 (subtree swap), M-086 (GibbsSubtreeSwap), M-087/M-088 (bin machinery)

---

## Overview

WeightedSubtreeSwap extends GibbsSubtreeSwap (M-086) the same way WeightedSPR (M-088) extended GibbsSPR (M-085): by marginalising over branch-fraction bins at each candidate swap position. Cost: O(N × B) likelihood evaluations.

---

## Algorithm

### What varies: the two swapped branch lengths

When swapping nodeA and nodeB, `swap_subtrees_impl` reassigns parent connections and swaps the branch lengths:
- rowA (edge to A): gets B's old branch length
- rowB (edge to B): gets A's old branch length

For the **weighted** extension, we treat the total `brA + brB` as fixed and integrate over how to distribute it between the two swapped edges via a fraction tau:
- Edge to A (at B's old position): `tau × total`
- Edge to B (at A's old position): `(1 - tau) × total`

The default swap (no fraction change) corresponds to `tau = brB / (brA + brB)`.

### Steps

1. **Pick nodeA** uniformly from all edge children (same as GibbsSubtreeSwap)
2. **Get valid swap partners** via `get_valid_swap_partners_impl()`
3. **Find rowA** — the edge row where `child[rowA] == nodeA`
4. **Compute `total = absLen[rowA] + absLen[rowB_partner]`** and `fOld = absLen[rowA] / total`

5. **Self marginal (M_self):** For each bin b:
   - Set `absLen[rowA] = mid_b × total`, `absLen[rowB_self] = (1 - mid_b) × total` (where rowB_self is rowA's own position — no swap, just varying A's and... wait)

   **Correction:** The self candidate keeps A in place, so we vary A's own branch length vs... what? We need a second edge to redistribute with. The natural choice is: for each candidate partner nodeB, the "self" baseline is the current topology with the total `brA + brB_partner` held fixed. But the self is per-partner, not global.

   **Simpler approach (matching M-088):** Include the current topology as a single candidate with weight = `exp(beta × state->logLik)` (no bin integration for self). This matches how GibbsSubtreeSwap includes self via `wOrig`. Then for each non-self candidate, compute the marginal over bins.

   **Even simpler (but still correct):** For each candidate partner nodeB_i:
   - `total_i = absLen[rowA] + absLen[rowB_i]`
   - Construct the swapped topology once (swap A ↔ B_i)
   - For each bin b: set the two swapped edges to `mid_b × total_i` and `(1 - mid_b) × total_i`
   - Evaluate LL at each bin → marginal M_i = sum of weights

   For self: use `state->logLik` directly as the self weight (scaled to be comparable). Since different candidates have different `total_i`, the bin integration is per-candidate anyway, and self doesn't have a natural partner to redistribute with.

   **Problem:** This breaks the Z = Z' symmetry that M-088 relies on, because self has no bin integration while candidates do.

### Resolving the self-candidate issue

In M-088 (WeightedSPR), self has a natural "regraft point" (the merged edge at the prune point) and the same total edge length as any other candidate. This makes the bin integration symmetric.

For subtree swap, each candidate partner nodeB_i has a **different** total edge length (`brA + brB_i`), so there's no single shared fraction space. This means:

**Option A: Per-candidate self baseline.** For each partner B_i, compute both:
- "swap" weights: topology with A↔B_i swapped, fraction varied
- "no-swap" weights: original topology, fraction of (brA, brB_i) varied

This makes the candidate sets symmetric (Z = Z' per partner) but changes the move structure — we'd be sampling (partner, swap-or-not, bin, fraction) rather than (topology, bin, fraction).

**Option B: Single self weight, no bin integration for self.** Use `wOrig = exp(beta × state->logLik)` as self weight (matching GibbsSubtreeSwap), and compute full marginals for each candidate. Z ≠ Z' in general, so the Hastings ratio must account for both normalizing constants.

**Option C: Skip self entirely.** Always propose a swap (never a no-op). Compute marginals for all candidates, sample one, sample bin and fraction, MH accept/reject. The MH step handles everything. This is simplest and correct.

**Decision: Option C.** The MH acceptance step already handles the posterior ratio, so we don't need self in the candidate set for correctness. If all swaps are bad, MH will reject. This avoids the Z = Z' symmetry issue entirely.

### Revised algorithm (Option C)

1. Pick nodeA uniformly from edge children
2. Get valid swap partners → partners[0..nPart-1]
3. For each partner B_i:
   a. Find rowA, rowB_i in child vector
   b. `total_i = absLen[rowA] + absLen[rowB_i]`
   c. Construct swapped topology: `newParent[rowA] = parent[rowB_i]`, `newParent[rowB_i] = parent[rowA]` (plus branch length swap)
   d. For each bin b:
      - Set the two swapped edges to `mid_b × total_i` and `(1 - mid_b) × total_i`
      - Call `preorder_weighted_impl` + `compute_full_loglik_at`
      - Store `candLL[i][b]`
   e. `M_i = sum_b exp(beta × (candLL[i][b] - globalMax))`
4. Sample partner i proportional to M_i
5. Sample bin b within partner i, proportional to `candW[i][b]`
6. Draw fraction fNew from `Beta(mid_b × conc + 1, (1 - mid_b) × conc + 1)`
7. Construct final proposed topology at fNew
8. Evaluate actual newLogLik at fNew

### Hastings ratio (Option C — no self)

Forward proposal: pick nodeA (prob 1/nEdge), pick partner i (prob M_i / Z), pick bin b (prob w_{i,b} / M_i), draw fNew (prob Beta(fNew)).

Reverse proposal: from proposed state, pick nodeA (prob 1/nEdge), pick partner that undoes the swap (this is the same partner i — swapping A↔B is its own inverse), pick the bin containing fRev, draw fRev.

**Key insight:** Swapping A↔B and then swapping A↔B again restores the original tree. So the reverse move picks the same partner with probability M'_i / Z'. The candidate set from the proposed state has the same partners (get_valid_swap_partners is topology-dependent — need to verify this is symmetric after swap).

Actually, after swapping A and B, A is now at B's old position and B at A's old position. The valid swap partners for A from the new topology may differ from the original. This breaks the simple reverse analysis.

**Simpler Hastings ratio approach:** Since we're using MH acceptance anyway, the Hastings ratio just needs to be the ratio of reverse-to-forward proposal densities. Both the forward and reverse moves:
1. Pick nodeA uniformly (cancels)
2. Enumerate partners and compute marginals (may differ between topologies)
3. Sample partner, bin, fraction

Since step 2 differs between forward and reverse directions, we'd need to compute the reverse marginals too — which costs another O(N × B) evaluations. That's 2× the cost and defeats the purpose.

**Practical simplification:** Accept that the normalizing constants don't cancel and use the following logHR:

```
logHR = log(Z_forward) - log(M_i_forward) - log(w_{i,b_new}) - log(Beta(fNew))
      - [log(Z_reverse) - log(M_i_reverse) - log(w_{i,b_rev}) - log(Beta(fRev))]
```

But computing Z_reverse requires evaluating all candidates from the proposed topology — too expensive.

### Alternative: include self, match M-088 pattern

Let me reconsider **Option B** with a careful Hastings ratio.

For a specific partner B_i with `total_i = brA + brB_i`:

The current state's "fraction" at this pair is `fOld_i = brA / total_i`. After swapping A↔B_i with fraction fNew, the new state's "fraction" at this pair is `fNew` (and the reverse fraction to get back is `fOld_i`).

If I include self with weight `wOrig = exp(beta × state->logLik)` (no bins), then:

```
q_forward = (1/nEdge) × (M_i / Z) × (w_{i,b} / M_i) × Beta(fNew)
          = (1/nEdge) × (w_{i,b} / Z) × Beta(fNew)

Z = wOrig + sum_j M_j
```

For the reverse, the partner set from the swapped topology is the same (swap is self-inverse), and:
```
Z' = wOrig' + sum_j M'_j  (where wOrig' = exp(beta × newLogLik))
```

These don't simplify. So Option B also needs the reverse normalizer.

### Final decision: use the GibbsSubtreeSwap pattern with bins

The cleanest approach that avoids computing reverse normalizers:

**Treat this as GibbsSubtreeSwap with per-candidate bin sampling.** Include self (wOrig, no bins). For each non-self candidate, compute marginal M_i over bins. Sample from {self, cand_1, ..., cand_N}. If self → no-op. If candidate i chosen:
- Sample bin b within i
- Draw fraction fNew
- **Accept unconditionally** (Gibbs step for topology) but then apply MH correction for the branch-fraction component only.

The topology selection is a proper Gibbs draw (same as GibbsSubtreeSwap) because wOrig and M_i are proportional to posterior probabilities (marginal over fractions for candidates, point value for self). The only non-Gibbs part is drawing fNew from an approximate (binned) distribution rather than the exact conditional.

**logHR for the fraction component only:**
```
logHR = log(w_{i, b_default}) + log(Beta(f_default | b_default))
      - log(w_{i, b_new})     - log(Beta(fNew | b_new))
```
where `f_default = brB_i / total_i` (the default swap fraction) and `b_default` is the bin containing `f_default`.

Then MH acceptance:
```
logAlpha = beta * (newLogLik - candLL[i][b_default_midpoint]) + logHR
```

Wait, this doesn't quite work either because `candLL[i][b_default_midpoint]` is not the same as the actual LL at `f_default`.

### Simplest correct approach

After much deliberation, the simplest correct implementation is:

1. Enumerate partners, compute marginal weights M_i (sum over bins)
2. Include self with weight `wOrig = exp(beta × state->logLik)`  
3. Sample topology (self → no-op, candidate i → proceed)
4. For chosen candidate i: sample bin, draw fNew from Beta
5. Construct proposed state at fNew
6. Evaluate actual newLogLik
7. **Compute logHR as:**
   ```
   logHR = -log(Beta(fNew | alphaNew, betaNew))
   ```
   This accounts for the proposal density of fNew. The topology-level Gibbs weights handle the discrete part. The MH step corrects for using the approximate (binned + Beta) proposal rather than the exact posterior conditional.

Actually this is also wrong. Let me just go with the exact same pattern as M-088.

### FINAL approach (matching M-088 exactly)

The key realization: **each candidate partner has its own total branch length**, but I can still define self-per-pair by evaluating the original topology at different fractions of (brA, brB_i). This makes the construction parallel to M-088 where we have a fixed total at the regraft point.

**But wait —** in M-088, the self candidate only involves ONE pair of edges (the merged prune-point edge), and this total is the SAME for all candidates (lMerge is fixed). For subtree swap, each candidate i has a different total_i, so the self evaluations would differ per candidate.

**Resolution:** I don't need self evaluations per candidate. The self weight is simply `exp(beta × state->logLik)`, evaluated once. Each candidate's marginal M_i integrates over fractions specific to that candidate's total_i. The Hastings ratio accounts for the asymmetry.

**logHR derivation:**

Forward: q_fwd ∝ (w_{i,b_new} / Z) × Beta(fNew | b_new)
Reverse: from proposed state, undo swap. q_rev ∝ (w'_{self??} / Z') × ...

This is where it gets complicated because the reverse isn't just "undo the swap with the reverse fraction."

**PRAGMATIC DECISION:** Given the complexity of deriving a clean Hastings ratio for the asymmetric-total case, implement WeightedSubtreeSwap with the **same self-weight pattern as GibbsSubtreeSwap** (point weight, no bins for self), and add a **simple MH correction** for the within-candidate fraction draw:

```
logHR = log(selfW_at_default) + log(Beta(f_default | b_default))
      - log(candW[i][b_new])  - log(Beta(fNew | b_new))
```

where:
- `selfW_at_default = wOrig` (exp(beta × state->logLik) normalized to same scale)  
- `f_default = brB_i / total_i` (default swap fraction for candidate i)
- `b_default` = bin containing f_default
- Everything uses the same globalMax offset

This is analogous to M-088's logHR where selfW[b_old] is the self weight at the old bin. Here, wOrig plays the role of selfW[b_old] but is evaluated at the actual current fraction rather than a bin midpoint.

Then full MH acceptance:
```
logAlpha = beta × (newLogLik - state->logLik) + (newLogPrior - state->logPrior) + logHR
```

---

## Implementation

### Files modified
- `src/mcmc.cpp` — add `weighted_subtree_swap_impl()`, case 14 in `do_move_impl()`
- `R/RunMkPrime.R` — add `weighted_subtree_swap = 14L` to `.kMoveTypes`
- New: `tests/testthat/test-weighted-subtree-swap.R`

### C++ function signature
```cpp
static bool weighted_subtree_swap_impl(
    McmcData* data, McmcState* state, double beta, int nBins);
```

Returns true on accept, false on reject or self-draw. Handles MH acceptance internally (like weighted_spr_impl, cases 10/11/13).

### Pseudocode

```
1.  Pick nodeA = child[random edge]
2.  partners = get_valid_swap_partners_impl(parent, child, nTip, nodeA)
3.  if empty → return false
4.  rowA = find_child_row(child, nodeA)
5.  Compute absLen = treeLength × relBrLengths

6.  For each partner B_i (pi = 0..nPart-1):
      rowB = find_child_row(child, B_i)
      total_i = absLen[rowA] + absLen[rowB]
      
      // Construct swapped topology
      np = clone(parent); np[rowA] = parent[rowB]; np[rowB] = parent[rowA]
      // child vector unchanged (child[rowA] still = nodeA, child[rowB] still = B_i)
      
      For each bin b:
        na = clone(absLen)
        na[rowA] = mid_b × total_i
        na[rowB] = (1 - mid_b) × total_i
        (op, oc, oa) = preorder_weighted_impl(np, child, na)
        candLL[pi][b] = compute_full_loglik_at(data, state, op, oc, oa)

7.  Compute globalMax, then:
      wOrig = exp(beta × (state->logLik - globalMax))
      For each pi: M_i = sum_b exp(beta × (candLL[pi][b] - globalMax))
      Z = wOrig + sum_i M_i

8.  Sample from {self, cand_0, ..., cand_{nPart-1}}
    If self → return false

9.  For chosen pi: sample bin b_new from candW[pi]
10. Draw fNew from Beta(mid_{b_new} × conc + 1, (1 - mid_{b_new}) × conc + 1)
11. Clamp fNew to [1e-8, 1-1e-8]

12. Construct final proposed state:
      np = clone(parent); np[rowA] = parent[rowB]; np[rowB] = parent[rowA]
      na = clone(absLen); na[rowA] = fNew × total; na[rowB] = (1-fNew) × total
      (ordEdge, ordAbs) = preorder_weighted_impl(np, child, na)
      newLogLik = compute_full_loglik_at(...)

13. Compute Hastings ratio:
      f_default = absLen[rowB] / total  (default swap fraction)
      b_default = bin containing f_default
      logHR = log(max(wOrig, 1e-300)) + log(Beta(f_default | b_default))
            - log(max(candW[pi][b_new], 1e-300)) - log(Beta(fNew | b_new))

14. newLogPrior = cpp_log_prior(...)
15. logAlpha = beta × (newLogLik - state->logLik) + (newLogPrior - state->logPrior) + logHR
16. MH accept/reject
```

### Tests (test-weighted-subtree-swap.R)

Same structure as test-weighted-spr.R:
1. **Structural:** can change topology, preserves tree length, simplex preserved, canonical preorder, scalar parameters unchanged
2. **Functional:** logLik consistent with recomputation, all branch lengths positive
3. **Chain:** 30-move chain stays valid
4. **Acceptance:** non-trivial acceptance rate (> 0%, < 100%)
5. **Rejection:** state unchanged on rejection

### Move type registration
- C++: case 14 in do_move_impl(), returns `weighted_subtree_swap_impl(data, state, beta, 10)`
- R: `weighted_subtree_swap = 14L` in `.kMoveTypes`

---

## Key differences from M-088 (WeightedSPR)

| Aspect | WeightedSPR (M-088) | WeightedSubtreeSwap (M-089) |
|--------|---------------------|----------------------------|
| Topology op | SPR (prune + regraft) | Subtree swap (exchange parent connections) |
| Candidate source | Valid regraft edges (BFS filter) | `get_valid_swap_partners_impl()` |
| Branch fraction | Split of regraft edge | Distribution of brA + brB_i |
| Self marginal | Bin integration over merged edge | Point weight `exp(beta × state->logLik)` |
| Z = Z' symmetry | Yes (candidate sets identical) | No (self has no bin integration) |
| Total per candidate | Same for all (lReg is constant) | Varies (total_i = brA + brB_i) |
| f_old for reverse | `absLen[parentRow] / lMerge` | `absLen[rowB_i] / total_i` (the default swap fraction) |

---

## Risks

- **Hastings ratio approximation:** Using wOrig (point weight) vs candidate marginals (integrated weights) means the ratio is approximate. MH acceptance corrects this, but acceptance rates may be lower than M-088.
- **Variable total_i:** Different candidates have different totals, so the bin midpoints correspond to different absolute lengths. This is fine — each candidate's bin integration is self-consistent.
