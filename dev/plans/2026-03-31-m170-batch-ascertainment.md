# M-170: Batch ascertainment correction (`constant_site_prob_jc`)

**Date:** 2026-03-31  
**Author:** Agent D  
**Status:** Ready to implement

---

## Problem

`constant_site_prob_jc(k)` accounts for **18.9% of total warmup CPU** (23.1s / 122.6s,
Sun2018 54 taxa, nCat=6, 3000 warmup iters). VTune run: `vtune-warmup/`.

The function computes P(constant site) for JC(k) by running a full Felsenstein
traversal with `nChar = kStates = k` pseudo-characters — one per constant state.

**The redundancy:** JC is state-symmetric, so all k constant patterns have exactly
equal probability. Currently we compute P(all tips = s) for s = 0..k-1 and sum,
but by symmetry they are all equal:

```
P(constant) = sum_{s=0}^{k-1} P(all tips = s) = k × P(all tips = 0)
```

We need only **1 pseudo-character** (all tips in state 0) and multiply by k.
This is a **k× speedup** within the function (k = 2..8 for Sun2018; average ~4×).

---

## Root cause

```cpp
// In ascertainment.cpp, line 39
const int nChar = kStates;        // ← k traversals per call
const int stride = nChar * kStates;
```

With nCat=6 and k=5 (e.g., kObs=3 + ko=2), the current code runs:
- 5 pseudo-characters × 5 states × 6 rate cats × 106 edges (Sun2018)

With the fix:
- 1 pseudo-character × 5 states × 6 rate cats × 106 edges

---

## Fix: single pseudo-character in `constant_site_prob_jc`

**Only `src/ascertainment.cpp` changes.** Function signature and return value are
identical; all call sites are unaffected.

### Changes in `constant_site_prob_jc`

**1. Reduce `nChar` and `stride`:**

```cpp
// OLD
const int nChar = kStates;
const int stride = nChar * kStates;

// NEW — JC symmetry: 1 pseudo-char suffices
// P(all tips = s) is identical for all s; use s=0, multiply by kStates at end.
const int nChar = 1;
const int stride = kStates;
```

**2. Simplify tip initialisation:**

```cpp
// OLD
for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = cl_flat.data() + tip * stride;
    for (int s = 0; s < kStates; ++s)
        cl[s * kStates + s] = 1.0;  // k sparse writes
    cl_init[tip] = 1;
}

// NEW — only state 0 is 1; rest already 0 from vector initialisation
for (int tip = 1; tip <= nTip; ++tip) {
    cl_flat[tip * stride] = 1.0;  // cl[0] = 1, rest already 0
    cl_init[tip] = 1;
}
```

**3. Simplify root accumulation and return value:**

```cpp
// OLD — sum over all k characters
const double* clRoot = cl_flat.data() + root * stride;
for (int s = 0; s < kStates; ++s) {
    const int offset = s * kStates;
    double site_lik = 0.0;
    for (int i = 0; i < kStates; ++i)
        site_lik += root_freqs[i] * clRoot[offset + i];
    total_const_prob += site_lik;
}
...
return total_const_prob / nCat;

// NEW — single character, compensate by × kStates
{
    const double* clRoot = cl_flat.data() + root * stride;
    double site_lik = 0.0;
    for (int i = 0; i < kStates; ++i)
        site_lik += inv_k * clRoot[i];  // root_freqs[i] = inv_k = 1/kStates
    total_const_prob += site_lik;
}
...
return total_const_prob * kStates / nCat;  // × kStates restores the k-fold sum
```

The inner edge-traversal loop is unchanged except that `nChar = 1` makes the
inner `for (int c = 0; c < nChar; ++c)` loop degenerate to a single iteration.

---

## Out of scope for M-170

- `singleton_site_prob_jc`: not used yet (`coding = "variable"`, singleton
  correction is Phase 7). Unchanged.
- `het_constant_site_prob`: F81 path; state-symmetry doesn't hold under
  arbitrary base-frequency rotation. Unchanged.
- `constant_site_prob_mkn`: MkN binary asymmetric model; two patterns differ
  in probability. Unchanged.
- Fusing CSP into `pruning_jc_acrv_persite` for the Gibbs sweep: worthwhile
  follow-up (could eliminate the separate traversal entirely), but more invasive.
  File as M-172 if the 1-char fix doesn't close the gap sufficiently.

---

## Expected impact

| Function | Before | After (est.) |
|----------|--------|-------------|
| `constant_site_prob_jc` | 18.9% | ~4–5% |
| `pruning_jc_acrv_persite` | 68.2% | ~73% (share increases) |
| **Total CPU** | 100% | **~86%** (est. 14% saving) |

Speedup within the function ≈ k_avg / 1 ≈ 4×. Net wall-time saving across the
full MCMC ≈ 14%.

---

## Tests

New test file: `tests/testthat/test-m170-csp-symmetry.R`

1. **Numerical agreement:** Call `constant_site_prob_jc()` (exported) for
   k = 2, 3, 4, 5 with a fixed tree and branch lengths. Compare before/after
   (run against installed reference via `.vtune-lib` or use the prior-build
   snapshot). Actually: check P(constant) ≤ 1, > 0, and consistent with a
   hand-computed value for k=2 on a simple 3-tip tree.

2. **Invariance to kObs:** CSP depends only on k and tree, not on observed
   character data. Verify that two partitions with different observed data but
   same k, tree, and edge lengths give identical CSP.

3. **No MCMC regression:** Add a short fixed-seed MCMC run and compare
   log-posterior trace vs. a stored snapshot (or just verify it doesn't crash
   and converges to the same ballpark).

Existing 4000+ tests cover the remaining correctness burden.

---

## Files changed

| File | Change |
|------|--------|
| `src/ascertainment.cpp` | `constant_site_prob_jc`: nChar=1, stride=kStates, simplified tip init and root, ×kStates in return |
| `tests/testthat/test-m170-csp-symmetry.R` | New regression tests |
| `to-do.md` | M-170: ASSIGNED (D) |
| `completed-tasks.md` | M-170 entry on completion |
| `AGENTS.md` | Performance notes update |
