# M-169: Fix `het_constant_site_prob` per-edge-product bug

**Date:** 2026-03-31
**Priority:** P2
**Assignee:** Agent B

## Problem

`het_constant_site_prob()` (mcmc_likelihood.cpp:1348–1411) computes
P(constant site) using a per-edge-product formula:

```
P(const in state s) = π_s × Π_e P_ss(t_e)
```

This ignores tree topology entirely — the `parent`/`child` args are
accepted but never used. The formula is equivalent to a star tree
(all tips connected to root), not the actual tree.

M-157 fixed this for the **main hot-path** by fusing F81 Felsenstein
pruning of k pseudo-characters into the flat pruning functions
(`pruning_f81_het_acrv_flat`, etc.) via the `outConstProb` parameter.

However, `const_site_prob_for_k()` (line 1522) — called from the
**Gibbs kPrime sweep** (mcmc.cpp:3221) — still dispatches to the buggy
standalone `het_constant_site_prob`. This biases k' inference whenever
Q-het + variable coding + transformational characters are active.

## Scope of impact

- Only affects the Q-heterogeneity (F81) path — JC path uses
  `constant_site_prob_jc()` (correct Felsenstein pruning).
- Only affects the Gibbs kPrime sweep, which calls `const_site_prob_for_k()`
  to compute ascertainment corrections for candidate k' values.
- The main likelihood evaluation (MH proposals, log-posterior) uses
  the fused ascertainment (M-157), which is correct.

## Fix

Replace the body of `het_constant_site_prob()` with proper Felsenstein
pruning of k pseudo-characters over the tree, identical in structure to
the M-157 fused ascertainment code but as a standalone function.

The approach:
1. Allocate a CL buffer for k pseudo-characters (each with k states):
   stride = k × k.
2. Init tip CLs: pseudo-char s has CL = e_s (identity — all tips in
   state s for the constant-state-s pattern).
3. For each (cat × bin × rot) combination:
   a. Build pi vector (F81 frequencies).
   b. Compute mu = 1/(1 - Σπ²).
   c. Traverse tree in postorder (using parent/child):
      `P_ij(t) = π_j(1-d) + δ_ij·d` where `d = exp(-μt)`.
      Hoist Σ_j π_j·cl_j for O(k) per pseudo-char per edge (OPP-1).
   d. At root: sum_s π_s · clRoot[s] for each pseudo-char, accumulate.
4. Return average over all components.

### Performance considerations

This function is called once per unique candidate k' value during the
Gibbs sweep (results cached in `cspCache`). Typically 5–20 unique k'
values. Each call does a full tree traversal with k pseudo-characters —
same cost as the M-157 fused ascertainment. This is fast enough since
it's called O(K_distinct) times per sweep, not O(nChar) times.

## Changes

### `src/mcmc_likelihood.cpp`

1. **Rewrite `het_constant_site_prob()`** (lines 1348–1411):
   Replace per-edge-product with Felsenstein pruning using k pseudo-chars.
   Signature stays the same (static, same args) so `const_site_prob_for_k()`
   needs no changes.

### `tests/testthat/test-fused-ascertainment.R`

2. **Add regression test**: Compare `const_site_prob_for_k()` in het mode
   against a brute-force reference on a small tree (e.g. 4–6 tips),
   verifying the standalone het ascertainment matches the fused result
   from `pruning_f81_het_acrv_flat(..., outConstProb = &p)`.

   The test calls `const_site_prob_for_k()` directly (exported as
   non-exported C++ accessible via `.Call` or through the R-side
   `eval_full_loglik_cpp` path), and compares against the fused result.

   Alternatively: verify that for het + variable coding, the correction
   from the standalone matches the correction applied by the full
   likelihood function (ll_var - ll_none = -nChar * log(1 - p_standalone)).

### `to-do.md`

3. Update M-169 status to DONE.

## Validation

- All existing tests pass (3960+).
- New test confirms het ascertainment is correct.
- Verify the fused-ascertainment test still passes (it uses the main
  path which was already correct).

## NOT in scope

- `het_singleton_site_prob()` — still returns 0.0, correct because
  informative coding is Phase 7.
- Performance optimization of the standalone function — it's called
  O(K_distinct) ≈ 5–20 times per Gibbs sweep, not a bottleneck.
