# M-114: Partial CL for Gibbs Moves under Q-Heterogeneity

## Problem

When `qHeterogeneity = TRUE`, both `gibbs_spr_impl` and
`gibbs_subtree_swap_impl` fall back to `_full` variants that run a
complete Felsenstein pruning (O(nEdge)) per candidate. The partial CL
optimization (O(depth) per candidate) is bypassed, losing the ~3-4×
speedup measured in M-105/M-111.

The root cause: partial CL caches JC/MkN conditional likelihoods, but
Q-het uses F81 transitions with per-component frequency vectors
π(betaBin, rotation). The cached CLs are invalid for other components.

## Approach: Streaming Partial CL

Rather than caching all components simultaneously (memory-prohibitive:
nAcrvCat × nBetaCat × nRot sets of I/F buffers), **stream** over
(betaCat, rotation) components. For each component:

1. Set the CLGroup's F81 parameters (π vector, μ)
2. Run `caching_downpass` with F81 transitions → fills I/F buffers
3. Compute residual CLs (SPR only) with F81 transitions
4. For each candidate, evaluate via partial CL path → **accumulate**
   per-site raw likelihoods into an external array
5. Repeat for next component (CLGroup buffers overwritten)

After all components, convert accumulated per-site likelihoods to
log-likelihoods by dividing by totalComp = nAcrvCat × nBetaCat × nRot
and taking log.

**Memory:** Same as current partial CL (one CLGroup per group, reused
across components). Extra: nCand × nChar doubles for accumulator per
group (~360 KB for Sun2018). Negligible.

**Expected speedup:** Same factor as homogeneous case (≈ nEdge/depth),
because the component count cancels out (both full and partial CL scale
linearly in nComponents). For Sun2018: ~8× per candidate evaluation.

## Implementation Steps

### Step 1: Add F81 transition support to CLGroup

In `gibbs_partial_cl.h`:

- Add `modelType` enum (`JC`, `MKN`, `F81`) and F81 fields to
  `CLGroup`:
  ```cpp
  enum CLModelType { CL_JC, CL_MKN, CL_F81 };
  CLModelType modelType = CL_JC;
  double f81Pi[16];  // frequency vector (max k=16)
  double f81Mu;      // 1 / (1 - Σπ²)
  ```
  Change `isMkN` to be derived: `bool isMkN() const { return modelType == CL_MKN; }`

  **Breaking change note:** `isMkN` is used in `caching_downpass`,
  `evaluate_candidate`, `evaluate_swap_impl`, `compute_residual_cl`,
  and the ascertainment functions. All must be updated.

- Add `f81_transition()` helper:
  ```cpp
  // (P × cl)_i = (1 - e^{-μt}) × dot(π, cl) + e^{-μt} × cl_i
  inline void f81_transition(const double* cl, double* result,
                             int nChar, int kStates,
                             const double* pi, double mu, double t);
  ```

### Step 2: Update existing partial CL functions

Modify the transition dispatch in:
- `caching_downpass`: add F81 branch alongside JC/MKN
- `evaluate_candidate` (applyTransition lambda): add F81 branch
- `compute_residual_cl`: add F81 branch (if it has inline transitions)
- `evaluate_swap_impl`: add F81 branch
- Root frequency: use `f81Pi[s]` when modelType == CL_F81

These are small, localized changes — add a third branch to existing
`if (grp.isMkN) { ... } else { ... }` blocks.

### Step 3: Create accumulating evaluate variants

New functions that accumulate per-site raw likelihoods into an external
array instead of returning a scalar log-likelihood:

```cpp
// Accumulate per-site likelihoods (sum over ACRV cats) into siteLikAccum.
// Does NOT divide by nCat or take log.
static void evaluate_candidate_accum(
    const CLGroup& grp, const TreeNav& topo,
    const ResidualCL& res, const NumericVector& rates,
    int v, int u, int sibNode, double lMerge,
    int a, int b, double lHalfReg, double lPrune,
    double* siteLikAccum);  // [nChar], accumulated into

// Same for swap
static void evaluate_swap_accum(
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int nodeA, int nodeB,
    int pA, int slotA, double lenA,
    const std::vector<int>& pathA,
    const std::vector<int>& pathAIdx,
    double* siteLikAccum);

// Same for const_prob accumulation (ascertainment correction)
static void evaluate_const_prob_accum(..., double* constProbAccum);
static void evaluate_swap_const_prob_accum(..., double* constProbAccum);
```

These are thin wrappers around the existing evaluation logic — same path
propagation, but writing to an external accumulator. The inner loop body
(per-site root likelihood) changes from:

```cpp
siteLikSum[c] += sl;  // internal array, then log(avg)
```

to:

```cpp
siteLikAccum[c] += sl;  // external array, no log yet
```

### Step 4: Build F81 component parameters helper

```cpp
// Compute π and μ for one (betaBin, rotation) component.
// Mirrors the logic in pruning_f81_het_acrv_flat (mcmc_likelihood.cpp:634-653).
static void build_f81_component(
    double betaVal, int rot, int kStates, double baseRL,
    double* pi, double& mu);
```

### Step 5: Q-het streaming wrapper for SPR

In `mcmc.cpp`, modify `gibbs_spr_impl`:

```cpp
if (data->qHeterogeneity) {
  // Instead of: return gibbs_spr_impl_full(data, state, beta);
  // Use streaming partial CL:
  return gibbs_spr_impl_het(data, state, beta);
}
```

New `gibbs_spr_impl_het()`:
- Same prune edge selection and candidate enumeration as current
- Same CLGroup setup, but set `modelType = CL_F81`
- Allocate per-candidate, per-group siteLikAccum arrays
- Compute betaBins from `state->betaScale` and each group's kStates
- Stream over (betaCat, rot):
  - Set each group's f81Pi, f81Mu
  - caching_downpass (overwrites I/F)
  - compute_residual_cl
  - For each candidate: evaluate_candidate_accum
  - For pseudo-groups (ascertainment): same streaming
- After all components: convert accumulators to log-likelihoods
- Sampling and commit: identical to current

### Step 6: Q-het streaming wrapper for subtree swap

Same pattern as Step 5, but for `gibbs_subtree_swap_impl_het()`:
- No residual CL (swap doesn't use it)
- Uses evaluate_swap_accum instead
- Otherwise identical streaming structure

### Step 7: Validation tests

Following the M-105/M-111 validation pattern:

1. **F81 transition unit test**: verify `f81_transition()` against
   known analytic values for small k.

2. **Partial CL vs full agreement (Q-het)**: For a small tree (6-8
   tips) with qHeterogeneity enabled, enumerate all SPR/swap candidates.
   For each, compare the partial CL F81 log-likelihood against
   `compute_full_loglik_at()`. Require agreement to ~1e-12.

3. **End-to-end MCMC with Q-het + Gibbs**: Short run (1000 iter) on a
   small dataset with qHeterogeneity. Verify Gibbs moves fire (nonzero
   acceptance) and logPosterior is reasonable.

### Step 8: Remove _full fallback guard

Delete the `if (data->qHeterogeneity) return _full(...)` guards in
both functions. The partial CL path now handles Q-het natively.

(Keep the `_full` functions for potential future validation use.)

## Complexity

Most of the complexity is in Step 3 (accumulating variants) and Step 5/6
(streaming wrappers). Steps 1-2 are straightforward mechanical changes.
Step 7 reuses the established validation pattern.

The `evaluate_candidate_accum` functions duplicate the path propagation
logic of `evaluate_candidate`. An alternative would be to refactor
`evaluate_candidate` to support both modes (return scalar vs accumulate),
but that adds parameter complexity to an already-long function. The
duplication is acceptable given the functions are in the same header
file and change together.

## Files Modified

| File | Changes |
|------|---------|
| `src/gibbs_partial_cl.h` | F81 transitions, modelType enum, accum variants |
| `src/mcmc.cpp` | `gibbs_spr_impl_het`, `gibbs_subtree_swap_impl_het`, remove fallback |
| `src/mcmc_likelihood.cpp` | `compute_het_bins` → make non-static (shared) or duplicate |
| `tests/testthat/test-gibbs-het.R` (new) | Q-het partial CL validation |

## Risks

- **Correctness**: F81 math must match `pruning_f81_het_acrv_flat`
  exactly. Validated via candidate-by-candidate comparison.
- **Memory**: Per-candidate accumulators scale as O(nCand × max_nChar).
  For 200 cands × 225 chars = ~360 KB. Acceptable.
- **Performance regression for non-Q-het**: Adding modelType branches
  to the hot path is negligible (branch predictor handles it). Verify
  with a quick before/after benchmark.
