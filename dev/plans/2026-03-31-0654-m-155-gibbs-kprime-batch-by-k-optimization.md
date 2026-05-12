# M-155: Gibbs kPrime Sweep — Batch-by-k' Optimization

**Task:** M-155 (P1)  
**Status:** Plan  
**Date:** 2026-03-31  

## Problem

S-PROF round 3 found the Gibbs kPrime sweep (moveType 25) is **10× slower**
than equivalent int_walk coverage: 340 ms vs 33 ms for 99 trans chars on
Sun2018 (54 taxa, nCat=6).

**Root cause:** `gibbs_kprime_sweep_impl()` calls `single_char_loglik_jc()`
in a tight loop — once per character per candidate k' value (~300 calls
total). Each call does a full tree traversal with per-call heap allocation
(`jc1_acrv()` allocates CL + flag vectors on every invocation). The
partition-level flat-buffer pruning used by int_walk batches all characters
in one traversal with pre-allocated workspace, amortizing cost to
0.007–0.07 ms/char.

## Key Structural Insight

Transformational partitions are **grouped by kObs** (see `partition.R`
line 32–47). All characters within a partition share the same kObs. This
means we can directly reuse each partition's `tipStates` matrix at any
candidate k ≥ kObs with **zero tip-data copy**: states 0..kObs-1 are valid
indices for a kStates=k JC model, and missing data (-1) is handled
identically regardless of kStates.

## Solution: Precompute All Per-Site Likelihoods in Batched Passes

Replace ~300 individual tree traversals with ~30-40 batched traversals
(one per distinct (partition, candidate_k) combination), then sample
from precomputed log-weights.

### Algorithm

```
1. Pre-compute shared resources (edge lengths, ACRV rates, prior params)

2. Pre-compute per-site log-likelihoods for all (partition, k) combos:
     For k = kMin_global to kMax_global:          // ~18 distinct k values
       CSP_k = getCSP(k)
       For each trans partition p with kObs_p ≤ k: // 1-3 partitions per k
         pruning_*_persite(tree, p.tipStates, kStates=k, ...)
           → fills siteLL[0..nCharPart-1]
         Apply corrections (ascertainment, relabeling) per-char
         Store: charLL[ti][k - kObs_i] for each char in partition

3. Sample k' for each character (random order):
     For each trans char i:
       Look up precomputed charLL[ti][0..K_PRECOMP-1]
       Compute log-weights = β × logLik + logPrior
       Early terminate when weight drops below peak + LOG_CUTOFF
       If not terminated within K_PRECOMP: fall back to single_char_loglik_jc
       Sample from categorical → state->kPrime[gi]

4. Rebuild partition likelihoods + prior (unchanged from current)
```

### Expected Performance (Sun2018, 54 taxa, 99 trans, nCat=6)

| What                    | Count | Est. cost each | Total   |
|-------------------------|-------|----------------|---------|
| Batched pruning calls   | ~40   | ~0.3 ms        | ~12 ms  |
| CSP cache computation   | ~18   | ~0.1 ms        | ~2 ms   |
| Sampling (99 chars)     | 99    | ~0.001 ms      | ~0.1 ms |
| Final partition rebuild | 5     | ~0.4 ms        | ~2 ms   |
| **Total**               |       |                | **~16 ms** |

Compare to current 340 ms → **~20× speedup** (conservative estimate; actual
may be higher due to cache benefits).

## Implementation Steps

### Step 1: Per-site pruning functions (`src/mcmc_likelihood.cpp`)

Create two new static functions by copying the flat-buffer variants and
modifying only the root-extraction / return logic:

**`pruning_jc_acrv_persite()`**
- Signature: same as `pruning_jc_acrv_flat()` plus `double* siteLL` output
  array (length nChar).
- Returns `void` (no summed log-likelihood needed).
- Identical tree traversal (tip init, category loop, edge loop, OPP-1
  JC symmetry).
- Root extraction: instead of accumulating into `logLik`, fills
  `siteLL[c] = log(site_lik_sum[c] / nCat)` for each character.
  Returns `R_NegInf` for non-positive site likelihoods.

**`pruning_f81_het_acrv_persite()`**
- Same pattern for the Het/F81 path.
- Root extraction fills `siteLL[c] = log(site_lik_sum[c] / totalComp)`.

Both functions share 95%+ code with their `_flat` counterparts; the only
change is the last ~10 lines.

### Step 2: Add Gibbs workspace to McmcState (`src/mcmc_state.h`)

Add a persistent `ClWorkspace gibbsWs` field to `McmcState`:

```cpp
struct McmcState {
  // ... existing fields ...
  ClWorkspace gibbsWs;   // Dedicated workspace for Gibbs kPrime sweep
};
```

Allocated once on first use (lazy init in the sweep function). Sized for
the max stride across all (trans partition, candidate k) combinations:

```
maxGibbsStride = max over trans partitions p of (p.nChar × (max_kObs + K_PRECOMP))
```

For Sun2018: max(68×22, 23×23, 8×24) = max(1496, 529, 192) = 1496.
Buffer: 108 × 1496 × 8 ≈ 1.3 MB. Negligible.

### Step 3: Build lookup maps (`src/mcmc.cpp`, in sweep function)

Build once at top of sweep:

- **`globalToTransIdx[nChar]`**: maps global char index → position in
  `transIdxGlobal` (for indexing into `charLL`). -1 for non-trans chars.
  (Same spirit as existing `charToLocalIdx`, but for trans indexing.)

- **`transPartitions`**: list of (partIdx, kObs, nChar) for trans partitions
  only, sorted by kObs. Cache from `data->parts` on first call.

### Step 4: Rewrite `gibbs_kprime_sweep_impl()` (`src/mcmc.cpp`)

Replace the inner loop (lines 3118-3174) with:

**Phase A: Determine k range**
```cpp
const int K_PRECOMP = 20;  // covers ~99.99% of characters
int kMinGlobal = INT_MAX, kMaxObs = 0;
for (trans partitions) { kMinGlobal = min(kObs_p); kMaxObs = max(kObs_p); }
int kMaxGlobal = kMaxObs + K_PRECOMP;
```

**Phase B: Allocate charLL array**
```cpp
// charLL[ti * K_PRECOMP + k_offset] = corrected log-lik for char ti at k = kObs_i + k_offset
std::vector<double> charLL(nTrans * K_PRECOMP, R_NegInf);
```

**Phase C: Ensure Gibbs workspace fits**
```cpp
int maxGibbsStride = compute_max_gibbs_stride(data, kMaxGlobal);
if (!state->gibbsWs.fits(maxNode, maxGibbsStride))
  state->gibbsWs.allocate(maxNode, maxGibbsStride);
```

**Phase D: Batched likelihood precomputation**
```cpp
std::vector<double> siteLL(maxNCharPart);  // temp output buffer

for (int k = kMinGlobal; k <= kMaxGlobal; ++k) {
  double csp = getCSP(k);
  double logAscCorr = (coding != 0) ? -std::log(1.0 - csp) : 0.0;

  for (each trans partition p with kObs_p <= k) {
    int k_offset = k - kObs_p;
    if (k_offset >= K_PRECOMP) continue;

    // Run batched per-site pruning (JC or Het)
    if (useHet) {
      compute_het_bins(betaScale, k, nBC, hetBins);
      pruning_f81_het_acrv_persite(..., p.tipStates, k, ...,
        gibbsWs.buf, gibbsWs.init, gibbsWs.strideMax, siteLL.data());
    } else {
      pruning_jc_acrv_persite(..., p.tipStates, k, acrvRates,
        gibbsWs.buf, gibbsWs.init, gibbsWs.strideMax, siteLL.data());
    }

    // Store corrected log-likelihoods
    for (int c = 0; c < p.nChar; ++c) {
      int gi = p.globalCharIdx[c];
      int ti = globalToTransIdx[gi];
      double ll = siteLL[c] + logAscCorr;
      if (data->relabel)
        ll += mk_prime_relabel_log(k, kObs_p);
      charLL[ti * K_PRECOMP + k_offset] = ll;
    }
  }
}
```

**Phase E: Sampling (mostly unchanged)**
```cpp
// Random permutation (same as current)
for (int si = 0; si < nTrans; ++si) {
  int ti = perm[si];
  int gi = data->transIdxGlobal[ti];
  int kObs_i = data->kObs[gi];

  double logW[K_PRECOMP + K_FALLBACK_MAX];
  int nCand = 0;
  double maxLogW = R_NegInf;

  // Use precomputed likelihoods
  for (int ko = 0; ko < K_PRECOMP; ++ko) {
    int k = kObs_i + ko;
    double logLik_k = charLL[ti * K_PRECOMP + ko];
    if (!R_FINITE(logLik_k)) break;  // R_NegInf → no more precomputed

    double logPrior_k = /* same as current */;
    double w = beta * logLik_k + logPrior_k;
    logW[nCand++] = w;
    if (w > maxLogW) maxLogW = w;
    if (w < maxLogW + LOG_CUTOFF) break;
  }

  // Fallback: if not terminated, use per-char evaluation (rare)
  if (nCand == K_PRECOMP && logW[nCand-1] >= maxLogW + LOG_CUTOFF) {
    for (int ko = K_PRECOMP; ko < K_MAX_CAND; ++ko) {
      int k = kObs_i + ko;
      double csp = getCSP(k);
      double logLik_k = single_char_loglik_jc(...);  // fallback
      // ... same weight computation ...
    }
  }

  // Sample from categorical (same as current)
  // ...
  state->kPrime[gi] = kObs_i + chosen;
}
```

**Phase F: Rebuild** (unchanged from current lines 3176-3199)

### Step 5: Tests (`tests/testthat/`)

1. **Numerical equivalence test:** For a small test case (e.g., 10 taxa,
   15 trans chars), verify that the batched Gibbs sweep produces the
   same per-character log-likelihoods as `single_char_loglik_jc()` for
   each (char, k) pair. Use a fixed seed and compare log-weight vectors.

2. **Per-site pruning unit tests:** Verify `pruning_jc_acrv_persite()`
   produces per-site values that sum to match `pruning_jc_acrv_flat()`.
   Test with kStates ∈ {2, 3, 5}, nCat ∈ {1, 4, 6}, including missing
   data and boundary cases.

3. **Full MCMC integration test:** Run a short MCMC with Gibbs sweep
   enabled and verify logLik drift check passes (existing consistency
   infrastructure).

4. **Het path test:** If Het is enabled, verify per-site Het pruning
   matches single-char Het evaluation.

### Step 6: Cleanup

- Keep `jc1_acrv()` and `single_char_loglik_jc()` for the fallback path
  and potential use by other code. Consider removing in a follow-up task
  if dead code analysis confirms they're only used by Gibbs.
- Remove the `charToLocalIdx` rebuild from the sweep (it was needed for
  the per-char approach; the batched approach uses `globalToTransIdx` and
  partition's `globalCharIdx` directly).

## Risks and Mitigations

| Risk | Mitigation |
|------|------------|
| Numerical divergence between per-site and summed pruning | Step 5.2: per-site values must sum to match `_flat` return value exactly (same FP operations) |
| K_PRECOMP too small for some datasets | Fallback path (Phase E) handles overflow transparently. K_PRECOMP=20 is very conservative. |
| Gibbs workspace memory for very large datasets | Max ~10 MB even for 200 taxa × 500 chars × kMax=30. Negligible. |
| Het path adds complexity | Structure mirrors JC path exactly. Same code pattern with `pruning_f81_het_acrv_persite`. |

## Files Modified

| File | Changes |
|------|---------|
| `src/mcmc_likelihood.cpp` | Add `pruning_jc_acrv_persite()`, `pruning_f81_het_acrv_persite()` |
| `src/mcmc_state.h` | Add `ClWorkspace gibbsWs` to `McmcState` |
| `src/mcmc.cpp` | Rewrite `gibbs_kprime_sweep_impl()` |
| `tests/testthat/test-gibbs-kprime.R` | Numerical equivalence + per-site pruning tests |

Estimated: ~350–400 lines changed/added.
