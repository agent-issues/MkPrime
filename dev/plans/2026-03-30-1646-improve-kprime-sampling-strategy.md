# Improve kPrime Sampling Strategy

**Created:** 2026-03-30
**Status:** IMPLEMENTED (M-154)

## Problem

kPrime (k'_i, the true number of character states for each transformational
character) mixes poorly. The current proposal is a **univariate integer random
walk** (±window, typically window=1) applied to one randomly-chosen character
per proposal. This has several compounding weaknesses:

1. **Random walk inefficiency:** Moving from k'=kObs to k'=kObs+3 requires
   three consecutive accepted +1 steps. Each step can be rejected, so reaching
   distant values takes many iterations.

2. **Correlation through the hyperparameter p:** Under the geometric prior,
   all k' values are coupled through p. When many characters have k' = kObs,
   p is pulled high (favoring k' ≈ kObs), which makes it even harder for
   individual characters to explore higher k' values. Moving between the
   "all low k'" and "some high k'" collective modes requires many characters
   to shift simultaneously — exponentially unlikely under one-at-a-time proposals.

3. **Low weight relative to other moves:** kPrime proposals get weight
   `2 × nTrans`, but with 40+ transformational characters, each individual
   character is only proposed ~2 times per MCMC cycle. With a modest
   acceptance rate, effective updates per character per cycle are very low.

## Proposed Solution: Two New Move Types

### A. Gibbs kPrime Sweep (Primary — highest impact)

**Concept:** Analogous to the existing block Gibbs branch sweep (moveType 15),
sample each k'_i from its exact **full conditional distribution** in a single
random-order scan of all transformational characters.

**Full conditional for character i:**

```
P(k'_i = v | rest) ∝ exp(β × [loglik_i(v) + relabel_i(v)] + prior_i(v))
```

where:
- `loglik_i(v)` = Felsenstein log-likelihood for character i under JC(v),
  with ascertainment correction
- `relabel_i(v)` = lgamma(v+1) - lgamma(v-kObs_i+1)
- `prior_i(v)` = log(p) + (v-kObs_i)·log(1-p) [geometric]
  or v·log(c) - log(v) - log(-log(1-c)) [logseries]
- `β` = chain temperature (for parallel tempering compatibility)

The likelihood+relabelling terms are tempered; the prior is not (matching the
existing MH convention).

**Enumeration range:** For each character, enumerate v ∈ {kObs_i, ...,
kObs_i + K_max}. Use adaptive K_max: start from kObs_i and stop when the
log-weight drops below `max_weight - 25` (≈ exp(-25) ≈ 10^{-11} relative
probability). Cap at a hard maximum (e.g., 30) for safety.

**Cost per sweep:** O(nTrans × K_mean × nNode × k_mean) where K_mean is
the average number of candidates evaluated (typically 5–15) and k_mean is
the average k' value. For hyoliths (40 trans chars, 107 nodes, k≈3–5,
K_mean≈10): ~40 × 10 × 107 × 5 ≈ 214,000 operations — roughly equivalent
to one full partition likelihood evaluation.

**Why this works:**
- **Always accepts** (Gibbs update) — no wasted proposals
- **Can jump directly** from kObs to kObs+5 in one step
- **Full sweep** decorrelates all k' values in one "move"
- Followed by the existing **Gibbs p update**, the full (k', p) block is refreshed

**Key advantage over the integer walk:** The integer walk at window=1 is
essentially a Markov chain that must walk between modes one step at a time.
The Gibbs sampler samples directly from the target distribution, so it
converges in one step regardless of how far the mode is.

### B. Block Uniform Shift (Secondary — complements Gibbs sweep)

**Concept:** Propose to shift ALL transformational characters by the same
integer δ, addressing the collective mode between "all k' ≈ kObs" and
"all k' somewhat above kObs."

**Proposal:**
1. Draw δ ~ Uniform({-W, ..., -1, +1, ..., +W})  (exclude 0)
2. Compute k'_new_i = k'_old_i + δ for all transformational characters
3. If any k'_new_i < kObs_i, reject immediately (O(1) check: just test
   whether δ < -min(k'_i - kObs_i))
4. Evaluate full likelihood (all transformational partitions recomputed)
5. Prior ratio: nTrans × δ × log(1-p) [geometric]; for logseries,
   sum of individual prior ratios
6. MH accept/reject

**Cost:** One full likelihood evaluation (comparable to a topology move).

**Why include this in addition to Gibbs sweep?**
The Gibbs sweep samples each character conditionally on the others, which is
excellent for marginal exploration but doesn't address collective modes. If
the posterior has ridges where many characters are simultaneously high or low
(mediated by p), the block shift moves along these ridges. The Gibbs sweep
+ block shift combination handles both marginal and collective mixing.

**Hastings ratio:** Symmetric (δ and -δ equally likely), but boundary
rejections create implicit asymmetry. The boundary check means the effective
proposal is truncated: for δ < 0, the proposal is only valid if
min(k' - kObs) ≥ |δ|. The Hastings ratio remains 0 because the reverse
move (from the proposed state with -δ) faces the same boundary check at the
same min-gap (all characters shifted by the same amount, so the minimum gap
is unchanged). **Symmetric and correct.**

## Implementation Plan

### Step 1: Per-character likelihood function

New C++ function in `mcmc_likelihood.cpp`:

```cpp
double cpp_single_char_loglik(
    const McmcData& data, int partIdx, int localCharIdx,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, int kCandidate,
    double rateLogSd, double betaScale,
    ClWorkspace* ws);
```

This runs Felsenstein pruning on a single column of the partition's tipStates
matrix at JC(kCandidate). Includes ascertainment correction (constant site
probability depends only on k and tree, so can be precomputed per k and
shared across characters with the same kObs).

**Optimization opportunity:** Precompute the constant-site probability for
each candidate k' value once per partition per sweep, rather than once per
character. Characters in the same partition share kObs, so they share the
ascertainment correction at each candidate k'.

### Step 2: Gibbs kPrime sweep implementation

New C++ function `gibbs_kprime_sweep_impl()` as a standalone (like
`block_gibbs_branch_sweep_impl`), not inline in `do_move_impl()`.

Pseudocode:
```
gibbs_kprime_sweep_impl(data, state, beta):
  edgeLen = treeLength * relBrLengths
  
  // Random-order scan
  perm = random permutation of transIdxGlobal
  
  for each char i in perm:
    partIdx = charToPartition[i]
    localIdx = charToLocalIdx[i]  // index within partition
    kObs_i = data.kObs[i]
    
    // Precompute ascertainment correction per candidate k
    // (can be cached across characters in same partition)
    
    logWeights = []
    maxW = -Inf
    for v in {kObs_i, kObs_i+1, ..., kObs_i+K_max}:
      ll = single_char_felsenstein(i, v, edgeLen, ...)
      if (data.relabel):
        ll += lgamma(v+1) - lgamma(v - kObs_i + 1)
      if (coding != 0):
        ll -= log(1 - constant_site_prob(v, edgeLen, ...))
      
      // Prior contribution for this character
      if (geometric):
        lp = log(p) + (v - kObs_i) * log1p(-p)
      else:  // logseries
        lp = v * log(c) - log(v) - log(-log(1-c))
      
      w = beta * ll + lp
      logWeights.push_back(w)
      maxW = max(maxW, w)
      
      // Early termination
      if (v > kObs_i + 2 && w < maxW - 25): break
    
    // Sample from categorical distribution
    // (log-sum-exp for numerical stability)
    probs = exp(logWeights - maxW)
    probs /= sum(probs)
    state->kPrime[i] = kObs_i + sample(probs)
  
  // Recompute full likelihood and prior (many characters changed)
  state->logLik = compute_full_loglik(...)
  state->logPrior = cpp_log_prior(...)
  state->partLogLik.clear()
  state->nodeCL.valid = false
  
  return true  // always accepts (Gibbs)
```

### Step 3: Block shift implementation

New case in `do_move_impl()` (or standalone function):

```
block_kprime_shift_impl(data, state, beta, window):
  nTrans = data.transIdxGlobal.size()
  
  // Draw delta, excluding 0
  delta = random_int(-window, window)  // excluding 0
  
  // Boundary check
  minGap = INT_MAX
  for each trans char i:
    gap = state->kPrime[i] - data->kObs[i]
    minGap = min(minGap, gap)
  if (delta < -minGap): return false  // reject
  
  // Save old values for rollback
  oldKPrime = clone(state->kPrime)  // or just snapshot trans chars
  
  // Apply shift
  for each trans char i:
    state->kPrime[i] += delta
  
  // Evaluate (full likelihood recompute needed — all trans partitions change)
  newLogLik = cpp_log_likelihood(...)
  newLogPrior = cpp_log_prior(...)
  
  // MH acceptance (Hastings = 0, symmetric)
  logAlpha = beta * (newLogLik - state->logLik) + (newLogPrior - state->logPrior)
  if (log(U) < logAlpha):
    accept, update state
  else:
    rollback kPrime from snapshot
```

### Step 4: R-side move registration

In `.BuildMoves()`:
```r
if (nTrans > 0) {
  kPrimeMoves <- list(
    # Existing integer walk (reduced weight)
    list(name = "kPrime", type = "int_walk", target = "kPrime",
         weight = max(1, nTrans), dim = 1L),
    # NEW: Gibbs kPrime sweep
    list(name = "gibbs_kPrime", type = "gibbs_kprime_sweep",
         target = "kPrime", weight = nTrans, dim = nTrans),
    # NEW: Block uniform shift
    list(name = "block_kPrime", type = "block_kprime_shift",
         target = "kPrime", weight = 2, dim = 1L)
  )
  # ... existing p Gibbs move
}
```

### Step 5: Move type codes and dispatch

- Assign new integer codes (25 = gibbs_kprime_sweep, 26 = block_kprime_shift)
- Add cases to `do_move_impl()` dispatcher
- Add to `.kMoveTypes` in `RunMkPrime.R`

### Step 6: Acceptance tracking and adaptation

- Gibbs kPrime sweep: acceptance is always 1 (like Gibbs p and Gibbs branch).
  Track wall time for adaptive scheduler.
- Block shift: track acceptance rate for MH tuning. Window parameter adapts
  during warmup (like int_walk_window).

### Step 7: Testing

- **Correctness test:** On a small dataset (5 taxa, 3 characters), verify that
  the Gibbs sweep produces the correct marginal distribution of k' by
  comparing to brute-force enumeration.
- **Mixing test:** Compare ESS/second for kPrime parameters between baseline
  (int_walk only) and new moves on hyoliths dataset.
- **Reversibility test for block shift:** Verify that the MH ratio is correct
  by checking detailed balance empirically (ratio of forward/reverse transition
  probabilities matches the acceptance formula).

### Step 8: charToLocalIdx mapping

The existing `charToPartition` map gives the partition index for each global
character index. We also need the **local index within the partition**
(i.e., which column of `part.tipStates` corresponds to this character).
Add `charToLocalIdx` to `McmcData` during `prepare_mcmc_data()`.

## Priority and Phasing

1. **Gibbs kPrime sweep** (Steps 1–2, 4–5, 7–8): This is the primary
   improvement and should be implemented first. Expected to dramatically
   improve kPrime ESS.

2. **Block uniform shift** (Step 3, 4–5, 6–7): Implement after verifying
   the Gibbs sweep helps. May be unnecessary if the Gibbs sweep + existing
   Gibbs p fully resolves mixing.

3. **Keep existing integer walk** at reduced weight as a fallback/complement.

## Risks and Mitigations

- **Cost concern:** The Gibbs sweep is O(nTrans × K × nNode), which could
  be expensive for datasets with many transformational characters (>100) or
  large trees (>100 tips). Mitigation: the adaptive scheduler will naturally
  down-weight the sweep if its ESS/second is poor. Also, K_max is adaptive
  (early termination when weights are negligible).

- **ACRV interaction:** When ACRV is active, each candidate k' must be
  evaluated across all rate categories. This multiplies the cost by nCat
  (typically 4). Still manageable for morphological datasets.

- **Node CL cache invalidation:** The Gibbs sweep changes many kPrime values
  simultaneously. The node CL cache stores conditional likelihoods that depend
  on kPrime (through the P(t) matrices). Must invalidate after the sweep.
  This is already handled by `.valid = false`.

## Expected Impact

For the hyoliths dataset (40 transformational characters):
- Current: ~2 int_walk proposals per character per cycle, ~30-50% acceptance
  → ~0.6–1.0 effective updates per character per cycle
- Gibbs sweep: 1 sweep per invocation, all 40 characters sampled from full
  conditional → 40 effective updates per invocation. Even if the sweep costs
  10× more wall time than a single int_walk proposal, the ESS/second
  improvement should be substantial (order of magnitude or more).
