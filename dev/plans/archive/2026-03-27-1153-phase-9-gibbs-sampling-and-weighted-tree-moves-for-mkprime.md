# Phase 9: Gibbs Sampling and Adaptive Moves for MkPrime

## Motivation

The current MCMC engine uses Metropolis-Hastings (MH) for all parameters.
Two parameters have analytically tractable full conditionals that can be
sampled exactly via Gibbs, eliminating proposal tuning overhead and
achieving acceptance rate 1.0. Additionally, static move weights miss
opportunities to adapt the move mixture based on observed acceptance rates.

## Model recap (relevant parameters)

The prior on k'_i (transformational characters) is a shifted geometric:
```
k'_i | p ~ p * (1-p)^(k'_i - kObs_i),  k'_i >= kObs_i
p ~ Beta(a, b)
```

**p does not appear in the tree likelihood** — it only enters through the
prior on k'. This makes p a textbook conjugate Gibbs target.

k'_i enters the tree likelihood through JC(k'_i) but is discrete with a
geometrically decaying prior, making it amenable to enumeration Gibbs.

---

## Tasks

### M-082: Gibbs update for p (P1)

**Conjugate Gibbs draw.** The full conditional for p given all k'_i is:

```
p | k', ... ~ Beta(a + nTrans, b + Σ(k'_i - kObs_i))
```

No likelihood evaluation needed. Implementation:

1. **C++ (`mcmc.cpp`):** Add moveType 9 = `gibbs_p` to `do_move_impl()`.
   Instead of propose/evaluate/accept, directly:
   - Compute `shape1 = data.kprimeHyperA + nTrans`
   - Compute `shape2 = data.kprimeHyperB + sumU` where
     `sumU = Σ(kPrime[gi] - kObs[gi])` over transformational indices
   - Draw `p_new = R::rbeta(shape1, shape2)`
   - Update `state->p = p_new`
   - Recompute `state->logPrior` (p changed, affects prior on k' and p itself)
   - Return `true` always (Gibbs move, acceptance = 1.0)
   - No likelihood recomputation needed (p doesn't affect likelihood)

2. **R (`.BuildMoves()`):** Change the `p` move from `type = "scale"` to
   `type = "gibbs_p"`. Keep the same weight (1). Remove `scale_p` tuning
   parameter from adaptation (Gibbs needs no tuning).

3. **R (`.kMoveTypes`):** Add `p_gibbs = 9L` mapping.

4. **Acceptance tracking:** Gibbs moves always accept. The
   `accept_counts` / `propose_counts` bookkeeping already works — just
   always increment both. The printed acceptance rate will be 100%,
   which is correct and informative.

5. **Tests:**
   - Unit test: draw many p samples from known k' vector, verify empirical
     distribution matches Beta(a + n, b + sumU)
   - Integration test: run short MCMC with Gibbs p vs. MH p, verify
     comparable posterior (matching marginal distribution of p)
   - Verify that removing scale_p from adaptation doesn't break warmup

**Impact:** Eliminates ~30-40% rejection rate on p proposals. Saves one
`cpp_log_prior()` call per p move (though prior is cheap). Primary value
is correctness and simplicity — one fewer tuning parameter.

---

### M-083: Per-character log-likelihood function (P1)

**Prerequisite for M-084 (Gibbs k'_i).** Add a C++ function that computes
the log-likelihood for a single transformational character under a
specified k value. Currently the code only computes partition-level or
full-tree likelihoods.

1. **C++ (`mcmc_likelihood.cpp`):** New function:
   ```cpp
   double cpp_single_char_log_likelihood(
       const McmcData& data, int globalCharIdx, int kPrime,
       IntegerVector parent, IntegerVector child,
       NumericVector edgeLen,
       double rateLogSd, ClWorkspace* ws);
   ```
   This extracts the tip states for the single character from its partition,
   runs JC(kPrime) pruning for that one column, applies ACRV if
   `rateLogSd > 0`, and applies ascertainment correction (per-site portion).
   If `data.relabel`, adds the relabelling correction for that character.

2. **Declare in `mcmc_state.h`** for use by `mcmc.cpp`.

3. **Tests:** Verify that summing per-character likelihoods for all
   characters in a transformational partition matches the partition
   likelihood from `cpp_partition_log_likelihood()`.

**Note:** This function processes one character at a time, so it's `O(nNode)`
per call. For Gibbs enumeration of ~10 candidate k values per character,
this means ~10 pruning passes per character per Gibbs step — cheap
relative to a full-tree evaluation.

---

### M-084: Gibbs update for k'_i (P2)

**Discrete enumeration Gibbs for each transformational character's k'.**

The full conditional for k'_i = k is:

```
P(k'_i = k | rest) ∝ L_i(data | k, tree, ...) × p × (1-p)^(k - kObs_i) × C(k, kObs_i)
```

where L_i is the per-character likelihood and C is the relabelling
correction. The geometric prior decays exponentially in k, so only a
finite range of k values has non-negligible probability.

1. **C++ (`mcmc.cpp`):** Add moveType 10 = `gibbs_kprime` to
   `do_move_impl()`:
   - Select a random transformational character index
   - Enumerate k from `kObs_i` to `kObs_i + maxEnum` (default maxEnum = 20)
   - For each candidate k, compute:
     `logW[k] = single_char_loglik(i, k) + log(p) + (k - kObs_i) * log(1-p) + relabel_log(k, kObs_i)`
   - Log-sum-exp to normalize, then sample from the discrete distribution
   - Update `state->kPrime[gi] = k_sampled`
   - Recompute `state->logPrior` (k' changed)
   - Update `state->logLik` via partition cache: recompute the affected
     partition and update `partLogLik`
   - Return `true` always

2. **Truncation:** The geometric prior ensures rapid decay. Use a
   dynamic cutoff: stop enumerating when `logW[k] - logW_max < -30`
   (probability ratio < 1e-13). This naturally handles the heavy tail.

3. **CL workspace:** The single-character pruning calls are small (stride = 1 × k
   per call). The existing workspace should handle this easily. However,
   after k' changes, the affected partition's cached log-likelihood must
   be refreshed — call `cpp_partition_log_likelihood()` for the affected
   partition and update `partLogLik`.

4. **R side:** Replace the `int_walk` entry in `.BuildMoves()` with
   `gibbs_kprime`. Same weight (`max(1, 2 * nTrans)`). Remove
   `int_walk_window` from adaptation.

5. **Tests:**
   - Unit: verify sampled k' marginal matches known posterior on a small
     hand-crafted example (2-3 taxa, 1 character)
   - Verify that Gibbs k' produces the same posterior as int_walk on a
     short run (KS test on marginals)
   - Verify partition cache consistency after Gibbs k' update

**Impact:** The current int_walk proposes k' ± 1 with some rejection.
Gibbs jumps directly to the correct conditional distribution, potentially
making larger k' changes in one step and eliminating all rejections.
Especially valuable when the posterior on k' is multimodal or when the
current k' is far from the mode.

---

### M-085: Adaptive move weights during warmup (P3)

**Online adaptation of move weights based on acceptance rates.**

The current move weights are static (set in `.BuildMoves()`, never
updated). During warmup, acceptance rates are tracked and used to tune
proposal widths — but the relative frequency of each move type is fixed.

1. **Design:** At each adaptation interval (every 200 iterations during
   warmup), update move weights proportional to:
   `w_new[m] = w_init[m] × f(accept_rate[m])`
   where f is a function that increases weight for moves with moderate
   acceptance (~25-45%) and decreases weight for moves with very low
   (<5%) or very high (>90%) acceptance. Very high acceptance suggests
   the move is too conservative (small steps); very low suggests it's
   wasteful.

   For Gibbs moves (acceptance = 100%), keep weight constant — they're
   always productive.

2. **Implementation:** In the R adaptation callback
   (`.AdaptProposals()`), after tuning proposal widths, compute
   updated weights and pass them to the next `run_mcmc_batch_cpp()`
   call. The C++ side already accepts `moveWeights` as a parameter.

3. **Post-warmup:** Freeze weights at the end of warmup (same as
   proposal widths are frozen). This maintains detailed balance during
   sampling, since the move mixture doesn't change.

4. **Tests:**
   - Verify weights change during warmup and freeze after
   - Verify that adapted weights improve ESS/iter on a test case
   - Regression test: adapted weights don't break convergence on the
     hyoliths dataset

**Impact:** This is lower priority (P3) because the static weights already
work reasonably well. The main scenario where this helps is when some
move types have persistently low acceptance and waste iterations.

---

## Implementation order

```
M-082 (Gibbs p)  ←  no dependencies, simplest change
    ↓
M-083 (per-char lik) ← standalone utility needed by M-084
    ↓
M-084 (Gibbs k')  ← depends on M-083
    ↓
M-085 (adaptive weights) ← independent, lowest priority
```

M-082 and M-083 can be done in parallel since they're independent.
M-085 is independent of all others.

## Validation strategy

- All tasks include unit tests
- After M-082 + M-084, run a full MCMC comparison:
  - Same dataset, same seed, Gibbs vs. MH for p and k'
  - Compare ESS/iter, ESS/second, posterior marginals
  - This validates that the Gibbs moves produce the correct posterior
- GHA CI for regression testing

## Files touched

| Task | New files | Modified files |
|------|-----------|----------------|
| M-082 | — | `src/mcmc.cpp`, `R/RunMkPrime.R` |
| M-083 | — | `src/mcmc_likelihood.cpp`, `src/mcmc_state.h` |
| M-084 | — | `src/mcmc.cpp`, `R/RunMkPrime.R` |
| M-085 | — | `R/RunMkPrime.R` |
