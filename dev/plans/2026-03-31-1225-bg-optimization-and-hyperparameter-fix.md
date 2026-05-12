# Plan: BG Benchmark Documentation, Gibbs Sweep Optimization, Hyperparameter Fix

**Created:** 2026-03-31 12:25  
**Scope:** AGENTS.md update + 2 optimization tasks  
**Estimated tasks:** M-163 (hyperparameter fix), M-164 (Gibbs sweep optimization)

---

## Execution Order

1. **Update AGENTS.md** — Document BG benchmark results (quick, standalone)
2. **Fix hyperparameter freezing (M-163)** — Root-cause bug fix + slice samplers
3. **Optimize Gibbs kPrime sweep (M-164)** — Tighten early termination
4. **Validate** — Combined benchmark on Sun2018

Order rationale: Fixing hyperparameters first is essential because (a) it's
a correctness bug, (b) properly-functioning α,β affect the BG prior shape
which changes the Gibbs sweep termination profile, and (c) the sweep
optimization should be validated against a correctly-functioning sampler.

---

## Step 1: Update AGENTS.md Performance Notes

Add a new subsection under "Performance notes" documenting:
- BG vs Geometric benchmark results (wall time, ESS/s, per-iteration mixing)
- The kPrime freezing pathology in geometric and how BG fixes it
- The 10× per-iteration cost tradeoff

This is a straightforward documentation edit.

---

## Step 2: Fix Hyperparameter Freezing (M-163)

### Root Cause

**Bug:** `kprime_alpha` and `kprime_beta` are absent from the `targets`
vector in `.AdaptTuning()` (line 3503). When the adaptation loop runs:

```r
target <- targets["kprime_alpha"]  # → NA_real_
adj <- exp(0.5 * (rate - NA))      # → NA
tuning[["scale_kprime_alpha"]] <- tuning[["scale_kprime_alpha"]] * NA  # → NA
```

After 20 proposals (~2 batches), the tuning value is corrupted to NA.
`%||%` doesn't rescue this (NA is not NULL). In C++, `std::exp(NaN * ...)` →
NaN → proposal rejected by positivity check. **All subsequent proposals fail.**

The reported "1.7–2% acceptance" is entirely from the ~20 proposals before
the first adaptation round corrupts the tuning. Post-corruption: 0% acceptance.

### Fix (3 sub-tasks)

#### 2a. Add adaptation targets

In `.AdaptTuning()`, add entries to the `targets` vector:
```r
kprime_alpha = 0.35, kprime_beta = 0.35
```

This enables the standard scale adaptation loop to adjust the proposal
width during warmup. Target 0.35 matches other scale proposals.

#### 2b. Fix initial scale defaults

In `MkPrimeMCMC.R`, reduce the default scale tuning:
```r
scale_kprime_alpha = 0.05,  # was 0.5
scale_kprime_beta = 0.05,   # was 0.5
```

With N=225 characters, the posterior of (α,β) is concentrated (SD ~ 1/√N ≈ 0.07).
A scale of 0.5 on the log scale overshoots massively. Starting at 0.05
gives the adaptation a better starting point, while the adaptation loop
(now enabled by 2a) will fine-tune further.

#### 2c. Add prior-only slice samplers for α and β

**Motivation:** Scale proposals + adaptation may still struggle because
(α,β) are correlated in the posterior (both appear in B(α,β) terms).
Slice samplers don't need tuning and are guaranteed to explore the support.

**Implementation:**

Extend the existing slice sampler infrastructure. The key advantage: these
are **prior-only** parameters — they don't affect the likelihood. The slice
target function only needs:

```cpp
double eval_kprime_hyper_target(McmcData* data, McmcState* state, double beta) {
    return cpp_log_prior(...);  // No likelihood evaluation needed
}
```

This is O(N) beta function evaluations — orders of magnitude cheaper than
the existing slice samplers (which evaluate the full likelihood).

Concrete changes:

1. **`src/mcmc.cpp`:** Add `slice_kprime_alpha_impl` and `slice_kprime_beta_impl`
   functions, following the `slice_scalar_impl` pattern but calling only
   `cpp_log_prior()` (no `cpp_log_likelihood()`). Add new move type codes
   (29 = slice_kprime_alpha, 30 = slice_kprime_beta).

2. **`R/RunMkPrime.R`:**
   - Add `slice_kprime_alpha` and `slice_kprime_beta` moves in `.BuildMoves()`
     when `kPrimePrior == "beta_geometric"`, with appropriate weights
   - Add to `alwaysAcceptTypes` (slice samplers are auto-pinned)
   - Add slice width defaults and adaptation in `.AdaptSliceWidths()`

3. **`R/MkPrimeMCMC.R`:** Add `slice_width_kprime_alpha` and
   `slice_width_kprime_beta` defaults (initial width = 1.0).

4. **Keep scale proposals as complement.** Don't remove moveTypes 27/28;
   they're cheap and provide variety. The slice samplers are the primary
   movers; scale proposals help with mixing when the slice gets stuck.

### Tests for Step 2

- Unit test: `.AdaptTuning()` with kprime_alpha/beta returns non-NA tuning
- Unit test: slice samplers for α and β explore correctly (given fixed k' values)
- Integration test: short BG run produces non-constant α,β traces

---

## Step 3: Optimize Gibbs kPrime Sweep (M-164)

### Current Bottleneck

`LOG_CUTOFF = -57.5` is extremely permissive. A candidate contributing
exp(-57.5) ≈ 1.3e-25 relative probability is below double precision
(~1e-16) and cannot affect the categorical sampling outcome. The BG prior's
slow (logarithmic) tail decay means characters survive to much higher ko
values than necessary, each requiring an expensive tree traversal.

### Optimization: Tighten LOG_CUTOFF

**Change:** `LOG_CUTOFF = -57.5` → `LOG_CUTOFF = -25.0`

**Safety analysis:**

| LOG_CUTOFF | Relative prob | Effect on sampling |
|------------|--------------|-------------------|
| -57.5 | 1.3e-25 | Current; far below double precision |
| -35.0 | 6.3e-16 | Matches double precision; theoretically tight |
| -25.0 | 1.4e-11 | 1e-11 relative error; well below sampling noise |
| -20.0 | 2.1e-9 | Marginal for edge cases with many candidates |

At -25.0, even 50 such candidates contribute 7e-10 total, which is below
uniform RNG precision (~1e-16 per draw). This is safe.

**Expected impact:** For BG with α≈1, β≈1, the per-step weight decline is
~1-3 log-units (likelihood dilution + slow prior decay). Tightening by
32.5 log-units should eliminate ~10-20 unnecessary ko values per character,
translating to ~30-50% fewer tree traversals in the sweep.

### Additional: Prior-ceiling pre-termination check

Before the expensive tree traversal at ko, for each active character,
check if the prior alone exceeds any theoretical possibility of the weight
meeting the cutoff. Specifically:

For a character where the corrected likelihood has been monotonically
decreasing over the last 2 ko steps, and the current weight is already
10 log-units below the cutoff, mark it as terminable with high confidence.

Actually, a simpler and safe optimization: **track the maximum corrected
likelihood seen across all ko values for each character**. If this maximum
is `llMax`, then the best possible weight at ko is `beta * llMax + bgLogPrior[ko]`.
When this drops below `charMaxLogW + LOG_CUTOFF`, terminate early without
computing the actual likelihood. This is safe because the likelihood
generally decreases with k (dilution dominates corrections).

However, this bound is not rigorous (the likelihood could theoretically
increase at higher k if the relabelling correction compensates). So I'll
implement this as an **optimistic pre-filter** that skips individual characters
before the traversal, but never skips the traversal itself (which runs for
all active characters in the partition).

Implementation:
- Add a `charMaxLL` vector tracking the best corrected ll per character
- Before tree traversal at ko, remove characters from the active set where
  `beta * charMaxLL[ti] + bgLogPrior[ko] < charMaxLogW[ti] + LOG_CUTOFF`
- This is cheap (array comparison) and can avoid some sub-matrix columns

### Tests for Step 3

- Regression test: Gibbs sweep produces identical distributions with
  tighter cutoff (compare sampling distributions with -57.5 vs -25.0
  using a fixed seed and large sample)
- Performance benchmark: measure sweep wall time before/after

---

## Step 4: Combined Validation

Run the Sun2018 benchmark (54 taxa, 225 chars, nCat=6) with both fixes
applied. Compare:
- Per-iteration cost (ms/iter) — expect significant reduction from step 3
- α,β mixing (ESS) — expect non-zero ESS from step 2
- kPrime mixing — should remain good (inherited from BG)
- Overall ESS/s — goal is to close the 10× gap while maintaining BG's
  mixing advantages

---

## File Change Summary

| File | Changes |
|------|---------|
| `AGENTS.md` | Add BG benchmark subsection to Performance notes |
| `R/RunMkPrime.R` | Fix `.AdaptTuning()` targets; add slice sampler moves |
| `R/MkPrimeMCMC.R` | Reduce initial scale defaults; add slice width defaults |
| `src/mcmc.cpp` | Add slice_kprime_{alpha,beta}_impl; tighten LOG_CUTOFF; add pre-filter |
| `tests/testthat/` | New test files for adaptation fix + Gibbs cutoff regression |

## Task IDs

- **M-163:** Fix BG hyperparameter (α,β) freezing — adaptation bug + slice samplers
- **M-164:** Optimize Gibbs kPrime sweep for BG prior — tighter early termination
