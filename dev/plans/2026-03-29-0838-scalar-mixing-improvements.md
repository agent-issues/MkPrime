# Plan: Scalar Parameter Mixing Improvements

**Problem:** On Sun2018 (54 taxa, 225 chars), `rate_loss` and `rate_neo`
have ESS ~40 after 280k iterations (PSRF ~1.9). Root cause: initial move
weights give them <0.4% of proposals (1.5/400 and 1.0/400), starved by
kPrime (49%) and branch_lengths (26%).

**Four fixes, in dependency order:**

---

## 1. Scalar parameter weight floor in `.BuildMoves()`

**File:** `R/RunMkPrime.R` — `.BuildMoves()`

After constructing all moves, apply a post-hoc weight floor: each scalar
parameter move (dim=1, not topology/branch) gets at least
`scalarFloorFraction` of the total weight budget.

```
scalarFloorFraction = 0.02  (2% of total)
```

For Sun2018 (total weight ~400): each scalar gets at least 8 (up from 1–1.5).
This guarantees ~2% of proposals per scalar parameter.

**Implementation:**
- Identify "scalar" moves: those with `dim == 1L` and `type %in%
  c("scale", "int_walk", "gibbs_p")`. Exclude topology moves
  (nni/spr/tbr/gibbs_spr/etc.) which have dim=1 but are structural.
- Compute `totalWeight = sum(all move weights)`.
- For each scalar move: `weight <- max(weight, totalWeight * 0.02)`.
- This inflates total weight slightly; that's fine — it just means slightly
  more proposals per thinning cycle, not a normalisation issue (the C++
  loop uses raw weights for weighted sampling).

**Scalar moves affected:** `tree_length`, `rate_loss`, `rate_neo`,
`rate_log_sd`, `beta_scale`, `p`, `kPrime`.

Wait — `kPrime` already has weight ~198 (2×nTrans). The floor wouldn't
touch it. And `tree_length` has weight 1, which would go to 8. That's
reasonable: tree_length is important but cheap.

Actually, let me reconsider which moves get the floor. The problem is
specifically that `rate_loss` (1.5) and `rate_neo` (1.0) are swamped.
The fix should apply broadly to all scalar moves so future parameters
get the same protection, but the actual impact will be on the low-weight
ones.

**Edge case:** `kPrime` with `weight = 2 * nTrans = 198` — floor of
`400 * 0.02 = 8` doesn't touch it. Good.

**Edge case:** Very small datasets (5 taxa, 10 chars) — total weight
might be ~30. Floor = 0.6, which is less than current 1.0–1.5 for
rate_loss. Floor is a no-op. Good.

---

## 2. Joint (rate_loss, rate_neo) scale proposal

**Files:** `src/mcmc.cpp`, `R/RunMkPrime.R`

Add move type 18: `neo_joint_scale`. Proposes a single multiplicative
perturbation applied to **both** `rate_loss` and `rate_neo`
simultaneously.

**Motivation:** These parameters jointly control neomorphic character
likelihoods. If they're positively correlated in the posterior
(likely — both scale neomorphic rates), a joint move along the
correlation ridge is much more efficient than two independent moves.

**C++ implementation (case 18 in do_move_impl):**
```cpp
case 18: { // neo_joint_scale — scale rate_loss and rate_neo together
  double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
  state->rateLoss *= mult;
  state->rateNeo  *= mult;
  logHastings = 2.0 * std::log(mult);  // Jacobian for 2 parameters
  break;
}
```

**R-side wiring:**
- Add to `.kMoveTypes`: `neo_joint_scale = 18L`
- Add to `.BuildMoves()` when `hasNeo`:
  ```r
  list(name = "neo_joint", type = "neo_joint_scale",
       target = c("rate_loss", "rate_neo"),
       weight = <floor-adjusted>, dim = 2L)
  ```
  Weight: same as individual moves after floor, or `totalWeight * 0.02`.
- Add to `.BuildScaleTuningMatrix()`: map `neo_joint` to
  `tun$scale_rate_loss` (shared tuning parameter; both are on same
  scale).
- Add to `.AdaptTuning()`: target acceptance 0.35 (same as rate_loss),
  tuning key `scale_neo_joint` (new).
- Add to rollback in `do_move_impl`: restore both `rateLoss` and
  `rateNeo` from `oldRL`/`oldRN` (already saved at top of function).
- Partial likelihood cache: affects same partitions as case 1/3
  (neoPartIndices only).

**Rollback:** Already handled — `oldRL` and `oldRN` are saved at
the top of `do_move_impl` and restored in the rejection block.

**Partial cache:** Add case 18 alongside cases 1, 3 in the
cache-aware likelihood switch.

---

## 3. Scale tuning budget with total move weight

**File:** `R/RunMkPrime.R` — tuning phase setup

**Problem:** Default `tuningBudget = 10000`. With `thin = auto ≈ 60`
(Sun2018), each tuning window gets ~167 samples. But the bandit runs
4 candidates per round (current + 3 perturbations), so each candidate
gets ~42 samples. Computing min-ESS from 42 samples is unreliable,
especially for parameters proposed in only ~1% of iterations.

**Fix:** After `.BuildMoves()` resolves `thin`, scale tuning budget:

```r
# In the Warmup → Tuning transition:
effectiveTuningBudget <- max(
  mcmc$tuningBudget,
  resolvedThin * 100L * (1L + 3L)  # 100 samples × 4 candidates
)
```

For Sun2018: `max(10000, 60 * 100 * 4) = max(10000, 24000) = 24000`.
Each candidate window gets ~100 samples — enough for meaningful ESS.

**Where:** The scaling happens at the point where tuning phase begins
(around line 672 in RunMkPrime.R), not in `MkPrimeMCMC()` — because
`thin` isn't resolved until `.BuildMoves()` runs.

The user-specified `tuningBudget` acts as the explicit floor; auto-
scaling can only increase it. If the user set `tuningBudget = 5000`
explicitly, we should respect that. So: only auto-scale when
`tuningBudget` was left at its default (10000).

Actually, simpler: just always take the max. If the user explicitly
set `tuningBudget = 50000`, the auto-scale for Sun2018 (24000) won't
touch it. If they set 5000, we still honour the auto-scale because
5000 is genuinely too few for this dataset. The user can set
`autoTune = FALSE` to skip tuning entirely if they don't want it.

---

## 4. Slice sampling for scalar parameters

**Files:** `src/mcmc.cpp`, `R/RunMkPrime.R`

### Where slice sampling helps

Slice sampling is most effective for parameters where:
- The optimal MH step size is hard to determine
- The conditional posterior is unimodal but may have varying curvature
- The parameter is expensive to tune but cheap to evaluate

**Best candidates in MkPrime:**

| Parameter | Why slice helps | Current proposal |
|-----------|----------------|-----------------|
| `rate_loss` | Only affects neo partitions (cheap eval); hard to tune; often poorly mixed | scale (MH) |
| `rate_neo` | Same as rate_loss | scale (MH) |
| `rate_log_sd` | Affects all partitions but is a single scalar; curvature varies with data | scale (MH) |
| `tree_length` | Affects all partitions; already mixes OK but could benefit | scale (MH) |
| `beta_scale` | Only when qHeterogeneity=TRUE; similar profile to rate_log_sd | scale (MH) |

**Not good candidates:**
- `kPrime` — discrete (integer), not suitable for standard slice
- Topology moves — discrete, complex proposal mechanism
- Branch lengths — high-dimensional simplex, slice not applicable
- `p` — already Gibbs-sampled (conjugate)

### Implementation: univariate stepping-out slice sampler

Add move type 19: `slice_scalar`. The slice sampler replaces the
MH accept/reject with an internal loop that always produces an
accepted sample.

**C++ implementation — new function:**

```cpp
// Slice-sample a single scalar parameter.
// Returns true (always "accepted" — slice sampling is exact).
static bool slice_scalar_impl(
    McmcData* data, McmcState* state,
    int paramIdx,   // 0=treeLength, 1=rateLoss, 2=rateLogSd, 3=rateNeo, 16=betaScale
    double width,   // initial bracket width (tunable)
    double beta,    // temperature
    int maxSteps    // max stepping-out steps (default 10)
) {
  // 1. Get current value and log-target
  double x0 = get_scalar(state, paramIdx);
  double logY0 = beta * state->logLik + state->logPrior;
  
  // 2. Draw slice height
  double logZ = logY0 + std::log(R::unif_rand());  // log(y * U)
  
  // 3. Stepping out: find bracket [L, R]
  double L = x0 - width * R::unif_rand();
  double R = L + width;
  // Enforce positivity (all scalar params are > 0)
  if (L <= 0.0) L = 1e-12;
  
  for (int j = 0; j < maxSteps; ++j) {
    set_scalar(state, paramIdx, L);
    double logTargetL = eval_log_target(data, state, paramIdx, beta);
    if (logTargetL <= logZ) break;
    L = std::max(L - width, 1e-12);
  }
  for (int j = 0; j < maxSteps; ++j) {
    set_scalar(state, paramIdx, R);
    double logTargetR = eval_log_target(data, state, paramIdx, beta);
    if (logTargetR <= logZ) break;
    R += width;
  }
  
  // 4. Shrink in: sample uniformly from [L, R], shrink on rejection
  for (int iter = 0; iter < 100; ++iter) {
    double x1 = L + R::unif_rand() * (R - L);
    if (x1 <= 0.0) x1 = 1e-12;
    set_scalar(state, paramIdx, x1);
    double logTarget1 = eval_log_target(data, state, paramIdx, beta);
    if (logTarget1 >= logZ) {
      // Accept: state already updated
      state->logLik = ...; state->logPrior = ...;
      return true;
    }
    // Shrink bracket
    if (x1 < x0) L = x1; else R = x1;
  }
  
  // Fallback: restore original (shouldn't happen in practice)
  set_scalar(state, paramIdx, x0);
  return false;
}
```

**Helper functions needed:**
- `get_scalar(state, paramIdx)` / `set_scalar(state, paramIdx, val)`:
  switch on paramIdx to read/write the appropriate field
- `eval_log_target(data, state, paramIdx, beta)`: compute
  `beta * logLik + logPrior` for the current state. Use partial cache
  when available (rate_loss/rate_neo → neo partitions only).

**Key design decisions:**

1. **Slice width as tuning parameter.** Initial width = current MH
   scale tuning value. Adaptive: track accepted slice widths and
   update (wider bracket → more evals but better exploration;
   narrower → fewer evals but may miss modes). Or just use a fixed
   generous width and let stepping-out handle it.

2. **Partial likelihood evaluation.** For `rate_loss` and `rate_neo`,
   only neo partitions change — reuse existing partial cache logic.
   This makes each slice evaluation cheap (neo chars only, not full
   tree). For `rate_log_sd`, all partitions change.

3. **Multiple evaluations per "move".** Slice sampling typically
   requires 3–8 likelihood evaluations per update (stepping out +
   shrinking). The `moveTimeNs` counter captures total wall time,
   so the adaptive weight scheduler naturally accounts for this cost.
   The `acceptance` counter records 1.0 (always accepts), but the
   ESS-per-second metric in the tuning phase handles the tradeoff
   correctly.

4. **Temperature handling.** The slice height is drawn from
   `beta * logLik + logPrior`, so the tempered distribution is
   correctly sampled. This preserves detailed balance for parallel
   tempering.

**R-side wiring:**
- Add to `.kMoveTypes`: `slice_scalar = 19L`
- New move entries in `.BuildMoves()`:
  ```r
  list(name = "slice_rate_loss", type = "slice_scalar",
       target = "rate_loss", weight = <floor>, dim = 1L,
       paramIdx = 1L)
  list(name = "slice_rate_neo", type = "slice_scalar",
       target = "rate_neo", weight = <floor>, dim = 1L,
       paramIdx = 3L)
  ```
- Keep existing MH scale proposals as well — the adaptive scheduler
  will allocate weight between MH and slice based on ESS/s.
- Width tuning: add `slice_width_rate_loss`, `slice_width_rate_neo`,
  etc. to the tuning list. Initialize to 1.0. Adapt during warmup
  based on average bracket width used.

**Passing paramIdx to C++:** The `do_move_impl` signature already
takes `moveType` and `charIdx`. For slice moves, repurpose `charIdx`
as `paramIdx` (0=treeLength, 1=rateLoss, etc.). This avoids signature
changes.

**Integration with do_move_impl:** The slice sampler needs a different
control flow (no external accept/reject). Two options:

**Option A:** Implement slice as a separate code path in
`run_mcmc_batch_cpp`, before the `do_move_impl` call:
```cpp
if (moveType == 19) {
  // Slice sampling — handles everything internally
  bool ok = slice_scalar_impl(data, states[ch], charIdx, sliceWidth, betas[ch]);
  if (ok) acceptCounts(ch, moveIdx)++;
} else {
  bool accepted = do_move_impl(...);
  ...
}
```

**Option B:** Handle inside `do_move_impl` with an early return after
the slice loop. The slice sampler sets `state->logLik`, `state->logPrior`
directly and returns `true`, bypassing the MH evaluation block.

Option A is cleaner — avoids complicating `do_move_impl` with a fundamentally
different control flow. Go with Option A.

**Passing slice width:** Add a `sliceWidths` parameter to
`run_mcmc_batch_cpp` (NumericMatrix, nChains × nMoves, like
scaleTunings). For non-slice moves, the value is ignored.

---

## Implementation order

1. **Weight floor** (R only, quick) — immediate impact on rate_loss/rate_neo
2. **Joint neo proposal** (C++ case + R wiring) — ~30 min
3. **Tuning budget scaling** (R only, quick) — ensures bandit can detect issues
4. **Slice sampling** (C++ new function + R wiring) — largest change

## Testing

- Existing tests: must all pass (weight floor changes initial weights, so
  any tests asserting exact weight values may need updating).
- New unit tests:
  - Weight floor: verify scalar moves get >= 2% of total weight for large datasets
  - Joint neo: verify Hastings ratio = 2*log(mult), verify rollback on reject
  - Tuning budget: verify scaling with total move weight
  - Slice: verify bracket finding on a simple unimodal target; verify
    temperature handling; verify partial cache used for neo params
- Integration: short Sun2018 run (5k iter) to verify rate_loss/rate_neo
  acceptance rates improved
