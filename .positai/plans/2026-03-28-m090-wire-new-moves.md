# M-090: Wire Gibbs/Weighted moves into dispatch + MkPrimeMCMC API

**Task:** Make the five new move types (GibbsSPR, GibbsSubtreeSwap,
WeightedBranchScale, WeightedSPR, WeightedSubtreeSwap) user-configurable
and include them in the move pool when enabled.

**Date:** 2026-03-28

---

## Current state

- **C++ dispatch is complete.** Cases 10–14 are already in `do_move_impl()`
  and functional (tested via direct `do_move_cpp()` calls in unit tests).
- **R-side wiring is missing.** `.BuildMoves()` does not create move entries
  for any of the five new types. `.kMoveTypes` is missing entries for
  `gibbs_spr` (10) and `gibbs_subtree_swap` (11). `MkPrimeMCMC()` has no
  parameters to enable/disable these moves.
- **`nBranchBins` is hardcoded** as literal `10` in `do_move_impl()` cases
  12, 13, 14. No way for the user to control it.

---

## Design

### 1. New MkPrimeMCMC parameters

Add six new parameters to `MkPrimeMCMC()`:

| Parameter | Type | Default | Rationale |
|-----------|------|---------|-----------|
| `gibbsSpr` | logical | `TRUE` | Cheap O(N) move; on by default |
| `gibbsSubtreeSwap` | logical | `TRUE` | Cheap O(N) move; on by default |
| `weightedBranchScale` | logical | `FALSE` | O(B) per call; off by default |
| `weightedSpr` | logical | `FALSE` | O(N×B) per call; expensive, off by default |
| `weightedSubtreeSwap` | logical | `FALSE` | O(N×B) per call; expensive, off by default |
| `nBranchBins` | integer | `10L` | Number of bins for weighted moves |

These are stored on the MkPrimeMCMC S3 object and passed to `.BuildMoves()`.

### 2. Move pool construction (`.BuildMoves()`)

New moves are inserted in the topology moves block (requires `!fixTopology
&& nEdge >= 5L`). Order within the block:

```
nni             (existing)
spr             (existing)
gibbs_spr       (if mcmc$gibbsSpr)
gibbs_subtree_swap  (if mcmc$gibbsSubtreeSwap)
weighted_branch_scale (if mcmc$weightedBranchScale)   -- branch move, not topology
weighted_spr    (if mcmc$weightedSpr)
weighted_subtree_swap (if mcmc$weightedSubtreeSwap)
```

Note: `weighted_branch_scale` modifies branch fractions, not topology, so
it could live outside the topology guard. However, it only makes sense when
there are enough edges to be interesting, and grouping it here keeps the
code simpler. The guard is `nEdge >= 5L`, which is always true for real
analyses.

**Default weights:**

| Move | Default weight | Notes |
|------|---------------|-------|
| `gibbs_spr` | `max(1, nEdge / 4)` | Same as plain SPR (both O(N)) |
| `gibbs_subtree_swap` | `max(1, nEdge / 6)` | Slightly less than SPR |
| `weighted_branch_scale` | `max(1, nEdge / 6)` | O(B) evals |
| `weighted_spr` | `max(1, nEdge / 8)` | O(N×B), expensive |
| `weighted_subtree_swap` | `max(1, nEdge / 8)` | O(N×B), expensive |

These weights are intentionally conservative starting points. M-092
(adaptive scheduler) will tune them based on ESS/wall-time during warmup.

### 3. `.kMoveTypes` update

Add the two missing entries:

```r
gibbs_spr = 10L,
gibbs_subtree_swap = 11L
```

(Codes 12, 13, 14 already present as `weighted_branch_lengths`,
`weighted_spr`, `weighted_subtree_swap`.)

### 4. `nBranchBins` → C++ via `McmcData`

Add `int nBranchBins = 10;` field to the `McmcData` struct (in
`mcmc_state.h`). Default is 10 so all existing tests and code paths
work unchanged.

Add a C++ setter:

```cpp
// [[Rcpp::export]]
void set_branch_bins(SEXP dataPtr, int nBins) {
  Rcpp::XPtr<McmcData>(dataPtr)->nBranchBins = nBins;
}
```

In `RunMkPrime()` and `ResumeMkPrime()`, call `set_branch_bins(mcmcData,
mcmc$nBranchBins)` right after `.InitMcmcData()` — before any batch calls.

In `do_move_impl()`, replace the hardcoded `10` with `data->nBranchBins`
in cases 12, 13, 14.

This approach:
- Avoids changing the 19-parameter `prepare_mcmc_data()` signature
- Keeps nBranchBins accessible everywhere via `data->` (no threading)
- Default of 10 means existing direct-call tests work unchanged

### 5. Adaptation (`.AdaptTuning()`)

The five new moves have **no tunable parameters** — they are Gibbs-like
(no step size to adapt). Add them to the `targets` and `tuningKeys`
lookups:

```r
targets <- c(
  ...,
  gibbs_spr = NA_real_,       # no target needed
  gibbs_subtree_swap = NA_real_,
  weighted_branch_scale = NA_real_,
  weighted_spr = NA_real_,
  weighted_subtree_swap = NA_real_
)

tuningKeys <- c(
  ...,
  gibbs_spr = NA_character_,
  gibbs_subtree_swap = NA_character_,
  weighted_branch_scale = NA_character_,
  weighted_spr = NA_character_,
  weighted_subtree_swap = NA_character_
)
```

The adaptation loop already skips moves where `tuningKeys[nm]` is `NA`.

### 6. `.BuildScaleTuningMatrix()`

The new moves don't use scale tuning. The existing `default` case
already returns 0.5. No changes needed — the switch statement's default
covers unknown names.

### 7. Documentation

- **`MkPrimeMCMC()`**: Add `@param` entries for the six new parameters.
  Add a `## Gibbs and weighted moves` section to `@details` explaining
  costs and when to enable weighted moves.
- **`RunMkPrime()`**: No signature change; no doc change needed.

---

## Files changed

| File | What changes |
|------|-------------|
| `R/MkPrimeMCMC.R` | 6 new params, validation, `@param` docs, `@details` section |
| `R/RunMkPrime.R` | `.BuildMoves()` gains new moves; `.kMoveTypes` gains codes 10, 11; `.AdaptTuning()` gains entries; `set_branch_bins()` call in `RunMkPrime()` and `ResumeMkPrime()` |
| `src/mcmc_state.h` | `int nBranchBins = 10;` field on `McmcData` |
| `src/mcmc.cpp` | Replace hardcoded `10` with `data->nBranchBins` in cases 12–14 |
| `src/mcmc_likelihood.cpp` | `set_branch_bins()` function (or a new `mcmc_config.cpp`) |
| `R/RcppExports.R` | Auto-generated by `Rcpp::compileAttributes()` |
| `src/RcppExports.cpp` | Auto-generated |
| `tests/testthat/test-m090-move-wiring.R` | New tests (see below) |

---

## Correctness criteria (tests)

### test-m090-move-wiring.R

1. **MkPrimeMCMC param storage.** Create config with various enable/disable
   combos; verify the returned object stores all six fields correctly.

2. **nBranchBins validation.** `nBranchBins = 0` and `nBranchBins = -1`
   error; `nBranchBins = 5L` stores correctly.

3. **Move pool construction — defaults.** With default `MkPrimeMCMC()`,
   `.BuildMoves()` includes `gibbs_spr` and `gibbs_subtree_swap` but NOT
   the three weighted moves.

4. **Move pool construction — all enabled.** With all five enabled,
   `.BuildMoves()` includes all five new moves with correct names.

5. **Move pool construction — all disabled.** With all five disabled
   (including `gibbsSpr = FALSE`), `.BuildMoves()` returns only the
   original moves (tree_length, branch_lengths, nni, spr, kPrime, etc.).

6. **Move pool construction — fixTopology.** With `fixTopology = TRUE`,
   none of the four topology moves (gibbs_spr, gibbs_subtree_swap,
   weighted_spr, weighted_subtree_swap) appear.
   `weighted_branch_scale` is included inside the topology guard, so
   it also doesn't appear when fixTopology is TRUE. This is fine — the
   move barely matters on a 4-edge tree.

7. **Move type codes.** All five new move names resolve to the correct
   integer codes via `.kMoveTypes`.

8. **nBranchBins reaches C++.** Create McmcData, call
   `set_branch_bins(ptr, 5L)`, run a weighted move via `do_move_cpp()`
   with moveType 12 — verify it doesn't error (functional proof that the
   config propagates).

9. **Full MCMC run with Gibbs moves.** `RunMkPrime()` with
   `gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE` on a small dataset;
   verify it produces a valid MkPosterior. (Integration test; guard
   with `skip_slow_tests()`.)

10. **Full MCMC run with all moves.** Same but with all five enabled.
    (Integration test; guard with `skip_slow_tests()`.)

---

## Implementation order

1. `src/mcmc_state.h` — add `nBranchBins` field
2. `src/mcmc_likelihood.cpp` — add `set_branch_bins()` export
3. `src/mcmc.cpp` — replace hardcoded 10 with `data->nBranchBins`
4. `Rcpp::compileAttributes()` to regenerate exports
5. `R/MkPrimeMCMC.R` — add parameters, validation, docs
6. `R/RunMkPrime.R` — `.kMoveTypes`, `.BuildMoves()`, `.AdaptTuning()`,
   `set_branch_bins()` calls in RunMkPrime/ResumeMkPrime
7. `tests/testthat/test-m090-move-wiring.R` — all tests
8. Build + full test suite pass
