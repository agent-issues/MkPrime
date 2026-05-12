# M-092: Adaptive Move Scheduler

**Date:** 2026-03-28
**Branch:** `feature/gibbs-weighted-moves` (mkp-gibbs worktree)
**Depends on:** M-090 (complete)

## Goal

During warmup, learn optimal move frequencies by tracking per-move
acceptance rates and wall-clock cost, then reweight moves via softmax
with temperature annealing. Freeze weights at warmup end to preserve
detailed balance. Allow user-pinned weights via `moveWeights` parameter
in `MkPrimeMCMC()`.

---

## 1. Mathematical Specification

### Score formula

For each move `m` in the pool with at least `minProposals` (20)
proposals from the cold chain of run 1:

```
accept_rate[m] = n_accepted[m] / n_proposed[m]
mean_cost_s[m] = total_time_ns[m] / (n_proposed[m] × 1e9)
score[m]       = accept_rate[m] / max(mean_cost_s[m], 1e-9)
```

Score is "acceptances per second" — a proxy for ESS/wall-time
efficiency.

### Softmax with temperature annealing

```
progress = batchEnd / warmup          # 0 → 1 over warmup
T        = T_start + (T_end - T_start) × progress
         = 2.0 − 1.5 × progress       # 2.0 → 0.5
```

For free (non-pinned) moves with sufficient proposals:

```
raw[m]  = exp(score[m] / T)
prob[m] = raw[m] / sum(raw) × budget
```

where `budget = 1 − sum(pinned_weights)`.

### Floor

After softmax, apply floor:

```
prob[m] = max(w_min, prob[m])   for each free move
```

where `w_min = 0.05 / n_free_moves` (scales with pool size so total
floor allocation stays reasonable). Then renormalize free moves to fill
`budget`.

### Pinned weights

User supplies `moveWeights = list(nni = 0.3, spr = 0.2, ...)`.
Pinned moves have fixed probabilities; remaining `budget` goes to
free moves. Validate: `sum(pinned) <= 1.0`, names are valid move
names. Warn if a pinned move is not in the pool (e.g., NNI with
`fixTopology = TRUE`).

### Freezing

At the first batch where `batchEnd > warmup`, the current moveWeights
are frozen. All subsequent batches use the frozen weights. This
preserves detailed balance in the sampling phase.

### Moves with insufficient proposals

Moves that haven't accumulated `minProposals` (20) proposals by a given
adaptation window keep their current weight (unchanged). This avoids
noisy early adaptation for rare moves.

---

## 2. C++ Changes

### File: `src/mcmc.cpp`

**Change:** Add per-move wall-time tracking to `run_mcmc_batch_cpp()`.

1. Add `#include <chrono>` at the top of the file.

2. After accept/propose counter initialization (line ~1503), add:

   ```cpp
   NumericMatrix moveTimeNs(nChains, nMoves);
   ```

3. Wrap the `do_move_impl()` call (lines 1540–1547) with timing:

   ```cpp
   auto t0 = std::chrono::steady_clock::now();
   bool accepted = do_move_impl(
     data, states[ch], moveType, charIdx,
     chainScaleTunings(ch, moveIdx),
     chainBsmpTunings[ch], chainIntWalkWins[ch],
     betas[ch]
   );
   auto t1 = std::chrono::steady_clock::now();
   moveTimeNs(ch, moveIdx) +=
     (double)std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count();
   if (accepted) acceptCounts(ch, moveIdx)++;
   ```

4. Add `move_time_ns` to the return List:

   ```cpp
   return List::create(
     _["accept_counts"]  = acceptCounts,
     _["propose_counts"] = proposeCounts,
     _["move_time_ns"]   = moveTimeNs,       // NEW
     _["swap_accept"]    = swapAccept,
     // ... rest unchanged
   );
   ```

**Impact:** Minimal. `std::chrono::steady_clock::now()` is ~20ns on
modern CPUs; with ~10 moves per iteration, overhead is ~200ns/iter
(negligible vs. likelihood evals).

No signature change — the return value gains one new element.
Existing R code that accesses `result$accept_counts` etc. is unaffected.

---

## 3. R Changes

### File: `R/MkPrimeMCMC.R`

Add `moveWeights = NULL` parameter:

```r
MkPrimeMCMC <- function(..., moveWeights = NULL, ...) {
```

Validation (in the body):
- If not NULL, must be a named numeric vector or list coercible to one.
- All names must be in the set of valid move names (`.kMoveTypes` names).
- All values must be positive (`> 0`).
- Sum must be `<= 1.0`. If exactly 1.0, adaptive scheduling is fully
  disabled.
- Warn if any named move is not in the current pool (defer to
  `.BuildMoves()` to check; here just validate types).

Store in the returned mcmc list as `moveWeights`.

### File: `R/RunMkPrime.R`

#### New function: `.AdaptMoveWeights()`

```r
.AdaptMoveWeights <- function(currentWeights, acceptCount, proposeCount,
                               moveTimeNs, moveNames, pinnedWeights,
                               warmupProgress, tStart = 2.0, tEnd = 0.5,
                               wMin = 0.05, minProposals = 20L)
```

- `currentWeights`: numeric vector (current move probabilities)
- `acceptCount`: named integer vector (cold chain, cumulative)
- `proposeCount`: named integer vector (cold chain, cumulative)
- `moveTimeNs`: named numeric vector (cold chain, cumulative, in ns)
- `moveNames`: character vector (move names, same order as weights)
- `pinnedWeights`: named numeric vector or NULL
- `warmupProgress`: fraction of warmup completed (0 to 1)

Algorithm:
1. Identify free vs. pinned moves.
2. For free moves with `proposeCount >= minProposals`, compute score.
3. Compute softmax temperature from warmupProgress.
4. Apply softmax → raw weights for scoreable free moves.
5. Free moves below minProposals keep current weight, scaled into
   the free budget.
6. Apply per-move floor.
7. Normalize free moves to fill `budget = 1 - sum(pinned)`.
8. Return updated weight vector.

#### New function: `.LogMoveWeights()`

```r
.LogMoveWeights <- function(moveWeights, moveNames, logFilePaths)
```

Writes a `#` comment line to each log file:
```
# Adapted move weights: nni=0.12 spr=0.08 gibbs_spr=0.15 ...
```

Tracer ignores `#` lines, so this is safe. Also emits a
`cli_alert_info` to the console.

#### Modifications to `.InitRun()`

Add `move_time_ns` tracking (analogous to `chain_accept`/`chain_propose`):

```r
chainTimeNs <- vector("list", nChains)
for (ch in seq_len(nChains)) {
  chainTimeNs[[ch]] <- numeric(length(moves))
  names(chainTimeNs[[ch]]) <- moveNames
}
```

Add to the returned list as `chain_time_ns`.

#### Modifications to `RunMkPrime()`

1. **Normalize initial weights** (after line 152):

   ```r
   moveWeights <- vapply(moves, `[[`, numeric(1), "weight")
   moveWeights <- moveWeights / sum(moveWeights)  # normalize to probs
   ```

2. **Apply user-pinned weights** (overlay pinned values):

   ```r
   pinnedWeights <- .ResolvePinnedWeights(mcmc$moveWeights, moveNames)
   if (!is.null(pinnedWeights)) {
     for (nm in names(pinnedWeights)) {
       idx <- match(nm, moveNames)
       if (!is.na(idx)) moveWeights[idx] <- pinnedWeights[nm]
     }
     # Renormalize free moves to fill remaining budget
     moveWeights <- .NormalizeMoveWeights(moveWeights, pinnedWeights)
   }
   ```

3. **Accumulate timing data** (after line 222, alongside accept/propose):

   ```r
   for (ch in seq_len(mcmc$nChains)) {
     r$chain_time_ns[[ch]] <- r$chain_time_ns[[ch]] +
       as.numeric(result$move_time_ns[ch, ])
   }
   ```

4. **Move weight adaptation** (inside `if (batchEnd <= mcmc$warmup)`
   block, after `.AdaptTuning()` calls, before `runs[[run]] <- r`):

   ```r
   # Adapt move weights using run 1's cold chain stats
   if (run == 1L) {
     moveWeights <- .AdaptMoveWeights(
       moveWeights, r$chain_accept[[1L]], r$chain_propose[[1L]],
       r$chain_time_ns[[1L]], moveNames, pinnedWeights,
       warmupProgress = batchEnd / mcmc$warmup
     )
   }
   ```

5. **Log final weights** (after warmup ends — detect transition):

   Track `weightsLogged <- FALSE` before the loop. Inside the loop,
   after the adaptation block:

   ```r
   if (batchEnd > mcmc$warmup && !weightsLogged) {
     if (isStreaming) .LogMoveWeights(moveWeights, moveNames, logFilePaths)
     cli::cli_alert_info("Move weights frozen: {.LogMoveWeightsStr(moveWeights, moveNames)}")
     weightsLogged <- TRUE
   }
   ```

#### Modifications to `ResumeMkPrime()`

Mirror all changes from `RunMkPrime()`:

1. Load `moveWeights` from checkpoint if available (else compute fresh).
2. Restore `chain_time_ns` from checkpoint.
3. Accumulate timing each batch.
4. Adapt during warmup (if still in warmup phase).
5. Log final weights on warmup→sampling transition.

#### Modifications to `.SaveCheckpoint()`

Add `moveWeights` to the checkpoint payload:

```r
payload <- list(runs = serialRuns, mcmc = mcmc, iter = iter,
                timestamp = Sys.time(), version = version,
                moveWeights = moveWeights)  # NEW
```

Timing accumulators (`chain_time_ns`) are already stored in `runs[[r]]`
since they're part of the run list.

---

## 4. File Change Summary

| File | Change |
|------|--------|
| `src/mcmc.cpp` | `#include <chrono>`, timing around `do_move_impl()`, return `move_time_ns` |
| `R/MkPrimeMCMC.R` | `moveWeights` parameter + validation |
| `R/RunMkPrime.R` | `.AdaptMoveWeights()`, `.LogMoveWeights()`, timing accumulation in `RunMkPrime()` + `ResumeMkPrime()`, checkpoint saves moveWeights |
| `R/RcppExports.R` | Auto-regenerated (no manual edit) |
| `src/RcppExports.cpp` | Auto-regenerated (no manual edit) |
| `tests/testthat/test-m092-adaptive-scheduler.R` | New test file |

---

## 5. Test Plan

### Unit tests for `.AdaptMoveWeights()`

1. **Score computation:** Given known accept/propose/timing, verify
   scores are `accept_rate / mean_cost_s`.
2. **Softmax:** At T=very_large, weights should be nearly uniform.
   At T=very_small, weights should concentrate on highest-score move.
3. **Temperature annealing:** `warmupProgress=0` → T=2.0;
   `warmupProgress=1` → T=0.5.
4. **Floor:** No free move falls below `wMin / nFreeMoves`.
5. **Pinned weights:** Pinned moves keep exact specified values;
   free moves sum to `1 - sum(pinned)`.
6. **All pinned:** If all moves are pinned (sum ≈ 1), return weights
   unchanged.
7. **Insufficient proposals:** Moves below `minProposals` keep
   current weight (proportionally rescaled within free budget).
8. **Normalization:** Output always sums to 1.0 (within tolerance).

### Integration tests

9. **C++ timing returned:** Run a short batch, verify
   `result$move_time_ns` is a numeric matrix of correct dimension
   with all values ≥ 0.
10. **Weights adapt during warmup:** Run ~400 iterations with
    `warmup = 200`. Verify moveWeights at end differ from initial.
11. **Weights frozen after warmup:** Verify moveWeights don't change
    in post-warmup batches.
12. **User-pinned weights respected:** Set `moveWeights = list(nni = 0.4)`
    in `MkPrimeMCMC()`. Verify NNI weight stays at 0.4 throughout.
13. **Checkpoint preserves weights:** Save checkpoint during warmup,
    resume, verify adapted weights are restored.
14. **Log file comment:** With streaming enabled, verify `#` comment
    line appears in log file after warmup.
15. **MkPrimeMCMC validation:** Invalid moveWeights (negative, sum > 1,
    bad names) produce appropriate errors.

### Slow integration test (guarded with `skip_slow_tests()`)

16. **Full MCMC with adaptive scheduler:** Run a small analysis to
    completion, verify valid MkPosterior returned and no errors.

---

## 6. Correctness Criteria

- **Detailed balance preserved:** Weights are frozen after warmup.
  During warmup, adaptation doesn't affect correctness since warmup
  samples are discarded.
- **No regression:** All existing tests pass unchanged. The only C++
  change is adding timing around existing code and a new return field;
  existing return fields unchanged.
- **Score formula is sensible:** Moves that accept more often and cost
  less per proposal get higher scores. Gibbs moves (100% acceptance)
  are weighted by their cost. Expensive weighted moves are naturally
  down-weighted if their cost exceeds the acceptance benefit.
- **Numerics:** Softmax overflow guarded by subtracting max(score) 
  before exponentiation. Floor prevents any move from being starved.
