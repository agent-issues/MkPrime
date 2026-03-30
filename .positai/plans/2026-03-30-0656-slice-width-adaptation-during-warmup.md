# Plan: Slice Width Adaptation During Warmup

## Goal

Adapt slice sampler initial widths during warmup to match the posterior's
scale on the log-transformed parameter space. Currently all five slice
widths are fixed at 1.0 throughout the run. An ill-matched width wastes
likelihood evaluations: too narrow → many stepping-out expansions; too
wide → many shrink-in iterations before acceptance.

## Design

### Metric: average stepping-out expansions

The natural diagnostic for slice width quality is the number of
stepping-out steps per call. With a well-tuned width:

- **Target: ~2–4 total expansions** (left + right combined). This means
  the initial bracket is slightly narrower than the slice, and 1–2
  expansions on each side find the boundary efficiently.
- 0 expansions → width too wide (wasting shrink-in evaluations).
- 8+ expansions → width too narrow (many expensive steps to find the boundary).

### Adaptation rule

After each warmup batch, for each slice move with ≥ 10 proposals:

```
avg_expansions = total_expansions / n_proposals
ratio = avg_expansions / target_expansions   # target = 3.0
width_new = width_old * ratio
width_new = clamp(width_new, 0.05, 10.0)    # guard rails on log scale
```

This is a simple proportional controller: if we're doing 6× more
expansions than the target of 3, double the width; if we're doing 1.5×,
halve it. Convergence is fast because the relationship between width and
expansion count is roughly linear.

### Where the data comes from

The C++ `slice_scalar_impl()` already counts stepping-out steps in its
loops (`j` iterates in the two stepping-out for-loops). We need to:

1. **Accumulate** the total step count per slice move across the batch.
2. **Return** it to R alongside `accept_counts` / `propose_counts`.
3. **Use it** in a new `.AdaptSliceWidths()` helper called from warmup.

## Implementation Steps

### 1. C++ changes (`src/mcmc.cpp`)

#### 1a. Add an expansion-count accumulator to `run_mcmc_batch_cpp`

Add a `NumericMatrix sliceExpansions(nChains, nMoves)` alongside
`acceptCounts` / `proposeCounts`. Initialize to zero.

#### 1b. Modify `slice_scalar_impl` to return expansion count

Change the return type from `bool` (accepted) to a struct or just use an
output parameter. Simplest: add an `int& nExpansions` output parameter.
Count `j` iterations in both stepping-out loops and write the total to
the output parameter.

Signature becomes:
```cpp
static bool slice_scalar_impl(McmcData* data, McmcState* state,
                               int paramIdx, double width,
                               double beta, int maxSteps = 10,
                               int* nExpansionsOut = nullptr);
```

In the stepping-out loops, count iterations:
```cpp
int nExp = 0;
for (int j = 0; j < maxSteps; ++j) {
    // ... existing code ...
    ++nExp;
}
for (int j = 0; j < maxSteps; ++j) {
    // ... existing code ...
    ++nExp;
}
if (nExpansionsOut) *nExpansionsOut = nExp;
```

#### 1c. Accumulate in `run_mcmc_batch_cpp`

At the slice call site:
```cpp
int nExp = 0;
accepted = slice_scalar_impl(data, states[ch], charIdx,
                              sliceWidths(ch, moveIdx), betas[ch],
                              10, &nExp);
sliceExpansions(ch, moveIdx) += nExp;
```

#### 1d. Return the matrix

Add to the return `List::create(...)`:
```cpp
_["slice_expansions"] = sliceExpansions
```

### 2. R changes (`R/RunMkPrime.R`)

#### 2a. Accumulate expansion counts in warmup loop

After each batch result (near lines 738–749 where accept/propose counts
are accumulated), also accumulate:
```r
r$chain_slice_exp[[ch]] <- r$chain_slice_exp[[ch]] +
    as.numeric(result$slice_expansions[ch, ])
```

Initialize `chain_slice_exp` in `.InitRun()` as a list of zero-vectors
(same shape as `chain_accept`).

#### 2b. New `.AdaptSliceWidths()` helper

```r
.AdaptSliceWidths <- function(tuning, proposeCount, sliceExpCount,
                               moves, target = 3.0) {
  sliceKeys <- c(
    slice_rate_loss = "slice_width_rate_loss",
    slice_rate_neo = "slice_width_rate_neo",
    slice_rate_log_sd = "slice_width_rate_log_sd",
    slice_tree_length = "slice_width_tree_length",
    slice_beta_scale = "slice_width_beta_scale"
  )
  for (move in moves) {
    nm <- move$name
    tk <- sliceKeys[nm]
    if (is.na(tk) || is.null(tuning[[tk]])) next
    nProp <- proposeCount[nm]
    if (nProp < 10) next
    avgExp <- sliceExpCount[nm] / nProp
    ratio <- avgExp / target
    # Clamp ratio to avoid extreme jumps
    ratio <- max(0.25, min(ratio, 4.0))
    tuning[[tk]] <- tuning[[tk]] * ratio
    tuning[[tk]] <- max(0.05, min(tuning[[tk]], 10.0))
  }
  tuning
}
```

#### 2c. Call from warmup adaptation block

In the warmup block (around line 813), after `.AdaptTuning()`:
```r
r$chain_tuning[[ch]] <- .AdaptSliceWidths(
    r$chain_tuning[[ch]],
    r$chain_propose[[ch]],
    r$chain_slice_exp[[ch]],
    moves
)
```

The expansion counts and proposal counts are cumulative across warmup
(same convention as `chain_accept` / `chain_propose`). The adaptation
divides cumulative expansions by cumulative proposals to get average
expansions per call across the full warmup window. This dampens noise
from individual batches.

#### 2d. Initialise in `.InitRun()`

Add `chain_slice_exp` as a named numeric vector of zeros, same length and
names as `chain_accept`.

### 3. No test changes needed

The existing slow MCMC tests exercise slice moves implicitly. The
adaptation is warmup-only and doesn't affect correctness. The width
guards (0.05–10.0) prevent degenerate behaviour.

## Verification

- Build and run full fast test suite (`devtools::test()`) — should pass
  unchanged since adaptation only fires in warmup.
- Run a short MCMC with `MKPRIME_SLOW_TESTS=true` to confirm slice
  widths actually adapt (add a diagnostic message or inspect the returned
  tuning parameters).

## Notes

- The expansion count is cheap to track (one integer increment per
  stepping-out step, which we already iterate through).
- The target of 3.0 expansions is a soft heuristic. On the log scale
  with width 1.0 (≈ factor of 2.7×), 3 expansions explore ≈ 3 units of
  log-space (factor of ~20×), which is a reasonable initial bracket for
  most rate parameters.
- Guard rails at 0.05 and 10.0 on the log scale correspond to factors of
  ~1.05× and ~22,000× respectively — very permissive bounds.
