# Plan: Auto-thin default based on number of active moves

**Task:** M-116 — Change `thin` default from hardcoded `10L` to `"auto"`,
where auto resolves to the number of distinct move types at run time.

## Motivation

Storing every single MCMC iteration is wasteful: consecutive samples differ by
at most one parameter. RevBayes defaults to one "full cycle" (one attempt per
move type) between stored samples. MkPrime already has `thin = 10L` which
roughly matches a typical ~10-move configuration, but when optional moves are
enabled/disabled the actual move count can range from 5 to 14+, making 10 a
poor fit.

## Design

### `MkPrimeMCMC()` change

- Change default from `thin = 10L` to `thin = "auto"`.
- Validation: accept `"auto"` (character) or a positive integer. Defer
  resolution of `"auto"` to run time (when the move list is known).

### Resolution in `RunMkPrime()` / `ResumeMkPrime()`

- After `.BuildMoves()` returns the move list, resolve `thin`:
  ```r
  if (identical(mcmc$thin, "auto")) {
    mcmc$thin <- length(moves)
  }
  ```
- This counts *distinct move types*, not total weight. A config with
  {tree_length, branch_lengths, nni, spr, gibbs_spr, gibbs_subtree_swap,
  tbr, kPrime, p, rate_log_sd} = 10 moves → `thin = 10`.
- Mutate `mcmc$thin` in place so downstream code (C++ call, streaming,
  convergence window) sees a concrete integer with no further changes.

### Stepping-stone path

`mkp_stepping_stone()` also calls `.BuildMoves()`. Apply the same resolution
there.

### Resume path

`ResumeMkPrime()` also calls `.BuildMoves()`. Apply the same resolution. For
checkpointed runs, the resolved integer is already stored in the checkpoint's
`mcmc$thin`, so auto-resolution only fires on fresh starts.

### Documentation

- Update `@param thin` roxygen to document `"auto"` and explain the
  move-count heuristic.
- Mention that users can override with any positive integer.

## Files to modify

| File | Change |
|------|--------|
| `R/MkPrimeMCMC.R` | Default `thin = "auto"`, validation accepts character or integer |
| `R/RunMkPrime.R` | Resolve `"auto"` → `length(moves)` after `.BuildMoves()` in `RunMkPrime()`, `ResumeMkPrime()`, and `mkp_stepping_stone()` |
| `tests/testthat/test-validation.R` or new `test-thin.R` | Test that auto resolves correctly, that explicit integers still work |

## What stays the same

- C++ `run_mcmc_batch_cpp` — already receives `thin` as an integer; no change.
- Streaming, convergence window, checkpoint — all use `mcmc$thin` after
  resolution; no change.
- `ReadMkLog()` — unaffected.

## Not doing

- Adaptive thinning during the run (ACT-based). Complexity outweighs benefit;
  the move-count heuristic captures the main win.
