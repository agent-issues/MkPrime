# MkPrimeRecover from log file

## Problem

When `RunMkPrime()` is interrupted and the user supplied a named `logFile`
in `MkPrimeMCMC()`, the current `MkPrimeRecover()` cannot build an
`MkPosterior` because:

1. The named log files are not "temp" files — so `.mkp_env$recovery` may
   not reference them (the interrupt handler stores recovery metadata, but
   it checks `active_temp_logs`, and user-supplied log files are not temp).
2. Even when recovery metadata exists, it's session-local and lost on restart.

Meanwhile the log files and checkpoint are sitting on disk with everything
needed to construct a posterior.

## Design

Add a `logFile` parameter to `MkPrimeRecover()`:

```r
MkPrimeRecover(logFile = NULL, checkpointFile = NULL)
```

### Behaviour

1. **`logFile` is NULL** (default): existing behaviour — look for
   `.mkp_env$recovery` metadata from the current session.

2. **`logFile` is provided**: read samples from disk and build an
   `MkPosterior`:
   - Discover per-run files using the `_N.log` naming convention
     (e.g. `"hyoliths.log"` → look for `hyoliths_1.log`, `hyoliths_2.log`,
     …). If those don't exist, try the exact path.
   - Read samples via `ReadMkLog()`.
   - Look for a checkpoint file: explicit `checkpointFile` arg, or
     auto-derive as `sub("\\.[^.]+$", ".ckp", logFile)`.
   - If checkpoint exists, extract `model` and `mcmc` from it.
   - Build `MkPosterior` with available metadata. Fields not recoverable
     from the checkpoint (`data`, `acceptance`, `tuning`, trees) are set to
     `NULL` / empty.
   - Set `result$partial <- TRUE`, `result$stop_reason <- "recovered"`.

### Handling missing `data`

The checkpoint stores `model` and `mcmc` but not `data` (MkPrimeData).
The `print.MkPosterior` method accesses `x$data$nChar` and `x$data$type`
to show character counts. Two options:

- **Option A (chosen):** Guard `print.MkPosterior` so it skips the
  character-count line when `x$data` is NULL. Same for `summary` and any
  other method that accesses `$data`. Minimal, non-breaking change.
- Option B: Start saving `data` in checkpoints. Good idea but separate
  task — the checkpoint would grow, and existing `.ckp` files wouldn't
  have it.

### Multi-run file discovery

When `logFile = "hyoliths.log"`, try numbered files first:
1. Check if `hyoliths_1.log` exists → multi-run pattern.
2. Enumerate `hyoliths_N.log` for N = 1, 2, … until one is missing.
3. If `hyoliths_1.log` doesn't exist, check if `hyoliths.log` itself
   exists → single-run file.

This matches `.LogFilePaths()` and `.OpenLogFiles()` conventions.

## Changes

### `R/streaming.R` — `MkPrimeRecover()`

- Add `logFile = NULL, checkpointFile = NULL` parameters.
- When `logFile` is non-NULL, add the log-file recovery path described
  above, before the existing `.mkp_env$recovery` path.
- Factor out the `MkPosterior`-construction into a small helper so both
  paths share it.

### `R/MkPosterior.R` — `print.MkPosterior()`

- Guard the character-count `cli_ul` item: skip when `x$data` is NULL.
- Guard `x$mcmc$nChains` (fallback to 1L).
- Similarly guard `summary.MkPosterior` if it accesses `$data`.

### `R/burnin.R` — `.PostBurninData()`

- Line 172 accesses `posterior$mcmc$thin` unconditionally. Guard: if
  `mcmc` is NULL, default `thin` to 1L. (This only matters if no
  checkpoint is found — rare.)

### Tests

- `tests/testthat/test-streaming.R` (or new `test-recover.R`):
  - Run a short `RunMkPrime()` with `logFile`, interrupt-safe test via
    `nIter` cap, then call `MkPrimeRecover(logFile = ...)` on the
    resulting logs.
  - Verify the returned object is class `MkPosterior`, has correct
    sample dimensions, and `$partial` is TRUE.
  - Test discovery of `_1.log` / `_2.log` pattern.
  - Test fallback when no checkpoint exists (minimal posterior).

### Documentation

- Update `MkPrimeRecover` roxygen with the new parameters and examples.

## Out of scope

- Saving `data` (MkPrimeData) in checkpoints — useful but separate task.
- Tree recovery from `.nwk` file — could be added later.
- Extending `ResumeMkPrime()` — that's a different entry point.
