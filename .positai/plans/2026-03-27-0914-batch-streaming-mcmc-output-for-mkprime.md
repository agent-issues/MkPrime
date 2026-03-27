# Plan: Batch-Streaming MCMC Output for MkPrime

**Date:** 2026-03-27  
**Status:** Proposed

---

## Motivation

The current design accumulates all MCMC samples in RAM and stores the full
run history in checkpoints. This has two consequences:

1. **No real-time monitoring.** Users cannot open a `.log` file in Tracer
   mid-run to assess mixing and ESS — a core part of the Bayesian phylogenetics
   workflow.
2. **Checkpoints grow without bound.** Late checkpoints are expensive to
   write and read, and redundantly re-serialize the entire run history.

The goal is **batch streaming**: buffer a configurable number of thinned
samples in memory, then flush to a Tracer-compatible `.log` file. Checkpoints
become state-only, storing just enough to resume the MCMC chain. The log
file is the persistent record of the run.

---

## Design overview

```
R outer loop
  ↓ C++ batch (200 iterations)
  ↓ R processes: store in flush buffer (per run)
  ↓ if buf_idx == bufferSize → flush to .log file, reset buffer
  ↓ maintain rolling convergence window (separate, in memory)
  ↓ convergence check uses window → PSRF / ESS
  ↓ checkpoint: state only (no samples)
  ↓ end of run: flush remainder, populate MkPosterior
```

Two modes are supported:

| | **In-memory mode** (current) | **Streaming mode** (new) |
|---|---|---|
| Trigger | `logFile = NULL` (default) | `logFile = "run.log"` |
| Sample storage | Full matrix in RAM | Rolling buffer; flushed to disk |
| Checkpoint content | Full run (state + samples) | State only |
| Tracer compatibility | No (until run ends) | Yes (after each flush) |
| Memory growth | Linear with run length | Bounded by `bufferSize + convWindow` |

In-memory mode is preserved unchanged for short/interactive runs. Streaming
mode activates when `logFile` is set in `MkPrimeMCMC`.

---

## New `MkPrimeMCMC` parameters

Two new parameters are added to `MkPrimeMCMC()`:

```r
logFile    = NULL   # character: path to tab-separated scalar log.
                    # NULL → in-memory mode (existing behaviour).
                    # Non-NULL → streaming mode.
                    # Multi-run: suffixed automatically, e.g.
                    #   "run.log" → "run_1.log", "run_2.log"

bufferSize = 500L   # integer: thinned samples to hold per run before
                    # flushing to logFile. User-tunable.
                    # Smaller → more frequent flushes (lower latency
                    #   for Tracer, slightly more I/O).
                    # Larger  → fewer flushes (better I/O efficiency,
                    #   higher memory use per run).
                    # Ignored in in-memory mode.
```

Validation: `logFile` must be a length-1 character or NULL; `bufferSize`
must be a positive integer.

---

## Log file format (Tracer-compatible)

```
Sample\tlog_posterior\tlog_likelihood\ttree_length\trate_loss\t...
10\t-312.45\t-298.1\t1.83\t0.021\t...
20\t-309.88\t-296.4\t1.91\t0.019\t...
```

- First column is `Sample` (the actual MCMC iteration number, not the
  thinned-sample index). Tracer requires this.
- Tab-separated, no quotes.
- Header is written once on file creation.
- Subsequent flushes append rows without re-writing the header.
- On resume, rows are appended to the existing file.

Tree files (`treeFile`) are unaffected: the existing per-sample Newick
append already achieves low-latency output and is unchanged.

---

## In-memory structures (streaming mode, per run)

```
runs[[r]]$flush_buf     # matrix(NA_real_, bufferSize, nParams)
runs[[r]]$flush_idx     # integer: rows filled in flush_buf

runs[[r]]$conv_window   # matrix(NA_real_, convWindowSize, nParams)
                        # rolling window for PSRF/ESS
runs[[r]]$conv_head     # integer: next write position (circular)
runs[[r]]$conv_filled   # logical: has the window wrapped around?

runs[[r]]$saved_idx     # integer: total thinned samples saved so far
                        # (flush count × bufferSize + flush_idx)
```

`convWindowSize` is computed once at initialization:

```r
convWindowSize <- max(bufferSize, 4L * (checkEvery %/% thin))
```

This ensures the convergence window always spans at least four full
`checkEvery` intervals, giving a reasonable ESS estimate.

The `samples` and `tree_samples` fields present in in-memory mode are
**not allocated** in streaming mode. Trees continue to be appended to
`treeFile` per-sample (existing logic, unchanged).

---

## Files changed

### `R/MkPrimeMCMC.R`
- Add `logFile = NULL` and `bufferSize = 500L` parameters.
- Validate both in the constructor body.
- Document in roxygen.

### `R/RunMkPrime.R`

**Initialization block** (`RunMkPrime`, before the main loop):
- If streaming mode: call `.OpenLogFiles()` to create file(s) and write
  headers; initialize flush buffers and convergence windows.
- If in-memory mode: allocate `samples` and `tree_samples` as before.

**Main batch loop** (inside `for (run in seq_len(nRuns))`):
- After processing each C++ result, route sample storage through
  `.AddToFlushBuffer()`.
- `.AddToFlushBuffer()` also pushes the row into the circular
  convergence window.
- When `flush_idx == bufferSize`: call `.FlushBuffer()` and reset.
- Convergence check (every `checkEvery`): use `conv_window` in streaming
  mode instead of `samples[1:saved_idx, ]`.

**`.SaveCheckpoint()`**:
- In streaming mode: omit `samples` and `tree_samples` from the
  checkpoint list; store `saved_idx` and `logFile` path instead.
- In in-memory mode: unchanged.

**`ResumeMkPrime()`**:
- Detect streaming mode from checkpoint's `logFile` field.
- Open log file(s) in append mode.
- Initialize flush buffers and convergence windows as above.
- Start batching from `checkpoint$iter + 1L`.

**`.BuildResult()`**:
- In streaming mode:
  - Flush any remaining buffer rows.
  - Set `$samples = NULL` in the returned `MkPosterior` object.
  - Set `$logFile` field to the path(s) of the written log(s).
  - Print an informative message directing the user to `ReadMkLog()`.
- In in-memory mode: unchanged.

### `R/streaming.R` (new internal file)

Internal helpers, not exported:

```r
.OpenLogFiles(logFile, paramNames, nRuns)
  # Creates file path(s), writes header, returns path vector.

.FlushBuffer(buffer, nRows, logFile, iterNums)
  # Appends nRows rows of buffer to logFile (write.table, append=TRUE,
  # col.names=FALSE, sep="\t", quote=FALSE).
  # iterNums: integer vector of MCMC iteration numbers for the Sample col.

.AddToFlushBuffer(run, scalarRow, iterNum)
  # Adds one row to flush buffer and convergence window.
  # Returns invisibly; side-effects run list in caller via <<-.
  # (Or: pass run by reference via environment.)

.InitStreamBuffers(nParams, paramNames, bufferSize, convWindowSize)
  # Returns a named list: flush_buf, flush_idx, conv_window,
  # conv_head, conv_filled, saved_idx.
```

### `R/MkPosterior.R` (or `R/streaming.R`)

New exported function:

```r
ReadMkLog(logFile)
```

- Reads a tab-separated `.log` file written by MkPrime.
- Returns a named numeric matrix (same structure as `MkPosterior$samples`).
- Handles multi-run logs (vector of paths → list of matrices, or rbind
  with a `run` column — TBD).
- Validates that the `Sample` column is present and strictly increasing.

---

## Checkpoint format change

Streaming-mode checkpoint omits the heavy fields:

```r
# Streaming checkpoint
list(
  version  = 2L,          # bumped to distinguish from v1
  timestamp,
  iter,
  mcmc,                   # MkPrimeMCMC config (for verification)
  logFile,                # character vector of log file path(s)
  runs = list(
    per_run = list(
      chains     = ...,   # MCMC state (unchanged)
      saved_idx  = ...,   # thinned sample count only (no sample matrix)
      chain_tuning, betas, chain_accept, chain_propose,
      swap_accept, swap_propose
    )
  )
)

# In-memory checkpoint (unchanged)
list(
  version = 1L,
  ...                     # full samples + tree_samples as before
)
```

`ResumeMkPrime` dispatches on `checkpoint$version` (1 = in-memory,
2 = streaming). Version 1 checkpoints remain fully functional.

---

## Open questions (decisions needed before or during implementation)

1. **Multi-run log naming.** If `logFile = "run.log"` and `nRuns = 2`,
   produce `run_1.log` + `run_2.log`? Or a single file with a `Run`
   column? Separate files are simpler for Tracer (each chain opened
   independently); a merged file is simpler for post-processing.
   **Recommendation**: separate files.

2. **`$samples` at end of run.** Should `MkPosterior$samples` be
   auto-populated from `ReadMkLog()` at the end of `RunMkPrime()`?
   Convenient, but defeats the memory-saving purpose for long runs.
   **Recommendation**: leave `$samples = NULL` in streaming mode;
   let the user call `ReadMkLog()` explicitly. Print a clear message.

3. **`convWindowSize` adequacy.** Using `max(bufferSize, 4 × checkEvery/thin)`
   means PSRF is computed on a rolling window, not the full chain. This
   is an approximation. Document the limitation; add a note that
   convergence-based stopping in streaming mode is heuristic.

---

## Implementation order

1. `MkPrimeMCMC.R` — add + validate parameters  
2. `R/streaming.R` — internal helpers  
3. `RunMkPrime.R` — initialization  
4. `RunMkPrime.R` — batch loop routing  
5. `RunMkPrime.R` — `.SaveCheckpoint()` version 2  
6. `RunMkPrime.R` — `ResumeMkPrime()` streaming path  
7. `RunMkPrime.R` — `.BuildResult()` streaming path  
8. `ReadMkLog()` exported function  
9. Tests  

---

## Tests

- `test-streaming.R`:
  - Short run with `logFile` set: verify `.log` file exists, header is
    correct, row count matches `nIter / thin`, Sample column matches
    iteration numbers.
  - Flush boundary: run with `bufferSize = 3`, `nIter = 10`, `thin = 1`
    → verify multiple flushes all land in file.
  - Multi-run: verify separate `_1.log`, `_2.log` files created.
  - Resume: run to iter 500, stop, resume to 1000, verify log is
    contiguous (no duplicate rows, no gap).
  - In-memory mode unchanged: `logFile = NULL` still returns
    `$samples` matrix as before.
  - `ReadMkLog()`: round-trip — write log, read back, compare to
    in-memory result.
  - Checkpoint version dispatch: v1 checkpoint resumes in in-memory
    mode; v2 checkpoint resumes in streaming mode.
