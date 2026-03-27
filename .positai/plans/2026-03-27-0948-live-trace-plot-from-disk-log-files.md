# Plan: Live trace plot from disk log files — M-076

**Date:** 2026-03-27  
**Author:** Agent C  
**Status:** Draft

---

## Background

`RunMkPrime()` now writes a Tracer-compatible TSV log file to disk when
`logFile` is set in `MkPrimeMCMC()` (M-075). This enables a new, cleaner
approach to live monitoring:

- The existing `progressFn = "default"` / `MkpTracePlot()` mechanism
  fires inside the MCMC loop, passing in-memory data. It works for a
  single R process but can't monitor truly parallel runs (separate R
  processes writing to separate log files).
- A **standalone disk-based watcher** can be called independently —
  even in a different R session — as long as it can read the log files.
  It polls, reads the latest data, and redraws the plot window.

The user's request is for this simpler "live plot in the plot window"
pattern, enabled by the log files now being on disk.

---

## Design

### New exported function: `MkpWatchLog()`

```r
MkpWatchLog(logFiles,
            interval  = 2,
            params    = NULL,
            maxSamples = 2000,
            warmup    = NULL)
```

**Parameters:**

| Arg | Type | Default | Meaning |
|-----|------|---------|---------|
| `logFiles` | `character` | — | One or more log file paths. One file per run; use the `$logFile` field of `MkPrimeMCMC()` directly. For a multi-run analysis, e.g. `"run.log"` → `c("run_1.log", "run_2.log")`. |
| `interval` | `numeric` | `2` | Seconds between poll/redraw cycles. |
| `params` | `character` or `NULL` | `NULL` | Column names to plot. Default (`NULL`) auto-selects key scalar parameters (see below). |
| `maxSamples` | `integer` | `2000` | Maximum number of rows read per file (tail-subsample for speed when the log is large). |
| `warmup` | `integer` or `NULL` | `NULL` | If provided, draws a vertical dashed line at the warmup/sampling boundary on each trace panel. Not required — log files only contain post-warmup samples anyway. |

**Behaviour:**

1. Validates that at least one `logFiles` path exists. If none exist yet,
   waits up to 30 s in 1-second increments with a `cli_alert_info` message,
   then aborts if still absent.
2. Enters a `repeat` loop:
   a. Read each log file with `ReadMkLog()`, capping at `maxSamples` rows
      (tail of file = most recent samples).
   b. Call `.WatchLogPlot()` to draw trace panels.
   c. `Sys.sleep(interval)`.
3. The loop is wrapped in `tryCatch(..., interrupt = ...)` so Ctrl+C gives
   a clean message and returns the last-read data invisibly.
4. Returns the last sample matrix (or named list of matrices for multi-run)
   **invisibly**, so the user can capture it:
   ```r
   last <- MkpWatchLog("run.log")   # Ctrl+C to stop
   summary(last)
   ```

**Default `params` selection** (mirrors `MkpTracePlot`):
```
log_posterior, tree_length, rate_loss, rate_log_sd, p
```
Plus any `kPrime_*` columns present, up to a maximum of 3 (to avoid
overwhelming the panel layout). Branch-length (`br_*`) and
`log_likelihood` columns are excluded by default.

---

### New internal function: `.WatchLogPlot(matList, params, nSamples)`

Separated from `MkpWatchLog` for testability. Takes a named list of
matrices (one per run, the output of `ReadMkLog()`) and draws the trace
plot.

Shares layout logic with the existing `.PlotTracePanel()`:
- Multi-panel grid: up to 3 columns, rows = `ceiling(nParams / 3)`
- Multiple runs → `hcl.colors(nRuns, "Set 2")`, one colour per matrix
- X-axis: MCMC iteration number from `rownames(mat)` (the `Sample` column
  stored by `ReadMkLog()`), **not** sample index — shows true iteration
  numbers on the x-axis
- Title line (via `mtext`): total samples read, last iteration number,
  wall-clock time of last read (via `format(Sys.time(), "%H:%M:%S")`)

---

### File changes

| File | Change |
|------|--------|
| `R/PlotDuringMCMC.R` | Add `MkpWatchLog()` (exported) and `.WatchLogPlot()` (internal) |
| `tests/testthat/test-progress.R` | Add tests for `MkpWatchLog` and `.WatchLogPlot` |
| `to-do.md` | M-076 → DONE (C) |
| `completed-tasks.md` | Append M-076 row |

No changes to `RunMkPrime.R`, `streaming.R`, or the existing
`MkpTracePlot` / `MkpPngProgress` functions.

---

## Testing strategy

The polling loop itself (the `repeat { Sys.sleep }` part) is not
testable in unit tests. Test the pieces separately:

1. **`.WatchLogPlot()` renders without error** — construct a named list of
   mock matrices (matching `ReadMkLog()` output format), call
   `.WatchLogPlot()`, expect no error.

2. **`.WatchLogPlot()` handles single-run and multi-run** — one-element
   list vs two-element list.

3. **`MkpWatchLog()` errors on missing file** — expect `cli_abort` when
   `logFiles` doesn't exist and `maxWait = 0` (expose an internal
   parameter for tests).

4. **`MkpWatchLog()` with `maxUpdates = 1`** — add an internal
   `maxUpdates = Inf` parameter; call with `maxUpdates = 1` and a
   pre-existing log file. Should read once, plot, and return cleanly.

5. **Return value structure** — after one update, check the returned
   matrix has the right column names and row count.

---

## Usage example

```r
# Start a long run with logFile
result <- RunMkPrime(pd, tree,
  mcmc = MkPrimeMCMC(nRuns = 2L, logFile = "hyoliths.log"))

# OR: start in the background, then watch in the console
# (background job or callr::r_bg())
MkpWatchLog(c("hyoliths_1.log", "hyoliths_2.log"), interval = 3)
```

---

## What this does NOT include (out of scope)

- **Shiny dashboard**: `MkpPngProgress` + a Shiny app already covers
  that use case. `MkpWatchLog` is for the interactive R console.
- **Auto-stopping the watcher**: When the MCMC finishes, the log file
  stops growing. The watcher just keeps polling (each cycle shows the
  same final state). The user presses Ctrl+C when done. A future
  enhancement could check for a sentinel file, but not needed now.
- **Changes to `progressFn`**: The callback infrastructure is unchanged.
  Both approaches coexist.

---

## Task ID

**M-076**: `MkpWatchLog()` — live trace plot from disk log files  
**Priority**: P2  
**Estimate**: ~100 lines of code, ~5 tests
