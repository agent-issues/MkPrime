# Plan: Pass-through `...` from `RunMkPrime()` to `MkPrimeMCMC()`

## Goal

Allow users to pass MCMC configuration arguments directly to `RunMkPrime()`
instead of constructing a separate `MkPrimeMCMC` object.

**Before:**
```r
mcmc <- MkPrimeMCMC(nIter = 50000, logFile = "out.log")
RunMkPrime(data, model = model, mcmc = mcmc)
```

**After (also works):**
```r
RunMkPrime(data, model = model, nIter = 50000, logFile = "out.log")
```

The explicit `mcmc =` path remains for power users and complex configs.

## Design decisions

1. **Error on conflict.** If the user passes both `mcmc =` and any `...` arg
   that belongs to `MkPrimeMCMC()`, error immediately. No silent merging —
   it's ambiguous which should win.

2. **Only unknown `...` args go to `MkPrimeMCMC()`.** The `...` must not
   contain any argument name that matches a named parameter of `RunMkPrime()`
   itself (`data`, `tree`, `neomorphic`, `knownStates`, `model`, `mcmc`,
   `fixTopology`, `overwrite`). R handles this naturally since those are named
   formals, but we should validate to give a clear error for typos.

3. **Empty `...` preserves existing behaviour.** When neither `mcmc` nor `...`
   are supplied, `MkPrimeMCMC()` with defaults is constructed (status quo).

4. **`match.call()` in `MkPrimeMCMC()` works correctly.** The internal
   `"autoTune" %in% names(match.call())` check (line 299) functions properly
   under `do.call()` — only args actually passed appear in the call.

## Changes

### 1. `R/RunMkPrime.R` — signature and preamble

**Signature** (line 62):
```r
RunMkPrime <- function(data, tree = NULL,
                       neomorphic = integer(0),
                       knownStates = integer(0),
                       model = NULL,
                       mcmc = NULL,
                       fixTopology = FALSE,
                       overwrite = FALSE,
                       ...) {
```

**Preamble** — replace the current `if (is.null(mcmc)) mcmc <- MkPrimeMCMC()`
(line 71) with:

```r
dots <- list(...)
if (length(dots) > 0L && !is.null(mcmc)) {
  cli::cli_abort(c(
    "Supply MCMC options via {.arg mcmc} or {.code ...}, not both.",
    "i" = "Either pass {.code mcmc = MkPrimeMCMC(...)}, or pass \\
          MCMC arguments directly (e.g. {.code nIter = 50000})."
  ))
}
if (is.null(mcmc)) {
  mcmc <- do.call(MkPrimeMCMC, dots)
}
```

### 2. `R/RunMkPrime.R` — roxygen

Add to the existing `@param` block:

```r
#' @param ... Additional arguments forwarded to [MkPrimeMCMC()]. Allows
#'   passing MCMC configuration inline (e.g. `nIter`, `logFile`, `nChains`)
#'   without constructing a separate object. Cannot be combined with an
#'   explicit `mcmc` argument.
```

Update the examples section to show both styles.

### 3. Vignette — optional simplification

The demo run in `vignettes/hyoliths.qmd` could be simplified:

```r
posteriorDemo <- RunMkPrime(
  data       = pd,
  neomorphic = neomorphicChars,
  model      = model,
  nIter      = 5000L,
  thin       = 10L,
  maxWarmup  = 2000L,
  autoTune   = FALSE,
  nRuns      = 2L,
  nChains    = 1L
)
```

The production config is complex enough that a named `MkPrimeMCMC` object
is still clearer, so keep that as-is.

### 4. Tests

Add a handful of targeted tests in a new section of an existing test file
(or a small new file `test-runmkprime-dots.R`):

- `RunMkPrime(d$pd, d$tree, nIter = 500L, maxWarmup = 100L)` works
  (constructs mcmc from dots)
- `RunMkPrime(d$pd, d$tree, mcmc = MkPrimeMCMC(), nIter = 500L)` errors
  (conflict)
- `RunMkPrime(d$pd, d$tree)` works (empty dots, defaults — status quo)

### 5. Out of scope

- Same treatment for `MkPrimeModel` (could do later if desired)
- No changes to `ResumeMkPrime()` (gets mcmc from checkpoint)
- No changes to `MkPrimeMCMC()` itself

## Risk

Low. Purely additive — no existing code paths change. The only new code is
the 6-line preamble and the roxygen update.
