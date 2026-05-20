# HARNESS-001 investigation: `maxTime=7.5h` graceful exit didn't fire

## Summary

`run_one.R` passes `maxTime = 7.5 * 3600` and `nRuns = 2L, maxRhat = 1.1` to
`MkPrimeMCMC()`. That config routes the dispatcher through `.RunSerialRuns`
(R/RunMkPrime.R:385-394). In serial mode, the **inner**
`.RunMkPrimeSingleRun` resets its own `startTime` clock at the top of each
run — so each of the two runs is allowed its own 7.5 h budget, giving a
worst-case total wall budget of 15 h on top of warmup overhead. SLURM kills
at 8 h, well before the second run's `maxTime` check can fire.

## (a) Where `maxTime` checks fire

Three sites in `R/RunMkPrime.R`:

1. **Inner batch loop** in `.RunMkPrimeSingleRun`, line 1300-1317. Check
   fires *between batches*, after `run_mcmc_batch_cpp()` returns. On the
   `max_time` break it flushes the stream buffer and writes a checkpoint.
   `startTime` for this clock is set fresh at **line 757**, on entry to the
   function.
2. **`.RunSerialRuns` Phase 2 loop** (R-hat resumption), line 1556-1567.
   Uses a `startTime` set at line 1482 (entry to `.RunSerialRuns`). This is
   the only `maxTime` check that spans multiple serial runs, but Phase 2
   only runs after Phase 1 finishes all `nRuns` runs.
3. **`.RunParallelRuns` polling loop**, line 1705. Not on the hamilton code
   path (`parallel = FALSE`, default in `MkPrimeMCMC`).

The inner check uses `samplingBatch = 5000L` (R/RunMkPrime.R:693) as the
batch unit, so the granularity of the wall-clock check is one Sample-phase
batch.

## (b) Why it didn't fire in production

**Primary cause — per-run `startTime` reset.** `run_one.R` requests
`nRuns = 2L` + `maxRhat = 1.1`, which selects the `.RunSerialRuns` branch
(R/RunMkPrime.R:385). Inside Phase 1 (lines 1487-1512) each run is invoked
through `.RunMkPrimeSingleRun`, which resets `startTime` at line 757. The
Phase-1 loop only breaks early on `stop_reason == "cancelled"` (line 1506),
**not on `"max_time"`** — so a run that hits its 7.5 h ceiling returns
cleanly, and run 2 starts a fresh 7.5 h budget. SLURM SIGKILL's at 8 h
before run 2 can hit its own `max_time` check.

**Contributing cause — coarse check granularity.** Even within a single
run, the check only fires between batches. With `samplingBatch = 5000L`
and 4 PT chains on real datasets (large trees, 17+ taxa, hundreds of
characters), a single batch can take many minutes. If the in-flight batch
straddles the 7.5 h boundary, the gap until next check eats further into
the 30-min SLURM cushion. The C++ batch call (`run_mcmc_batch_cpp`) does
not poll for elapsed time and does not yield to R.

**Tertiary — Warmup phase is not bounded by `maxTime` either**, but it
*is* bounded by `mcmc$warmup` (capped at `maxWarmup = 5000L` in `run_one`),
so warmup can't run away on its own. It still consumes the inner
`startTime` budget, but the check fires at every batch boundary in Warmup
too (the `if (!is.null(mcmc$maxTime) …)` at line 1300 is not gated by
`phase`). So warmup *will* exit on `maxTime` — only the per-run reset is
the real harness bug.

## (c) Recommended fix

Two layered fixes in `R/RunMkPrime.R`:

1. **Hard fix — share the wall-clock budget across serial runs.** In
   `.RunSerialRuns` (around line 1482), record a single shared
   `startTime`. Pass a derived per-run budget down to each
   `.RunMkPrimeSingleRun` call: e.g. inject
   `innerMcmc$maxTime <- mcmc$maxTime - (proc.time()["elapsed"] - startTime)`
   before each `.RunMkPrimeSingleRun()` call (lines 1493-1505 and the
   analogous Phase-2 epoch calls around line 1574). Also break the Phase-1
   loop on `stop_reason == "max_time"` (line 1506) — currently it only
   breaks on `"cancelled"`.

2. **Defence in depth — keep the SLURM cushion realistic.** Reduce
   `samplingBatch` from `5000L` to e.g. `1000L` for hamilton-class runs,
   or expose it via `MkPrimeMCMC()` so `run_one.R` can override. With
   5000-iter batches the inner check has multi-minute jitter; 1000-iter
   batches let `maxTime` land within ~60 s of the requested budget.

A separate small fix in `data-raw/hamilton/run_one.R`: wrap the
`RunMkPrime()` call so that even on an interrupt or shorter cushion, a
`partial` summary `.rds` is written from the streaming log file (via
`MkPrimeRecover()` / `ReadMkLog()`) before the script exits — so future
SIGKILL near-misses still leave a non-trivial summary.

## (d) Scope: which phases are affected

- **Sample phase**: affected — convergence/checkpoint AND `maxTime` checks
  are batch-granular and the per-run reset wastes 50% of the budget on
  retry.
- **Warmup**: per-run reset still applies (each run starts its own
  warmup), but `maxWarmup = 5000L` caps it independently, so warmup itself
  can't blow the wall.
- **Tuning**: same as Warmup; bounded by `tuningBudget` / `tuningRounds`.
- **`.RunParallelRuns`** (line 1693): not used by this harness, but it
  *does* enforce `maxTime` correctly because it polls in R after a single
  shared `startTime`. Worth referencing in the fix as the model.

## Recommended regression test

Add to `tests/testthat/test-stopping.R`:

```r
test_that("serial nRuns=2 + maxRhat respects total maxTime budget", {
  # Total wall time must be within maxTime + one-batch slack,
  # not 2 * maxTime.
  start <- proc.time()["elapsed"]
  res <- RunMkPrime(pd, tree,
    mcmc = MkPrimeMCMC(nRuns = 2L, thin = 5L,
                       maxWarmup = 200L, minWarmup = 200L,
                       autoTune = FALSE,
                       maxRhat = 1.001,   # unreachable -> forces maxTime
                       checkEvery = 100L,
                       maxTime = 2))      # 2 s budget
  elapsed <- proc.time()["elapsed"] - start
  expect_lt(elapsed, 4)                   # not 2 * maxTime
  expect_equal(res$stop_reason, "max_time")
})
```

The existing test at `test-stopping.R:352` accepts either `converged` or
`max_time` and does not time the run — it would have passed even with the
current bug.
