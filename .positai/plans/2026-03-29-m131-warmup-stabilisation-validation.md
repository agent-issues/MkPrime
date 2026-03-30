# M-131: Warmup Stabilisation Detector — Validation Study & Adaptive Tuning

## Problem

The warmup stabilisation detector (`.CheckStabilisation()` in `R/RunMkPrime.R`,
line ~2703) uses a Geweke-style z-test comparing two windows of `windowSize`
consecutive log-posterior snapshots. Current defaults are all hardcoded:

| Parameter | Default | Effect |
|-----------|---------|--------|
| `windowSize` | 10 | Snapshots per Geweke window |
| `zThreshold` | 1.5 | |z| below which a check passes |
| `nStableRequired` | 3 | Consecutive passes to declare stable |
| (batch size) | 500 iter | Iterations between snapshots |

Since the detector needs `2 × windowSize` snapshots before its first check, and
then `nStableRequired` consecutive passes, the **minimum possible warmup** is:

    (2 × windowSize + nStableRequired) × batchSize
    = (20 + 3) × 500 = 11,500 iterations

On the hyoliths dataset (29 tips, 53 chars), the chain appeared visually
stationary by ~5,000 iterations but the detector didn't fire until ~19,000.
This suggests the detector is working (catching genuine slow drift), but the
11,500-iteration floor may be unnecessarily high for small/easy datasets, and
unnecessarily low for large/hard ones.

### Proposed fix: scale `windowSize` with tree size

Larger trees have longer autocorrelation times in logP, so each 500-iteration
snapshot is more correlated with its neighbours. A fixed `windowSize` gives
fewer effective samples per Geweke window on large trees and more than needed
on small trees. Scaling `windowSize` with `nEdge` addresses this directly:

```r
windowSize <- max(5L, min(20L, nEdge %/% 10L))
```

Implications by tree size:

| nTip | nEdge | windowSize | Min warmup | Window span |
|:----:|:-----:|:----------:|:----------:|:-----------:|
| 10 | 17 | 5 | 6,500 | 2,500 iter |
| 29 | 55 | 5 | 6,500 | 2,500 iter |
| 50 | 97 | 9 | 10,500 | 4,500 iter |
| 75 | 147 | 14 | 15,500 | 7,000 iter |
| 100 | 197 | 19 | 20,500 | 9,500 iter |
| 120+ | 237+ | 20 (cap) | 21,500 | 10,000 iter |

We may also want to adjust `nStableRequired` (e.g. increase it for small
`windowSize` to compensate for noisier individual checks).

**This must be validated empirically** — we need to confirm that the smaller
windows on small trees don't cause premature termination (false positives),
and that the larger windows on big trees aren't wastefully conservative.

---

## Validation Study Design

### Overview

Run MkPrime on ~8 benchmark datasets spanning a range of tree sizes, with
**warmup termination disabled** (very high `maxWarmup`), logging the logP
snapshot trace. Then **replay the detector offline** with many `windowSize`
configurations to find where each would have terminated warmup. Compare
against ground-truth stationarity (determined from the long traces).

This is efficient: the Hamilton runs need only happen once; all detector
tuning is done in post-hoc R analysis.

### Datasets

Selected from `../TreeSearch-a/data-raw/*.nex` to cover the nTip range.
All are real morphological matrices used in published phylogenetic analyses.

| Dataset | nTip | nChar | nEdge | Why included |
|---------|:----:|:-----:|:-----:|--------------|
| Longrich2010 | 20 | 93 | 37 | Smallest available |
| Vinther2008 | 23 | 57 | 43 | Small, few chars |
| DeAssis2011 | 33 | 50 | 63 | Medium-small |
| Wortley2006 | 37 | 105 | 71 | Medium, more chars |
| Sun2018 | 54 | 225 | 105 | Hyoliths-scale (our primary user dataset) |
| Wilson2003 | 61 | 165 | 119 | Medium-large |
| Zhu2013 | 75 | 253 | 147 | Large, many chars |
| Dikow2009 | 88 | 220 | 173 | Largest available |

### Hamilton job specification

**Per dataset:** 4 independent seeds, each running for a fixed 200,000
iterations with `maxWarmup = 200000L` (effectively no early termination).
The run should log the logP snapshot trace at every warmup batch endpoint.

**Modifications needed before the Hamilton run:**

1. `.CheckStabilisation()` currently does NOT save its logP trace to the
   output — it's a local variable (`logPostHistory`) that gets discarded.
   **Add a mechanism to export the logP snapshot vector** from the warmup
   phase. Options:
   - Easiest: save `logPostHistory` into the checkpoint/result object
     (e.g. `posterior$warmup_trace`).
   - Alternative: write it to a small CSV alongside the main log file.

2. The `RunMkPrime()` progress callback already receives `coldLogpost` at
   each batch. Confirm this is the same value appended to `logPostHistory`
   (line 689). If so, an alternative is to have the Hamilton script capture
   progress callbacks and log the trace externally. But modifying
   `RunMkPrime()` to save `logPostHistory` in the result is cleaner and
   useful beyond this study.

**MCMC configuration for the validation runs:**

```r
mcmc_val <- MkPrimeMCMC(
  nIter       = 200000L,
  thin        = 50L,
  maxWarmup   = 200000L,   # never terminate warmup
  minWarmup   = 200000L,   # never terminate warmup
  autoTune    = FALSE,     # stay in warmup phase throughout
  nRuns       = 1L,
  nChains     = 1L         # no tempering — isolate single-chain behaviour
)
```

**Output per run:** An RDS file containing:
- The `logPostHistory` vector (one logP per 500-iteration batch = 400 values)
- Dataset name, seed, nTip, nEdge, nChar
- Wall time

**Resource estimate:** 32 runs × ~30–90 min each (larger datasets slower) =
~16–48 core-hours. A single Hamilton node (40 cores) can run all 32 in
parallel in ~90 min wall time.

### Post-hoc analysis (local R, no Hamilton needed)

Load all 32 RDS files and, for each trace, replay the detector with a grid of
configurations:

```r
windowSizes     <- c(3, 5, 7, 10, 12, 15, 20)
nStableOptions  <- c(2, 3, 4, 5)
```

For each `(windowSize, nStableRequired)` combination:

1. **Find the iteration at which the detector would have declared stability.**
   This is a simple loop over `logPostHistory` calling `.CheckStabilisation()`.

2. **Assess ground-truth stationarity** for each trace:
   - Visual: plot the full 200k logP trace and mark the detector's chosen
     cutpoint for each configuration.
   - Quantitative: use the Heidelberger–Welch stationarity test (from `coda`)
     on the full trace as a reference. Or simply compare the mean logP from
     the last 50k iterations to the detector's chosen start-of-stationarity.

3. **Classify each (dataset × seed × config) as:**
   - **Too early:** detector fired before genuine stationarity (mean logP
     in post-warmup window is > 2 logP units below the final-50k mean)
   - **On time:** fired within a reasonable window after stationarity
   - **Too late:** fired > 2× the ground-truth stationarity iteration

4. **Summarise** in a table: for each `windowSize` formula (fixed or
   nEdge-scaled), fraction of runs that are too-early / on-time / too-late,
   and median overhead (extra warmup beyond ground truth).

### Key analysis questions

1. Does `windowSize = max(5, nEdge %/% 10)` produce **zero** "too early"
   results across all 32 runs? (Hard requirement — false positives are the
   main risk.)

2. What is the **median overhead** of the nEdge-scaled formula vs the
   current fixed `windowSize = 10`, stratified by tree size?

3. Would a different formula (e.g. `nEdge %/% 8`, or `max(5, nEdge %/% 12)`)
   be better?

4. Does increasing `nStableRequired` from 3 to 4 or 5 help protect small
   `windowSize` values from false positives without excessive overhead?

5. Is there a dataset where the detector fires much later than ground truth
   regardless of configuration? (Would indicate a fundamentally different
   problem, e.g. multi-modal logP requiring a different test.)

---

## Implementation Plan (after validation results)

### Step 1: Expose `logPostHistory` in the result object

In `RunMkPrime()` (R/RunMkPrime.R), after the warmup→tuning/sample
transition, save `logPostHistory` into the result. This is a small vector
(≤ 400 entries for 200k warmup at batch=500) and is useful diagnostically.

**File:** `R/RunMkPrime.R`
- After line ~704 (`"Chain stabilised at iteration {batchEnd}."`), or in the
  checkpoint-save logic, store `logPostHistory` in `r$logPostHistory`.
- In `MkPosterior` construction, pass it through so `posterior$warmup_trace`
  is available. (Or keep it internal — just needs to be in the RDS for the
  validation study.)

### Step 2: Hamilton validation runs

Write `inst/hamilton/m131-warmup-validation.R`:
- Loops over datasets × seeds
- Reads each `.nex` from a specified directory
- Runs `RunMkPrime()` with the validation MCMC config above
- Saves one RDS per run

Write `inst/hamilton/m131-warmup-validation.slurm`:
- SLURM array job, one task per (dataset × seed)
- Estimated: 32 tasks, 1 core each, 2 hours walltime, 4 GB RAM

### Step 3: Post-hoc analysis

Write `inst/hamilton/m131-analyse-warmup.R`:
- Loads all 32 RDS files
- Replays detector with the grid of configurations
- Produces summary table + diagnostic plots
- Outputs a recommended `windowSize` formula

### Step 4: Implement the chosen formula

In `.CheckStabilisation()`, replace the hardcoded `windowSize = 10L` default
with the validated formula. The function already takes `windowSize` as a
parameter, so the change is at the **call site** (line ~691):

```r
# Current:
stabResult <- .CheckStabilisation(logPostHistory, nStableConsecutive)

# After:
stabResult <- .CheckStabilisation(
  logPostHistory, nStableConsecutive,
  windowSize = max(5L, min(20L, nEdge %/% 10L))
)
```

If the validation study suggests a different formula or adjusted
`nStableRequired`, update accordingly.

### Step 5: Test

- Existing warmup tests (`test-burnin.R`, `test-stopping.R`) should still
  pass — they use small fixed warmup so the detector is irrelevant.
- Add a new test that verifies `windowSize` scales correctly with `nEdge`.
- Re-run the hyoliths vignette to confirm warmup ends earlier than 19k.

---

## Code locations

| What | Where |
|------|-------|
| `.CheckStabilisation()` | `R/RunMkPrime.R` line ~2703 |
| Call site (warmup loop) | `R/RunMkPrime.R` line ~691 |
| `logPostHistory` accumulation | `R/RunMkPrime.R` line ~689 |
| `warmupBatch` (500L) | `R/RunMkPrime.R` line ~432 |
| `minWarmup` / `maxWarmup` defaults | `R/MkPrimeMCMC.R` (search `minWarmup`) |
| Benchmark datasets | `../TreeSearch-a/data-raw/*.nex` |

## Risk assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| False-positive stability (premature warmup end) | Medium for small windowSize | High — biased posterior | Validation study; conservative `nStableRequired` |
| Excessive warmup on large trees | Low | Medium — wasted compute | Cap at windowSize=20; maxWarmup provides hard ceiling |
| Formula doesn't generalise beyond test datasets | Low | Medium | 8 diverse datasets + 4 seeds = 32 runs; real published matrices |
| logPostHistory storage bloats checkpoint files | Very low | Low | ≤ 400 doubles = 3.2 KB |
