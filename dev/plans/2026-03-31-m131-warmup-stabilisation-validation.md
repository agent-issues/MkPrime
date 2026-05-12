# M-131: Warmup Stabilisation Validation Study & Adaptive windowSize

## Summary

Implement adaptive `windowSize` scaling in `.CheckStabilisation()` and validate
empirically on Hamilton HPC. The existing briefing
(`.positai/plans/2026-03-29-m131-warmup-stabilisation-validation.md`) is the
authoritative design document; this plan covers the concrete implementation
steps.

---

## Step 1: Expose `logPostHistory` in the result object

**Goal:** Make the warmup logP snapshot trace available in the `MkPosterior`
object so Hamilton validation runs can save it in an RDS.

### Changes

**`R/MkPosterior.R`** — Add `warmup_trace` field:

```r
MkPosterior <- function(samples, trees, acceptance, model, data, mcmc,
                        warmup, tuning, warmup_trace = NULL) {
  structure(
    list(
      samples = samples, trees = trees, acceptance = acceptance,
      model = model, data = data, mcmc = mcmc,
      warmup = warmup, tuning = tuning,
      warmup_trace = warmup_trace
    ),
    class = "MkPosterior"
  )
}
```

**`R/RunMkPrime.R`** — At `.RunMkPrimeSingleRun()` return, save the trace:

The `logPostHistory` vector is already persisted in checkpoint state
(`r$logPostHistory`, line 1005). We need to pass it through to
`MkPosterior()`. Find the 3 `MkPosterior(...)` construction sites
(lines ~2007, ~2040, ~2060) and add `warmup_trace = runs[[1]]$logPostHistory`
(or equivalent per-run traces).

For the single-run codepath, the state `r$logPostHistory` is available
when `.RunMkPrimeSingleRun()` returns (it's in the checkpoint state).
Capture it from the returned state and pass to MkPosterior.

### Tests

- Verify `posterior$warmup_trace` is a numeric vector of length > 0 after
  a short `RunMkPrime()` call.

---

## Step 2: Write Hamilton validation scripts

### `inst/hamilton/m131-warmup-validation.R`

Single-dataset, single-seed validation run script. Called from SLURM with
two env vars: `DATASET` and `SEED`.

```r
# inst/hamilton/m131-warmup-validation.R
# Usage: DATASET=Sun2018 SEED=1234 Rscript inst/hamilton/m131-warmup-validation.R

library(MkPrime)
library(TreeTools)

dataset_name <- Sys.getenv("DATASET")
seed         <- as.integer(Sys.getenv("SEED"))
nex_file     <- file.path("data-raw", paste0(dataset_name, ".nex"))

dat <- ReadCharacters(nex_file)
mkd <- MkPrimeData(dat)

mcmc_val <- MkPrimeMCMC(
  nIter       = 200000L,
  thin        = 50L,
  maxWarmup   = 200000L,
  minWarmup   = 200000L,
  autoTune    = FALSE,
  nRuns       = 1L,
  nChains     = 1L
)

set.seed(seed)
posterior <- RunMkPrime(mkd, mcmc = mcmc_val)

result <- list(
  dataset       = dataset_name,
  seed          = seed,
  nTip          = length(mkd$tipLabels),
  nEdge         = 2L * length(mkd$tipLabels) - 3L,
  nChar         = sum(lengths(mkd$partitions)),
  warmup_trace  = posterior$warmup_trace,
  wall_time     = as.numeric(difftime(Sys.time(), .start, units = "secs"))
)

outdir <- Sys.getenv("OUTDIR", "results")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
saveRDS(result, file.path(outdir, sprintf("m131_%s_seed%d.rds",
                                           dataset_name, seed)))
```

### `inst/hamilton/m131-warmup-validation.slurm`

SLURM array job. 8 datasets × 4 seeds = 32 tasks.

```bash
#!/bin/bash
#SBATCH --job-name=m131-warmup
#SBATCH --partition=shared
#SBATCH --array=0-31
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=04:00:00
#SBATCH --output=logs/m131_%A_%a.out

DATASETS=(Longrich2010 Vinther2008 DeAssis2011 Wortley2006
          Sun2018 Wilson2003 Zhu2013 Dikow2009)
SEEDS=(2847 5193 7621 9034)

IDX=$SLURM_ARRAY_TASK_ID
DS_IDX=$((IDX / 4))
SEED_IDX=$((IDX % 4))

export DATASET=${DATASETS[$DS_IDX]}
export SEED=${SEEDS[$SEED_IDX]}
export OUTDIR=results

module load R/4.4.1

Rscript inst/hamilton/m131-warmup-validation.R
```

### Upload & submit

Use the `hamilton-hpc` skill workflow:
1. Upload the package tarball + validation scripts via `ssh::scp_upload()`
2. Install MkPrime on Hamilton (`R CMD INSTALL`)
3. Upload the relevant `.nex` files from `../TreeSearch-a/data-raw/`
4. Submit with `sbatch`
5. Record the job in `remote-jobs.md`

---

## Step 3: Post-hoc analysis script

### `inst/hamilton/m131-analyse-warmup.R`

After downloading the 32 RDS files, replay the stabilisation detector
offline with a grid of `(windowSize, nStableRequired)` configurations.

Key design points:

- **Detector replay function:** Extract `.CheckStabilisation()` logic into
  a standalone function that takes a logP trace and returns the iteration
  at which it would have declared stability for a given
  `(windowSize, zThreshold, nStableRequired)`.

- **Grid:**
  ```r
  windowSizes    <- c(3, 5, 7, 10, 12, 15, 20)
  nStableOptions <- c(2, 3, 4, 5)
  ```

- **Ground-truth stationarity:** Use the last 25% of each 200k trace as
  "equilibrium reference." Compare post-warmup mean logP to this reference.
  A detection is "too early" if the post-detection mean is > 2 logP units
  below equilibrium mean.

- **Adaptive formula replay:** For each trace, also replay with the
  proposed `windowSize = max(5, min(20, nEdge %/% 10))` and compare against
  the fixed `windowSize = 10` baseline.

- **Output:** Summary table (dataset × seed × config → iteration detected,
  classification: too-early / on-time / too-late), diagnostic plots.

---

## Step 4: Implement the validated formula

After analysing the Hamilton results (Step 3), update the call site in
`R/RunMkPrime.R` (line ~1007):

```r
# Current:
stabResult <- .CheckStabilisation(
  logPostHistory, nStableConsecutive
)

# After:
stabResult <- .CheckStabilisation(
  logPostHistory, nStableConsecutive,
  windowSize = max(5L, min(20L, nEdge %/% 10L))
)
```

The exact formula may be adjusted based on validation results.

---

## Step 5: Tests

- Existing warmup tests in `test-burnin.R` and `test-stopping.R` should
  pass unchanged (they use small fixed warmup / `minWarmup`).
- New test: verify `windowSize` scales correctly with `nEdge` for a range
  of tree sizes (unit test of the formula, no MCMC needed).
- New test: verify `posterior$warmup_trace` is populated.

---

## Execution order

Steps 1 and 2 can be done now (code changes + Hamilton submission).
Step 3 happens after Hamilton results return (~1-4 hours walltime).
Steps 4 and 5 depend on Step 3 analysis.

**This session:** Implement Steps 1–2, submit Hamilton job, record in
`remote-jobs.md`. Steps 3–5 will be picked up in a later session when
results are ready.

---

## Risk notes

- `logPostHistory` storage is negligible (≤400 doubles = 3.2 KB per run).
- The Hamilton runs use `nChains = 1` (no tempering) to isolate
  single-chain warmup behaviour. This is intentional.
- `autoTune = FALSE` keeps the chain in warmup phase for the full 200k
  iterations, ensuring the full logP trace is captured.
