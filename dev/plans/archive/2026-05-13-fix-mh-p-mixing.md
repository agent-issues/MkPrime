# Plan: Diagnose & fix poor mixing under `empirical_geometric` prior

## Current Hamilton state (READ FIRST)

- **Branch on GitHub**: `feature/empirical-geometric-prior` (commits `9332a54` + `50a3c15`).
- **Hamilton repo**: `/nobackup/pjjg18/mkp-study/mkp` checked out on the feature branch.
- **Hamilton installed lib**: `/nobackup/pjjg18/mkp-study/lib` (rebuilt this session; new prior loads correctly — `MkPrimeModel()$kPrimePrior == "empirical_geometric"` confirmed).
- **Active SLURM array**: **`17135710`** (260 tasks, all RUNNING, hitting `maxTime=6h` without converging). **Cancel before re-submitting.**
- **Per-task data**: `tree_NN/rep_MM` under `/nobackup/pjjg18/mkprime-files/tree-inference`.
- **Per-task results dir**: `/nobackup/pjjg18/mkp-study/results/t<NN>_r<MM>/` (checkpoints live here).
- **Logs**: `/nobackup/pjjg18/mkp-study/logs/mkp_eg_17135710_<idx>.{out,err}`.
- **SSH**: must use Windows-native `/c/Windows/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk` (Git-Bash ssh can't see the Windows ssh-agent service). Key loaded via `Start-Service ssh-agent` + `ssh-add C:/Users/pjjg18/id_ed25519`.

## What happened

The Hamilton smoke test (15-min single task) ran cleanly: 232k iterations, no errors. But the full 260-task array shows:

- **logP swings ±30** every few thousand samples (e.g. `-349 → -314 → -340`) — chain bouncing between distant modes.
- **minESS = 6-54** after 395k samples — extreme autocorrelation.
- **ETA estimates jump 7h–85h** between adjacent windows — no convergence trajectory.

At this rate all 260 tasks will exhaust the 6h `maxTime` budget producing posteriors too noisy to compare against the `mkp` arm.

## Why the smoke test missed it

The smoke test (15 min, single task) checked:
1. Package installs ✓
2. `run_one.R` parses arguments ✓
3. MCMC engine doesn't crash and iterates ✓

It did **not** check:
- **Per-move acceptance rates** — `mh_p` could be at 0.1% acceptance and the smoke test wouldn't notice; it just reports `Sample N │ logP:X │ minESS: ?` and that's fine for a too-short chain.
- **Effective sample size growth rate** — 15 min is shorter than the warmup phase, so ESS is always `?`.
- **logP convergence** — short chains haven't moved far from initialisation.

**Lesson for next smoke test**: run for ≥ `2 × maxWarmup × thin` iterations (here ≥ 100k samples) and check (a) per-move acceptance rates, (b) minESS rate-of-rise, (c) logP variance after warmup. A 30–45 min single-task smoke would have caught this.

## Hypothesis (most likely root cause)

The `mh_p` move I added is multiplicative on `p ∈ (0, 1)` (reuses moveType 8 `scale_p`: `p_new = p · exp(σ · z)` where `z` is a Bactrian draw). When `p` is near 1, any positive perturbation pushes `p > 1` → prior boundary returns `-Inf` → rejected. Under the empirical-geometric prior the posterior on `p` typically concentrates high (large empirical mass at low `kObs` → small `u` → high `p`). So acceptance rate is likely <5% and `p` is effectively stuck. The Gibbs `kPrime` sweep then conditions on a frozen `p`, producing the long-range autocorrelation we see.

## Plan

### Step 1 — Stop wasting cluster time

```
ssh pjjg18@hamilton8.dur.ac.uk "scancel 17135710"
```

### Step 2 — Diagnose locally (confirm hypothesis)

Run a single-dataset analysis locally with the same `MkPrimeMCMC` config as the Hamilton job (but short, ~5–10 min) and inspect per-move acceptance rates from `result$acceptance`. Compare `mh_p` acceptance vs other moves. Expected finding: `mh_p` acceptance << 10%.

```r
# In R, in mkp/ project root
suppressMessages(devtools::load_all('.', quiet = TRUE))
library(ape)
# Use one of the Hamilton datasets, or simulate similarly
set.seed(1)
tree <- ape::rtree(20); tree$edge.length <- runif(nrow(tree$edge), 0.05, 0.25)
# simulate ~50 chars, kTrue mixed 2-5
# (use code from data-raw/smoke_empirical_geometric.R)
mkd <- MkPrimeData(pd)
res <- RunMkPrime(mkd, tree,
                   model = MkPrimeModel(kPrimePrior = "empirical_geometric"),
                   mcmc = MkPrimeMCMC(nIter = 20000L, thin = 10L,
                                       maxWarmup = 5000L, minWarmup = 5000L,
                                       autoTune = FALSE, nRuns = 2L,
                                       nChains = 2L))
print(res$acceptance)   # this is the diagnostic — look for mh_p row
```

If `mh_p` acceptance is low, hypothesis confirmed.

### Step 3 — Fix the `p` sampler

Replace the multiplicative MH move with one of these (in order of simplicity / robustness):

**Option A — Logit-scale MH (recommended).** Propose on the real line so the move never overshoots `(0, 1)`:

```
logit_p_new = logit(p) + σ · Bactrian
p_new = sigmoid(logit_p_new)
log_hastings = log(p_new · (1 - p_new)) - log(p · (1 - p))  # Jacobian
```

Implementation: add a new moveType code in `src/mcmc.cpp` for `mh_logit_p`, mirroring case 8 but on the logit scale. Update `R/RunMkPrime.R` move scheduling to use it under `empirical_geometric`.

**Option B — Slice sampler.** The package already has a slice-sampling framework (`slice_kprime_hyper`, `slice_rate_loss`). Add a `slice_p` variant. Robust, tuning-free, but more invasive.

**Option C — Proposal scale autotuning.** Keep the multiplicative move but engage the autotuning path (`autoTune = TRUE` in `MkPrimeMCMC`) so the proposal scale adapts down. Easiest, but `run_one.R` deliberately uses `autoTune = FALSE` to freeze move weights for reproducibility, so this is a config change with side effects.

**Recommended**: Option A. Smallest code change with the largest impact, and matches the existing pattern of using Bactrian perturbations.

### Step 4 — Re-smoke locally (longer)

Run a ~30-min smoke under the fix:

```r
# As above but nIter = 100000L (4× warmup), check:
#   - res$acceptance["mh_logit_p", ] > 0.20 ideally
#   - minESS rises monotonically with samples
#   - logP variance shrinks after warmup
```

Also keep `data-raw/smoke_empirical_geometric.R` updated to do this check by default in future.

### Step 5 — Update run_one.R if needed

If the move name changes (e.g. `mh_p` → `mh_logit_p`), no caller changes required — the move is scheduled internally by `.BuildMoves`. But verify `R/RunMkPrime.R` references are updated.

### Step 6 — Commit, push to feature branch, redeploy on Hamilton

```
git add ...
git commit -m "Fix p sampler under empirical_geometric (logit-MH)"
git push origin main:feature/empirical-geometric-prior
```

On Hamilton:
```
cd /nobackup/pjjg18/mkp-study/mkp && git pull --ff-only
module load r/4.5.1 gcc/14.2
rm -f src/*.o src/*.so
R CMD build --no-build-vignettes --no-manual --no-resave-data .
mv MkPrime_0.0.0.9000.tar.gz /nobackup/pjjg18/mkp-study/
cp data-raw/hamilton/run_one.R /nobackup/pjjg18/mkp-study/run_one.R
cd /nobackup/pjjg18/mkp-study && sbatch install_mkp.sh
# wait for install to finish (~1.5 min)
```

### Step 7 — Cleanup stale checkpoints from job 17135710

The previous failed array left per-task `t<NN>_r<MM>/mkp_eg_*.{rds,log,nwk}` files. Without cleanup, resubmitted tasks will try to resume bad checkpoints. Purge:

```
ssh pjjg18@hamilton8.dur.ac.uk "find /nobackup/pjjg18/mkp-study/results -name 'mkp_eg_*' -delete && find /nobackup/pjjg18/mkp-study/results -name '.slurm_job_id' -delete"
```

### Step 8 — One-task wallclock-realistic smoke on Hamilton

Don't relaunch all 260 immediately. Run **one task with full 6h `maxTime`** to confirm the fix works at scale:

```
ssh pjjg18@hamilton8.dur.ac.uk \
  "cd /nobackup/pjjg18/mkp-study && sbatch --array=0 --job-name=mkp-eg-test mkp_eg_array.slurm"
```

Watch `minESS`/logP in the .err log; expect minESS to climb past 200 within 1–2 hours if the fix worked.

### Step 9 — Re-submit the full array

If the single-task smoke converges, fire the full array:

```
ssh pjjg18@hamilton8.dur.ac.uk "cd /nobackup/pjjg18/mkp-study && sbatch mkp_eg_array.slurm"
```

## Stretch goal

Investigate the `[DIAG] fopen failed iter=...` spam every 500 iters — it predates this work (present in the existing MCMC engine) but is noise that pollutes diagnostic logs. Probably a leftover diagnostic file path that doesn't exist on Hamilton's filesystem layout. Out of scope for this fix but worth opening as a separate issue.

## Files in play

| Path | Role |
|---|---|
| `src/mcmc.cpp` | Add new `mh_logit_p` move type (case N) + helper for logit/sigmoid transform with Jacobian |
| `R/RunMkPrime.R` | Update `.BuildMoves` empirical_geometric branch + `.kMoveTypes` map + R-side `.DoMove` switch + tuning matrix |
| `data-raw/smoke_empirical_geometric.R` | Extend smoke to check `res$acceptance` and minESS rate of rise |
| `data-raw/hamilton/run_one.R` | No change expected (move scheduling is internal) |
| `tests/testthat/test-empirical-geometric-prior.R` | Add acceptance-rate test for the new move |
