# Three-Phase MCMC: Warmup / Tuning / Sampling

**Goal:** Replace the fixed warmup→sampling boundary with an
automatic three-phase architecture: **Warmup** (burn-in + tuning
adaptation), **Tuning** (ESS-aware move-weight optimisation), and
**Sample** (frozen weights, convergence-targeted sampling). Progress
display shows `Warmup`, `Tuning`, or `Sample` as the first word on
the ticker line (6 chars each, matching existing layout).

## Motivation

Runs terminate on a minESS / maxPSRF threshold, so the relevant
objective is *maximising min-ESS across parameters*. The current
warmup bandit scores moves by `acceptance_rate × dim / cost` — a
reasonable warm-start heuristic, but blind to actual ESS. A dedicated
tuning phase, run after the chain has stabilised, can optimise move
weights against the true objective (min-ESS/s) using stationary
samples.

Asking the user to set a warmup length is also undesirable: too short
risks non-stationarity, too long wastes compute. Auto-detection of
stabilisation removes this guesswork.

---

## Design

### Phase 1: Warmup

**What happens now, plus auto-termination.**

- Tuning parameters (scale, window, β-simplex concentration) adapt
  toward target acceptance rates — no change.
- Move weights adapt via the existing acceptance-rate/cost softmax
  heuristic — no change.
- Temperature ladder adapts (parallel tempering) — no change.
- **New:** A stabilisation detector monitors the cold-chain
  log-posterior and ends warmup automatically.

#### Stabilisation detector

Every `warmupBatch` (500 iter), after at least `minWarmup` iterations
(default 2000), compare the mean log-posterior in the most recent
window to the previous window using a Geweke-style z-score:

```
z = (mean_recent - mean_prev) / sqrt(var_recent/n_recent + var_prev/n_prev)
```

Declare stable when `|z| < 1.5` for `nStableRequired` consecutive
checks (default 3). This guards against false positives from
plateaus.

**Fallback:** A `maxWarmup` ceiling (default 50 000) hard-stops
warmup if stabilisation is never detected, with a cli warning.

**Edge case — parallel tempering:** Use the cold chain (β = 1) only
for stabilisation detection. Heated chains exist to help the cold
chain mix; their log-posteriors are not directly informative.

#### Display

```
Warmup  1234 │ logP: -567.8 │ minESS: ?
```

Denominator is omitted or shown as `~N` (estimated from
`maxWarmup`). The ticker already handles `NA` totals for `nIter =
Inf`.

### Phase 2: Tuning

**New phase.** Chain is approximately stationary; optimise move
weights via a min-ESS/s bandit.

At warmup→tuning transition:
1. Freeze tuning parameters (scale, window, etc.) — no further
   adaptation of step sizes.
2. Freeze temperature ladder.
3. Reset acceptance/timing counters (cumulative stats from warmup
   are stale w.r.t. the stationary regime).
4. Begin saving samples to a **tuning buffer** (circular, not
   written to the log file). These are discarded after tuning.

#### Bandit: perturbation-based min-ESS/s

The move pool is small (4–8 moves), so a simple perturbation scheme
suffices:

1. Run an **evaluation window** of `tuningWindowSize` iterations
   (default 500) under the current weight vector.
2. Compute per-parameter ESS on the tuning buffer using
   `coda::effectiveSize`. Take min-ESS. Divide by wall-time for
   the window → **min-ESS/s**.
3. Generate `nPerturbations` (default 3) candidate weight vectors
   by perturbing the current weights: pick a random free (unpinned)
   move, shift its weight by ±δ (drawn from a small distribution),
   renormalise.
4. Run an evaluation window under each candidate. Measure
   min-ESS/s.
5. Adopt the best-performing weight vector (current or candidate).
6. Repeat for `nTuningRounds` (default 5) or until min-ESS/s
   improvement < ε for 2 consecutive rounds (early stopping).

**Why perturbation, not a proper bandit (EXP3, Thompson)?** The
search space is continuous (weight simplex) and low-dimensional.
Perturbation with greedy selection is simpler to implement, debug,
and reason about. The evaluation windows are long enough (500 iter)
that min-ESS/s estimates have acceptable variance for comparative
ranking, even if individual ESS values are noisy.

**Interactions are captured automatically:** Because min-ESS/s is a
system-level outcome, it reflects indirect effects (e.g. better
topology mixing improving branch-length ESS). No per-parameter move
attribution is needed.

#### Tuning budget

Default: `nTuningRounds × (1 + nPerturbations) × tuningWindowSize`
= 5 × 4 × 500 = 10 000 iterations. This is configurable via
`tuningBudget` (a hard cap) in `MkPrimeMCMC()`.

**Tuning samples are discarded.** They are not written to the log
file or included in the posterior. This avoids contaminating the
posterior with samples drawn under shifting move weights (which
technically breaks detailed balance, even though the stationary
distribution is the same — the mixing properties change, so ESS
interpretation would be muddied).

#### Display

```
Tuning  1234 │ logP: -567.8 │ minESS/s: 0.42
```

The right-hand ticker shows `minESS/s` (the optimisation target)
during this phase rather than cumulative `minESS`, since the tuning
samples are discarded. Could also show `round 2/5` or similar.

### Phase 3: Sample

**Current post-warmup behaviour, unchanged except for the entry
point.**

At tuning→sample transition:
1. Freeze move weights (already done by the bandit).
2. Reset the tuning buffer. Begin writing samples to the log file
   / in-memory storage.
3. Log the frozen weights to the log file (existing
   `.LogMoveWeights()`).
4. Convergence checks proceed as today: ESS + PSRF at `checkEvery`
   intervals.

#### Display

```
Sample  5678 │ logP: -567.8 │ minESS: 234  PSRF: 1.02
```

Identical to the current post-warmup display, with `Sample` instead
of `iter`.

---

## `MkPrimeMCMC()` parameter changes

| Parameter | Current | Proposed | Notes |
|-----------|---------|----------|-------|
| `warmup` | Integer (default `nIter/2` or 5000) | **Deprecated but still accepted.** If supplied, used as `maxWarmup`; auto-detection still runs. `NULL` (default) → `maxWarmup = 50000`. | Backward-compatible. |
| `minWarmup` | — | New. Default 2000. Minimum iterations before stabilisation checks begin. | |
| `maxWarmup` | — | New. Default 50 000. Hard ceiling on warmup. | |
| `tuningBudget` | — | New. Default 10 000. Max iterations in the tuning phase. | |
| `tuningRounds` | — | New. Default 5. Number of perturbation rounds. | |

The `warmup` parameter is retained for backward compatibility: if a
user passes `warmup = 3000`, we set `maxWarmup = 3000` and
`minWarmup = min(minWarmup, maxWarmup / 2)`. A deprecation message
is emitted.

Users who want to skip automatic warmup detection entirely can pass
`warmup = 0` to go straight to tuning (or `tuning = FALSE` to skip
tuning too, preserving current behaviour for debugging).

---

## Implementation plan

### Step 1: Stabilisation detector

Add `.CheckStabilisation()` — a pure function taking a numeric
vector of log-posterior values (one per batch endpoint, accumulated
during warmup) and returning `TRUE`/`FALSE`.

```r
.CheckStabilisation <- function(logPostHistory, windowSize = 10L,
                                 zThreshold = 1.5,
                                 nStableRequired = 3L)
```

Unit-testable in isolation with synthetic log-posterior trajectories
(trending, plateau-then-jump, already-stationary).

**Tests:**
- Trending series → not stable.
- Stationary series → stable after `minWarmup`.
- Plateau then jump → resets stability counter.

### Step 2: Phase state machine in `.RunMkPrimeSingleRun()`

Replace the binary `warmup`/`iter` phase label with a three-state
machine:

```r
phase <- "Warmup"   # or "Tuning" or "Sample"
```

Current code structure:

```
if (batchEnd <= mcmc$warmup)  → adapt
if (batchEnd > mcmc$warmup && !weightsLogged)  → freeze weights
```

New structure:

```
if (phase == "Warmup") {
  # Adapt tuning, temperatures, move weights (acceptance heuristic)
  # Check stabilisation → transition to "Tuning"
}
if (phase == "Tuning") {
  # Run bandit evaluation windows
  # On completion → transition to "Sample"
}
if (phase == "Sample") {
  # Save samples, check convergence (existing code)
}
```

The batch-size logic also changes:
- Warmup: 500 (unchanged)
- Tuning: `tuningWindowSize` (500, aligned to evaluation rounds)
- Sample: 5000 (unchanged)

### Step 3: Tuning-phase bandit

Add `.RunTuningRound()` — takes current weights, runs an evaluation
window, returns min-ESS/s. The outer loop in the main batch loop
manages rounds and perturbations.

Tuning buffer: a small circular matrix (like the convergence window)
that accumulates cold-chain samples during tuning. Reused from the
existing streaming infrastructure (`.AddToStreamBuffer()` pattern)
but not flushed to disk.

### Step 4: Progress display updates

- `phaseLabel` becomes `phase` (already `"Warmup"` / `"Tuning"` /
  `"Sample"`).
- `.BuildProgressInfo()` gains a `phase` field (replacing
  `inWarmup`).
- `MkpTracePlot` and `MkpProgressJson` updated to use `phase`
  instead of `inWarmup`. The `inWarmup` field is retained for
  backward compat (TRUE when `phase == "Warmup"`).
- During tuning, the ticker shows `minESS/s` as the optimisation
  target.

### Step 5: `MkPrimeMCMC()` parameter updates

- Add `minWarmup`, `maxWarmup`, `tuningBudget`, `tuningRounds`.
- Deprecate `warmup` (accept with message, map to `maxWarmup`).
- Validate: `minWarmup <= maxWarmup`, `tuningBudget > 0`.
- Add `tuning = TRUE` flag to allow disabling the tuning phase
  entirely.

*Name collision:* The existing `tuning` parameter is a list of
step-size tuning values. The new "skip tuning phase" toggle needs a
different name. Options: `autoTune` (default `TRUE`), or
`moveTuning` (default `TRUE`). I'd suggest `autoTune` since `tuning`
is already taken.

### Step 6: Checkpoint compatibility

Checkpoints must record the current phase so that a resumed run
re-enters the correct phase. Add a `phase` field to the checkpoint
payload. Old checkpoints without `phase` default to the current
warmup/sampling logic (backward compat).

### Step 7: Documentation and tests

- Update `MkPrimeMCMC()` roxygen (section "Adaptive move
  scheduling" → "Three-phase warmup / tuning / sampling").
- Add integration test: run a short analysis with `minWarmup = 200`,
  `maxWarmup = 500`, `tuningBudget = 1000`, verify all three phases
  execute and the final weights differ from initial.
- Test stabilisation detector edge cases.
- Test checkpoint round-trip through all three phases.

---

## What this does NOT change

- The C++ hot loop (`run_mcmc_batch_cpp`) — no changes needed. It
  already accepts move weights and doesn't know about phases.
- The convergence check machinery (`.CheckConvergence()`) — used
  only in the Sample phase, unchanged.
- The parallel orchestrator — it polls log files for convergence.
  The warmup/tuning phases happen inside each worker. The
  orchestrator only cares about the Sample phase.
- Streaming to log files — samples are only streamed during the
  Sample phase, same as today (warmup samples are not logged).

---

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| Stabilisation detector false-positives (plateau before mode jump) | `nStableRequired = 3` consecutive passes; `minWarmup` floor; parallel tempering helps find modes faster. |
| Tuning phase too long on easy problems | Early stopping when improvement < ε; hard `tuningBudget` cap. |
| Tuning phase too short on hard problems | User can increase `tuningBudget`. Default of 10k iters is 5 rounds × 4 candidates × 500 = adequate for 4–8 moves. |
| Noisy ESS estimates in short windows | 500-iter windows produce rough ESS estimates, but we only need *comparative ranking* (is A better than B?), not precise values. |
| Backward compat breakage | `warmup` still accepted (mapped to `maxWarmup`). Old checkpoints handled via default phase. `inWarmup` retained in progress info. |

---

## Open questions

1. **Should tuning samples contribute to the posterior?** Current
   design discards them. Alternative: keep them, since the chain is
   stationary and the target distribution hasn't changed. The move
   weights shift, but that only affects mixing efficiency, not the
   target. Discarding is conservative; keeping would reduce total
   runtime. Leaning toward discard for simplicity — revisit if
   tuning budget becomes a significant fraction of total runtime.

2. **Perturbation magnitude.** The δ for weight perturbations needs
   tuning itself. Start with δ ~ Uniform(0.02, 0.10) per move,
   which produces modest shifts. Could adapt δ based on
   improvement magnitudes.

3. **Multi-run coordination.** In parallel mode, each run does its
   own warmup/tuning independently. Should they share tuning
   results? Probably not in v1 — each run may find different local
   optima, and independent tuning is simpler. Could revisit for v2.
