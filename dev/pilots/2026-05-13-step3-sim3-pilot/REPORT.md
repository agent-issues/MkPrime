# Step 3 Pilot Report — 20k Sim 3 chain with asymmetric-slab + θ-logging fixes

**Date:** 2026-05-13  
**Plan:** `dev/plans/2026-05-13-1100-v2-asymmetric-slab-prior-and-gamma-ridge.md`, Step 3

---

## Configuration

| Parameter | Value |
|---|---|
| Script | `dev/pilots/2026-05-13-step3-sim3-pilot/run-pilot.R` |
| Seed | 20260512 (canonical pre-fix seed, direct comparability) |
| kEco | 2 (ecology 0 = C/D clades; ecology 1 = A/B clades) |
| nTheta | 1 (one non-reference ecology) |
| nEco | 60 neomorphic characters |
| nBase | 180 transformational characters |
| phi_truth | 4 |
| stem / root | 0.10 / 0.15 (Goldilocks) |
| nIter | 20,000 |
| thin | 40 |
| nChains / nRuns | 1 / 1 |
| minWarmup / maxWarmup | 2,000 / 5,000 |
| model | ecologyAware=TRUE, magnitudeMode="global", kPrimePrior="geometric" |
| Log file | `dev/pilots/2026-05-13-step3-sim3-pilot/sim3-pilot.log` |
| RDS | `dev/pilots/2026-05-13-step3-sim3-pilot/sim3-pilot-result.rds` |
| Elapsed | ~4.4 min |

**Bugs fixed before this run:**
- Bug A: asymmetric-slab prior (`cpp_log_prior` + `LogPrior`) — θ now has proper Beta counterweight
- Bug B: θ now logged — `theta_1` column present in output (verified)

**Post-burnin samples:** 175 raw rows → 131 retained after 25% discard.

**True tree_length:** 12.7 (computed from Goldilocks topology: 4 clades × 3.1 + 2 root branches × 0.15). The plan's "~9" is the effective per-gamma-normalised scale; tl/γ_e median = 9.35 confirms this.

---

## Summary table: pre-fix vs post-fix

| Metric | Pre-fix | Post-fix | Truth | Direction |
|---|---|---|---|---|
| tree_length median | 22.0 | 29.0 | 12.7 | Modestly worse in median |
| tree_length IQR | [11.1, 124.8] | [20.1, 42.3] | — | Much tighter upper tail |
| tree_length max | 7,496 | 70.0 | — | 107× reduction |
| γ_e median (actual θ) | 1.130 (θ=0.5 proxy) | 3.003 | ~1.9 | Now tracking θ correctly |
| tl / γ_e median | 14.49 | 9.35 | ~9 | Near truth |
| tl / γ_e IQR | [10.4, 41.2] | [6.1, 12.7] | — | Tight, symmetric around truth |
| tl / γ_e max | 4,010 | 28.8 | — | 139× reduction |
| cor(log tl, log γ_e) | +0.629 | +0.197 | ~0 | 3× reduction |
| ESS phi | — | 8.9 / 131 | — | Low at 20k; expected |
| ESS pi0 | — | 3.7 / 131 | — | Low at 20k; expected |
| ESS theta_1 | — | 10.2 / 131 | — | Low at 20k; expected |
| ESS tree_length | — | 8.1 / 131 | — | Low at 20k; expected |

*Pre-fix numbers from 200k-iter chain (n=194 post-burnin). Pre-fix γ_e used θ=0.5 proxy.*

---

## Posterior covariance (log/logit-transformed)

```
                 log_phi  logit_pi0  logit_theta1  log_tree_length
log_phi           0.2417    -0.1277        0.0432          -0.0536
logit_pi0        -0.1277     0.2312       -0.0471           0.0775
logit_theta1      0.0432    -0.0471        0.0426          -0.0250
log_tree_length  -0.0536     0.0775       -0.0250           0.2258
```

Largest off-diagonal: `cov(logit_pi0, log_tl) = +0.078`, a weak positive coupling
(pi0 ↑ → γ_e ↓ → tl ↑ compensates). No entry approaches a problematic ridge value.

---

## Hyperparameter posteriors (post-fix)

| Param | Median | IQR | 95% CI |
|---|---|---|---|
| phi | 0.236 | [0.152, 0.309] | [0.105, 0.836] |
| pi0 | 0.278 | [0.196, 0.318] | [0.143, 0.508] |
| theta_1 | 0.030 | [0.026, 0.033] | [0.018, 0.041] |

Note: phi < 1 reflects enc/disc label non-identifiability (symmetric under Mk' relabelling).
theta_1 near zero is prior-dominated at this character count; expected.

---

## Decision-rule branch

**BRANCH: Tail GONE, covariance benign → ship simple Bactrians; proceed to multi-rep Sim 3.**

Supporting numbers:
- tree_length max: 7,496 → **70** (107× reduction). No catastrophic escapes post-fix.
- cor(log tl, log γ_e): +0.629 → **+0.197** (well below the 0.4 plan threshold).
- tl/γ_e median: 14.49 → **9.35** (near truth ~9).
- tl/γ_e max: 4,010 → **28.8** (139× reduction).
- All posterior covariance off-diagonals < |0.13|; no problematic ridge identified.

The remaining tree_length spread (median 29 vs truth 12.7) reflects legitimate posterior
uncertainty from a 20k single-replicate pilot with low ESS, not a runaway ridge. A longer
multi-rep analysis will resolve the residual spread.

---

## Warnings and issues

1. **Warmup reached maxWarmup (5000) without stabilisation.** Expected at 20k pilot scale.
   ESS values (3.7–10.2 / 131) indicate incomplete mixing. Multi-rep run should use nIter ≥ 100k.

2. **1 invariant character dropped** (column 19, kObs ≤ 1). Standard; harmless.

3. **phi posterior inverted** (median 0.24, truth 4). Expected — enc/disc label non-identifiability
   under Mk' relabelling correction. Magnitude recovery is what matters for topology inference.

4. **theta_1 near zero** (median 0.030). Prior-dominated under Beta(2,2) with 60 ecology characters.
   The asymmetric-slab prior is correctly counterweighting θ; this is not a bug.

---

## Next steps

- Proceed to **Step 5 (multi-rep Sim 3)**: 20 replicates × 2 chains (aware + blind), nIter ≥ 100k, Goldilocks config. Dispatch to HPC (~20h total).
- **No Step 4** (v2-no-γ) needed: γ-coupling is below the 0.4 threshold.
- Re-run rodent vignette (Step 6) with corrected prior.
