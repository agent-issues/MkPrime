# Plan: M-120 — 2D Joint Bactrian Proposals for Correlated Parameter Pairs

**Task:** M-120 from `to-do.md`

## Motivation

MkPrime's scalar parameters (tree_length, rate_loss, rate_log_sd, rate_neo,
beta_scale) are currently proposed one-at-a-time via independent 1D Bactrian
scale moves. When two parameters are posteriorly correlated—e.g., tree_length
and rate_log_sd often trade off because increasing ACRV variance requires
compensating tree length—independent proposals explore the correlated ridge
inefficiently.

The existing `neo_joint` move (moveType 18) already handles rate_loss ×
rate_neo by scaling both identically, but that's a 1D move along the
diagonal. A true 2D Bactrian proposes correlated perturbations whose
correlation matches the posterior, moving along the ridge rather than across
it.

`research-mcmc-optimizations.md` §3 recommends this as a lightweight
alternative to full AVMVN (which isn't justified for MkPrime's ~5 scalar
parameters). Potential benefit: up to 2–5× ESS improvement on correlated
pairs, at near-zero computational cost.

## Design

### 2D Bactrian kernel

Generate correlated perturbations `(z₁, z₂)` with the bimodal structure of
the 1D Bactrian, extended to 2D:

```
bactrian_2d(ρ) → (z₁, z₂):
  M  = 0.95                 // standard Bactrian mode
  σ  = sqrt(1 - M²)         // ≈ 0.312
  s  = random sign ±1        // both components go to same mode

  // Correlated Gaussian noise via Cholesky
  n₁ ~ N(0, 1)
  n₂ ~ N(0, 1)
  e₁ = σ * n₁
  e₂ = σ * (ρ * n₁ + sqrt(1 - ρ²) * n₂)

  z₁ = (s * M + e₁) * BACTRIAN_SCALE
  z₂ = (s * M + e₂) * BACTRIAN_SCALE
```

where `BACTRIAN_SCALE = 1/√12` (M-118 normalization).

Properties:
- Bimodal: modes at `(±M, ±M) × BACTRIAN_SCALE ≈ (±0.274, ±0.274)`
- Marginally each z has the same distribution as the 1D Bactrian
- Correlation between z₁ and z₂ ≈ ρ (not exact due to the bimodal
  structure, but close for |ρ| < 0.95)
- Symmetric: P(z₁, z₂) = P(-z₁, -z₂) → no Hastings correction beyond
  the Jacobian

### Move mechanics

For pair (param_A, param_B):

```
mult_A = exp(tuning * z₁)
mult_B = exp(tuning * z₂)
param_A_new = param_A_old * mult_A
param_B_new = param_B_old * mult_B
logHastings = log(mult_A) + log(mult_B) = tuning * (z₁ + z₂)
```

Uses a single `scaleTuning` value adapted via the standard acceptance-rate
targeting. The individual scale moves already handle per-parameter
exploration; the joint move's value is the *correlation*.

### Target pairs

| Pair | When | Rationale |
|------|------|-----------|
| tree_length × rate_log_sd | Always | ACRV variance and tree length trade off |
| tree_length × rate_loss | Neo chars present | Loss rate and tree length trade off |

The `rate_loss × rate_neo` pair is already handled by `neo_joint` (moveType
18) which scales both identically—appropriate since they're in the same model
component and are roughly proportional.

### Adaptive ρ estimation

During warmup, estimate ρ from the chain's posterior samples:

1. After the first ~200 accepted iterations, compute Pearson correlation
   of `(log tree_length, log rate_log_sd)` from the samples matrix
   returned by `run_mcmc_batch_cpp`.
2. Update ρ at each adaptation step (every `checkEvery` iterations during
   warmup).
3. Cap |ρ| at 0.95 for numerical stability.
4. Start with ρ = 0 (independent proposals) — the move degenerates to
   two independent 1D Bactrians, which is harmless.
5. Freeze ρ when tuning freezes (end of tuning phase), along with step
   sizes.

Storage: per-chain, per-joint-move. Passed to C++ via a new
`NumericMatrix jointRhos` (nChains × nMoves) parameter to
`run_mcmc_batch_cpp`. Non-joint moves have rho = 0 (ignored).

### Partial likelihood cache

Both tree_length and rate_log_sd affect all partitions, so the joint move
requires full likelihood recomputation (the `default` branch of the PLC
switch). No new PLC logic needed.

For `tree_length × rate_loss`: tree_length affects all partitions, so
full recomputation regardless.

## Files to modify

### C++ (`src/mcmc.cpp`)

1. **New kernel function** `bactrian_2d_perturbation(double rho, double& z1, double& z2)`:
   - Inline, next to existing `bactrian_perturbation()`
   - Cholesky decomposition for correlated noise

2. **Exported test helper** `bactrian_2d_draws(int n, double rho)`:
   - Returns Nx2 matrix for property tests

3. **New move type codes** in `do_move_impl`:
   - `case 21:` joint_tl_rls (tree_length × rate_log_sd)
   - `case 22:` joint_tl_rl (tree_length × rate_loss)
   - Each reads `jointRho` from a new per-move array

4. **`run_mcmc_batch_cpp` signature**: add `NumericMatrix jointRhos`
   (nChains × nMoves). Read `jointRhos(ch, moveIdx)` for joint moves,
   ignored for all others.

5. **Rollback**: Both parameters are already snapshotted at the top of
   `do_move_impl` (oldTL, oldRL, oldRLSD), so rollback is automatic.

6. **PLC switch**: Joint moves fall through to `default` (full
   recomputation). No changes needed.

### R side

7. **`R/MkPrimeMCMC.R`**:
   - New parameter `joint2d = TRUE` (default on).
   - Validation: logical.
   - New tuning defaults: `scale_joint_tl_rls = 0.5`,
     `scale_joint_tl_rl = 0.5`.

8. **`R/RunMkPrime.R` — `.BuildMoves()`**:
   - Add `joint_tl_rls` move (always, weight = 1, dim = 1L).
   - Add `joint_tl_rl` move (only when `hasNeo`, weight = 1, dim = 1L).

9. **`R/RunMkPrime.R` — `.kMoveTypes`**:
   - `joint_tl_rls = 21L`, `joint_tl_rl = 22L`.

10. **`R/RunMkPrime.R` — `.BuildScaleTuningMatrix()`**:
    - Add entries for joint_tl_rls, joint_tl_rl reading from tuning list.

11. **`R/RunMkPrime.R` — New `.BuildJointRhoMatrix()`**:
    - Returns nChains × nMoves matrix; 0.0 for non-joint moves.
    - Called before each batch; updated during warmup adaptation.

12. **`R/RunMkPrime.R` — `.AdaptTuning()`**:
    - Add target acceptance rate for joint moves (0.25 — lower than 1D
      because 2D proposals are harder to accept).
    - Add tuning keys: `scale_joint_tl_rls`, `scale_joint_tl_rl`.

13. **`R/RunMkPrime.R` — Warmup adaptation loop**:
    - After each batch during warmup, compute ρ from cold chain's
      recent samples. Update the chain's rho values.
    - Use `cor(log(tree_length), log(rate_log_sd))` from the samples
      matrix (available from `run_mcmc_batch_cpp` return).
    - Store in `r$chain_rhos[[ch]]` (named list per chain).

14. **`R/RunMkPrime.R` — `run_mcmc_batch_cpp` call site**:
    - Build and pass `jointRhos` matrix.

15. **`R/RunMkPrime.R` — `ResumeMkPrime()` and `mkp_stepping_stone()`**:
    - Same treatment as other move types.

### Tests

16. **`tests/testthat/test-joint-2d.R`** (new):
    - 2D kernel properties: marginal distributions match 1D Bactrian,
      correlation ≈ ρ (within statistical tolerance), symmetry.
    - Joint move integration test: short run with joint2d = TRUE
      produces valid posterior, acceptance > 0.

17. **`tests/testthat/test-m092-adaptive-scheduler.R`**:
    - Add "joint_tl_rls" and "joint_tl_rl" to valid move names.

### Documentation

18. **Roxygen**: `@param joint2d` in `MkPrimeMCMC()`, mention of the
    adaptive correlation.

## What stays the same

- 1D Bactrian scale moves remain unchanged — joint moves supplement,
  not replace.
- `neo_joint` (moveType 18) unchanged.
- PLC logic unchanged (joint moves hit the default branch).
- Checkpoint/resume: rho values are derived from samples during warmup,
  so they're re-estimated on resume. No need to checkpoint ρ separately
  (but the tuning step sizes are checkpointed as usual).

## Risks and mitigations

| Risk | Mitigation |
|------|------------|
| ρ estimate noisy early in warmup | Start with ρ=0; require ≥200 samples before estimating |
| Joint acceptance too low on small trees | Single scaleTuning adapts down; weight adapts via scheduler |
| Overhead of ρ estimation | O(nSamples) per adaptation step; negligible vs likelihood cost |
| Interaction with adaptive scheduler | Joint moves have separate weights; scheduler can down-weight if unhelpful |

## Estimated scope

~200 lines C++ (kernel + two case blocks + batch parameter), ~100 lines R
(move registration + ρ estimation + tuning), ~80 lines tests.
Medium effort, straightforward extension of existing infrastructure.
