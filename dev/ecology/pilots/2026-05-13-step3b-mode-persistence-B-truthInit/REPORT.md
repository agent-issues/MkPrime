# Step 3b Mode-Persistence Report — Chain B (truth-init, 100k)

**Date:** 2026-05-14
**Plan:** `dev/plans/2026-05-13-1100-v2-asymmetric-slab-prior-and-gamma-ridge.md`, Open question: label-switching

---

## Configuration

| Parameter | Value |
|---|---|
| Script | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/run.R` |
| Seed (data) | 20260512 (same as 20k pilot — identical dataset) |
| Seed (chain) | 20260514 (independent of chain A and 20k pilot) |
| kEco | 2 |
| nTheta | 1 |
| nEco | 60 neomorphic characters |
| nBase | 180 transformational characters |
| phi_truth | 4 |
| stem / root | 0.10 / 0.15 (Goldilocks) |
| **phi_init** | **4.0 (truth-mode override)** |
| **theta_init** | **0.97 (truth-mode override)** |
| nIter | 100,000 |
| thin | 40 |
| nChains / nRuns | 1 / 1 |
| minWarmup / maxWarmup | 2,000 / 5,000 |
| model | ecologyAware=TRUE, magnitudeMode="global", kPrimePrior="geometric" |
| Log file | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/chain.log` |
| RDS | `dev/pilots/2026-05-13-step3b-mode-persistence-B-truthInit/result.rds` |

**Post-burnin samples:** 2125 raw rows -> 1593 retained after 25% discard.

---

## 0. Init verification (critical)

The truth-mode init is injected via a monkey-patch of `.InitState` that overrides
`phi` and `theta` and recomputes `log_lik`, `log_prior`, `log_post` before the
first MCMC iteration.

**Stdout evidence:** `[truth-init] phi=4.0000 theta=0.9700`

**INIT VERIFIED.** phi_init=4.0000, theta_init=0.9700 confirmed in stdout. The override ran correctly.

> Note: this confirms the override took at **iter 0**. Whether the chain then
> stayed in the phi>1 mode (across warmup and sampling) is the question answered
> by Section 1 below.

---

## 1. Mode-crossing check (critical)

| Metric | Value |
|---|---|
| phi < 1 (reflected mode) | 0 samples |
| phi >= 1 (true mode) | 1,593 samples |
| Zero-crossings of log(phi) | **0** |
| First phi<1 crossing | — |
| First sample: phi / theta_1 | 3.4859 / 0.9898 |
| Last sample: phi / theta_1 | 6.9422 / 0.8871 |
| phi range: min / max | 1.0018 / 77.9822 |

**No mode crossings detected.** The chain, initialised at the truth-mode (phi=4, theta=0.97), never crossed log(phi) = 0. Mode persistence is confirmed for both init directions.

Trace plot: `phi-trace.png`

---

## 2. Settled mode

| Parameter | Median | Truth | Reflected |
|---|---|---|---|
| phi | 3.6445 | 4 | 0.25 |
| theta_1 | 0.9898 | ~0.97 | ~0.03 |

**Chain settled in the **true mode** (phi >= 1, theta_1 near 1).**

---

## 3. Mixing health (ESS)

| Parameter | ESS | n retained | ESS/n |
|---|---|---|---|
| phi | 43.3 | 1593 | 0.027 |
| pi0 | 8.3 | 1593 | 0.005 |
| theta_1 | 8.2 | 1593 | 0.005 |
| tree_length | 99.2 | 1593 | 0.062 |

---

## 4. Tail check (tree_length)

| Metric | Chain B | 20k pilot | Truth |
|---|---|---|---|
| median | 23.280 | ~29 | 12.7 |
| IQR | [14.387, 41.060] | [20.1, 42.3] | — |
| max | 440.908 | ~70 | — |

---

## 5. Comparison to sibling chains

Both chain A (default init, phi<1 mode) and chain B (truth init, phi>1 mode) showed zero mode crossings over 100k iterations. Combined with the sibling chain's result, this confirms bilateral mode persistence: neither mode is escapable on this chain length. Post-hoc relabelling (Option 1) is sufficient.

---

## Conclusion

Truth-init chain stayed in the phi>1 mode for the full 100k iter run (Option 1 sufficient).

