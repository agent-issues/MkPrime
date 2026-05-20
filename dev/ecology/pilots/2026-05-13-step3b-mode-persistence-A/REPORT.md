# Step 3b Mode-Persistence Report — 100k Sim 3 chain

**Date:** 2026-05-14
**Plan:** `dev/plans/2026-05-13-1100-v2-asymmetric-slab-prior-and-gamma-ridge.md`, Open question: label-switching

---

## Configuration

| Parameter | Value |
|---|---|
| Script | `dev/pilots/2026-05-13-step3b-mode-persistence-A/run.R` |
| Seed | 20260513 (one more than the 20k pilot — independent chain) |
| kEco | 2 (ecology 0 = C/D clades; ecology 1 = A/B clades) |
| nTheta | 1 (one non-reference ecology) |
| nEco | 60 neomorphic characters |
| nBase | 180 transformational characters |
| phi_truth | 4 |
| stem / root | 0.10 / 0.15 (Goldilocks) |
| nIter | 100,000 |
| thin | 40 |
| nChains / nRuns | 1 / 1 |
| minWarmup / maxWarmup | 2,000 / 5,000 |
| Initial state | default (.InitState) |
| model | ecologyAware=TRUE, magnitudeMode="global", kPrimePrior="geometric" |
| Log file | `dev/pilots/2026-05-13-step3b-mode-persistence-A/chain.log` |
| RDS | `dev/pilots/2026-05-13-step3b-mode-persistence-A/result.rds` |

**Post-burnin samples:** 2125 raw rows → 1593 retained after 25% discard.

---

## 1. Mode-crossing check (critical)

| Metric | Value |
|---|---|
| phi < 1 (reflected mode) | 0 samples |
| phi >= 1 (true mode) | 1,593 samples |
| Zero-crossings of log(phi) | **0** |
| First sample: phi / theta_1 | 3.1915 / 0.9187 |
| Last sample: phi / theta_1 | 2.3947 / 0.9908 |
| phi range: min / max | 1.0750 / 44.1676 |

**No mode crossings.** The chain never crossed log(phi) = 0.

Trace plot: `phi-trace.png`

---

## 2. Settled mode

| Parameter | Median | Truth | Reflected |
|---|---|---|---|
| phi | 3.5999 | 4 | 0.25 |
| theta_1 | 0.9405 | ~0.97 | ~0.03 |

**Chain settled in the **true mode** (phi >= 1, theta_1 near 1).**

---

## 3. Mixing health (ESS)

| Parameter | ESS | n retained | ESS/n | 20k pilot ESS/n |
|---|---|---|---|---|
| phi | 47.1 | 1593 | 0.030 | 8.9/131 = 0.068 |
| pi0 | 19.1 | 1593 | 0.012 | 3.7/131 = 0.028 |
| theta_1 | 12.0 | 1593 | 0.008 | 10.2/131 = 0.078 |
| tree_length | 71.4 | 1593 | 0.045 | 8.1/131 = 0.062 |

---

## 4. Tail check (tree_length)

| Metric | 100k chain | 20k pilot | Truth |
|---|---|---|---|
| median | 20.433 | ~29 | 12.7 |
| IQR | [13.849, 33.954] | [20.1, 42.3] | — |
| max | 182.818 | ~70 | — |

---

## 5. Comparison to 20k pilot

100k chain: **zero mode crossings**. Single-chain mode persistence is confirmed. The 20k pilot result (phi inverted to reflected mode, theta_1 near zero) is replicated: the chain consistently finds one of the two symmetric modes and stays there. Post-hoc relabelling (Option 1) is sufficient.

---

## Conclusion

A single 100k chain shows zero log(phi) zero-crossings. Mode persistence is confirmed. Post-hoc relabelling (Option 1 from the plan) is sufficient to handle the phi <-> 1/phi symmetry. No prior symmetry-breaking is required for single-chain inference; multi-rep analyses should apply relabelling before cross-replicate summary.

