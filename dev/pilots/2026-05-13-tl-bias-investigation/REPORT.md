# tree_length bias investigation (read-only)

**Date:** 2026-05-13 (analysis), persisted 2026-05-14
**Source data:** `dev/pilots/2026-05-13-step3b-mode-persistence-A/result.rds` and
`chain.log` (100k chain A, default init, truth mode: phi=3.6, theta=0.94,
pi0=0.235).

**Verdict:** No bug located in C++ likelihood or R prior. Most likely
explanation is a unit/parameterisation mismatch between simulator and model.
A discriminating re-sim is needed before "no bug" is airtight.

## Summary table

| Hypothesis | Verdict | Confidence |
|---|---|---|
| 1. gamma-normalisation direction error | NOT a bug | High |
| 2. Reference-ecology dummy coding miss | NOT a bug | High |
| 3. Mixing / burn-in artefact | NOT the cause | High |
| 4. Prior on tree_length pushing up | NOT the cause | High |
| 5. Rate composition / unit mismatch | Likely root cause | Medium |
| (new) pi0 bias | Concern, possibly contributing | Medium |

## Hypothesis 1 — gamma-normalisation direction

Checked `src/mcmc_ecology.cpp`:
- `gamma_e_compute` (line 249): `pi0 + (1-pi0)*(theta_e*phi + (1-theta_e)/phi)`.
  Equals `E[mu | z]` over slab/spike. Correct prior-expected rate factor.
- `trans_rate_factor` (line 260, 266): `return mu / gamma_e;` — **divided** by
  gamma_e. Correct direction.
- `mkn_rates_for_state` (lines 490-491): `r01 = rate01_base * mu01 / gamma_e;
  r10 = rate10_base * mu10 / gamma_e;` — also divided. Correct.
- Reference ecology short-circuits with factor=1 in both paths (lines 263,
  366-367, 476-480). Correct.

**Not a bug.**

## Hypothesis 2 — reference-ecology dummy coding miss

Checked JC path (lines 392-398) and MkN path: zMat indexed at `kEco-1` columns
with explicit skip past `refEcology`; refEcology branches short-circuit before
consulting zMat. Code is unambiguous.

**Not a bug.**

## Hypothesis 3 — mixing / burn-in artefact

Quartile breakdown of post-burn-in tree_length (n=1593 retained):

| Quartile | tl median | tl mean |
|---|---|---|
| Q1 | 18.86 | 29.81 |
| Q2 | 22.30 | 28.36 |
| Q3 | 18.51 | 23.71 |
| Q4 | 21.10 | 26.81 |

Stationary, no monotone drift. ESS low (~71-99) but chain is sampling a
genuine posterior with median ~20.

**Not a convergence issue.**

## Hypothesis 4 — prior on tree_length

`R/MkPrimeModel.R:525-528`: `dgamma(tree_length, shape=2, rate=2/expSteps =
2/10 = 0.2)`. Gamma(2, 0.2) has mean 10, mode 5, 95% quantile 23.7. Truth
(12.7) and posterior median (20.4) both inside the prior's bulk. Prior is
diffuse — not pulling tree_length up.

Relative branch lengths use Dirichlet(1,...,1), uniform on simplex.
Marginal `E[br_i] = 1/30`.

**Not the cause** (but Dirichlet matters for Hyp 5).

## Hypothesis 5 — rate composition / unit mismatch (likely root cause)

This is the most likely explanation. It is **not a bug** — it is a mismatch
between simulator parameterisation and model parameterisation.

**Simulator** (`inst/simulations/ecology/sim3-simulate.R`, lines 154/178/239):
for trans chars, `r = baseRate * mult[z+1] / gE`, with `gE = 1` when
`normalize=FALSE` (the default; chain A's run does not pass `normalize=TRUE`).
With `baseRate=0.5`, trans rate per char per edge = 0.5 regardless of ecology.
Per-edge expected substitutions for trans char on a 0.5-time edge = 0.25.

**Model** (`src/mcmc_likelihood.cpp:210-213`): JC kernel uses `t = elPtr[e]`
directly. No separate `baseRate` for trans chars — rate fully absorbed into
edge_length. For neo chars there is a separate `rate_neo` scalar
(`mcmc_ecology.cpp:965`).

The user's framing "truth tree_length = 12.7" implicitly equates simulator
time units with model substitution units. They are not equal:

- Trans on ref edges: 0.5 * 6.7 = **3.35** model units
- Trans on non-ref edges: 0.5 * 6.0 = **3.0** model units (z=0 truth, mult=1, gE_sim=1)
- Trans-equivalent total: **6.35** model units
- Neo on ref edges: 1 * 6.7 = **6.7** model units
- Neo on non-ref edges (z=1): 4 * 6.0 = **24** model units
- Neo-equivalent total: **30.7** model units

The simulator's rate-time product is character-class-dependent; the model has
a single `tree_length` (with `rate_neo` / `rate_loss` adjusters) and finds a
compromise. Posterior median 20.4 sits between the trans-view (~6) and
neo-view (~30).

**Per-edge diagnostic is inconclusive.** Computing `br_i * tree_length` per
sample then taking column medians: positions 1, 16 (roots) at ~0.55
(truth 0.15, 3.7x); positions 2, 9, 17, 24 (clade stems) at ~0.50
(truth 0.10, 5x); positions 3-8, 10-15, 18-23, 25-30 (within-clade) at ~0.50
(truth 0.50, 1.0x). This looks like a stems-inflate-but-internal-recovered
story — but every column is suspiciously near `tree_length / 30 ≈ 0.67`, the
Dirichlet(1,...,1) marginal mean × posterior tree_length. With topology
mixing (most common topo_hash holds 0.5% of samples), `br_i` does not
consistently track one biological edge. The apparent pattern is most
plausibly Dirichlet-prior reversion under topology averaging, not a real
signal.

## New concern — pi0 bias

Truth empirical pi0 = 180/240 (all trans chars z=0) = **0.75**. Posterior
pi0 median = **0.235**. That is 3x downward bias.

Truth-aware gamma_e using pi0_truth = 0.75:
`0.75 + 0.25*(0.94*3.6 + 0.06/3.6) ≈ 0.75 + 0.85 = 1.60`. Posterior gamma_e
≈ 2.65 (1.65x inflated). The chain compensates with larger non-ref edge
lengths.

Probably not a bug — likely the slab prior is too weak to constrain pi0 when
most chars are weakly informative. Worth tracking through multi-rep.

## Posterior summary (chain A, post 25%-burnin, n=1593)

```
phi   median: 3.600 IQR: [2.397, 5.931]    truth 4
pi0   median: 0.235 IQR: [0.171, 0.344]    truth empirical 0.75
theta median: 0.941 IQR: [0.919, 0.960]    truth 1.0
tl    median: 20.43 IQR: [13.85, 33.95]    truth-time 12.7, model-units ambiguous
gE    median: 2.65  IQR: [1.85, 4.50]      truth-pi0-aware 1.60
rate_loss   median: 0.434 (truth=1)
rate_neo    median: 2.10
cor(log tl, log gE): +0.344
```

## Recommendations

1. **Do not patch.** No bug located.

2. **Run the discriminating test** before declaring "no bug" airtight.
   Re-simulate with simulator parameterisation aligned to the model:
   - `.SimulateMkPrimeEcology(..., baseRate = 1.0, normalize = TRUE,
     pi0 = 0.75, theta = 1.0, refEcology = 0)`
   - Both `normalize` and `baseRate` are already exposed (lines 151-156).
   - Then re-run a 100k chain. If posterior tl median collapses to ~12.7 and
     pi0 to ~0.75, the unit-mismatch hypothesis is confirmed. If bias
     persists, Hypothesis 5 needs deeper investigation.

3. **Flag pi0 bias for multi-rep monitoring.** Even after parameterisation
   alignment, if pi0 systematically under-shoots, that is a separate
   inference concern (not the gamma fix). Add to the multi-rep summary
   script.

4. **Document in vignette** that `tree_length` in v2 is in model
   substitution-per-site units, not simulator time units. Conversion depends
   on `baseRate`, `rate_neo`, and per-character ecology assignments — no
   single scalar. For Sim 3 reporting either re-simulate with aligned
   parameterisation, or report tl in model units without comparing to a
   "truth" scalar.

## Key files referenced

- `src/mcmc_ecology.cpp` (lines 249, 260-267, 360-378, 470-492)
- `src/mcmc_likelihood.cpp` (lines 207-225)
- `R/MkPrimeModel.R` (lines 520-534)
- `inst/simulations/ecology/sim3-simulate.R` (lines 150-185, 210-249)
- `dev/pilots/2026-05-13-step3b-mode-persistence-A/chain.log`
- `dev/pilots/2026-05-13-step3b-mode-persistence-A/run.R`
