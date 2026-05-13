# Empirical-geometric k' prior: over-recovery diagnosis (2026-05-13)

## Status

Diagnosed; **Fix B implemented and committed** with `Beta(15, 1)` as the
default empirical-Bayes hyperprior on `p` for `kPrimePrior =
"empirical_geometric"`. Smoke test on the original 8-tip, 65-char,
truth-Σu = 50 dataset now lands at Σu posterior ≈ 67 (134 % of truth, in
the target [25, 75]); pre-fix was 172.6 (345 %).

A second smoke (truth Σu = 0, 7-tip / 60 binary chars) revealed a
**structural limit** of Fix B: at the prior asymptote `p → 1` the
per-character k' = kObs → kObs + 1 prior log-ratio approaches
`log(P_emp(kObs + 1) / P_emp(kObs))` ≈ −1.14 for kObs = 2, while the
relabel correction at the same step contributes +1.10. The two nearly
cancel, so under truth = 0 the posterior settles at Σu ≈ 28–32 (≈ 0.5
unseen states per character) regardless of how tight the Beta on `p`
is — verified by sweeping `(α, 1)` up to (200, 1). Pre-fix on the same
dataset gives Σu ≈ 42 (Beta(1, 1)), so Fix B does provide a partial
counteraction (~25 % reduction), but cannot drive Σu close to zero
under truth = 0 alone.

Driving Σu below ~5 under truth = 0 would require either Fix A on top
of Fix B (only helps for kObs > 2 characters, since Z_kObs(p) = 1 when
kObs = 2 — the support is unrestricted), or reconsideration of the
relabel correction itself.

## Symptom

[data-raw/smoke_empirical_geometric.R](../../data-raw/smoke_empirical_geometric.R)
on simulated JC data (8 tips, 65 variable chars, k_true ∈ {2,3,4,5},
true Σu = 50) produced:

| arm | Σ posterior u | % of truth | Mean CID |
|---|---|---|---|
| `empirical_geometric` | 172.6 | **345 %** | 0.84 |
| `geometric` | 5.5 | 11 % | 0.81 |

`block_kPrime` acceptance 0.000 in both arms (symptom, not bug — see
below); both arms hit `maxWarmup` without stabilising.

Logs preserved:
- [data-raw/smoke_eg_2026-05-13.log](../../data-raw/smoke_eg_2026-05-13.log)
- [data-raw/smoke_eg_fixedtree_2026-05-13.log](../../data-raw/smoke_eg_fixedtree_2026-05-13.log)
  — fixed-tree replication, Σu = 154.9 (310 %), p posterior mean
  ≈ 0.32. Topology is innocent.

## Root cause

The chain is **correctly** sampling the posterior under the current
model. Verified line-by-line:

- R and C++ EG prior implementations match
  ([R/MkPrimeModel.R:348–377](../../R/MkPrimeModel.R#L348),
   [src/mcmc.cpp:269–315](../../src/mcmc.cpp#L269),
   [src/mcmc.cpp:3403–3429](../../src/mcmc.cpp#L3403)).
- `mh_logit_p` (case 30) and `joint_p_kprime` (case 31) MH ratios
  including Jacobians are correct
  ([src/mcmc.cpp:4274–4370](../../src/mcmc.cpp#L4274)).
- `K_MAX_CAND = 50`, `LOG_CUTOFF = −25` candidate-set asymmetry
  between `logZ(p)` and `logZ(p')` is bounded ≤ e⁻²⁵ per character —
  negligible.
- Analytical posterior at p ≈ 0.32 gives E[Σu] ≈ 173, matching the
  chain to 3 sig figs.

The **prior model** is the problem. The EG construction
`π(k'|p) = Σ_{j=2}^{k'} P_emp(j) · p·(1−p)^{k'−j}` is a marginal on
**k′**, not on **u given kObs**. Its peak at u = 0 caps at
`P_emp(kObs)`:

| kObs | max π(k' = kObs \| p → 1) |
|---|---|
| 2 | 0.671 |
| 3 | 0.215 |
| 4 | 0.067 |

Meanwhile the relabel correction in the likelihood
`lgamma(k'+1) − lgamma(k'−kObs+1)`
([src/corrections.cpp:35](../../src/corrections.cpp#L35)) pulls k′ upward
by a factor ~10 between k' = kObs and k' = kObs + 2 (kObs = 3). The
EG prior cannot counter this pull regardless of p, so the chain
lowers p to widen the geometric tail, accommodating the relabel-pull
and giving Σu ≈ 173.

Geometric beats this trivially because `(1−p)^u` at high p
concentrates per-character at u = 0 regardless of kObs.

## What the previous (wrong) diagnosis suggested

A subagent first proposed truncating the convolution at `j ≥ kObs_i`
per character — this is **circular** (uses the focal character's
own data to build its own prior; the prior should act *before*
observation). See memory file
`feedback_prior_modelling.md` for the principle.

The unconditional formula is the correct empirical-Bayes marginal on
k′. Conditioning on kObs would be a model re-specification, not a
bug fix.

## Two candidate fixes

### Fix A — per-character prior renormalisation by `Z_{kObs}(p)`

Replace `π(k'|p)` with `π(k'|p, kObs_i) = π(k'|p) · 1[k'≥kObs_i] /
Z_{kObs_i}(p)`. Mathematically not a no-op (changes the marginal
posterior on p via factor `∏_i 1/Z_{kObs_i}(p)`, upweighting high p).
But it uses kObs_i in the renormaliser — defensible (the likelihood
already enforces the same support) but in tension with the
"prior-must-not-depend-on-focal-data" principle.

### Fix B — empirical-Bayes informative Beta(α, β) prior on `p` *(chosen first)*

Keep `π(k'|p)` unconditional. Replace `Beta(1, 1)` on `p` with
`Beta(α, β)` calibrated from `empiricalNObs` so that the *marginal-
on-N_obs* under the convolution roughly matches the corpus pmf at
its mean. Pushes p toward a regime where the geometric tail is
short — bounding the relabel correction's pull.

No new per-character normalisation; no use of focal kObs in the
prior; pure prior on p tuned from corpus (independent population).
Cleaner Bayes.

### Implementation pointers for Fix B

- Locate `kprimeHyperA` / `kprimeHyperB` defaults
  (probably in [R/MkPrimeModel.R](../../R/MkPrimeModel.R) or
  [R/MkPrimePriors.R](../../R/MkPrimePriors.R) — grep first).
- Decide calibration target: simplest is to match
  `E[N_obs | p_obs = ?]` to the corpus mean. Solve for the (α, β)
  whose marginal on k′ under the convolution gives a corpus-like
  marginal on N_obs. Could do this numerically with `optim`.
- The C++ side uses these hyperparameters in `cpp_log_prior` for the
  `p` term — verify the same Beta is applied there. Grep for
  `kprimeHyper`.

## Test gap

[tests/testthat/test-empirical-geometric-prior.R:213–286](../../tests/testthat/test-empirical-geometric-prior.R#L213)
asserts only `expect_gt(sum(uMedEmp), sum(uMedGeo))` — one-sided.
A 345 % overshoot passes. Needs:

- **Two-sided** assertion: `|Σu_post − Σu_true| < tol · Σu_true`,
  e.g. tol = 0.5, on the smoke seed/dataset.
- **Σu_true = 0 case**: synthetic dataset where every char has
  k_true = kObs (e.g. 7 tips, 60 chars all binary, JC, short branches)
  — assert `Σu_post < 5`. Catches the relabel-overruns-prior pathology.
- **R ↔ C++ parity**: pin `p` (via `moveWeights` zeroing all p
  moves, or a freeze flag) and verify `π(k' = kObs)` matches
  `P_emp(kObs)` to 1e-10 on both sides.

## `block_kPrime` zero acceptance

Confirmed **symptom only**. The move proposes a uniform integer
shift in k′; under EG the prior is heterogeneous in k′ (depends on
k', not on u), so a uniform shift always lands in a wildly different
prior region. Under geometric this collapses to `(1−p)^{n·δ}` and
works. The move design predates EG. Will recover once Fix B brings
the EG posterior into a sensible regime; otherwise leave alone for
now ([src/mcmc.cpp:3737–3823](../../src/mcmc.cpp#L3737)).

## Pickup checklist for remote machine

1. Pull this branch.
2. Read this plan + the two smoke logs in `data-raw/`.
3. Implement Fix B:
   - locate `kprimeHyper{A,B}` defaults;
   - calibrate (α, β) from `empiricalNObs` (numeric optimisation);
   - update default in R model constructor and verify C++ honours it;
   - rerun `data-raw/smoke_empirical_geometric.R` — target Σu within
     [25, 75] of the true 50.
4. Add the three tests in
   [tests/testthat/test-empirical-geometric-prior.R](../../tests/testthat/test-empirical-geometric-prior.R)
   listed under "Test gap".
5. If Fix B doesn't bring Σu within tolerance, fall back to Fix A or
   reconsider the relabel correction itself.
