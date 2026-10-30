# Does the empirical k′ prior beat a blind geometric one? (#179)

`empirical_geometric` writes k′ = N_obs + N_unobs, with N_obs drawn from the
corpus pmf of observed state counts (`empiricalNObs`) and N_unobs ~ Geo(p).
`geometric` writes k′ = 2 + Geo(p). Both learn p under the same Beta(1, 1).
The study asks whether the corpus shape buys more accurate k′ on data whose
true k′ follow it, and what it costs when they do not.

## Design

| | |
|---|---|
| Truth | `empirical`: k′ from `empirical_geometric`; `geometric`: k′ from `geometric`. Both truncated to [2, 30], as inference (`kprimeTruncK = 30`). |
| Data | 30 tips, 100 transformational characters, JC(k′), no rate variation; tree from the inference prior (Gamma(2) length, mean 3; Dirichlet split). An invariant (k′, y) pair is redrawn whole (joint ascertainment). |
| Fits | Both priors on every dataset; scores are paired by dataset. |
| k′ scores | Mean per-character log score log Pr(k′ = kTrue) (primary, proper); 90% HPD coverage; bias and RMSE of the posterior mean. |
| Tree scores | Mean normalized clustering information distance and RF from posterior trees to the true tree; both split-based, so blind to rooting. |

`exact.R` scores the exact k′ posterior on the true tree (per-k′ likelihoods
enumerated, p integrated on a grid): no sampler, no topology uncertainty,
seconds per dataset. `mcmc.R` fits with free topology and also records the
exact score, so sampler error is visible per dataset.

The `empirical` configs with p = 0.6 and expected tree length 3 give recorded
kObs that match the corpus (two thirds binary); that is the headline.

## Expected direction

Under the `empirical` truth, `empirical_geometric` should score higher and
carry less bias. The effect is small: the head of the corpus pmf is itself
close to geometric (consecutive ratios 0.32, 0.31, 0.39, 0.43), and a free p
lets `geometric` mimic most of it. Under a mostly-binary `geometric` truth
(p = 0.9) the corpus floor should make `empirical_geometric` overshoot k′ and
lose by more than it gains on its own truth. Topology should not differ.

## Running

From the repository root:

```bash
Rscript dev/red-team/heavy-tests/empirical-benefit/test-forward.R
Rscript dev/red-team/heavy-tests/empirical-benefit/exact.R 50
Rscript dev/red-team/heavy-tests/empirical-benefit/mcmc.R 1 --pilot
```

The full MCMC study is `submit-empirical-benefit.sh` (Hamilton, 200 tasks);
`aggregate.R` summarizes either arm.
