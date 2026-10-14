# MkPrimeEcology

The ecology-biased equilibrium (EBE) model, rebuilt as an MkPrime engine
extension (#441). It replaces the archived `worktree-ecology-aware` branch
(tag `ecology-aware-2026-05-29`, tip `fa20e85`), which extended MkPrime's
engine directly. MkPrime itself carries no ecology code: this package
claims MkPrime's neomorphic characters (extension ABI 3) and scores them
with an edge-aware kernel on MkPrime's `Pruner` (header API 2).

```r
ext <- EcologyExtension(data, ecology, neomorphic = neo)
res <- MkPrime::RunMkPrime(data, rootedTree, neomorphic = neo,
                           model = MkPrime::MkPrimeModel(extensions = list(ext)))
samples <- CanonicaliseEcology(res$samples, ext)
```

## Model

As `dev/ecology/ebe-spec.md` on the archived branch, sections 1-3 and 7:

- an equilibrium tilt `logit(pi1) += s_z log(phi)` per (character, non-reference
  ecology), with the decay rate fixed at 2, so `z = 0` is MkN exactly;
- edge ecology weights: the mean of the JC(K) marginals of tip ecology at the
  edge's two ends (raw edge lengths, uniform root);
- root frequencies `(pi0_base, pi1_base)`; rooted trees, with the root sampled
  by MkPrime's topology moves;
- priors `phi ~ LogNormal(0, sigmaPhi)`, `pi0 ~ Beta(360, 120)`,
  `theta ~ Beta(2, 2)`, `P(z = 0) = pi0`, `P(z = 1) = (1 - pi0) theta`.

## Differences from the archived branch

- Neomorphic edge lengths use MkPrime's normalised `neo_scale` (the RB-style
  rate_neo split now on `main`), so `z = 0` reproduces `main`'s neomorphic
  likelihood.
- Under variable coding a character is conditioned on varying among its
  observed tips, as on `main` (#213), rather than among all tips.
- `z` is updated one character at a time by an exact Gibbs draw over its
  `3^(K-1)` rows (weight `zWeight`), not by a full sweep every
  `gibbsZEvery` generations.
- Not ported: per-chain RNG and OpenMP (T-017), the ecology partial-CL caches,
  the fixed-root option, informative coding, `phi`-tempo variants. Each
  likelihood evaluation recomputes the ecology reconstruction and prunes
  all neomorphic characters: correct, not yet fast.

## Tests

`tests/testthat/test-ebe.R` follows the spec's contracts: T1 (all `z = 0`
equals MkPrime's MkN, through MkPrime's own data and state pointers), T2
(an independent R oracle, with and without ACRV and ascertainment), T3
(non-stationarity under re-rooting), T4 (pattern probabilities sum to one),
the exact Gibbs identity, and a short `RunMkPrime()` run.
.github/workflows/ecology-extension.yml builds MkPrime and this package and
runs them.
