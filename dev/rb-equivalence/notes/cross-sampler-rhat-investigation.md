# Cross-sampler R-hat investigation: MkPrime vs RevBayes

**Status:** RESEARCH ONLY — no code changes. Read-only audit of model parameterisation differences likely to drive the residual cross-sampler R-hat observed in the smoke run (pid 635 × by_nt_9v).

## TL;DR

Two specific source-level asymmetries between MkPrime's and RevBayes's implementations of the "shared" Mk + ACRV + asymmetric-Mk2 model parsimoniously explain the entire signature of the observed cross-sampler drift (`tree_length` rhat 1.037, `rate_log_sd` rhat 1.027, `rate_neo` rhat 1.030; `rate_loss` and tree-CID agree). **(1) Partition-rate parameterisation:** RB applies `partition_rate := [rate_neo/(1+rate_neo), 1/(1+rate_neo)] * sum(nChar)/nChar` to both partitions, so the nChar-weighted mean rate is 1 by construction; MkPrime applies `rate_neo` as a one-sided neomorphic branch-length multiplier (`neoEl = edgeLen * rate_neo`) and leaves the transformational partition at rate 1, so the weighted mean rate is `(n_neo · rate_neo + n_trans) / (n_neo + n_trans)`. That makes `rate_neo` a fundamentally different parameter on the two sides and breaks the implicit `tree_length` × `rate_neo` identification, explaining the exact pair of params (`tree_length`, `rate_neo`) that drift. **(2) ACRV discretisation normalisation:** both samplers use `qnorm((i+0.5)/nCat, μ=−σ²/2, σ=rate_log_sd)` for the bin medians, but MkPrime then renormalises so the arithmetic mean of the discrete rates is exactly 1, while RB's `fnDiscretizeDistribution` does NOT renormalise. At σ=1 this is a ~13 % multiplicative bias on the per-site rate, absorbed by `tree_length` and `rate_log_sd`. These two asymmetries are sufficient to produce a 1.03–1.04 rhat gap on exactly the observed parameters; `rate_loss` and the tree CID are correctly insensitive (rate_loss enters the Q-matrix only, not the time axis).

> **Correction (2026-09-24, agent-issues/MkPrime#212, #213):** two of this note's AGREED verdicts were wrong. (a) MkPrime's neomorphic Q was not normalised: its π-weighted mean rate was 4r/(1+r)², whereas RevBayes' `fnFreeK` rescales to 1, so `rate_loss` *did* scale the neomorphic time axis in MkPrime. MkPrime now normalises to mean rate 1. (b) §4: MkPrime pinned every tip to the constant state, whereas RevBayes marginalises each character's ?/- tips (one correction per missing-data mask). MkPrime now does the same. Cross-sampler comparisons recorded before the fix compare different models on `rate_loss`, `rate_neo` and `tree_length`, and on any matrix with missing data.

---

## Context

Smoke run: pid 635, by_nt_9v, 2 runs × 4 PT chains per sampler. Within-sampler maxRhat = 1.001. Cross-sampler rhat (treating both samplers' cold chains as four sources):

| param | cross-sampler rhat | pooled ESS |
|---|---:|---:|
| tree_length   | 1.037 | 445  |
| rate_log_sd   | 1.027 | 794  |
| rate_loss     | **1.004** | 3131 |
| rate_neo      | 1.030 | 1099 |
| cid_to_median | **1.002** | 1983 |

So `rate_loss` (Q-matrix asymmetry, no time-axis role) and the tree CID (topology) agree; the three drifting parameters are all on the time-axis / rate scale.

---

## Candidate-by-candidate audit

### 1. Branch-length parameterisation — AGREED

**RB** (`templates/by_nt_9v.template.Rev:29–34`):
```
rel_br_lengths ~ dnDirichlet( rep(1.0, nEdge) )
br_lengths := rel_br_lengths * tree_length
phylogeny := treeAssembly(topology, br_lengths)
```

**MkPrime** carries the same `tree_length × rel_br_lengths` parameterisation internally: `R/RunMkPrime.R:1098`, `:3137-3144`, `:3340-3341`, `:3884`, `:4085` all reconstruct `edge.length = tree_length * rel_br_lengths`. The Dirichlet(1,…,1) prior on the simplex is included as `lfactorial(nEdges - 1)` in `R/MkPrimeModel.R:434–436`, equivalent to RB's `dnDirichlet(rep(1, nEdge))` log-density (an additive constant on the prior that does not bias posteriors).

Verdict AGREED. No predicted bias.

### 2. ACRV discretisation — DIFFER (normalisation)

**RB** (`src/core/functions/distribution/DiscretizeDistributionFunction.cpp:99–108`):
```cpp
for (int i=0; i<n_cats; ++i) {
    double p = (i+0.5)/n_cats;
    (*value)[i] = dist->quantile(p);
}
```
That is the median of each equal-mass bin under `dnLognormal(0, σ)`. No renormalisation.

**MkPrime** (R fallback `R/acrv.R:33–39`):
```r
midpoints <- (seq_len(nCat) - 0.5) / nCat
logRates <- qnorm(midpoints, mean = mu, sd = rateLogSd)
rates <- exp(logRates)
rates / mean(rates)             # <-- normalisation absent in RB
```
C++ hot path (`src/mcmc_likelihood.cpp:75–87`, with `mu = -rateLogSd^2/2`):
```cpp
for (int i = 0; i < nCat; ++i) {
  rates[i] = std::exp(mu + rateLogSd * acrvZ[i]);
  total   += rates[i];
}
for (int i = 0; i < nCat; ++i) rates[i] *= nCat / total;  // <-- ditto
```

Median-quantile midpoints of `LogNormal(−σ²/2, σ)` do not have arithmetic mean 1 in finite-bin discretisations; for σ=1 the empirical mean of six bins is ~0.885 (verified `Rscript`). RB therefore applies an effective per-site rate that is ~10–15 % lower than MkPrime at the same `rate_log_sd`.

Predicted bias direction: at matched σ, RB's likelihood "thinks" the average rate is ~0.88, so the posterior compensates by inflating `tree_length` and (to a lesser extent) `rate_log_sd` (since increasing σ further shifts the median-quantile mean even further from 1, providing a non-trivial Jacobian-like coupling). Predicted: `tree_length_RB > tree_length_MkP`, with the gap scaling with `rate_log_sd`.

Verdict DIFFER. Effect on tree_length and rate_log_sd.

### 3. Partition rate scaling — DIFFER (this is the smoking gun for `rate_neo`)

**RB** (`templates/by_nt_9v.template.Rev:47`):
```
partition_rate := [ rate_neo / (1 + rate_neo), (1 / (1 + rate_neo)) ] / nChar * sum(nChar)
```
Both partitions get a scalar, with normaliser `sum(nChar)/nChar[i]`, so the nChar-weighted mean rate is `(rate_neo + 1)/(1+rate_neo) = 1`. `tree_length` retains its interpretation as "expected number of state changes per site averaged over partitions weighted by partition size".

**MkPrime** (`src/mcmc_likelihood.cpp:2371–2372` and `R/likelihood.R:102`):
```cpp
NumericVector neoEl(edgeLen.size());
for (int i = 0; i < edgeLen.size(); ++i) neoEl[i] = edgeLen[i] * rateNeo;
```
The transformational partition uses `edgeLen` unscaled (no symmetric counter-scaling). The weighted mean rate is therefore `(n_neo · rate_neo + n_trans) / (n_neo + n_trans)`, which is 1 only at `rate_neo = 1`. So on the MkPrime side, the joint `(tree_length, rate_neo)` parameterisation is reparameterised to a different basis than on the RB side. Equivalently, a fixed posterior on "expected total changes" maps to:
- RB: identifiable `tree_length` (it IS that quantity, weighted mean over partitions = 1)
- MkPrime: `tree_length × (n_neo · rate_neo + n_trans)/(n_neo + n_trans)` (an over-parameterised pair: TL and rate_neo trade off linearly)

Predicted bias direction: posterior marginals on `rate_neo` differ in mean and variance between samplers; `tree_length` posteriors differ in proportion to the difference in the rate_neo posterior, in a way that is consistent across replicates (hence systematic, not within-sampler noise). This explains why exactly `tree_length` and `rate_neo` drift in lockstep, while `rate_loss` (which enters the Q-matrix asymmetry, not the time axis) is invariant across samplers.

Verdict DIFFER. Predicted to dominate the `rate_neo` drift; mediates a fraction of the `tree_length` drift.

### 4. Variable-coding ascertainment correction — DIFFERED (missing data; fixed in #213)

**Mkn (neomorphic):** MkPrime's `constant_site_prob_mkn` (`src/ascertainment.cpp:428–514`) constructs 2 pseudo-characters with `cl[s][s]=1` (i.e. all-0 and all-1 patterns), prunes under the same Q + rate_loss + ACRV, sums by `root_freqs`, and the caller subtracts `nChar * log(1 - p_const)` (`src/mcmc_likelihood.cpp:2431`). RB's `coding="variable"` excludes the same two constant patterns under the asymmetric Q + stationary root frequencies. ✓

**JC(9) (transformational):** `constant_site_prob_jc` (`src/ascertainment.cpp:32+`) and `pruning_jc_acrv` symmetrically subtract all 9 monomorphic patterns. ✓

Verdict AGREED. No predicted bias.

### 5. Compound-Dirichlet vs decoupled branch lengths — AGREED

See §1. MkPrime mirrors RB's compound-Dirichlet parameterisation; the Dirichlet log-density constant is included.

Verdict AGREED.

### 6. Tree-topology prior — AGREED (in expectation)

RB: `dnUniformTopology(taxa)`. MkPrime has no explicit topology prior term; with symmetric NNI/SPR/TBR proposals (`R/RunMkPrime.R:3416,3454,3461`) the implied prior is uniform on labelled topologies. The tree-CID rhat of 1.002 corroborates that topology posteriors agree.

Verdict AGREED. Could matter if MkPrime adds Gibbs topology moves with subtly asymmetric Hastings ratios, but the CID rhat says it doesn't in practice.

### 7. Rooted vs unrooted — AGREED for likelihood, NEEDS-EMPIRICAL-CHECK for proposal symmetry

MkPrime stores trees in `TreeTools::Preorder` form (rooted-by-convention edge ordering), but the Felsenstein pruning with JC root_freqs = (1/k,…,1/k) is invariant to root location for JC partitions (`src/acrv.cpp:30+`). For the MkN partition the root frequencies are the stationary distribution (`src/mcmc_likelihood.cpp:89–93`), so the likelihood is also root-invariant for MkN.

RB's `dnUniformTopology(taxa)` is the unrooted-binary topology space.

Verdict AGREED. Tree-CID rhat of 1.002 says this is empirically clean.

### 8. Float / numerical precision — AGREED in expectation

MkPrime uses `double` throughout; expm1 and inv-k tricks are mathematically equivalent to RB's matrix-exponential approach modulo ULP-level rounding. The `expm1` form *improves* small-`λt` accuracy relative to the naive `1 - exp(arg)` formulation. Cannot produce a 1.03 rhat at this sample count.

Verdict AGREED.

---

## Why `rate_loss` does not drift

`rate_loss` enters MkPrime *only* through the Mk2 Q-matrix construction (`src/mcmc_likelihood.cpp:2290–2293`, `pruning_mkn`/`pruning_mkn_acrv` at `src/acrv.cpp:290–293`) and the corresponding root stationary frequencies. It does NOT multiply branch lengths and does NOT interact with the partition_rate formula on either side. Both samplers therefore agree on its likelihood and the posterior is identical up to MC noise — exactly the observed rhat = 1.004.

`rate_neo` *does* multiply branch lengths, and that multiplication differs between the two samplers (§3). The asymmetry between `rate_loss` and `rate_neo` posteriors is the strongest single piece of evidence that the partition-rate formula is the culprit.

---

## Ranked causes (most likely first)

| Rank | Cause | Predicted effect | Magnitude |
|---|---|---|---|
| 1 | §3 partition-rate parameterisation | `rate_neo` drift; coupled `tree_length` drift | dominant (explains 1.030 on `rate_neo`) |
| 2 | §2 ACRV bin normalisation | `tree_length` ↑ in RB, `rate_log_sd` weakly | ~10–15 % multiplicative on TL at σ=1 |
| 3 | §1, §5 BL prior constants | none on posteriors, only on absolute log-priors | nil |
| 4 | §6 topology prior | nil under symmetric proposals | nil |
| 5 | §4 ascertainment | differed on characters with ?/- cells (#213) | unmeasured |
| 6 | §7 rooting / §8 float | nil at this sample count | nil |

Notably, **(1) and (2) operate on different parameters but both impinge on `tree_length`.** Their combination cleanly reproduces the observed pattern:
- `tree_length` (rhat 1.037): both causes contribute (sum of ACRV mean shift + rate_neo reparam absorption).
- `rate_neo` (rhat 1.030): partition-rate reparam (cause 1).
- `rate_log_sd` (rhat 1.027): ACRV mean shift (cause 2); σ posterior shifts so that the discrete-rate mean approximates 1 in RB.
- `rate_loss` (rhat 1.004): clean (no time-axis role).
- tree CID (rhat 1.002): clean (topology is invariant to all of the above).

---

## Empirical checks I did not run (research-only constraint)

1. Re-fit MkPrime with `partition_rate`-style scaling: replace `neoEl = edgeLen * rate_neo` with the RB-equivalent symmetric formula and confirm `rate_neo` rhat drops to ~1.005. The right place to test this is `src/mcmc_likelihood.cpp:2371–2372` plus the analogous transformational path (which currently uses `edgeLen` unscaled).
2. Re-fit MkPrime without the `rates / mean(rates)` normalisation in `R/acrv.R:39` and `src/mcmc_likelihood.cpp:85` and confirm `tree_length` and `rate_log_sd` rhats drop.
3. Alternatively, re-fit RB *with* the median-quantile-then-renormalise scheme (custom deterministic node) and confirm symmetric convergence.

Each fix is a one-line change on one side; the two together should drop all three cross-sampler rhats to within-sampler noise.

---

## Two related but lower-priority discrepancies surfaced during the audit

- **Mk' relabelling correction** (`src/corrections.cpp:35–45`) is added on the MkPrime side for transformational characters but is absent in RB's `fnJC(9)` model. When `kPrime == kObs` (always true under by_nt_9v because `knownStates` pins k=9 = max possible observed states for 9-state chars), the correction reduces to `log(kObs!)` per character — a constant that cancels in MH ratios on a single dataset, so it does not bias posteriors at fixed `kPrime`. *But* if any transformational character has `kPrime > kObs` floating across the chain, that correction term will systematically pull MkPrime's posterior. For by_nt_9v with `kPrime = 9 = kObs`, this is inert; for by_nt_kv it could matter.
- **Per-partition rate normalisation for the transformational partition under MkPrime** is implicitly 1, but RB's formula gives `1/(1+rate_neo) * sum(nChar)/nChar[trans]`, which for the by_nt_9v Casali numbers (mostly trans-heavy) is close to 1 only at `rate_neo ≈ 1`. So even the trans partition sees a different effective rate between samplers when `rate_neo ≠ 1`. This is a subset of cause (3) but worth flagging since it means the bias on `tree_length` is not solely mediated through the neomorphic partition.

---

## Files referenced

- `dev/rb-equivalence/templates/by_nt_9v.template.Rev:1–87`
- `dev/rb-equivalence/run_mkprime.R:113–125`, `:206–211`
- `dev/rb-equivalence/R/utils.R:84–97`
- `R/MkPrimeModel.R:104–229`, `:393–572` (priors)
- `R/acrv.R:15–40` (R-fallback discretisation)
- `R/likelihood.R:43–251` (orchestrator)
- `src/mcmc_likelihood.cpp:73–87` (ACRV bins in C++), `:2325–2433` (per-partition likelihood, including the `neoEl = edgeLen * rateNeo` one-sided rate scalar)
- `src/acrv.cpp:24–139` (JC ACRV pruning), `:263–370` (MkN ACRV pruning)
- `src/ascertainment.cpp:32+` (JC constant-site prob), `:428–514` (MkN constant-site prob)
- `src/corrections.cpp:35–72` (relabelling correction)
- RB upstream: `revbayes/src/core/functions/distribution/DiscretizeDistributionFunction.cpp:99–108`
