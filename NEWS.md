# MkPrime (development version)

## Pooled half-normal hyperprior on per-class `class_rate_log_sd`

Per-class ACRV-shape parameters `σ_c = class_rate_log_sd[c]` now share
information across user classes via a half-normal hyperprior on a
population scale τ. The structure (default for `unlink = "shape"` with
two or more classes) is

```
σ_c | τ  ~ HalfNormal(τ)    # per-class scale
τ        ~ HalfNormal(1)    # population scale
```

implemented in **non-centred** parameterisation `σ_c = τ · z_c` with
`z_c ~ HalfNormal(1)` i.i.d. to avoid the small-τ funnel that the
centred form induces. `σ_c` is what the likelihood consumes; `z_c` and
`τ` are the parameters sampled by the chain. Trace plots and ESS now
include `hyper_tau` and `class<c>_rate_log_sd_z` columns in
`result$samples` (appended after the existing `class<c>_rate_log_sd`
columns; downstream readers index by name and so are unaffected).

**Rationale.** AutoPart Casali pilot runs (95 cells in the first batch)
showed `class*_rate_log_sd` as the dominant ESS bottleneck on multi-class
cells. The smallest class in a cell carries little likelihood signal on
within-class rate dispersion, so its independent `Gamma` prior leaves
σ_c diffusing across orders of magnitude. Pooling σ_c through τ uses the
larger classes to anchor the population scale, dramatically lifting
mixing on the small class.

**Deviation from the literal `σ_c | τ ~ HN(τ)` prior.** The
non-centred parameterisation samples `z_c ~ HN(1)` and reconstructs
`σ_c = τ · z_c`. The resulting marginal on σ_c IS half-normal with scale
τ; the deviation is purely in the sampling representation, not in the
prior support or density.

**Selection / opt-out.** The new prior is the default. The legacy
independent `Gamma(rateLogSdShape, rateLogSdRate)` prior on each σ_c
remains selectable via
`MkPrimeModel(priorOnClassRateLogSd = "gamma_independent")` for
backward-compatibility or prior-sensitivity analyses.

**Move set under `unlink = "shape"`.** The legacy scalar `rate_log_sd`
moves (`rate_log_sd`, `slice_rate_log_sd`, `joint_tl_rls`) are dropped
from the partitioned move list — they would break the lockstep between
the scalar `state->rateLogSd` and `state->classRateLogSd[0]` that the
partial-CL fallback paths rely on. σ_0 is moved via case 31 with
`classIdx = 1` (i.e. `scale_class_rate_log_sd_1`). Under the new prior
a single global `scale_hyper_tau` Bactrian move on τ is also added.

The symmetric concern in the **shape-linked** partitioned regime
(`partition != NULL` but `"shape" %notin% unlink`) — where the scalar
`state->rateLogSd` is the only σ but the partition-aware likelihood
path still reads from `state->classRateLogSd[0]` — predates this change
and is not addressed here. A future tidy-up should either fold the
shape-linked case through the same move filter or update case 2 / case
19-slice to keep both fields in lockstep when `state->usePartitioned`
is true.

**Casali ESS benchmark.**
`dev/red-team/heavy-tests/casali-ess-hyperprior-vs-gamma.R` runs the
seven Casali pilot cells (`auto-part/dev/benchmarks/casali`) whose
`class*_rate_log_sd` were the dominant ESS bottleneck in the first
production batch. Each cell runs under both priors at the production
config (`nGen = 5e6`, `thin = 1000`, `nChains = 4` PT, `nCat = 4`),
with seeds paired across the prior pair so the only difference between
the two runs is the σ_c prior structure (SLURM arrays `17304194` and
`17306284`, Hamilton, 2026-05-28).

Per-class σ_c ESS (γ = `gamma_independent`; H = `hyperprior_pooled`).
Five cells are reported at N = 1 paired seed; the two K = 3 cells were
replicated at N = 3 (`17306284`) after a tail finding on the first
batch — see the **N = 3 replication** subsection below.

| cell | K | σ-min γ | σ-min H | ratio | σ-mean γ | σ-mean H | ratio | TL γ | TL H | τ ESS | wall γ (h) | wall H (h) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Allain2012/T2a         | 11 |  870 |  918 | 1.06 | 1028 | 1368 | **1.33** | 1849 |  997 | 1241 | 5.96 | 5.62 |
| Brochu2010/T1          |  2 | 1431 | 1453 | 1.02 | 1899 | 2111 | 1.11 | 3574 | 3826 | 2910 | 2.11 | 1.97 |
| Allain2012/T1          |  2 | 1773 | 1923 | 1.08 | 2122 | 2435 | 1.15 | 1600 | 1360 | 2850 | 3.54 | 3.36 |
| AllainAquesbi2008/T1   |  2 | 1946 | 2337 | **1.20** | 2946 | 3193 | 1.08 | 4789 | 4789 | 3073 | 4.10 | 4.50 |
| CarranoSampson2008/T2a |  3 | *see below* |  |  |  |  |  |  |  |  |  |  |
| Burns2011/T4           |  3 | *see below* |  |  |  |  |  |  |  |  |  |  |
| Godefroit2008/T1       |  2 | 2803 | 2682 | 0.96 | 3070 | 2839 | 0.92 | 4320 | 3804 | 2597 | 2.21 | 2.14 |

**N = 3 replication on the K = 3 cells.** The first single-seed pair on
Burns2011/T4 showed an alarming σ-min ratio of 0.38 (one class dropped
to 836 ESS under the hyperprior while the other two and τ mixed fine).
Two additional seeds per prior on Burns2011/T4 and the K = 3 control
CarranoSampson2008/T2a (SLURM array `17306284`) resolve the picture:

| cell | rep | γ per class | H per class | γ σ-mean | H σ-mean |
|---|---|---|---|---:|---:|
| Burns2011/T4 (20/19/21)     | rep1 | 2182 / 2912 / 2390 | 1652 / 2159 / **836**  | 2495 | 1549 |
| Burns2011/T4                | rep2 | 1788 / 2128 / 1685 | 2283 / 1949 / 2346 | 1867 | 2193 |
| Burns2011/T4                | rep3 | 1276 / 1886 / 2000 | 1809 / 2720 / 2007 | 1721 | 2179 |
| **Burns2011/T4 — N = 3 mean** |    |    |    | **2028** | **1974** |
| CarranoSampson/T2a (28/22/22) | rep1 | 2312 / 2034 / 2120 | 2135 / 1957 / 2236 | 2155 | 2109 |
| CarranoSampson/T2a            | rep2 | 2530 / 2124 / 2671 | 1931 / 2233 / 2635 | 2442 | 2267 |
| CarranoSampson/T2a            | rep3 | 2075 / 2186 / 2218 | 2380 / 2507 / **1005** | 2159 | 1964 |
| **CarranoSampson/T2a — N = 3 mean** |    |    |    | **2252** | **2113** |

Averaged over three seeds, the pooled prior is **essentially neutral**
on both K = 3 cells (σ-mean ratios 0.97 and 0.94). The class-stuck
behaviour observed on Burns2011/T4 rep1 (class3) and CarranoSampson
rep3 (class3) reflects a **stochastic tail hazard of the non-centred
parameterisation**: one z_c can occasionally remain near its
initialisation when the population scale τ moves quickly past it, while
τ and the other z_c continue to mix. In our 6-seed budget under the
hyperprior on K = 3 cells this happened on 2 of 6 seeds, never affected
more than one class at a time, and never affected τ or tree-length
mixing. It is not data-specific.

**Reading the table.** The pooled prior:
* Lifts mean per-class σ ESS on the highest-K cell (Allain2012/T2a,
  K = 11): +33 %, the regime the change was designed for. The N = 1
  result is consistent with a real gain, but the seed-to-seed noise
  measured on the K = 3 cells (σ-mean varies ± ~25 % across seeds
  under either prior) means a single +33 % point should be read as
  suggestive, not confirmed.
* Gives a +20 % σ-min lift on AllainAquesbi2008/T1 (N = 1; same
  caveat).
* Is essentially neutral on K ∈ {2, 3} once seed variance is
  controlled for (N = 3 means on the two replicated K = 3 cells;
  the four other K ∈ {2} cells run at N = 1 with per-cell ratios in
  [0.92, 1.20]).
* Carries a per-seed stuck-z_c tail risk (~1/3 frequency in our
  budget) where one of the K classes drops to ~50 % of its
  γ-equivalent ESS. Other classes, τ, tree-length and log-posterior
  mixing are unaffected.
* Tree-length ESS is mostly flat-to-slightly-worse under the
  hyperprior on the N = 1 cells; this is plausibly the dropped
  `joint_tl_rls` move (acknowledged in the prior-section above).
  Wallclock is comparable (± 10 % cell-by-cell). `hyper_tau` itself
  mixes well (ESS 1241–3073 across all cells).

Net read: the pooled prior is a reasonable default at every K — the
Casali pilot's worst-case (Burns2011/T4) is **neutral, not regressive**
once seed noise is controlled for, and the K = 11 regime shows the
intended pooling benefit. A future tidy-up should consider either a
centred fallback or a τ-aware joint move on z_c to suppress the stuck-
z_c tail; in the meantime, users running a small number of seeds on
K ∈ {2, 3} balanced partitions should be aware that any single seed
may show one class with ~50 % of expected ESS without the chain being
in trouble.

The earlier smoke test
(`dev/red-team/heavy-tests/funnel-stress-hyperprior-sigma.R`, 5-class
random-binary fixture) is preserved as a plumbing check — it confirms
both priors run end-to-end but its ratios are not representative.

Raw per-cell RDS:
`dev/red-team/heavy-tests/casali-ess-results/<matrix>__<treatment>__<prior>[__<tag>].rds`;
tabulated summary: `dev/red-team/heavy-tests/casali-ess-results.{csv,md}`;
per-rep diagnostic: `dev/red-team/heavy-tests/casali-ess-rep-summary.R`.

**Move-sampler validation.** `tests/testthat/test-partition-hyperprior.R`
group (F) runs the per-class and `scale_hyper_tau` MH moves at β = 0
(likelihood disabled) and verifies the empirical (τ, z_c) moments match
`HalfNormal(1)` to within 15 %/20 % (mean/sd). This catches Hastings-
ratio bugs or prior-aware acceptance regressions in the move dispatcher
that the analytic R↔C++ density tests would miss.

**Degenerate cases.**
* `unlink` without `"shape"`, or `nClasses == 1` (silently coerced by
  `.ValidatePartitionArgs` when the user passes `unlink = "shape"` at
  `nClasses == 1`): the legacy single-σ Gamma prior remains in effect.
  Numerical equivalence with `eval_log_prior_cpp` to ~1e-10 is
  guaranteed (`test-partition-hyperprior.R` (A)).
* `partition = NULL` (§7a): bit-identical to the unchanged legacy code
  path. Reference test `test-partition-bitcompat-null.R` still green.

## Partition API — Layer 1 complete

The `feature/partition-api` branch adds per-character user-class
partitioning to `RunMkPrime` (plan v4 in `NOTES/partition-api-plan.md`).
Layer 1 is now fully implemented and end-to-end tested.

### Shipped in Layer 1

* `RunMkPrime()` accepts new arguments `partition` (integer vector
  assigning each character to a 1:nClasses user class) and `unlink`
  (character vector of model components to unlink across classes).
  All five AutoPart / Casali production treatments are dispatchable
  end-to-end:

  | Treatment | Description | `partition` | `unlink` |
  |-----------|-------------|-------------|----------|
  | T0 | Unpartitioned | `NULL` | — |
  | T1 | Anatomical K-class | integer vector | `c("shape", "ratemultiplier")` |
  | T2a | AutoPart isolated | integer vector | `c("shape", "ratemultiplier")` |
  | T2b | AutoPart merged | integer vector | `c("shape", "ratemultiplier")` |
  | T4 | Random control | integer vector | `c("shape", "ratemultiplier")` |

* **§7a contract honoured:** `partition = NULL` routes through the
  unchanged legacy code path bit-for-bit. A regression test locks the
  reference sample from `Lobo.phy`.
* **§7b numeric-equivalence honoured:** at the trivial spec
  (`nClasses = 1`, `classRate = 1`, `etaNeo = 1`) and at `t = 0` of any
  spec, the partitioned likelihood surface matches the legacy
  `cpp_log_likelihood` to ~1e-10.
* **Per-class moves implemented:**
  * `scale_class_rate_log_sd` (case 31): Bactrian MH on each
    `classRateLogSd[c]` independently (unlink `"shape"`).
  * `dirichlet_simplex_class_w` (case 32): Dirichlet simplex proposal on
    the class-weight vector `(w_1, …, w_K)`, maintaining the simplex
    invariant `sum(w_c) = 1` (unlink `"ratemultiplier"`).
* **Per-class priors implemented:**
  * `class_w ~ Dirichlet(alpha)` with default `alpha = 1` (flat);
    tunable via `MkPrimeModel(classRateConcentration = ...)`.
  * `class_rate_log_sd[c] ~ Gamma(shape, rate)` i.i.d. across classes.
  * `eta_neo ~ LogNormal(0, rateNeoSdlog)` — single global scalar
    for neomorphic gain/loss asymmetry; **frozen at 1.0 in Layer 1**
    (Casali data are all transformational so `nNeo = 0` and this never
    bites; semantics for `eta_neo ≠ 1` with neomorphic data to be
    settled in a follow-up).
* `unlink` token resolution is case-insensitive with partial-prefix
  matching (warns on prefix; errors on ambiguous prefix or unknown
  token, with an `agrep`-driven "did you mean" suggestion). `unlink` is
  silently coerced to `character(0)` with an info alert when
  `partition = NULL`, so legacy callers passing `unlink` accidentally
  get a single uniform behaviour.

### Deferred to Layer 2

* `unlink = "brlens"` (per-class branch lengths under a shared topology,
  MrBayes subParam idiom) is not yet implemented. Passing `"brlens"`
  errors cleanly with an informative message.

### Determinism note (by design)

The MCMC engine is intentionally non-deterministic across separate R
processes at long chain lengths: the adaptive move scheduler tunes move
weights against wall-clock cost, which varies between runs. Within a
single process, identical seeds produce identical samples (validated by
`tests/testthat/test-determinism-baseline.R` and the §7a bit-identity
reference, which is generated and consumed in the same process). The
§7a reference uses `gibbsSubtreeSwap = FALSE, joint2d = FALSE,
autoTune = FALSE, nIter = 100, maxWarmup = 50` to keep the
bit-comparison robust against wall-clock-driven adaptation.

---

* `RunMkPrime(tree = NULL)` now starts from a greedy parsimony
  stepwise-addition tree (`TreeSearch::AdditionTree`) when `TreeSearch`
  is installed, falling back to `TreeTools::NJTree` otherwise.  The
  parsimony tree is a methodologically preferred MCMC start for
  morphological data, and `AdditionTree` is essentially instant
  (no iterative search).  Behaviour is unchanged when a `tree`
  argument is supplied.

* Multiple MCMC runs now execute in parallel when `nCore > 1`
  (defaults to `getOption("mc.cores", 1L)`). Backend: `callr::r_bg()`.
  Mirrors the TreeDist option-driven pattern.
* **Breaking (pre-release):** `parallel` argument to `MkPrimeMCMC()`
  removed; use `nCore` instead. The `future` package is no longer used.

* New `kPrimePrior = "empirical_geometric"` prior (now the default).  Decomposes
  `k' = N_obs + N_unobs`, with `N_obs` drawn from an empirical pmf tabulated
  from real morphological matrices (`empiricalNObs`) and `N_unobs` from a
  `Geometric(p)` with `Beta(a, b)` hyperprior on `p`.  The convolution of the
  two distributions defines the prior on `k'` and counters the previous
  tendency of inference to collapse the number of unobserved states to zero.
* New `MkPrimeEmpiricalPrior()` constructor and `empiricalNObs` package
  dataset, allowing users to override the empirical component.
* `gibbs_p` Gibbs draw on `p` is disabled under `"empirical_geometric"`
  (no Beta conjugacy under the convolution); replaced with a logit-scale
  Metropolis-Hastings move `mh_logit_p` (move type 30).  The multiplicative
  `mh_p` move (type 8) is rejected for most proposals when the posterior on
  `p` concentrates near 1 (typical under the empirical_geometric prior),
  causing the adaptive scheduler to crush its weight to the floor and `p`
  to stay stuck; proposing on the unbounded logit scale fixes this.
