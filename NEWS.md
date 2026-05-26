# MkPrime (development version)

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
