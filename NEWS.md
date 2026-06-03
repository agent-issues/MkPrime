# MkPrime (development version)

## `empirical_geometric` honours `priorVariant = "unconditional"` (Model A)

The `empirical_geometric` arm now respects the model's `priorVariant`. Under
`"unconditional"` (Model A, used by the SBC harness) the prior is placed on the
full support `k' >= 2` with no per-character `Z_i(p)` truncation correction; the
`k' >= kObs` floor is enforced by the likelihood. Under `"conditional"` (Model B,
the default production prior) the EG-001 `Z_i(p)` normaliser is applied as before
(behaviour unchanged). The two agree whenever every `kObs == 2`, and differ by
exactly `sum_i log Z_i(p)` otherwise. Reuses the existing `unconditionalPrior`
McmcData field already honoured by the plain-`geometric` arm; both
`.LogPriorEmpiricalGeometric` (R) and the `empirical_geometric` branch of
`cpp_log_prior` (C++) gate the correction on it (R<->C++ parity tested). A
`.SampleNObsEmpirical` forward-draw helper is added for Model A SBC.

## `empirical_geometric` prior: per-character truncation normaliser (EG-001 fix)

`LogPrior` enforces `k'_i >= kObs_i` for transformational characters, so the
`empirical_geometric` k'-prior is now renormalised over that truncated support.
Both `.LogPriorEmpiricalGeometric` (R) and the `empirical_geometric` branch of
`cpp_log_prior` (C++) subtract the per-character normaliser

```
log Z_i(p) = log( sum_{k >= kObs_i} P(k | p) ) = log( 1 - sum_{m=2}^{kObs_i-1} P(m | p) )
```

Previously the untruncated convolution log-pmf was returned, biasing the
marginal on `p` by tens of nats on a 50-character dataset (`|log Z_i|` reaches
~4.7 nats at kObs=8, p=0.7). The term is a no-op at kObs=2, depends on `p` (so
it enters `mh_logit_p` acceptance and does **not** cancel there), and cancels in
fixed-`p` Gibbs k'-sweeps. Verified by independent re-derivation (renormalises
to 1 to ~1e-14; R<->C++ parity to ~1e-9) and a regression block in
`test-empirical-geometric-prior.R`.

**Posterior samples on `p` (and downstream summaries) under
`kPrimePrior = "empirical_geometric"` are not comparable across this fix.** The
section 7a partition bit-identity reference (`partition-bitcompat-null-ref.rds`)
was regenerated because the default model uses this prior; the legacy
`partition = NULL` code path itself is unchanged.

## `likelihoodMode = "marginal_k"`: Rao-Blackwellised k'_i (geometric arm only)

`MkPrimeModel()` gains a `likelihoodMode = c("sampled_k", "marginal_k")`
option. Under `"marginal_k"` the per-character `k'_i` state is
analytically marginalised out of the likelihood at every evaluation
(`logSumExp` over `u ∈ {0, .., u_max}` weighted by `P(u | p)`), so the
chain carries only `(tree, mu, sigma, p)` — no slow discrete coordinate.
The geometric arm is the only one supported in v1; the
`empirical_geometric` / `beta_geometric` / `logseries` arms abort with a
"§11 follow-up" message. Het + marginal-k and partition-API + marginal-k
abort similarly (plan §13). See
`dev/notes/2026-05-28-marginal-k-plan.md` for the full design.

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

**Funnel-stress benchmark.** `dev/red-team/heavy-tests/funnel-stress-hyperprior-sigma.R`
is a smoke test of the plumbing: it runs a 5-class fixture
(sizes 5, 50, 50, 50, 50) under both priors at nGen = 20000 and prints
per-class ESS on σ_c. The fixture uses random-binary data with weak
per-class signal so the pooling effect is modest; representative ESS
gains need to be measured against the Casali production cells that
motivated this change, where σ_small mixes catastrophically under the
old prior.

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

## Partition-rate parameterisation fix (audit Issue 1, 2026-05-27)

* **RB-style partition-rate normalisation.** `rate_neo` used to enter the
  likelihood as a one-sided neo-only multiplier (`neoEl = edgeLen *
  rate_neo`), giving an nChar-weighted mean partition rate `(n_neo·r +
  n_trans) / nChar` that equals 1 only when `r = 1`. `tree_length` thus
  lost its "expected substitutions per character on the average edge"
  interpretation on mixed (neo + trans) datasets — implicit inflation
  factor ~1.31 at the 635 `by_nt_9v` posterior median. Now uses the
  symmetric two-sided RevBayes formula

      neoScale   = r / (1+r) · (n_neo + n_trans) / n_neo
      transScale = 1 / (1+r) · (n_neo + n_trans) / n_trans

  so `n_neo·neoScale + n_trans·transScale ≡ n_neo + n_trans` identically.
  When `n_neo == 0` or `n_trans == 0`, rate_neo has no identifiable effect
  and both scales collapse to 1.0 — **trans-only and neo-only datasets
  are bit-identical pre-/post-patch**, and the existing §7a bit-identity
  reference (Casali production workload, `hasNeo == false`) is unchanged.
* **Mixed-dataset posteriors will shift.** On any dataset with both neo and
  trans characters, posterior samples of `tree_length` and `rate_neo` are
  not directly comparable pre- and post-patch (the marginal on each
  shifts by a deterministic rescaling). The shape of the posterior on
  `rate_loss`, `kPrime`, topology, and branch lengths is unaffected to
  first order; check before reusing any pre-patch chains.
* **F1: partial-cache invalidation widened.** `do_move_impl` case 3
  (rate_neo scale) and case 18 (neo_joint = rate_loss + rate_neo), plus
  the `paramIdx == 3` branch of `eval_slice_target` / `slice_scalar_impl`,
  used to refresh only `data->neoPartIndices` after a rate_neo change.
  Under the new normalisation rate_neo also shifts `transScale`, so all
  partitions go stale — these paths now fall through to full
  recomputation. Node-CL cache invalidation widened analogously
  (`invalidate_all_cls` instead of `invalidate_neo_cls`).

## Numerical & ascertainment corrections (red-team campaign 2026-05-26)

* **FAST-EXP-001 fix.** Closed-form JC / MkN / F81 transition probabilities
  now use `std::expm1` to compute `1 − exp(arg)`, avoiding catastrophic
  cancellation when `rate × branch` is tiny. Worst observed relative error
  on the audit grid dropped from ~0.39 (at rt = 1e-15) to ≤ 2.4e-16.
  Touches 27 call sites across `rate_matrix.cpp`, `likelihood.cpp`,
  `acrv.cpp`, `ascertainment.cpp`, `mcmc_likelihood.cpp`,
  `gibbs_partial_cl.h`, and `node_cl_cache.h`. Sampled values change at
  ULP scale on every iteration — the §7a bit-identity reference was
  regenerated (commit accompanying this entry).
* **LIKE-001 fix.** Under `coding = "informative"`, the singleton-site
  ascertainment term is now added in `const_site_prob_for_k`
  (`mcmc_likelihood.cpp`) and at the three partial-CL sites in
  `node_cl_cache.h::cache_total_loglik` (type 0 / 1 / 2 branches).
  Previously the Gibbs k′ sweep and cache-accelerated MH moves dropped
  the singleton term, producing a chain that did not target the desired
  posterior under informative coding. Default `coding = "variable"`
  unaffected.
* **CONV-002 fix.** `R/Convergence.R::.ComputeRhat` now tail-equalises
  per-run sample matrices to the minimum nrow before `cbind`,
  eliminating a latent silent-recycle bug exposed by adaptive truncation
  or partial run drops.

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
