# Beta-Geometric Per-Character kPrime Prior

**Created:** 2026-03-31
**Status:** DRAFT

## Problem

The current hierarchical geometric prior on k' (the true number of character
states) uses a **shared hyperparameter p** across all transformational characters:

```
k'_i | p  ~  Geometric(p)  truncated at kObs_i
p         ~  Beta(a, b)
```

This creates a structural over-shrinkage problem. The conjugate Gibbs update
for p concentrates it based on the **collective** behavior of all characters:

```
p | k'  ~  Beta(a + nTrans, b + Σ(k'_i - kObs_i))
```

When most/all k' = kObs (the common case), the posterior for p becomes
Beta(n+1, 1), giving p ≈ 1. This imposes a per-character prior penalty of
log(1-p) ≈ -9 log-units for kObs → kObs+1, making it nearly impossible for
any character to explore hidden states.

**Empirical evidence (benchmarks/02_kprime_ess_benchmark.R):**
- p = 0.9999 at equilibrium on Sun2018 (225 trans chars)
- int_walk kPrime achieves 0% effective acceptance (MH accept ≈ 0.01%)
- Even the Gibbs sweep visits kObs+1 only 0.1–1.2% of samples
- The data shows ΔlogLik of -0.02 to -2.7 for kObs+1, which is far smaller
  than the prior penalty of -9.08

The feedback loop:
1. Most k' = kObs → p concentrates near 1
2. p ≈ 1 → prior massively penalizes k' > kObs
3. No k' can leave kObs → back to step 1

This is a **model structure** problem, not an MCMC mixing problem.

## Proposed Solution: Per-Character p with Shared Hyperparameters

Replace the shared p with per-character p_i, marginalized analytically:

```
k'_i | p_i  ~  Geometric(p_i)
p_i         ~  Beta(α, β)        # per-character, marginalized out
(α, β)      ~  hyperprior        # shared, estimated from data
```

Marginalizing p_i yields the **Beta-Geometric** prior per character:

```
P(k'_i = kObs + u | α, β) = B(α+1, β+u) / B(α, β)
```

### Why this fixes the problem

| Model | Penalty (kObs → kObs+1) | Depends on n? |
|-------|-------------------------|---------------|
| Shared p, point (current) | -9.08 (n=225) | Yes — worsens with more chars |
| Shared p, collapsed | -5.42 (n=225) | Yes |
| **Per-character p_i** | **-1.10** (α=β=1) | **No** |

The per-character penalty is constant, determined only by (α, β). Characters
are coupled softly through the shared hyperparameters, not rigidly through a
single point p.

MH acceptance for kObs → kObs+1 improves from 0.01% to **5–33%** (measured
on Sun2018 with actual likelihood ratios).

### Relationship to ACRV

This is structurally analogous to how ACRV models rate heterogeneity:

| | ACRV | Beta-Geometric kPrime |
|--|------|----------------------|
| Per-character parameter | rate_i | p_i |
| Distribution | LogNormal(0, σ) | Beta(α, β) |
| Shared hyperparameter | σ (= rate_log_sd) | (α, β) |
| Marginalized? | Yes (discretized) | Yes (analytical) |

### Mathematical details

The Beta-Geometric(α, β) prior on u = k' - kObs:

```
P(u=0) = α / (α + β)
P(u=1) = α / (α + β) × β / (α + β + 1)
P(u=k) = α / (α + β) × ∏_{j=0}^{k-1} (β + j) / (α + β + 1 + j)
```

For the Gibbs sweep, the log-prior is computed **incrementally**:

```
logPrior(u=0) = log(α) - log(α + β)
logPrior(u=k) = logPrior(u=k-1) + log(β + k - 1) - log(α + β + k)
```

No `lgamma()` or `lbeta()` calls needed in the hot loop — just one addition
and two `log()` calls per candidate k-offset. This is cheaper than the current
geometric formula (which requires `log(p)` and `log(1-p)` lookups).

### Hyperparameter estimation

With p marginalized, (α, β) are estimated from the data. The log-posterior
for (α, β) given all k' is:

```
log P(α, β | k') ∝ log P(α, β) + Σ_i [lbeta(α+1, β+u_i) - lbeta(α, β)]
```

This requires MH updates for α and β (no conjugate form). But:
- Only 2 scalar parameters → negligible cost
- Can use Bactrian or scale proposals
- Updated once per cycle (like rate_log_sd)

**Default hyperprior:** α ~ Exp(1), β ~ Exp(1), giving a broad prior that
lets the data determine whether hidden states are common or rare. With
α = β = 1, the marginal prior puts 50% at kObs, 17% at kObs+1, etc.
(heavy-tailed relative to the current model).

### Identifiability note

The model math spec notes: *"Per-character u is not identifiable
(Spearman r < 0.31 across 15 priors at realistic tree lengths)."* This applies
to both the current model and the proposed one. The Beta-Geometric doesn't
improve identifiability — it fixes the **over-shrinkage** that prevents
the MCMC from exploring the (admittedly uncertain) posterior on k'.

## Implementation Plan

### Step 1: Add `kPrimePrior = "beta_geometric"` option

**R/MkPrimeModel.R:**
- Add `"beta_geometric"` to the `kPrimePrior` choices
- New parameters: `kprimeAlpha = 1`, `kprimeBeta = 1` (hyperprior on α, β)
- Store in model object: `kprimeAlpha`, `kprimeBeta`

### Step 2: R-side LogPrior

**R/MkPrimeModel.R (`LogPrior`):**

Replace the geometric block with a Beta-Geometric block when
`kPrimePrior == "beta_geometric"`:

```r
# Beta-Geometric: per-character p_i marginalized
# log P(k'_i | α, β) = lbeta(α+1, β+u_i) - lbeta(α, β)
u <- state$kPrime[transIdx] - mkd$kObs[transIdx]
alpha <- state$kprime_alpha
beta_  <- state$kprime_beta
lp <- lp + sum(lbeta(alpha + 1, beta_ + u) - lbeta(alpha, beta_))

# Hyperprior on (α, β): Exponential(1)
lp <- lp + dexp(alpha, rate = 1, log = TRUE)
lp <- lp + dexp(beta_, rate = 1, log = TRUE)
```

No `p` parameter, no Beta hyperprior on p.

### Step 3: C++ log prior (`cpp_log_prior`)

**src/mcmc.cpp:**

Add a third branch alongside `kPriorLogseries` / geometric:

```cpp
} else if (data.kPriorBetaGeometric) {
  double a = state->kprimeAlpha;
  double b = state->kprimeBeta;
  double lbAB = R::lbeta(a, b);
  for (int i = 0; i < nTrans; ++i) {
    int gi = data.transIdxGlobal[i];
    int u = kPrime[gi] - data.kObs[gi];
    lp += R::lbeta(a + 1.0, b + (double)u) - lbAB;
  }
  // Hyperprior: Exp(1) on α and β
  lp += R::dexp(a, 1.0, 1);
  lp += R::dexp(b, 1.0, 1);
}
```

### Step 4: McmcState / McmcData updates

**src/mcmc_state.h:**
- Add `bool kPriorBetaGeometric` to McmcData
- Add `double kprimeAlpha, kprimeBeta` to McmcState
- Remove `double p` from McmcState (or keep for backward compat; set unused)

**src/mcmc.cpp (`init_mcmc_state`, `get_mcmc_state`):**
- Pass and retrieve `kprime_alpha`, `kprime_beta`
- Skip `p` when using beta_geometric

### Step 5: Gibbs kPrime sweep update

**src/mcmc.cpp (`gibbs_kprime_sweep_impl`):**

The prior computation in the inner loop (currently lines 3258–3264) changes.
Currently:
```cpp
double logPrior_k;
if (isGeometric) {
  logPrior_k = logP + ko * log1mP;
}
```

Replace with incremental Beta-Geometric:
```cpp
double logPrior_k;
if (isBetaGeometric) {
  // Computed incrementally: see precomputation below
  logPrior_k = bgLogPrior[ko];
}
```

Before the k-offset loop, precompute the incremental log-prior:
```cpp
double bgAlpha = state->kprimeAlpha;
double bgBeta  = state->kprimeBeta;
std::vector<double> bgLogPrior(K_MAX_CAND);
bgLogPrior[0] = std::log(bgAlpha) - std::log(bgAlpha + bgBeta);
for (int ko = 1; ko < K_MAX_CAND; ++ko) {
  bgLogPrior[ko] = bgLogPrior[ko - 1]
                  + std::log(bgBeta + ko - 1)
                  - std::log(bgAlpha + bgBeta + ko);
}
```

This is O(K_MAX_CAND) precomputation, then O(1) per character per candidate.
The M-155 batched architecture and early termination are unchanged.

### Step 6: MH updates for (α, β) hyperparameters

**New move types** (27 = scale_kprime_alpha, 28 = scale_kprime_beta):

Simple Bactrian scale proposals, identical structure to `tree_length` or
`rate_log_sd`. Each proposes:
```
α_new = α_old × exp(λ × bactrian())
```

Accept/reject via the full log-prior ratio (all kPrime terms + hyperprior).
This is O(nTrans) per proposal — cheap, since it only evaluates the prior
(no likelihood change).

**Weight:** Low (1–2 per cycle each). These hyperparameters change slowly
and don't affect the likelihood directly.

### Step 7: Remove Gibbs p update for beta_geometric

When `kPrimePrior == "beta_geometric"`:
- Do not include `gibbs_p` (moveType 9) in the move pool
- Do not include `p` in paramNames / sample columns
- Include `kprime_alpha` and `kprime_beta` instead

### Step 8: int_walk kPrime — now viable

With the prior penalty reduced from -9 to -1.1, int_walk kPrime (moveType 7)
becomes effective. No code changes needed — the MH ratio naturally improves.
The adaptive scheduler will detect the improved acceptance and upweight it.

### Step 9: Block kPrime shift — reassess

May contribute now that the prior penalty is weaker. Monitor after the
change; the adaptive scheduler handles weight allocation automatically.

### Step 10: Parameter names and output

**R/RunMkPrime.R (`.ParamNames`):**
- When `kPrimePrior == "beta_geometric"`: replace `"p"` with
  `"kprime_alpha"` and `"kprime_beta"`

**R/RunMkPrime.R (`.BuildMoves`):**
- Add `scale_kprime_alpha` and `scale_kprime_beta` moves
- Condition `gibbs_p` on `kPrimePrior != "beta_geometric"`

### Step 11: Checkpoint compatibility

Existing checkpoints have `p` in the state. `ResumeMkPrime` should handle:
- Old checkpoint (geometric) → resume with geometric (unchanged)
- New checkpoint (beta_geometric) → resume with beta_geometric
- Mismatch → error with clear message

### Step 12: Tests

1. **Prior correctness:** Verify R and C++ LogPrior agree for beta_geometric
   with various (α, β, k') values.
2. **Gibbs sweep correctness:** On a small dataset, verify that the Gibbs
   sweep samples from the correct full conditional by comparing to brute-force
   enumeration. Mirror existing test-gibbs-kprime.R structure.
3. **Hyperparameter MH:** Verify acceptance rates are reasonable and that
   (α, β) converge to sensible values.
4. **int_walk kPrime acceptance:** Verify that int_walk now achieves non-zero
   acceptance on Sun2018 (the whole point).
5. **Benchmark regression:** Re-run benchmarks/02_kprime_ess_benchmark.R
   with beta_geometric and compare kPrime ESS/second across move configs.
6. **Backward compatibility:** Verify that `kPrimePrior = "geometric"` still
   works identically to before.

## Migration Strategy

- **Default prior:** Keep `"geometric"` as the default for now (backward
  compat). Add `"beta_geometric"` as an opt-in choice.
- **Once validated:** Consider making `"beta_geometric"` the default in a
  future release, with `"geometric"` retained for reproducibility.
- **Vignette:** Update hyoliths.qmd to demonstrate both priors and discuss
  the trade-off (shrinkage vs. flexibility).

## Files Modified

| File | Changes |
|------|---------|
| R/MkPrimeModel.R | Add `"beta_geometric"` option, `kprimeAlpha/Beta` params, LogPrior branch |
| R/RunMkPrime.R | `.ParamNames`, `.BuildMoves`, `.InitState` for new prior; new move types |
| src/mcmc_state.h | `kPriorBetaGeometric` flag, `kprimeAlpha/kprimeBeta` in McmcState |
| src/mcmc.cpp | `cpp_log_prior` branch, `gibbs_kprime_sweep_impl` prior formula, MH moves for α/β, `init_mcmc_state`/`get_mcmc_state` |
| src/mcmc_likelihood.cpp | `prepare_mcmc_data` — pass new flag/params |
| tests/testthat/ | New test file(s) for beta_geometric prior, Gibbs sweep, MH moves |

## Risks

- **Too permissive:** With α = β = 1, the Beta-Geometric prior puts 50% mass
  on k' = kObs. If the data doesn't support hidden states, the posterior may
  still favor k' = kObs but waste MCMC effort exploring alternatives. This is
  arguably correct (the data speaks) but slower to converge than the current
  model. Mitigation: (α, β) are estimated from data, so if the posterior
  favors k' = kObs, α will increase relative to β, naturally concentrating.

- **Identifiability of (α, β):** With per-character k' only weakly identified,
  (α, β) may be poorly constrained. This is fine for inference (marginalizing
  over uncertainty) but may cause slow mixing of α and β themselves.
  Mitigation: monitor ESS for α and β; consider stronger hyperprior if needed.

- **Interaction with logseries:** The beta_geometric replaces only the
  geometric prior. The logseries option is unaffected. If users want
  a "beta-logseries" variant, that's a separate extension.

## Expected Impact

On Sun2018 (or similar datasets):
- **int_walk becomes viable:** MH acceptance 5–33% (from 0.01%)
- **Gibbs sweep explores more:** Should visit kObs+1 at ~10–30% rate (from 0.5%)
- **Scalar ESS improves:** The 19% budget currently wasted on int_walk kPrime
  now contributes real mixing, indirectly improving overall ESS/second.
- **More accurate posterior on k':** Characters that genuinely have hidden
  states can explore them without being suppressed by the collective.
