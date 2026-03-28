# M-052: Discretized Dirichlet-marginal Q-matrix heterogeneity

**Status**: Planning (revised)
**Agent**: A
**Date**: 2026-03-28

## Objective

Implement a lightweight Q-matrix heterogeneity option for **all character
state counts** (not just binary), so users can compare topologies under
homogeneous vs heterogeneous models, following the manuscript recommendation
to "apply multiple models and report only relationships that are recovered
consistently."

## Evidence from the neotrans manuscript

The Smith & Yang manuscript tested heterogeneous Q-matrix models across
69 datasets. Summary of findings relevant to our design:

### Model fit
- Het1 (Beta mixture across ALL 2-state characters) outperforms
  Het (neomorphic only) in 28/4 datasets. Most of the fit
  improvement reflects latent rate variation across all binary
  characters, not processes specific to neomorphic traits.
- Further partitioning (separate Beta distributions for neo vs trans)
  offers negligible improvement.

### Tree topology (6 datasets with reference trees)
- In 4/6 datasets, heterogeneous models produced **worse** trees
  despite better model fit.
- In 1/6 (Halliday complete), Het models substantially improved trees.
- In 1/6, results were mixed.
- "No single class of models consistently produces trees that are
  closer to well-corroborated topologies."

### Computational cost
- 2–4× more time, 2× more memory.

### Conclusion
Het models are worth offering for model comparison but should **not** be
the default. The primary value is letting users check whether their
topology is sensitive to model choice.

## Reference

Wright, Lloyd & Hillis (2016, Syst. Biol. 65:602–611) introduced the
discretized Beta approach for binary characters. The RevBayes tutorial
`mkv_model_gamma_discretized_multistate.Rev` extends this to multistate
characters using Dirichlet-drawn frequency vectors as MCMC parameters
with reversed copies for symmetry.

Our design takes a different approach: instead of MCMC-sampled frequency
vectors (which add many parameters and mix slowly), we **discretize the
Dirichlet marginal** to keep a single scalar parameter for any k.

## Design

### Mathematical framework

Each character c with k_c states evolves under an F81 model with
equilibrium frequencies π_c drawn from a symmetric Dirichlet:

```
π_c ~ Dirichlet(α, ..., α)    [k_c components]
```

The per-character likelihood marginalizes over π_c:

```
L_c = ∫ L_c(π) × Dir(π | α,...,α) dπ
```

We approximate this integral using the **Dirichlet marginal
discretization** described below.

### Key insight: Dirichlet marginal

The marginal distribution of any single component of Dirichlet(α,...,α)
with k components is:

```
π_j ~ Beta(α, (k−1)α)
```

with mean 1/k and variance proportional to 1/α. This is a well-known
property of the Dirichlet distribution.

For k = 2, Beta(α, (2−1)α) = Beta(α, α) — the symmetric Beta. This
is exactly the distribution used in Wright et al. (2016) and our
previous binary-only design.

### Discretization

Discretize Beta(α, (k−1)α) into B equal-probability bins:

```
β_b = E[ X | X ∈ I_b ],   b = 1, ..., B
```

where I_b is the b-th quantile interval. This is the standard
fnDiscretizeBeta / fnDiscretizeGamma scheme.

### Rotation for state-labelling symmetry

Each bin β_b gives the frequency of one "special" state. To avoid
breaking the arbitrary state labelling, we cycle which state is
elevated. For each bin b and each state s ∈ {0, ..., k−1}:

```
π^(b,s)_j = β_b                if j = s
          = (1 − β_b)/(k − 1)  if j ≠ s
```

The approximate marginal likelihood is:

```
L_c ≈ (1/B) × (1/k) × Σ_b Σ_s L_c(π^(b,s))
```

This uses B × k equal-weight mixture components.

### Verified properties

1. **Uniform average root frequency.**
   For any α and any b: (1/k) × Σ_s π^(b,s) = (1/k, ..., 1/k).
   Each state appears as "special" exactly once, and the contributions
   from all rotations cancel to uniform. The model preserves the
   "no preferred state" property.

2. **Homogeneous limit (α → ∞).**
   Beta(α, (k−1)α) concentrates at 1/k. All bins → 1/k, all frequency
   vectors → (1/k,...,1/k), and F81 → JC(k). The Het model recovers
   the homogeneous Mk model.

3. **k = 2 backward compatibility.**
   Beta(α, α) is symmetric around 0.5, so bin b and bin B+1−b
   satisfy β_b = 1 − β_{B+1−b}. The rotation over s ∈ {0, 1}
   produces the same set of frequency vectors as the B bins alone
   (each appearing twice with half the weight). The sum reduces to:
   ```
   L_c ≈ (1/B) × Σ_b L_c(π^(b,0))
   ```
   which is exactly the current binary implementation. No code change
   needed for k = 2.

4. **Distinct from ACRV.**
   ACRV scales the overall substitution rate but keeps all off-diagonal
   transition probabilities equal (under JC). Het changes the
   *equilibrium frequencies*, producing F81 Q-matrices where transitions
   TO different states have genuinely different probabilities.
   These are orthogonal model components.

### F81 transition probabilities (closed form)

For F81 with frequencies π, the transition probability matrix has an
analytical solution:

```
P_ij(t) = π_j + (δ_ij − π_j) × exp(−μt)
```

where

```
μ = 1 / (1 − Σ_l π_l²)
```

For our "one-special-state" structure π = (β, r, ..., r) with
r = (1−β)/(k−1):

```
Σ π_l² = β² + (1 − β)² / (k − 1)
μ = 1 / (1 − β² − (1 − β)² / (k − 1))
```

This is O(1) per matrix entry — no eigendecomposition or matrix
exponentiation needed. The same efficiency as the JC analytical formula.

### Composition with `rate_loss` (neomorphic characters, k = 2)

For neomorphic characters, `rate_loss` sets the base gain/loss asymmetry.
The Beta bins modulate around this base:

```
gain_b = [1 / (1 + ρ)] × 2β_b
loss_b = [ρ / (1 + ρ)] × 2(1 − β_b)
```

This is algebraically equivalent to F81 with:
```
π₁ = gain_b / (gain_b + loss_b)
π₀ = 1 − π₁
```

For k ≥ 3 characters (transformational, known-k), there is no `rate_loss`
analog. The F81 model is parameterized purely by the frequency vector.

### Scope: all characters

The Het mixture applies to **all characters**, regardless of state count:

| Character type | k | Base frequencies | Rotation |
|---|---|---|---|
| Neomorphic | 2 | Asymmetric (from `rate_loss`) | Redundant (Beta symmetry) |
| Transformational k'=2 | 2 | Symmetric (0.5, 0.5) | Redundant |
| Transformational k'≥3 | k' | Uniform (1/k') | k' rotations |
| Known-k, k=2 | 2 | Symmetric | Redundant |
| Known-k, k≥3 | k | Uniform (1/k) | k rotations |

### Single parameter: `beta_scale`

The entire mixture is controlled by **one scalar** α (= `beta_scale`).

- Prior: `beta_scale ~ Gamma(shape, rate)`
- Move: scale proposal
- Default: OFF (`qHeterogeneity = FALSE`)

The discretized bins are a deterministic function of `beta_scale` and k.
When `beta_scale` changes, all bins for all state counts are recomputed.

### Equal bin weights — no Dirichlet on mixture probabilities

We fix the component weights at 1/(B × k) each (or equivalently, 1/B
per bin with a 1/k rotation average). This eliminates the
multi-dimensional Dirichlet on `matrix_probs` used by RevBayes. The
trade-off is slightly less flexibility for substantially fewer
parameters and better MCMC mixing.

### User interface

```r
MkPrimeModel(
  data,
  qHeterogeneity = FALSE,    # enable Dirichlet-marginal Q-matrix mixture
  nBetaCat = 4L,             # number of discretization bins per state count
  betaScaleShape = 1,        # Gamma prior shape for beta_scale
  betaScaleRate = 1,         # Gamma prior rate for beta_scale
  # ... existing parameters unchanged
)
```

### Likelihood computation

For a k-state character, the per-character likelihood becomes:

```
L_c = (1/R) × (1/B) × (1/k) × Σ_j Σ_b Σ_s L_c(rate_j, Q^(b,s))
```

where j ∈ {1..R} (ACRV), b ∈ {1..B} (frequency bins), s ∈ {0..k−1}
(rotations).

For k = 2, the rotation sum is redundant (see property 3 above), so
the effective computation is (1/R) × (1/B) × Σ_j Σ_b L_c(rate_j, Q^(b,0)).

### Computational cost

| k | Distinct Q's | Pruning passes (B=4, R=6) | Overhead vs JC |
|---|---|---|---|
| 2 | B = 4 | 24 | 4× |
| 3 | 3B = 12 | 72 | 12× |
| 4 | 4B = 16 | 96 | 16× |
| 5 | 5B = 20 | 120 | 20× |

Typical morphological matrices are ~80% binary, ~15% 3-state, ~5% 4+
state. The overall overhead is dominated by binary characters. Expected
wall-clock impact: 4–6× for typical datasets (consistent with the
neotrans finding of 2–4× for binary-only Het).

## Implementation steps

### 1. R layer: `MkPrimeModel()` and priors

- Add `qHeterogeneity`, `nBetaCat`, `betaScaleShape`, `betaScaleRate`
  parameters to `MkPrimeModel()`.
- Store in model object. Add `beta_scale` to the MCMC state when enabled.
- Add prior: `beta_scale ~ Gamma(shape, rate)` in `LogPrior()`.
- Add Scale proposal for `beta_scale` in `.BuildMoves()`.
- Update `print.MkPrimeModel()` to show Het status.

### 2. C++ bin computation: `compute_het_bins()`

For a given `beta_scale` (= α) and state count k, compute the B bins
of Beta(α, (k−1)α) using quantile midpoints:

```cpp
void compute_het_bins(double alpha, int k, int B, double* bins) {
    double a = alpha;
    double b = (k - 1) * alpha;
    for (int i = 0; i < B; i++) {
        double lo = R::qbeta((double)i / B, a, b, 1, 0);
        double hi = R::qbeta((double)(i + 1) / B, a, b, 1, 0);
        double p_lo = R::pbeta(lo, a + 1, b, 1, 0);
        double p_hi = R::pbeta(hi, a + 1, b, 1, 0);
        double denom = R::pbeta(hi, a, b, 1, 0) - R::pbeta(lo, a, b, 1, 0);
        bins[i] = (a / (a + b)) * (p_hi - p_lo) / denom;
    }
}
```

This generalizes the existing `compute_beta_bins()` (which hardcodes
k = 2, i.e. b = α) to arbitrary k.

### 3. C++ pruning: `pruning_f81_acrv_flat()`

New function that handles **any k** using the F81 closed-form P(t):

1. Takes precomputed bins for the relevant k.
2. For each (rate_j, bin_b, rotation_s) combination:
   - Constructs π^(b,s) from the bin value and rotation.
   - Computes μ = 1/(1 − β² − (1−β)²/(k−1)).
   - Computes P_ij(t) = π_j + (δ_ij − π_j) exp(−μ × rate_j × t).
   - Runs Felsenstein pruning.
3. Averages over all combinations and logs.

For k = 2, the rotation loop runs once (redundant by symmetry).
For k ≥ 3, it runs k times.

### 4. Dispatch in `mcmc_likelihood.cpp`

Update `cpp_partition_log_likelihood()`:

- For **all** character types when Het enabled: call
  `pruning_f81_acrv_flat()` with the appropriate k and base frequencies.
- For type 0 (neo): pass `rateLoss` to compose with bin frequencies.
- For type 1 (trans): for each k' sub-group, use k'-specific bins.
- For type 2 (known-k): for each k, use k-specific bins.
- When Het is disabled: existing `pruning_jc_acrv_flat()` unchanged.

Bin computation is per-unique-k. Cache bins for each distinct k present
in the dataset (typically k ∈ {2, 3, 4, maybe 5}).

### 5. MCMC engine: proposal and state management

- Add `beta_scale` to the flat state buffer.
- Add a Scale move for `beta_scale` (reuse existing scale proposal
  infrastructure).
- When `beta_scale` changes, recompute bins for **all** k values and
  recompute **all** partition likelihoods (they all depend on beta_scale
  through their k-specific bins).
- Cache bins per k in `McmcState` (small: B doubles per distinct k).

### 6. Tests

- **Unit**: Verify `pruning_f81_acrv_flat()` matches hand-computed
  likelihood for a 3-taxon tree with known bins, for k = 2, 3, 4.
- **Homogeneous limit**: beta_scale → ∞ should recover JC(k) log-
  likelihood to within tolerance, for all k.
- **k=2 backward compat**: Het on binary characters should produce
  identical log-likelihoods under the old and new code paths.
- **Rotation symmetry**: For k=3, relabelling states 0↔1 in the data
  should produce the same marginal likelihood.
- **Reversibility**: MCMC with Het enabled should not alter `rate_loss`
  recovery on simulated data generated under homogeneous model.
- **Integration**: Full `RunMkPrime()` with `qHeterogeneity = TRUE`
  completes without error on Vinther dataset (which has multistate chars).

### 7. Documentation

- Document `qHeterogeneity` in `MkPrimeModel()` man page.
- Add a note that this is for model comparison, not recommended as
  default, citing the neotrans findings.
- Document the Dirichlet-marginal discretization with mathematical
  derivation in a vignette or package-level documentation.

## Approximation quality

The "one-special-state" frequency vector π^(b,s) cannot represent
arbitrary Dirichlet draws. For k = 3, a true draw might be (0.5, 0.3, 0.2),
which is not expressible as any rotation of (β, (1−β)/2, (1−β)/2).

This is inherent to any finite discretization of a multivariate
distribution via its marginal. The approximation captures:
- The **marginal** variation in individual state frequencies (exactly).
- The **mean** frequency vector (exactly: always 1/k by rotation).
- The dominant mode of variation: one state elevated, rest depressed.

It does **not** capture simultaneous elevation of two states at the
expense of others. For the Mk model (where state labels are arbitrary
and carry no biological meaning), this is an acceptable trade-off:
the primary signal is "some states are more common than others," not
"states 0 and 1 are both common while state 2 is rare."

## What this does NOT include

- No Dirichlet on matrix_probs (equal weights only).
- No separate Beta distributions for neo vs trans (no HetB/HetM/HetBM).
- No MCMC-sampled frequency vectors (RevBayes multistate approach).
- No stepping-stone marginal likelihood estimation (separate task).

## Comparison with RevBayes multistate approach

| Aspect | MkPrime (this design) | RevBayes multistate |
|---|---|---|
| Parameters | 1 scalar (`beta_scale`) | n_cats × k frequency vectors + mixture weights per k |
| Mixture components | B × k per k (deterministic) | 2 × n_cats per k (stochastic) |
| Symmetry mechanism | Rotation averaging | Reversed copies |
| MCMC mixing | Fast (1D scale move) | Slower (simplex moves in high dimensions) |
| Flexibility | Moderate (marginal discretization) | High (arbitrary frequency vectors) |
| Extension to any k | Natural (Beta(α,(k−1)α)) | Natural (Dirichlet(α,...,α)) |
