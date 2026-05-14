# v2 asymmetric-slab prior fix + γ-ridge diagnosis

## Status (2026-05-14, identifiability resolved)

Steps 1, 2, 3 are **complete**. Bugs A+B fixed; 20k Sim 3 pilot confirms the
γ-ridge tail is resolved (max tree_length 7,496 → 70; cor(log tl, log γ_e)
+0.629 → +0.197). Step 4 (v2-no-γ) is **not needed**.

The label-switching identifiability (φ ↔ 1/φ, θ ↔ 1-θ) surfaced by the
pilot is **resolved via Option 1 (post-hoc relabel)**, verified by:

- **Mode persistence:** two independent 100k Sim 3 chains (default and truth
  init) show **zero** within-run mode crossings.
- **`RelabelEcology()`** function implemented in `R/RelabelEcology.R` with
  27/27 tests passing. Smoke test on the 20k pilot (which was entirely in
  the reflected mode) recovers truth-scale parameters: phi 0.28 → 3.58,
  theta_1 0.03 → 0.97, z labels enc↔disc swap correct, z=none unchanged,
  function idempotent.

Step 5 (multi-rep Sim 3 HPC dispatch) is now **unblocked**. Step 6 (rodent
re-run + vignette) is next after the model spec stops moving.

## Original status (2026-05-13, pre-pilot)

v2 ecology-aware MkPrime is *functionally* running on Sim 3 but exhibits a
heavy tail in `tree_length` (median ≈ 22, max ≈ 7500, truth ≈ 9). Investigation
showed this is partly explained by γ-normalisation coupling tree_length to
(φ, π₀, θ) — and partly by two outright **bugs** in the v2 implementation that
must be fixed before the γ-coupling is re-evaluated.

## What v2 was supposed to be

Three changes over v1:

1. **γ-normalisation** of per-edge per-character rate factors so that the
   prior-expected rate on each non-reference ecology equals
   `γ_e = π₀ + (1-π₀)[θ_e·φ + (1-θ_e)/φ]`.
   Goal: close the v1 rate-time identifiability ridge.
2. **Reference-ecology dummy coding** — z-matrix has `(K-1)` columns, not `K`.
   Goal: remove K-vs-K-1 redundancy.
3. **Asymmetric slab** — Beta-distributed θ controls the relative probability
   of z = encouraged vs discouraged.
   Goal: avoid the unrealistic symmetric-slab assumption.

## What v2 actually is, today

| change | status | location | notes |
|---|---|---|---|
| γ-normalisation | ✅ implemented | `mcmc_ecology.cpp:249` (`gamma_e_compute`) | used in every per-edge rate factor |
| refEcology coding | ✅ implemented | zMatrix is nChar × (kEco-1) | |
| asymmetric slab (likelihood) | ✅ implemented | θ feeds into γ_e | |
| asymmetric slab (prior) | ❌ **broken** | `cpp_log_prior` mcmc.cpp:410–413 | uses `(1-π₀)/2` for slab regardless of θ |
| θ sampling | ✅ exists | move 33, `scale_theta` in mcmc.cpp:4562 | |
| θ logged | ❌ **missing** | `.ParamNames` RunMkPrime.R:3389–3396 stops at pi0 | also `.StateToRow:3438` |
| (rate × tree_length) joint scaler | ✅ exists | RunMkPrime.R:2997, 3002 (joint2d) | the pre-existing Mk ridge IS handled |
| (φ, tree_length) joint scaler | ❌ not present | — | candidate if pilot evidence demands it |

### Bug A — asymmetric-slab prior is symmetric

`cpp_log_prior` (mcmc.cpp:410–413):

```cpp
lp += R::dbeta(pi0, data.rho0Alpha, data.rho0Beta, 1);
if (nNone > 0) lp += nNone * std::log(pi0);
if (nSlab > 0) lp += nSlab * (std::log1p(-pi0) - std::log(2.0));
```

The `nSlab` term lumps encouraged and discouraged with prior `(1-π₀)/2` each.
The correct asymmetric form (matching what γ_e expects):

- `P(z = enc)  = (1-π₀) θ_e`
- `P(z = disc) = (1-π₀) (1-θ_e)`
- `P(z = none) = π₀`

Plus the Beta(α_θ, β_θ) prior on each θ_e.

**Consequence:** θ enters the likelihood (via γ_e) but has no proper Beta
prior counterweight on the z-cell labels. θ wanders freely → γ_e moves → drags
tree_length. Posterior interpretation of θ is also broken: it's just a
γ-knob, not a slab-shape parameter.

### Bug B — θ not in the MCMC log

`.ParamNames` (RunMkPrime.R:3389–3396):

```r
if (isTRUE(ecologyAware) && nPhi >= 1L) {
  if (nPhi == 1L) nms <- c(nms, "phi") else nms <- c(nms, paste0("phi_", ...))
  nms <- c(nms, "pi0")
}
```

θ has length `kEco - 1`; should append `paste0("theta_", seq_len(kEco - 1))`
after `pi0`. `.StateToRow:3438–3442` likewise. The C++ batch sample-matrix
path may also need extension — to verify.

## γ-ridge diagnosis (pre-fix)

On the existing v2 Sim 3 aware chain (200k iter, post-burnin n=194):

```
phi          median=0.862  IQR=[0.327, 1.272]
pi0          median=0.414  IQR=[0.216, 0.632]
tree_length  median=22.00  IQR=[11.06, 124.78]  max=7496.19
gamma_e(θ=0.5 proxy)  median=1.130  IQR=[1.006, 2.174]

tree_length / gamma_e   median=14.49  IQR=[10.43, 41.21]  max=4009.94

cor(log tl, log phi) = -0.027
cor(log tl, log gE)  = +0.629
```

γ-normalisation **does** couple tree_length to (φ, π₀, θ) — the +0.63
correlation is the smoking gun — but dividing tl by γ_e doesn't collapse the
tail to truth, so γ is a major contributor, not the sole cause. The unlogged θ
likely contributes the rest by drifting under the symmetric prior.

## Plan

Sequence chosen with the advisor: fix bugs first, then re-evaluate the ridge.

### Step 1 — fix the asymmetric-slab prior (model correctness)

In `cpp_log_prior` (src/mcmc.cpp:410–413), and the mirror in R (`LogPrior` in
R/MkPrimeModel.R:391):

- Count z-cells per ecology e: `nEnc_e`, `nDisc_e`, `nNone_e`.
- Slab contribution: `nEnc_e · log[(1-π₀) θ_e] + nDisc_e · log[(1-π₀)(1-θ_e)]`.
- None contribution: `Σ_e nNone_e · log π₀`.
- Beta prior on θ_e: `R::dbeta(theta_e, alpha_theta, beta_theta, 1)`, summed
  over e.

New hyperparameters `thetaAlpha`, `thetaBeta` (already declared in
MkPrimeModel.R per existing model spec; check actual defaults `2, 2`).

Verify R-side `LogPrior` matches C++. Add test pinning numerical equality
between cpp_log_prior and LogPrior for a fixed state.

### Step 2 — log θ (instrumentation)

- `.ParamNames`: add `nTheta` argument; emit `theta_1..theta_{nTheta}` columns
  after `pi0`. Wire through call site at RunMkPrime.R:234.
- `.StateToRow`: include `as.numeric(state$theta)` in the `ecoVal` block.
- Verify C++ batch sample-matrix layout matches (search for where samples
  matrix columns are filled per iteration).
- Update `sim3-v1v2-compare.R` post-processing to use the logged θ instead of
  the 0.5 proxy.

### Step 3 — 20k pilot Sim 3 with fixes

Single chain, simple Bactrians (no new joint scalers yet). Look at:

- `tree_length` trace: is the heavy tail gone, attenuated, or unchanged?
- Empirical posterior covariance of (log φ, logit π₀, logit θ, log tree_length).
- ESS on each scalar.
- `tree_length / γ_e` should be tightly concentrated near truth if the fixes
  are sufficient.

**Decision rule:**

- Tail gone, cov benign → ship simple Bactrians; multi-rep Sim 3 next.
- Tail attenuated but γ-tl cov still > 0.4 → implement v2-no-γ option as an
  A/B comparator. Compare on Sim 3 and pick the better-mixing parameterisation
  before the empirical run.
- Tail unchanged → something more fundamental is wrong; reconsider model.

### Step 4 — conditional: v2-no-γ option

Only if Step 3 shows γ-coupling persists post-prior-fix. Spec:

- Model flag `MkPrimeModel(gammaNormalise = TRUE/FALSE)`.
- In `gamma_e_compute`, return `1.0` when disabled.
- Keep refEcology coding + asymmetric slab. Reference ecology with rate=1
  anchors the timescale; asymmetric slab carries the directionality.
- Risk: re-opens the rate-time ridge if reference-ecology anchor is weak
  (i.e. too few reference edges). Diagnostic: `cor(log φ, log tree_length)`
  on the new chain. If >0.4 → add a joint `(φ, tree_length)` Bactrian.

### Step 5 — multi-rep Sim 3

Once the model is solid: 20 replicate datasets at the Goldilocks config
(nEco=60, nBase=180, φ=4, stem=0.10, root=0.15), each with blind + aware
chains. Headline metrics: P(true split), P(wrong split), CID gap, logL gap
(should be ≥ 0 under v2).

Dispatch decision: laptop if total time < 6h, HPC otherwise. At ~30 min per
chain × 2 chains × 20 reps = 20h total → HPC.

### Step 6 — rodent vignette

Background empirical run (bnqikq7l4) is using the *broken* prior. When that
completes, results will be partially wrong (the θ posterior is meaningless;
the z posterior is biased toward symmetric calls). Decide whether to re-run
post-fix or to add a caveat. Re-run is cleaner — empirical chain is ~3h.

## Background runs in progress

- `bnqikq7l4` — rodent empirical, 500k iter, ~3h. **Will need re-running**
  after Step 1.
- `ba0lf2yrk` — long Sim 3 v2 chain, 500k iter, ~85 min. **Also uses broken
  prior**; may be discarded or kept as a "before" snapshot for the diagnosis
  writeup.

## Critical files

To modify:

- [src/mcmc.cpp](src/mcmc.cpp) — `cpp_log_prior` asymmetric slab fix
- [R/MkPrimeModel.R](R/MkPrimeModel.R) — `LogPrior` mirror; `thetaAlpha`,
  `thetaBeta` defaults
- [R/RunMkPrime.R](R/RunMkPrime.R) — `.ParamNames`, `.StateToRow`, call sites
- [src/RcppExports.cpp](src/RcppExports.cpp) + [R/RcppExports.R](R/RcppExports.R) —
  if cpp_log_prior signature changes

To create:

- `tests/testthat/test-ecology-prior.R` — pin numerical equality of
  cpp_log_prior and LogPrior under asymmetric slab.

## Reused machinery (do not rewrite)

- 2D joint Bactrian: RunMkPrime.R:2997, 3002 — pattern to follow if a
  (φ, tree_length) scaler becomes necessary in Step 4.
- `gamma_e_compute`: mcmc_ecology.cpp:249 — single function; wrap in a
  `useGamma` flag for Step 4.
- Streaming log buffer: streaming.R `.AddToStreamBuffer` — flushes
  paramNames-shaped rows; just extend nParams.

## Step 3 pilot result (2026-05-13)

Full report: `dev/pilots/2026-05-13-step3-sim3-pilot/REPORT.md`.

20k iter, seed=20260512, kEco=2, Goldilocks (nEco=60, nBase=180, φ_truth=4,
stem=0.10, root=0.15), 175 raw samples → 131 retained, ~4.4 min.

| Metric | Pre-fix | Post-fix | Truth |
|---|---|---|---|
| tree_length median | 22.0 | 29.0 | 12.7 |
| tree_length max | 7,496 | **70** | — |
| tl / γ_e median | 14.49 | **9.35** | ~9 |
| tl / γ_e max | 4,010 | **28.8** | — |
| cor(log tl, log γ_e) | +0.629 | **+0.197** | ~0 |

All posterior covariance off-diagonals < |0.13|; no ridge detected. ESS
3.7–10.2 per 131 retained — low at 20k but adequate for diagnosis. Multi-rep
must use nIter ≥ 100k.

**Decision-rule branch: tail GONE, cov benign.** Bug A + Bug B were sufficient
to resolve the γ-ridge. Step 4 (v2-no-γ) is not required.

## Open question: label-switching (φ ↔ 1/φ, θ ↔ 1-θ)

Pilot recovered `phi ≈ 0.236` (truth 4) and `theta_1 ≈ 0.030` (truth ≈ 0.97
under the convention "encouraged"=high-rate). This is the reflected mode of
the likelihood symmetry

  (φ, θ_e, enc-cells, disc-cells)  ↔  (1/φ, 1 − θ_e, disc-cells, enc-cells).

The asymmetric-slab construction *names* one direction "encouraged" but does
not assign it a privileged prior weight: Beta(2, 2) on θ_e and LogNormal(0,
0.5) on log φ are both symmetric, so the posterior is bimodal with equal
posterior mass at the two reflections.

Consequences:

- **Topology/tree-shape inference:** unaffected. The mixture transition
  probabilities under (φ, θ) and (1/φ, 1-θ) are identical, so likelihood and
  branch lengths are recovered correctly. The 107× tree_length-tail
  improvement holds regardless of which mode the chain visits.
- **Per-character z interpretation:** flipped under the reflected mode. If
  the chain settles in the φ < 1 mode, posterior "encouraged" labels mean
  *biologically discouraged*. For the methods-paper narrative this is a
  presentation problem, not an inference problem.
- **Mixing across modes:** the pilot stayed in one mode; multi-rep will see
  different replicates settle in different modes unless we break the
  symmetry. This shows up as bimodality in cross-rep summary plots of φ.

### Options, in increasing intrusiveness

1. **Post-hoc relabel** every chain to the φ ≥ 1 representative. Cheap;
   preserves the symmetric prior; documented as a presentation choice in the
   vignette. Per-character z labels are inverted after relabelling if the
   chain was in the φ < 1 mode.
2. **Constrain `phi ≥ 1`** in the proposal / prior (e.g. half-LogNormal on
   log φ ≥ 0; or `phi = 1 + exp(log_phi_excess)` with log_phi_excess
   unconstrained). Eliminates the bimodality at the cost of a less natural
   prior shape; θ_e then carries the "what fraction of slab cells are
   strongly affected" interpretation unambiguously.
3. **Asymmetric Beta prior on θ_e** with α_θ > β_θ (e.g. Beta(3, 1.5),
   E[θ] = 0.67), reflecting the modelling claim that ecologies more often
   encourage than discourage convergent change. Soft symmetry-break: the
   reflected mode still exists but has lower posterior mass.
4. **Per-character orientation indicator** — overkill for a global slab
   shape parameter; would explode parameter count.

Initial recommendation was **option 2** (constrain phi≥1). Advisor pushed
back: option 1 is strictly less work, standard practice for label-switching
in mixture models, and doesn't require a model-spec rework. Decision
discriminator: does a single chain stay in one mode for the whole run? If
yes, option 1 is sufficient.

### Verification (2026-05-14)

Two 100k-iter Sim 3 chains run as a single-chain mode-persistence check:

| Chain | Init | Settled mode | phi median | theta_1 median | Zero-crossings |
|---|---|---|---|---|---|
| A (seed 20260513) | default | φ ≥ 1 (truth) | 3.60 | 0.94 | **0** |
| B (seed 20260514) | truth (phi=4, θ=0.97) | φ ≥ 1 (truth) | 3.64 | 0.99 | **0** |
| 20k pilot (seed 20260512) | default | φ < 1 (reflected) | 0.28 | 0.03 | **0** |

Different default-init seeds find different modes (RNG-dependent), but
**within a chain, modes do not mix** — exactly the condition where option 1
suffices.

### Option 1 implementation (2026-05-14)

`R/RelabelEcology.R` exposes `RelabelEcology(result, magnitudeMode)`. Per
sample, if `phi < 1` (global) or `phi_e < 1` (per-ecology), applies
`(phi → 1/phi, theta_e → 1-theta_e, swap z=1 ↔ z=2)`. Streaming-mode results
read samples from `result$logFile` on demand. Idempotent (post-flip every
phi ≥ 1, so re-application is a no-op). 27/27 unit tests pass; real-data
smoke on the 20k reflected-mode pilot recovers phi 0.28 → 3.58,
theta_1 0.03 → 0.97, z=enc/dis totals swap correctly (836 ↔ 31025), z=none
unchanged (9964).

Pipeline guidance for Step 5: after each chain finishes, pass the result
through `RelabelEcology()` before aggregating across replicates. This
guarantees consistent (φ≥1, θ near truth) interpretation regardless of
which mode each individual chain settled in.

## Out of scope

- Per-character θ (would explode parameter count).
- Per-ecology φ_e *and* per-ecology θ_e simultaneously (identifiability).
- Hierarchical character dependencies (deferred extension).
- Tip ecology uncertainty.
