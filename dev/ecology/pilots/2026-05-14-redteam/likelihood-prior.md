# Likelihood + prior red-team audit (2026-05-14)

**Summary.** I'm confident the core ecology mixture math (`gamma_e`,
`trans_rate_factor`, per-edge per-character mixing, reference-ecology
skip, refEcology indexing in zMat) is internally consistent between
`src/mcmc_ecology.cpp` and `src/mcmc.cpp`'s ecology prior block. I have
concerns about (a) the simulator/model normalisation contract when
kEco > 2, (b) the silent dependence on a caller convention that
`z[, refEcology + 1] == 0`, and (c) the **wEdge root-edge handling
departs from the vignette spec** — code uses `0.5 * (margRoot +
margRootChild)`, spec says `1/K` at the root edge. I found the
previously-flagged R/C++ prior mismatches plus one new prior bug (R
rejects `pi0 = 0`/`pi0 = 1` strict; that's symmetric with theta) and one
simulator/model alignment bug (scalar theta in simulator vs per-ecology
theta_e in model). None of these obviously explain the +307 nat
posterior gap on their own; the wEdge root-edge departure and the
ecology rate hard-coded to 1.0 are the most likely culprits worth
chasing if accumulator-side debugging clears.

---

## Bugs / discrepancies

### B1. wEdge root-edge mismatch with spec
`src/mcmc_ecology.cpp:104-106, 234-241, 936-942`;
`vignettes/ecology-details.qmd:58-62`.

Spec: "with the root edge using the stationary distribution 1/K".
Code: every edge — including the root edge — uses
`0.5 * (margParent + margChild)` where `margParent` for the root edge is
the *posterior* marginal at the root node, not `1/K`. `bwdRoot` is
seeded as `1/K` (line 106) but then multiplied by `fwd[root]` and
normalised (lines 168-174), so `marg[root]` is the data-conditional
posterior, not the prior `1/K`. Off-spec, and probably biases wEdge at
the deepest internal edge — the very edge most likely to drive
phi/pi0 over-recovery.

**Verify:** for a simple 4-tip tree, hand-compute `marg[root]` vs `1/K`
and compare to `EcologyEdgeWeights` output on the root edge.

**Fix:** in `compute_ecology_node_marginals` (or in
`EcologyEdgeWeights`/the inline loop at line 940), special-case the
root edge to use `1/K` for the parent side. Alternative: update the
spec to match the code (the comment at lines 9-18 already documents
"mean of parent + child" as an intentional cheap symmetric blend; just
remove the "root edge uses 1/K" sentence from the vignette).

### B2. Simulator uses **scalar** theta; model uses per-ecology theta_e
`inst/ecology/simulations/sim3-simulate.R:180-184`.

Line 180: `theta <- mean(z == 1L) / max(1, mean(z != 0L))` is a scalar
averaged over all non-reference ecology columns; line 183 then loops
ecologies but always plugs in the same scalar `theta`. The model
(`gamma_e_compute`, `cpp_log_prior` lines 426-432) uses **per-ecology
theta_e**. For Sim 3 with `kEco = 2` (one non-ref ecology) this is
harmless. For any sim with kEco > 2 and non-symmetric z distribution
across columns, simulator and model disagree on the normalisation
constant — silent rate-scale bias.

**Verify:** call `.SimulateMkPrimeEcology(normalize = TRUE)` with
`kEco = 3` and z columns having different (enc/disc) proportions, then
compute `gamma_e` from the simulator's scalar and the model's per-column
values — they differ.

**Fix:** in sim3-simulate.R, allow `theta` to be a vector of length
`kEco - 1` and index by ecology when computing `gammaE[e]`. Mirror the
model's column-index convention (`j = (s < refEcology) ? s : s - 1`).

### B3. Simulator's `z` has refEcology column; model's `zMatrix` does not
`inst/ecology/simulations/sim3-simulate.R:32` (caller) and `:211`
(simulator read); `src/mcmc_ecology.cpp:692, 762, 947, 1154`.

The simulator takes `z` with `ncol(z) == kEco` (all ecology columns,
including refEcology). The model's `zMatrix` has `ncol == kEco - 1`. All
callers happen to zero out the refEcology column (`zFull[, 1] = 0`),
which makes refEcology behave as "no deviation" — matching the model's
"refEcology forced to factor 1". But this is an unenforced convention.
A caller who sets `z[c, refEcology + 1] = 1L` would silently produce
data inconsistent with any fit.

**Verify:** set `z[1, 1] = 1L` (refEcology = 0 → column 1) in a Sim 3
caller; the simulator will apply `mult[2] = phi` on refEcology edges,
producing a dataset the model cannot represent at all.

**Fix:** in `.SimulateMkPrimeEcology`, require `z` to have `kEco - 1`
columns and apply the same `j = (s < refEcology) ? s : s - 1` mapping
the C++ uses. Or assert `all(z[, refEcology + 1] == 0L)` at entry.

### B4. R `LogPrior` rejects `theta ∈ {0, 1}` and produces NaN at boundary
`R/MkPrimeModel.R:656, 668-677`.

Already flagged in the brief. Two related issues:
- Line 656: `any(theta <= 0) || any(theta >= 1)` returns `-Inf` at the
  boundary; the C++ guard (line 408) does the same, but the
  *simulator* allows `theta = 1` and the vignette (lines 122-124)
  explicitly says "admits the boundary cases".
- Lines 669-677: `nDiscCol * log1p(-theta) + nEncCol * log(theta)`
  produces `0 * -Inf = NaN` at `theta == 1` and `nDiscCol == 0`. The
  C++ guard at lines 428-430 of `src/mcmc.cpp` is correct (skips when
  the count is zero).

**Fix:** mirror the C++ guard in R:
```r
contrib <- 0
if (any(nNoneCol > 0)) contrib <- contrib + sum(nNoneCol * log(pi0))
if (any(nEncCol  > 0)) contrib <- contrib + sum(nEncCol  * (log1p(-pi0) + log(theta)))
if (any(nDiscCol > 0)) contrib <- contrib + sum(nDiscCol * (log1p(-pi0) + log1p(-theta)))
lp <- lp + contrib
```
For the strict-open bound on `theta`: either relax to `theta < 0 || theta > 1`
(matching the Beta(1, 1) default's closed support) or document that the
boundary is excluded by both prior and proposal — and remove the 0.999
nudge in the simulator-comparison harness. The current half-fix
(R rejects, C++ allows) means the R posterior reported during a chain
that briefly visits `theta = 1` from a tuned proposal will silently
become `-Inf` whereas C++ won't.

### B5. Ecology rate fixed to 1.0 in `compute_ecology_node_marginals`
`src/mcmc_ecology.cpp:205, 933, 1244` (all three call sites pass
`rateMultiplier = 1.0`).

Not strictly wrong — the ecology covariate is an exogenous "fixed"
trait — but the marginal-reconstruction smearing of `w_{p->c}(e)`
depends on this rate. With a long stem branch (`stemBr = 0.10` in
Sim 3) and the JC-K rate at 1.0, the marginal at the root will be
already substantially smeared toward `1/K`, which weakens whatever
ecology signal exists. **Suspicious as a contributor to the +307 nat
log-post gap** because the truth-init chain's wEdge will not match the
simulator's hard edge labels at all.

**Verify:** print `wEdge` at truth-init and compare to the
hard-`edgeEcology` matrix that the simulator used. If they are far
apart on deep edges, this is a model-vs-data mismatch baked into the
likelihood spec, not an accumulator bug.

**Fix (longer-term):** expose the ecology rate as a parameter, or use a
narrower prior that the data can update (data-augmentation MCMC over
edge labels is the principled route, mentioned in the file header).

---

## Looks correct

- **`gamma_e_compute` formula** (`mcmc_ecology.cpp:249-251`) matches
  the simulator's normalize block (line 183) and the vignette (line
  137) exactly.
- **`trans_rate_factor`** (lines 260-267) and `mkn_rates_for_state`
  (lines 470-492) implement the vignette's μ → μ/γ_e prescription
  correctly. Reference ecology returns factor = 1 unconditionally on z.
- **Mixture is applied per-(edge, character)**: lines 384-409 (JC) and
  599-619 (MkN) build `psMix`, `pdMix` (or P-matrix mix) inside the
  per-character loop, using the per-character z-row and per-edge w.
  Not pulled outside.
- **z-column indexing** when `s == refEcology` vs not is correct in
  all five places that do it (`pruning_jc_acrv_flat_ecology` lines
  391-398, 417-424; `pruning_mkn_acrv_flat_ecology` lines 603-609,
  625-631; per-char wrapper). The `(s < refEcology) ? s : s - 1`
  pattern is consistent.
- **zMatrix shape** is enforced `nChar × (kEco - 1)` at all entry
  points (`PruningJcEcology:761`, `PruningMknEcology:692`); the per-
  char copy at `cpp_log_likelihood_ecology:954-958` and the Gibbs
  scaffold at `:1308-1310` use the same shape.
- **C++ prior empty-count guard** (`src/mcmc.cpp:428-430`) correctly
  short-circuits the `0 * log(0)` case for boundary theta/pi0 when a
  column has zero counts in some category.
- **Ecology forward/backward** in
  `compute_ecology_node_marginals`: standard postorder Felsenstein +
  preorder upward-message; folds root prior `1/K` into `bwdRoot` so
  `marg = fwd * bwd / Z` is the posterior at every node. The sibling
  lookup (lines 109-115) correctly assumes a strictly binary tree.
- **JC mixture preserves JC structure**: each `(p_same, p_diff)`
  contribution is averaged with the same `w`, so the mixed matrix is
  still doubly-stochastic with a single diff-coefficient → the
  reduced O(k) propagation is valid. Good — the comment at lines
  272-274 calls this out.
- **MkN mixture is NOT collapsible to a single (pi0, pi1)** and the
  code correctly does a full 2×2 matvec (lines 617-618), not the
  short-cut MkN form.
- **Variable-coding correction** computes per-character constant-site
  probability with the same per-character z row, matching what the
  pruning pass actually used.
- **kEco = 1 degenerate case**: zMat has 0 columns, prior loop skipped
  in both R and C++; likelihood loop runs with w(refEco) = 1, factor =
  1 — equivalent to baseline Mk. ✓.

## Could not verify

- Whether the +307 nat log-post drift is *caused* by any of B1–B5
  rather than the suspected accumulator bug — I'd need to run the
  truth-init script with a `wEdge` audit print and compare against
  `0.5 * (margParent + margChild)` vs the hard `edgeEcology` matrix.
- Whether the R `LogPrior`'s `dbeta(theta, ...)` Beta(thetaAlpha,
  thetaBeta) hyperprior matches the C++ exactly when those defaults
  differ between fit-time and truth-init; I checked the math, not the
  default-resolution code path.
- `data.relabel` interaction with z under transformational characters:
  vignette identifiability note (lines 92-110 of sim3-simulate.R)
  applies to neomorphic only, but the relabel correction is added for
  transformational only (mcmc_ecology.cpp:1055-1059) — checked the
  symmetry argument briefly; nothing leapt out as wrong but I didn't
  exhaust it.
