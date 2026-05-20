# v2 ecology-aware redesign — red-team notes

Red-team review of the v2 spec (rate normalisation by gamma, reference
ecology, asymmetric slab) and of existing v1 simulation artefacts, while
the v2 implementation is blocked by sub-agent permissions.

## Verified-correct items

### γ normalisation math

For a non-reference ecology e and a single cell (c, e) under the
prior on (z, π₀, θ_e):

```
E[μ(z, φ) / γ_e]
  = (1/γ_e) · {π₀·1 + (1-π₀)·θ_e·φ + (1-π₀)·(1-θ_e)·(1/φ)}
  = (1/γ_e) · γ_e = 1     ✓
```

So the prior-expected per-cell rate factor on non-reference ecology is
exactly 1, and `tree_length` recovers its standard Mk interpretation
on reference-ecology edges. Identifiability ridge is closed.

### Reference-ecology arbitrariness

Swapping the reference ecology gives an equivalent posterior up to
relabelling: ((π₀, θ, φ), z) → ((π₀, 1−θ, 1/φ), 1−z) maps every
realised rate ratio to its reciprocal, which is just relabelling the
reference. Tree posterior, P(z ≠ none) per cell, and |log φ| are all
invariant.

### Identifiability of z direction for neomorphic

Neomorphic characters are never relabelled (0 = absent, 1 = present
biologically). So z = enc (more gains, fewer losses) and z = disc
(fewer gains, more losses) are genuinely distinguishable from the
data via per-clade state-1 frequencies. The v1 vignette's caveat
about "magnitude recovered, sign not" was wrong category-wide
(stale relic from a different model form).

### Identifiability of z direction for transformational

The off-diagonal rate matrix pattern is symmetric in state labels,
but total transition probability over branch length t is monotone
in the rate factor φ. So z = enc (factor φ) and z = disc (factor
1/φ) give different total amounts of change and are distinguishable
data-side.

### Sim 2's φ pull diagnosis

Re-checked: tip ecology was sampled `c(0L, 1L)` with replacement
from 16 tips and only constrained at ≥5 of each. With ecology-1
common at tip level, Fitch-assigned edge ecology gave 26/30 edges
in ecology 1 — yes, this is the rate-vs-baseline degeneracy.
Under v2 with γ normalisation, this should improve because
tree_length is anchored by the few edges with reference (=ecology 0)
state, but identifiability of φ separately from baseline rate
remains weak when one ecology dominates the tree.

### Sim 3 mixing diagnosis

Aware-chain logL was 76 nats WORSE than blind despite the model
nesting blind (set all z = none → recovers blind exactly). This is
a definitive signature of mixing failure, not model fit. The
diagnosis (φ–tree_length identifiability ridge) is correct because:

- Without normalisation: doubling all branch lengths and halving the
  factor on z=enc edges gives the same probabilities.
- With normalisation: same scaling on reference edges (factor 1)
  fixes baseline × tree_length, breaking the ridge.

v2 directly targets this.

## Concerns surfaced by red-team

### Concern 1: simulator/model mismatch under γ normalisation

The forward simulator (`sim3-simulate.R`) generates data with rate
factor `μ(z, φ)` un-normalised. The v2 model fits with rate factor
`μ(z, φ)/γ_e` on non-reference ecology. So the model parameters
under which the simulator generated data and the model parameters
that best fit those data are related by reparameterisation:

- True φ_sim = 4, true π₀_sim = 0.75 (180/240 cells with z=none in
  ecology 1)
- Model wants: γ_e · (rate ratio on z=none cells) = 1 — so γ_e ≈ 1
  if the z=none cells really do show rate 1 over reference. That
  pins γ_e = 1, which forces either π₀ → 1 or φ → 1 (contradicting
  true φ = 4).

**Conclusion:** under v2, simulating un-normalised and fitting
normalised produces a model misspecification. The chain has to
compromise on (π₀, θ, φ) such that the rate ratios fit per-cell
even though no single setting matches simulator truth exactly.

**Tree topology is unaffected** — rate ratios are right per-cell,
and topology is identified by those ratios. But the *interpretation*
of the (π₀, θ, φ) posterior under v2 from un-normalised simulator
data is not "the true generating parameters."

**Options:**

1. **Document and accept.** The model's reparameterisation is
   different from the simulator's, but the tree inference works.
   Add a vignette note.
2. **Normalise the simulator too** (add `normalize = TRUE` flag,
   default off but tested on). Then simulator-truth and model-posterior
   live on the same scale; recovery sims directly test "can we
   recover (π₀, θ, φ) = simulated values."
3. **Both:** keep `normalize = FALSE` as the realistic biology, add
   `normalize = TRUE` as a controlled recovery test.

Recommendation: option 3. The realistic-biology test (option 1) is
the headline argument; the normalised recovery test (option 2)
provides clean parameter-recovery evidence.

### Concern 2: ascertainment correction under γ

`const_site_prob_jc_eco_single` and the corresponding neomorphic
ascertainment helper compute P(all states equal | character variable).
This must use the γ-normalised rates too, otherwise the marginal
likelihood is internally inconsistent. The v2 implementation plan
notes this (in `mcmc_ecology.cpp:756` area).

If the implementation agent forgets, the symptom is: posterior on
tree_length drifts systematically off-truth in a way the parameter
posteriors don't explain. Add a regression test that checks the
ratio `L_v2 / L_blind → 1` as `(π₀ → 1, all z → none)` to catch
this.

### Concern 3: θ_e identifiability with small effect counts

θ_e is the slab balance: P(z = enc | non-none) = θ_e. If only ~60
cells have non-none z (as in sim3), and they're all z = enc (true
direction), the posterior on θ should concentrate at ≈ 1 — but the
Beta(2, 2) prior pulls toward 0.5. Effective posterior ≈
Beta(60 + 2, 0 + 2) → mean ≈ 0.97.

If the number of non-none cells is much smaller (real data with
diffuse signal), the prior may dominate and θ stays near 0.5,
which inflates γ_e in ways that compete with φ. Worth a sensitivity
analysis.

### Concern 4: per-ecology magnitudeMode with K = 2

Under `magnitudeMode = "per_ecology"`, there's one φ_e per ecology
state. With v2's reference convention, φ_ref ≡ 1 inert. For K = 2,
that means only ONE free φ — identical to `magnitudeMode = "global"`.
The two modes coincide for binary ecology. Worth documenting so users
understand it's not a free choice when K = 2.

### Concern 5: vignette `ecology-details.qmd` simulation evidence

The simulation evidence section in the math vignette describes Sim 1,
Sim 2, Sim 3 results obtained under v1 parameterisation. The headline
findings (π₀ ≈ 0.38 under sim3 aware) are valid for v1; under v2
they will change. After v2 sims complete, re-run and update the
vignette numerics.

## Open questions / follow-ups

1. After v2 implementation lands: re-run Sim 3 and compare aware
   logL to blind. v2 prediction: aware logL ≥ blind (rate-time ridge
   closed → no mixing failure).
2. Check whether `const_site_prob_jc_eco_single` has been updated to
   use γ — add a small unit test for the all-z = none case (should
   reduce to v1 ascertainment exactly when π₀ = 1).
3. Refine `refEcology` choice: edge-mass argmax is good for fixed
   tree, but during MCMC the tree changes. Should refEcology be
   re-chosen each iteration? Probably no — keeping it fixed avoids
   label switching, and the choice is mathematically arbitrary.
4. Sub-agents are silently failing in this session (0-byte output
   files from `general-purpose` Opus agents). Implementation needs
   either a different agent type, manual execution, or different
   permission setup.
