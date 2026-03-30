# MCMC Optimizations from MrBayes, BEAST, and RevBayes

**Agent D research report — 2026-03-29**

Surveyed: MrBayes (NBISweden/MrBayes `develop`), BEAST X (beast-dev/beast-mcmc `main`),
RevBayes (revbayes/revbayes `development`), and published literature on each.

---

## Executive Summary

MkPrime already has a strong MCMC engine — Gibbs SPR/swap with partial CL,
weighted branch moves, adaptive move weights via bandit tuning, parallel tempering
with adaptive temperatures, and slice sampling. Most of the "big wins" from these
codebases are already implemented or close analogues exist. But there are several
medium-value ideas worth considering:

| Priority | Idea | Source | Estimated gain | Effort |
|----------|------|--------|----------------|--------|
| **HIGH** | Bactrian proposal kernels for scalar moves | Yang & Rodríguez 2013 / BEAST2 ORC | ~15–50% ESS/iter for scalars | Low |
| **HIGH** | Parsimony-guided SPR weighting | MrBayes 3.2 (Ronquist et al. 2020) | Up to 10× convergence speedup | Medium |
| **MEDIUM** | AVMVN (Adaptive Variance Multivariate Normal) for correlated scalars | BEAST (Baele et al. 2017) | Up to 14× for correlated params | Medium |
| **MEDIUM** | Conditional likelihood node-level dirty flags | MrBayes | Avoids full tree LL for parameter-only moves | Low-Medium |
| **MEDIUM** | Extending SPR/TBR with geometric stopping | MrBayes (Lakner et al. 2008) | Better than uniform SPR, worse than Gibbs SPR | Already have analogues |
| **LOW** | SubTreeLeap operator | BEAST | Replaces NNI+SPR+node-slider as single move | Medium |
| **LOW** | SSE/SIMD vectorization of inner pruning loop | MrBayes 3.2 | ~2× raw speed for likelihood | High |
| **LOW** | BEAGLE library integration | MrBayes + BEAST | GPU acceleration for large datasets | Very High |

---

## 1. Bactrian Proposal Kernels (HIGH priority, LOW effort)

### What they are

The standard Gaussian random walk proposal q(x'|x) = N(x, σ²) spends too much time
proposing values very close to x (small steps → high autocorrelation) or very far from x
(rejected → wasted). The **Bactrian kernel** (Yang & Rodríguez 2013, PNAS 110:19307)
uses a mixture of two Gaussians displaced from the centre:

```
q(Δ | σ, m) = ½ N(+m, σ²(1-m²)) + ½ N(-m, σ²(1-m²))
```

with m ≈ 0.95 by default. This has a bimodal shape (like a Bactrian camel's humps) that
suppresses near-zero steps while keeping the same variance.

### Evidence

- Yang & Rodríguez showed ≥50% efficiency gain over Gaussian for all target distributions tested.
- Douglas et al. (2021, PLOS Comp Bio) confirmed 15–20% ESS improvement in phylogenetic
  relaxed clock context with essentially zero computational overhead.
- Now standard in BEAST2's ORC package and recommended by the developers.

### What MkPrime should do

Replace the uniform perturbation `U - 0.5` in all scale proposals (`exp(λ × (U - 0.5))`)
with a Bactrian kernel. The implementation is trivial:

```cpp
// Bactrian kernel (m = 0.95 by default)
double bactrian(RNG& rng, double m = 0.95) {
    double s = sqrt(1.0 - m * m);
    double z = R::rnorm(0.0, s);
    return (rng.uniform01() < 0.5) ? z + m : z - m;
}
// Then in scale proposal:
// newVal = oldVal * exp(lambda * bactrian(rng));
```

This applies to: `tree_length`, `rate_loss`, `rate_neo`, `rate_log_sd`, `beta_scale`,
and the `beta_simplex` concentration parameter.

---

## 2. Parsimony-Guided SPR Proposals (HIGH priority, MEDIUM effort)

### What they are

MrBayes 3.2 introduced parsimony-biased SPR (pSPR) and TBR (pTBR) proposals
(Ronquist et al. 2020, Syst Biol 69:1016). Instead of choosing a regraft position
uniformly at random (standard SPR) or with a geometric stopping probability (extending SPR),
the pSPR computes the **parsimony score** of each candidate tree and weights the
proposal probability proportionally to exp(-β × parsimony_score).

### How it works

1. Prune subtree at random edge (same as standard SPR)
2. Evaluate parsimony score at each candidate regraft position (O(n) via Fitch algorithm)
3. Weight candidates by exp(-β × score), normalize to probabilities
4. Sample regraft position from this distribution
5. Correct Hastings ratio for the asymmetric proposal

### Evidence

- Ronquist et al. 2020 showed pSPR can improve convergence speed by **an order of magnitude**
  on some empirical datasets.
- The parsimony calculation is O(n_char × n_tips) per candidate — much cheaper than
  a full likelihood evaluation. For morphological data with ~100 characters, this is trivial.

### Relevance to MkPrime

MkPrime already has **Gibbs SPR** which weights candidates by the **full likelihood**,
not just parsimony. This is strictly better but much more expensive per proposal.
Parsimony-guided SPR sits in a sweet spot:

- Cheaper than Gibbs SPR (no likelihood eval, just Fitch parsimony)
- Better guided than extending/geometric SPR
- Could serve as a **complement** to Gibbs SPR, especially in early warmup when
  likelihood is expensive but parsimony already guides well.

**Recommendation**: Consider adding pSPR as an option. During warmup (when we want
cheap, exploratory moves), pSPR could replace some of the NNI weight. The Fitch
downpass for Mk data is extremely cheap (binary characters → bitwise OR).

---

## 3. Adaptive Variance Multivariate Normal (AVMVN) Kernel (MEDIUM priority)

### What it is

BEAST implements an **adaptive multivariate normal** proposal for all continuous
parameters simultaneously (Baele et al. 2017, Bioinformatics 33:1798). Instead of
updating parameters one at a time, the AVMVN learns the posterior covariance matrix
during MCMC and proposes joint updates along correlated directions.

### Evidence

- Up to 14× ESS improvement for highly correlated parameters.
- Most effective when many continuous parameters have strong posterior correlations.

### Relevance to MkPrime

MkPrime has relatively few correlated scalar parameters (tree_length, rate_loss,
rate_neo, rate_log_sd — maybe 4-5 scalars). The AVMVN payoff scales with the
number of correlated parameters, so the benefit may be modest. But:

- The `rate_loss` × `rate_neo` correlation is already handled by a joint move.
- `tree_length` × `rate_log_sd` may also be correlated.

**Recommendation**: Low-hanging fruit would be a 2D Bactrian joint proposal for the
most correlated pairs, before investing in full AVMVN machinery. Monitor correlations
in early runs to see if this warrants a full multivariate operator.

---

## 4. Conditional Likelihood Dirty Flags (MEDIUM priority, LOW effort)

### What MrBayes does

MrBayes uses per-node `upDateCl` and `upDateTi` flags. When a scalar parameter
changes (e.g., substitution model parameter), only the nodes whose transition
probability matrices (`Ti`) need recomputing get flagged. The CL recomputation
propagates only from flagged nodes up to the root, skipping unchanged subtrees.

For parameter-only moves (no topology change), this means:
- Changing a rate parameter: flag all `Ti` that use that rate, recompute CLs from tips to root
- Changing one branch length: flag only that node's `Ti`, recompute CLs from that node up

### What MkPrime currently does

MkPrime has **partition-level caching** (M-064): when a move only affects neomorphic
partitions (e.g., `rate_loss`), only those partition likelihoods are recomputed. This is
the right level of granularity for the current model.

However, for **topology moves that don't change all CLs** (like NNI affecting only 2 nodes),
the full tree is currently recomputed for all partitions. Node-level dirty tracking could
save significant work: an NNI only invalidates CLs along the path from the modified nodes
to the root (O(depth) instead of O(n_nodes)).

### Recommendation

This is essentially what the partial CL framework (M-105, M-111) already does for Gibbs
moves. Extending it to standard NNI/SPR would provide the same benefit. This is worth
doing but requires careful engineering of the CL rollback on rejection.

---

## 5. Extending SPR with Geometric Stopping (Already implemented — MkPrime has analogues)

### What it is

MrBayes's `Move_ExtSPR` and `Move_ExtTBR` move the regraft point away from the prune
point with a geometric stopping probability (extension probability ~0.5). This makes
larger rearrangements less likely, keeping proposals in a "neighborhood" of the current tree.

### Relevance

MkPrime's TBR (M-053) already includes subtree re-rooting on the pruned side, and the
Gibbs SPR evaluates all candidates weighted by likelihood. The extending SPR is strictly
inferior to Gibbs SPR for our use case. **No action needed.**

---

## 6. Node Slider Move (Low priority, but notable omission)

### What it is

MrBayes has `Move_NodeSlider` which slides an internal node along its branch, changing
the relative split of branch length between parent and children. This is different from
the Beta Simplex move (which redistributes mass between two randomly chosen branches).

### Relevance

For unrooted trees with our tree_length × relative_branch_lengths parameterisation,
the Beta Simplex already serves this role. **No action needed.**

---

## 7. SSE/SIMD Vectorization (LOW priority for morphological data)

### What MrBayes does

MrBayes 3.2 uses SSE (Streaming SIMD Extensions) intrinsics for the inner pruning loop,
processing 2 or 4 states simultaneously. For DNA (k=4), this gives a ~2× speedup.

### Relevance

MkPrime uses the JC symmetry optimization (OPP-1) which already reduces the inner loop
from O(k²) to O(k). For morphological data with k=2-10 states, SIMD vectorization would
have minimal benefit because the inner loops are already tiny. The bottleneck is the
O(n_chars × n_nodes) traversal, not the O(k) per-site computation.

**Recommendation**: Not worth pursuing for morphological data. Focus optimization
effort on reducing the number of full tree traversals instead.

---

## 8. BEAGLE Library Integration (LOW priority)

### What it is

BEAGLE is a shared library for GPU-accelerated phylogenetic likelihood computation,
used by both MrBayes and BEAST. It can provide 50× speedup for codon models and
~10× for DNA.

### Relevance

For morphological data with k ≤ 10 states and typically < 500 characters, the
likelihood computation is already fast (~ms per evaluation). GPU overhead would
likely dominate. **Not worth pursuing.**

---

## 9. Adaptive Operator Weighting During Sampling (Already implemented)

Both BEAST2 (Douglas et al. 2021) and MkPrime implement adaptive operator weighting.
MkPrime's bandit tuning (Phase 2) is arguably more sophisticated than BEAST2's approach.
**No action needed.**

---

## 10. Adaptive Parallel Tempering (Already implemented)

BEAST2's CoupledMCMC package (Müller & Bouckaert 2020) adapts temperature spacing
to target ~23.4% swap acceptance. MkPrime already adapts the heat parameter during
warmup to target 25% swap acceptance. **No action needed.**

---

## 11. Other Notable Techniques from the Literature

### Stepping Stone Sampling (Model Comparison)
Both MrBayes 3.2 and BEAST implement stepping-stone / path-sampling for marginal
likelihood estimation (Bayes factors). This is for model comparison, not for improving
MCMC mixing per se. Not relevant to the current speed/ESS goal.

### Thermodynamic Integration
Similar to stepping stone; used for Bayes factors. Not relevant here.

### Stochastic Variational Inference (SVI)
Dang & Kishino 2019 explored SVI for phylogenetics. Still experimental and not
adopted by any major package. Not recommended.

---

## Prioritized Recommendations

### Quick Wins (implement now)
1. **Bactrian kernels** for all scale proposals — trivial to implement, guaranteed ~15-50% ESS improvement per scalar parameter iteration.

### Medium-term (next phase)
2. **Parsimony-guided SPR** — add as a complement to Gibbs SPR, especially useful during warmup. The Fitch downpass for Mk data is essentially free (bitwise operations on {0,1,...,k-1} state sets).

### Worth investigating
3. **2D joint Bactrian proposals** for correlated parameter pairs (tree_length × rate_log_sd, etc.)
4. **Node-level CL dirty flags** for standard NNI/SPR (extending the partial CL framework already used by Gibbs moves)

### Skip
- SSE/SIMD for morphological data (k too small)
- BEAGLE integration (data too small)
- Extending SPR (Gibbs SPR already better)
- Full AVMVN (too few correlated scalars to justify)

---

## Key References

- Yang Z, Rodríguez CE (2013). Searching for efficient MCMC proposal kernels. PNAS 110:19307-12.
- Ronquist F, Kudlicka J, et al. (2020). Using parsimony-guided tree proposals to accelerate convergence in Bayesian phylogenetic inference. Syst Biol 69:1016-32.
- Baele G, Lemey P, Rambaut A, Suchard MA (2017). Adaptive MCMC in Bayesian phylogenetics. Bioinformatics 33:1798-805.
- Douglas J, Zhang R, Bouckaert R (2021). Adaptive dating and fast proposals: revisiting the phylogenetic relaxed clock model. PLOS Comp Biol 17:e1008322.
- Müller NF, Bouckaert R (2020). Adaptive Metropolis-coupled MCMC for BEAST 2. PeerJ 8:e9473.
- Lakner C, et al. (2008). Efficiency of MCMC tree proposals in Bayesian phylogenetics. Syst Biol 57:86-103.
- Höhna S, Drummond AJ (2012). Guided tree topology proposals for Bayesian phylogenetic inference. Syst Biol 61:1-11.
