# Mk' Relabelling Correction: Proof of the Correct Formula

## Summary

The relabelling correction in the Mk' (MkPrime) model is **wrong** in both
the MkPrime R package and in RevBayes (`PhyloCTMCSiteHomogeneousMkPrime`).
The error causes the correction to **penalise** higher k' when it should
**favour** it, pushing k' toward kObs and effectively disabling hidden-state
inference.

**Wrong formula** (currently in both codebases):

```
log C = log(kObs!) − log(k'!) + log((k' − kObs)!) + kObs × log(k')
      = log( k'^kObs / C(k', kObs) )          [C = binomial coefficient]
```

This *decreases* with k'.

**Correct formula:**

```
log C = log(k'!) − log((k' − kObs)!)
      = log P(k', kObs)                        [P = falling factorial]
```

This *increases* with k'.

The error is a factor of `k'^kObs / (kObs! × C(k', kObs)²)` in the
likelihood, or equivalently `2 × log C(k', kObs) − kObs × log(k') +
log(kObs!)` nats of log-likelihood error per character per MCMC step.  For
the k' = 2 → 3 transition with kObs = 2, the per-character error is exactly
**log 4 ≈ 1.386 nats**.

---

## Setup

Consider a single morphological character observed at *n* tips of a
phylogeny *T*.  The character has **kObs** distinct states in the observed
data.  Under the Mk' model, the true number of states is **k' ≥ kObs**;
the remaining *u = k' − kObs* states are unobserved (no extant tip occupies
them).

The substitution model is Jukes–Cantor with k' states: all off-diagonal
rates equal, normalised so that Q_ii = −1.  Root frequencies are uniform:
π_s = 1/k' for s = 0, …, k' − 1.

---

## The Felsenstein likelihood computes one specific labelling

The standard Felsenstein pruning algorithm computes:

```
L_fels(data | T, k') = Σ_{s=0}^{k'-1}  (1/k') × CL(root, s)
```

where CL(root, s) is the conditional likelihood at the root given root
state *s*, computed by the usual post-order traversal.

Crucially, this computation assumes a **specific assignment** of the
observed data labels {0, 1, …, kObs − 1} to model states {0, 1, …,
kObs − 1}.  That is: observed label *i* ↔ model state *i*, and model
states kObs, …, k' − 1 are unobserved (no tip has partial likelihood 1
for those states).

---

## Counting equivalent labellings

Under JC(k'), all k' states are exchangeable.  The observed labels could
equally well have been assigned to **any** kObs of the k' model states, in
**any** order.  The number of such injective maps is the **falling
factorial**:

```
P(k', kObs) = k'! / (k' − kObs)!
            = C(k', kObs) × kObs!
```

where C(k', kObs) is the binomial coefficient (choosing *which* kObs states
are observed) and kObs! is the number of label permutations.

By JC symmetry, every such map produces the **same** Felsenstein
likelihood.  Therefore the marginal probability of the observed data,
summed over all equivalent labellings, is:

```
P(data | T, k') = P(k', kObs) × L_fels(data | T, k')
```

The log relabelling correction to be **added** to log L_fels is:

```
log P(k', kObs) = log(k'!) − log((k' − kObs)!)
```

### Value at k' = kObs

When k' = kObs, the correction is log(kObs!).  This is a constant
independent of k' and cancels in all MH ratios, but it is needed for
correct absolute log-posterior values (e.g. for stepping-stone marginal
likelihoods).

### Behaviour

| k' | kObs | P(k', kObs) | log P(k', kObs) |
|----|------|-------------|-----------------|
| 2  | 2    | 2           | 0.693           |
| 3  | 2    | 6           | 1.792           |
| 4  | 2    | 12          | 2.485           |
| 5  | 2    | 20          | 2.996           |
| 3  | 3    | 6           | 1.792           |
| 4  | 3    | 24          | 3.178           |
| 5  | 3    | 60          | 4.094           |

The correction **increases** with k': more model states means more ways
the observed labels could have been assigned.

---

## Proof by enumeration

For a star tree with 3 tips, branch length *t* = 0.3, and data pattern
(0, 0, 1) (kObs = 2):

We enumerate **all** injective maps from {0, 1} → {0, …, k' − 1} and,
for each map, compute the Felsenstein likelihood of the resulting
model-state tip pattern.  All maps produce the same L by JC symmetry.

The total (summed) probability equals P(k', 2) × L_fels exactly:

```r
# Direct enumeration (all injective maps from kObs labels to k model states)
direct_enum_L <- function(k, obs_pattern, t) {
  total_L <- 0;  n_maps <- 0
  for (a in 0:(k-1))
    for (b in 0:(k-1)) {
      if (a == b) next
      model_tips <- ifelse(obs_pattern == 0, a, b)
      total_L <- total_L + brute_force_L(k, model_tips, t)
      n_maps <- n_maps + 1
    }
  list(total = total_L, n_maps = n_maps)
}

# Results:
#   k=2: 2 maps,  total P = 0.17470145   = 2  × 0.08735072  ✓
#   k=3: 6 maps,  total P = 0.16461328   = 6  × 0.02743555  ✓
#   k=4: 12 maps, total P = 0.15880135   = 12 × 0.01323345  ✓
#   k=5: 20 maps, total P = 0.15532411   = 20 × 0.00776621  ✓
#   k=6: 30 maps, total P = 0.15303815   = 30 × 0.00510127  ✓
```

---

## What the current code computes

### MkPrime R package (`src/corrections.cpp`)

```cpp
double log_correction =
    std::lgamma(kObs + 1.0) -       // log(kObs!)
    std::lgamma(kPrime + 1.0) +     // −log(k'!)
    std::lgamma(kPrime - kObs + 1.0) + // log((k'−kObs)!)
    kObs * std::log(kPrime);         // kObs × log(k')
```

This computes `log(kObs!) − log(k'!) + log((k'−kObs)!) + kObs × log(k')`
= `log( k'^kObs / C(k', kObs) )`.

### RevBayes (`PhyloCTMCSiteHomogeneousMkPrime.h`, line 268)

```cpp
return logFactorial_[kObs] -
       static_cast<double>(kObs) * std::log(static_cast<double>(kObs)) -
       logFactorial_[k] +
       logFactorial_[k - kObs] +
       static_cast<double>(kObs) * std::log(static_cast<double>(k));
```

This computes the same thing with an extra `−kObs × log(kObs)` term that
normalises C to 0 at k' = kObs.  **The MH-ratio delta is identical in both
implementations.**

### Comparison table (kObs = 2)

| k' | Correct log P(k',2) | Code delta | Correct delta | Error/char |
|----|---------------------|------------|---------------|------------|
| 2  | 0.693               | 0.000      | 0.000         | 0.000      |
| 3  | 1.792               | −0.288     | +1.099        | **1.386**  |
| 4  | 2.485               | −0.406     | +1.792        | **2.197**  |
| 5  | 2.996               | −0.470     | +2.303        | **2.773**  |

The code penalises higher k' when it should favour it.  The error grows
with k'.

---

## Impact

The incorrect correction adds ~1.4 nats of bias **per character** against
each unit increase in k'.  In a dataset with *n* transformational
characters, the bias against moving even a single character from k' = kObs
to k' = kObs + 1 is ~1.4 × n nats (overwhelmingly against).  This
effectively prevents the model from ever inferring hidden states, making
k' ≡ kObs a near-certain outcome regardless of the data.

Observed behaviour: in the MkPrime R package, 44 of 54 transformational
characters **never** leave k' = kObs in the posterior, and the remaining 10
leave only ~0.5% of the time.  The hyperparameter *p* converges to ~0.98,
reinforcing the lock-in.

---

## The fix

Replace the formula in both codebases with:

```cpp
// Correct: log of the falling factorial P(k', kObs) = k'! / (k'−kObs)!
return std::lgamma(kPrime + 1.0) - std::lgamma(kPrime - kObs + 1.0);
```

### RevBayes patch

In `PhyloCTMCSiteHomogeneousMkPrime.h`, replace
`computeLogRelabellingCorrection`:

```cpp
inline double RevBayesCore::PhyloCTMCSiteHomogeneousMkPrime
    ::computeLogRelabellingCorrection(size_t k, size_t kObs) const {
    if (k == 0 || k > kMax_ || kObs > kMax_ || kObs > k) {
        return RbConstants::Double::neginf;
    }
    // Falling factorial: P(k, kObs) = k! / (k − kObs)!
    return logFactorial_[k] - logFactorial_[k - kObs];
}
```

---

## Test strategy

1. **Unit test**: verify `mk_prime_relabel_log(k, kObs)` returns
   `lfactorial(k) − lfactorial(k − kObs)` for a range of k and kObs
   values.

2. **Brute-force enumeration test**: for a small star tree, verify that
   `P(k', kObs) × L_fels` equals the sum over all injective maps.

3. **Posterior test**: run a short MCMC on a dataset where hidden states
   are expected (e.g. short tree length) and verify that k' > kObs appears
   with nontrivial posterior probability.

4. **RevBayes cross-validation**: compare log-likelihoods from the fixed
   MkPrime R package with RevBayes (after patching) for a reference dataset
   at known parameter values.

---

## Origin of the error

The incorrect formula `k'^kObs / C(k', kObs)` likely arose from a
derivation that attempted to combine the binomial coefficient with a
`1/k'^kObs` term (perhaps from the root frequencies), but the root
frequencies are already included in the Felsenstein likelihood.  The
comment in `corrections.cpp` even flags the confusion: *"Wait — let's be
precise."*  The result double-counts (and inverts) the root-frequency
contribution.
