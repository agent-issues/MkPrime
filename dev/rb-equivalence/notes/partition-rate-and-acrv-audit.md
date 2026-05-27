# Partition-rate and ACRV discretisation audit: MkPrime vs RevBayes

**Status:** Mathematical audit, no code changes. Follow-up to
`dev/rb-equivalence/notes/cross-sampler-rhat-investigation.md`; supersedes that
note's §2 (ACRV) framing, which understated the discrepancy.

**TL;DR.** Both implementations have one defensible parameterisation choice and
one indefensible one. RB gets the **partition-rate** convention right
(nChar-weighted mean rate = 1, edge lengths interpretable across Mk and
Mk\_n; user's design intent is therefore correctly realised by RB).
MkPrime gets the **ACRV discretisation** right (μ on log-scale = −σ²/2 *and*
explicit renormalisation to mean 1; this is Harrison & Larsson 2014 eq. 2, the
canonical published convention). Each side should adopt the other's
parameterisation. The user's worry "do we need n\_neo / n\_trans weighting?"
is answered **yes** — and that weighting is already correctly applied in RB,
not in MkPrime.

---

## Issue 1 — Partition rate parameterisation

### Notation

Let $n_1 = n_{\text{neo}}$, $n_2 = n_{\text{trans}}$, $n = n_1 + n_2$, and let
the partition rate vector be $\rho = (\rho_1, \rho_2)$. Let $r := r_{\text{neo}}$,
$T :=$ `tree_length`, $b_e :=$ `rel_br_lengths[e]` (with $\sum_e b_e = 1$).
The per-edge expected substitutions for a character in partition $i$ is

$$\lambda_{i,e} \;=\; T \cdot b_e \cdot \rho_i.$$

### What the user means by "edge lengths comparable between Mk and Mk\_n"

Of the three candidate definitions in the brief, candidate (C) — derived from
the user's stated semantics for `tree_length` — is the only one that survives
scrutiny.

- **(A) Marginal posterior on `tree_length` is invariant when neo data is
  removed.** This is too strong: dropping a partition changes the likelihood,
  which changes the posterior on every parameter that the partition was
  informing. No proposal scaling can make posteriors identical across datasets.

- **(B) `tree_length` is "expected substitutions per character on the average
  edge".** This is the right *interpretation* but requires a precise
  mathematical statement to be checkable.

- **(C, defensible).** Adopt (B) literally: define `tree_length` such that
  *averaging over the empirical character distribution* (i.e. weighting by
  partition size), the expected number of substitutions per character on
  edge $e$ equals $T \cdot b_e$. Formally,

  $$\frac{1}{n} \sum_{i=1}^{2} n_i \cdot \lambda_{i,e} \;=\; T \cdot b_e
    \quad \Longleftrightarrow \quad \boxed{\;\frac{1}{n}\sum_i n_i \rho_i = 1\;}$$

  This is the unique condition that gives `tree_length` a sampler-independent
  meaning: the total number of substitutions summed over the dataset and over
  the tree is $T \cdot \sum_i n_i \rho_i = T \cdot n$, regardless of
  $r_{\text{neo}}$. Drop one partition and the same definition applies to the
  remaining one (with weight $n_i / n_i = 1$). Hence "edge lengths comparable
  between Mk and Mk\_n" means: the constraint that the data-weighted mean
  partition rate equals 1.

### Checking each implementation against (C)

**RB** (`dev/rb-equivalence/templates/by_nt_9v.template.Rev:47`):
```rev
partition_rate := [ rate_neo / (1 + rate_neo), 1 / (1 + rate_neo) ] / nChar * sum(nChar)
```
Elementwise, with `nChar = (n_1, n_2)` and `sum(nChar) = n`:

$$\rho_1 = \frac{r}{1+r} \cdot \frac{n}{n_1}, \qquad
  \rho_2 = \frac{1}{1+r} \cdot \frac{n}{n_2}.$$

Weighted mean:

$$\frac{1}{n}\sum_i n_i \rho_i
  \;=\; \frac{1}{n}\!\left( n_1 \cdot \frac{r}{1+r} \cdot \frac{n}{n_1}
                          + n_2 \cdot \frac{1}{1+r} \cdot \frac{n}{n_2}\right)
  \;=\; \frac{r}{1+r} + \frac{1}{1+r} \;=\; 1. \quad\checkmark$$

So RB satisfies (C) identically, for any $r > 0$, any $n_1, n_2 > 0$.
Algebraically verified for $n_1=5, n_2=12, r=2$: $\rho = (5.111, 0.383)$,
weighted mean $= 1$.

**MkPrime** (`src/mcmc_likelihood.cpp:2371-2372` for the neo path;
`src/mcmc_likelihood.cpp:2464,2474,2499,2511,2515,2520,2522` for the trans
paths which all pass `edgeLen` unscaled):

$$\rho_1 = r, \qquad \rho_2 = 1.$$

Weighted mean:

$$\frac{n_1 r + n_2}{n}.$$

This equals 1 only when $r = 1$. For the smoke-run dataset
(`out/mkprime_635_by_nt_9v.rds`, $n_1 = 5$ neo, $n_2 = 12$ trans), the
posterior median of $r_{\text{neo}}$ is **0.204**, giving a weighted mean of
$(5 \cdot 0.204 + 12)/17 = 0.766$. So in MkPrime's parameterisation,
`tree_length` is implicitly inflated by $1/0.766 \approx 1.31$ relative to the
"expected substitutions per character" interpretation it has under RB —
mediated through the $(T, r_{\text{neo}})$ ridge.

### The user's worry: is n\_neo / n\_trans weighting needed?

**Yes.** Without the $n/n_i$ weighting, the rate scalars are dimensionally
meaningful but **not** dataset-invariant: e.g. RB's
$\rho_1 / \rho_2 = r \cdot n_2 / n_1$, so adding more trans characters (with
fixed $r$) increases the neo:trans rate ratio applied to edges. This is the
right behaviour, because what's being preserved is the *aggregate* expected
number of substitutions, and aggregating over a partition with more characters
puts more characters at the lower rate, requiring the partition's per-character
rate to compensate to hold the total constant.

Equivalent, more transparent formulation: parameterise the "raw" rate of
partition $i$ as $w_i$ (in MkPrime: $w_1 = r$, $w_2 = 1$), then set

$$\rho_i \;=\; \frac{w_i}{\bar w}, \qquad \bar w \;:=\; \frac{1}{n}\sum_j n_j w_j.$$

This automatically enforces $\frac{1}{n}\sum_i n_i \rho_i = 1$ for **any**
choice of raw rates and any number of partitions, and generalises RB's formula
(which is the special case $w = (r, 1)$).

The user's *worry* that the weighting might already be present and shouldn't
be added in MkPrime is unfounded: the weighting is *not* present in MkPrime,
which is the cause of the (`tree_length`, `rate_neo`) drift signature.
Without the weighting, $T$ and $r$ are jointly under-identified in a
data-asymmetry-dependent way, so RB and MkPrime see slightly different
posterior support on the marginal of each.

### Verdict (Issue 1)

**RB is correct; MkPrime should be patched.** The user's design intent —
that edge lengths be directly comparable between Mk and Mk\_n — is realised
by RB and not by MkPrime.

### Recommended canonical form (Issue 1)

For $K$ partitions with character counts $n_1, \dots, n_K$ and unnormalised
relative rates $w_1, \dots, w_K$ (e.g. $w = (r_{\text{neo}}, 1)$ for the 2-
partition case),

$$\boxed{\;\rho_i \;=\; \frac{w_i \cdot \sum_j n_j}{\sum_j n_j w_j}
                  \;=\; \frac{w_i}{\bar w}, \quad
       \bar w := \frac{1}{n} \sum_j n_j w_j.\;}$$

Identity: $\frac{1}{n}\sum_i n_i \rho_i = 1$ for any $w$.

For $K=2$ partitions with $w = (r, 1)$, equivalently
$\rho_1 = r \cdot n / (n_1 r + n_2)$ and $\rho_2 = n / (n_1 r + n_2)$. RB's
template uses a different (but algebraically equivalent) factorisation
$\rho_i = (w_i / \sum w) \cdot (n / n_i)$, which yields identical $\rho_i$:

$$\rho_1 = \frac{r}{1+r} \cdot \frac{n}{n_1}, \qquad
  \rho_2 = \frac{1}{1+r} \cdot \frac{n}{n_2}.$$

Substituting $n_1 r + n_2$ for the denominator in the first form recovers RB's
expression after multiplying numerator and denominator by $1/(1+r)$.

### Patch suggestion (Issue 1)

Apply to MkPrime, *not* RB. In `src/mcmc_likelihood.cpp`:

- Compute `wbar = (n_neo * rate_neo + n_trans) / n_total` once at the top of
  the partition loop (or pass it in alongside `rateNeo`).
- Replace `neoEl[i] = edgeLen[i] * rateNeo` (line 2372) with
  `neoEl[i] = edgeLen[i] * rateNeo / wbar`.
- Construct a parallel `transEl[i] = edgeLen[i] / wbar` and pass it to *all*
  trans-partition (type==2) pruning sites: lines 2464, 2474, 2499, 2504, 2511,
  2515, 2520, 2522, and 2525 (the `constant_site_prob_jc` call). Also the
  `type==0`/transformational branch at ≥2536.

This is a multi-site edit (≈15 call-sites accept `edgeLen`); not a trivial
one-line patch. No patch attached.

Caveat (Issue 1, item 1): the fix changes the *meaning* of `tree_length` on
the MkPrime side. Existing posterior samples are not directly comparable to
post-patch samples for `tree_length` and `rate_neo` (a deterministic
rescaling, but not the identity). Lane MCMC-diagnostician should re-validate
SBC and posterior-predictive coverage after the patch.

Caveat (Issue 1, item 2): there is no per-`kPrime` group weighting analogous
to the trans/neo split. If transformational characters in the future are
split into multiple Q-matrices (e.g. by `kPrime`), the same weighting formula
extends naturally with $K > 2$ but must be added explicitly. Should be picked
up by **math-prover** lane during any future model-splitting work.

---

## Issue 2 — ACRV (across-character rate variation) discretisation

### Two independent discrepancies, not one

The prior investigation note (line 7 and §2) framed this as a
renormalisation-only difference, claiming "both samplers use
`qnorm((i+0.5)/nCat, μ=−σ²/2, σ=rate_log_sd)`, but MkPrime then renormalises".
That description is **incorrect for RB**. There are two independent
discrepancies, and they partially cancel:

1. **Log-scale mean $\mu$.** RB's template (`by_nt_9v.template.Rev:41`) uses
   `dnLognormal( 0, rate_log_sd )`. From the RevBayes documentation and
   `LognormalDistribution.cpp`, the first argument is the log-scale mean.
   So RB samples categories from a $\text{LogNormal}(\mu=0,\,\sigma)$, which
   has *true* mean $e^{\sigma^2/2}$ — not 1. MkPrime uses $\mu = -\sigma^2/2$
   (`src/mcmc_likelihood.cpp:78`, `R/acrv.R:21`), so its underlying continuous
   distribution has true mean 1.

2. **Discrete renormalisation.** MkPrime divides by the empirical mean of the
   six bin medians (`R/acrv.R:39`, `src/mcmc_likelihood.cpp:85`). RB's
   `fnDiscretizeDistribution` does not
   (`src/core/functions/distribution/DiscretizeDistributionFunction.cpp:99-108`).

The combined effect: at the same $\sigma$, MkPrime's six discrete rates have
arithmetic mean exactly 1 (renormalised); RB's have arithmetic mean
$\approx e^{\sigma^2/2}$ minus the finite-bin "quantile median ≤ continuous
mean" Jensen-like correction.

### Numerical mean of the discretised lognormal

Computed with `qnorm((i-0.5)/6, μ, σ)`, $i = 1, \dots, 6$, for $n_{\text{Cat}} = 6$:

| $\sigma$ | RB ($\mu = 0$): mean | MkPrime ($\mu = -\sigma^2/2$): mean (pre-renorm) | MkPrime mean (post-renorm) |
|---:|---:|---:|---:|
| 0.10 | 1.0040 | 0.9990 | **1.0000** |
| 0.30 | 1.0413 | 0.9910 | **1.0000** |
| 0.50 | 1.1039 | 0.9742 | **1.0000** |
| 1.00 | 1.4591 | 0.8850 | **1.0000** |
| 1.50 | 1.9869 | 0.7197 | **1.0000** |
| 2.00 | 3.7084 | 0.5019 | **1.0000** |
| 3.00 | 30.34 | 0.1361 | **1.0000** |

(R numerical verification:
`exp(qnorm((seq_len(6)-0.5)/6, μ, σ))`, see Bash invocations in this audit's
working session.)

At $\sigma = 1$, RB applies effective rates with mean **1.459**, MkPrime
applies effective rates with mean **1.000**. RB's posterior on `tree_length`
must therefore deflate by a factor of ≈1.46 relative to MkPrime's to match the
likelihood — a much larger shift than the prior note's "10–15%". The prior
note's mistake was an algebra/notation slip on the value of $\mu$ used in
`dnLognormal`.

### The smoke-run posterior context

Smoke run pid 635, by\_nt\_9v ($n_1 = 5$ neo, $n_2 = 12$ trans, 19 taxa):
posterior median $\sigma = 0.55$ (MkPrime side, run 1 of 2). At that $\sigma$:

- RB's effective mean rate: ≈ 1.13 (from interpolation of the table above).
- MkPrime's: 1.00.

A 13% inflation of `tree_length` is consistent with the observed cross-sampler
rhat of 1.037 on `tree_length`, especially when partially absorbed by the
$(T, r_{\text{neo}})$ ridge from Issue 1. The two issues compound: Issue 1
drives the rate-neo drift directly, and Issue 2 drives an additional ~13%
tree-length inflation on the RB side. Both posteriors of `rate_log_sd` also
shift, because at smaller $\sigma$ the RB-side mean-rate inflation shrinks
(and so the posterior $\sigma$ is mildly pulled downwards in RB to compensate
for the over-multiplication of tree length). This is consistent with the
observed `rate_log_sd` rhat of 1.027.

### Mathematical check: what is the analytic limit?

Let $r_i = e^{\mu + \sigma \Phi^{-1}(p_i)}$ with $p_i = (i - 1/2) / n$. Then

$$\bar r \;=\; \frac{1}{n}\sum_{i=1}^n r_i.$$

As $n \to \infty$, $\bar r \to E[\text{LogNormal}(\mu, \sigma)] = e^{\mu + \sigma^2/2}$
(this is a Riemann-sum / Monte-Carlo argument: the $p_i$ are exactly the
midpoint quadrature nodes of the unit interval, and the integrand
$e^{\mu + \sigma \Phi^{-1}(p)}$ has integral $e^{\mu + \sigma^2/2}$).

For finite $n$, the midpoint quadrature underestimates this integral because
the integrand is convex in $\Phi^{-1}(p)$ near $p = 1$ (heavy upper tail of
the lognormal). Hence:

- With $\mu = -\sigma^2/2$ (MkPrime): $\bar r \to 1$ as $n \to \infty$, but
  $\bar r < 1$ for finite $n$. The deficit grows monotonically with $\sigma$
  (≈ 1% at σ=0.3, ≈ 12% at σ=1, ≈ 50% at σ=2).
- With $\mu = 0$ (RB): $\bar r \to e^{\sigma^2/2}$ as $n \to \infty$.

In neither case does $\bar r = 1$ hold exactly for finite $n$ without explicit
renormalisation. **The user's recollection that "the maths should check out
to a normalisation on 1" is correct only with explicit Σ-rescaling.**

### Published canon

Harrison & Larsson 2014 (Syst Biol, doi:10.1093/sysbio/syu098,
"Among-Character Rate Variation Distributions in Phylogenetic Analysis of
Discrete Morphological Characters") explicitly prescribes:

$$r_i \;=\; \frac{e^{\mu + \sigma \Phi^{-1}(p_i)}}{(1/n)\sum_j e^{\mu + \sigma \Phi^{-1}(p_j)}}$$

with $\mu = -\sigma^2/2$ (their eq. 2; verified via WebFetch on Oxford
Academic abstract+methods, which reproduces the formula and states "rates
[were] rescaled so that the mean of the discrete distribution is 1.0").
This is the formula MkPrime implements; it is the canonical method.

Yang 1994 (J Mol Evol, doi:10.1007/BF00160154), the foundational reference
for discretised rate categories, also normalises: Yang's discrete-gamma method
prescribes that the discrete category rates have arithmetic mean 1 by
construction, and provides explicit normalising factors for the mean-based
and median-based discretisations. RB's `fnDiscretizeGamma` follows Yang
(rescales by `factor = a / b * nCats`; see
`src/core/functions/distribution/DiscretizeGammaFunction.cpp`); but RB's
**generic** `fnDiscretizeDistribution` does *not* rescale.

Wagner 2012 (Biol Lett, doi:10.1098/rsbl.2011.0523) introduces the lognormal
ACRV model for morphological data, with discretisation by midpoint of equal-
mass bins, expressed "relative to $\delta$" (the mean rate parameter). Wagner
does not give an explicit renormalisation formula in the paper, but the
"relative to $\delta$" framing implies rescaling — the formal calculation
becomes Harrison & Larsson 2014's eq. 2.

Wright & Wynd 2024 (bioRxiv 2024.06.26.600858) is cited in the template at
line 41 as a justification for $n_{\text{Cat}} = 6$ categories; it does not
prescribe a renormalisation rule and inherits whatever the underlying
implementation (RB's `fnDiscretizeDistribution`) provides.

### Is `fnDiscretizeDistribution` correct *as a generic discretiser*?

Yes. As inspected at
`src/core/functions/distribution/DiscretizeDistributionFunction.cpp:99-108`,
it is a generic equal-mass bin median tool. The contract is: the user
supplies a distribution whose mean is whatever they want it to be, and gets
back $n$ representative quantile medians. If the user passes a distribution
with mean $\neq 1$, the output has mean $\neq 1$.

But that contract is violated in the RB by\_nt\_9v template: the user passes
`dnLognormal(0, σ)`, whose true mean is $e^{\sigma^2/2}$, not 1. To make the
template correct under the generic-discretiser contract, either:

(a) Parameterise the lognormal with $\mu_{\log} = -\sigma^2/2$ so the
    continuous mean is 1:
    `dnLognormal(-rate_log_sd^2/2, rate_log_sd)`. Caveat: the discrete
    midpoint mean is then still slightly < 1 for finite $n$ (the Harrison &
    Larsson 2014 deficit), so a posterior shift remains at large $\sigma$.

(b) Both adjust $\mu$ *and* renormalise. This is the H&L 2014 prescription
    and exactly matches MkPrime.

(c) Leave $\mu = 0$ and renormalise: this gives $r_i$'s with mean exactly 1
    but with a different shape (more right-skew) than (b). At small $\sigma$
    the two are nearly identical; at $\sigma \gtrsim 1$ they diverge.
    Defensible but not the published convention.

### Verdict (Issue 2)

**MkPrime follows the published convention (Harrison & Larsson 2014 eq. 2;
Yang 1994 style); RB does not.** RB's template silently breaks the
discretiser's "supply a mean-1 distribution" contract, by passing in a
lognormal with mean $e^{\sigma^2/2}$ and not renormalising. The user's
recollection that they had lifted the formula from a published demo is
consistent — the **MkPrime** ACRV code matches the demo (Harrison & Larsson
2014). The **RB template** does not.

### Recommended canonical form (Issue 2)

For $n$ rate categories with shape parameter $\sigma > 0$:

$$\boxed{\;r_i \;=\;
  \frac{\exp\!\big(-\sigma^2/2 \;+\; \sigma \, \Phi^{-1}((i-1/2)/n)\big)}
       {\tfrac{1}{n}\sum_{j=1}^{n}\exp\!\big(-\sigma^2/2 \;+\; \sigma \, \Phi^{-1}((j-1/2)/n)\big)},
  \quad i = 1,\dots,n.\;}$$

with $\Phi^{-1}$ the standard normal quantile. Properties:

1. $\frac{1}{n}\sum_i r_i = 1$ exactly, by construction.
2. As $n \to \infty$, $r_i$ approaches the corresponding continuous-lognormal
   quantile (its underlying lognormal has mean 1 by $\mu = -\sigma^2/2$).
3. At $\sigma = 0$: $r_i = 1$ for all $i$ (degenerate; MkPrime handles this
   explicitly at `R/acrv.R:16-18` and `src/mcmc_likelihood.cpp:77`).

The $\mu = -\sigma^2/2$ shift in the exponent and the $\Sigma^{-1}$ rescaling
are *both* required to match Harrison & Larsson 2014 eq. 2. (Either alone
gives a slightly wrong distribution shape; together they are exact.)

### Patch suggestion (Issue 2)

Apply to RB, *not* MkPrime. In `dev/rb-equivalence/templates/by_nt_9v.template.Rev:41`:

Replace
```rev
rate_categories := fnDiscretizeDistribution( dnLognormal( 0, rate_log_sd ), 6)
```
with
```rev
# Harrison & Larsson 2014 eq. 2: μ = −σ²/2, then renormalise to discrete mean 1
raw_log_rates := fnDiscretizeDistribution( dnLognormal( -rate_log_sd^2/2, rate_log_sd ), 6 )
rate_categories := raw_log_rates * 6 / sum(raw_log_rates)
```

(Verify RB Rev syntax for `sum()` on a deterministic simplex/vector and for
scalar/vector arithmetic; if `simplex` typing complicates the division, a
small `for`-loop deterministic node may be needed instead. The mathematical
content is what the patch is delivering.)

Apply analogous fix to `by_nt_kv.template.Rev` and any other RB templates that
use `fnDiscretizeDistribution(dnLognormal(...))` directly.

The fix is single-site on the RB side (one template line per template). No
patch attached because changes are to the dev/rb-equivalence/ templates, not
to MkPrime source, and the user has not asked for a Rev-side trivial-fix
patch.

Caveat (Issue 2): the recommended fix changes the meaning of `rate_log_sd` on
the RB side at large $\sigma$ (because the new RB rates are no longer
proportional to the old ones), so post-patch RB posteriors on `rate_log_sd`
are not directly comparable to pre-patch RB samples. Should be picked up by
**mcmc-diagnostician** lane during the re-validation sweep.

---

## Reconciling with the prior investigation note

The prior note (`cross-sampler-rhat-investigation.md`) has the *direction* of
all four implicated effects correct (`tree_length`, `rate_log_sd`,
`rate_neo` drift; `rate_loss` and tree-CID clean), but two specific quantitative
claims should be corrected:

- Line 7, §2: "both samplers use `qnorm((i+0.5)/nCat, μ=−σ²/2, σ=rate_log_sd)`,
  but MkPrime then renormalises so the arithmetic mean of the discrete rates
  is exactly 1, while RB's `fnDiscretizeDistribution` does NOT renormalise" —
  this **understates** the RB-side discrepancy. RB uses $\mu = 0$, not
  $-\sigma^2/2$. The RB-side mean is $\approx e^{\sigma^2/2}$ (1.46 at σ=1),
  not 0.88. The 13% smoke-run gap is consistent with σ≈0.55, but at higher
  σ the discrepancy grows much faster than the prior note predicted.
- §2 line 69: "for σ=1 the empirical mean of six bins is ~0.885 (verified
  Rscript)" — that number is for MkPrime *before* renormalisation; the
  comparable number for RB is ~1.459. The number 0.885 is the RB pre-fix
  mean only if one (mistakenly) imagines RB uses $\mu = -\sigma^2/2$.

The §3 partition-rate analysis is correct as written. Ranking:
issue 1 dominates `rate_neo`; issue 2 dominates `tree_length` (worse than
the prior note thought); both contribute to `rate_log_sd`.

---

## Summary table

| | partition-rate | ACRV discretisation |
|---|---|---|
| Definition of "comparable units" | $\frac{1}{n}\sum_i n_i \rho_i = 1$ | $\frac{1}{n}\sum_i r_i = 1$ |
| RB satisfies it? | **Yes** (algebraic identity, any $r$, any $n_i$) | **No** (mean $\approx e^{\sigma^2/2}$ for $\sigma > 0$) |
| MkPrime satisfies it? | **No** (mean $= (n_1 r + n_2)/n$, = 1 only at $r=1$) | **Yes** (explicit Σ-rescaling) |
| Published convention | — (RB's formula is itself the convention; no upstream cite) | Harrison & Larsson 2014 eq. 2; Yang 1994 |
| Recommended patch site | MkPrime: `src/mcmc_likelihood.cpp:2371-2372` + ~12 trans pruning call-sites | RB: `templates/*.Rev` `fnDiscretizeDistribution` lines |

---

## Files referenced

- `dev/rb-equivalence/templates/by_nt_9v.template.Rev:1-87`
- `dev/rb-equivalence/notes/cross-sampler-rhat-investigation.md:7,42-73,75-95`
- `dev/rb-equivalence/out/mkprime_635_by_nt_9v.rds` (smoke run, posterior median σ=0.55, $n_1=5$, $n_2=12$)
- `R/acrv.R:15-40` (R-fallback discretisation with renormalisation)
- `src/mcmc_likelihood.cpp:73-87` (C++ discretisation with renormalisation)
- `src/mcmc_likelihood.cpp:2371-2372` (neo edge length scaling; one-sided)
- `src/mcmc_likelihood.cpp:2434-2534` (trans/known-k partition; `edgeLen` passed unscaled)
- RB upstream (verified via WebFetch):
  - `src/core/functions/distribution/DiscretizeDistributionFunction.cpp:99-108` (generic discretiser, no rescaling)
  - `src/core/functions/distribution/DiscretizeGammaFunction.cpp` (gamma-specific, *does* rescale)
  - `src/core/distributions/math/LognormalDistribution.cpp` (first arg is log-scale mean)
- Harrison & Larsson 2014, Syst Biol, doi:10.1093/sysbio/syu098 (canonical ACRV renormalisation formula)
- Yang 1994, J Mol Evol, doi:10.1007/BF00160154 (discrete rate categories; mean-1 convention)
- Wagner 2012, Biol Lett, doi:10.1098/rsbl.2011.0523 (lognormal ACRV for morphology)
- Wright & Wynd 2024, bioRxiv 2024.06.26.600858 (motivates $n_{\text{Cat}} = 6$)
