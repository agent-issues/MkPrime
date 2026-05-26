# L7 — ACRV (Among-Character Rate Variation) discretisation

## Theorem (informal title)

The discrete-lognormal quantile-midpoint scheme used by MkPrime —
defined in `R/acrv.R` and replicated in `src/mcmc_likelihood.cpp`,
`src/gibbs_partial_cl.h`, `src/node_cl_cache.h` — produces, after the
final rescaling step, a finite distribution
$\{r_1, \dots, r_K\}$ whose arithmetic mean is exactly $1$.
The rescaling is algebraically necessary: the unrescaled
quantile-midpoint rates do **not** have mean $1$ in general; the
implementation's explicit `rates *= K / sum(rates)` (equivalently
`rates / mean(rates)`) makes the rescaled mean equal $1$ identically
(to floating-point precision).

We further bound the behaviour as $\sigma \to \infty$: the rescaled
top category $r_K$ approaches $K$ from below (a *hard ceiling*), the
bottom category $r_1$ decays as $\Theta\!\bigl(K \exp(-\sigma\,
\Phi^{-1}(1-1/(2K)))\bigr)$, and the ratio $r_K/r_1$ blows up
exponentially in $\sigma$. The ceiling matters: it caps the maximum
effective branch-length scaling at $K \times \text{(largest edge)}$,
which is the mechanism by which the likelihood saturates and can mask
state drift (per project memory
`feedback_acrv_saturation_masks_state_drift`).

## Assumptions

1. $K \ge 1$ is the integer number of rate categories
   (`nCat`, `data->nCat`, argument to `DiscreteLognormalRates`).
2. $\sigma \ge 0$ is `rateLogSd`. The code branches at
   $\sigma \le 0$ and returns $K$ equal rates of $1$
   (`R/acrv.R:16-18`, `src/mcmc_likelihood.cpp:77`,
   `src/gibbs_partial_cl.h:187`, `src/node_cl_cache.h:43`).
3. The underlying continuous distribution of rates across characters is
   $\operatorname{LN}(\mu, \sigma^2)$ with $\mu = -\sigma^2/2$, the
   unique choice making the continuous mean equal $1$
   (lognormal mean is $\exp(\mu + \sigma^2/2)$; this is set at
   `R/acrv.R:21` and `src/mcmc_likelihood.cpp:78`).
4. The quantile function $\Phi^{-1}$ is computed by R/Rmath
   (`qnorm`/`R::qnorm`); the code uses the **midpoint** quantiles
   $(i-\tfrac12)/K$ for $i=1,\dots,K$, not the bin boundaries
   $i/K$ — see `R/acrv.R:34` and the precomputation
   `src/mcmc_likelihood.cpp:2862`.
5. Arithmetic is double-precision IEEE 754; standard cancellation
   error of order $K\,\varepsilon_{\text{mach}}$ is tolerated.
   The empirical check at the end of this proof confirms the rescaled
   mean is $1$ to within $1\,\text{ulp}$ for $\sigma\in[10^{-2},
   20]$.

## Statement

Define
$$
z_i \;=\; \Phi^{-1}\!\bigl((i-\tfrac{1}{2})/K\bigr),
\qquad i = 1,\dots,K,
$$
$$
\tilde r_i \;=\; \exp\!\bigl(\mu + \sigma\, z_i\bigr)
      \;=\; \exp\!\bigl(-\sigma^2/2 + \sigma\, z_i\bigr),
\qquad
S \;=\; \sum_{i=1}^{K}\tilde r_i,
$$
$$
r_i \;=\; \frac{K}{S}\,\tilde r_i. \qquad\qquad (\text{rescaled rates})
$$

**Claim 1 (mean exactly $1$).** $\dfrac{1}{K}\sum_{i=1}^{K} r_i \;=\; 1$
identically for all $K\ge 1$ and all $\sigma\ge 0$.

**Claim 2 (ceiling).** For $\sigma \ge 0$ and $K \ge 2$, $\;r_K \le K$,
with $r_K \to K$ as $\sigma \to \infty$.

**Claim 3 (floor decay).** For fixed $K\ge 2$, as $\sigma\to\infty$,
$$
r_1 \;\sim\; K\,\exp\!\Bigl[-\sigma\bigl(z_K - z_1\bigr) + o(\sigma)\Bigr],
$$
so $r_1 \to 0$ super-exponentially fast (note $z_K = -z_1 > 0$ by the
symmetry of $\Phi^{-1}$ about $1/2$).

**Claim 4 (degeneracy).** $\sigma = 0$ ⇒ $r_i = 1$ for all $i$;
$K = 1$ ⇒ $r_1 = 1$.

## Proof

### Claim 1: mean exactly 1.

By construction,
$$
\frac{1}{K}\sum_{i=1}^{K} r_i
   \;=\; \frac{1}{K}\sum_{i=1}^{K}\frac{K}{S}\,\tilde r_i
   \;=\; \frac{1}{S}\sum_{i=1}^{K} \tilde r_i
   \;=\; \frac{S}{S}
   \;=\; 1,
$$
provided $S>0$. Since each $\tilde r_i = \exp(\cdot)>0$, we have
$S>0$. The result is therefore an algebraic identity in $K$ and the
$\tilde r_i$; it makes no use of the specific values of $z_i$ — in
particular, it holds for **any** $\sigma$ (large or small) and is
robust to the quadrature choice (midpoint, mean, anything). $\square$

Note the subtle point. Yang 1994's median scheme **without**
rescaling has mean different from $1$. The unrescaled mean is
$$
\bar{\tilde r} \;=\; \frac{1}{K}\sum_{i=1}^{K}
        \exp(-\sigma^2/2 + \sigma z_i),
$$
which (by Jensen, since $\exp$ is convex and $\sum z_i = 0$ by
symmetry of $z_i$ about $0$) is **strictly greater than**
$\exp(-\sigma^2/2)$ — i.e. it can be much less than $1$ when $\sigma$
is moderate, and grows without bound as $\sigma$ grows because the
top category dominates. So rescaling is the load-bearing step.

### Claim 2: ceiling $r_K \le K$.

By definition $r_K = K \tilde r_K / S$. Each $\tilde r_i > 0$, so
$S > \tilde r_K$ if $K \ge 2$, hence $r_K < K$. Equality is
approached as $\sigma \to \infty$: the categories $i<K$ shrink
super-exponentially relative to $\tilde r_K$ (their log-ratio is
$\sigma(z_i - z_K) < 0$ with $|z_i - z_K| \ge |z_{K-1}-z_K|>0$),
so $S/\tilde r_K \to 1$ and $r_K \to K$. $\square$

### Claim 3: floor decay.

For large $\sigma$, $S$ is dominated by $\tilde r_K$:
$S = \tilde r_K(1 + O(\exp(-\sigma\,(z_K-z_{K-1}))))$.
Therefore
$$
r_1 \;=\; \frac{K\,\tilde r_1}{S}
   \;\sim\; K\,\frac{\tilde r_1}{\tilde r_K}
   \;=\; K\,\exp(\sigma(z_1 - z_K)).
$$
With $z_K = -z_1 > 0$ this gives the stated $K\exp(-\sigma\cdot 2 z_K)$
behaviour. $\square$

### Claim 4: degeneracy.

- $\sigma=0$: handled by an explicit branch
  (`R/acrv.R:16`, `src/mcmc_likelihood.cpp:77`,
  `src/gibbs_partial_cl.h:187`, `src/node_cl_cache.h:43`). Even
  ignoring the branch, the formula gives
  $\tilde r_i = \exp(0 + 0\cdot z_i) = 1$ for all $i$, $S=K$, and
  $r_i = K\cdot 1 / K = 1$.
- $K=1$: $z_1 = \Phi^{-1}(0.5) = 0$, $\tilde r_1 = \exp(-\sigma^2/2)$,
  $S = \tilde r_1$, $r_1 = 1\cdot\tilde r_1/\tilde r_1 = 1$. $\square$

## Implementation cross-check

The midpoint–rescale recipe is implemented in **four** places that
must agree. They do.

| Location | Midpoint quantiles | Mean formula | Rescale |
|---|---|---|---|
| R | `R/acrv.R:34` `midpoints <- (seq_len(nCat) - 0.5) / nCat` then `qnorm(midpoints, mean = mu, sd = rateLogSd)` at line 35 | `R/acrv.R:21` `mu <- -rateLogSd^2 / 2` | `R/acrv.R:39` `rates / mean(rates)` |
| C++ likelihood | `src/mcmc_likelihood.cpp:2862` precomputes $z_i = $`R::qnorm((i+0.5)/nCat, 0, 1, 1, 0)`; consumed at `src/mcmc_likelihood.cpp:82` `rates[i] = std::exp(mu + rateLogSd * acrvZ[i])` | `src/mcmc_likelihood.cpp:78` `double mu = -rateLogSd * rateLogSd / 2.0` | `src/mcmc_likelihood.cpp:85` `rates[i] *= nCat / total` |
| Gibbs partial-CL | consumes `data->acrvZ` at `src/gibbs_partial_cl.h:192` | `src/gibbs_partial_cl.h:188` | `src/gibbs_partial_cl.h:195` |
| Node-CL cache | consumes `data.acrvZ` at `src/node_cl_cache.h:48` | `src/node_cl_cache.h:44` | `src/node_cl_cache.h:51` |

Notes on the cross-check:

1. **Index offset.** R uses `(i - 0.5)/K` for `i = 1..K`; C++ uses
   `(i + 0.5)/K` for `i = 0..K-1`. These are the same sequence,
   $\{(2j-1)/(2K) : j=1\dots K\}$. ✓
2. **`mu` matches** in both languages: $-\sigma^2/2$
   (`R/acrv.R:21`; `src/mcmc_likelihood.cpp:78`). ✓
3. **Rescale identity.** Three of four implementations
   (C++) use `rates[i] *= nCat / total`. The R implementation uses
   `rates / mean(rates)`. These are algebraically identical:
   $K/S = 1/\bar{\tilde r}$. ✓
4. **The pruning kernels themselves** (`src/acrv.cpp:30`,
   `src/acrv.cpp:152`, `src/acrv.cpp:264`) take `rate_multipliers`
   as input; they do not synthesise rates. The mixture is computed
   correctly: per-site likelihoods are summed across categories
   (`src/acrv.cpp:120`) and then divided by `nCat` before logging
   (`src/acrv.cpp:128`, `src/acrv.cpp:244`, `src/acrv.cpp:354`).
   This is the uniform-mixture-weight assumption: $p(\text{cat}=i) =
   1/K$. This is the standard treatment (Yang 1994) and is correct
   given the discretisation construction above.

### What `src/acrv.cpp` itself does

The file `src/acrv.cpp` contains only the **pruning kernels**
(`pruning_jc_acrv`, `pruning_jc_acrv_collapsed`, `pruning_mkn_acrv`).
It does not construct rate multipliers — those arrive as the
`rate_multipliers` parameter. The discretisation lives in `R/acrv.R`
and the three C++ helpers above. The lane task wording assumes the
discretisation is in `src/acrv.cpp`, which is not the case; the
spec's "src/acrv.cpp:LINE" requests are answered with the actual
locations above.

## Edge cases

| Case | Outcome | Source verification |
|---|---|---|
| $\sigma = 0$ | All `rates[i] = 1` via branch | `R/acrv.R:16–18`; `src/mcmc_likelihood.cpp:77`; `src/gibbs_partial_cl.h:187`; `src/node_cl_cache.h:43` |
| $K = 1$ | `rates = c(1)` (midpoint at $\Phi^{-1}(0.5)=0$) | algebra; one ACRV call always passes ≥1 |
| $K = 2$ | $z_1=-z_2 \approx -0.6745$; rates $(r_1, r_2)$ symmetric in log space and rescaled to mean $1$ | numerical check, see below |
| $\sigma$ huge ($=20$) | rates span 24 orders of magnitude; $r_1 \approx 5.66\times 10^{-24}$, $r_K \to 6.0$ for $K=6$; rescaled mean still $1$ | numerical check, see below |
| Negative $\sigma$ | Same branch as $\sigma=0$ (the `<= 0` check); legitimate because $\sigma$ is a standard deviation and the model prior is Gamma (positive) — see `R/MkPrimeModel.R:457–464`, where the prior log-density at $\sigma=0$ is $-\infty$ when `rateLogSdShape > 1` | enforced at prior level |

### Numerical confirmation (mean = 1 to machine precision)

A single R one-liner reproducing the implementation analytically:

```r
d <- function(sd, K = 6) {
  if (sd <= 0) return(rep(1, K))
  mu  <- -sd^2 / 2
  mid <- (seq_len(K) - 0.5) / K
  r   <- exp(qnorm(mid, mu, sd))
  r / mean(r)
}
```

Output across $\sigma \in \{0.01, 0.5, 1, 2, 5, 10, 20\}$ at $K=6$:

```
sd=  0.01 mean=1.0000000000000000 min=9.8623e-01 max=1.0139e+00
sd=  0.50 mean=0.9999999999999999 min=4.5369e-01 max=1.8088e+00
sd=  1.00 mean=1.0000000000000000 min=1.7191e-01 max=2.7324e+00
sd=  2.00 mean=1.0000000000000000 min=1.6965e-02 max=4.2861e+00
sd=  5.00 mean=1.0000000000000000 min=5.7297e-06 max=5.8130e+00
sd= 10.00 mean=0.9999999999999999 min=5.8244e-12 max=5.9949e+00
sd= 20.00 mean=1.0000000000000000 min=5.6635e-24 max=6.0000e+00
```

The two cases displaying `0.9999999999999999` are off by $1\,\text{ulp}$
($\varepsilon_{\text{mach}} \approx 2.22\times 10^{-16}$),
which is the expected double-precision rounding error of summing $K$
positive doubles and dividing.

The ceiling behaviour $r_K \to K$ is also clearly visible:
$r_K = 5.81, 5.99, 6.00\dots$ for $\sigma = 5, 10, 20$.

## Behaviour as $\sigma \to \infty$ (state-drift implications)

Per `Claim 2`, the top rate category is hard-bounded by $K$. The
**effective branch length** seen by the top rate category is at most
$K \times \text{(edge length)}$ — for $K=6$, never more than $6\times$
the actual edge. This is the mechanism by which the likelihood
saturates as $\sigma$ inflates: the top category cannot scale rates
arbitrarily; instead, the top category looks like a fixed multiplier
$\le K$ on the true branch length, while the bottom $K-1$ categories
contribute essentially constant-site likelihood (they have rates
$\le r_{K-1}$ which still go to $0$, just less rapidly than $r_1$).

Concretely, for large $\sigma$ and $K=6$:
- one category sees rate $\approx 6$ on every branch (highly
  saturating);
- five categories see rate $\to 0$ (nearly constant);

so each site likelihood saturates to roughly
$\tfrac{1}{K}\bigl(\text{lik at rate }K\bigr) +
 \tfrac{K-1}{K}\bigl(\text{lik at rate }0\bigr)$.
The second term equals the equilibrium-frequencies prior at the root,
which is data-independent. That is the **free-equilibrium floor**
referenced in `feedback_acrv_saturation_masks_state_drift.md`: the
model can place $(K-1)/K$ of its mixture weight on "this site is
constant at frequency $\pi$" without paying a likelihood penalty,
while the one "saturating" category does little to discriminate
states because it has scaled time $\approx K\cdot t$ which is already
near equilibrium. The bound `logL` is therefore not evidence that
state has stabilised — exactly the warning in the memory item.

This is a **structural feature** of the renormalised midpoint scheme
(not a bug). It is shared with discrete-gamma. The conclusion is
operational: use the prior $p(\sigma)$ to keep $\sigma$ in a sensible
range; do not use bounded `logL` as a proxy for "state has converged".

## Comparison to Yang 1994 truncated-mean

Yang (1994, *J Mol Evol* 39:306–314) discusses two
discretisation schemes for the gamma:

- **Mean (truncated-mean) approximation.** $r_i = E[R \mid R \in
  B_i]$ where $B_i = [q_{(i-1)/K}, q_{i/K}]$ is the $i$-th
  $K$-tile of the continuous distribution. **By construction**
  $\frac{1}{K}\sum r_i = E[R] = 1$; no rescaling is needed.
- **Median (midpoint) approximation.** $r_i = q_{(i-\tfrac12)/K}$.
  This does *not* have mean $1$, so Yang explicitly recommends
  rescaling: "we may rescale ... to make the mean of $r$ equal
  to $1$" (Yang 1994 §"Discrete approximation").

The MkPrime implementation does the **median** scheme on the lognormal
**with rescaling**. That matches Yang's median-with-rescaling
prescription, transported from gamma to lognormal. Felsenstein 2004
§16.3 (and Yang 2014 §4.3) confirm this is the standard practice.

**Tradeoff.** The mean (truncated-mean) scheme has lower
discretisation bias for moderate $\sigma$ but requires evaluating
$E[R \mid R\in B_i]$, which for the lognormal is the truncated-lognormal
mean
$$
E[R \mid R \in [a,b]]
   \;=\; \frac{\Phi\bigl(\tfrac{\ln b - \mu}{\sigma} - \sigma\bigr)
              - \Phi\bigl(\tfrac{\ln a - \mu}{\sigma} - \sigma\bigr)}
              {\Phi\bigl(\tfrac{\ln b - \mu}{\sigma}\bigr)
              - \Phi\bigl(\tfrac{\ln a - \mu}{\sigma}\bigr)}\,
        \exp(\mu+\sigma^2/2),
$$
which is more expensive than a single `qnorm`/`exp`. Yang (1994)
reports the median scheme is "almost as accurate" as the mean scheme
for $K\ge 4$ when applied to the gamma; the same intuition transports
to the lognormal because the construction is purely a quadrature
choice. The implementation choice (median + rescale) is the standard
fast option and does not bias any claim of the MkPrime model.

## Verdict

**Watertight.** The mean-1 property is an algebraic identity
($r_i = K\tilde r_i / S \Rightarrow \tfrac{1}{K}\sum r_i = 1$)
that holds for **any** $\sigma\ge 0$, any $K\ge 1$, and any choice of
quadrature points $\{z_i\}$. The implementation is correct in all
four locations (`R/acrv.R:39`, `src/mcmc_likelihood.cpp:85`,
`src/gibbs_partial_cl.h:195`, `src/node_cl_cache.h:51`). Behaviour at
large $\sigma$ is well-characterised (top rate saturates at $K$,
bottom rate decays super-exponentially), confirming the spirit of
`feedback_acrv_saturation_masks_state_drift`: the rescaled
midpoint scheme imposes a *finite* ceiling on the per-category
rate-multiplier, which is the mechanism by which `logL` can bound
while state drifts in latent dimensions. This is a feature of the
mathematics (and of the analogous discrete-gamma), not a bug. No
implementation patch attached.
