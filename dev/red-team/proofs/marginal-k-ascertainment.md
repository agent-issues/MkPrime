# Ascertainment correction for the marginal-k (Rao-Blackwellised) likelihood under `coding="variable"`

> **⚠ SUPERSEDED (2026-06-02).** This proof derives the **ratio-of-sums**
> $P = \frac{\sum_k w_k L(y\mid k)}{\sum_k w_k a_k}$ for the **pair-rejection**
> forward (`e321036`), where the simulator rejects-and-redraws *both* `k'` and
> the data jointly. HEAD uses the **redraw-data-only** forward (`d337750`:
> `k'_i` drawn once unconditionally, then only the *data* is redrawn until the
> character is variable). Under that forward the correct estimator is the
> **sum-of-ratios** $\sum_k w_k L(y\mid k)/a_k$ — which is what the production
> code computes — and the only missing factor was the truncation normaliser
> `Z(p)` (fixed in MARGINAL-K-TRUNC-001, proof
> `marginal-k-truncation-normaliser.md`). Kept for the record; do **not** treat
> the ratio-of-sums conclusion below as describing current code.

## Theorem (informal title)

Under Model A (`k'_i ~ 2 + Geo(p)`, unconditional) and `coding="variable"`, the
correct per-character marginal data likelihood for an observed *variable*
character, with the latent state count `k` marginalised out, is the
**ratio-of-sums**

$$P(y_i \mid p, \mathrm{tree}, \mathrm{var}) \;=\; \frac{\sum_{k} w_k(p)\, L(y_i\mid k)}{\sum_{k} w_k(p)\, a_k},$$

**not** the sum-of-ratios $\sum_k w_k(p)\,L(y_i\mid k)/a_k$. The denominator
$D(p,\mathrm{tree}) = \sum_k w_k(p)\,a_k$ is a per-character ascertainment
normaliser that depends only on $(p,\mathrm{tree})$, is shared across
characters, and must be applied **once per character, outside the k-sum**. The
implementation instead applies a per-$k$ factor $1/a_k$ *inside* the k-sum and
applies no outer normaliser. This is an **implementation bug** in both the
marginal likelihood and the sampled-k Gibbs full-conditional.

## Assumptions

1. Branch lengths are non-negative reals; the tree is fixed within a single
   likelihood evaluation.
2. The substitution model is symmetric JC/Mk on $k$ states (so
   `const_site_prob_for_k` returns the constant-site probability under uniform
   root frequencies; verified at `src/mcmc_likelihood.cpp:2300-2310`).
3. `coding="variable"` $\Leftrightarrow$ `data.codingType == 1`; the
   ascertainment event is $E = \{y\ \text{is variable}\} = \{y\ \text{shows} \ge 2$
   distinct states$\}$, the complement of the constant-site event.
4. $a_k := P(\mathrm{variable}\mid k,\mathrm{tree}) = 1 - \mathrm{csp}_k$, where
   $\mathrm{csp}_k = $ `const_site_prob_for_k(..., k, ...)` is the constant-site
   probability $P(\mathrm{invariant}\mid k,\mathrm{tree})$. **Verified below.**
5. $L(y_i\mid k) = P(y_i\mid k,\mathrm{tree})$ is the *raw* (unconditioned) Mk
   pruning likelihood of pattern $y_i$ given $k$ states.
6. Model A prior weight on the state count: $w_k(p) = p(1-p)^{k-2}$, $k\ge 2$.
   (Model B conditional weight: $w_k^{B}(p) = p(1-p)^{k-k_{\mathrm{obs}}}$.)
7. The SBC forward simulator draws each *observed* character by the
   resample-until-variable scheme: draw $k\sim w_k(p)$, draw $y\sim P(\cdot\mid
   k,\mathrm{tree})$, reject the **whole pair** $(k,y)$ if $y$ is invariant, and
   repeat. (Stated in the brief; matches commit `e321036`
   "resample-until-variable forward".)

### Verification of Assumption 4 (csp = P(invariant))

`const_site_prob_for_k` (`src/mcmc_likelihood.cpp:2280-2311`) returns, for
`codingType == 1` and no Q-heterogeneity,

```
constant_site_prob_jc(parent, child, edgeLen, nTip, kStates, rootFreqs, acrvRates)
```

with `rootFreqs = 1/k` uniform (`:2300-2302`). This is exactly
$\sum_{s} \pi_s\,P(\text{all tips} = s) = P(\mathrm{invariant}\mid k)$. The
`codingType == 2` (informative) branch *adds* the singleton-site probability
(`:2305-2308`); under `coding="variable"` (`codingType == 1`) that branch is not
taken, so $\mathrm{csp}_k = P(\mathrm{invariant}\mid k,\mathrm{tree})$ exactly,
hence $a_k = 1 - \mathrm{csp}_k = P(\mathrm{variable}\mid k)$. $\square$

## Statement

Fix $(p,\mathrm{tree})$. Let the generative model for a single observed
character be the rejection sampler of Assumption 7. Then:

**(I) Marginal likelihood.** The accepted pattern $y$ has density
$$P(y\mid p,\mathrm{tree},\mathrm{var}) = \frac{\sum_{k\ge 2} w_k(p)\,L(y\mid k)}{D(p,\mathrm{tree})}, \qquad D(p,\mathrm{tree}) = \sum_{k\ge 2} w_k(p)\,a_k.$$

**(II) Full conditional of k.** For an observed variable $y_i$,
$$P(k\mid y_i,p,\mathrm{tree},\mathrm{var}) = \frac{w_k(p)\,L(y_i\mid k)}{\sum_{k'\ge 2} w_{k'}(p)\,L(y_i\mid k')} \;\propto\; w_k(p)\,L(y_i\mid k).$$
The factor $a_k$ **cancels** and must not appear in the k-sampling weight.

**(III) Normaliser placement.** $D(p,\mathrm{tree})$ is independent of $y_i$.
With $n$ observed (independent) variable characters the total log-likelihood is
$$\log P(\mathbf{y}\mid p,\mathrm{tree},\mathrm{var}) = \sum_{i=1}^{n}\Big[\log\!\sum_{k} w_k(p) L(y_i\mid k)\Big] \;-\; n\log D(p,\mathrm{tree}).$$
(If $k_{\mathrm{obs},i}$ varies per character and one floors the k-sum at
$k=k_{\mathrm{obs},i}$, see Edge cases for the exact range bookkeeping; the
$-\log D$ correction is unchanged in form.)

## Proof

**Step 1 (joint of accepted pair).** The rejection sampler accepts $(k,y)$ iff
$y\in E$ (variable). Standard rejection sampling: the accepted draws have the
joint density equal to the proposal joint restricted to $E$ and renormalised. The
proposal joint is $w_k(p)\,P(y\mid k)$ (draw $k$, then $y$). Hence
$$P(k,y\mid \mathrm{var}) = \frac{w_k(p)\,P(y\mid k)\,\mathbb{1}[y\in E]}{Z}, \quad Z = \sum_{k'} w_{k'}(p)\sum_{y'\in E} P(y'\mid k').$$
Now $\sum_{y'\in E} P(y'\mid k') = P(\mathrm{variable}\mid k') = a_{k'}$, so
$Z = \sum_{k'} w_{k'}(p)\,a_{k'} = D(p,\mathrm{tree})$. (Ref: rejection-sampling
identity, e.g. Robert & Casella 2004, *Monte Carlo Statistical Methods*, §2.3;
ascertainment-bias treatment of variable-only morphological data, Lewis 2001,
*Syst. Biol.* 50:913, eq. for the Mkv conditional likelihood.)

**Step 2 (marginal over k — claim (I)).** For an observed (hence variable)
$y$, marginalise k from the Step-1 joint:
$$P(y\mid \mathrm{var}) = \sum_{k} P(k,y\mid\mathrm{var}) = \frac{\sum_{k} w_k(p)\,P(y\mid k)\,\mathbb{1}[y\in E]}{D} = \frac{\sum_{k} w_k(p)\,L(y\mid k)}{D},$$
using $\mathbb{1}[y\in E]=1$ for the observed variable $y$ and
$P(y\mid k)=L(y\mid k)$. This is the ratio-of-sums. $\square$

**Step 3 (sum-of-ratios is a different object).** The implementation's
per-character contribution is (in probability space)
$\sum_k w_k(p)\,L(y\mid k)/a_k$. Compare to (I): the $1/a_k$ sits *inside* the
sum and is $k$-dependent, so it does **not** factor out as the constant $1/D$.
The two agree only in the degenerate case where $a_k$ is independent of $k$ over
the support (e.g. a single candidate $k$, or $a_k\equiv a$), where both reduce
to $a^{-1}\sum_k w_k L$. In general they differ, and the sum-of-ratios is not a
normalised density in $y$: $\sum_{y\in E}\sum_k w_k L(y\mid k)/a_k = \sum_k w_k
\cdot (\sum_{y\in E}L(y\mid k))/a_k = \sum_k w_k\cdot a_k/a_k = \sum_k w_k \ne 1$
in general, and crucially $\ne D$, so no single outer constant repairs it.

**Step 4 (full conditional — claim (II)).** Condition the Step-1 joint on the
observed $y_i$:
$$P(k\mid y_i,\mathrm{var}) = \frac{P(k,y_i\mid\mathrm{var})}{\sum_{k'} P(k',y_i\mid\mathrm{var})} = \frac{w_k(p)L(y_i\mid k)/D}{\big(\sum_{k'} w_{k'}(p)L(y_i\mid k')\big)/D} = \frac{w_k(p)L(y_i\mid k)}{\sum_{k'} w_{k'}(p)L(y_i\mid k')}.$$
$D$ cancels (numerator and denominator carry the same $1/D$), and **no $a_k$
appears**: the only $k$-dependence is through $w_k$ and $L(\cdot\mid k)$. Adding
a $1/a_k$ factor to the k-sampling weight therefore tilts the conditional toward
large-$k$ candidates (where $a_k$ is larger... note $a_k$ *grows* with $k$ on a
typical tree, so $1/a_k$ down-weights large $k$ — see Step 6). $\square$

**Step 5 (numerical confirmation).** A toy rejection sampler ($k\in\{2,3,4\}$,
geometric $w_k$, three patterns one invariant) was simulated $4$-$6\times10^5$
times:

- Marginal over variable patterns: empirical $(0.6164, 0.3836)$ matched
  ratio-of-sums $(0.6156, 0.3844)$; renormalised sum-of-ratios gave
  $(0.6199, 0.3801)$ — outside Monte-Carlo error.
- Conditional $P(k\mid y\!=\!\text{pattern }2)$: empirical
  $(0.5755, 0.2866, 0.1379)$ matched $w_kL$ $(0.5747, 0.2874, 0.1379)$, not the
  $w_kL/a_k$ form $(0.4938, 0.3086, 0.1975)$.

(One-off check, not committed as a test.)

**Step 6 (sign of the bias — consistency with the observed SBC failure).** On a
non-degenerate tree, larger $k$ makes a randomly evolved character *more* likely
to be variable, so $a_k$ is increasing in $k$ and $1/a_k$ is *decreasing* in
$k$. The erroneous $1/a_k$ factor in the marginal therefore systematically
**down-weights large-$k$ candidates**. Under Model A, large $k$ is favoured by
small $p$ (since $w_k = p(1-p)^{k-2}$ puts more mass on large $k$ as $p\to 0$).
Suppressing large-$k$ likelihood contributions pushes the posterior on $p$
**upward** (toward the small-$k$ regime), producing exactly the reported
rank-spike-at-0 / posterior-biased-high signature on $p$. Because the missing
$-n\log D(p,\mathrm{tree})$ couples $p$ and the tree, `tree_length` inherits the
bias, while `rate_log_sd` — orthogonal to the k-marginalisation weighting — stays
calibrated. This is corroborating (not proof-grade) evidence that the derived
bug is the cause; the SBC rerun after the fix is the confirmation
(mcmc-diagnostician lane).

## Implementation cross-check

| Theorem element | Correct form | Source | Status |
|---|---|---|---|
| $a_k = 1-\mathrm{csp}_k$ | $P(\mathrm{var}\mid k)$ | `src/mcmc_likelihood.cpp:2300-2310` | OK |
| per-$k$ weight | $w_k L^\beta$ (log: $\beta\log L + \log w_k$) | `src/mcmc.cpp:3956-3959, 4087-4099` | **BUG: includes $-\log a_k$** |
| marginal over $k$ | $\log\sum_k w_k L$ | `src/mcmc.cpp:4270-4275` | logSumExp OK, but operates on buggy weights |
| outer normaliser | $-\log D$, $D=\sum_k w_k a_k$ | `src/mcmc.cpp:4276-4281` (uncond branch) | **BUG: absent** |
| Model A vs B factor | $(k_{\mathrm{obs}}-2)\log(1-p)$ | `src/mcmc.cpp:4276-4279, 4306-4309` | OK (orthogonal to this bug) |
| Gibbs k-sampling weight | $\propto w_k L$ | `src/mcmc.cpp:4374-4389` (uses `charLogW`) | **BUG: inherits $1/a_k$** |

Detailed citations:

1. **Per-$k$ ascertainment factor (the offending term).**
   `src/mcmc.cpp:3956-3959`:
   ```cpp
   double csp = getCSP(k);
   double logAscCorr = (coding != 0 && csp < 1.0) ? -std::log(1.0 - csp) : 0.0;
   if (coding != 0 && csp >= 1.0) logAscCorr = R_NegInf;
   ```
   `logAscCorr = -\log(1-\mathrm{csp}_k) = -\log a_k`. Applied at
   `src/mcmc.cpp:4087-4088` (`ll += logAscCorr`) and folded into the stored
   weight at `:4094` (`double w = beta * ll + logPrior_k`) → `:4099`
   (`charLogW[ti*kMaxKprimeCand + ko] = w`). So in log-space the stored per-$k$
   weight is
   $$w_{\text{stored}} = \beta\big(\log L(y_i\mid k) - \log a_k\big) + \log w_k(p),$$
   i.e. weight $\propto w_k(p)\,L(y_i\mid k)^\beta / a_k^{\beta}$. The $1/a_k$
   should not be here per Steps 2 and 4.

2. **Marginal logSumExp, no outer normaliser.**
   `src/mcmc.cpp:4270-4281` performs $\mathrm{charLL} = \mathrm{logSumExp}_{ko}
   (w_{\text{stored}})$ and, for `uncond`, adds $(k_{\mathrm{obs}}-2)\log(1-p)$
   (`:4276-4279`). There is no $-\log D$ term anywhere in the function (confirmed
   over `:4169-4313` and the cache fast-path `:4288-4311`). Hence
   $$\mathrm{charLL}_i = \log\!\sum_k w_k(p) L(y_i\mid k)^\beta / a_k^{\beta} \;+\;(k_{\mathrm{obs},i}-2)\log(1-p),$$
   which at $\beta=1$ is the **sum-of-ratios**, missing $-\log D$. This matches
   the brief's reading.

3. **Gibbs sweep reuses the same buggy weights.**
   `gibbs_kprime_sweep_impl` (`src/mcmc.cpp:4329-4390`) calls
   `compute_per_kprime_log_lik` (`:4347`) to fill `charLogW`, then samples $k$
   categorically from $\exp(\mathrm{charLogW})$ (`:4377-4389`). Since `charLogW`
   carries the $1/a_k^{\beta}$ factor, the realised k-conditional is
   $\propto w_k L^\beta/a_k^\beta$, contradicting claim (II)
   ($\propto w_k L$ at $\beta=1$). **The sampled-k Gibbs path is affected.**

4. **Post-sweep / block-shift relikelihood is consistent within itself but uses
   the same per-$k$ ascertainment convention.** After the Gibbs sweep,
   `cpp_partition_log_likelihood` (`:4397-4403`) recomputes the likelihood at the
   chosen fixed $k$; for a fixed $k$ the trans-partition path applies the same
   `-log(1-csp)` correction (the `coding!=0` ascertainment in the partition
   kernel). This is *internally* a per-character $L(y_i\mid k)/a_k$ — the
   standard Mkv conditional-on-variable likelihood for a *known* $k$ (Lewis
   2001), which is correct **when $k$ is fixed/observed**, but is *not* what the
   marginal-over-$k$ target wants. The two notions of "ascertainment" must not be
   conflated: the Mkv per-character $/a_k$ is correct for the conditional model
   $P(y\mid k,\mathrm{var})$ with *known* $k$; the marginal-$k$ model needs the
   single shared $/D$ instead.

## Edge cases

1. **Single candidate $k$ (e.g. $k$ pinned, or `kMaxKprimeCand` collapses to
   one).** Then $\sum_k w_k L/a_k = w_{k_0}L/a_{k_0}$ and $D = w_{k_0}a_{k_0}$;
   ratio-of-sums $= L/a_{k_0}$ while sum-of-ratios $= w_{k_0}L/a_{k_0}$. These
   differ by the constant $w_{k_0}$, absorbed into the prior normaliser — so the
   bug is *invisible* when only one $k$ is ever live. This is why a saturated /
   single-$k$ toy can pass while a realistic multi-$k$ regime fails, consistent
   with commit `40da2b1` ("move regime off saturated toy").

2. **$a_k$ constant across the live $k$-range.** If $a_k\equiv a$ for all live
   $k$, then $1/a$ factors out of the sum and equals $1/D$ up to the constant
   $\sum_k w_k$. Bug invisible. Occurs only on degenerate trees (e.g. very long
   branches where every $k\ge 2$ is almost surely variable, $a_k\approx 1$).

3. **$\mathrm{csp}_k \to 1$ (very short tree / large $k$ saturating constant
   sites).** Code sets `logAscCorr = R_NegInf` (`:3959`), killing that $k$. Under
   the correct formula such a $k$ has $a_k\to 0$ so $D$ gets a vanishing
   contribution but $w_k L(y_i\mid k)$ stays finite — that $k$ should **not** be
   killed in the numerator. The current `R_NegInf` guard is correct for the
   sum-of-ratios formulation but would need re-examination under the fix:
   the numerator term $w_k L(y_i\mid k)$ must survive even as $a_k\to 0$.
   See Fix note 4.

4. **$k_{\mathrm{obs}}=2$ (Model A floor).** Then $w_k = w_k^{B}$ and the
   $(k_{\mathrm{obs}}-2)\log(1-p)=0$ correction vanishes; Model A and Model B
   coincide. The $-\log D$ bug is present regardless.

5. **`coding="all"` (`codingType==0`).** `getCSP` returns 0, `logAscCorr=0`,
   $a_k=1$, $D=\sum_k w_k$ (constant in tree). No ascertainment, no bug — the
   ratio-of-sums and sum-of-ratios coincide. The bug is specific to
   `coding="variable"`/`"informative"`.

## Exact fix

The fix has two coordinated parts. Let $S_i := \sum_k w_k(p) L(y_i\mid k)$
(numerator) and $D := \sum_k w_k(p) a_k$ (shared denominator).

**Fix part 1 — remove the per-$k$ $1/a_k$ from the weights.** In
`compute_per_kprime_log_lik` (`src/mcmc.cpp:3956-3959, 4087-4090`), drop the
`logAscCorr = -log(1-csp)` contribution from the stored per-$k$ weight so that
$$w_{\text{stored}} = \beta\log L(y_i\mid k) + \log w_k(p).$$
Concretely: do not add `logAscCorr` to `ll` at `:4087-4088`. **Retain** the
`csp >= 1.0 → R_NegInf` guard *only* if the corresponding $k$ truly cannot
appear; but note Edge case 3 — under the correct numerator a $k$ with $a_k\to 0$
should still contribute $w_k L$, so the kill should be removed from the numerator
path (it belongs only to the $D$ accumulation). Safer: keep $L(y_i\mid k)$ finite
and let it enter both $S_i$ and the (separately accumulated) $D$.

This fix part 1 **simultaneously repairs the Gibbs k-conditional**: once
`charLogW` $= \beta\log L + \log w_k$, the categorical draw at
`src/mcmc.cpp:4377-4389` samples $\propto w_k L^\beta$, which at $\beta=1$ is
claim (II). No further change to `gibbs_kprime_sweep_impl` is needed.

**Fix part 2 — add the outer $-\log D$ per character in the marginal.** In
`cpp_log_likelihood_marginal` (`src/mcmc.cpp:4270-4281` and the cache fast-path
`:4302-4310`), after forming $\mathrm{charLL}_i = \log S_i$, subtract $\log D$:
$$\mathrm{charLL}_i \mathrel{-}= \beta\log D(p,\mathrm{tree}).$$
$D = \sum_k w_k(p)\,a_k$ must be computed once per evaluation (it is shared
across all $n$ characters of a partition; it does depend on $p$, the tree, ACRV
rates, and the per-partition $k$-range). Compute it from the cached
$a_k = 1-\mathrm{csp}_k$ (already available via `getCSP`) and the $w_k(p)$
series — a single $O(k_{\max})$ logSumExp:
$$\log D = \mathrm{logSumExp}_k\big(\log w_k(p) + \log a_k\big).$$
**Tempering:** since the impl tempers the whole log-likelihood by $\beta$
(per-$k$ weights carry $\beta\log L$), the ascertainment normaliser is part of
the likelihood and should be tempered consistently, i.e. subtract
$\beta\log D$, not $\log D$. (If the house convention is that ascertainment
normalisers are *not* tempered, subtract $\log D$; this should be reconciled with
how `cpp_partition_log_likelihood` tempers its own $/a_k$ for fixed-$k$ chars.
Flagged as a caveat below.)

**Range bookkeeping for $D$ (per-character $k_{\mathrm{obs}}$).** The k-sum in
$S_i$ runs over $k = k_{\mathrm{obs},i} + ko$, $ko=0..$, with prior weight
$w_k^B = p(1-p)^{ko}$ in the impl (Model A's extra $(k_{\mathrm{obs}}-2)\log(1-p)$
added separately). For $D$ to be the correct *unconditional* Model A normaliser,
it must use the **same Model A weights $w_k(p)=p(1-p)^{k-2}$ over the same support
$k\ge 2$** as the numerator's effective prior. Two consistent options:

  - (a) Compute $S_i$ and $D$ both in Model A weights over $k\ge 2$, then
    $\mathrm{charLL}_i = \log S_i - \beta\log D$ with $D$ truly $y$-independent
    and **shared** (one $D$ for all chars sharing the same candidate $k$-range
    and tree). This is the clean form.
  - (b) Keep the impl's per-character Model B factorisation
    ($S_i^B$ then $+(k_{\mathrm{obs},i}-2)\log(1-p)$). Then $D$ must be the Model
    A normaliser $\sum_{k\ge 2} w_k a_k$, **the same for every character**
    (it does not carry $k_{\mathrm{obs},i}$). Because $D$ is computed over the
    full $k\ge 2$ support, characters with different $k_{\mathrm{obs},i}$ share
    one $D$ only if their candidate $a_k$ values are tabulated over the common
    $k\ge 2$ range. `getCSP` already caches $a_k$ by absolute $k$
    (`src/mcmc.cpp:3763-3776`), so a single $\log D$ over $k=2..k_{\max}$ is the
    right shared quantity. **Caution:** the numerator's per-character $k$-range
    starts at $k_{\mathrm{obs},i}$, not 2; the Model A factor
    $(k_{\mathrm{obs},i}-2)\log(1-p)$ supplies the missing low-$k$ prior mass in
    $S_i$ but those low-$k$ states ($2\le k<k_{\mathrm{obs},i}$) are
    *impossible* for that character (it shows $k_{\mathrm{obs},i}$ states), so
    $L(y_i\mid k)=0$ there and they correctly contribute nothing to $S_i$ — but
    they **do** contribute to $D$ (a character *could* have been variable with
    fewer states). So $D$ must run over the full $k\ge 2$, while $S_i$ runs over
    $k\ge k_{\mathrm{obs},i}$. This asymmetry is correct and is the subtle part
    of the fix.

**Model A vs Model B.** Under **Model B** (production, conditional prior
$w_k^B = p(1-p)^{k-k_{\mathrm{obs}}}$ where the prior is *already* conditioned on
$k\ge k_{\mathrm{obs}}$), the ascertainment normaliser is
$D^B_i = \sum_{k\ge k_{\mathrm{obs},i}} w_k^B(p)\,a_k$, which **is**
$k_{\mathrm{obs},i}$-dependent and therefore *not* shared across characters with
different $k_{\mathrm{obs}}$. So:

  - Model A: one shared $D = \sum_{k\ge 2} w_k a_k$; subtract $\beta\log D$ from
    every character (plus the existing $(k_{\mathrm{obs}}-2)\log(1-p)$ Model A
    factor).
  - Model B: per-$k_{\mathrm{obs}}$ normaliser $D^B_i = \sum_{k\ge
    k_{\mathrm{obs},i}} w_k^B a_k$; subtract $\beta\log D^B_i$ per character.
    Characters sharing $k_{\mathrm{obs}}$ share $D^B$.

In **both** cases the per-$k$ $1/a_k$ inside the sum (Fix part 1) is wrong and
must be removed; the difference is only in how $D$ is aggregated.

A patch was **not** prepared: the fix is non-trivial (it changes the cached
weight semantics in `compute_per_kprime_log_lik`, requires a new $\log D$
accumulation with careful $S_i$-vs-$D$ range asymmetry, and touches the
$\beta$-tempering convention). This exceeds the trivial-fix threshold and is
reported in the body for the implementer.

> Caveat (1): Whether the ascertainment normaliser $\log D$ should be tempered
> by $\beta$ or not depends on the project's Metropolis-coupling convention,
> which I did not fully derive here. Closing this requires confirming how
> `cpp_partition_log_likelihood` tempers its fixed-$k$ $/a_k$ correction and
> matching the marginal path to it. Should be picked up by lane (marginal-k
> tempering) / role math-prover or numerical-auditor.

> Caveat (2): The SBC rerun is the empirical confirmation that this fix removes
> the $p$ rank-spike and the `tree_length` bias. Closing this requires running
> the honest-SBC harness post-fix. Should be picked up by role
> mcmc-diagnostician.

> Caveat (3): The interaction of the `csp >= 1.0 → R_NegInf` guard
> (`src/mcmc.cpp:3959`) with the corrected numerator (Edge case 3) needs a
> numerical check on short-tree regimes where $a_k\to 0$: the numerator term
> $w_k L$ must survive even as the $D$ contribution $w_k a_k\to 0$. Should be
> picked up by role numerical-auditor.

## Verdict

**Implementation bug.** The marginal-k likelihood under `coding="variable"`
computes the **sum-of-ratios** $\sum_k w_k(p)L(y_i\mid k)/a_k^{\beta}$ with no
outer ascertainment normaliser, whereas the correct marginal (derived from the
resample-until-variable generative model and confirmed numerically) is the
**ratio-of-sums** $\big(\sum_k w_k(p)L(y_i\mid k)\big)\big/ D(p,\mathrm{tree})$,
$D = \sum_k w_k(p)a_k$.

- The per-$k$ factor $-\log a_k$ (`src/mcmc.cpp:3956-3959, 4087-4088`) is wrong
  and must be removed.
- The shared outer normaliser $-\log D$ (Model A: one $D$ over $k\ge 2$; Model B:
  per-$k_{\mathrm{obs}}$ $D^B_i$) is missing from
  `cpp_log_likelihood_marginal` (`src/mcmc.cpp:4270-4311`) and must be added.
- **The sampled-k Gibbs full-conditional is also affected** (claim II): it
  currently samples $\propto w_k L/a_k$ but should sample $\propto w_k L$. This
  is repaired automatically by Fix part 1 (removing the per-$k$ $1/a_k$ from
  `charLogW`); no separate edit to `gibbs_kprime_sweep_impl` is required.

The bug's sign (down-weighting large $k$, hence biasing $p$ high and dragging
`tree_length`, while leaving `rate_log_sd` calibrated) is consistent with the
reported SBC failure; the post-fix SBC rerun is the confirmation (Caveat 2). No
patch attached — fix is non-trivial.
