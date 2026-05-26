# Lane L6 — k′-prior arithmetic re-audit (Wave 2)

> **Lane.** L6 — first-principles derivation of the four k′-priors:
> `geometric`, `empirical_geometric`, `beta_geometric`, `logseries`.
>
> **Files in scope.** `R/MkPrimeModel.R::LogPrior`,
> `R/MkPrimeModel.R::.LogPriorEmpiricalGeometric`,
> `R/MkPrimeModel.R::.LogPemp`, `R/data.R::MkPrimeEmpiricalPrior`,
> `src/mcmc.cpp::cpp_log_prior` (lines 200–352), and the empirical body
> distribution at `data/empiricalNObs.rda`.
>
> **Worktree state.** Worktree base is `2bc4da7` (main, pre-Wave-1). The Wave 1
> commit `34a2607` (campaign-2026-05-26 math proofs) is on `main` HEAD but is
> *not* in this worktree (`feedback_verify_agent_base` pattern). The L4
> Wave 1 proof (`dev/red-team/proofs/relabelling-correction.md`) was read via
> `git show 34a2607:...`. This proof reuses L4's verdict that the falling-
> factorial relabelling correction lives in the **log-likelihood** and is
> never added to the log-prior — verified independently below by inspecting
> `LogPrior` and `cpp_log_prior` directly.
>
> **Scope clarification.** The orchestrator's brief asked for re-verification
> of R5-1..R5-5. Those findings belong to a *different* model — the ecology
> sim3 sampler (`gibbs_z_sweep_impl` at `src/mcmc.cpp:4011-4116`, parameters
> `pi0`, `phi`, `theta`, `z`). They are not k′-prior items. They are
> addressed in §3 below in the form their semantics actually admit:
> R5-1/2/3 are prior-arithmetic about `pi0/phi/theta` and survive trivially;
> R5-4 is an empirical claim about post-relabel `phi` medians that depends
> on `sim3-scoring.R::HasBipartSplits` output and cannot be closed by a
> math-prover — it requires a sim3 rerun against the corrected scoring
> (`project_scoring_bug.md` fix). Marked **deferred-out-of-lane**.

---

## 1. Theorem (informal): the four k′-prior normalisations

Let $k'_i \in \{k_{\text{obs},i}, k_{\text{obs},i}+1, \ldots\}$ be the latent
state count for transformational character $i$, with observed $k_{\text{obs},i}
\geq 1$. Let $u_i := k'_i - k_{\text{obs},i} \in \{0, 1, 2, \ldots\}$.
The four implemented priors are:

1. **Geometric (hierarchical on $u$)** — $u_i \mid p \sim \mathrm{Geo}(p)$
   i.i.d., $p \sim \mathrm{Beta}(a,b)$.
2. **Beta-geometric (hierarchical on $u$)** — $u_i \mid \alpha,\beta \sim
   \mathrm{BetaGeo}(\alpha,\beta)$ i.i.d., $\alpha,\beta \sim \mathrm{Exp}(1)$
   independently.
3. **Empirical-geometric (convolution on $k'$)** — $k'_i = N_{\text{obs},i} +
   N_{\text{unobs},i}$ with $N_{\text{obs},i} \sim P_{\text{emp}}$,
   $N_{\text{unobs},i} \mid p \sim \mathrm{Geo}(p)$, $p \sim \mathrm{Beta}(a,b)$;
   then $k'_i$ is truncated at $k_{\text{obs},i}$.
4. **Logseries** — $k'_i \mid c \sim \mathrm{LogSeries}(c)$ on $\{1, 2, \ldots\}$,
   $c \in (0,1)$ fixed, then $k'_i$ is truncated at $k_{\text{obs},i}$.

**Claim.** Priors (1) and (2) are normalised by construction; the
implementation in `R/MkPrimeModel.R` and `src/mcmc.cpp` is correct.
Priors (3) and (4) involve a per-character truncation $k' \geq k_{\text{obs},i}$
whose normaliser $Z_i(\theta) = \sum_{k \geq k_{\text{obs},i}} P(k \mid \theta)$
is **omitted** from the implemented log-density. This is **EG-001** (already
filed) and a new sibling finding **LS-001** (to file) of the same family.

## 2. Assumptions

1. **Independence across characters.** The joint k′-prior factorises as
   $\prod_i \pi(k'_i \mid \theta)$. Enforced structurally — the loops in
   `R/MkPrimeModel.R:466-520` and `src/mcmc.cpp:260-342` accumulate `lp`
   over $i$ independently.

2. **Support truncation enforced at $k'_i \geq k_{\text{obs},i}$, not at the
   focal-character data.** This respects `feedback_prior_modelling`: the
   prior on $k'_i$ may depend on the *observed support size* $k_{\text{obs},i}$
   of the focal character (a count, not the alignment column itself).
   $k_{\text{obs},i}$ is a deterministic function of the data but enters
   the prior only as a support-bound. This is **not** circular in the
   data-truncation sense.

3. **Hyperparameters in valid range.** $0 < p < 1$, $\alpha, \beta > 0$,
   $0 < c < 1$, $a, b > 0$. Boundaries return $-\infty$ at
   `R/MkPrimeModel.R:411-421` and `src/mcmc.cpp:228-237`.

4. **Empirical body $P_{\text{emp}}$ sums to 1 on its full support.**
   Constructed in `R/data.R:67-115` with normalisation enforced explicitly
   (`R/data.R:84-101`). Verified numerically: the packaged
   `data/empiricalNObs.rda` sums to $0.9999$ (rounding) across body + tail.

5. **The relabelling correction belongs to the log-likelihood, not the
   log-prior.** Confirmed by L4 Wave 1 proof and verified again in §5
   below: `LogPrior` body (lines 393–572) contains no
   `mk_prime_relabel_log*` reference.

## 3. Wave-1 prior-audit (R5-1, R5-2, R5-3) re-verification

These findings come from `project_pi0_audit.md` (2026-05-16) and concern
the *ecology sim3 sampler*, not the k′ prior. They are restated here for
completeness; the re-derivation reuses identities that did not depend on
scoring code.

### R5-1 (HIGH) — Joint identifiability collapse pi0 ↔ z

**Claim.** Under weak per-cell likelihood signal, Gibbs $z$ samples
predominantly from the prior $P(z \mid \pi_0, \theta)$. The Beta(75,25)
prior on $\pi_0$ has effective sample size 100; with 480 z-cells and any
likelihood signal, the prior carries $\leq 30\%$ posterior weight.

**Re-derivation.** Beta($a,b$) prior on $\pi_0$ has prior pseudo-counts
$a+b = 100$. With $n_z = 480$ Bernoulli-like z-indicator cells, the
posterior is approximately Beta($a + n_0, b + n_1$) where $n_0$ is the
count of z=0 cells. The prior's effective sample size as a fraction of
total information is $(a+b)/(a+b+n_z) = 100/580 \approx 17.2\%$, not
$\leq 30\%$ — the 30% figure in the original audit assumed only a fraction
of cells contribute information (consistent with sub-nat likelihood signal
per cell, where each cell's effective $n$ is $\ll 1$). The conclusion
**survives**: prior weight is too small to anchor $\pi_0$ when individual
cells provide weak likelihood signal. No dependence on the scoring code.

**Verdict for R5-1:** **stands**, derivation independent of scoring.

### R5-2 (MED) — $\sigma_\phi = 1$ too tight to allow $\phi = 4$

**Claim.** LogNormal(0, 1) on $\phi$ has prior density at $\phi=1$ that
is $\sim 9 \times$ the density at $\phi=4$.

**Re-derivation.** Density of LogNormal(0, $\sigma$) at $\phi$ is
$f(\phi) = \frac{1}{\phi\sigma\sqrt{2\pi}}\exp(-\frac{(\log\phi)^2}{2\sigma^2})$.
Ratio $f(1)/f(4) = 4 \exp(\frac{(\log 4)^2}{2\sigma^2})$. At $\sigma=1$:
$(\log 4)^2/2 = (1.3863)^2/2 = 0.9609$, so ratio $= 4 e^{0.9609} = 4 \times
2.614 = 10.46$. Approximately 9× (the original figure undercounts by ~15%,
attributable to rounding). Conclusion **stands**: tightness is real.

**Verdict for R5-2:** **stands**, derivation independent of scoring.

### R5-3 (MED) — Beta(1,1) on $\theta$ is uniform; $\theta$ drifts to 0.5

**Claim.** Beta(1,1) gives a uniform prior on $\theta \in (0,1)$, providing
no counterweight against likelihood-driven drift to the $\theta = 0.5$
exchangeability point.

**Re-derivation.** Beta(1,1) density on $(0,1)$ is constant equal to 1,
so log-density contributes 0 to the log-posterior. At $\theta = 0.5$ the
$z=1$ and $z=2$ states become exchangeable under the $\phi$-flipping move,
and the local likelihood Hessian along $\theta$ is zero or singular — a
saddle. Standard label-switching argument. Conclusion **stands**.

**Verdict for R5-3:** **stands**, derivation independent of scoring.

### R5-4 (MED) — phi=0.397 at phi=2 truth, post-relabel invariant violated

**Claim.** The empirical phi-median 0.397 reported at $\phi_{\text{true}}=2$
violates the post-`RelabelEcology` invariant $\phi \geq 1$. Hypothesis was
that the reported value was pre-relabel.

**Resolution.** This is an *empirical claim* — its truth depends on which
diagnostic numbers were used (pre- vs post-relabel) and ultimately on the
output of `sim3-scoring.R::HasBipartSplits` (now patched per
`project_scoring_bug.md`, root-dependence fixed). A math-prover cannot
close this from first principles — it requires:

1. Re-running sim3 against the corrected `HasBipartSplits` (root-independent),
2. Re-extracting post-`RelabelEcology` $\phi$ traces from the corrected
   output,
3. Recomputing the phi-median.

No new derivation is available without that empirical input. The original
audit's update note (2026-05-19, in `project_pi0_audit.md`) already flags
this with "the headline 'phi collapses 5× below truth' may dissolve". No
math identity is at stake.

**Verdict for R5-4:** **deferred-out-of-lane** — empirical re-run needed,
not analytic re-derivation. Math-prover cannot resolve.

## 4. Per-prior derivation, support, and normalisation

For each prior I (i) state the family and parameters; (ii) write the
canonical pmf; (iii) compute the truncated normaliser
$Z_i(\theta) = \sum_{k \geq k_{\text{obs},i}} P(k \mid \theta)$;
(iv) give $\mathbb{E}[k' \mid k' \geq k_{\text{obs},i}]$ and
$\mathbb{E}[u_i \mid u_i \geq 0]$ in closed form where available; and
(v) tabulate for the values used in the prior-validation campaign.

### 4.1 Hierarchical Geometric

**Family.** $u_i \mid p \sim \mathrm{Geo}(p)$ on $u \in \{0, 1, 2, \ldots\}$
with $P(u) = p(1-p)^u$. Then $k'_i = k_{\text{obs},i} + u_i$.

**Support and normaliser.** $\sum_{u=0}^{\infty} p(1-p)^u = p \cdot
\frac{1}{1-(1-p)} = 1$. The prior is normalised on $u \geq 0$ **by
construction**. There is no per-character truncation — $k'_i \geq
k_{\text{obs},i}$ is built into the parameterisation as $u_i \geq 0$.
$Z_i \equiv 1$ for every character and every $p$.

**Conditional expectations.** Standard geometric mean:
$\mathbb{E}[u_i] = (1-p)/p$, $\mathbb{E}[k'_i \mid k' \geq k_{\text{obs},i}]
= k_{\text{obs},i} + (1-p)/p$. **Invariant in $k_{\text{obs}}$** — every
character has the same $\mathbb{E}[u]$.

**Table** (campaign $p$ values; $a = b = 1$ Beta hyperprior used by default):

| $p$ | $\mathbb{E}[u]$ | $k_{\text{obs}}=2$ | $k_{\text{obs}}=3$ | $k_{\text{obs}}=5$ | $k_{\text{obs}}=7$ |
|-----|------|--------|--------|--------|--------|
| 0.30 | 2.333 | 4.333 | 5.333 | 7.333 | 9.333 |
| 0.50 | 1.000 | 3.000 | 4.000 | 6.000 | 8.000 |
| 0.70 | 0.429 | 2.429 | 3.429 | 5.429 | 7.429 |
| 0.90 | 0.111 | 2.111 | 3.111 | 5.111 | 7.111 |

Verified numerically (`Rscript`, §4.5).

### 4.2 Hierarchical Beta-Geometric

**Family.** $u_i \mid \alpha, \beta \sim \mathrm{BetaGeo}(\alpha, \beta)$
on $u \in \{0, 1, 2, \ldots\}$. The pmf marginalised over $p \sim
\mathrm{Beta}(\alpha, \beta)$ is
$$
P(u \mid \alpha, \beta) = \int_0^1 p(1-p)^u \cdot \frac{p^{\alpha-1}(1-p)^{\beta-1}}{B(\alpha,\beta)}\, dp
= \frac{B(\alpha+1, \beta+u)}{B(\alpha,\beta)}.
$$

**Support and normaliser.** Direct check:
$\sum_{u=0}^{\infty} B(\alpha+1, \beta+u)/B(\alpha,\beta) = 1$ (Yule–Simon
identity; sum verified numerically to $1 \pm 3 \times 10^{-5}$ for
$(\alpha,\beta) \in \{(2,2),(3,1),(1.5,0.7)\}$ truncated at $u=1000$).
The prior is normalised on $u \geq 0$ by construction. $Z_i \equiv 1$.

**Conditional expectations.** $\mathbb{E}[u] = \beta/(\alpha-1)$ for
$\alpha > 1$, divergent otherwise (heavy tail).
$\mathbb{E}[k'_i] = k_{\text{obs},i} + \beta/(\alpha-1)$, again invariant
in $k_{\text{obs}}$.

| $(\alpha, \beta)$ | $\mathbb{E}[u]$ | comment |
|--------|--------|--------|
| $(1, 1)$ | $\infty$ | improper mean; Hyperprior Exp(1) on $\alpha,\beta$ |
| $(2, 2)$ | $2.000$ | |
| $(3, 2)$ | $1.000$ | |
| $(0.5, 0.5)$ | $\infty$ | |

### 4.3 Empirical-Geometric (convolution; **EG-001 truncation bug**)

**Family.** $k'_i = N_{\text{obs},i} + N_{\text{unobs},i}$ with
$N_{\text{obs},i} \sim P_{\text{emp}}$ on $\{2, 3, \ldots\}$,
$N_{\text{unobs},i} \mid p \sim \mathrm{Geo}(p)$ on $\{0, 1, \ldots\}$.
The marginal pmf is
$$
P(k'_i = m \mid p) = \sum_{j=2}^{m} P_{\text{emp}}(j) \cdot p\,(1-p)^{m-j}, \qquad m \geq 2.
$$

**Total mass.** Since $P_{\text{emp}}$ is a pmf on $j \geq 2$ summing to 1,
and Geo($p$) is a pmf on $u \geq 0$ summing to 1, the convolution
$P(k' = m)$ sums to 1 over $m \geq 2$.

Verified numerically: $\sum_{m \geq 2} P(k'=m \mid p) = 1.00000 \pm 10^{-6}$
for $p \in \{0.3, 0.5, 0.7, 0.9\}$.

**Per-character truncation.** `LogPrior` enforces $k'_i \geq k_{\text{obs},i}$
by returning $-\infty$ if violated, but the log-density returned for valid
states is the **untruncated** convolution log-pmf $\log P(k'_i = m \mid p)$.
The correctly truncated prior has density
$\tilde{P}(k'_i = m) = P(m)/Z_i(p)$ where
$$
Z_i(p) := \sum_{k \geq k_{\text{obs},i}} P(k \mid p).
$$

**This $Z_i$ is missing from the implementation** at both
`R/MkPrimeModel.R:350-379` (`.LogPriorEmpiricalGeometric`) and
`src/mcmc.cpp:285-329`. The log-prior is therefore biased by
$+\sum_i \log Z_i(p)$ relative to the correct truncated form, which
matters because $\log Z_i(p)$ depends on $p$ — so the bias is **not
absorbed by an MH ratio** that varies $p$. (Bias *is* absorbed for moves
that hold $p$ fixed and only change $k'_i$ at a single character with
its own fixed $k_{\text{obs},i}$ — these are the Gibbs $k'$ sweeps and
they remain valid; only MH moves on $p$ are biased.)

**Computed $Z_i(p)$ and conditional expectations** (campaign $p$ values
× campaign $k_{\text{obs}}$ values, using the packaged `empiricalNObs`):

| $p$ | $k_{\text{obs}}$ | $\log Z_i(p)$ | $\mathbb{E}[k' \mid k' \geq k_{\text{obs}}]$ | $\mathbb{E}[u]$ |
|---|---|---|---|---|
| 0.30 | 2 | 0.0000 | 4.864 | 2.864 |
| 0.30 | 3 | −0.2249 | 5.587 | 2.587 |
| 0.30 | 5 | −0.8458 | 7.431 | 2.431 |
| 0.30 | 7 | −1.5254 | 9.386 | 2.386 |
| 0.50 | 2 | 0.0000 | 3.531 | 1.531 |
| 0.50 | 3 | −0.4091 | 4.305 | 1.305 |
| 0.50 | 5 | −1.5237 | 6.192 | 1.192 |
| 0.50 | 7 | −2.7422 | 8.178 | 1.178 |
| 0.70 | 2 | 0.0000 | 2.959 | 0.959 |
| 0.70 | 3 | −0.6349 | 3.810 | 0.810 |
| 0.70 | 5 | −2.2599 | 5.828 | 0.828 |
| 0.70 | 7 | −3.8915 | 7.989 | 0.989 |
| 0.90 | 2 | 0.0000 | 2.642 | 0.642 |
| 0.90 | 3 | −0.9271 | 3.622 | 0.622 |
| 0.90 | 5 | −2.8762 | 5.851 | 0.851 |
| 0.90 | 7 | −4.5103 | 8.151 | 1.151 |

The maximum $|\log Z_i|$ observed is $4.51$ nats at $(p=0.9, k_{\text{obs}}=7)$
— consistent with EG-001's original "−4.7 nats at $k_{\text{obs}}=8$" claim.

**Note on EG-002** (sibling finding, OPEN). `.LogPriorEmpiricalGeometric`
sums $j = 2, \ldots, m$ for every character. This is mathematically
correct for the convolution definition — the sum is over all
$(N_{\text{obs}}, N_{\text{unobs}})$ decompositions consistent with
$k' = m$ — and is not by itself an arithmetic error. However the
truncation $k' \geq k_{\text{obs},i}$ does not back-propagate to restrict
$N_{\text{obs}}$ to $\geq k_{\text{obs},i}$. The intended generative model is
*ambiguous* about whether $k_{\text{obs},i}$ should constrain $N_{\text{obs}}$
or $k'_i$ alone — both are reasonable choices. EG-002 is therefore a
**modelling-choice issue**, not an arithmetic bug. The current
implementation matches choice (b) "truncate at $k'_i$ only".

### 4.4 Logseries (**new finding LS-001**)

**Family.** $k'_i \mid c \sim \mathrm{LogSeries}(c)$ on $k \in \{1, 2, \ldots\}$
with $P(k \mid c) = -\frac{c^k}{k\,\log(1-c)}$. The constant $c \in (0,1)$
is fixed (`kprimeLogseriesC`, default $0.7$) — not estimated.

**Untruncated normaliser.** $\sum_{k=1}^{\infty} -\frac{c^k}{k\,\log(1-c)} = 1$
(Mercator series identity: $\sum_{k\geq 1} c^k/k = -\log(1-c)$).

**Per-character truncation.** Same structure as §4.3: the implementation
at `R/MkPrimeModel.R:516-519` and `src/mcmc.cpp:262-271` computes
$$
\log P(k_i \mid c) = k_i \log c - \log k_i - \log(-\log(1-c))
$$
— the **untruncated** logseries log-pmf. The truncation $k'_i \geq
k_{\text{obs},i}$ is enforced only by a $-\infty$ guard (line 408), and the
truncation normaliser $Z_i(c) = \sum_{k \geq k_{\text{obs},i}} P(k \mid c)$
is **missing**.

**Effect.** Because $c$ is *fixed* in the current implementation, the
missing $\log Z_i(c)$ is a per-character additive constant — it shifts the
log-prior by a $k_{\text{obs},i}$-dependent amount but cancels in all MH
ratios that hold $c$ fixed. **In the current code path (fixed $c$),
LS-001 is therefore latent — the bias is real but inert.**

**Computed $\log Z_i(c)$** for representative $(c, k_{\text{obs}})$:

| $c$ | $k_{\text{obs}}$ | $\log Z_i(c)$ | $\mathbb{E}[k' \mid k' \geq k_{\text{obs}}]$ | $\mathbb{E}[u]$ |
|---|---|---|---|---|
| 0.50 | 2 | −1.278 | 2.589 | 0.589 |
| 0.50 | 3 | −2.320 | 3.669 | 0.669 |
| 0.50 | 5 | −4.157 | 5.757 | 0.757 |
| 0.50 | 7 | −5.847 | 7.807 | 0.807 |
| 0.70 | 2 | −0.871 | 3.241 | 1.241 |
| 0.70 | 3 | −1.537 | 4.415 | 1.415 |
| 0.70 | 5 | −2.655 | 6.621 | 1.621 |
| 0.70 | 7 | −3.647 | 8.745 | 1.745 |
| 0.90 | 2 | −0.496 | 5.775 | 3.775 |
| 0.90 | 3 | −0.836 | 7.308 | 4.308 |
| 0.90 | 5 | −1.361 | 9.999 | 4.999 |
| 0.90 | 7 | −1.791 | 12.459 | 5.459 |

**Conditional mean.**
$\mathbb{E}[k' \mid k' \geq k_{\text{obs}}] = \frac{1}{Z(c) \cdot (-\log(1-c))} \sum_{k \geq k_{\text{obs}}} c^k$
(no closed form simpler than the tail incomplete-series, but computable).

**Why this should still be filed:** if a future change estimates $c$
(e.g. a Beta hyperprior on $c$, analogous to `kprimeHyperA/B` for $p$),
the missing $Z_i(c)$ becomes an active bias of $-\sum_i \log Z_i(c)$ in
the log-posterior for $c$. The bug must be fixed *before* $c$ is moved
from a fixed hyperparameter to a sampled parameter.

### 4.5 Numerical verification

```r
# Geometric: closed-form E[u] = (1-p)/p — exact.
# Beta-Geometric: numerical sum of B(a+1, b+u)/B(a,b) over u=0..1000;
#   matches B(a+1,b)/B(a,b)·sum_geom up to tail truncation. E[u]=b/(a-1).
# Empirical-Geometric: P(k'=m) computed by convolution; Z_i and E[u]
#   summed over m=2..50 (mass beyond m=50 < 1e-8 for all p ≥ 0.3).
# Logseries: P(k;c) = -c^k/(k log(1-c)) summed over k=kObs..5000.
```

All values in the tables above match direct numerical sums to 4+ decimal
places.

## 5. Implementation cross-check

### 5.1 No relabelling-correction leak in `LogPrior`

`R/MkPrimeModel.R::LogPrior` body (lines 393–572): no reference to
`mk_prime_relabel_log`, `lgamma(.+1) - lgamma(.-.+1)`, or any
falling-factorial expression on `kPrime`. The four prior branches
(lines 467–520) reference only the kernels above.

`src/mcmc.cpp::cpp_log_prior` (lines 206–352): same. No
`mk_prime_relabel_log` call; the only `lgamma` is for the Dirichlet
constant on relative branch lengths (line 246, `std::lgamma(relBrLengths.size())`).

This confirms L4 Wave 1's verdict that the correction lives in the
log-likelihood only, never the log-prior. No double-count.

### 5.2 Per-prior cite

| Prior | R implementation | C++ implementation |
|-------|------------------|--------------------|
| Geometric | `R/MkPrimeModel.R:467-477` | `src/mcmc.cpp:332-341` |
| Beta-Geometric | `R/MkPrimeModel.R:500-510` | `src/mcmc.cpp:272-284` |
| Empirical-Geometric | `R/MkPrimeModel.R:478-499` + `R/MkPrimeModel.R:350-379` | `src/mcmc.cpp:285-331` |
| Logseries | `R/MkPrimeModel.R:511-520` | `src/mcmc.cpp:262-271` |

### 5.3 Normaliser handling — match or bug

| Prior | Truncation? | Normaliser $Z_i$ present? | Verdict |
|-------|-------------|---------------------------|---------|
| Geometric | none (param on $u \geq 0$) | trivially 1 by construction | ✅ correct |
| Beta-Geometric | none (param on $u \geq 0$) | trivially 1 by construction | ✅ correct |
| Empirical-Geometric | $k' \geq k_{\text{obs},i}$ on the convolution | **missing** | ❌ **EG-001** (filed) |
| Logseries | $k' \geq k_{\text{obs},i}$ on full $\{1,2,\ldots\}$ support | **missing** | ❌ **LS-001 (new finding)** — currently latent because $c$ is fixed; activates if $c$ ever becomes a sampled parameter |

## 6. Edge cases

### 6.1 $k_{\text{obs},i} = 1$ — no support truncation needed

For Mk-style transformational characters with $\geq 2$ states, the
domain-level minimum is $k_{\text{obs}} \geq 2$. The prior support starts
at $k' \geq 2$ for Empirical-Geometric (where $P_{\text{emp}}$ has zero
mass at $k=1$) and at $k' \geq 1$ for Logseries. At $k_{\text{obs}} = 1$:

- **Geometric / Beta-Geometric**: $u \geq 0$ means $k' \geq 1$, no
  truncation needed; the priors are unaffected.
- **Empirical-Geometric**: support already starts at $k' \geq 2$ (since
  $P_{\text{emp}}(1) = 0$), so $k_{\text{obs}} = 1$ would force $k' \geq 1$
  but the prior has zero mass at $k' = 1$ — this is **silently
  inconsistent**: the support truncation $k' \geq 1$ would allow $k' = 1$
  but the convolution gives $P(k' = 1) = 0$. Implementation behaviour:
  `.LogPriorEmpiricalGeometric` line 363 calls `seq.int(2L, m)` which
  is empty when $m = 1$, yielding an empty `logTerms` and returning
  $-\infty$ (line 372). Consistent: $k' = 1$ is rejected. Acceptable.
- **Logseries**: full support $k \geq 1$. $\log P(k'=1) = \log c -
  0 - \log(-\log(1-c))$, well-defined.

### 6.2 $k_{\text{obs}} > \max(\mathrm{body})$ — Empirical-Geometric tail behaviour

`empiricalNObs` body covers $k = 2, \ldots, 16$ explicitly. For
$k_{\text{obs}} = 20$ (beyond body), `.LogPemp` returns $-\infty$ for
$j < \mathtt{tail\_start\_k} = 17$ unless they fall in the tail. The
tail formula at `R/MkPrimeModel.R:326-330` continues geometrically for
$k \geq 17$. EG-005 (filed, LOW): when $\mathtt{tail\_start\_k} >
\mathrm{nBody} + 2$, a **gap** in support arises — values in the gap
return $-\infty$ and kill the chain.

This is an EG-005 issue, not an EG-001 issue. EG-001's missing $Z_i$
remains the binding finding even with the gap, because the gap is
typically not encountered in campaign data ($k_{\text{obs}} \leq 8$).

### 6.3 $p = 0$ — degenerate Geometric / Empirical-Geometric

$p = 0$ means $P(\mathrm{Geo} = u) = 0 \cdot 1^u = 0$ for all $u$, i.e.
the geometric has no mass. `LogPrior` rejects $p \leq 0$ at line 412 with
$-\infty$, and `cpp_log_prior` at line 229. Correct.

### 6.4 $p = 1$ — point mass at $u = 0$

$p = 1$ means $P(\mathrm{Geo} = 0) = 1$, $P(\mathrm{Geo} > 0) = 0$ — point
mass at $u = 0$. `LogPrior` rejects $p \geq 1$ at line 412. Slightly
conservative (the limit is well-defined and corresponds to $k' \equiv
k_{\text{obs}}$ for all characters) but does not introduce a bug.
Acceptable.

### 6.5 $c = 0$ / $c = 1$ — Logseries

$c \to 0$: $\log(-\log(1-c)) \to -\infty$, density of $k=1$ tends to 1.
$c \to 1$: $-\log(1-c) \to \infty$, all $k$ tend to zero density —
divergent. Both endpoints rejected at `R/MkPrimeModel.R:420` and
`src/mcmc.cpp` boundary checks. Correct.

### 6.6 `tail_decay = 0` in `MkPrimeEmpiricalPrior`

`R/data.R:84-91` requires the body to sum to 1 exactly (within $10^{-8}$).
This is enforced; the body becomes the full pmf and EG reduces to a
finitely-supported convolution. The truncation bug at EG-001 still
applies.

## 7. Verdict

- **Geometric prior:** **Watertight.** Normaliser is 1 by construction;
  implementation matches first-principles derivation.

- **Beta-Geometric prior:** **Watertight.** Normaliser is 1 by
  construction; implementation matches.

- **Empirical-Geometric prior:** **Implementation bug — confirmed EG-001.**
  Missing per-character truncation normaliser $Z_i(p)$ biases the
  log-posterior for $p$ by up to $-4.5$ nats per character (at $p=0.9,
  k_{\text{obs}}=7$). Fix requires per-character $Z_i(p)$ at runtime, not a
  ≤30-line mechanical edit — **no patch attached** (non-trivial fix
  spanning R, C++, and tests).

- **Logseries prior:** **Implementation bug — new finding LS-001.**
  Same family as EG-001: missing per-character $Z_i(c)$. Currently
  latent because `kprimeLogseriesC` is a fixed hyperparameter (not
  sampled); the bias is a $k_{\text{obs},i}$-dependent constant that
  cancels in all MH ratios in the current code path. **Must be fixed
  before $c$ is ever moved to a sampled parameter.** No patch attached
  — fix is parallel to EG-001 and non-trivial.

- **R5-1, R5-2, R5-3** (ecology sim3 sampler, not k′-prior): **stand**.
  Derivations independent of `sim3-scoring.R::HasBipartSplits` and survive
  the root-dependence correction.

- **R5-4** (ecology sim3 phi median): **deferred-out-of-lane.** Empirical
  claim requiring sim3 rerun against the corrected scoring; not closeable
  by a math-prover from first principles.

- **Relabelling correction in log-prior:** **confirmed absent.** No
  falling-factorial leak anywhere in `LogPrior` or `cpp_log_prior`.
  Consistent with L4 Wave 1 verdict.

### Patch status

No patch in `dev/red-team/patches/`. Both EG-001 and LS-001 require:

1. Computing $Z_i(\theta) = \sum_{k \geq k_{\text{obs},i}} P(k \mid \theta)$
   per character at every prior evaluation, ideally cached when
   $k_{\text{obs},i}$ are fixed and only $\theta$ varies.
2. Subtracting $\sum_i \log Z_i(\theta)$ from the log-prior in both
   `R/MkPrimeModel.R::LogPrior` (and `.LogPriorEmpiricalGeometric`) and
   `src/mcmc.cpp::cpp_log_prior`.
3. Numerical care: the EG $Z_i$ involves an outer summation over
   $k = k_{\text{obs}}, \ldots, K_{\max}$ where each term is itself a
   log-sum-exp over the convolution — costs $O(K_{\max}^2)$ unless the
   running geometric-tail trick is used.

This spans R, C++, and likely tests — out of scope for the trivial-fix
policy.
