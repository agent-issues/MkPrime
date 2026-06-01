# Marginal-k Mk′ — Rao–Blackwell consistency of `sampled_k` and `marginal_k` (Stage 2)

> **Lane.** Red-team analytic proof for the `marginal-k` feature on branch
> `feat/marginal-k` at `C:\Users\pjjg18\GitHub\worktrees\mkp\marginal-k`,
> HEAD `644cb4a`. This is the **Stage 2 RB-consistency re-proof** flagged as
> *Caveat (1)* / follow-up *lane (marginal-k RB-consistency) / role math-prover*
> in the sibling normaliser proof
> `dev/red-team/proofs/marginal-k-truncation-normaliser.md` (read that first;
> its §2 closed form for $Z(p)$ and §4.4 Model B normaliser are **reused, not
> re-derived**).
>
> **Claim proved.** For each transformational character $i$, summing
> (`logSumExp`) the *truncated* `sampled_k` per-character **joint**
> $\big(\log\pi(k'_i) + L_i(k'_i)\big)$ over $k'_i\in[k_{\mathrm{obs},i},K]$
> reproduces the `marginal_k` per-character value `charLL`$_i$, for **both**
> Model A (unconditional) and Model B (conditional). Lifting over characters
> (the likelihood factorises and the $k'$-independent prior terms are identical
> in both modes), the two `likelihoodMode`s target the **same** posterior over
> the shared parameters $\vartheta=(\text{tree},\mu,\sigma,p,\dots)$.
>
> **Do not** edit `R/`, `src/`, `tests/`. No patch attached: this proof finds
> **no discrepancy** between the stated formulas and the source — the verdict is
> *Watertight* (the implementation already matches; the existing regression
> test `tests/testthat/test-marginal-k-truncation.R:258–333` is its $n=1$
> empirical anchor and is the routine this proof justifies analytically).

---

## 0. TL;DR

* **Verdict: Watertight** (under the explicit assumptions of §2, chiefly the
  *empirical* likelihood-equivalence A6, and the *no-pruning* idealisation A7).
  Every formula in the brief was checked line-by-line against
  `src/mcmc.cpp` and `src/mcmc_state.h`; **all match**.
* The decisive per-character identity, for $k=k_{\mathrm{obs},i}+c$,
  $$\operatorname*{logSumExp}_{k=k_{\mathrm{obs},i}}^{K}
    \big[\log\pi(k) + L_i(k)\big]
    \;=\;\underbrace{\operatorname*{logSumExp}_{c=0}^{K-k_{\mathrm{obs},i}}
      \big[\mathrm{rawLL}_i(c) + \log p + c\log q\big]}_{\text{evaluator core, }\,\mathtt{logSumExp}(w)}
    \;+\;\mathrm{corr}_i
    \;=\;\mathrm{charLL}_i,$$
  closes because every $k$-**in**dependent term — the per-character constants
  $\mathrm{corr}_i$ ($=(k_{\mathrm{obs},i}-2)\log q-\log Z_A$ for Model A;
  $=-\log Z_{B,i}$ for Model B) — factors **out** of the `logSumExp`, leaving
  exactly the marginal evaluator's core sum.
* Verified to machine precision on synthetic $\mathrm{rawLL}$ for random
  $(p,K,k_{\mathrm{obs}})$, both models, including the $k_{\mathrm{obs}}=2$ and
  $k_{\mathrm{obs}}=K$ edge cases (§Edge cases).
* **Modes vs models.** `sampled_k`↔`marginal_k` agree **within each** of Model A
  and Model B. Model A and Model B do **not** agree with each other: their
  per-character constants differ by a function **of $p$**, so they are genuinely
  different priors / posteriors on $p$ (this is exactly the SBC asymmetry in the
  memory note `project_sbc_kprime_structural`). The "reparameterisation by a
  per-character constant" (point 3) is a statement about the **weights inside
  the sum**, not about the posterior on $p$.

---

## Theorem (Rao–Blackwell consistency of the two likelihood modes)

After the Stage 2 edit (commit `e5ebf90` + `e074109`; truncating and
renormalising the `sampled_k` geometric prior to match the `marginal_k`
evaluator), the `sampled_k` and `marginal_k` likelihood modes target the
**same** posterior over the shared parameters
$\vartheta=(\text{tree},\mu,\sigma,p,\beta_{\mathrm{scale}},\dots)$, for both
Model A (`unconditionalPrior = true`) and Model B (`false`). Concretely, for
each transformational character $i$,
$$\log\!\!\sum_{k'_i=k_{\mathrm{obs},i}}^{K}
   \exp\!\big(\log\pi(k'_i\mid p) + L_i(k'_i)\big)
   \;=\;\mathrm{charLL}_i(p),$$
where the left side is built from the **edited** `cpp_log_prior` plus the
shared per-$k'$ likelihood, and the right side is the per-character value
returned by `cpp_log_likelihood_marginal`. Summing over the $n_{\mathrm{trans}}$
transformational characters and adding the (mode-identical) $k'$-independent
prior terms gives
$$\log\pi_{\mathrm{marginal}}(\vartheta\mid\mathbf y)
  \;=\;\log\!\!\sum_{k'_{1:n}}\exp\!\big(\log\pi_{\mathrm{sampled}}(\vartheta,k'_{1:n}\mid\mathbf y)\big)
  \;+\;\text{const}(\mathbf y),$$
i.e. `marginal_k` is the exact marginalisation of the `sampled_k` joint over
the latent state counts $k'_{1:n}$.

---

## Assumptions

1. **Parameters.** $p\in(0,1)$, $q:=1-p\in(0,1)$. The `cpp_log_prior` boundary
   guard rejects $p\le 0$ and $p\ge 1$ for the geometric arm
   (`src/mcmc.cpp:295–296`), so $q\in(0,1)$ strictly throughout.
2. **Truncation cap.** $K:=$ `data.kprimeTruncK` is a single declared integer
   constant, **identical** between the two modes for a given model object
   (it is the same field `McmcData::kprimeTruncK`, default `30` at
   `src/mcmc_state.h:160`, wired from `MkPrimeModel(kprimeTruncK=…)` whose
   package default is 200). Require $K\ge k_{\mathrm{obs},i}$ for every $i$
   (else that character has empty support; see Edge cases).
3. **Observed state count.** $k_{\mathrm{obs},i}\ge 2$ is the number of distinct
   realised states in column $i$ (`data.kObs[gi]`), fixed by the data (not by
   any parameter).
4. **Per-character likelihood.** $L_i(k)$ is the per-character log-likelihood
   when character $i$ has exactly $k$ states. It is a well-defined function of
   $k$ for $k\ge k_{\mathrm{obs},i}$, and $L_i(k)=-\infty$ for
   $k<k_{\mathrm{obs},i}$ (the relabelling falling-factorial
   $R(k,k_{\mathrm{obs},i})$ vanishes — `relabelling-correction.md`). $L_i$ may
   embed the Lewis-Mkv ascertainment $-\log a_k$ and the relabelling $\log R$;
   the proof treats $L_i(\cdot)$ as a black box and never needs its internals.
5. **Conditional independence / factorisation.** Given $\vartheta$ and the
   $k'_{1:n}$, the characters are conditionally independent, so the total
   log-likelihood is $\sum_i L_i(k'_i)$ and the latent counts $k'_{1:n}$ are a
   priori independent across characters given $p$. (Standard for site/character
   models; Felsenstein 2004 ch. 16; matches the per-partition independent
   pruning in `cpp_partition_log_likelihood`.)
6. **(EMPIRICAL — load-bearing.) Likelihood-equivalence.**
   $$L_i(k_{\mathrm{obs},i}+c)\;=\;\mathrm{rawLL}_i(c)\qquad(c=0,1,2,\dots),$$
   where $\mathrm{rawLL}_i(c)$ is the marginal evaluator's stored raw
   log-likelihood (`state.charLLCache[ti*kMaxKprimeCand+c]` =
   `charLogW − logPriorByU[c]`, `src/mcmc.cpp:4348–4351`), i.e. the
   `sampled_k` per-character likelihood at fixed $k'_i=k_{\mathrm{obs},i}+c$
   equals the marginal evaluator's per-candidate likelihood. This is **not
   proved here**; it is established empirically by the **C-i guard**
   (`tests/testthat/test-marginal-k-truncation.R:107–152`), which evaluates the
   `sampled_k` likelihood at each pinned $k'$ via `eval_full_loglik_cpp(modS)`
   and matches it to the marginal evaluator's full reference at low $p$. Both
   paths call the same per-$k'$ pruning machinery
   (`compute_per_kprime_log_lik`; `compute_full_loglik` →
   `cpp_log_likelihood_marginal` vs `cpp_log_likelihood` at
   `src/mcmc.cpp:1056–1075`), which is *why* the equivalence holds, but the
   guarantee rests on the test, not on a separate analytic argument.
   > The C-i guard runs under **Model A only**; the per-$k'$ likelihood is
   > **prior-variant-independent** (the prior never enters the pruning kernel),
   > so A6 carries unchanged to Model B. The Stage-2 test
   > (`:258–333`) additionally exercises both variants end-to-end.
7. **(IDEALISATION.) No pruning of the candidate sum.** The marginal
   evaluator's `logSumExp` runs over $c=0,\dots,n_{\mathrm{Eff}}-1$ with
   $n_{\mathrm{Eff}}=\min(\mathtt{nCand},\,K-k_{\mathrm{obs},i}+1)$
   (`src/mcmc.cpp:4341`); $\mathtt{nCand}$ may be reduced below
   $K-k_{\mathrm{obs},i}+1$ by the M-164 log-cutoff
   $\mathtt{kKprimeLogCutoff}=-25$ and the absolute cap
   $\mathtt{kMaxKprimeCand}=256$ (`src/mcmc_state.h:397–398`). The analytic
   identity is stated for the **full** sum $c=0,\dots,K-k_{\mathrm{obs},i}$
   (i.e. $n_{\mathrm{Eff}}=K-k_{\mathrm{obs},i}+1$). Pruning is a separate,
   bounded numerical approximation: the omitted tail mass is
   $<(K-k_{\mathrm{obs},i})\,e^{-25}\approx 4\times10^{-10}$ in absolute
   likelihood-weight terms (see Edge cases / §Pruning). The `sampled_k`
   explicit sum has **no** cutoff, so the only non-bit-exactness between the two
   modes is this pruned tail.
8. **Shared-parameter prior block is mode-identical.** Every $k'$-independent
   prior term — tree-length $\mathrm{Gamma}$, Dirichlet branch-proportions,
   $\sigma=$`rate_log_sd` $\mathrm{Gamma}$, neomorphic rate log-normals,
   $\beta_{\mathrm{scale}}$, **and the $p$-hyperprior $\mathrm{dbeta}(p)$** — is
   added by the *same* `cpp_log_prior` in *both* modes, with the **same** value
   (verified in §Implementation cross-check, point 8).

---

## Statement (symbols)

| Symbol | Meaning | Code name / site |
|---|---|---|
| $p,\;q=1-p$ | geometric success prob. and complement | `state.p`; `logP=log(p)`, `log1mP=log1p(-p)` |
| $K$ | truncation cap on $k'$ | `data.kprimeTruncK` (`mcmc_state.h:160`) |
| $k_{\mathrm{obs},i}\ge2$ | observed distinct states, char $i$ | `data.kObs[gi]` |
| $k'_i\in[k_{\mathrm{obs},i},K]$ | latent state count | `kPrime[gi]` (sampled) / summed (marginal) |
| $c=k'_i-k_{\mathrm{obs},i}\ge0$ | offset (sum index) | candidate index `c` / `ko` |
| $L_i(k)$ | per-char log-lik at $k$ states ($-\infty$ if $k<k_{\mathrm{obs},i}$) | pruning kernel |
| $\mathrm{rawLL}_i(c)$ | evaluator raw LL at $k=k_{\mathrm{obs},i}+c$ | `charLLCache[ti,c]` (`mcmc.cpp:4351`) |
| $\log\pi_A,\log\pi_B$ | truncated sampled-$k$ log-pmf, Model A / B | `cpp_log_prior` (`mcmc.cpp:428–445`) |
| $Z_A,\;Z_{B,i}$ | truncation normalisers | `logZA` (`:431/4296`), `logZB` (`:443/4371`) |
| $w_{i,c}$ | evaluator weight $=\mathrm{rawLL}_i(c)+\log p+c\log q$ | `charLogW[ti,c]` (`mcmc.cpp:4345`) |
| $\mathrm{corr}_i$ | per-char additive constant | `:4369` (A) / `:4371` (B) |

**Truncated sampled-$k$ log-pmf (the Stage 2 edit).** For
$k\in[k_{\mathrm{obs},i},K]$ (the likelihood support; see point 2 of §Proof for
the Model A nuance below $k_{\mathrm{obs},i}$):

$$\log\pi_A(k\mid p) = \log p + (k-2)\log q - \log Z_A,\qquad
  \log Z_A=\log\!\big(1-q^{\,K-1}\big),\quad k\in[2,K]; \tag{A}$$
$$\log\pi_B(k\mid p) = \log p + (k-k_{\mathrm{obs},i})\log q - \log Z_{B,i},\qquad
  \log Z_{B,i}=\log\!\big(1-q^{\,K-k_{\mathrm{obs},i}+1}\big),\quad k\in[k_{\mathrm{obs},i},K]. \tag{B}$$

**Marginal-$k$ per-character value.** With $\mathrm{rawLL}_i(c)$ the per-$k'$
likelihood at $k=k_{\mathrm{obs},i}+c$ and (idealising pruning, A7) the full
range $c=0,\dots,K-k_{\mathrm{obs},i}$,

$$\mathrm{charLL}_i \;=\;
  \operatorname*{logSumExp}_{c=0}^{K-k_{\mathrm{obs},i}}
   \Big[\mathrm{rawLL}_i(c) + \underbrace{\log p + c\log q}_{\mathtt{logPriorByU}[c]}\Big]
   \;+\;\mathrm{corr}_i, \tag{M}$$
$$\mathrm{corr}_i = \begin{cases}
   (k_{\mathrm{obs},i}-2)\log q - \log Z_A, & \text{Model A (\texttt{mcmc.cpp:4369/4404})},\\[4pt]
   -\log Z_{B,i}, & \text{Model B (\texttt{mcmc.cpp:4371/4406})}.
  \end{cases}$$

---

## Proof

### Part 1 — Prior normalisation (point 1)

Reuse Lemma §2 of `marginal-k-truncation-normaliser.md`: for $j$ a counting
index and $q\in(0,1)$, $\sum_{j=0}^{m}q^{j}=\frac{1-q^{m+1}}{1-q}$.

**Model A.** The unnormalised weights are $w_k = p\,q^{k-2}$ on $k\in\{2,\dots,K\}$.
Substituting $j=k-2$ ($j=0,\dots,K-2$):
$$Z_A=\sum_{k=2}^{K}p\,q^{k-2}=p\sum_{j=0}^{K-2}q^{j}
   =p\cdot\frac{1-q^{K-1}}{1-q}=p\cdot\frac{1-q^{K-1}}{p}=1-q^{K-1},$$
so $\log Z_A=\log(1-q^{K-1})$, and (A) is a genuine normalised pmf on
$\{2,\dots,K\}$: $\sum_{k=2}^K \pi_A(k)=Z_A/Z_A=1$. $\square$

**Model B.** Weights $w_k=p\,q^{k-k_{\mathrm{obs},i}}$ on
$k\in\{k_{\mathrm{obs},i},\dots,K\}$. Substituting $u=k-k_{\mathrm{obs},i}$
($u=0,\dots,K-k_{\mathrm{obs},i}$):
$$Z_{B,i}=\sum_{u=0}^{K-k_{\mathrm{obs},i}}p\,q^{u}
   =p\cdot\frac{1-q^{K-k_{\mathrm{obs},i}+1}}{p}=1-q^{\,K-k_{\mathrm{obs},i}+1},$$
so $\log Z_{B,i}=\log(1-q^{K-k_{\mathrm{obs},i}+1})$, and (B) normalises on
$\{k_{\mathrm{obs},i},\dots,K\}$. $\square$

**Boundary / limit behaviour.**
* $K\to\infty$: $q^{K-1}\to0$, so $Z_A\to1$, $Z_{B,i}\to1$ — the untruncated
  geometric needs no normaliser (the $K=\infty$ model of the pre-Stage-1
  `sampled_k`); both modes reduce to it.
* $p\to1^-$ ($q\to0^+$): $Z_A\to1-0=1$, $Z_{B,i}\to1$; all mass on $k=2$ (A) or
  $k=k_{\mathrm{obs},i}$ (B). Well-defined.
* $p\to0^+$ ($q\to1^-$): $Z_A=1-q^{K-1}\to0$ and $\log Z_A\to-\infty$ (linearised
  $Z_A\approx(K-1)p$). The renormalisation $-\log Z_A\to+\infty$ per character;
  the posterior stays proper because the $\mathrm{dbeta}(p)$ hyperprior and the
  likelihood product dominate, and `cpp_log_prior` rejects $p=0$ exactly.
* $K=2$: $Z_A=1-q^{1}=p$, single state $k=2$ carries all mass
  ($\pi_A(2)=p/p=1$). Requires $k_{\mathrm{obs},i}\le2$, i.e. $k_{\mathrm{obs},i}=2$.
* **Empty support** ($K<k_{\mathrm{obs},i}$): the sum has no terms; the code
  returns $-\infty$ for that character (see Edge cases).

These match the closed forms used at `src/mcmc.cpp:431` (`logZA`) and
`:443/4371` (`logZB`), evaluated stably as
`log1p(-exp((K-1)*log1p(-p)))` and `log1p(-exp((K-kObs+1)*log1p(-p)))`.

### Part 2 — The per-character RB identity, both models (point 2)

Fix character $i$ and write $k=k_{\mathrm{obs},i}+c$, so $c$ ranges over
$0,\dots,K-k_{\mathrm{obs},i}$ as $k$ ranges over $k_{\mathrm{obs},i},\dots,K$
(a bijection). The `sampled_k` per-character **joint** at fixed $k'_i=k$ is
$$J_i(k):=\log\pi(k\mid p)+L_i(k).$$
We compute $\operatorname*{logSumExp}_{k} J_i(k)$ and show it equals (M).

Throughout we use the **shift-invariance of `logSumExp`**: for any constant $a$
(independent of the summation index),
$$\operatorname*{logSumExp}_{c}\big[x_c+a\big]
  =\log\sum_c e^{x_c+a}=\log\Big(e^{a}\sum_c e^{x_c}\Big)
  =a+\operatorname*{logSumExp}_{c}[x_c]. \tag{$\star$}$$

**Model B.** Insert (B) and A6 ($L_i(k)=\mathrm{rawLL}_i(c)$):
$$
J_i^B(k)=\Big[\log p + c\log q - \log Z_{B,i}\Big] + \mathrm{rawLL}_i(c)
        =\underbrace{\big[\mathrm{rawLL}_i(c)+\log p + c\log q\big]}_{w_{i,c}}
          \;-\;\log Z_{B,i}.
$$
The term $-\log Z_{B,i}$ is **independent of $c$** (it depends only on $p,K,
k_{\mathrm{obs},i}$). By ($\star$) with $a=-\log Z_{B,i}$,
$$
\operatorname*{logSumExp}_{k=k_{\mathrm{obs},i}}^{K} J_i^B(k)
  =\operatorname*{logSumExp}_{c=0}^{K-k_{\mathrm{obs},i}} w_{i,c}\;-\;\log Z_{B,i}
  =\mathrm{charLL}_i\quad(\text{Model B, since } \mathrm{corr}_i=-\log Z_{B,i}).
$$
$\blacksquare$ (Model B)

**Model A.** Insert (A). Note the exponent is $k-2$, *not* $c=k-k_{\mathrm{obs},i}$;
split it as $k-2=(k_{\mathrm{obs},i}-2)+c$:
$$
J_i^A(k)=\Big[\log p + (k-2)\log q - \log Z_A\Big] + \mathrm{rawLL}_i(c)
        =\underbrace{\big[\mathrm{rawLL}_i(c)+\log p + c\log q\big]}_{w_{i,c}}
          \;+\;\underbrace{(k_{\mathrm{obs},i}-2)\log q-\log Z_A}_{=\,\mathrm{corr}_i,\ \text{const in }c}.
$$
The bracketed constant is independent of $c$ (it depends only on
$p,K,k_{\mathrm{obs},i}$). By ($\star$) with $a=(k_{\mathrm{obs},i}-2)\log q-\log Z_A$,
$$
\operatorname*{logSumExp}_{k=k_{\mathrm{obs},i}}^{K} J_i^A(k)
  =\operatorname*{logSumExp}_{c=0}^{K-k_{\mathrm{obs},i}} w_{i,c}
   +(k_{\mathrm{obs},i}-2)\log q-\log Z_A
  =\mathrm{charLL}_i\quad(\text{Model A}).
$$
$\blacksquare$ (Model A)

In both cases the **only** thing that happened is: the per-character constant
($\mathrm{corr}_i$) pulled out of the `logSumExp` via ($\star$), leaving the
evaluator's core sum $\operatorname*{logSumExp}_c w_{i,c}$ — which is *literally*
what `cpp_log_likelihood_marginal` computes as
`charLL = mx + log(s)` over `w = charLogW[ti,c] = rawLL + logPriorByU[c]`
(`src/mcmc.cpp:4362–4367`), before adding `corr_i` at `:4369/4371`. The Model A
and Model B results coincide with the evaluator's **two branches** because the
two `corr_i` expressions are *identical* to the two branch additions in the
source.

**The Model A sub-$k_{\mathrm{obs}}$ nuance (point 2, "unconditional"
semantics).** Model A's prior pmf (A) has support $k\in[2,K]$, i.e. it places
prior mass on $k\in[2,k_{\mathrm{obs},i})$ too. But the **likelihood** kills
those states: $L_i(k)=-\infty$ for $k<k_{\mathrm{obs},i}$ (A4), so the joint
$J_i^A(k)=-\infty$ there and they contribute nothing to the `logSumExp`
($e^{-\infty}=0$). Equivalently: in the marginal evaluator the sum starts at
$c=0$ ($k=k_{\mathrm{obs},i}$), and the only trace of the $[2,k_{\mathrm{obs},i})$
prior mass is *inside the normaliser* $Z_A$ (whose support is $[2,K]$) and in the
shift $(k_{\mathrm{obs},i}-2)\log q$ that converts the Model-B-style weight base
$k_{\mathrm{obs},i}$ to the Model-A base $2$. This is the **intended generative
("unconditional") semantics**: $k'_i$ is drawn from a $\mathrm{TruncGeom}_{[2,K]}(p)$
*independent of the data*, then the data are generated and the column happens to
realise $k_{\mathrm{obs},i}\le k'_i$ distinct states. The prior mass on
$k<k_{\mathrm{obs},i}$ is the prior probability of drawing a $k'_i$ too small to
ever produce this column — it is correctly present in the normaliser and
correctly carries zero likelihood. This is **consistent between forward and
inference** (the SBC forward in `marginal-k-truncation-normaliser.md` §4.2 draws
$k'\sim\mathrm{TruncGeom}_{[2,K]}$ and the inference normalises over the same
$[2,K]$), which is what SBC calibration requires (cf. that proof's §5.1
theorem).

> Contrast Model B, whose support is exactly $[k_{\mathrm{obs},i},K]$ — the
> *conditional* ("given the data realised at least $k_{\mathrm{obs},i}$ states")
> semantics. There is no sub-$k_{\mathrm{obs}}$ mass in $Z_{B,i}$.

### Part 3 — Model A vs Model B relationship (point 3)

The marginal evaluator's own comment (`src/mcmc.cpp:4274–4282`) asserts: the
Model A weight equals the Model B weight times the per-character constant
$q^{\,k_{\mathrm{obs},i}-2}$. Check: at $k=k_{\mathrm{obs},i}+c$,
$$
w_k^A = p\,q^{k-2}=p\,q^{(k-k_{\mathrm{obs},i})+(k_{\mathrm{obs},i}-2)}
      = \big(p\,q^{k-k_{\mathrm{obs},i}}\big)\,q^{\,k_{\mathrm{obs},i}-2}
      = w_k^B\;q^{\,k_{\mathrm{obs},i}-2}.\quad\checkmark
$$
In log-space, $\log w_k^A-\log w_k^B=(k_{\mathrm{obs},i}-2)\log q$, a constant
**in the sum index $c$** (and in $k$), so it factors out of the `logSumExp` by
($\star$). This is precisely the $(k_{\mathrm{obs},i}-2)\log q$ that appears in
$\mathrm{corr}_i^A$ but not $\mathrm{corr}_i^B$. The evaluator exploits this by
storing the cache and `logPriorByU` in the **Model B** convention (located at
$k_{\mathrm{obs},i}$, `:4271, 4348`) and adding the constant only for Model A
(`:4369`). Both models therefore satisfy the per-character RB identity (Part 2),
and they are a **reparameterisation differing by a per-character constant** —
*in the summand*.

> **Crucial caveat (do not over-read point 3).** $q^{\,k_{\mathrm{obs},i}-2}$ is
> constant in $c$ but **not constant in $p$**. Hence Model A and Model B do
> **not** target the same posterior on $p$: the difference of their
> per-character contributions is
> $$\mathrm{corr}_i^A-\mathrm{corr}_i^B
>   =(k_{\mathrm{obs},i}-2)\log q-\log Z_A+\log Z_{B,i},$$
> a non-trivial **function of $p$** (and of $k_{\mathrm{obs},i},K$). Summed over
> characters it reweights the $p$-likelihood, so the $p$-posteriors differ.
> This is consistent with the project's SBC finding
> (`project_sbc_kprime_structural`): Model A (unconditional, generative joint)
> can pass SBC on $p$; Model B (conditional on $k_{\mathrm{obs}}$) cannot, being
> a non-generative reparameterisation. **The RB identity is strictly
> *within-model*:** sampled-A ↔ marginal-A, and sampled-B ↔ marginal-B. It is
> *not* a claim that A ↔ B.

### Part 4 — Lift to the total posterior (point 4)

Both modes share the same parameter vector $\vartheta$ and the same
$k'$-independent prior block (A8; verified at point 8 of the cross-check).
Write the `sampled_k` log-joint as
$$
\log\pi_{\mathrm{sampled}}(\vartheta,k'_{1:n}\mid\mathbf y)
  = \underbrace{G(\vartheta)}_{k'\text{-indep. prior}}
  + \sum_{i=1}^{n_{\mathrm{trans}}}\Big[\log\pi(k'_i\mid p)+L_i(k'_i\mid\vartheta)\Big]
  - \log m(\mathbf y),
$$
where $G(\vartheta)$ collects tree-length Gamma, Dirichlet, $\sigma$-Gamma,
neomorphic, $\beta_{\mathrm{scale}}$ **and** $\mathrm{dbeta}(p)$ (all
mode-identical), and $\log m(\mathbf y)$ is the (parameter-free) log evidence.
Marginalising the latent counts:
$$
\sum_{k'_{1:n}}\exp\!\Big(\sum_i\big[\log\pi(k'_i\mid p)+L_i(k'_i)\big]\Big)
  \;\overset{(\dagger)}{=}\;
  \prod_{i=1}^{n_{\mathrm{trans}}}\ \sum_{k'_i=k_{\mathrm{obs},i}}^{K}
     \exp\!\big(\log\pi(k'_i\mid p)+L_i(k'_i)\big),
$$
where $(\dagger)$ is the standard product-of-sums factorisation: the joint
exponent is a **sum over $i$** of terms each depending only on its **own** $k'_i$
(A5, conditional independence), so the multi-index sum factors into a product of
single-index sums. Taking $\log$ and applying the **per-character identity of
Part 2** to each factor,
$$
\log\sum_{k'_{1:n}}\exp\Big(\sum_i\big[\log\pi(k'_i\mid p)+L_i(k'_i)\big]\Big)
  =\sum_{i=1}^{n_{\mathrm{trans}}}\mathrm{charLL}_i(p)
  =\text{(marginal-}k\text{ transformational log-lik)}.
$$
Adding $G(\vartheta)-\log m(\mathbf y)$ to both sides:
$$
\boxed{\;\log\!\!\sum_{k'_{1:n}}\exp\big(\log\pi_{\mathrm{sampled}}(\vartheta,k'_{1:n}\mid\mathbf y)\big)
  \;=\;\log\pi_{\mathrm{marginal}}(\vartheta\mid\mathbf y)\;}
$$
which is exactly "`marginal_k` = the `sampled_k` joint with $k'_{1:n}$
analytically summed out". Therefore the marginal of the `sampled_k` posterior
over $(\vartheta,k'_{1:n})$ onto $\vartheta$ equals the `marginal_k` posterior on
$\vartheta$: **the two modes target the same posterior over the shared
parameters.** $\blacksquare$

> The neomorphic partitions are handled identically in both modes — in
> `marginal_k` they fall straight through to `cpp_partition_log_likelihood`
> with $k=k_{\mathrm{obs}}$ (`src/mcmc.cpp:4252–4259`), and in `sampled_k` the
> chain's `kPrime` for neo chars never moves off $k_{\mathrm{obs}}$ (case 25
> skips them). They contribute the *same* $k'$-independent term to $G(\vartheta)$
> and drop out of the RB sum.

---

## Implementation cross-check

Every formula in the brief was checked against `src/mcmc.cpp` @ `644cb4a`.
**No discrepancies found.**

| Theorem element | Source | Match |
|---|---|---|
| 1. $p\in(0,1)$ guard (geometric) | `mcmc.cpp:295–296` | ✓ rejects $p\le0,\,p\ge1$ |
| 1. $k'\in[k_{\mathrm{obs}},K]$ support (sampled) | `mcmc.cpp:309–314` | ✓ `kPrime<kObs`→`-Inf`; `isPlainGeom && !marginalK && kPrime>kprimeTruncK`→`-Inf` |
| (A) $\log\pi_A=\log p+(k-2)\log q-\log Z_A$ | `mcmc.cpp:431,434` | ✓ `logZA=log1p(-exp((K-1)*log1mP))`; `lp += logP+(kPrime-2)*log1mP-logZA` |
| (B) $\log\pi_B=\log p+u\log q-\log Z_{B,i}$ | `mcmc.cpp:439–443` | ✓ `u=kPrime-kObs`; `lp += logP+u*log1mP-log1p(-exp((K-kObs+1)*log1mP))` |
| Prior added only in sampled mode | `mcmc.cpp:424` | ✓ `if (!data.marginalK)` wraps the k′-prior block |
| (M) core $w_{i,c}=\mathrm{rawLL}_i(c)+\log p+c\log q$ | `mcmc.cpp:4271,4345` | ✓ `logPriorByU[ko]=logP+ko*log1mP`; `w=charLogW[ti,c]` |
| (M) $\mathrm{rawLL}_i(c)=w-\mathtt{logPriorByU}[c]$ (A6 cache) | `mcmc.cpp:4348–4351` | ✓ |
| (M) cap $n_{\mathrm{Eff}}=\min(\mathtt{nCand},K-k_{\mathrm{obs}}+1)$ | `mcmc.cpp:4341` | ✓ sum runs `c<nEff` (`:4363`) |
| (M) core `logSumExp` | `mcmc.cpp:4362–4367` | ✓ `mx+log(s)` |
| $\mathrm{corr}_i^A=(k_{\mathrm{obs}}-2)\log q-\log Z_A$ | `mcmc.cpp:4369` (non-cache), `:4404` (cache) | ✓ identical both branches |
| $\mathrm{corr}_i^B=-\log Z_{B,i}$ | `mcmc.cpp:4371` (non-cache), `:4406` (cache) | ✓ identical both branches |
| A↔B weight ratio $q^{k_{\mathrm{obs}}-2}$ | `mcmc.cpp:4274–4282` (comment) | ✓ algebra confirmed (Part 3) |
| Cache fast-path reuses same `logPriorByU`, inherits cap via `charLLNCand` | `mcmc.cpp:4382,4389,4398–4400` | ✓ `nCand=charLLNCand[ti]` (=`nEff`); `w=rawLL+logPriorByU[c]` |
| $K$ shared field | `mcmc_state.h:160`; `mcmc.cpp:425,4295` | ✓ same `data.kprimeTruncK` both functions |
| 8. dbeta(p) in **both** modes | `mcmc.cpp:447` | ✓ `lp += dbeta(p,...)` is **outside** the `if(!marginalK)` guard (closes at `:446`) |
| 8. tree/Dirichlet/σ priors in both | `mcmc.cpp:320–336` | ✓ precede the `hasTrans` block, unconditional on `marginalK` |
| Dispatch: `eval_log_prior_cpp`→`cpp_log_prior` | `mcmc.cpp:879–887` | ✓ |
| Dispatch: `eval_full_loglik_cpp`→`compute_full_loglik`→(marginalK)`cpp_log_likelihood_marginal` | `mcmc.cpp:1056–1062, 1078–1093` | ✓ |

**Critical detail — the joint, not the likelihood, is what marginalises.** A
reader might worry that the marginal evaluator's `corr_i` is a *likelihood*
correction while the `sampled_k` $-\log Z$ is a *prior* correction, so summing
"prior + likelihood" might double-count. It does not: under `marginal_k`,
`cpp_log_prior` adds **only** `dbeta(p)` for the geometric arm (the
`if(!data.marginalK)` guard at `:424` skips the entire $k'$-prior block), so the
$k'$-prior mass — *including* the normaliser — is carried **solely** inside
`cpp_log_likelihood_marginal`'s `corr_i`. Under `sampled_k`, that same mass is
carried **solely** inside `cpp_log_prior`'s per-character term. The Stage-2 test
compares **(prior + lik)** on each side (`test-marginal-k-truncation.R:297, 301`)
precisely to make this bookkeeping watertight, and it passes — confirming the
$-\log Z$ lives on exactly one side in each mode and the totals agree.

---

## Edge cases

Verified analytically and (where arithmetic) numerically in pure R on synthetic
$\mathrm{rawLL}$ (random $(p,K,k_{\mathrm{obs}})$, both models): the identity
$\operatorname*{logSumExp}_k J_i(k)=\mathrm{charLL}_i$ held to $\le4.4\times10^{-16}$.

1. **$k_{\mathrm{obs},i}=2$ (lowest, heaviest tail).** Then
   $K-k_{\mathrm{obs},i}+1=K-1$, so $Z_{B,i}=1-q^{K-1}=Z_A$, and
   $(k_{\mathrm{obs},i}-2)\log q=0$. Hence $\mathrm{corr}_i^A=\mathrm{corr}_i^B=-\log Z_A$:
   **Model A and Model B coincide exactly at $k_{\mathrm{obs}}=2$.** (Confirmed:
   `logZA=logZB=-0.256100` at $p=0.05,K=30$.) The RB identity holds trivially in
   both. This is the regime the C-i and Stage-2 tests stress most (`kObs2` case),
   correctly.
2. **$k_{\mathrm{obs},i}=K$ (single-point support).** The sum has one term
   ($c=0$, $k=K$). $\operatorname*{logSumExp}$ of one element is the element:
   $\mathrm{charLL}_i=\mathrm{rawLL}_i(0)+\log p+0+\mathrm{corr}_i$. Model B:
   $Z_{B,i}=1-q^{1}=p$, $\mathrm{corr}_i^B=-\log p$, so
   $\mathrm{charLL}_i^B=\mathrm{rawLL}_i(0)$ — the bare likelihood at $k=K$ (the
   prior is a point mass, log-pmf $=0$). Model A:
   $\mathrm{charLL}_i^A=\mathrm{rawLL}_i(0)+\log p+(K-2)\log q-\log Z_A=\mathrm{rawLL}_i(0)+\log\pi_A(K)$,
   i.e. likelihood + the (non-degenerate) Model A log-pmf at $k=K$. Both equal
   $\operatorname*{logSumExp}_{k=K}J_i(k)=J_i(K)$. ✓ (Confirmed numerically.)
3. **$K<k_{\mathrm{obs},i}$ (empty support).** No valid $k'$. Marginal:
   `nEff = min(nCand, K-kObs+1) <= 0` → `charLLNCand=0`, `mx` stays $-\infty$, the
   character contributes $-\infty$ and `totalLL=R_NegInf` (`mcmc.cpp:4341–4361`).
   Sampled: every pinned $k'\ge k_{\mathrm{obs}}>K$ is rejected by the support
   guard `kPrime>kprimeTruncK`→`-Inf` (`:312`), and $k'<k_{\mathrm{obs}}$ by
   `:311`. Both sides are $-\infty$ → identity holds vacuously ($-\infty=-\infty$).
   *Assumption 2 ($K\ge\max_i k_{\mathrm{obs},i}$) should be enforced at setup*
   (the `mcmc_state.h:158–159` comment says so); if violated, the model is simply
   unsupportable for that character, identically in both modes — a coherent (if
   useless) state, not an inconsistency.
4. **Equal rates / $\sigma=0$.** `rate_log_sd=0` collapses ACRV to a single rate
   (`mcmc.cpp:4309–4311`, `acrvRates=1`). This changes $\mathrm{rawLL}_i(c)$
   identically in both modes (it is a likelihood input, orthogonal to the $k'$
   weighting) and does not touch the RB algebra. The tests run at
   `rate_log_sd=0`. ✓
5. **$p\to0^+$.** $Z_A,Z_{B,i}\to0$; $-\log Z\to+\infty$ per character on **both**
   sides (same closed form), so they track exactly. The heavy geometric tail
   makes the high-$c$ candidates non-negligible — this is exactly why the cap and
   the normaliser must both be present and must match between modes (the
   `marginal-k-truncation-normaliser.md` §5.3 "cap + Z, not Z alone" point). The
   Stage-2 test runs $p\in\{0.02,0.08\}$. ✓

---

## Pruning (Assumption 7) — the one non-bit-exact spot

The identity is **exact in exact arithmetic when the marginal sum reaches $K$**,
i.e. $n_{\mathrm{Eff}}=K-k_{\mathrm{obs},i}+1$. Two devices can shorten the
marginal sum but not the `sampled_k` explicit sum:

* **Absolute cap $\mathtt{kMaxKprimeCand}=256$** (`mcmc_state.h:397`). Since the
  package default is $K=200<256$ and tests pin $K\in\{30,200\}$, the cap never
  binds for $K\le256$; the Stage-1b test (`:172–229`) explicitly confirms the
  numerator reaches the **full** $[2,200]$ support at $K=200$. For $K>256$ the
  cap *would* bind — out of scope for the current defaults, but a latent
  divergence if $K$ is ever raised past 256 without raising
  $\mathtt{kMaxKprimeCand}$.
* **M-164 log-cutoff $\mathtt{kKprimeLogCutoff}=-25$** (`mcmc_state.h:398`;
  applied `mcmc.cpp:4061, 4182`). This drops candidates whose optimistic weight
  is $>25$ nats below the per-character max. The omitted absolute likelihood
  mass is bounded by $(\text{\#dropped})\cdot e^{-25}<K\,e^{-25}\approx
  4\times10^{-10}$, so the marginal `logSumExp` is low by at most $\sim10^{-10}$
  nats. The `sampled_k` explicit sum (e.g. the test's `vapply` over `ks`,
  `:294`) has **no** cutoff and sums every $k'$.

Hence the two modes agree **exactly** in the no-pruning idealisation (A7) and
**up to $<\sim10^{-10}$ nats** in the implemented evaluator — far below the
$10^{-7}$ tolerance of the Stage-2 regression test (`:307`), and orders of
magnitude below any prior-bug signature ($\ge\sim0.01$ nats; a missing $-\log Z$
is $\sim0.26$ nats at $p=0.05,K=30$). This is the only source of
non-bit-identity between the modes; it is a benign, bounded numerical
approximation, not a modelling inconsistency.

> **Caveat (1).** The $n=1$ empirical anchor. The Stage-2 regression test
> (`test-marginal-k-truncation.R:258–333`) and the C-i guard (`:107–152`) both
> use single-character data (`mkChar`, `ncol=1`) with all `kPrime` pinned equal,
> so they verify the **per-character** identity (Part 2) and the pruning bound,
> but **not** the multi-character lift (Part 4 step $(\dagger)$). The lift is
> proved here analytically from conditional independence (A5); it relies on no
> code path beyond the per-character one already tested. A direct empirical
> check of $(\dagger)$ would require a multi-character `sampled_k` chain summed
> over the full $k'_{1:n}$ grid (combinatorially infeasible to enumerate for
> $n>$ a few), so the analytic argument is the appropriate closure. Should a
> belt-and-braces check be wanted, an mcmc-diagnostician overlap test
> (`T-OVL-sampled-vs-marginal`: do `sampled_k` and `marginal_k` chains produce
> matching $\vartheta$-posteriors on a small multi-character dataset?) is the
> right instrument — picked up by role **mcmc-diagnostician**.

> **Caveat (2).** Assumption A6 ($L_i=\mathrm{rawLL}_i$) is **empirical**, not
> analytic — it is the content of the C-i guard. Everything downstream (Parts
> 2–4) is exact algebra *given* A6. If a future refactor diverges the
> `sampled_k` likelihood path (`cpp_log_likelihood`) from the marginal
> per-$k'$ path (`compute_per_kprime_log_lik`), A6 — and hence this whole
> proof — would silently break; the C-i guard is the tripwire. No further
> work needed now; flagged so the dependency is explicit.

---

## Verdict

**Watertight** (under Assumptions 1–8; the load-bearing non-analytic ones are
A6 likelihood-equivalence, established by the C-i regression test, and A7
no-pruning, whose violation is bounded at $<\sim10^{-10}$ nats).

* **Prior normalisation (point 1):** proved; $Z_A=1-q^{K-1}$,
  $Z_{B,i}=1-q^{K-k_{\mathrm{obs},i}+1}$ are genuine truncated-geometric
  normalisers, matching `mcmc.cpp:431, 443` exactly. Limits ($K\to\infty$,
  $p\to0/1$, $K=2$, empty support) all behave.
* **Per-character RB identity (point 2):** proved for **both** models via
  `logSumExp` shift-invariance; the per-character constants
  $\mathrm{corr}_i^{A/B}$ factor out of the sum and are **bit-identical** to the
  evaluator's two branch additions (`mcmc.cpp:4369/4404` and `:4371/4406`). The
  Model A sub-$k_{\mathrm{obs}}$ prior mass is correctly accounted (in $Z_A$,
  zero likelihood) — the intended generative/unconditional semantics,
  forward-inference consistent.
* **Model A vs B (point 3):** the weight ratio $q^{k_{\mathrm{obs},i}-2}$ is
  constant **in the sum index** (factors out) but a **function of $p$** — so the
  RB identity is strictly *within-model* (sampled-A↔marginal-A,
  sampled-B↔marginal-B), and A and B are *different* posteriors on $p$ (matching
  the SBC asymmetry on record). The brief's point-3 claim is correct as a
  statement about the weights; I have flagged the over-reading that would turn
  it into "A and B target the same posterior" (they do not).
* **Total posterior (point 4):** lifted via the product-of-sums factorisation
  $(\dagger)$ (conditional independence, A5) and the mode-identical
  $k'$-independent prior block — including, critically, $\mathrm{dbeta}(p)$ at
  `mcmc.cpp:447` being **outside** the `if(!marginalK)` guard. Conclusion:
  $\log\pi_{\mathrm{marginal}}(\vartheta\mid\mathbf y)=\log\sum_{k'_{1:n}}
  \exp(\log\pi_{\mathrm{sampled}}(\vartheta,k'_{1:n}\mid\mathbf y))$ — the two
  modes target the same posterior over $\vartheta$.

**No discrepancy found between the brief's stated formulas and the source**, and
**no patch attached** (nothing to fix). The implementation at `cpp_log_prior`
(`mcmc.cpp:424–447`) and `cpp_log_likelihood_marginal` (`mcmc.cpp:4294–4410`)
realises the proved identity exactly.

**Caveats carried forward:** (1) the empirical anchor is $n=1$ — the
multi-character lift is analytic-only; an mcmc-diagnostician overlap test would
close it belt-and-braces. (2) A6 is empirical (C-i guard is the tripwire). (3)
For $K>256$, the absolute candidate cap $\mathtt{kMaxKprimeCand}=256$ would
bind and break exactness — latent, out of scope at current defaults
($K\le200$), but worth a setup assertion if $K$ is ever raised.

---

## References

* Felsenstein, J. (2004). *Inferring Phylogenies.* Sinauer. (Ch. 16,
  conditional independence of characters; pruning.)
* Lewis, P. O. (2001). A likelihood approach to estimating phylogeny from
  discrete morphological character data. *Syst. Biol.* 50:913–925. (Mkv
  ascertainment.)
* Talts, S. et al. (2018). Validating Bayesian inference algorithms with
  simulation-based calibration. arXiv:1804.06788. (SBC forward=inference
  requirement, used in the sibling normaliser proof.)
* Sibling proof: `dev/red-team/proofs/marginal-k-truncation-normaliser.md`
  (§2 $Z(p)$ closed form, §4.4 Model B normaliser, §5.3 cap+Z, §6.3 the RB
  caveat this proof closes).
* Regression anchors: `tests/testthat/test-marginal-k-truncation.R`
  (C-i guard `:107–152`; Stage-1b cap-coupling `:172–229`; Stage-2
  RB-consistency `:258–333`).
