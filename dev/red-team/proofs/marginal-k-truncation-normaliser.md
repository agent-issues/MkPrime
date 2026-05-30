# Marginal-k Mk′ — geometric prior truncation normaliser (Model A / Model B)

> **Lane.** Red-team proof for the `marginal-k` feature on branch
> `feat/marginal-k` at `C:\Users\pjjg18\GitHub\worktrees\mkp\marginal-k`.
> This is the **geometric analog of EG-001** (the `empirical_geometric`
> truncation-normaliser fix, `project_eg001_fix` memory;
> `dev/red-team/proofs/kprime-priors.md` §4.3).
>
> **Deliverable.** A derivation of the missing per-character truncation
> normaliser $-\log Z(p)$ in the marginal-k geometric likelihood, the exact
> formula to implement (inference C++ **and** the SBC forward R), the code
> sites, the Model A / Model B and `marginal_k` / `sampled_k` scope, and the
> numerical + SBC behaviour the fix predicts.
>
> **Do not** edit `R/`, `src/`, `tests/`. No patch attached (the fix is
> non-trivial: it changes the marginal evaluator's per-character arithmetic,
> requires a declared model constant $K$, and must be coordinated with the
> SBC forward and with the Rao-Blackwell equivalence to `sampled_k`).

---

## 0. TL;DR for the implementer

* **Closed-form normaliser (Model A, unconditional).** For a truncated
  geometric on $k \in \{2,\dots,K\}$ with weight $w_k(p)=p(1-p)^{k-2}$,
  $$\boxed{\,Z(p) \;=\; \sum_{k=2}^{K} p(1-p)^{k-2} \;=\; 1-(1-p)^{K-1}.\,}$$
  Verified to machine precision (§7.1).
* **Correction term.** Subtract **one** $-\log Z(p)$ **per transformational
  character** from the marginal log-likelihood (it is the same for all
  characters under Model A, so the total is $-n_{\mathrm{trans}}\log Z(p)$).
  $Z(p)$ depends only on $(p,K)$ — **separable from the ascertainment**
  $a_k$, which lives inside the $k$-sum and depends on $(k,\text{tree})$.
* **Where (inference, C++).** `cpp_log_likelihood_marginal`
  (`src/mcmc.cpp:4190`): subtract $\log Z(p)$ from `charLL` in **both** the
  non-cache branch (`src/mcmc.cpp:4316–4322`) and the cache fast-path
  (`src/mcmc.cpp:4346–4352`), right where the `uncond` Model-A shift
  $(k_{\mathrm{obs}}-2)\log(1-p)$ is already added.
* **Where (forward, R).** `dev/red-team/heavy-tests/marginal-k/`
  `T-SBC-marginal-geometric.R:172–173`: replace the `pmin` cap
  `kTrue <- pmin(2L + rgeom(N_CHAR, p), K_MAX_PRIOR)` (which piles a point
  mass at $K$) with a **rejection-redraw** of a truncated geometric on
  $[2,K]$.
* **$K$ must be a declared model constant** used **identically** by forward
  and inference, and $Z(p)$ must be the **analytic** $1-(1-p)^{K-1}$ for that
  fixed $K$ — **not** the sum of the actually-evaluated candidate weights
  (which floats per character via the M-164 cutoff). See §5.3.
* **Scope.** Bug is in `likelihoodMode = "marginal_k"`, **geometric arm**.
  `sampled_k` is internally consistent (untruncated prior + untruncated
  latent) **but** adopting the truncated model for `marginal_k` breaks the
  exact Rao-Blackwell equivalence to `sampled_k` unless `sampled_k` is also
  truncated (§6). Model B needs a **per-$k_{\mathrm{obs}}$** normaliser
  $Z_i^B(p)=1-(1-p)^{K-k_{\mathrm{obs},i}+1}$ (§4.4).

---

## 1. The exact generative model (Model A)

### 1.1 Symbols (using the code's names)

| Symbol | Meaning | Code name |
|---|---|---|
| $p$ | geometric success prob, $p\in(0,1)$, $p\sim\mathrm{Beta}(a,b)$ | `state.p`, `kprimeHyperA/B` |
| $k'_i$ | latent true state count of transformational char $i$ | `kPrime` / marginalised |
| $k_{\mathrm{obs},i}$ | number of **distinct realised** states in column $i$, $\ge 2$ | `data.kObs[...]` |
| $u_i$ | excess count $k'_i - k_{\mathrm{obs},i}\ge 0$ | offset `ko` |
| $K$ | truncation cap on $k'$ (proposed declared constant) | `K_MAX_PRIOR` (fwd), see §5.3 |
| $w_k(p)$ | Model A weight $p(1-p)^{k-2}$ | `logPriorByU` + `uncond` shift |
| $a_k$ | $P(\text{variable}\mid k,\text{tree})=1-\mathrm{csp}_k$ | `1 - getCSP(k)` |
| $L(y_i\mid k)$ | raw Felsenstein per-pattern likelihood at $k$ states | pruning kernel |
| $R(k,k_{\mathrm{obs}})$ | falling-factorial relabelling factor | `mk_prime_relabel_log` |

### 1.2 The intended model, stated precisely

**Model A (unconditional), `coding="variable"`, Lewis-Mkv.** Per transformational
character $i$, independently given $(\text{tree},\mu,\sigma,p)$:

1. **Latent state count** drawn from a geometric *located at 2*, truncated to
   the cap $K$:
   $$k'_i \;\sim\; \mathrm{TruncGeom}_{[2,K]}(p), \qquad
     P(k'_i = k) = \frac{w_k(p)}{Z(p)} = \frac{p(1-p)^{k-2}}{Z(p)},
     \quad k\in\{2,\dots,K\}. \tag{1}$$
   This is the **only** consistent finite-support definition; §2 derives
   $Z(p)$ and §2.3 contrasts it with the *censored / point-mass-at-$K$* model
   the current forward implements.

2. **Data given $k'_i$, conditioned on variability** (Lewis 2001 Mkv;
   `dev/red-team/proofs/ascertainment.md`). The character is generated under
   the symmetric $k'_i$-state JC/Mk model on the tree and **conditioned on
   being variable**:
   $$y_i \mid k'_i=k \;\sim\; P(y\mid k,\text{tree},\text{var})
       = \frac{L(y\mid k)\,R(k,k_{\mathrm{obs}})}{a_k},
       \qquad a_k = P(\text{variable}\mid k,\text{tree}). \tag{2}$$
   The relabelling factor $R$ accounts for the unobserved-state labellings
   (`dev/red-team/proofs/relabelling-correction.md`); it and $a_k$ both
   depend on $k$ and sit **inside** the marginal (cf.
   `marginal-k-geometric.md` §3).

3. **Data constraint.** $k'_i \ge k_{\mathrm{obs},i}$: one cannot realise
   more distinct states than exist, and $L(y_i\mid k)=0$ for
   $k<k_{\mathrm{obs},i}$ (the falling factorial $R(k,k_{\mathrm{obs}})$
   vanishes). So in the marginal the sum **effectively** starts at
   $k_{\mathrm{obs},i}$, while the **prior support** $[2,K]$ — and hence
   $Z(p)$ — is the same for every character.

The chain carries $(\text{tree},\mu,\sigma,p)$ only; $k'_{1:n}$ is
analytically marginalised (`marginalK = true`). The Rao-Blackwell
preservation of the $(\text{tree},\mu,\sigma,p)$-posterior is proved in
`marginal-k-geometric.md` §2 **for the untruncated model**; §6 here revisits
it under truncation.

### 1.3 The forward semantics decide which inference is correct (read this)

There are **two** corrections people conflate, attached to **two different
forward models**. The live forward (HEAD) fixes which one is needed:

* **Current forward (commit `d337750`, "Lewis-Mkv (Model I-a) — redraw data
  only"):** draw $k'_i$ **once** from the (capped) prior, then **redraw the
  data only** (k held fixed) until the column is variable
  (`T-SBC-marginal-geometric.R:172–189`). Then $y\mid k \sim L(y\mid k)/a_k$
  exactly, so the marginal is the **sum-of-ratios**
  $\frac{1}{Z}\sum_k w_k\,L(y\mid k)/a_k$ — i.e. **$a_k$ inside the sum is
  correct**, and the *only* missing piece is the prior normaliser $Z(p)$.
  **This is the defect this proof formalises.**

* **Superseded forward (commit `e321036`, "resample-until-variable"):** draw
  the **pair** $(k,y)$ and reject the *whole pair* if invariant. That model's
  correct inference is the **ratio-of-sums**
  $\big(\sum_k w_k L\big)\big/\big(\sum_k w_k a_k\big)$ with the shared
  ascertainment normaliser $D=\sum_k w_k a_k$ — the object of the **untracked**
  proof `dev/red-team/proofs/marginal-k-ascertainment.md` (written
  2026-05-29 09:23, between `e321036` 08:24 and `d337750` 09:44). That proof's
  Assumption 7 describes the pair-rejection forward and is therefore **stale
  with respect to HEAD**: under the current redraw-data-only forward, the
  $D=\sum_k w_k a_k$ ascertainment normaliser is **not** required, and
  sum-of-ratios is correct.

> **Reader caution.** If you have `marginal-k-ascertainment.md` open: it is
> correct **for the pair-rejection forward (`e321036`)** but its premise no
> longer matches the code. The SBC evidence in this proof's brief was measured
> on the current `d337750` forward, so the live defect is the **prior
> truncation normaliser $Z(p)$ derived here**, not the ascertainment $D$.
> Both proofs can coexist; they address different (forward-model) regimes.

---

## 2. Claim and derivation of $Z(p)$

### 2.1 Claim

For Model A with cap $K$, the truncated-renormalised geometric (1) has
normaliser
$$Z(p) \;=\; \sum_{k=2}^{K} p(1-p)^{k-2} \;=\; 1-(1-p)^{K-1},
  \qquad p\in(0,1),\ K\ge 2. \tag{3}$$

### 2.2 Proof

Substitute $j=k-2$, $j=0,\dots,K-2$:
$$\sum_{k=2}^{K} p(1-p)^{k-2}
  = p\sum_{j=0}^{K-2}(1-p)^{j}
  = p\cdot\frac{1-(1-p)^{K-1}}{1-(1-p)}
  = p\cdot\frac{1-(1-p)^{K-1}}{p}
  = 1-(1-p)^{K-1}. \qquad\square$$

Sanity: $K\to\infty\Rightarrow Z\to 1$ (the untruncated geometric needs no
normaliser — this is exactly `marginal-k-geometric.md` §1's "no $Z_i$" claim,
**valid in the $K=\infty$ limit**). $K=2\Rightarrow Z=p$ (the single state
count $k=2$ has all the mass, $w_2/Z = p/p = 1$). Both correct.

$Z(p)$ is **non-decreasing in $p$**, strictly increasing on $(0,\approx 0.3)$
and numerically $\equiv 1$ beyond (e.g. $Z(0.02)=0.443$, $Z(0.05)=0.774$,
$Z(0.10)=0.953$, $Z(0.30)\approx 1$ for $K=30$; §7.1). This monotonicity is
the engine of the bias direction (§3.2, §7.2).

### 2.3 Why truncate+renormalise, not censor (pmin)

The current forward does **not** implement (1). It implements the **censored**
model
$$k'_i = \min\!\big(2+\mathrm{Geom}(p),\,K\big)
  \;\Rightarrow\;
  P(k'_i=k)=\begin{cases} p(1-p)^{k-2}, & 2\le k<K,\\[2pt]
  (1-p)^{K-1}, & k=K\ (\text{tail piled at }K).\end{cases} \tag{4}$$
(`T-SBC-marginal-geometric.R:173`, `kTrue <- pmin(2L + u_true, K_MAX_PRIOR)`.)
The mass $\sum_{k\ge K}p(1-p)^{k-2}=(1-p)^{K-1}$ that the truncated model (1)
*discards* is instead **dumped onto the single point $k=K$**. Model (4) is a
legitimate distribution, but it is **a different distribution from (1)** and,
crucially, from the unrenormalised weights the inference uses. We recommend
(1) for the model definition because:

* **Cleanliness / separability.** Under (1) the normaliser is a single
  $-\log Z(p)$ shared across all characters (Model A), trivially separable
  from ascertainment and relabelling. Under (4) the point mass at $K$ makes
  $w_K = (1-p)^{K-1}$ a special case in the weight table, complicating the
  evaluator with no modelling benefit.
* **No spurious atom.** (4) places artificial probability on "exactly $K$
  states", an artefact of the computational cap, not a modelling belief. (1)
  simply asserts "at most $K$ states", which is what a finite cap *means*.
* **Negligible at moderate $p$, controlled at small $p$.** For $p\ge 0.26$ the
  discarded tail $(1-p)^{K-1}<10^{-6}$ ($K=30$;
  `marginal-k-geometric.md` §4.1), so (1), (4), and the $K=\infty$ ideal all
  coincide. They diverge **only at small $p$** — exactly where the bug bites
  (§3). Adopting (1) on **both** sides makes forward and inference identical
  for **all** $p$, closing the gap by construction.

(The advisor on this lane recommended truncate+renormalise on **both** sides
at a common $K$; this section justifies that recommendation.)

---

## 3. The correct marginal likelihood, and the bug

### 3.1 Correct per-character marginal (Model A)

Marginalising $k'_i$ out of (1)–(2) for an observed (hence variable) $y_i$:
$$\log L_i(p) \;=\;
  \log\!\Bigg(\sum_{k=\max(2,\,k_{\mathrm{obs},i})}^{K}
    \underbrace{p(1-p)^{k-2}}_{w_k(p)}\,
    \frac{L(y_i\mid k)\,R(k,k_{\mathrm{obs},i})}{a_k}\Bigg)
  \;-\; \log Z(p), \tag{5}$$
with $Z(p)=1-(1-p)^{K-1}$ shared across **all** $i$. In numerically stable
form (matching the code's logSumExp):
$$\log L_i(p) = \operatorname*{logSumExp}_{k}\!\Big[
  \log L(y_i\mid k) + \log R(k,k_{\mathrm{obs},i}) - \log a_k + \log w_k(p)
  \Big] \;-\; \log Z(p). \tag{5$'$}$$
The first three bracket terms are **inside** the sum (correct — they depend on
$k$); $-\log Z(p)$ is **outside** the sum and **outside** the per-character
$k_{\mathrm{obs}}$ floor (it is a property of the prior on $[2,K]$, the same
for every character). **Separability claim, proved:** $Z(p)=1-(1-p)^{K-1}$ is
a function of $(p,K)$ **only**; $a_k=1-\mathrm{csp}_k$ is a function of
$(k,\text{tree})$ only. They never multiply or interact — $Z(p)$ factors
cleanly out of the $k$-sum because it has no $k$ dependence. $\square$

Total transformational log-likelihood:
$$\log L(\mathbf y\mid p) = \sum_{i=1}^{n_{\mathrm{trans}}}\log L_i(p)
  = \underbrace{\sum_i \operatorname*{logSumExp}_k[\cdots]}_{\text{what the code computes}}
    \;-\; n_{\mathrm{trans}}\log Z(p). \tag{6}$$

### 3.2 What the code actually computes (the bug)

`cpp_log_likelihood_marginal` (`src/mcmc.cpp:4190`) forms, per character,
`charLL = logSumExp_ko(charLogW[ti,ko])` with
`charLogW = β·LL + logPrior_k`, where (geometric arm)
`logPrior_k = logP + ko*log1mP` $= \log\!\big(p(1-p)^{ko}\big)$
(`src/mcmc.cpp:4114–4115`, `4236–4237`), the ascertainment $-\log a_k$ is
folded into `LL` inside the sum (`src/mcmc.cpp:3997–3999, 4128–4133`), and for
`uncond` (Model A) the shift $(k_{\mathrm{obs}}-2)\log(1-p)$ is added
(`src/mcmc.cpp:4317–4320, 4347–4350`). **There is no $-\log Z(p)$ term
anywhere** in the function (confirmed over `:4229–4353`, both the non-cache and
cache branches). So the code computes
$$\widehat{\log L_i}(p) = \operatorname*{logSumExp}_k[\cdots]
  \;=\; \log L_i(p) + \log Z(p), \tag{7}$$
i.e. **each character is too large by $\log Z(p)$** ($\le 0$, so too *small in
magnitude*; the omission inflates the likelihood at small $p$ relative to the
correct value). Equivalently the code's posterior is
$$\pi_{\mathrm{code}}(p\mid\mathbf y) \;\propto\;
  \pi(p)\prod_i \widehat{L_i}(p)
  \;=\; \pi(p)\,Z(p)^{\,n_{\mathrm{trans}}}\prod_i L_i(p)
  \;=\; Z(p)^{\,n_{\mathrm{trans}}}\;\pi_{\mathrm{correct}}(p\mid\mathbf y). \tag{8}$$

**Direction (pure algebra — regime-independent).** (8) says the code's
$p$-posterior equals the correct one **multiplied by $Z(p)^{n}$**. Since
$Z(p)=1-(1-p)^{K-1}$ is non-decreasing in $p$ and $\to 0$ as $p\to 0$, the
factor $Z(p)^n$ **up-weights large $p$ and crushes small $p$**. Therefore:

* the code's $p$-posterior is **biased high**;
* at small true $p$, $Z(p_{\mathrm{true}})^n$ is tiny, so almost no posterior
  mass sits at or below $p_{\mathrm{true}}$ ⇒ the SBC rank of $p_{\mathrm{true}}$
  collapses to **0** (the true value at the far-left tail of the posterior);
* at moderate $p$, $Z(p)\approx 1$ so $Z(p)^n\approx 1$ and the bias vanishes
  ⇒ **calibrated**.

This is **exactly** the reported SBC signature: $p_{\mathrm{true}}<0.10$ →
rank mean $5.4/134\approx 0.04$, fraction rank $\le 2$ = 0.92 (broken, biased
high); $p_{\mathrm{true}}\ge 0.10$ → rank mean $66.2/134\approx 0.49$,
fraction $\le 2$ = 0.03 (calibrated). It is "latent at moderate $p$" because
$Z(p)\equiv 1$ there, and "bites at small $p$" because the heavy geometric tail
makes $Z(p)\ll 1$ and the forward's $K$-cap / inference's truncation diverge.
The omission of $-\log Z$ **and** the censored-vs-truncated forward mismatch
push the same direction; both are fixed by adopting (1) on both sides.

`rate_log_sd` is orthogonal to the $k$-marginalisation weighting (it enters
only through ACRV rates inside $L(y_i\mid k)$, not through $w_k$ or $Z$), so it
stays calibrated (AD $p=0.425$) — consistent with the report.

---

## 4. The exact formulae to implement

### 4.1 Inference (C++) — Model A

Replace, per transformational character $i$, the code's `charLL` by
$$\mathrm{charLL}_i \;\mathrel{-}=\; \log Z(p),
  \qquad \log Z(p)=\log\!\big(1-(1-p)^{K-1}\big)
        = \mathtt{log1p(-\,exp((K-1)\,log1p(-p)))}. \tag{9}$$
Numerically stable evaluation: `(K-1)*log1p(-p)` $=(K-1)\log(1-p)$, then
`std::expm1` / `log1p`:
```
double logOneMinusP = std::log1p(-p);              // log(1-p)
double logTail       = (K - 1) * logOneMinusP;      // log (1-p)^{K-1}
double logZ          = std::log1p(-std::exp(logTail)); // log(1 - (1-p)^{K-1})
```
Subtract `logZ` from each character's `charLL` (one shared scalar; compute
once per evaluation). At small $p$, $(1-p)^{K-1}$ is bounded away from 1 unless
$p$ is extremely small; guard $p\to 0$: as $p\to0$, $\log Z\to\log((K-1)p)\to
-\infty$, and `charLL -= logZ` $\to+\infty$ per character — but this is
correct (the truncated prior places vanishing mass and the renormalisation
blows up the per-state weight); the Beta hyperprior and the $\prod$ over
characters keep the posterior proper. (`LogPrior` already rejects $p=0$
exactly.)

### 4.2 Forward (R) — Model A

Replace `T-SBC-marginal-geometric.R:172–173`
```r
u_true <- stats::rgeom(N_CHAR, p_true)
kTrue  <- pmin(2L + u_true, K_MAX_PRIOR)            # WRONG: censoring / point mass at K
```
with a **truncated-geometric** draw on $\{2,\dots,K\}$. Two equivalent methods:

* **Rejection-redraw (simplest, exact):**
  ```r
  kTrue <- integer(N_CHAR)
  for (j in seq_len(N_CHAR)) {
    repeat { kk <- 2L + stats::rgeom(1L, p_true); if (kk <= K_MAX_PRIOR) break }
    kTrue[j] <- kk
  }
  ```
  Expected redraws per draw $= 1/Z(p) = 1/(1-(1-p)^{K-1})$, $\le$ a few even at
  small $p$ for $K=30$ (e.g. $p=0.02\Rightarrow 1/0.443\approx 2.3$).
* **Inverse-CDF (vectorised):** draw $V\sim\mathrm{U}(0,1)$, set
  $u=\big\lceil \log(1-V\,Z(p))/\log(1-p)\big\rceil-1$ truncated to
  $\{0,\dots,K-2\}$, $k=2+u$. Equivalent in distribution to (1).

The redraw-**data**-only Lewis-Mkv loop (`:178–189`) is **unchanged** — only
the $k$-draw changes from `pmin` to truncated.

### 4.3 Why `pmin`-capping is wrong (one line)

`pmin` realises the **censored** law (4), with $P(k=K)=(1-p)^{K-1}$, whereas
the inference weights $p(1-p)^{k-2}$ (renormalised by $Z$) realise the
**truncated** law (1), with $P(k=K)=p(1-p)^{K-2}/Z$. These differ by the
piled tail; SBC compares the forward $k$-law against the inference $k$-law, so
the mismatch directly miscalibrates $p$ (§5.1).

### 4.4 Model B (conditional) — per-character normaliser

Under **Model B** (`priorVariant = "conditional"`, the **default non-SBC
path**), the per-character weight is $p(1-p)^{k-k_{\mathrm{obs},i}}$, i.e. the
geometric is located at $k_{\mathrm{obs},i}$. Truncated at the same cap $K$,
its normaliser is
$$Z_i^B(p) = \sum_{k=k_{\mathrm{obs},i}}^{K} p(1-p)^{k-k_{\mathrm{obs},i}}
  = 1-(1-p)^{\,K-k_{\mathrm{obs},i}+1}, \tag{10}$$
which **depends on $k_{\mathrm{obs},i}$** and is therefore **not shared** across
characters (contrast Model A's single $Z(p)$). The fix for Model B subtracts a
**per-$k_{\mathrm{obs}}$** $\log Z_i^B(p)$; characters sharing $k_{\mathrm{obs}}$
share it (tabulate $\log Z^B$ by $k_{\mathrm{obs}}$ value, exactly as EG-001
tabulated `logZByKObs`; `project_eg001_fix` "Files touched"). Note the Model A
↔ Model B identity is preserved: $w_k^A = w_k^B\,(1-p)^{k_{\mathrm{obs}}-2}$
and $Z(p) = Z_i^B(p)\,(1-p)^{k_{\mathrm{obs},i}-2}/1\ldots$ — concretely,
$Z(p)/Z_i^B(p) = \big(1-(1-p)^{K-1}\big)/\big(1-(1-p)^{K-k_{\mathrm{obs},i}+1}\big)$,
and $-\log Z + (k_{\mathrm{obs},i}-2)\log(1-p)$ (Model A) versus
$-\log Z_i^B$ (Model B) leave the **same** posterior on $p$ **iff** the same
$K$ is used. (For SBC, the forward in §4.2 must match whichever variant
inference uses; the brief's SBC harness uses Model A.)

---

## 5. SBC validity: forward = inference $\Rightarrow$ calibration

### 5.1 The theorem

**Claim.** If the data-generating prior on $k'_{1:n}$ equals the inference
prior on $k'_{1:n}$ (both Model T with the **same** $K$), and the Lewis-Mkv
data step matches (redraw-data-only forward ↔ sum-of-ratios inference), then
the SBC rank of $p_{\mathrm{true}}$ is uniform.

**Proof.** SBC (Talts et al. 2018) is exact whenever the generative joint
$\pi(\vartheta)\,\pi(\text{data}\mid\vartheta)$ used to simulate equals the
joint the inference targets, for $\vartheta=(\text{tree},\mu,\sigma,p)$ with
$k'_{1:n}$ marginalised on **both** sides. With forward (1)–(2) and inference
(5), the marginal data law is identical:
$$\pi_{\mathrm{fwd}}(y_i\mid\vartheta)
  = \sum_{k=2}^K \frac{w_k(p)}{Z(p)}\,\frac{L(y_i\mid k)R}{a_k}
  = \pi_{\mathrm{inf}}(y_i\mid\vartheta),$$
because $k'_i$ is drawn from the **same** TruncGeom and the data from the
**same** Mkv conditional. Hence the data-averaged posterior equals the prior
and the rank is uniform. $\square$

The current code violates the hypothesis **twice**: (i) forward is censored
(4) not truncated (1); (ii) inference omits $-\log Z(p)$ (and uses a floating
cap, §5.3). Both are repaired by §4.1 + §4.2. The bias direction follows from
(8), **not** from any simulation — see the honesty note in §7.3.

### 5.2 Must the inference cap equal the forward cap?

**For SBC validity: yes.** The forward $k$-law (1) has support $[2,K_{\rm fwd}]$;
the inference renormaliser must use the **same** $K=K_{\rm fwd}$ or the two
$k$-laws differ at the boundary and $p$ miscalibrates (the discrepancy is
$O((1-p)^{\min(K_{\rm fwd},K_{\rm inf})-1})$ — negligible at moderate $p$,
divergent at small $p$, i.e. exactly where SBC must pass). The brief's harness
uses `K_MAX_PRIOR = 30` for the forward; the inference must adopt $K=30$ for
$Z(p)$ (and for the candidate sum — see §5.3).

**For real-data inference (no forward): $K$ is a declared model constant.** It
should be chosen so truncation is negligible over the posterior bulk of $p$:
the discarded tail is $(1-p)^{K-1}$, so for the smallest plausible $p$ in the
posterior, pick $K$ with $(1-p_{\min})^{K-1}\ll$ the LL tolerance. For
$p_{\min}=0.05$, $K=30\Rightarrow(0.95)^{29}\approx 0.23$ (NOT negligible!);
$K=200\Rightarrow(0.95)^{199}\approx 4\times10^{-5}$. **Implication:** for
real data with appreciable posterior mass at small $p$, $K=30$ is **too small**
and the truncation is a genuine (not just numerical) modelling choice — $K$
must either be large (hundreds) or the truncated model (1) must be embraced as
*the* model with $-\log Z$ applied (which makes any $K$ self-consistent). The
latter is why (1)+(9) is the robust fix: it is correct for **any** declared
$K$, whereas "make $K$ large enough to ignore $Z$" silently re-enters the bug
whenever the posterior wanders to small $p$.

**Behaviour as $p\to 0$.** $Z(p)\to 0$; the truncated prior (1) concentrates
its (renormalised) mass near $k=2$ but with heavy spread to $K$; $-\log Z\to
+\infty$ per character. The posterior remains proper because $\prod_i L_i$ and
the Beta hyperprior dominate; `LogPrior` rejects $p=0$ exactly
(`kprime-priors.md` §6.4). The fix is well-defined on $(0,1)$.

### 5.3 Critical spec subtlety: $Z(p)$ is analytic, the candidate sum floats

> **CORRECTION (2026-05-30, after an R spec-check; supersedes the "fine for the
> numerator" claim below).** The fix is **two-part**, not one. Subtracting
> $-\log Z(p)$ **alone is wrong** and over-corrects: with the numerator still
> summing the (non-negligible at small $p$) candidates $k\in(K,\,k_{\rm obs}+49]$
> while $Z$ normalises only $[2,K]$, the $p$-posterior is biased **low** at small
> $p$ (toy SBC: AD $=0.000$, low-$p$ mean rank $0.86$). The numerator **must be
> capped at $k\le K$ per character** ($c\le K-k_{\rm obs,i}$) in **both** the
> $\max$ scan and the $\sum\exp$, in **both** code branches, **and** $-\log Z$
> subtracted. Only the two together calibrate (toy SBC, CAP+Z: low-$p$ mean rank
> $0.37$, no rank-0 piling, vs current bug $0.035$/77 % piling). Eq. (5)'s upper
> limit $K$ is the authority; the "candidate cap stays as-is" sentence below is
> **wrong** for the $(K,\,k_{\rm obs}+49]$ terms (it is right only for the truly
> negligible tail beyond $k_{\rm obs}+49$). The beyond-$K$ mass is already in the
> current numerator (`kMaxKprimeCand=50 > K=30`, $-25$ cutoff doesn't fire at
> small $p$), so capping is a real change, not a no-op.

The inference sum in (5) is evaluated over a **floating** per-character range:
$k=k_{\mathrm{obs},i}+ko$ for $ko=0,\dots,\mathrm{charNCand}[i]-1$, where
`charNCand` is bounded by `kMaxKprimeCand = 50` (`src/mcmc_state.h:378`) **and**
early-terminated by the M-164 prior-ceiling cutoff `kKprimeLogCutoff = -25.0`
(`src/mcmc_state.h:379`; `src/mcmc.cpp:4007–4037, 4148–4154`). This floating
cap is **fine for the numerator $S_i$** (it drops only negligible terms,
`marginal-k-geometric.md` §4.2) — **but the renormaliser $Z(p)$ must NOT be
the sum of the actually-evaluated candidate weights.** That sum is
*likelihood-dependent* (the cutoff depends on $L$), so renormalising by it
would (a) reintroduce a data-dependent, non-separable normaliser and (b) differ
from the forward's fixed-$K$ law. **$Z(p)$ must be the analytic
$1-(1-p)^{K-1}$ for a single declared constant $K$** (e.g. $K=$ `K_MAX_PRIOR`
shared with the forward, exposed as a model constant). The candidate-evaluation
cap (`kMaxKprimeCand`, the cutoff) is an orthogonal numerical-truncation device
and stays as-is; only the **prior** normaliser is added.

> This is the single most likely implementation error: do **not** compute
> `logZ = logSumExp(logPriorByU[0..nCand-1])`. Compute the closed form (9)
> with the fixed $K$.

---

## 6. Scope — `marginal_k` vs `sampled_k`, Model A vs B, and the RB caveat

### 6.1 `marginal_k` geometric (Model A) — **buggy, fix per §4.1/§4.2**

The marginal evaluator omits $-\log Z(p)$. Confirmed at
`src/mcmc.cpp:4316–4322` (non-cache) and `4346–4352` (cache). **Both branches**
must subtract `logZ`.

### 6.2 `marginal_k` geometric (Model B) — **same defect, per-$k_{\mathrm{obs}}$
normaliser**

Model B (`unconditionalPrior = false`, default) takes the cache as Model B
weights and never adds an $-\log Z_i^B$. The fix is (10): subtract
$\log Z_i^B(p)=\log(1-(1-p)^{K-k_{\mathrm{obs},i}+1})$ per character (shared by
$k_{\mathrm{obs}}$). Same two code sites.

### 6.3 `sampled_k` (`cpp_log_prior`, NOT the marginal evaluator) — **internally
consistent, but see RB caveat**

The `sampled_k` geometric prior (`src/mcmc.cpp:399–414`) uses
`nTrans*log(p) + sumU*log1p(-p)` $=\prod_i p(1-p)^{u_i}$, the **untruncated**
geometric on $u\ge 0$, which **sums to 1 by construction** — no $Z$, no bug.
The latent $k'_i$ is a sampled state with no analytic truncation (move bounds
only), and the conjugate $p$-draw (case 9, `src/mcmc.cpp:4855–4874`) relies on
this untruncated pmf: $p\mid k' \sim \mathrm{Beta}(a+n_{\rm trans},
b+\sum_i u_i)$. So `sampled_k` is the $K=\infty$ model, self-consistent.
**No edit to `cpp_log_prior` or case 9 is needed for `sampled_k` in
isolation.**

> **Caveat (1) — Rao-Blackwell equivalence breaks under the fix.**
> `marginal-k-geometric.md` §2 proves the marginal-k and sampled-k chains
> target the **same** $(\text{tree},\mu,\sigma,p)$-posterior — **under the
> assumption that both use the untruncated $p(1-p)^u$ on $u\ge 0$**. If
> `marginal_k` adopts the truncated+renormalised Model T (this proof's fix)
> while `sampled_k` stays untruncated ($K=\infty$, Beta-conjugate via case 9),
> the **two modes now target different posteriors**, agreeing to
> $O((1-p)^{K-1})$ — i.e. negligibly at moderate $p$ but **divergent exactly
> at small $p$** where the truncation matters. Three coherent resolutions, in
> decreasing order of cleanliness:
>
> 1. **Truncate `sampled_k` too** (TruncGeom on $k'\in[2,K]$). Then the modes
>    re-agree exactly. **Cost:** case-9 Beta conjugacy dies (the truncated
>    geometric is not Beta-conjugate); $p$ must move via the MH
>    `mh_logit_p` (case 30) instead — **which already exists and is already
>    the route for `empirical_geometric`** (case 9 returns `false` for
>    `kPriorEmpiricalGeometric` at `src/mcmc.cpp:4861`). Mirror that guard for
>    the truncated geometric.
> 2. **Accept approximate-RB.** Keep `sampled_k` untruncated; document that
>    `marginal_k` (truncated) and `sampled_k` (untruncated) coincide only to
>    $O((1-p)^{K-1})$ and may disagree on $p$ at small $p$. Acceptable **iff**
>    $K$ is large enough that this regime is outside the posterior bulk for
>    the target datasets — but see §5.2 ($K=30$ is **not** large enough at
>    $p\lesssim 0.1$).
> 3. **Embrace $K=\infty$ everywhere** (no truncation in either mode; raise
>    `kMaxKprimeCand` so the marginal sum's omitted tail is below tolerance,
>    and **drop the forward `pmin` in favour of an uncapped $2+\mathrm{Geom}$**).
>    Then $Z\equiv 1$ and **no $-\log Z$ is needed** — but the uncapped forward
>    can draw $k'$ in the hundreds at small $p$, with pruning cost $O(k'^2)$
>    or worse, so this is impractical for the SBC simulator. Rejected on
>    performance grounds.
>
> **Recommendation:** resolution 1 (truncate both, MH on $p$). Closing the
> choice requires confirming the project's stance on `marginal_k`/`sampled_k`
> bit-equivalence vs. the cost of dropping case-9 conjugacy. Should be picked
> up by lane **(marginal-k RB-consistency)** / role **math-prover** for the
> equivalence re-proof, and **mcmc-diagnostician** for the post-fix
> sampled-vs-marginal overlap test
> (`T-OVL-sampled-vs-marginal.R`, `marginal-k-geometric.md` Caveat 1).

### 6.4 Other arms

`empirical_geometric` (EG-001, fixed on `main`), `logseries` (LS-001, latent),
`beta_geometric` are **out of scope** for the marginal evaluator (it `stop()`s
on them, `src/mcmc.cpp:4199–4204`). Their truncation normalisers are tracked
separately (`marginal-k-geometric.md` §7).

---

## 7. Verification

### 7.1 $Z(p)$ closed form (machine-precision check)

`Rscript` direct sum vs closed form, $K=30$ (one-off, not committed):

| $p$ | $1-(1-p)^{29}$ | $\sum_{k=2}^{30}p(1-p)^{k-2}$ | abs. diff |
|---|---|---|---|
| 0.02 | 0.44338335 | 0.44338335 | $3.9\times10^{-16}$ |
| 0.05 | 0.77406446 | 0.77406446 | $6.7\times10^{-16}$ |
| 0.10 | 0.95289871 | 0.95289871 | $3.3\times10^{-16}$ |
| 0.20 | 0.99845257 | 0.99845257 | $3.3\times10^{-16}$ |
| 0.50 | 1.00000000 | 1.00000000 | 0 |
| 0.90 | 1.00000000 | 1.00000000 | 0 |

Formula (3) confirmed. Note $Z$ rising sharply on $(0,0.3)$ then $\equiv 1$ —
the locality engine.

### 7.2 The $Z(p)^n$ bias mechanism (no phylo needed)

From (8), $\pi_{\mathrm{code}}\propto Z(p)^n\,\pi_{\mathrm{correct}}$. Taking
the cleanest exhibition — the **prior-dominated regime** (data uninformative
about $k$, the EG-003 regime of `marginal-k-geometric.md` §4.2), where
$\pi_{\mathrm{correct}}(p)\approx\mathrm{U}(0,1)$ so
$\pi_{\mathrm{code}}(p)\propto Z(p)^n$ — at $n=100$, $K=30$ (one-off Rscript):

* Correct $p$-posterior $=\mathrm{U}(0,1)$: mean 0.500, $P(p<0.10)=0.100$.
* **Code** $p$-posterior $\propto Z(p)^{100}$: **mean 0.580, mode $\approx
  0.725$, $P(p<0.10)\approx 5\times10^{-5}$.**
* SBC rank of $p_{\mathrm{true}}$ (correct vs code):
  $0.03\!\to\!(0.030,\,0.0000)$, $0.05\!\to\!(0.050,\,0.0000)$,
  $0.08\!\to\!(0.080,\,0.0000)$, $0.15\!\to\!(0.150,\,0.009)$,
  $0.30\!\to\!(0.300,\,0.164)$, $0.50\!\to\!(0.500,\,0.403)$.

The code rank collapses to 0 for small $p_{\mathrm{true}}$ and converges to the
correct (uniform) rank as $p$ grows — **posterior-$p$-high, localised to small
$p$**, matching the reported $p<0.10$ failure and $p\ge0.10$ pass to within
Monte-Carlo error. The mechanism is **regime-independent algebra** (8); the
flat-prior case is merely the cleanest illustration.

### 7.3 Honesty note on a failed toy

A from-scratch 5-tip phylo SBC (real Felsenstein pruning, Lewis-Mkv,
forward=inference Model T) **failed to calibrate even in the FIXED
configuration** (mean rank $\approx 0.04$ at $N=800$). Diagnosis: the toy
**omitted the relabelling factor $R(k,k_{\mathrm{obs}})$** from $L(y\mid k)$,
which biases the marginal over $k$ and breaks calibration **globally** —
orthogonal to the $Z(p)$ defect. It is therefore **not cited as evidence**.
The load-bearing evidence is §7.1 (exact $Z$) + §7.2 (the $Z^n$ algebra). The
**definitive** empirical confirmation is the post-fix SBC rerun on the real
implementation.

> **Caveat (2) — empirical confirmation.** The post-fix SBC rerun (forward
> §4.2 + inference §4.1, $K=30$ both sides) must show $p$ calibrated at
> $p_{\mathrm{true}}<0.10$ (AD $p>0.4$ overall). Closing this requires running
> `T-SBC-marginal-geometric.R` post-fix. Should be picked up by role
> **mcmc-diagnostician**.

---

## 8. Implementation spec (concise, for the coder)

**Code sites (all in the worktree `feat/marginal-k`):**

| What | File:line | Change |
|---|---|---|
| Declare model constant $K$ | new (e.g. `McmcData` field `kprimeTruncK`, default = forward cap) | wire from `MkPrimeModel()`; reuse `K_MAX_PRIOR` semantics |
| Model A $-\log Z(p)$, non-cache | `src/mcmc.cpp:4316–4322` | after `uncond` shift, `charLL -= logZ` with `logZ = log1p(-exp((K-1)*log1p(-p)))` |
| Model A $-\log Z(p)$, cache path | `src/mcmc.cpp:4346–4352` | identical `charLL -= logZ` |
| Model B $-\log Z_i^B(p)$ | same two sites, B branch | `charLL -= log1p(-exp((K - kObs_i + 1)*log1p(-p)))` (tabulate by $k_{\mathrm{obs}}$) |
| Forward $k$-draw | `T-SBC-marginal-geometric.R:172–173` | replace `pmin(2+rgeom,K)` with rejection-redraw truncated geometric (§4.2) |
| `sampled_k` (case 9 conj.) | `src/mcmc.cpp:399–414, 4855–4874` | **only if** adopting RB-resolution 1 (§6.3): truncate prior + guard case 9 off like `empirical_geometric` (`:4861`), route $p$ via case 30 |

**Per-character correction term (the one number):**
$$\Delta_i = -\log Z(p),\quad Z(p)=1-(1-p)^{K-1}\ \text{(Model A, shared)};
  \qquad \Delta_i = -\log\!\big(1-(1-p)^{K-k_{\mathrm{obs},i}+1}\big)\ \text{(Model B)}.$$
$Z$ is **analytic at fixed declared $K$** — never the floating candidate sum
(§5.3). Separable from ascertainment $a_k$ (proved §3.1).

**Tests to add (mirroring EG-001's regression block,
`project_eg001_fix` "Tests"):** (i) no-op at $K\to\infty$ ($Z\to1$, bit-equal
to current); (ii) hand-verified $-\log Z$ at small $K$ and several $p$;
(iii) R↔C++ agreement of the marginal log-lik to $10^{-9}$ across mixed
$k_{\mathrm{obs}}\in\{2,3\}$, three $p$, Model A and B; (iv) the §6 brute-force
bit-identity test (`marginal-k-geometric.md` §6) re-run **with $-\log Z$**.

---

## 9. Verdict

**Implementation bug (prior-normalisation), confirmed and localised.** The
`marginal_k` geometric evaluator (`cpp_log_likelihood_marginal`,
`src/mcmc.cpp:4190`) computes the marginal log-likelihood **without** the
truncation normaliser $-\log Z(p)$ that the finite candidate sum requires,
and the SBC forward (`T-SBC-marginal-geometric.R:173`) draws a **censored**
(`pmin`) rather than **truncated** geometric. Under the current
redraw-data-only Lewis-Mkv forward (commit `d337750`), the ascertainment
$a_k$-inside-the-sum is **correct** (the competing untracked proof
`marginal-k-ascertainment.md` addresses the **superseded** pair-rejection
forward `e321036` and is stale w.r.t. HEAD); the **only** live defect is the
prior normaliser.

* **Closed form:** $Z(p)=1-(1-p)^{K-1}$ (Model A, §2; verified §7.1);
  $Z_i^B(p)=1-(1-p)^{K-k_{\mathrm{obs},i}+1}$ (Model B, §4.4).
* **Fix:** subtract $\log Z(p)$ per character in **both** marginal branches
  (`src/mcmc.cpp:4316–4322` and `4346–4352`), with $Z$ **analytic at a
  declared fixed $K$** (§5.3); change the forward to a truncated-geometric
  draw (§4.2); $K$ must match between forward and inference for SBC (§5.2).
* **Direction & locality:** the code's posterior $=Z(p)^{n}\times$ correct
  (§3.2, eq. 8), pushing $p$ **high** and crushing small $p$ — SBC rank of
  $p_{\mathrm{true}}$ collapses to 0 for $p_{\mathrm{true}}<0.10$ and is
  uniform for $p_{\mathrm{true}}\ge0.10$, **exactly the reported signature**
  (§7.2). `rate_log_sd` stays calibrated (orthogonal).
* **Scope:** `marginal_k` Models A and B both need the fix; `sampled_k` is
  internally consistent (untruncated) and needs **no** edit in isolation, but
  the fix **breaks the exact Rao-Blackwell equivalence** to `sampled_k`
  unless `sampled_k` is also truncated (Caveat 1, §6.3) — recommend truncating
  both and routing $p$ via the existing `mh_logit_p` (case 30), mirroring
  `empirical_geometric`.

**Caveats:** (1) RB equivalence to `sampled_k` — §6.3, math-prover +
mcmc-diagnostician follow-up. (2) Post-fix SBC rerun is the empirical
confirmation — §7.3, mcmc-diagnostician. (3) For real-data inference with
posterior mass at small $p$, $K=30$ is **not** large enough to ignore $Z$
(§5.2); the analytic $-\log Z$ fix is robust for any $K$, "make $K$ huge" is
not.

No patch attached: the fix is non-trivial (new declared $K$ constant, two
evaluator branches, forward change, and the RB-consistency decision for
`sampled_k`), exceeding the trivial-fix threshold.

---

## SUMMARY (one paragraph a coder can act on)

The `marginal_k` geometric likelihood is missing a per-character prior
truncation normaliser. **Add $-\log Z(p)$ to each transformational
character's `charLL`** in `cpp_log_likelihood_marginal`
(`src/mcmc.cpp:4316–4322` non-cache **and** `4346–4352` cache), where for
**Model A** (`unconditionalPrior=true`) $Z(p)=1-(1-p)^{K-1}$ is shared across
characters, and for **Model B** (default) $Z_i^B(p)=1-(1-p)^{K-k_{\mathrm{obs},i}+1}$
is keyed by $k_{\mathrm{obs}}$ (tabulate like EG-001's `logZByKObs`). Compute
$Z$ as the **closed form at a single declared constant $K$**
(`log1p(-exp((K-1)*log1p(-p)))`), **not** as the sum of the floating
candidate weights. $Z(p)$ depends only on $(p,K)$ and is fully separable from
the ascertainment $a_k=1-\mathrm{csp}_k$ (which stays inside the $k$-sum). For
SBC, also fix the forward (`T-SBC-marginal-geometric.R:172–173`): replace
`kTrue <- pmin(2L + rgeom(N_CHAR, p), K_MAX_PRIOR)` (censoring — point mass at
$K$) with a rejection-redraw of a truncated geometric on $[2,K]$, using the
**same** $K$ as inference. Without the fix, the $p$-posterior is the correct
one multiplied by $Z(p)^{n_{\rm trans}}$, biased high and crushing small $p$,
which is why the SBC rank of $p_{\mathrm{true}}$ collapses to 0 only for
$p_{\mathrm{true}}<0.10$ (heavy-tail regime where $Z\ll1$) and is calibrated
above (where $Z\approx1$). Heads-up: this fix breaks the exact
Rao-Blackwell equivalence with `sampled_k` (which stays untruncated) at small
$p$ — either truncate `sampled_k` too and move $p$ via the existing
`mh_logit_p` (case 30, as already done for `empirical_geometric`), or document
the small-$p$ discrepancy. `rate_log_sd` is unaffected; `sampled_k` in
isolation has no bug.
