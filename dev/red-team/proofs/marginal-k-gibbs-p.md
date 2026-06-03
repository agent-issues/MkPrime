# Marginal-k Mk′ — data-augmentation Gibbs-p update (Metropolis-within-Gibbs)

> **Lane.** Red-team derivation for the `marginal-k` feature on branch
> `feat/marginal-k` at `C:\Users\pjjg18\GitHub\worktrees\mkp\marginal-k`.
>
> **Problem.** Under `likelihoodMode="marginal_k"` the latent state counts
> $k'_i$ are integrated out, so the chain holds no $u_i = k'_i - k_{\mathrm{obs},i}$
> to feed the conjugate Beta full-conditional that `sampled_k` (case 9) used.
> `p` therefore moves **only** via the random-walk `mh_logit_p` (case 30), and
> mixes poorly (measured p-ESS 45–256 at 80–200k vs `sampled_k`'s ~1250). This
> note derives a **data-augmentation Metropolis-within-Gibbs** `p`-update that
> recovers near-Gibbs mixing for `marginal_k` **without changing the target**.
>
> **Deliverable.** The exact augmented target, the imputation distribution, the
> proposal, the closed-form Metropolis accept ratio (Model A and Model B), the
> reduction to the `sampled_k` conjugate Gibbs in the untruncated limit, an
> invariance / detailed-balance argument, the numerical guards, and the code
> sites. **No patch in this note** — it is the pre-implementation correctness
> artefact; implementation follows once math-prover-checked.

---

## 0. TL;DR for the implementer

A single `p`-update sweep, run **only when the per-(char,k′) cache is valid**
(`charLLCacheReady`, support `nEff_i = charLLNCand[i]`, raw LL in
`charLLCache[i,·]`):

1. **Impute** $u_i$ for each transformational character $i$ independently:
   $$u_i \sim \mathrm{Categorical}\big(\{0,\dots,n\mathrm{Eff}_i-1\}\big),\qquad
     P(u_i = u)\ \propto\ \exp\!\big(\underbrace{\mathrm{charLLCache}[i,u] + \log p + u\log(1-p)}_{=\ \mathrm{charLogW}[i,u]}\big).$$
   (The Model-A factor $(1-p)^{k_{\mathrm{obs},i}-2}$ and the normaliser
   $Z_i(p)$ are constant in $u$, so they cancel out of the categorical.)
   Accumulate $S \equiv \sum_i u_i$ and $n \equiv |\text{trans}|$.

2. **Propose** $p^\star$ from the **untruncated** conjugate Beta (independence
   proposal):
   $$p^\star \sim \mathrm{Beta}\big(a+n,\ b + S + c_A\big),\qquad
     c_A = \begin{cases}\sum_i (k_{\mathrm{obs},i}-2) & \text{Model A (unconditional)}\\[2pt] 0 & \text{Model B (conditional)}\end{cases}$$
   with $a=\texttt{kprimeHyperA}$, $b=\texttt{kprimeHyperB}$.

3. **Accept** with probability $\min(1,e^{\log\alpha})$,
   $$\boxed{\ \log\alpha\ =\ \sum_{i\in\text{trans}}\big[\log Z_i(p)-\log Z_i(p^\star)\big]\ }$$
   $$Z_i(p) = \begin{cases} 1-(1-p)^{K-1} & \text{Model A (shared} \Rightarrow \log\alpha = n[\log Z^A(p)-\log Z^A(p^\star)])\\[2pt] 1-(1-p)^{K-k_{\mathrm{obs},i}+1} & \text{Model B (per-char sum)}\end{cases}$$
   $K=\texttt{kprimeTruncK}$ (default 200). On accept set `state->p = p*`, recompute
   `logLik` (cache-valid re-marginalise, cheap) and `logPrior`; on reject leave all three.

**Why it is fast.** For the default $K=200$ and any non-tiny $p$, $(1-p)^{K-1}$
underflows to $0$, so $Z_i\approx 1$, $\log\alpha\approx 0$, and the move is a
near-pure Gibbs draw. The Metropolis correction only bites in the deep small-$p$
tail ($p\lesssim 1/K$), where the truncation actually distorts the geometric.

**Keep `mh_logit_p` in the move set.** The new move is added, not substituted;
`mh_logit_p` guarantees irreducibility across the whole $(0,1)$ range (incl. the
small-$p$ tail where this independence proposal can stall), and the Gibbs-p just
accelerates the bulk.

---

## 1. The marginal-k target as a function of `p`

Hold the tree $T$, branch lengths, and the other scalars $\theta=(\text{rateLoss},
\text{rateLogSd},\text{rateNeo},\dots)$ fixed. The shipped marginal evaluator
`cpp_log_likelihood_marginal` (`src/mcmc.cpp:4445–4476`, cache fast-path) computes,
per transformational character $i$,
$$\log m_i(p) \;=\; \operatorname*{logSumExp}_{u\in S_i}\big(\mathrm{charLLCache}[i,u] + \log p + u\log(1-p)\big)\;+\;\Delta^A_i(p)\;-\;\log Z_i(p),$$
with support $S_i=\{0,\dots,n\mathrm{Eff}_i-1\}$ (set at cache-fill, **fixed**
during a `p`-move), raw per-(char,$k'$) log-likelihood $\mathrm{charLLCache}[i,u]
=\log L_i(u)$, the Model-A shift $\Delta^A_i(p)=(k_{\mathrm{obs},i}-2)\log(1-p)$
under `uncond` (else $0$; `src/mcmc.cpp:4469–4473`), and the truncation normaliser
$Z_i(p)$ as in §0.3.

`compute_log_prior` under `marginalK` adds **only** the bare Beta prior on $p$
(`src/mcmc.cpp:447`, gated by `if (!data.marginalK)` at `:424` which drops every
$k'$-prior and $Z$ term); all $k'$ / truncation mass lives in the likelihood.
Hence the `p`-conditional target is
$$\pi(p\mid\theta,T,\text{data}) \;\propto\; \mathrm{Beta}(p;a,b)\,\prod_{i\in\text{trans}} m_i(p),\qquad
m_i(p)=\frac{g_i(p)\sum_{u\in S_i} L_i(u)\,p(1-p)^{u}}{Z_i(p)},$$
where $g_i(p)=(1-p)^{k_{\mathrm{obs},i}-2}$ (Model A) or $1$ (Model B) — i.e.
$\Delta^A_i = \log g_i$.

**Note on exactness.** $Z_i(p)$ is the analytic full-$K$-truncation normaliser
($1-(1-p)^{K-1}$ etc.), whereas the numerator sums only over the cached/
early-terminated support $S_i$, which can be a strict subset of $\{0,\dots,K-k_{\mathrm{obs},i}\}$.
This mismatch is **the evaluator's own definition of $m_i(p)$** — the derivation
below leaves *that* $m_i(p)$ invariant exactly, whatever $S_i$ and $Z_i$ are, so
the Gibbs-p targets precisely the posterior the rest of the marginal sampler
already targets (and hence, by the RB identity in
`marginal-k-sampled-rb-consistency.md`, the `sampled_k` posterior up to the same
controlled early-termination approximation that is already shared between modes).

---

## 2. Data augmentation

Introduce one latent $u_i\in S_i$ per transformational character with augmented
joint (as a function of $p$ and $\{u_i\}$, other params fixed)
$$\pi_{\mathrm{aug}}(p,\{u_i\}) \;\propto\; \mathrm{Beta}(p;a,b)\,\prod_{i\in\text{trans}}\frac{g_i(p)\,L_i(u_i)\,p(1-p)^{u_i}}{Z_i(p)}.$$

**Marginalisation identity.** $Z_i(p)$ and $g_i(p)$ do not depend on $u_i$, so
$$\sum_{\{u_i\}} \pi_{\mathrm{aug}}(p,\{u_i\}) = \mathrm{Beta}(p;a,b)\prod_i \frac{g_i(p)\sum_{u\in S_i}L_i(u)p(1-p)^u}{Z_i(p)} = \pi(p),$$
**exactly**, for any $S_i$ / $Z_i$. So a sampler that leaves $\pi_{\mathrm{aug}}$
invariant leaves the $p$-marginal $\pi(p)$ invariant.

---

## 3. The two Gibbs sub-steps

### 3.1 Impute $u_i \mid p,\theta,T,\text{data}$ (exact conditional)
From $\pi_{\mathrm{aug}}$, the $u_i$ are conditionally independent with
$$P(u_i=u\mid p,\dots) = \frac{L_i(u)\,p(1-p)^u}{\sum_{u'\in S_i}L_i(u')\,p(1-p)^{u'}},\quad u\in S_i,$$
since $g_i(p)/Z_i(p)$ cancels between numerator and the normalising sum. In log
weights this is $\mathrm{charLLCache}[i,u] + \log p + u\log(1-p) = \mathrm{charLogW}[i,u]$,
the **exact** quantity the evaluator already forms in its per-char logSumExp
(`src/mcmc.cpp:4452–4458`). A standard Gumbel-max / inverse-CDF categorical draw
from these weights is an exact conditional sample. $\Rightarrow$ leaves
$\pi_{\mathrm{aug}}$ invariant.

### 3.2 Draw $p \mid \{u_i\},\theta,T,\text{data}$ (Metropolis-within-Gibbs)
Collect the $p$-dependence of $\pi_{\mathrm{aug}}$ at fixed $\{u_i\}$:
$$\pi(p\mid\{u_i\}) \;\propto\; \underbrace{p^{a-1}(1-p)^{b-1}}_{\text{Beta prior}}\;\cdot\; \prod_i \Big[p\,(1-p)^{u_i}\,g_i(p)\Big]\;\cdot\;\prod_i Z_i(p)^{-1}.$$
Using $\prod_i p = p^{n}$, $\prod_i(1-p)^{u_i}=(1-p)^{S}$, and $\prod_i g_i(p)=
(1-p)^{c_A}$ with $c_A=\sum_i(k_{\mathrm{obs},i}-2)$ (Model A) or $0$ (Model B),
$$\pi(p\mid\{u_i\}) \;\propto\; p^{\,a+n-1}(1-p)^{\,b+S+c_A-1}\cdot\prod_i Z_i(p)^{-1}.$$

Let the **independence proposal** be the $Z\!\equiv\!1$ conjugate density
$$q(p^\star) = \mathrm{Beta}\big(p^\star;\,a+n,\ b+S+c_A\big)\ \propto\ (p^\star)^{a+n-1}(1-p^\star)^{b+S+c_A-1}.$$
Because $q$ is exactly the Beta-kernel factor of $\pi(\cdot\mid\{u_i\})$, the
Metropolis–Hastings ratio for an independence proposal,
$$\alpha = \min\!\Big(1,\ \frac{\pi(p^\star\mid u)\,q(p)}{\pi(p\mid u)\,q(p^\star)}\Big),$$
has its entire Beta-kernel cancel — up to the proposal normaliser
$B(a+n,b+S+c_A)$, which is **identical** in $q(p)$ and $q(p^\star)$ because
$S$ and $c_A$ are fixed during the $p$-step (the imputation precedes the draw),
and up to the unknown $\pi(\cdot\mid u)$ normaliser; both cancel in the product:
$$\frac{\pi(p^\star\mid u)}{q(p^\star)} \propto \prod_i Z_i(p^\star)^{-1},\qquad
\frac{q(p)}{\pi(p\mid u)} \propto \prod_i Z_i(p),$$
$$\boxed{\ \log\alpha = \sum_i\big[\log Z_i(p) - \log Z_i(p^\star)\big].\ }$$
(The proportionality constants in the two factors are reciprocals of one another
and cancel, so the boxed result is **exact**, not merely proportional.)
This is a valid MH step targeting $\pi(p\mid\{u_i\})$ $\Rightarrow$ leaves
$\pi_{\mathrm{aug}}$ invariant.

### 3.3 Closed forms for $\sum_i\log Z_i$
* **Model A (unconditional).** $Z_i(p)=1-(1-p)^{K-1}$ is shared $\Rightarrow$
  $\sum_i\log Z_i(p) = n\,\log\!\big(1-(1-p)^{K-1}\big)$, computed as
  `n * log1p(-exp((K-1)*log1p(-p)))`.
* **Model B (conditional).** $Z_i(p)=1-(1-p)^{K-k_{\mathrm{obs},i}+1}$ per char
  $\Rightarrow$ $\sum_i\log Z_i(p)=\sum_i\log1p(-\exp((K-k_{\mathrm{obs},i}+1)\log1p(-p)))$.

Both forms are **identical** to the ones the marginal evaluator subtracts
(`src/mcmc.cpp:4353` and `:4469–4473`) and the `sampled_k` prior adds
(`src/mcmc.cpp:431, 443`), so the Gibbs-p uses the same arithmetic the rest of
the code is validated against.

---

## 4. Composite move and invariance

One Gibbs-p sweep = (3.1) then (3.2). Each sub-step leaves $\pi_{\mathrm{aug}}$
invariant; their composition therefore leaves $\pi_{\mathrm{aug}}$ invariant, and
by §2 the $p$-marginal $\pi(p)$ is preserved. Detailed balance is not required of
the composite (Gibbs sweeps need only invariance of each kernel); sub-step (3.1)
is an exact conditional draw and (3.2) is a $\pi(\cdot\mid u)$-reversible MH
kernel. Irreducibility/aperiodicity of the **whole** sampler is supplied by the
retained `mh_logit_p` move and the rest of the move set; the Gibbs-p need not be
irreducible on its own.

---

## 5. Reduction to the `sampled_k` conjugate Gibbs (untruncated limit)

As $K\to\infty$ (or, numerically, whenever $(1-p)^{K-1}$ and $(1-p^\star)^{K-1}$
underflow), $Z_i(p),Z_i(p^\star)\to 1$, so $\log\alpha\to 0$ and the move accepts
the Beta draw with probability 1:
$$p \sim \mathrm{Beta}(a+n,\ b+S+c_A).$$
* **Model B** ($c_A=0$): $p\sim\mathrm{Beta}(a+n,b+S)$ with $S=\sum_i u_i$. With
  $u_i = k'_i - k_{\mathrm{obs},i}$ imputed, this is **identical** to the
  `sampled_k` case-9 draw (`src/mcmc.cpp:5036–5038`: `shape1=a+nTrans`,
  `shape2=b+sumU`) and the pure-R reference `gibbs_p` (`R/RunMkPrime.R:4029–4031`).
* **Model A** ($c_A=\sum_i(k_{\mathrm{obs},i}-2)$): the corrected conjugate
  $p\sim\mathrm{Beta}(a+n,\ b+S+c_A)$. **N.B.** the disabled case-9 used
  `shape2=b+sumU` even under Model A, i.e. it omitted $c_A$. This is a **latent,
  not live**, error: case 9 now executes `return false` for the plain geometric
  arm (`src/mcmc.cpp:5027`), so the `b+sumU` path (`:5037`) is **unreachable** in
  shipping `sampled_k`. The omission therefore guards only a *future*
  re-enablement of the `sampled_k` Gibbs (deferred, §7) and must not be copied
  into this new move. The $c_A$ term is mandatory under Model A and follows
  directly from the $(1-p)^{k_{\mathrm{obs},i}-2}$ factor the marginal evaluator
  carries. (Numerically confirmed: dropping $c_A$ from the proposal while keeping
  the boxed accept ratio biases $E[p]$ by 0.35 — `gibbs-p-identity-check.R`.)

This is the unit-test oracle: with `kprimeTruncK` set large enough that the
truncation underflows, repeated Gibbs-p draws at fixed imputed $u$ must match
`rbeta(a+n, b+S+c_A)` in distribution, and (Model B, $k_{\mathrm{obs}}\equiv2$)
must match the `sampled_k`/R-reference draw exactly.

---

## 6. Numerical guards (implementation must include)

1. **Cache precondition.** Run the move only if `data.marginalK` and
   `state.charLLCacheReady` and the cache sizes match (`charLLNCand.size()==nTrans`,
   `charLLCache.size()==nTrans*kMaxKprimeCand`). Otherwise fall back: do one
   cache-filling marginal eval first, or return `false` (let `mh_logit_p` move
   `p` this iteration). Do **not** impute from a stale cache.
2. **Empty support.** If any `nEff_i = charLLNCand[i] <= 0` (i.e. $k_{\mathrm{obs},i}>K$,
   $-\infty$ character) the target is already degenerate — return `false`.
3. **Proposal in the open interval.** Both proposal shapes exceed 1
   ($a+n>1$ for $n\ge1$ and any $a>0$; $b+S+c_A>1$), so $\mathrm{Beta}$ has no
   boundary atom and the chain never *occupies* $0$/$1$; still, `R::rbeta` can
   round to exactly $0$/$1$ for extreme shapes, so if `p* <= 0 || p* >= 1`
   return `false` (rejecting only the boundary *proposal*, never the reverse).
4. **$\log\alpha$ finiteness.** In the deep tail $Z_i\to0$, $\log Z_i\to-\infty$;
   $\log\alpha$ may be $-\infty-(-\infty)$. Guard: if either $\sum\log Z_i(p)$ or
   $\sum\log Z_i(p^\star)$ is non-finite, return `false` (reject). `mh_logit_p`
   covers that region.
5. **logSumExp imputation.** Use the max-shift trick on `charLogW[i,·]` over
   $u\in S_i$ when forming the categorical, mirroring the evaluator, to avoid
   over/underflow; sample by inverse-CDF on the exponentiated shifted weights.
6. **State updates — FORCE COLD on accept (not the warm fast-path).** On accept:
   `state->p=p*`; **`state->charLLCacheReady=false`** then
   `state->logLik = compute_full_loglik` (which now takes the COLD path: re-prune
   the raw LLs and re-run the M-164 early-termination at $p^\star$);
   `state->logPrior = compute_log_prior`. On reject: leave `p`, `logLik`,
   `logPrior`, and the cache untouched. The new move id IS added to the entry/
   reject carve-outs (`src/mcmc.cpp:4770,5815`, `moveType != 30 && != 35`) so
   that on *entry* the cache (filled by the preceding move at the current $p$) is
   available for the imputation — but on *accept* the move deliberately refills
   it cold. **Why cold and not warm** (see §6.7): the cached candidate *support*
   is $p$-dependent, so the warm fast-path silently undercounts after a large
   downward $p$-jump.

---

## 6.7. The p-dependent cached support (why accept forces a cold refill)

The persisted cache stores, per character, the **raw** per-$(char,k')$
log-likelihoods $L_i(u)$ (p-independent) over a candidate set
$S_i=\{0,\dots,n\mathrm{Eff}_i-1\}$ whose *size* $n\mathrm{Eff}_i$ is set by the
M-164 prior-ceiling early-termination (`src/mcmc.cpp:4118,4239`, cutoff
$\log = -25$) **at the $p$ in force when the cache was filled**. The geometric
weight $p(1-p)^u$ has a heavier upper tail at smaller $p$, so candidates that are
negligible (dropped) at a fill-$p$ can carry real mass at a much smaller $p$.

Consequence: the warm fast-path (`:4445–4476`) re-weights only the *stored*
candidates. After a **large downward $p$-jump** it cannot recover the now-relevant
high-$u$ tail it never stored, and **undercounts** the marginal (empirically
$\sim10^{-6}$ nat on the 8-tip fixture; unbounded in principle for an extreme
jump). `mh_logit_p` (case 30) never hits this because its random-walk steps keep
$p^\star\approx p$, so $S_i(p^\star)\approx S_i(p)$ and warm $\approx$ cold to
$10^{-11}$. This Gibbs move makes **independence** jumps and routinely lands far
from $p$, so it MUST refill cold at $p^\star$ (re-deriving $S_i(p^\star)$).

What this means for invariance:

* The **imputation** (§3.1) always reads the cache at the *current* $p$ (the
  preceding move, or this move's own prior cold refill on accept, leaves the
  cache filled at the current $p$), so $S_i$ is the early-termination-appropriate
  support for $p$ and the imputed $u_i$ capture all mass down to $e^{-25}$.
* The **accept ratio** (§3.2) is the analytic $Z$-ratio and does *not* reference
  the cached support at $p^\star$ at all — so it is exact regardless.
* The move therefore leaves invariant $\pi_{S_i(p)}$ (the target with support
  frozen at the imputation-$p$). This differs from the ideal full-support
  $\pi_{\mathrm{full}}$ only by the early-termination tail ($\le e^{-25}$ per
  dropped candidate at the appropriate $p$; the $\sim10^{-6}$ figure above is the
  *warm-cache* artefact, which the cold refill removes from the committed state).
  This residual is **the same controlled early-termination approximation the
  entire marginal_k sampler already lives with** (every move's MH ratio compares
  a committed LL against a later cold eval), and is explicitly out of scope here
  per `marginal-k-truncation-normaliser.md`. An *exact* variant would store the
  full $p$-independent candidate set (no early termination) in the marginal
  cache; deferred as an efficiency/exactness lever, not a correctness blocker.

Verified: `tests/testthat/test-marginal-k-free-topology.R` asserts
committed $==$ cold and warm $==$ cold to $10^{-8}$ after every case-35 fire
(both Model A and B), and that the case-35 chain's $E[p]$ matches the
grid-tabulated analytic $\pi(p\mid\theta,\text{tree})$ to $0.02$.

---

## 7. Scope and what is *not* claimed

* **Plain geometric arm only.** Like the marginal evaluator (v1), this is for
  the geometric prior; `logseries` / `beta_geometric` / `empirical_geometric`
  carry no scalar conjugate `p` and must not schedule the move.
* **Neomorphic characters** carry $k=k_{\mathrm{obs}}$, no `p` dependence, and do
  not enter $n$, $S$, or $c_A$ (sum is over `transIdxGlobal` only) — matching
  case 9 and the prior.
* **Not a new target.** The move changes mixing, not the stationary
  distribution. The pre-registered acceptance gates are: (i) untruncated-limit
  reduction to `rbeta(a+n,b+S+c_A)` (and to case-9 under Model B, $k_{\mathrm{obs}}=2$);
  (ii) RB/posterior consistency on `p` between `marginal_k`+Gibbs-p and
  `sampled_k` (the overlap harness, now expected to show the p-ESS jump);
  (iii) full testthat green incl. a free-topology coherence guard that the move
  leaves `logLik`/cache coherent.
* **`sampled_k` is out of scope here** but the same MwG (with $u_i$ read from
  `state->kPrime` instead of imputed) would correctly re-enable a fast `p`-move
  there too (replacing the truncation-disabled case 9). Deferred.

---

## 8. Code sites (for the eventual patch — not applied in this note)

* New move impl + dispatch case in `src/mcmc.cpp` near the case-30 `mh_logit_p`
  block (`:5142`); reuse the cache fast-path arithmetic (`:4445–4476`) for the
  per-char weights and the $Z_i$ forms (`:4353`, `:4469–4473`).
* Add the new move id to the $p$-only carve-outs: entry invalidation
  (`:4705`, `moveType != 30`) and the post-move cache note (`:4753`, `:5800`).
* Schedule the move in `R/RunMkPrime.R` `.BuildMoves` (geometric arm,
  `marginalK` only), alongside — not replacing — `mh_logit_p`.
* Oracle/unit test: `tests/testthat/` — untruncated-limit Beta match + a
  free-topology coherence assertion; extend the overlap harness verdict.
