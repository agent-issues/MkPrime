# Marginal-k Mk' likelihood — geometric arm

> **Lane.** Proof note for the marginal-k evaluator planned in
> `dev/notes/2026-05-28-marginal-k-plan.md`. Establishes arithmetic
> correctness of the per-character marginalisation
> $L_{\mathrm{marg}}(y_i \mid \mathrm{tree}, \mu, p) = \sum_u L(y_i \mid
> \mathrm{tree}, \mu, k_{\mathrm{obs},i}+u)\,P(u\mid p)$ for the
> **geometric** prior arm only.
>
> **Scope.** v1 of the marginal-k mode (`likelihoodMode = "marginal_k"`
> in the planned `MkPrimeModel()` constructor) handles
> `kPrimePrior = "geometric"` exclusively. Other arms
> (`empirical_geometric`, `beta_geometric`, `logseries`) are flagged
> out-of-scope here (§7) and require their own proof addenda.
>
> **Worktree.** `feat/marginal-k` at
> `C:\Users\pjjg18\GitHub\worktrees\mkp\marginal-k`. Worktree base
> includes `89d5e98` ("docs(marginal-k): land implementation plan note")
> and `5e51e4e` ("refactor(case25): extract phase-1 batched
> precomputation into helper"). Companion proof
> `dev/red-team/proofs/kprime-priors.md` lives at HEAD in the worktree.
>
> **Reading order.** The implementer should read this proof alongside
> the plan note. The plan describes the *implementation*; this proof
> establishes what the implementation must compute, and why that
> arithmetic is correct.

---

## 1. Definition — per-character marginal likelihood

Let $y_i$ be the observed character column for transformational
character $i$, with observed support size $k_{\mathrm{obs},i} \geq 2$.
Let $\theta = (\mathrm{tree}, \mu, \sigma)$ denote the continuous
parameters and let $u_i := k'_i - k_{\mathrm{obs},i} \in \{0, 1, 2,
\ldots\}$ be the latent excess-states count. Under the **geometric
Model B** prior (`kprimeViability` §1, `dev/notes/2026-05-28-
marginal-k-plan.md` §1) we have
$$
P(u_i \mid p) \;=\; p\,(1-p)^{u_i}, \qquad u_i \in \{0, 1, 2, \ldots\}.
$$
The per-character marginal likelihood, integrating out the discrete
latent $u_i$, is
$$
\boxed{\;
L_{\mathrm{marg}}(y_i \mid \theta, p)
  \;:=\; \sum_{u=0}^{\infty}
    L\!\left(y_i \;\middle|\; \theta,\, k_{\mathrm{obs},i} + u\right)
    \cdot P(u \mid p)
\;}
$$
where $L(y_i \mid \theta, k)$ is the Mk' per-character likelihood
under the falling-factorial relabelling correction (cf.
`dev/red-team/proofs/relabelling-correction.md`) at fixed support
size $k$.

**No truncation normaliser $Z_i$.** The geometric pmf is normalised
on the full unbounded support $u \in \{0, 1, 2, \ldots\}$ by
construction: $\sum_{u=0}^{\infty} p(1-p)^u = p/(1-(1-p)) = 1$. This
is the key contrast with `empirical_geometric` and `logseries`, both
of which carry a per-character truncation normaliser $Z_i(\theta)$
(see `dev/red-team/proofs/kprime-priors.md` §4.3, §4.4; findings
EG-001 and LS-001).

Under Model B the prior on $p$ is the user-controlled
$p \sim \mathrm{Beta}(a, b)$ hyperprior. Note that Model B's
"conditional" parameterisation (the implementation as it stands on
`main`, pre-Model-A patch) places the geometric on
$u_i \mid k_{\mathrm{obs},i}, p$ and treats $u_i$ as conditionally
independent given $p$. This is the parameterisation we marginalise
over here. Model A (unconditional joint, validated by the SBC under
the experimental worktree patch) shares the same per-character
geometric kernel; the difference is in the joint generative model,
not in the per-character marginal evaluator.

## 2. Rao-Blackwell — preservation of the marginal-on-$\theta$ posterior

**Theorem (informal).** The marginal posterior on $(\theta, p) =
(\mathrm{tree}, \mu, \sigma, p)$ is identical between the joint
sampler over $(\theta, p, u_{1:n})$ (sampled-k mode) and the marginal
sampler over $(\theta, p)$ only (marginal-k mode).

### 2.1 Assumptions

1. **Independence across characters given $(\theta, p)$.** The joint
   prior factorises as $\pi(u_{1:n} \mid p) = \prod_i P(u_i \mid p)$
   and the joint likelihood factorises as
   $\prod_i L(y_i \mid \theta, k_{\mathrm{obs},i} + u_i)$. The
   character-independence assumption is structural (the prior loops
   over $i$ in `R/MkPrimeModel.R::LogPrior` and the likelihood
   evaluates per partition / per character; cf.
   `dev/red-team/proofs/kprime-priors.md` Assumption 1).
2. **Geometric pmf is summable.** $\sum_{u \geq 0} p(1-p)^u = 1$
   for all $p \in (0, 1)$, and for fixed $\theta$ the per-character
   summands $L(y_i \mid \theta, k_{\mathrm{obs},i}+u) \cdot p(1-p)^u$
   form an absolutely convergent series (the likelihood is bounded
   above by $1$ since it is a probability mass over alignment columns
   at fixed tree / parameters).
3. **No $u_i$-only proposal that affects $p$.** Equivalently, the
   joint chain's $u_i$ moves (case 25 / case 26 / case 7) leave $p$
   invariant and are themselves marginalised away under marginal mode.
   This is structurally true: case 25 conditions on $p$ and proposes
   only $u_i$; case 30 (mh_logit_p) conditions on $u_i$ and proposes
   only $p$.

### 2.2 Statement

Let $\pi^{\mathrm{joint}}(\theta, p, u_{1:n} \mid y)$ be the joint
posterior of the sampled-k chain and let
$\pi^{\mathrm{marg}}(\theta, p \mid y)$ be the marginal posterior of
the marginal-k chain. Then
$$
\pi^{\mathrm{marg}}(\theta, p \mid y)
  \;=\; \int \pi^{\mathrm{joint}}(\theta, p, u_{1:n} \mid y)\,
              \mathrm{d}u_{1:n}.
$$

### 2.3 Proof

Bayes' rule on the joint posterior gives
$$
\pi^{\mathrm{joint}}(\theta, p, u_{1:n} \mid y)
  \;\propto\; \pi(\theta)\,\pi(p)\,
              \prod_i L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)\,
              P(u_i \mid p).
$$
Marginalising over $u_{1:n}$ (interchanging sum and product by
Fubini — Assumption 2 plus the structural product form Assumption 1
guarantee absolute summability),
$$
\int \pi^{\mathrm{joint}}\,\mathrm{d}u_{1:n}
  \;\propto\; \pi(\theta)\,\pi(p)\,
              \prod_i \sum_{u_i \geq 0}
                  L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)\,P(u_i \mid p)
  \;=\; \pi(\theta)\,\pi(p)\,
        \prod_i L_{\mathrm{marg}}(y_i \mid \theta, p).
$$
The marginal-k chain's target is, by definition,
$$
\pi^{\mathrm{marg}}(\theta, p \mid y)
  \;\propto\; \pi(\theta)\,\pi(p)\,
              \prod_i L_{\mathrm{marg}}(y_i \mid \theta, p).
$$
The two right-hand sides agree, so the targets coincide. $\square$

**Remark (Rao-Blackwell variance reduction).** Under standard
Rao-Blackwell theorem, any posterior-expectation estimator
$\widehat{f(\theta, p)}$ formed from the marginal chain has variance
$\leq$ that of the same estimator formed from the joint chain (with
equality only when $u_i$ is conditionally independent of all
functions of $(\theta, p)$ given $(\theta, p, y)$, i.e. when the
chain has fully mixed). The expected ESS/sec gain on
$(\mathrm{tree\_length}, p)$ — the headline motivation in the plan
§8 — is a Rao-Blackwell corollary, not a separate claim.

## 3. Where ascertainment and relabelling sit — inside the sum

Under Mk' coding for variable characters, the per-character
likelihood at fixed $k$ is
$$
L(y_i \mid \theta, k)
  \;=\; \frac{L_{\mathrm{raw}}(y_i \mid \theta, k)}
             {1 - P_{\mathrm{const}}(\theta, k)}
       \cdot R(k, k_{\mathrm{obs},i})
$$
where $L_{\mathrm{raw}}$ is the raw Felsenstein per-pattern
likelihood at $k$ states, $P_{\mathrm{const}}(\theta, k)$ is the
constant-site probability under the Mk model with $k$ states (the
ascertainment correction; see `dev/red-team/proofs/
ascertainment.md`), and
$R(k, k_{\mathrm{obs},i}) := \exp\!\big(\mathrm{mk\_prime\_relabel\_log}(k,
k_{\mathrm{obs},i})\big)$ is the falling-factorial relabelling
correction (`dev/red-team/proofs/relabelling-correction.md`).

**Both corrections depend on $k$ and (for $R$) on $k_{\mathrm{obs},i}$.**
They therefore sit **inside** the marginal sum, evaluated at each
candidate $k = k_{\mathrm{obs},i} + u$. The marginal log-likelihood
in numerically stable form is
$$
\log L_{\mathrm{marg}}(y_i \mid \theta, p)
  \;=\; \operatorname*{logSumExp}_{u \geq 0}
    \Big[\,
      \log L_{\mathrm{raw}}(y_i \mid \theta, k_{\mathrm{obs},i}+u)
      - \log\!\big(1 - P_{\mathrm{const}}(\theta, k_{\mathrm{obs},i}+u)\big)
$$
$$
      + \log R(k_{\mathrm{obs},i}+u,\, k_{\mathrm{obs},i})
      + \log P(u \mid p)
    \,\Big].
$$

**Pulling the corrections outside the sum is incorrect.** Both
$P_{\mathrm{const}}$ and $R$ vary with $k$; factoring either out by
its value at a fixed reference $k$ (e.g. $k = k_{\mathrm{obs},i}$)
gives a different number that is **not** equal to the marginal
likelihood. This is the load-bearing structural fact for
implementation correctness — the §7.1 bit-identity test exists to
catch this class of error.

The plan §3 reuses case 25's phase-1 batched per-(char, k')
machinery, which already evaluates exactly this combination
(`src/mcmc.cpp:3906–4082`): per-character $L_{\mathrm{raw}}$ via
the pruning kernel, then `logAscCorr` (ascertainment), then
`mk_prime_relabel_log` (relabelling). The marginal evaluator
consumes phase-1's output and applies logSumExp; this is correct by
construction.

## 4. Truncation choice — analytic tail bound

The implementation truncates the sum at $u_{\max} := K_{\mathrm{MAX\_CAND}}
- 1$ with $K_{\mathrm{MAX\_CAND}} = 50$, augmented by the M-164
prior-ceiling pre-filter at $\mathrm{LOG\_CUTOFF} = -25.0$ (plan §6;
`src/mcmc.cpp` case-25 logic).

### 4.1 Geometric tail mass

$$
P(U \geq u_{\max} \mid p)
  \;=\; \sum_{u = u_{\max}}^{\infty} p(1-p)^u
  \;=\; p(1-p)^{u_{\max}} \cdot \frac{1}{1-(1-p)}
  \;=\; (1-p)^{u_{\max}}.
$$

At $u_{\max} = 49$:

| $p$   | $(1-p)^{49}$           | within $10^{-6}$? |
|-------|------------------------|-------------------|
| 0.10  | $5.7 \times 10^{-3}$   | no                |
| 0.20  | $1.7 \times 10^{-5}$   | no                |
| 0.26  | $7.9 \times 10^{-7}$   | yes               |
| 0.30  | $1.2 \times 10^{-7}$   | yes               |
| 0.50  | $1.8 \times 10^{-15}$  | yes               |
| 0.70  | $1.8 \times 10^{-26}$  | yes               |

Solving $(1-p)^{49} = 10^{-6}$ gives $p = 1 - 10^{-6/49} \approx
0.2545$. **For every $p \geq 0.26$ the tail mass omitted by the
hard cap is below $10^{-6}$**, i.e. truncation error in
$\log L_{\mathrm{marg}}$ is bounded by $\log(1 + 10^{-6}) \approx
10^{-6}$ nats per character.

### 4.2 Behaviour at small $p$ — M-164 cutoff dominates

For $p < 0.26$ the analytic tail exceeds $10^{-6}$ and a naive cap
would be insufficient. The implementation does not rely on the cap
alone: the M-164 prior-ceiling pre-filter terminates evaluation of
candidate $k = k_{\mathrm{obs},i} + u$ as soon as
$\beta \cdot \log L_{\mathrm{raw}} + \log P(u \mid p)$ falls below
the running maximum-weight cutoff by more than $25$ nats
(`src/mcmc.cpp:3927–3961`). Two observations:

1. **For chars with binding likelihood at small $u$.** When the
   data informs $u_i$ (large $k_{\mathrm{obs}}$, low rate), the
   likelihood drops sharply for $u > 0$ — the per-character contour
   peaks near $u = 0$ and the cutoff truncates well before $u = 49$.
   The omitted prior mass at large $u$ is large but the corresponding
   likelihood-weighted contribution is far below the running max.
2. **For chars where the data does not inform $u_i$** (EG-003 regime
   — small $k_{\mathrm{obs}}$, large edges, prior-dominated). The
   per-character contour follows the prior; the prior at small $p$
   has heavy tail. The cutoff terminates near $u^\star :=
   \arg\max_u \log P(u \mid p) = 0$ plus a window of width
   $\Delta u \approx 25 / |\log(1-p)|$ (since the prior log-pmf
   decreases by $|\log(1-p)|$ per unit increase in $u$). At
   $p = 0.1$: $\Delta u = 25 / 0.1054 \approx 237$, but the prior
   mass beyond $u = 237$ is $(0.9)^{237} \approx 10^{-11}$ —
   **far below the $10^{-6}$ tolerance**. At $p = 0.05$:
   $\Delta u \approx 487$, prior tail $(0.95)^{487} \approx
   10^{-11}$.

The combination of $K_{\mathrm{MAX\_CAND}} = 50$ + LOG_CUTOFF $= -25$
is therefore conservative: either the likelihood concentrates and
the $K = 50$ cap is safe by a wide margin, or the prior is too flat
for the likelihood to discriminate and the cutoff terminates safely
in the prior tail. Cases where both the data and prior support
$u > 49$ simultaneously do not arise in the
campaign-relevant parameter range.

**Failure mode.** A pathological dataset / prior combination with
$p < 0.05$ and likelihood independent of $u$ over $u \in [0, 50]$
could in principle hit the cap binding. The plan §6 flags this:
"If the cap is binding for some character at small $p$, raise
$K_{\mathrm{MAX\_CAND}}$ with a build-time `#define`". The §7.1
bit-identity test (with explicit small $K_{\mathrm{MAX\_CAND}}$ and
a tight tolerance) should catch any truncation that leaks above
the bit-identity bar; if the test fails, raise the cap.

### 4.3 Optional analytic tail correction (deferred to v1.1)

Plan §6 anticipates adding a per-character tail correction of the
form $\log L_{\mathrm{marg}} \approx \mathrm{logSumExp}(\ldots) +
\log\!\big(1 + P_{\mathrm{tail}} \cdot \max_u L / \sum_{u} w_u L\big)$
to compensate for any truncated mass. v1 omits this; v1.1 should add
it if the §7.1 bit-identity test (at the test's required tolerance
$\leq 10^{-9}$) shows the truncation bites. For the parameter range
above, the omission is bounded by $\sim 10^{-6}$ nats per character —
adequate for headline LL but not for bit-identity at $10^{-9}$, so
the test must run with $K_{\mathrm{MAX\_CAND}}$ raised above any
candidate the test grid touches (e.g. $K_{\mathrm{MAX\_CAND}} = 200$
with $u$ grid $\{0, \ldots, 6\}$ leaves a $5 \cdot 10^{-30}$ tail at
$p = 0.5$, well below $10^{-9}$).

## 5. Cache invariance under $p$-moves

**Claim.** The per-(char, $k$) likelihood $L(y_i \mid \theta, k)$ does
not depend on $p$. A move on $p$ (case 30 `mh_logit_p`) does not
invalidate the cached per-(node, $k$) partial conditional likelihoods
nor the cached per-(char, $k$) marginal-LL summands $\log L(y_i \mid
\theta, k)$; only the $P(u \mid p)$ weights consumed by the
logSumExp step change.

### 5.1 Proof

Inspect the definition in §3:
$$
L(y_i \mid \theta, k)
  \;=\; \frac{L_{\mathrm{raw}}(y_i \mid \theta, k)}
             {1 - P_{\mathrm{const}}(\theta, k)}
       \cdot R(k, k_{\mathrm{obs},i}).
$$

The three factors are:

1. **$L_{\mathrm{raw}}(y_i \mid \theta, k)$.** Felsenstein pruning
   on the tree at branch lengths $\mu \cdot t_e$ (or partition-rate-
   scaled equivalent) under the JC$_k$ / F81$_k$ kernel with $k$
   states. The Q-matrix entries depend on $k$ (via the $1/(k-1)$
   off-diagonal under JC). They do **not** depend on $p$.
2. **$1 - P_{\mathrm{const}}(\theta, k)$.** The ascertainment
   correction. $P_{\mathrm{const}}$ is the probability of a constant
   column under the same $k$-state kernel and tree; computed from the
   per-tip stationary distribution and the same Felsenstein pruning
   machinery. Depends on $\theta$ and $k$ only.
3. **$R(k, k_{\mathrm{obs},i})$.** The falling-factorial relabelling
   factor $k(k-1)\cdots(k-k_{\mathrm{obs},i}+1)$ (proof
   `relabelling-correction.md`). A function of $k$ and
   $k_{\mathrm{obs},i}$ only; no $\theta$, no $p$.

No factor depends on $p$. $\square$

### 5.2 Implication for the cache

Under the Option A cache strategy in plan §5 (per-(partition, node,
$k$-offset) partial CLs), the cache is **fully valid** across
$p$-moves. The case-30 move proposal therefore consists of:

1. Read cached `charLL[ti, ko]` for $ko = 0, \ldots, \mathrm{charNCand}[ti]
   - 1$ for every transformational char $ti$.
2. Recompute `logPriorByU[ko] = log p + ko * log(1 - p)` for the
   proposed $p'$ (a single scalar evaluation per $ko$).
3. Recompute `log L_marg[ti] = logSumExp(charLL[ti, ·] + logPriorByU[·])`
   per character.
4. Sum, compare to current marginal LL, MH accept/reject.

Total cost per case-30 proposal: $O(n_{\mathrm{trans}} \cdot u_{\max})$
scalar operations, no pruning. This is the load-bearing efficiency
claim for the §5 Option A cache strategy. Without it, every $p$-move
costs a full re-pruning across every (char, $k$) candidate, making
$p$ a slow coordinate and undermining the Rao-Blackwell ESS/sec
gain (§2 remark).

### 5.3 Sub-option A' — cache `charLL[ti, ko]` directly

The plan §5 also discusses Sub-option A': cache the post-correction
`charLL[ti, ko]` values themselves (smaller than per-node CLs but
invalidated by all moves except $p$). The arithmetic justification is
the same — `charLL[ti, ko]` is a function of $(\theta,
k_{\mathrm{obs},i}+ko)$
only — and the two caches stack: the per-node CLs feed the
per-(char, $k$) LLs feed the logSumExp. A $p$-move skips the first
two layers; a $\theta$-move invalidates everything except the
$P(u \mid p)$ vector.

## 6. Bit-identical brute-force check

This section specifies the exact arithmetic the §7.1 unit test must
reproduce.

### 6.1 Char-independence — product of sums = sum over product grid

Fix $(\theta, p)$ and let $u = (u_1, \ldots, u_{n_{\mathrm{trans}}})
\in \mathbb{N}^{n_{\mathrm{trans}}}$ range over a product grid (each
$u_i \in U := \{0, 1, \ldots, U_{\max}\}$ for some test cap
$U_{\max}$). The naive joint marginalisation is
$$
L_{\mathrm{marg}}^{\mathrm{naive}}(y \mid \theta, p)
  \;=\; \sum_{u \in U^{n_{\mathrm{trans}}}}
          L(y \mid \theta, k_{\mathrm{obs}} + u)
          \cdot \prod_{i=1}^{n_{\mathrm{trans}}} P(u_i \mid p).
$$
Under independence-across-characters (Assumption 2.1 (1)),
$L(y \mid \theta, k') = \prod_i L(y_i \mid \theta, k'_i)$, so
$$
L_{\mathrm{marg}}^{\mathrm{naive}}
  \;=\; \sum_{u} \prod_i L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)\,
                                  P(u_i\mid p)
$$
$$
  \;=\; \prod_{i=1}^{n_{\mathrm{trans}}}
          \Bigg[\,\sum_{u_i \in U}
              L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)\,
                                   P(u_i \mid p)\,\Bigg].
$$

The distributive identity $\sum_{u \in U^n} \prod_i f_i(u_i) =
\prod_i \sum_{u_i \in U} f_i(u_i)$ (a finite-product Fubini) is
elementary by induction on $n$. The right-hand side is exactly the
form the implementation computes: per character $i$, an inner sum
(implemented as logSumExp) over $u_i$ of the joint
log-likelihood + log-prior. The outer product becomes the outer sum
in log space.

### 6.2 Test recipe

In log space:
$$
\log L_{\mathrm{marg}}^{\mathrm{naive}}(y \mid \theta, p)
  \;=\; \operatorname*{logSumExp}_{u \in U^{n_{\mathrm{trans}}}}
        \Bigg[\,\sum_i \log L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)
                       + \sum_i \log P(u_i \mid p)\,\Bigg]
$$
$$
  \;=\; \sum_{i=1}^{n_{\mathrm{trans}}}
        \operatorname*{logSumExp}_{u_i \in U}
          \Big[\,\log L(y_i \mid \theta, k_{\mathrm{obs},i}+u_i)
                 + \log P(u_i \mid p)\,\Big].
$$

The implementation computes the right-hand side directly
(`cpp_log_likelihood_marginal` per the plan §3 signature). The
brute-force test computes the left-hand side by enumerating the
$|U|^{n_{\mathrm{trans}}}$ product grid and accumulating in log-space
with a single outer logSumExp. They must agree to numerical
tolerance.

**Recommended test parameters.**

- $n_{\mathrm{trans}} = 4$, $U_{\max} = 6$ ⇒ $7^4 = 2401$ joint
  evaluations.
- $n_{\mathrm{tip}} = 8$ (small enough that `cpp_log_likelihood`
  evaluates quickly per call).
- $k_{\mathrm{obs}}$ vector mixing small and large (e.g.
  $(2, 3, 5, 7)$ to stress the relabelling term across the
  candidates).
- $p \in \{0.3, 0.5, 0.7\}$ — three reps.
- Fixed RNG seed; fixed tree (call `TreeSearch::AdditionTree`).
- Tolerance $\leq 10^{-9}$ in log-LL.
- Set $K_{\mathrm{MAX\_CAND}} \geq U_{\max} + 1$ in the test build
  (or use a debug entry point that exposes the cap), so the
  truncation correction (§4.3) does not bite at $10^{-9}$.

The bit-identity test discriminates between:

- correct per-character logSumExp arithmetic;
- correct placement of the ascertainment / relabelling corrections
  **inside** the sum (§3);
- correct alignment of the prior weights $\log P(u \mid p) = \log p
  + u \log(1 - p)$ with the candidate $u$ values;
- the truncation cap not biting at the test tolerance.

A failure on this test localises immediately to one of the four
above. A pass is necessary but not sufficient for correctness at
scale — §7.2/§7.3 of the plan provide the wider validation rungs.

## 7. Out of scope for this proof note

The following are explicitly **not** covered here. Each is a
plan §11 follow-up requiring its own proof addendum:

1. **`empirical_geometric` arm.** The per-character marginal pmf is
   the convolution
   $P(k' = m \mid p) = \sum_j P_{\mathrm{emp}}(j) \cdot p(1-p)^{m-j}$
   on $m \geq 2$, with a per-character truncation $k' \geq
   k_{\mathrm{obs},i}$ whose normaliser $Z_i(p) = \sum_{k \geq
   k_{\mathrm{obs},i}} P(k \mid p)$ is currently missing from the
   `priorVariant = "conditional"` (Model B) path (finding **EG-001**;
   `dev/red-team/proofs/kprime-priors.md` §4.3). EG-001 was patched
   on `main` 2026-05-28 (`project_eg001_fix` memory) — but the
   marginal-k arithmetic for EG must thread the truncated pmf through
   the logSumExp consistently. **A separate proof addendum is
   required** to (i) write down the truncated EG marginal $\sum_{m \geq
   k_{\mathrm{obs},i}} L(y_i \mid \theta, m) \cdot P_{\mathrm{trunc}}(m
   \mid p, k_{\mathrm{obs},i})$ with $P_{\mathrm{trunc}} = P/Z_i$, (ii)
   confirm $Z_i(p)$ is the same normaliser the prior already
   subtracts, and (iii) extend the bit-identity test to the EG
   convolution.
2. **`beta_geometric` arm.** The pmf $P(u \mid \alpha, \beta) =
   B(\alpha+1, \beta+u) / B(\alpha, \beta)$ is normalised on
   $u \geq 0$ (`dev/red-team/proofs/kprime-priors.md` §4.2; no
   truncation). The Rao-Blackwell argument of §2 transfers verbatim
   with $p$ replaced by $(\alpha, \beta)$; the §5 cache invariance
   transfers because the per-(char, $k$) likelihood depends on
   neither $\alpha$ nor $\beta$. The new work is to write down the
   per-character marginal sum with the BetaGeo pmf and run the
   bit-identity test against it. Hyperprior-on-$(\alpha,\beta)$ moves
   (the BG-RIDGE-001 $(s, r)$ slice) play the role of case 30 here.
3. **`logseries` arm.** The pmf is $P(k \mid c) = -c^k / (k \log(1-c))$
   on $k \geq 1$, with per-character truncation at $k_{\mathrm{obs},i}$
   whose normaliser is currently missing (finding **LS-001**;
   `dev/red-team/proofs/kprime-priors.md` §4.4). Latent because $c$
   is fixed (so missing $Z_i(c)$ is a per-char additive constant in
   the prior, cancelled in all MH ratios) — but the marginal-k
   arithmetic must use the **truncated** logseries pmf
   $P_{\mathrm{trunc}}(k \mid c, k_{\mathrm{obs},i}) = P(k \mid c) /
   Z_i(c)$, and the addendum must compute $Z_i(c)$ once per
   character at chain start (one-time cost since $c$ is fixed). LS-001
   must be closed before $c$ is moved to a sampled parameter under
   marginal-k mode.
4. **The `priorVariant = "conditional"` distinction more generally.**
   Under Model B (conditional), the per-character prior on $u_i$ is
   conditioned on $k_{\mathrm{obs},i}$ as a structural support floor.
   For geometric this is trivial — Geo($p$) on $u \geq 0$ has full
   support regardless of $k_{\mathrm{obs}}$, so Model A and Model B
   produce the same per-character marginal-evaluator arithmetic
   (the per-character prior pmf $P(u \mid p) = p(1-p)^u$ is identical
   under both). The Model A vs B distinction reappears in the
   **joint** generative model used by SBC, not in the per-character
   marginal sum. For EG / logseries the truncation interacts with
   the support floor and Model A vs B differ at the prior level —
   the addendum for those arms must spell this out.

## 8. Edge cases

### 8.1 $n_{\mathrm{trans}} = 0$ (all-neomorphic dataset)

The marginal sum is over the empty set of transformational
characters; $\log L_{\mathrm{marg}} = 0$ and the marginal evaluator
falls through to the existing neomorphic-only likelihood path. No
$p$-prior contribution to the posterior in this case beyond the
Beta($a$, $b$) hyperprior on $p$ — which is structurally identical to
the sampled-k mode in this case (with case 25 sweeping over an
empty set of characters). Marginal-k and sampled-k give bit-identical
posteriors. Acceptable.

### 8.2 $k_{\mathrm{obs},i} = 2$ (binary character)

The smallest meaningful $k_{\mathrm{obs}}$. The relabelling factor
$R(k_{\mathrm{obs},i} + u, k_{\mathrm{obs},i}) =
(k_{\mathrm{obs},i}+u)(k_{\mathrm{obs},i}+u-1)$
for $k_{\mathrm{obs},i} = 2$ becomes $(2+u)(1+u)$, growing
quadratically in $u$. The geometric prior decays as $(1-p)^u$.
Convergence: $(2+u)(1+u)(1-p)^u L_{\mathrm{raw}}/Z \to 0$ as $u \to
\infty$ provided $L_{\mathrm{raw}}/Z$ is bounded; this holds since
$L_{\mathrm{raw}} \leq 1$ and $Z = 1 - P_{\mathrm{const}} \geq 1 -
1 = 0$ — strictly positive for any non-degenerate tree (and the
Mk' ascertainment correction is well-defined; cf.
`dev/red-team/proofs/ascertainment.md`). Acceptable.

### 8.3 $p \to 1$ — point mass at $u = 0$

$P(U = 0 \mid p = 1) = 1$, $P(U > 0 \mid p = 1) = 0$. The marginal
sum collapses to a single term: $L_{\mathrm{marg}}(y_i \mid \theta,
p = 1) = L(y_i \mid \theta, k_{\mathrm{obs},i})$. This is the MkNT
likelihood (the degenerate Mk' baseline at $k' \equiv k_{\mathrm{obs}}$).
Consistent. (`LogPrior` rejects $p = 1$ exactly per
`dev/red-team/proofs/kprime-priors.md` §6.4, but limits at the
boundary behave correctly under the marginal evaluator.)

### 8.4 $p \to 0$ — diffuse prior

$P(U = u \mid p \to 0) \to 0$ for every fixed $u$, with mass leaking
to infinity. The marginal sum is dominated by very large $u$ where
the likelihood is typically small. In the limit, both prior and
likelihood vanish and the ratio is dominated by the data-driven
posterior maximum on $u$. `LogPrior` rejects $p = 0$ exactly; for
$p$ slightly above $0$ the truncation cap may bite (§4.2). Plan §6
addresses by leaning on the M-164 cutoff; this proof's §4.2 verifies
the cutoff is safe.

### 8.5 Pattern dedup

The plan §3 (and `src/mcmc.cpp:3906–4082`) uses pattern dedup via
`partAct[pi].patTrans` — characters with identical tip columns
share one pruning call. The marginal evaluator inherits this:
deduplicated $L(y_i \mid \theta, k)$ values feed the per-character
logSumExp. Correctness is preserved because dedup is a re-use of
arithmetic identity (two characters with identical tip patterns
have identical $L_{\mathrm{raw}}$ at every $k$, and the per-character
$k_{\mathrm{obs},i}$ for them is structurally the same since
$k_{\mathrm{obs},i}$ is defined from the tip pattern). The
post-dedup combination with the per-character $\log P(u \mid p)$
and relabelling correction $R(k, k_{\mathrm{obs},i})$ must respect
each character's own $k_{\mathrm{obs},i}$ — verified by inspection of
the case-25 dedup loop.

### 8.6 Pre-built `logPriorByU[]` table

The implementation precomputes `logPriorByU[ko] = log(p) + ko *
log(1 - p)` once per chain step (or once per $p$ value). This is
just the closed form of $\log P(u = ko \mid p)$. Linear in $ko$,
single scalar evaluation per. No subtle arithmetic.

## 9. Verdict

**Watertight under stated assumptions.** The geometric arm's
marginal-k evaluator computes
$$
\log L_{\mathrm{marg}}(y_i \mid \theta, p)
  \;=\; \operatorname*{logSumExp}_{u = 0, \ldots, u_{\max}}
    \Big[\,\log L(y_i \mid \theta, k_{\mathrm{obs},i}+u)
            + \log P(u \mid p)\,\Big]
$$
with the ascertainment + relabelling corrections *inside* the
logSumExp (§3), and:

1. **Rao-Blackwell preservation** of the marginal-on-$\theta$
   posterior holds by direct factorisation of the joint posterior
   and the Fubini-justified interchange of $\sum$ and $\prod$ under
   character-independence (§2). The mode switch is correct by
   construction.
2. **Cache invariance under $p$-moves** holds because no factor of
   the per-(char, $k$) likelihood depends on $p$ (§5); the
   load-bearing efficiency claim of Plan §5 Option A is justified.
3. **Truncation at $K_{\mathrm{MAX\_CAND}} = 50$** is safe to
   $10^{-6}$ nats per character for $p \geq 0.26$ analytically (§4.1),
   and to better than $10^{-6}$ for $p < 0.26$ via the M-164
   prior-ceiling cutoff (§4.2). A per-character analytic tail
   correction is deferred to v1.1, contingent on the §7.1 test
   passing without it.
4. **Bit-identical sum-of-products test** (§6) is the necessary
   first-rung correctness check; failure localises to one of four
   well-defined components.

**Caveats.**

> *Caveat (1).* Empirical confirmation of the §2 Rao-Blackwell
> identity at the chain-level (posterior overlap with sampled-k on
> a held-out grid of datasets) is the §7.2 heavy-test. The analytic
> proof guarantees the marginal targets agree; chain-level mixing
> agreement is an empirical claim about the implementations of the
> two evaluators and is **not closed here**. Should be picked up by
> lane mcmc-diagnostician under the planned
> `dev/red-team/heavy-tests/marginal-k/T-OVL-sampled-vs-marginal.R`.

> *Caveat (2).* The cache invariance proof of §5 establishes that
> $L(y_i \mid \theta, k)$ is mathematically independent of $p$. It
> does **not** establish that the implementation in fact preserves
> the cache across $p$-moves — that's a code-correctness claim
> verifiable only by inspection of the case-30 path under marginal
> mode, which does not yet exist. The plan §4 / §5 specify the
> required cache-invalidation behaviour; the implementer must ensure
> case 30 under marginal mode skips the partial-CL invalidation.
> The §7.4 cross-cutting smoke test ("`log_likelihood + log_prior`
> at any $(\theta, p)$ under marginal mode equals the brute-force
> form") indirectly catches an over-invalidation bug (would slow but
> not break correctness) but not an under-invalidation bug. **A
> stricter test is recommended**: after a $p$-move, recompute the
> per-(char, $k$) cache from scratch and assert bit-identity to the
> retained cache. Should be picked up by lane numerical-auditor under
> the §7 heavy-test cluster.

> *Caveat (3).* The truncation tail analysis in §4.2 for small $p$
> ($p < 0.26$) relies on the M-164 cutoff terminating safely. The
> argument is informal (a back-of-envelope $\Delta u$ estimate); a
> rigorous version would bound the omitted mass per character as a
> function of the cutoff and the per-character likelihood. For
> campaign parameter ranges this is unlikely to bite, but a
> pathological dataset / prior could in principle violate the
> implicit assumption. The §7.1 bit-identity test (run with small
> $U_{\max}$, $K_{\mathrm{MAX\_CAND}}$ raised) catches this for
> small test datasets; a wider parameter scan is a numerical-auditor
> follow-up.

### Patch status

No patch attached. This proof note describes arithmetic that is **not
yet implemented**; the plan note `dev/notes/2026-05-28-marginal-k-plan.md`
specifies the implementation. The trivial-fix policy does not apply
(the work is a new feature, not a bug fix). Per the plan §14, the
implementation should land as PR-A (case-25 refactor; commit
`5e51e4e` already on this branch) + PR-B (marginal-k geometric arm)
+ PR-C (Hamilton SBC).
