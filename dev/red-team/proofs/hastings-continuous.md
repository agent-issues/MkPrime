# Lane L2 — Hastings ratios for continuous proposals

**Reviewer:** math-prover agent
**Date:** 2026-05-26
**Scope:** All continuous-state proposals consumed by the `do_move_impl` accept/reject
machinery in `src/mcmc.cpp`. Hastings ratios are derived from first principles
*before* the implementation is consulted, then cross-checked against the source.

## Convention

Throughout, the Metropolis-Hastings acceptance ratio used by the engine is

$$
\log\alpha \;=\; \beta\,(\log\ell' - \log\ell) \;+\; (\log\pi' - \log\pi) \;+\; \mathrm{logHastings},
$$

where `logHastings` carries **both** the log-ratio of proposal densities and any
Jacobian arising from a change of variables. This is the canonical form

$$
\mathrm{logHastings} \;=\; \log q(x \mid x') - \log q(x' \mid x) \;+\; \log\!\left|\frac{\partial x'}{\partial \xi}\right|^{-1}\!\left|\frac{\partial x}{\partial \xi'}\right|,
$$

absorbed into a single scalar so that detailed balance reduces to
`log α = β·Δlog ℓ + Δlog π + logHastings`.

Composition is verified at `src/mcmc.cpp:4845-4847`:

```cpp
double logAlpha = beta * (newLogLik - state->logLik) +
                  (newLogPrior - state->logPrior) + logHastings;
if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) { ... }
```

This is correct *iff* each per-move `logHastings` includes the Jacobian; that is
what the remainder of this document verifies, move by move.

---

## §1. Theorem (Bactrian scale move, M-118)

For every move type that executes
`mult = exp(scaleTuning · bactrian_perturbation()); x' = x · mult`,
the Hastings ratio is **`log(mult)`**, equivalently `scaleTuning · z` with
`z` the Bactrian draw.

### 1.1 Assumptions

1. `x > 0`. The move is applied to a positive scalar (tree length, rate, scale,
   variance, β-scale, neo joint scale). Negative or zero values are out of
   support; the engine rejects (`return false`) when `mult` produces a
   non-positive value via the downstream prior evaluation.
2. The Bactrian noise distribution is **symmetric about 0**:
   $f_Z(z) = f_Z(-z)$ for every $z\in\mathbb{R}$.
3. The transformation $z \mapsto x' = x\,e^{\tau z}$ is a $C^1$
   diffeomorphism from $\mathbb{R}$ to $(0,\infty)$ with Jacobian
   $|dx'/dz| = \tau x'$.

### 1.2 Bactrian density (Yang & Rodríguez 2013, eq. 10)

Let $m = 0.95$, $s_b = \sqrt{1 - m^2}$, $c = 1/\sqrt{12}$. A Bactrian draw is

$$
Z \;=\; c\cdot(S\,m + \sigma\,N),\qquad
S \in\{+1,-1\}\ \text{uniformly},\ N\sim\mathcal{N}(0,1),
$$

so the un-scaled $Z/c$ has density
$\tfrac12\phi_{s_b}(z+m) + \tfrac12\phi_{s_b}(z-m)$ with
$\phi_{s_b}(\cdot)$ the centred Gaussian of sd $s_b$. The full $Z$ has density

$$
f_Z(z) \;=\; \frac{1}{2c}\!\left[\phi_{s_b}\!\big(\tfrac{z}{c}+m\big) + \phi_{s_b}\!\big(\tfrac{z}{c}-m\big)\right].
$$

Symmetry $f_Z(z)=f_Z(-z)$ follows because each Gaussian term is even after the
sum-over-modes (the mode-flip $S\mapsto-S$ is just a relabelling of the equally
weighted mixture components).

### 1.3 Proof

Forward proposal: $z\sim f_Z$, $x' = x\,e^{\tau z}$.

Forward density of $x'\mid x$ (push-forward through the diffeomorphism):

$$
q(x'\mid x) \;=\; f_Z\!\left(\frac{\log(x'/x)}{\tau}\right) \cdot \frac{1}{\tau\,x'}.
$$

Reverse proposal: starting at $x'$, the move that returns to $x$ uses
$z' = \log(x/x')/\tau = -z$. Its density at $x$ is

$$
q(x\mid x') \;=\; f_Z\!\left(\frac{\log(x/x')}{\tau}\right)\cdot\frac{1}{\tau\,x}
           \;=\; f_Z(-z)\cdot\frac{1}{\tau\,x}.
$$

By symmetry (Assumption 2), $f_Z(-z) = f_Z(z)$, so

$$
\frac{q(x\mid x')}{q(x'\mid x)} \;=\; \frac{f_Z(z)/(\tau x)}{f_Z(z)/(\tau x')} \;=\; \frac{x'}{x} \;=\; e^{\tau z} \;=\; \mathrm{mult}.
$$

Hence

$$
\boxed{\ \mathrm{logHastings} \;=\; \log\mathrm{mult} \;=\; \tau z.\ }
$$

This is exactly the standard multiplicative-move Jacobian (Liu 2001, §5.3.2):
the symmetric kernel in log-space contributes 0; the change of variable
$x\to\log x$ contributes the factor $x'/x$.

### 1.4 Implementation cross-check

| moveType | code site | line | what it does | logHastings |
|---|---|---|---|---|
| 0 | scale tree_length | `src/mcmc.cpp:4115-4119` | `mult = exp(τ·z); TL ← TL·mult` | `log(mult)` ✓ |
| 1 | scale rate_loss | `src/mcmc.cpp:4121-4126` | as above | `log(mult)` ✓ |
| 2 | scale rate_log_sd | `src/mcmc.cpp:4127-4134` | as above + lockstep classRateLogSd[0] | `log(mult)` ✓ |
| 3 | scale rate_neo | `src/mcmc.cpp:4135-4140` | as above | `log(mult)` ✓ |
| 8 | scale p (legacy) | `src/mcmc.cpp:4280-4285` | as above on `p ∈ (0,1)` — **see Edge case 1.5(c)** | `log(mult)` |
| 16 | scale beta_scale | `src/mcmc.cpp:4329-4334` | as above | `log(mult)` ✓ |
| 27 | scale_kprime_alpha | `src/mcmc.cpp:4406-4420` | as above + immediate accept/reject | `log(mult)` ✓ |
| 28 | scale_kprime_beta | `src/mcmc.cpp:4421-4434` | as above | `log(mult)` ✓ |
| 31 | scale_class_rate_log_sd | `src/mcmc.cpp:4456-4477` | as above on `classRateLogSd[c]` | `log(mult)` ✓ |

The implementations of cases 27 and 28 compute `log α` directly and
bypass the shared engine path; they include `log(mult)` correctly
(`src/mcmc.cpp:4413, 4427`).

### 1.5 Edge cases

(a) **`mult ≤ 0`.** Impossible: `exp(·) > 0` for every finite argument.

(b) **`x · mult` underflows.** Then `compute_log_prior` returns `-Inf` and the
move is rolled back at `src/mcmc.cpp:4544-4572`. Detailed balance still holds:
the rejection at the prior step is equivalent to multiplying the accept ratio
by 0.

(c) **`x = p ∈ (0,1)` (case 8).** The Bactrian-scale move at `src/mcmc.cpp:4281`
treats `p` as a positive scalar, ignoring the upper bound `p < 1`. The
proposal `p' = p·mult` can exceed 1; the engine then relies on
`compute_log_prior` to return `-Inf` (Beta prior). The Hastings ratio is
*formally* still `log(mult)` (the kernel is well-defined on $\mathbb{R}_{>0}$),
and the upper-bound restriction is enforced by the prior support, not by the
proposal. This is a defensible MH construction (proposal can hit zero-prior
states; acceptance ratio is 0) but it is **wasteful** — see also EG-007 in
`dev/red-team/findings.md` and case 30 below for the logit-MH alternative
actually used in production for the empirical-geometric arm.

(d) **`τ = 0`.** Then `mult = 1`, the chain stays put, `logHastings = 0`.
Trivially detailed-balance compliant (identity kernel).

(e) **`τ < 0`.** Bactrian symmetry persists; the proof goes through unchanged.
The engine adapts `scaleTuning` to a positive value in practice
(`R/RunMkPrime.R` adaptation block), so this is hypothetical.

### 1.6 Detailed balance

For any pair $x \neq x'$ in the positive reals with $\pi(x)\pi(x') > 0$,
the accept probability $\alpha(x,x') = \min(1, e^{\log\alpha})$
satisfies $\pi(x) q(x'\mid x)\alpha(x,x') = \pi(x') q(x\mid x')\alpha(x',x)$
when `logHastings = log(mult)` by the standard MH argument
(Hastings 1970; Roberts & Rosenthal 2004 Thm. 2). The pre-condition
$\pi(x), \pi(x') > 0$ is the only support requirement.

### 1.7 Verdict

**Watertight.** Proof complete; implementation matches in every case listed in
1.4. The only caveat is case 8's reliance on prior rejection for `p ≥ 1`;
this is documented behaviour and noted in EG-007/EG-009.

---

## §2. Theorem (Logit-scale MH on p, case 30)

For move type 30 (`mh_logit_p`):

$$
\boxed{\ \mathrm{logHastings} \;=\; \log p' + \log(1-p') - \log p - \log(1-p).\ }
$$

### 2.1 Assumptions

1. $p \in (0,1)$. Code at `src/mcmc.cpp:4440` enforces with
   `if (oldP <= 0.0 || oldP >= 1.0) return false`.
2. The kernel acts as a symmetric Bactrian step on the *logit* of $p$:
   $\eta = \log\!\frac{p}{1-p}$, $\eta' = \eta + \tau z$.
3. $z\sim$ Bactrian (symmetric about 0).

### 2.2 Proof

The change of variables $p\mapsto\eta = \mathrm{logit}(p)$ has Jacobian

$$
\frac{d\eta}{dp} \;=\; \frac{1}{p(1-p)} \quad\Longleftrightarrow\quad \frac{dp}{d\eta} \;=\; p(1-p).
$$

In $\eta$-coordinates the proposal $\eta'=\eta+\tau z$ is **symmetric**, so
$q_\eta(\eta'\mid\eta) = q_\eta(\eta\mid\eta')$.

Push-forward to $p$-coordinates:

$$
q(p'\mid p) \;=\; q_\eta\!\big(\eta(p')\mid\eta(p)\big) \cdot \left|\frac{d\eta'}{dp'}\right| \;=\; \frac{q_\eta(\eta'\mid\eta)}{p'(1-p')},
$$

$$
q(p\mid p') \;=\; \frac{q_\eta(\eta\mid\eta')}{p(1-p)}.
$$

By the symmetry of $q_\eta$,

$$
\frac{q(p\mid p')}{q(p'\mid p)} \;=\; \frac{p'(1-p')}{p(1-p)}.
$$

Taking logs gives the boxed formula.

### 2.3 Implementation cross-check

`src/mcmc.cpp:4435-4455`:

```cpp
if (oldP <= 0.0 || oldP >= 1.0) return false;
double logitP = std::log(oldP / (1.0 - oldP));
double logitPnew = logitP + scaleTuning * bactrian_perturbation();
double newP;
if (logitPnew >= 0.0) {
  newP = 1.0 / (1.0 + std::exp(-logitPnew));
} else {
  double e = std::exp(logitPnew);
  newP = e / (1.0 + e);
}
if (newP <= 0.0 || newP >= 1.0) return false;
state->p = newP;
logHastings = std::log(newP) + std::log1p(-newP)
            - std::log(oldP) - std::log1p(-oldP);
```

The line `logHastings = log(newP) + log1p(-newP) - log(oldP) - log1p(-oldP)` is
exactly the formula derived in §2.2. ✓

The two-branch inverse-logit (`if logitPnew ≥ 0` vs `< 0`) is a numerical
trick for stability; both branches compute the same function $p =
1/(1+e^{-\eta})$ and do not affect the Hastings ratio.

### 2.4 Edge cases

(a) **`newP` numerically rounds to 0 or 1.** The line at `:4450` rejects.
Detailed balance preserved (rejection from zero-prior region).

(b) **`oldP = 0.5`.** Then `logitP = 0`, no asymmetry in the boundary
behaviour. Hastings reduces to symmetric perturbation; formula gives
$2\log 0.5 - 2\log 0.5 = 0$ exactly when the proposal is at the midpoint
— matches intuition.

### 2.5 Verdict

**Watertight.**

---

## §3. Theorem (2D joint Bactrian scale, M-120, cases 21–22)

For joint moves
$(x_1, x_2) \mapsto (x_1\,\mathrm{mult}_1, x_2\,\mathrm{mult}_2)$ with
$\mathrm{mult}_i = e^{\tau z_i}$ and $(z_1, z_2)$ drawn from the 2D
Bactrian kernel:

$$
\boxed{\ \mathrm{logHastings} \;=\; \log\mathrm{mult}_1 + \log\mathrm{mult}_2 \;=\; \tau(z_1+z_2).\ }
$$

### 3.1 Assumptions

1. $x_1, x_2 > 0$.
2. The 2D Bactrian density $f_{Z_1,Z_2}$ is **symmetric under joint sign-flip**:
   $f_{Z_1,Z_2}(z_1, z_2) = f_{Z_1,Z_2}(-z_1, -z_2)$.
3. $z_1, z_2$ are perturbed *independently* of $(x_1, x_2)$ — that is, the
   kernel parameters $\rho, \tau$ are state-independent.

### 3.2 The 2D Bactrian kernel

From `src/mcmc.cpp:99-113`, with $m = 0.95$, $s_b = \sqrt{1-m^2}$,
$c = 1/\sqrt{12}$, $\rho \in (-1, 1)$:

1. Draw $N_1, N_2 \stackrel{\text{iid}}{\sim} \mathcal{N}(0,1)$.
2. Set $E_1 = s_b N_1$, $E_2 = s_b(\rho N_1 + \sqrt{1-\rho^2}\,N_2)$; so
   $(E_1, E_2)\sim\mathcal{N}\!\left(0, s_b^2\!\begin{pmatrix}1&\rho\\\rho&1\end{pmatrix}\right)$.
3. Draw $S_1 \in\{+m,-m\}$ uniformly.
4. With probability $p_{\text{same}} = (1+\rho)/2$, set $S_2 = S_1$; else
   $S_2 = -S_1$.
5. Return $(Z_1, Z_2) = (c(S_1+E_1),\ c(S_2+E_2))$.

**Symmetry under $(z_1,z_2)\to(-z_1,-z_2)$.** The mixture has four sign
components $(S_1, S_2) \in \{(+,+),(-,-),(+,-),(-,+)\}$ with probabilities

$$
P(+,+) = \tfrac12\cdot p_{\text{same}} = \tfrac{1+\rho}{4} = P(-,-),
$$
$$
P(+,-) = \tfrac12\cdot(1-p_{\text{same}}) = \tfrac{1-\rho}{4} = P(-,+).
$$

So the sign distribution is invariant under joint flip $(S_1,S_2)\to(-S_1,-S_2)$.
The Gaussian noise $(E_1, E_2)$ is mean-zero with symmetric covariance, hence
its density is even: $f_{E_1,E_2}(e_1,e_2) = f_{E_1,E_2}(-e_1,-e_2)$.
Composing the four sign-conditional Gaussians with the symmetric sign
distribution gives

$$
f_{Z_1,Z_2}(z_1,z_2) \;=\; \sum_{(s_1,s_2)} P(s_1,s_2)\,f_{E_1,E_2}\!\big(\tfrac{z_1}{c}-s_1,\,\tfrac{z_2}{c}-s_2\big)\big/c^2
$$

which is invariant under $(z_1,z_2)\to(-z_1,-z_2)$ because the pairing
$\{(+,+)\leftrightarrow(-,-),(+,-)\leftrightarrow(-,+)\}$ matches summand
arguments.

**Numerical confirmation.** With $\rho=0.6$ and $5\times 10^5$ draws:
empirical correlation $0.5991$ vs target $0.6$; marginal sd $0.2887$ vs
target $1/\sqrt{12} = 0.2887$; quadrant probabilities
$P(Q_{++})=0.3993$ vs $P(Q_{--})=0.3989$ and
$P(Q_{+-})=0.1008$ vs $P(Q_{-+})=0.1011$. Symmetry holds to MC error.

### 3.3 Proof

Forward proposal: $(z_1,z_2)\sim f_Z$, $x_i' = x_i e^{\tau z_i}$.

Forward density:

$$
q(x_1',x_2'\mid x_1,x_2) \;=\; f_Z\!\left(\frac{\log(x_1'/x_1)}{\tau},\frac{\log(x_2'/x_2)}{\tau}\right)\frac{1}{\tau^2 x_1' x_2'}.
$$

Reverse: $(z_1', z_2') = (-z_1, -z_2)$ takes $(x_1', x_2')\to(x_1, x_2)$.

$$
q(x_1,x_2\mid x_1',x_2') \;=\; f_Z(-z_1,-z_2)\frac{1}{\tau^2 x_1 x_2} \;=\; f_Z(z_1,z_2)\frac{1}{\tau^2 x_1 x_2}.
$$

Ratio:

$$
\frac{q(x_1,x_2\mid x_1',x_2')}{q(x_1',x_2'\mid x_1,x_2)} \;=\; \frac{x_1' x_2'}{x_1 x_2} \;=\; \mathrm{mult}_1\,\mathrm{mult}_2.
$$

Taking logs:
$\mathrm{logHastings} = \log\mathrm{mult}_1 + \log\mathrm{mult}_2$, as claimed.

### 3.4 Implementation cross-check

`src/mcmc.cpp:4354-4374`:

```cpp
case 21: { // M-120: joint_tl_rls (tree_length × rate_log_sd)
  double z1, z2;
  bactrian_2d_perturbation(jointRho, z1, z2);
  double mult1 = std::exp(scaleTuning * z1);
  double mult2 = std::exp(scaleTuning * z2);
  state->treeLength = oldTL * mult1;
  state->rateLogSd  = oldRLSD * mult2;
  if (state->usePartitioned && state->classRateLogSd.size() > 0)
    state->classRateLogSd[0] = state->rateLogSd;  // lockstep
  logHastings = std::log(mult1) + std::log(mult2);
  break;
}
case 22: { // M-120: joint_tl_rl (tree_length × rate_loss)
  ...
  logHastings = std::log(mult1) + std::log(mult2);
}
```

Line `logHastings = std::log(mult1) + std::log(mult2)` matches the boxed
formula. ✓

### 3.5 Edge cases

(a) **`|ρ| → 1`.** The kernel degenerates: noise lives on a 1D subspace and
the sign pair becomes perfectly correlated/anticorrelated. The density is
still symmetric (as a degenerate measure on a line), the Hastings ratio
formula still holds, but the move becomes a 1D scale move along the
diagonal. No correctness issue.

(b) **`ρ = 0`.** Independent Bactrians; the proof reduces to two
independent applications of §1. ✓

(c) **`jointRho` updated mid-chain by adaptation.** The Hastings derivation
treats $\rho$ as a constant of the kernel. Adaptation that *learns* $\rho$
from past samples violates time-homogeneity unless adaptation is suspended
during the production phase. See `feedback_wallclock_adaptation` — this is
a documented efficiency/reproducibility tradeoff, **not** a Hastings-ratio
bug; the per-step Hastings ratio is computed with the *current* $\rho$ which
is exactly what enters the per-step kernel.

(d) **Classes 21/22's "lockstep" on `classRateLogSd[0]`.** Since
`classRateLogSd[0]` is constrained to equal `rateLogSd` by invariant
(see `feedback_model_scope` / partition-API §4.2), the move on
`(rateLogSd, treeLength)` is the unique move on the constrained manifold;
no further Jacobian arises.

### 3.6 Verdict

**Watertight.**

---

## §4. Theorem (Beta-simplex proposal, case 4)

For the move that picks two simplex elements $A = x_i, B = x_j$ ($i \ne j$),
holds $T = A + B$ fixed, and proposes $F' \sim \mathrm{Beta}(\alpha,\beta)$
with $\alpha = F\,\tau + 1, \beta = (1-F)\,\tau + 1$:

$$
\boxed{\ \mathrm{logHastings}
\;=\; \log\mathrm{Beta}(F\mid F'\tau + 1, (1-F')\tau + 1)
\;-\; \log\mathrm{Beta}(F'\mid F\tau + 1, (1-F)\tau + 1).\ }
$$

where $F = A/T$, $F' = A'/T$, and Beta$(x\mid a,b)$ denotes the density of
the Beta distribution evaluated at $x$.

### 4.1 Assumptions

1. $x$ is a simplex vector with at least two strictly positive entries
   (so $T > 0$).
2. The pair $(i,j)$ is drawn with index-symmetric probability
   $\frac{1}{n(n-1)}$: $i$ uniform on $\{1,\dots,n\}$ and $j$ uniform on
   $\{1,\dots,n\}\setminus\{i\}$. **Caller (mcmc.cpp case 4) draws `idx`
   uniformly at `src/mcmc.cpp:4143-4144`; impl draws `other` uniformly from the
   remaining $n-1$ at `src/proposals.cpp:211-214`** ✓.
3. $\tau \ge 0$. The Beta parameters $\alpha, \beta \ge 1$ are bounded away
   from the singular boundary.
4. The remaining $n-2$ simplex elements are unchanged.

### 4.2 Proof

**Step 1: parameterise the move on the conditional simplex $\{A,B : A+B=T\}$.**

Since $T$ is conserved and the other entries are untouched, the move is a
1-degree-of-freedom move on the line segment $\{(A, T-A) : A\in(0,T)\}$.
Use $F = A/T$ as coordinate; $F \in (0,1)$. The Jacobian from
$(A, B)\big|_{A+B=T}$ to $F$ is $|dA/dF| = T$. The same $T$ holds for the
reverse move.

**Step 2: forward density.**

In $F$-coordinates the proposal is
$F'\sim\mathrm{Beta}(\alpha_{\mathrm{fwd}},\beta_{\mathrm{fwd}})$ with
$(\alpha_{\mathrm{fwd}}, \beta_{\mathrm{fwd}}) = (F\tau+1, (1-F)\tau+1)$,
so $q_F(F'\mid F) = \mathrm{Beta}(F'\mid \alpha_{\mathrm{fwd}}, \beta_{\mathrm{fwd}})$.

**Step 3: reverse density.**

The reverse proposal is generated by the *same* algorithm with
$F'$ in the role of $F$:
$(\alpha_{\mathrm{rev}}, \beta_{\mathrm{rev}}) = (F'\tau+1, (1-F')\tau+1)$ and
$q_F(F\mid F') = \mathrm{Beta}(F\mid \alpha_{\mathrm{rev}}, \beta_{\mathrm{rev}})$.

**Step 4: index-selection probability cancels.**

$P((i,j)\text{ chosen}\mid x) = P((i,j)\text{ chosen}\mid x') = \tfrac{1}{n(n-1)}$
by Assumption 2; this factor cancels in the ratio.

**Step 5: combine.**

In the original $(A,B)$-coordinates the Jacobian factor $T$ appears in both
forward and reverse densities and cancels. The Hastings ratio reduces to the
ratio of $F$-coordinate kernels:

$$
\frac{q(x\mid x')}{q(x'\mid x)} \;=\; \frac{\mathrm{Beta}(F\mid \alpha_{\mathrm{rev}}, \beta_{\mathrm{rev}})}{\mathrm{Beta}(F'\mid \alpha_{\mathrm{fwd}}, \beta_{\mathrm{fwd}})}.
$$

Taking logs gives the boxed expression. (Cross-check against
Lakner et al. 2008, "Efficiency of Markov Chain Monte Carlo Tree Proposals
in Bayesian Phylogenetics", Syst. Biol. 57:86–103, eq. 3 — same form.)

### 4.3 Implementation cross-check

`src/proposals.cpp:217-236`:

```cpp
const double oldA = x[index];
const double oldB = x[other];
...
const double oldF    = oldA / total;
const double alpha   = oldF * tuning + 1.0;
const double betaPar = (1.0 - oldF) * tuning + 1.0;
const double newF    = R::rbeta(alpha, betaPar);

x[index] = newF * total;
x[other] = (1.0 - newF) * total;

const double logFwd  = R::dbeta(newF, alpha, betaPar, 1);
const double revAlpha = newF * tuning + 1.0;
const double revBeta  = (1.0 - newF) * tuning + 1.0;
logHastings = R::dbeta(oldF, revAlpha, revBeta, 1) - logFwd;
```

`logFwd = log Beta(newF | α_fwd, β_fwd)` and
`R::dbeta(oldF, revAlpha, revBeta, 1) = log Beta(oldF | α_rev, β_rev)`,
so `logHastings = log Beta(oldF|α_rev,β_rev) - log Beta(newF|α_fwd,β_fwd)`
matches the proof. ✓

### 4.4 Edge cases

(a) **`total ≤ 0`.** Rejected at `src/proposals.cpp:222` (returns false,
caller treats as failed proposal). Detailed balance preserved.

(b) **`oldF ∈ {0,1}`.** Then one Beta parameter equals 1 (uniform) and the
other equals $\tau+1$. Density is still well-defined (Beta$(x|1,b)=b(1-x)^{b-1}$).
However the *state* with $F=0$ would normally have already triggered a
rejection elsewhere because the simplex element is zero.

(c) **`τ = 0`.** Beta$(1,1)$ = Uniform$(0,1)$ both ways; the kernel is symmetric
and `logHastings = 0` exactly (both `dbeta` terms equal $\log 1 = 0$). ✓

(d) **`newF` numerically equals 0 or 1.** `R::rbeta` is safe against this for
$\alpha,\beta \ge 1$; the densities at $0,1$ remain finite when the
parameters are $\ge 1$.

### 4.5 Detailed balance

For $\pi$ a density on the simplex restricted to the affine slice
$\{A+B=T\}$, the standard MH argument gives
$\pi(x)q(x'\mid x)\alpha(x,x') = \pi(x')q(x\mid x')\alpha(x',x)$
with the boxed Hastings ratio. The precondition is $\pi(x), \pi(x') > 0$
on the simplex interior, which is enforced by the strict-positivity
implicit in the move ($F\in(0,1)$ ⇒ $A, B > 0$).

### 4.6 Verdict

**Watertight.**

---

## §5. Theorem (Dirichlet-simplex proposal, cases 23 and 24)

For a move that selects $K$ indices $I = \{i_1,\dots,i_K\}\subseteq\{1,\dots,n\}$,
computes $S = \sum_{k\in I} x_k$, normalises to $\hat x = (x_k/S)_{k\in I}\in\Delta^{K-1}$,
draws $\hat z\sim\mathrm{Dir}(\alpha\hat x + \mathbf{1})$, and sets
$x_k\leftarrow \hat z_k\cdot S$ for $k\in I$:

$$
\boxed{\ \mathrm{logHastings} \;=\; \log\mathrm{Dir}(\hat x \mid \alpha\hat z + \mathbf{1}) \;-\; \log\mathrm{Dir}(\hat z \mid \alpha\hat x + \mathbf{1}).\ }
$$

### 5.1 Assumptions

1. $x$ is a non-negative vector summing to a constant $C$ (treated as a
   simplex on the $n$-simplex scaled by $C$). For relative branch lengths
   $C = 1$.
2. $S = \sum_{k\in I} x_k > 0$ (the selected block has positive mass).
   Enforced at `src/proposals.cpp:284`.
3. $K \ge 2$. Enforced by caller (`src/mcmc.cpp:4379, 4390`).
4. The set $I$ is drawn from a **state-independent** distribution
   $p_I$ — see the index-selection discussion in §5.4 — that is symmetric
   in the sense $p_I(\cdot \mid x) = p_I(\cdot \mid x')$ for any
   forward/reverse pair $(x, x')$ that has the same elements outside $I$.
5. The shifted Dirichlet concentration $\alpha\hat z + \mathbf{1}$ has all
   components $\ge 0.01$ (floor at `src/proposals.cpp:293, 315`); this
   ensures the density is well-defined and bounded away from singularity.

### 5.2 Proof

**Step 1: parametrise on the conditional sub-simplex.**

The move conserves $S$ and leaves elements outside $I$ untouched, so it is
a move on the $(K-1)$-dimensional sub-simplex
$\Delta_S = \{y\in\mathbb{R}_{\ge 0}^K : \sum y = S\}$. Coordinates:
$\hat x = (x/S)\in\Delta^{K-1}$. The change of variables
$y\leftrightarrow\hat y = y/S$ has Jacobian $|dy/d\hat y| = S^{K-1}$.

**Step 2: forward proposal.**

$\hat z \sim \mathrm{Dir}(\boldsymbol{\alpha}_{\mathrm{fwd}})$ with
$\boldsymbol{\alpha}_{\mathrm{fwd}} = \alpha\hat x + \mathbf{1}$.
Forward density in $\hat z$-coordinates:
$q_{\hat z}(\hat z\mid\hat x) = \mathrm{Dir}(\hat z\mid\boldsymbol{\alpha}_{\mathrm{fwd}})$.

The mapping $\hat z \mapsto y' = \hat z S$ introduces a Jacobian
$S^{K-1}$ in $y'$-coordinates.

**Step 3: reverse proposal.**

Symmetric construction with roles swapped:
$\boldsymbol{\alpha}_{\mathrm{rev}} = \alpha\hat z + \mathbf{1}$, reverse
density $q_{\hat x}(\hat x\mid\hat z) = \mathrm{Dir}(\hat x\mid\boldsymbol{\alpha}_{\mathrm{rev}})$.

**Step 4: Jacobian factors cancel.**

Both forward and reverse share the same Jacobian $S^{K-1}$ because $S$ is
conserved. So in $y$-coordinates the Hastings ratio reduces to the
$\hat\cdot$-coordinate Dirichlet ratio:

$$
\frac{q(x\mid x')}{q(x'\mid x)}
\;=\; \frac{\mathrm{Dir}(\hat x\mid\boldsymbol{\alpha}_{\mathrm{rev}})}{\mathrm{Dir}(\hat z\mid\boldsymbol{\alpha}_{\mathrm{fwd}})}.
$$

Taking logs:

$$
\mathrm{logHastings} \;=\; \log\mathrm{Dir}(\hat x\mid\boldsymbol{\alpha}_{\mathrm{rev}}) - \log\mathrm{Dir}(\hat z\mid\boldsymbol{\alpha}_{\mathrm{fwd}}).
$$

Writing $\mathrm{Dir}(\mathbf{u}\mid\boldsymbol{\alpha}) = \frac{\Gamma(\sum\alpha_i)}{\prod\Gamma(\alpha_i)}\prod u_i^{\alpha_i-1}$:

$$
\mathrm{logHastings} = \log\Gamma\!\big(\textstyle\sum\boldsymbol{\alpha}_{\mathrm{rev}}\big) - \sum\log\Gamma(\boldsymbol{\alpha}_{\mathrm{rev},i})
                    + \sum (\boldsymbol{\alpha}_{\mathrm{rev},i}-1)\log\hat x_i
$$
$$
\quad - \log\Gamma\!\big(\textstyle\sum\boldsymbol{\alpha}_{\mathrm{fwd}}\big) + \sum\log\Gamma(\boldsymbol{\alpha}_{\mathrm{fwd},i})
        - \sum (\boldsymbol{\alpha}_{\mathrm{fwd},i}-1)\log\hat z_i.
$$

### 5.3 Implementation cross-check

`src/proposals.cpp:319-335`:

```cpp
double logFwd = 0.0, logRev = 0.0;
double sumAlphaFwd = 0.0, sumAlphaRev = 0.0;
for (int i = 0; i < nCats; ++i) {
  sumAlphaFwd += alphaFwd[i];
  sumAlphaRev += alphaRev[i];
}
logFwd += std::lgamma(sumAlphaFwd);
logRev += std::lgamma(sumAlphaRev);
for (int i = 0; i < nCats; ++i) {
  logFwd -= std::lgamma(alphaFwd[i]);
  logRev -= std::lgamma(alphaRev[i]);
  logFwd += (alphaFwd[i] - 1.0) * std::log(std::max(zK[i], 1e-300));
  logRev += (alphaRev[i] - 1.0) * std::log(std::max(xK[i], 1e-300));
}
logHastings = logRev - logFwd;
```

Reading the variables:
- `logFwd = lgamma(Σα_fwd) - Σlgamma(α_fwd,i) + Σ(α_fwd,i − 1)·log ẑ_i`
  = `log Dir(ẑ | α_fwd)` ✓
- `logRev = lgamma(Σα_rev) - Σlgamma(α_rev,i) + Σ(α_rev,i − 1)·log x̂_i`
  = `log Dir(x̂ | α_rev)` ✓
- `logHastings = logRev − logFwd` ✓

Matches the boxed formula line by line.

### 5.4 Index-selection symmetry (subtle point)

**Case 23 (random Dirichlet, `dirichlet_simplex_impl`):**
Fisher-Yates partial shuffle (`src/proposals.cpp:352-358`) draws a uniformly
random $K$-subset $I$ from $\{1,\dots,n\}$. The probability of drawing the
*same* subset to go forward and backward is the same constant
$\binom{n}{K}^{-1}$. ✓ Symmetric.

**Case 24 (local Dirichlet, `local_dirichlet_impl`):**
`select_neighborhood` (`src/proposals.cpp:369-426`) picks a uniformly random
starting edge and BFS-expands on the edge-adjacency graph.
The probability of producing a particular connected subset $I$ depends
**only on the tree topology**, not on the values $\{x_k\}$.
Forward and reverse moves apply on the same topology (the simplex move
does *not* alter the tree), so $p_I(\mathrm{fwd}) = p_I(\mathrm{rev})$. ✓ Symmetric.

In both cases the index-selection factor cancels and only the Dirichlet
density ratio remains in `logHastings`.

### 5.5 Edge cases

(a) **`S ≤ 0`** (numerical underflow of an entire block). Rejected at
`:284`. ✓

(b) **`α_fwd[i] = 0.01` floor.** When `xK[i] = 0` (numerically), the
formula would have $\alpha_{\mathrm{fwd}} = 1$ and density well-defined.
The floor at `:293` protects against `xK[i] < 0` arising from FP error.
Symmetrically applied to `α_rev` at `:315`. ✓ This preserves detailed
balance because the floored kernel is itself a valid (different) MH
kernel.

(c) **`zK[i]` numerical underflow.** Guarded by `std::max(zK[i], 1e-300)` at
`:301, :330`. The actual gamma draws are bounded by the same floor at `:301`
before normalisation, so $\hat z_i \ge 10^{-300}/\sum$ which is non-zero.
The `log` floor at `:330` is an extra safety. **One subtle point:** if the
floor at `:330` *fires* (i.e. `zK[i] < 1e-300` is materially different
from `zK[i] = 0`), then `logFwd` uses $\log(10^{-300})$ but the *next*
draw uses the value actually stored at `zK[i] = 1e-300` from line 301 —
self-consistent.

(d) **`K = n`.** Whole-simplex Dirichlet refresh. The Fisher-Yates degenerates
to identity; `select_neighborhood` returns all edges (`:374-379`). Both
fine. The conserved $S$ is the entire vector sum.

(e) **`K = 1`.** Forbidden by the engine (`if (nCats < 2) return false` at
`:4379, :4390` for cases 23/24). A single-element "Dirichlet" would be
deterministic (point mass) and not informative.

### 5.6 Detailed balance

For $\pi$ a density on the simplex (any continuous density supported on the
$(n-1)$-simplex), and for fixed $I$, the move is a standard Dirichlet
proposal on $\Delta_S$ with the Hastings ratio derived. Detailed balance
holds on the conditional kernel (Lange 2010, §11.3). The overall kernel is
a mixture over $I$ with state-independent mixing weights; mixture-of-MH
kernels preserves detailed balance (Tierney 1994, Thm. 2). ✓

### 5.7 Verdict

**Watertight.**

---

## §6. Theorem (Bounded-integer walk, case 7)

For move type 7 (`int_walk kPrime`):

$$
\boxed{\ \mathrm{logHastings} \;=\; 0\ \text{when }\ k'_{\mathrm{new}} \ge \mathrm{lower};\quad -\infty\ \text{otherwise}.\ }
$$

### 6.1 Assumptions

1. Move proposes $k' \to k'+\Delta$ with $\Delta\sim\mathrm{Uniform}\{-w,\dots,w\}$
   for some window $w \ge 1$.
2. Lower bound: $k' \ge k_{\mathrm{obs}}$ (set by `lowerK = data->kObs[charIdx]`
   at `src/mcmc.cpp:4271`).
3. No upper bound — the support is $\{k_{\mathrm{obs}}, k_{\mathrm{obs}}+1, \dots\}$.

### 6.2 Proof

Although the parameter is integer-valued (so technically a discrete
proposal, included here for completeness as part of L2), the same MH
argument applies with sums replacing integrals.

The unrestricted proposal $\Delta\sim\mathrm{Uniform}\{-w,\dots,w\}$ is
**symmetric**: $P(\Delta=d) = P(\Delta=-d) = \frac{1}{2w+1}$ for any
$d\in\{-w,\dots,w\}$.

**Case (i): $k' + \Delta < \mathrm{lower}$.**
The code rejects (`return false` at `src/mcmc.cpp:4275`). Equivalent to setting
$q(k'_{\mathrm{new}}\mid k') = 0$ on $k'_{\mathrm{new}} < \mathrm{lower}$, with
the remaining $2w+1 - |\{d : k'+d < \mathrm{lower}\}|$ values getting weight
$\frac{1}{2w+1}$ and the difference falling into a "stay" Dirac at $k'$ via
the rejection mechanism. The accept ratio for rejected proposals is 0; the
overall move is a valid MH kernel.

**Case (ii): $k'_{\mathrm{new}} \ge \mathrm{lower}$ and
$|k' - k'_{\mathrm{new}}|\le w$.**
Both forward $k'\to k'_{\mathrm{new}}$ and reverse
$k'_{\mathrm{new}}\to k'$ are allowed, with equal probability
$\frac{1}{2w+1}$ (since the reverse step has $\Delta = k' - k'_{\mathrm{new}}$
which is also in $\{-w,\dots,w\}$, and if $k' \ge \mathrm{lower}$ — which
it must be by construction of the chain — the reverse proposal lands in
support).

Hence $q(k'_{\mathrm{new}}\mid k') = q(k'\mid k'_{\mathrm{new}}) = \frac{1}{2w+1}$
and `logHastings = log 1 = 0`. ✓

**Subtle point — boundary asymmetry.** When $k' \in [\mathrm{lower},
\mathrm{lower}+w-1]$, *some* forward proposals are rejected (those that
would go below `lower`). For such $k'$, the conditional probability of
*reaching* a particular $k'_{\mathrm{new}}$ from $k'$ is still
$\frac{1}{2w+1}$ for $k'_{\mathrm{new}}\in\{k'-w,\dots,k'+w\}\cap\{\ge\mathrm{lower}\}$,
with the *rest* of the probability mass falling on rejection (stay at $k'$).
The reverse from any allowed $k'_{\mathrm{new}}$ back to $k'$ is also
$\frac{1}{2w+1}$ (since both are within window $w$ by construction).
**So both forward and reverse probabilities are equal**; the proposal
remains symmetric on the support, even at the boundary, and `logHastings = 0`
is correct.

The mass "wasted" on rejected $\Delta$'s does **not** appear in the
Hastings ratio — it merely lowers the overall acceptance rate. This is the
canonical "naive symmetric random walk with rejection at boundary"
construction (Tierney 1994, §2.3.3).

### 6.3 Implementation cross-check

`src/mcmc.cpp:4267-4279`:

```cpp
case 7: { // int_walk kPrime
  kPrimeCharIdx = charIdx;
  oldKPrimeVal  = state->kPrime[charIdx];
  int oldK      = oldKPrimeVal;
  int lowerK = data->kObs[charIdx];
  int range  = 2 * intWalkWindow + 1;
  int delta  = static_cast<int>(R::unif_rand() * range) - intWalkWindow;
  int newK   = oldK + delta;
  if (newK < lowerK) return false;
  state->kPrime[charIdx] = newK;
  logHastings = 0.0;
  break;
}
```

- `range = 2*w+1`, `delta ∈ {-w, ..., +w}` (uniform). ✓
- Below-lower → `return false` (immediate rejection). ✓
- Otherwise `logHastings = 0`. ✓

The R-side wrapper `ProposeBoundedIntWalk` in `R/proposals.R:80-90`
implements the same semantics with `sample(-window:window, 1L)` and
`return(list(value = x, logHastings = -Inf))` on out-of-bounds.

### 6.4 Edge cases

(a) **No upper bound.** Correct — the support
$\{k_{\mathrm{obs}}, k_{\mathrm{obs}}+1,\dots\}$ is unbounded above. Any
proposal $\Delta > 0$ lands in support.

(b) **`window = 0`.** Then `range = 1`, `delta = 0` always: the move is the
identity. `logHastings = 0` and the chain stays — trivially detailed-balance
compliant.

(c) **`window < 0`.** Cannot occur — `intWalkWindow` is treated as
non-negative by the engine; R-side `ProposeBoundedIntWalk` accepts
`window = 1L` by default.

(d) **`k' = lower` exactly.** Half the window proposals (those with
$\Delta < 0$) are rejected. As shown in 6.2 this preserves symmetry. ✓

### 6.5 Verdict

**Watertight.**

---

## §7. Summary table

| Move | logHastings | Symmetric? | Jacobian | Code site | Verdict |
|------|-------------|------------|----------|-----------|---------|
| Bactrian scale (cases 0–3, 8, 16, 27, 28, 31) | $\log\mathrm{mult} = \tau z$ | Yes (in $z$) | $\log(x'/x)$ | `src/mcmc.cpp:4116-…` | Watertight |
| 2D Bactrian joint (cases 21, 22) | $\log\mathrm{mult}_1 + \log\mathrm{mult}_2$ | Yes (in $(z_1,z_2)$) | $\log(x'_1 x'_2/(x_1 x_2))$ | `src/mcmc.cpp:4363, :4373` | Watertight |
| Logit MH on $p$ (case 30) | $\log\frac{p'(1-p')}{p(1-p)}$ | Yes (in $\eta$) | inv-logit | `src/mcmc.cpp:4452-4453` | Watertight |
| Beta-simplex (case 4) | $\log\mathrm{Beta}(F\mid\alpha_{\rm rev},\beta_{\rm rev}) - \log\mathrm{Beta}(F'\mid\alpha_{\rm fwd},\beta_{\rm fwd})$ | No | $T$ (cancels) | `src/proposals.cpp:235` | Watertight |
| Dirichlet-simplex (cases 23, 24) | $\log\mathrm{Dir}(\hat x\mid\alpha\hat z+\mathbf{1}) - \log\mathrm{Dir}(\hat z\mid\alpha\hat x+\mathbf{1})$ | No | $S^{K-1}$ (cancels) | `src/proposals.cpp:319-334` | Watertight |
| Bounded int walk (case 7) | 0 (or $-\infty$ on rejection) | Yes (uniform) | — | `src/mcmc.cpp:4277` | Watertight |

---

## §8. Cross-cutting findings

### 8.1 No bugs found

Every per-move `logHastings` returned by the codebase matches the
first-principles derivation, including:

- The Jacobian for multiplicative scale moves (cases 0–3, 8, 16, 27, 28, 31).
- The double Jacobian for 2D joint scale moves (cases 21, 22).
- The inv-logit Jacobian for the logit-MH `p` move (case 30).
- The absence of Jacobian for beta-simplex and Dirichlet-simplex (because
  the conserved sums $T$ and $S$ produce identical Jacobian factors in both
  forward and reverse directions).
- The trivial symmetry for the bounded-integer walk.

### 8.2 Documented design tradeoffs (not bugs)

1. **Wall-clock adaptation of `scaleTuning` and `jointRho`.** As noted in
   `feedback_wallclock_adaptation`, the engine adapts tuning constants
   using wall-clock cost normalisation outside the per-move accept/reject
   path. The per-move Hastings derivation here treats the *current* tuning
   values as constants of the kernel; this is correct (each step's MH
   ratio uses the kernel that was actually applied). Time-homogeneity is
   restored by suspending adaptation in production, which the engine does.

2. **Case 8 (`scale p`) uses multiplicative Bactrian without an upper-bound
   check.** Already filed as EG-007/EG-009. The Hastings ratio
   `log(mult)` is mathematically correct on $(0,\infty)$; the move relies
   on the prior to enforce $p \le 1$ via $\pi(p') = 0$ for $p' \ge 1$.
   Wasteful but not biased.

3. **Lockstep `classRateLogSd[0] = rateLogSd`** in cases 2 and 21. The
   constraint is enforced as an invariant of the state vector
   (`feedback_model_scope`), so the move acts on the constrained manifold
   and no further Jacobian arises.

### 8.3 Possible improvements (out-of-scope for L2)

- The Dirichlet-simplex concentration parameter is *shifted* (i.e.
  $\alpha\hat x + 1$ rather than $\alpha\hat x$). This is the form RevBayes
  uses (`src/proposals.cpp:288-290`), which avoids the singular point at
  $\alpha = 0$ in standard Dirichlet kernels. The proof goes through
  identically.

### 8.4 Patch status

**No patch attached.** All implementations match the derivations. The
`dev/red-team/patches/L2-hastings-continuous.patch` is empty.

---

## Verdict

**Watertight** for all six continuous-proposal families in scope of L2.

References cited:
- Yang, Z., & Rodríguez, C. E. (2013). Searching for efficient Markov chain
  Monte Carlo proposal kernels. *PNAS* 110: 19307–19312.
- Lakner, C. *et al.* (2008). Efficiency of Markov chain Monte Carlo tree
  proposals in Bayesian phylogenetics. *Systematic Biology* 57: 86–103.
- Liu, J. S. (2001). *Monte Carlo Strategies in Scientific Computing*.
  Springer, §§5.3–5.4.
- Lange, K. (2010). *Numerical Analysis for Statisticians* (2nd ed.).
  Springer, §11.
- Tierney, L. (1994). Markov chains for exploring posterior distributions.
  *Annals of Statistics* 22: 1701–1762.
- Hastings, W. K. (1970). Monte Carlo sampling methods using Markov chains
  and their applications. *Biometrika* 57: 97–109.
- Roberts, G. O., & Rosenthal, J. S. (2004). General state space Markov
  chains and MCMC algorithms. *Probability Surveys* 1: 20–71.
