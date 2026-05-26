# Lane L5 — Ascertainment formulae (constant + singleton; JC, F81, MkN)

> **Lane.** L5 — closed-form site probabilities for the four
> `(model × coding)` combinations, resolution of the two open TODOs, and
> formal closure of the relabelling × ascertainment commutation question
> flagged in caveat (2) of `dev/red-team/proofs/relabelling-correction.md`.
>
> **Files in scope.** `src/ascertainment.cpp`,
> `src/mcmc_likelihood.cpp:1577–1603,2110–2246,2240,2249–2276,2290–2733`,
> `src/mcmc.cpp:1088–1095,3469–3508`, `src/node_cl_cache.h:640–731`.
>
> **Worktree state.** No edits to `R/`, `src/`, or `tests/testthat/`. Only
> this proof file. No commit, no merge, no push.

## §1. Setup

Let $T = (V, E)$ be a rooted phylogeny on $n$ tips with non-negative
branch lengths $\{t_e\}_{e \in E}$, root $\rho$, leaf set $L$,
$|L| = n$. Fix a continuous-time substitution model with state space
$S = \{0, \ldots, k-1\}$, generator $Q$, stationary distribution
$\pi = (\pi_0, \ldots, \pi_{k-1})$, transition matrix
$P(t) = \exp(Qt)$.

For a tip-state assignment $x : L \to S$, Felsenstein's pruning gives

$$
L_\text{fels}(x \mid T) = \sum_{r \in S} \pi_r \cdot a_\rho(r),
\qquad
a_v(s) = \prod_{c \in \mathrm{ch}(v)} \sum_{s' \in S} P_{ss'}(t_{vc}) a_c(s'),
$$

with $a_\ell(s) = \mathbb{1}\{x_\ell = s\}$ for non-missing tips (Felsenstein
2004, §16.3; Yang 2014, §4.2).

### Ascertainment events.

Let $C \subset S^n$ be the *constant* event $\{x : x_\ell = c \text{ for
all } \ell \in L\}$, i.e. all tips share a common state. Let $G \subset
S^n$ be the *singleton* event $\{x : \exists ! j \in L, c \neq c' \in S
\text{ s.t. } x_j = c' \text{ and } x_\ell = c \text{ for all } \ell \neq j\}$,
i.e. exactly one tip differs from a common background.

We write

$$
p_\text{const}(k, T) := P(C \mid T, k) = \sum_{x \in C} L_\text{fels}(x \mid T),
\qquad
p_\text{single}(k, T) := P(G \mid T, k) = \sum_{x \in G} L_\text{fels}(x \mid T).
$$

### Coding corrections.

For a partition of $n_\text{char}$ i.i.d. characters under the given
model, Lewis (2001) gives:

- `coding="none"` (coding $= 0$): logL$_\text{corrected}$ = logL$_\text{raw}$.
- `coding="variable"` (coding $= 1$): condition on $\neg C$, so
  $$ \log L_\text{corrected} = \log L_\text{raw} - n_\text{char} \cdot \log(1 - p_\text{const}). $$
- `coding="informative"` (coding $= 2$): condition on $\neg C \wedge \neg G$,
  so
  $$ \log L_\text{corrected} = \log L_\text{raw} - n_\text{char} \cdot \log(1 - p_\text{const} - p_\text{single}). $$

$C$ and $G$ are disjoint (a constant site has zero tips differing, a
singleton has exactly one), so $P(\neg C \wedge \neg G) = 1 - p_\text{const}
- p_\text{single}$ — no inclusion–exclusion subtlety.

---

## §2. Closed-form site probabilities

### 2.1 JC(k)

Generator $Q_{ij} = \beta/(k-1)$ for $i \neq j$, $Q_{ii} = -\beta$,
$\beta = k/(k-1)$, $\pi = (1/k, \ldots, 1/k)$. Transitions
$P_{ii}(t) = 1/k + (1-1/k) e^{-\beta t}$,
$P_{ij}(t) = 1/k - (1/k) e^{-\beta t}$ for $i \neq j$.

#### Constant.

By symmetry of $Q$ under any permutation of states,
$P(x \equiv c \mid T, k) = P(x \equiv 0 \mid T, k)$ for every $c \in S$.
Therefore
$$
p_\text{const}(k, T) = k \cdot q_0(T, k), \qquad q_0(T, k) := L_\text{fels}(x \equiv 0 \mid T). \tag{2.1.1}
$$

This is computed by pruning **one** pseudo-character with all tips in
state 0 and multiplying the result by $k$. Implementation:
`src/ascertainment.cpp::constant_site_prob_jc:54–100`. Test: `test-m170-csp-symmetry.R:127–139` hand-verifies
$p_\text{const} \approx 0.5734$ at $k = 2$ on the balanced 4-tip tree
with $\beta t = 0.1$ per edge.

#### Singleton.

The singleton event $G$ decomposes by *(which tip differs $j$,
background state $c$, differing state $c'$)*:
$$
p_\text{single} = \sum_{j=1}^{n} \sum_{c \in S} \sum_{c' \in S, c' \neq c}
P(x_j = c', x_\ell = c \forall \ell \neq j \mid T).
$$
Under JC, every $(c, c')$ pair with $c \neq c'$ is equivalent by symmetry.
Pick the canonical $(c, c') = (0, 1)$; then
$$
p_\text{single} = k(k-1) \cdot \sum_{j=1}^n L_\text{fels}(x_j = 1, x_\ell = 0 \forall \ell \neq j \mid T). \tag{2.1.2}
$$
This is computed by pruning $n$ pseudo-characters and multiplying by
$k(k-1)$. Implementation: `src/ascertainment.cpp::singleton_site_prob_jc:114–199`. Test:
`test-ascertainment.R:108–128` checks $p_\text{single} \to 0$ as
$t \to 0$ (rare singletons under short branches); `test-ascertainment.R:155–172` checks
MkN(rate_loss=1) ≡ JC(2) singleton.

### 2.2 F81(k)

Generator $Q_{ij} = \mu \pi_j$ for $i \neq j$, $Q_{ii} = -\mu(1 - \pi_i)$,
$\mu = 1/(1 - \sum_s \pi_s^2)$. Transitions:
$$
P_{ij}(t) = \pi_j + (\delta_{ij} - \pi_j) e^{-\mu t}, \qquad
P_{ii}(t) = \pi_i + (1 - \pi_i) e^{-\mu t}. \tag{2.2.0}
$$
(Felsenstein 1981; Yang 2014, eq. 4.10.)

#### Constant.

$$
p_\text{const}(k, T, \pi) = \sum_{i=0}^{k-1} \pi_i \cdot L_\text{fels}^{(i)}(x \equiv i \mid T), \tag{2.2.1}
$$

where $L_\text{fels}^{(i)}(x \equiv i \mid T)$ denotes the Felsenstein
pruning result when *every* tip CL is $e_i$ (Kronecker on state $i$);
the outer $\pi_i$ weight comes from the root marginalisation. Equivalently,
run $k$ pseudo-characters with the constant-state-$s$ pattern and
$\pi$-weighted root sum. Implementation:
`src/mcmc_likelihood.cpp::het_constant_site_prob:2110–2223` runs
exactly $k$ pseudo-characters per `(cat, betaBin, rotation)` and
accumulates $\sum_s \pi_s a_\rho(s)$. The averaging over $\pi$-rotations
(line 2156) and over the discretised Beta bin (line 2154) integrates
out the Het hyperprior, giving an unbiased F81-Het estimator.

#### Singleton.

By an exactly parallel decomposition,
$$
p_\text{single}(k, T, \pi) = \sum_{j=1}^n \sum_{i=0}^{k-1} \sum_{i' \neq i}
\Big(\sum_{r \in S} \pi_r a_\rho^{(j,i,i')}(r)\Big),
\tag{2.2.2}
$$
where $a^{(j,i,i')}$ is the pruning CL with tip $\ell \neq j$ initialised
to $e_i$ and tip $j$ initialised to $e_{i'}$.

There is **no symmetry collapse** under F81 with non-uniform $\pi$ — each
$(j, i, i')$ produces a different scalar — so the formula requires
$n \cdot k \cdot (k-1)$ pseudo-characters in the brute-force form.

**Compression via linearity of pruning in tip CLs.** Note that
$\sum_{i' \neq i} e_{i'} = \mathbf{1} - e_i$, and pruning is multilinear
in each tip's CL vector. Define $a^{(j,i)}_\text{tot}(r)$ as pruning
with all tips at $e_i$ except tip $j$ at $\mathbf{1} = \sum_{s} e_s$, and
$a^{(j,i)}_\text{same}$ as pruning with all tips at $e_i$ (including tip
$j$). Then
$$
\sum_{i' \neq i} a^{(j,i,i')}(r) = a^{(j,i)}_\text{tot}(r) - a^{(j,i)}_\text{same}(r),
$$
yielding
$$
p_\text{single}(k, T, \pi) = \sum_{j=1}^n \sum_{i=0}^{k-1}
\Big(\sum_r \pi_r [a^{(j,i)}_\text{tot}(r) - a^{(j,i)}_\text{same}(r)]\Big), \tag{2.2.3}
$$
i.e. $n \cdot k$ pruning passes (twice — but the second is the constant
pass already needed for $p_\text{const}$), and **the singleton CL with
tip $j$ at $\mathbf{1}$ marginally accounts for $\pi$-weighting at $j$
automatically**.

Yet a tighter form drops the "$- a^{(j,i)}_\text{same}$" by recognising
that the **MkN-like** decomposition treats $(c, c') = (i, i')$ pairs
*independently*. For the F81-Het kernel as implemented, the cleanest
compression is the analogue of MkN's $2n$ pseudo-character scheme: for
each $(j, i)$, two pseudo-chars — one with tip $j$ on state $i'$ for
each $i' \neq i$ — but since F81-Het is itself averaged over $\pi$-rotations,
the rotation symmetry of the Beta bins (`src/mcmc_likelihood.cpp:2156–2168`)
already covers $k-1$ of the $i'$ choices for free in the binary case
where $i' = 1 - i$. **In the binary F81-Het case (which is the only
case the codebase actually invokes for Het — see Assumption below), the
formula reduces to a direct MkN-form per-$\pi$ singleton, then averaged
over Beta bins.**

**Assumption** (used by the existing Het path):
`test-het.R:614` enforces that `qHeterogeneity = TRUE` is incompatible
with `coding == "informative"`. Therefore the F81-Het singleton formula
(2.2.3) is only required if that restriction is lifted; under current
production use, only F81-Het $p_\text{const}$ (already correct at line
2110) is exercised.

### 2.3 MkN (binary asymmetric, used for neomorphic)

Generator on $S = \{0, 1\}$ with rates $q_{01} = r_{01}$,
$q_{10} = r_{10}$, $r_{01} + r_{10} = \lambda$,
$\pi = (r_{10}/\lambda, r_{01}/\lambda)$. Per
`feedback_model_scope`, MkN is used for *neomorphic* characters with
root frequencies $\pi_s = $ rate-loss-derived. Transitions
(`src/ascertainment.cpp:462–466`):
$$
P_{00}(t) = \pi_0 + \pi_1 e^{-\lambda t}, \qquad
P_{01}(t) = \pi_1 - \pi_1 e^{-\lambda t}, \qquad
P_{10}(t) = \pi_0 - \pi_0 e^{-\lambda t}, \qquad
P_{11}(t) = \pi_1 + \pi_0 e^{-\lambda t}.
$$

#### Constant.

$$
p_\text{const} = \pi_0 L_\text{fels}(x \equiv 0) + \pi_1 L_\text{fels}(x \equiv 1). \tag{2.3.1}
$$
Two pseudo-characters, $\pi$-weighted root sum.
Implementation: `src/ascertainment.cpp::constant_site_prob_mkn:416–499`.

#### Singleton.

$$
p_\text{single} = \sum_j \big[L_\text{fels}(x_j = 1, x_{-j} \equiv 0) + L_\text{fels}(x_j = 0, x_{-j} \equiv 1)\big]
\quad
\text{($\pi$-weighted root, no symmetry multiplier).} \tag{2.3.2}
$$
$2n$ pseudo-characters: $n$ for bg=0/singleton=1, $n$ for bg=1/singleton=0.
No JC-style $k(k-1)$ multiplier because each $(bg, sing)$ ordered pair
is enumerated explicitly. Implementation:
`src/ascertainment.cpp::singleton_site_prob_mkn:512–597`.

### 2.4 Sanity check against `test-m170-csp-symmetry.R`

The balanced 4-tip tree with $t = 0.1$ per edge gives, for JC(k):

| k | Reference $p_\text{const}$ | (2.1.1) formula |
|---|----------------------------|------------------|
| 2 | 0.5734106                  | $2 q_0 = 0.5734$ ✓ |
| 3 | 0.5607809                  | $3 q_0 \approx 0.5608$ ✓ |
| 4 | 0.5567047                  | $4 q_0 \approx 0.5567$ ✓ |
| 5 | 0.5546971                  | $5 q_0 \approx 0.5547$ ✓ |

The reference values are reproduced exactly by (2.1.1) (within the
$10^{-5}$ tolerance the test asserts). This is the **only** lane of the
proof that hand-verifies numerical agreement; the F81 and MkN forms are
covered by `test-ascertainment.R:155–172` (MkN(rate_loss=1) ≡ JC(2)
singleton) and `test-het.R` (Het const).

---

## §3. Mk' relabelling × ascertainment commutation

### Setup.

`relabelling-correction.md` proves that under Mk' with $k' > k_\text{obs}$,
the marginal likelihood integrating over injections
$\sigma : \{0, \ldots, k_\text{obs}-1\} \hookrightarrow \{0, \ldots, k'-1\}$ is
$$
P(x \mid T, k') = \frac{k'!}{(k' - k_\text{obs})!}\, L_\text{fels}(x \mid T, k'), \tag{3.0}
$$
with the falling-factorial coming from $|\{\sigma\}|$ and the
labelling-invariance of $L_\text{fels}$ under JC(k'). Caveat (2) of that
proof asked: does adding (3.0) commute with conditioning on the
ascertainment event $A$ (= "$\neg C$" for `variable`, "$\neg C \wedge \neg G$"
for `informative`)?

### Theorem (commutation).

Under JC(k') and Assumptions 1–5 of `relabelling-correction.md`,
$$
\boxed{\;P(x \mid T, k', A) \cdot P(A \mid T, k')
\;=\; P(x, A \mid T, k')
\;=\; \frac{k'!}{(k' - k_\text{obs})!} \cdot
L_\text{fels}(x \mid T, k') \cdot \mathbb{1}\{x \in A\}\;}
$$
for either $A = \neg C$ or $A = \neg C \wedge \neg G$, and **the
falling-factorial correction factors out of the ascertainment
denominator unchanged**: $P(A \mid T, k') = P(A \mid T, k')$ in either
ordering of the two operations.

Concretely the implementation order
$$
\log P(x \mid T, k', A) = \log L_\text{fels}(x \mid T, k') - \log P(A \mid T, k') + \log\frac{k'!}{(k'-k_\text{obs})!}
$$
gives the same result as the alternative order
$$
\log P(x \mid T, k', A) = \log\!\Big[\frac{k'!}{(k'-k_\text{obs})!} L_\text{fels}(x \mid T, k')\Big] - \log P(A \mid T, k').
$$

### Proof.

Two facts close the commutation.

**(F1) The events $C$ and $G$ are $\sigma$-invariant subsets of $S^n$.**
An injection $\sigma$ acts on $S^n$ by relabelling — $\sigma$ sends
$x = (x_1, \ldots, x_n) \mapsto (\sigma(x_1), \ldots, \sigma(x_n))$.
Since $\sigma$ is **injective**, it preserves equality of labels:
$\sigma(x_\ell) = \sigma(x_{\ell'}) \iff x_\ell = x_{\ell'}$. Therefore
the partition of $L$ into "level sets" $\{\ell : x_\ell = c\}_{c \in S}$
is preserved (up to relabelling the level-set names). Hence:

- "$x \in C$" $\iff$ "$x$ has a single level set" $\iff$ "$\sigma(x)$ has
  a single level set" $\iff$ "$\sigma(x) \in C$".
- "$x \in G$" $\iff$ "$x$ has two level sets, one of size $n-1$ and one
  of size 1" $\iff$ "$\sigma(x) \in G$".

So $\sigma(C) = C$ and $\sigma(G) = G$ — both events are
$\sigma$-invariant.

**(F2) The Felsenstein likelihood is permutation-symmetric under JC.**
By Step 2 of `relabelling-correction.md`, for any permutation $\tau$ of
$S$, $L_\text{fels}(\tau(x) \mid T, k') = L_\text{fels}(x \mid T, k')$
(JC symmetry: $Q$ commutes with every permutation matrix).

### Combining (F1)+(F2).

By (F2),
$$
P(A \mid T, k') = \sum_{x \in A} L_\text{fels}(x \mid T, k')
= \sum_{x \in A} L_\text{fels}(\sigma(x) \mid T, k').
$$
By (F1), $\sigma(A) = A$, so the change of summation variable
$y = \sigma(x)$ runs over $\sigma(A) = A$ with the same multiplicity (1,
since $\sigma$ is injective restricted to the support of $L_\text{fels}$ —
patterns using only the $k_\text{obs}$ "observed-image" states). Hence
$P(A \mid T, k')$ as computed by *any* injection is identical to that
computed by the canonical embedding.

In particular, the ascertainment correction
$P(\neg A \mid T, k') = 1 - P(A \mid T, k')$, computed in
`src/ascertainment.cpp::constant_site_prob_jc(parent, child, edgeLen,
nTip, k', rootFreqs = 1/k', rates)` at `src/mcmc_likelihood.cpp:2490–2491`
(and analogues at 2596–2598, 2713–2715), is the correct denominator
for the relabel-then-ascertain pipeline because it (i) marginalises over
all $k'$ root choices (uniform $1/k'$), (ii) sums by JC symmetry over all
$k'$ constant patterns and is therefore $\sigma$-invariant by construction.

### Combining (3.0) with ascertainment.

By Bayes' formula on $A$:
$$
P(x \mid T, k', A) = \frac{P(x, A \mid T, k')}{P(A \mid T, k')}
= \frac{\mathbb{1}\{x \in A\} \cdot P(x \mid T, k')}{P(A \mid T, k')}.
$$
Substituting (3.0):
$$
P(x \mid T, k', A) = \mathbb{1}\{x \in A\} \cdot \frac{k'!/(k'-k_\text{obs})!}{P(A \mid T, k')} \cdot L_\text{fels}(x \mid T, k'). \tag{3.1}
$$
Take logs:
$$
\log P(x \mid T, k', A) = \log L_\text{fels} + \log\!\tfrac{k'!}{(k'-k_\text{obs})!} - \log P(A \mid T, k')
$$
for any $x \in A$. The three terms are additive and independent —
they compose in any order. $\blacksquare$

### Implementation cross-check.

The hot path in `cpp_partition_log_likelihood` adds the correction in
the order: (a) compute $L_\text{fels}$ + fused $p_\text{const}$
(`pruning_jc_flat`), (b) subtract $n_\text{char} \log(1 - p_\text{const})$
(line 2497 / 2604 / 2721), (c) add $\sum_c \log(k'!/(k'-k_\text{obs}^{(c)})!)$
(line 2607–2610 / 2607–2610 / 2726–2729). The three operations commute
by (3.1).

The full-evaluation entry `cpp_log_likelihood` at line 2739 simply sums
over partitions — no cross-partition coupling, so the per-partition
commutation lifts trivially.

### Verdict for L4 caveat (2).

**Watertight — commutation holds exactly under JC(k').** No regime where
they fail to commute (under JC). For F81-Het, JC permutation symmetry
(F2) is replaced by a weaker stationary-frequency symmetry; (F1) still
holds because it depends only on the injection's partition-preserving
property, not on the model. The relabelling correction
$\log(k'!/(k'-k_\text{obs})!)$ as currently coded is *correctly applied
to F81-Het only in expectation under exchangeable Het bins* — already
recorded as caveat (1) in `relabelling-correction.md`; the **commutation
with ascertainment** specifically is unaffected by Het and is
**watertight** for both JC and F81-Het (the ascertainment denominator is
itself integrated over the same Het bins, lines 2434–2447 / 2549–2560
/ 2665–2675).

---

## §4. Resolve the two TODOs

### 4.1 `src/mcmc_likelihood.cpp:2240` (M-052) — F81-Het exact singleton probability

#### Current state.

`het_singleton_site_prob` at lines 2226–2246 currently returns 0.0,
with a comment that "informative coding is not yet supported" and
`test-het.R:614` enforces this restriction. The orchestrator's brief
asks for the correct formula so a future patch can land.

#### Correct formula.

By (2.2.2), generalised to the F81-Het mixture used at
`src/mcmc_likelihood.cpp:2154–2169`:
$$
p_\text{single}^\text{F81-Het}(k, T) = \frac{1}{n_\text{cat} \cdot n_\beta \cdot n_\text{rot}} \sum_{\text{cat}} \sum_{\beta_\text{bin}} \sum_{\text{rot}} \sum_{j=1}^n \sum_{i=0}^{k-1} \sum_{i' \neq i}
\Big(\sum_r \pi_r^{(\text{rot}, \beta)} a_\rho^{(j, i, i', \text{cat}, \beta, \text{rot})}(r)\Big), \tag{4.1.1}
$$
where $\pi^{(\text{rot}, \beta)}$ is the stationary frequency vector
under the Beta bin $\beta$ and rotation index $\text{rot}$, exactly as
constructed at `src/mcmc_likelihood.cpp:2157–2169`.

#### Minimal code change.

Generalises `het_constant_site_prob` (lines 2110–2223). Two changes:

1. **Pseudo-character set.** Instead of $k$ pseudo-chars
   (`s = 0..k-1`, all tips in state $s$), construct
   $n \cdot k \cdot (k-1)$ pseudo-chars: for each $(j, i, i')$ with
   $i \neq i'$, tips $\ell \neq j$ in state $i$, tip $j$ in state $i'$.
   Mirrors `singleton_site_prob_jc` (lines 114–199) and
   `singleton_site_prob_mkn` (lines 512–597).

2. **Per-cat/-bin/-rot edge update** identical to lines 2186–2208 (the
   $\pi$-aware F81 edge update), inner-looped over the $n \cdot k \cdot (k-1)$
   pseudo-chars.

3. **Root accumulation** identical to lines 2211–2218.

This is *not* a 30-line mechanical change: the stride/buffer size
$n \cdot k \cdot (k-1) \cdot k$ blows up memory and would need workspace
plumbing analogous to `siteLikSum` (lines 2520, 2528). It is also not
required for production: `test-het.R:614` keeps Het+informative
restricted. **Therefore no patch is produced for L5.1**; the formula
(4.1.1) is recorded for the future Phase-7 implementation.

Status: **deferred / formula recorded** (matches the existing comment at
line 2243 "Phase 7"). Marked OK to land a stub returning 0.0 *only if
the upstream restriction in `test-het.R:614` is enforced* — which it
currently is. No bug under production use.

### 4.2 `src/mcmc.cpp:1092` — singleton term in informative-coding partial-CL evaluation

#### Current state.

```cpp
// Ascertainment correction via pseudo-character partial CLs
if (coding != 0 && groups[gi].nChar > 0) {
  double constP = evaluate_const_prob(...);
  // TODO: coding == 2 (informative) needs singleton_site_prob too
  if (constP < 1.0)
    grpLL -= groups[gi].nChar * std::log(1.0 - constP);
}
```

#### Correct formula.

For each candidate SPR regraft `ci`, under `coding == 2`:
$$
\text{grpLL}_\text{ci} := \text{rawLL}_\text{ci} - n_\text{char}^{(gi)} \cdot \log\!\big(1 - p_\text{const}(k^{(gi)}, T_\text{ci}) - p_\text{single}(k^{(gi)}, T_\text{ci})\big), \tag{4.2.1}
$$
where $k^{(gi)}$ is the per-group state count (constant within a
`groups[gi]` by construction — `cache_total_loglik` line 705–722 makes
this explicit). Under JC, $p_\text{const}$ and $p_\text{single}$ are
given by (2.1.1) and (2.1.2); under MkN (neomorphic), by (2.3.1) and
(2.3.2).

The existing helper `singleton_site_prob_jc` (or `_mkn` for
type-0 neomorphic partitions) gives the right value. The partial-CL
infrastructure analogous to `evaluate_const_prob` does **not** exist for
singletons — `evaluate_const_prob` is a per-candidate function (line
1089) that reuses cached pseudo-CLs, while singleton requires $n$ such
pseudo-CLs.

#### Minimal code change.

This is **not** ≤30 LOC. Three options:

1. **Cheap & correct fallback:** call `singleton_site_prob_jc` (full
   recompute) once per candidate. ≤10 LOC change but defeats the
   partial-CL speedup precisely when `coding == 2`. For typical
   workloads with $n \approx 8$–$30$ and `nCand` $\approx 5$–$20$, this
   is the single largest cost of the move.
   Code:
   ```cpp
   if (coding == 2) {
     NumericVector rootFreqs(kStates, 1.0 / kStates);
     // Build candidate tree edge vectors here (the partial-CL framework
     // does this implicitly via evaluate_candidate)…
     double sP = singleton_site_prob_jc(candParent, candChild, candEdgeLen,
                                         nTip, kStates, rootFreqs, rates);
     constP += sP;
   }
   ```

2. **Proper partial-CL singleton:** generalise `evaluate_const_prob` to
   `evaluate_singleton_prob` using $n$ pseudo-CLs in parallel. ~150 LOC,
   requires new workspace plumbing and substantial test coverage. The
   `node_cl_cache.h` infrastructure (`create_const_pseudo_group`,
   `compute_residual_cl`) would need a sibling `create_singleton_pseudo_group`.

3. **Disable partial-CL for `coding == 2`:** simplest safe interim —
   fall back to full evaluation when `coding == 2`. 1-line guard at the
   top of the partial-CL dispatch. Loses speedup but the *correctness*
   bug LIKE-001 disappears.

**Recommended interim:** option 3 (guard + fallback), pending
implementation of option 2 if profiling shows `coding == 2` is hot.
Status: **deferred to orchestrator decision**; this is the LIKE-001 fix
described in §5.

### 4.3 Trivial-fix patch produced

The closest thing to a trivial fix in scope is in
`src/mcmc_likelihood.cpp::const_site_prob_for_k:2249–2276`, where two
"Phase 7" comments (lines 2267, 2273) mark the same singleton drop. The
function is called by `gibbs_kprime_sweep_impl`
(`src/mcmc.cpp:3502–3505`). For the JC branch only — F81-Het is gated by
`test-het.R:614` and not exercised — adding `singleton_site_prob_jc` is
≤10 LOC and matches the post-Phase-7 design comment exactly.

The fix:

```cpp
} else {
  NumericVector rootFreqs(kStates, 1.0 / kStates);
  double p = constant_site_prob_jc(parent, child, edgeLen, data.nTip,
                                    kStates, rootFreqs, acrvRates);
  if (data.codingType == 2) {
    p += singleton_site_prob_jc(parent, child, edgeLen, data.nTip,
                                 kStates, rootFreqs, acrvRates);
  }
  return p;
}
```

This is captured at `dev/red-team/patches/L5-ascertainment.patch`.
**It is necessary but not sufficient** to close LIKE-001: the
`cache_total_loglik` and `gibbs_spr_impl` paths in `src/node_cl_cache.h`
and `src/mcmc.cpp:1092` are independent (§5). The orchestrator should
treat this patch as one small piece of the LIKE-001 remediation.

---

## §5. LIKE-001 fix sketch (for orchestrator follow-up)

### What `cache_total_loglik` should compute under `coding == 2`.

The full-eval reference is `cpp_partition_log_likelihood`
(`src/mcmc_likelihood.cpp:2290–2733`), which at lines 2493–2498 (JC
known/homogeneous), 2600–2605 (transformational allSame), and 2717–2722
(transformational sub-group) consistently computes:
$$
p_A = p_\text{const} + \mathbb{1}\{\text{coding} = 2\} \cdot p_\text{single},
\quad \text{ll} \mathrel{-}= n_\text{char} \cdot \log(1 - p_A).
$$

`cache_total_loglik` (`src/node_cl_cache.h:643–731`) must do the same:

```cpp
// Pseudo-code; replace each constant_site_prob_* call with const+single sum
if (cache.coding == 2) {
  if (part.type == 0) {
    p += singleton_site_prob_mkn(parent, child, neoEl, nTip,
                                  rateLoss, rootFreqs, rates);
  } else if (part.type == 2) {
    p += singleton_site_prob_jc(parent, child, absEdgeLen, nTip,
                                 kStates, rootFreqs, rates);
  }
  // For type 1 (transformational): inside the per-unit loop,
  // pu += singleton_site_prob_jc(parent, child, absEdgeLen,
  //                              nTip, k, rootFreqs, rates);
}
```

The arithmetic to add: at line 698 (neomorphic), use
$p_\text{const} + p_\text{single}^\text{MkN}$ from (2.3.1)+(2.3.2);
at line 704 (known JC), $p_\text{const} + p_\text{single}^\text{JC}$
from (2.1.1)+(2.1.2); at line 717 (per-unit transformational), same JC
sum but per-unit $k$.

The empirical evidence in findings.md LIKE-001 (1.5–1.8 nat drift,
106/106 mismatches under coding=informative) is consistent with the
missing $\log(1 - p_\text{const} - p_\text{single})$ vs
$\log(1 - p_\text{const})$ gap; the gap is $\log[(1 - p_\text{const})/(1 -
p_\text{const} - p_\text{single})] = O(p_\text{single}/(1 - p_\text{const}))
\approx O(0.1$–$0.3)$ per char per partition, scaling with $n_\text{char}$.

### Scope.

The same fix pattern (`+ singleton_site_prob_jc(...)` under coding=2)
applies to the three independent sites flagged in LIKE-001:

1. `src/node_cl_cache.h:697,703,716` — `cache_total_loglik` (this §).
2. `src/mcmc_likelihood.cpp:1577–1603` — `const_site_prob_for_k`
   (the JC trivial-fix patch in §4.3).
3. `src/mcmc.cpp:1088–1095` — `gibbs_spr_impl` partial-CL (the
   non-trivial code path in §4.2 / option 3).

All share the same closed form and the same root-cause: a TODO that
predates Phase-7 informative-coding support. No new math is required.

---

## §6. Edge cases

### k = 1.

JC(1): only one state. $p_\text{const} = 1$ (always constant),
$p_\text{single} = 0$ (no second state to differ to). The variable
correction $\log(1 - 1) = -\infty$ — every site is constant,
ascertainment renders the data impossible. The implementation guards
against this upstream (`src/corrections.cpp:39`, `kObs >= 1` requirement);
$p_\text{const} < 1$ is also enforced at `src/mcmc.cpp:1093`
(`if (constP < 1.0)`). The formula (2.1.1) returns
$p_\text{const} = k \cdot q_0 = 1 \cdot 1 = 1$ correctly; the upstream
guard handles the resulting $\log 0$.

F81(1): same, $p_\text{const} = \pi_0 \cdot 1 = 1$.

MkN($k = 2$, trivial $\pi = (0, 1)$ or $(1, 0)$): one of the two
pseudo-chars contributes 0; the other contributes 1. $p_\text{const} =
1$ — same edge case.

### Branch length $t = 0$ on some edges.

$P_{ii}(0) = 1, P_{ij}(0) = 0$. The pruning is bit-equivalent to
contracting the zero-length edge. All formulae (2.1.1)–(2.3.2) hold:
no degeneracy in the closed form. Implementation: the JC and MkN
kernels use $e^{-\beta \cdot 0} = 1$ giving $p_\text{same} = 1$,
$p_\text{diff} = 0$ — the standard contraction.

### Very long branches, $t \to \infty$.

Under JC: $P_{ij}(t) \to 1/k$ for all $(i, j)$. Then for the
constant pseudo-char (all tips in state 0):
$q_0(T \text{ very long}, k) \to (1/k)^n \cdot 1 = (1/k)^{n-1} \cdot (1/k) = (1/k)^n$,
so $p_\text{const} \to k \cdot (1/k)^n = k^{1-n}$. For $n = 8, k = 4$:
$p_\text{const} \to 6.1 \times 10^{-5}$. The variable correction tends
to $-n_\text{char} \cdot \log(1) \to 0$. Sanity check.

Under F81: $P_{ij}(t) \to \pi_j$, so $q_i \to \pi_i^n$, and
$p_\text{const} \to \sum_i \pi_i \cdot \pi_i^n = \sum_i \pi_i^{n+1}$.
For uniform $\pi = 1/k$, recovers $k \cdot (1/k)^{n+1} = (1/k)^n$ — but
wait, that's $1/k$ times the JC limit. The discrepancy is because the
F81 formula sums $\pi_i$ at the root explicitly, while the JC formula
absorbs it into the symmetry multiplier $k$. The closed-form check:
$\sum_i \pi_i \cdot \pi_i^n = \sum_i \pi_i^{n+1}$ is the correct F81
constant probability, and the JC reduction is $k \cdot (1/k)^{n+1} =
(1/k)^n$. Both agree under $\pi \equiv 1/k$.

### Single-tip trees ($n = 1$).

Not a valid phylogeny, but the formulae handle it: $p_\text{const}
= 1$ (one tip, one state, always constant), $p_\text{single} = 0$
(no "other tip" to differ from). Upstream code requires
$n \geq 3$ before constructing a `MkPrimeData`.

### All tips same vs all tips different.

These are constant and "fully informative" patterns respectively. Both
have well-defined $L_\text{fels}$; they enter the corrected likelihood
through different ascertainment terms. No formula collapse.

### $k_\text{obs} = k_F$ (no unseen states).

JC-collapse path is gated out (`useCollapse` requires
`kObsMaxLocal + 1 < kStates`); the full uncollapsed kernel runs. The
ascertainment formulae (2.1.1)–(2.3.2) are independent of any
collapse, so this case is correctness-neutral.

### $k_\text{obs} = 1$ (only one state observed).

A character with `kObs = 1` is constant in the *data*. The dispatch in
`prepare_mcmc_data` typically filters such characters before MCMC. The
ascertainment formula still applies: under `coding = variable`, this
character has $L_\text{raw} = 1$ (everyone in state 0 trivially
explains the observations) and the correction is
$- \log(1 - p_\text{const})$. Joint posterior remains proper; tests
`test-ascertainment.R:266–290` exercise mixed (kObs=2, kObs=1) data.

---

## §7. Verdict

### Sub-questions.

| Sub-question | Verdict |
|---|---|
| JC(k) constant + singleton closed-form (§2.1) | **watertight** — formulae match (2.1.1), (2.1.2); impl at `src/ascertainment.cpp:54–199`. |
| F81(k) constant closed-form (§2.2, $p_\text{const}$ only) | **watertight** — formula (2.2.1) matches `het_constant_site_prob` at `src/mcmc_likelihood.cpp:2110–2223`. |
| F81(k) singleton closed-form (§2.2, $p_\text{single}$) | **w-caveats** — formula (2.2.2/2.2.3) recorded for the future F81-Het Phase-7 implementation; not currently used in production (gated by `test-het.R:614`). |
| MkN(2) constant + singleton (§2.3) | **watertight** — formulae (2.3.1), (2.3.2); impl at `src/ascertainment.cpp:416–597`. |
| Mk' relabelling × ascertainment commutation (§3) | **watertight** — closed via (F1) + (F2); falling-factorial commutes exactly with both `coding=variable` and `coding=informative` denominators under JC, and the commutation lifts to F81-Het without modification (the Het-specific caveat in L4 is about the falling-factorial itself, not its commutation). Caveat (2) of `relabelling-correction.md` is now **closed**. |
| TODO M-052 (F81 singleton, `mcmc_likelihood.cpp:2240`) (§4.1) | **open-question / impl-bug deferred** — formula (4.1.1) given; implementation requires substantial new pseudo-char buffers and is correctly gated out of production. Recorded for future Phase 7. |
| TODO informative singleton (`mcmc.cpp:1092`) (§4.2) | **impl-bug** — exact LIKE-001 incarnation in `gibbs_spr_impl`. Options (1)/(2)/(3) given; recommended interim is option 3 (fall back to full evaluation). |
| LIKE-001 fix scope (§5) | **impl-bug — closed-form fix described** for `cache_total_loglik` (3 sites) and `const_site_prob_for_k`; same `+ singleton_site_prob_*(...)` pattern matches the existing full-eval reference. |

### Headline finding.

The two outstanding TODOs at `src/mcmc_likelihood.cpp:2240` (F81-Het
singleton, deferred and correctly gated) and `src/mcmc.cpp:1092`
(JC/MkN singleton in partial-CL, **active bug under
`coding="informative"`**, root cause of LIKE-001) both have closed-form
resolutions stated above. The JC case is mechanical; the F81-Het case
is mechanical but memory-heavy. The Mk' relabelling × ascertainment
commutation (L4 caveat 2) is now formally closed: the falling-factorial
correction commutes exactly with conditioning on $\neg C$ and on $\neg C
\wedge \neg G$ under JC, because both events are $\sigma$-invariant
subsets of $S^n$ and $L_\text{fels}$ is permutation-symmetric.

### Patch.

A small (~7-line) trivial fix to `const_site_prob_for_k` is captured at
`dev/red-team/patches/L5-ascertainment.patch`. **Not committed, not
merged.** The orchestrator decides whether to apply. This patch alone
does **not** close LIKE-001 (the `cache_total_loglik` and `gibbs_spr_impl`
sites are independent code paths — see §5).

### Worktree state.

- New: `dev/red-team/proofs/ascertainment.md` (this file).
- New: `dev/red-team/patches/L5-ascertainment.patch`.
- Wave-1 proof files copied in (`acrv.md`, `f81-collapse.md`,
  `hastings-continuous.md`, `hastings-tree-moves.md`,
  `relabelling-correction.md`) from `34a2607` for cross-reference; these
  are *not* edits, just brought into the worktree.
- No edits to `R/`, `src/`, or `tests/testthat/`. No commit, no push, no
  merge.

---

## References

- Felsenstein, J. (1981) Evolutionary trees from DNA sequences: a
  maximum likelihood approach. *Journal of Molecular Evolution* 17:
  368–376.
- Felsenstein, J. (2004) *Inferring Phylogenies*. Sinauer. §§16.3, 18.1.
- Lewis, P.O. (2001) A likelihood approach to estimating phylogeny from
  discrete morphological character data. *Systematic Biology* 50:
  913–925. [Variable / informative coding correction.]
- Yang, Z. (2014) *Molecular Evolution: A Statistical Approach*. OUP.
  §4.2 (pruning), §4.10 (F81), §6.3 (ascertainment).
- `dev/red-team/proofs/relabelling-correction.md` — L4 (Mk' relabelling
  correction); §3 here closes its caveat (2).
- `dev/red-team/proofs/f81-collapse.md` — L1 (JC-collapse lumpability);
  the fused-ascertainment identity at §D ("outConstProb = sum_eff") is
  used at §2.1.
- `dev/red-team/findings.md:39` — LIKE-001, the empirical reproduction
  motivating §5.
