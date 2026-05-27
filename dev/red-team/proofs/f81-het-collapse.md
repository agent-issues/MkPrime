# F81-Het state-collapse design — proof

**Lane:** L? (math-prover, follow-up to L1 / `f81-collapse.md`).
**Scope:** Determine whether the JC(k) lumpability collapse landed in
PR #2 (`b99b0ae` / `8d4638e` / `a966d8b`) generalises to the F81-Het
kernel `pruning_f81_het_acrv_flat` (`src/mcmc_likelihood.cpp:1742–1983`)
and its constant-site companion `het_constant_site_prob`
(`src/mcmc_likelihood.cpp:2138–2253`).
**Author:** math-prover agent, 2026-05-27.
**Brief:** "Close (or definitively defer) the F81-Het collapse design"
(orchestrator wave-4).

> **Reference artefacts.**
> - JC proof `dev/red-team/proofs/f81-collapse.md` (especially §6).
> - Design stub on branch `wip/f81-het-collapse-design`, file
>   `dev/plans/2026-05-21-f81-het-collapse-design.md` (commit
>   `0efd7e7`).
> - Mk' relabelling proof
>   `dev/red-team/proofs/relabelling-correction.md`.

---

## Verdict (one-line, up front)

**HOLDS-PARTIAL — the *likelihood* collapse holds, by a pruning-algebraic
argument (not strong-lumpability); the *fused-ascertainment* collapse
holds only partially and needs the design choice spelled out below.**
Expected per-call speedup `≈ kFull / kEff` on the likelihood pass (4× at
kFull=12, kObs=2) plus ascertainment savings of a similar ratio on the
spike-on-O rotations only. Total kernel speedup ≈ 3×–4× in the
likelihood-dominated regime, contingent on profile evidence that the
F81-Het kernel is a bottleneck.

---

## Theorem (informal)

Under the F81-Het mixture used in MkPrime, the per-component pi-vector
takes the form
$\pi^{(\mathrm{rot})}_s = \beta\,\mathbb{1}\{s = \mathrm{rot}\} + r\,\mathbb{1}\{s \neq \mathrm{rot}\}$
with $r = (1-\beta)/(k_F-1)$. Fix the partition
$\Pi = \{\{0\},\ldots,\{k_O-1\},U\}$ where $U = \{k_O, \ldots, k_F-1\}$,
$n_U = k_F - k_O \geq 2$.

Then for every fixed (ACRV-cat, Beta-bin, rotation) component:

1. If $\mathrm{rot} \in O = \{0,\ldots,k_O-1\}$: pi is constant on $U$
   (every $s \in U$ has $\pi_s = r$). The chain is **strongly lumpable**
   on $\Pi$ (Kemeny–Snell), the pruning identity proven in
   `f81-collapse.md` Lemmas 1–5 carries over with the new lumped
   transition kernel below, and the lumped site likelihood equals the
   full-chain site likelihood exactly.
2. If $\mathrm{rot} \in U$: pi is *not* constant on $U$ (one element has
   $\pi = \beta$, the rest have $\pi = r$). Strong lumpability **fails**.
   However, by label-symmetry of unobserved states, all $n_U$ rotations
   in $U$ produce identical site likelihoods, so they can be replaced
   by a *single* rotation × $n_U$ weight. The single representative
   computation still requires a $k_F$-wide chain (no collapse), but it
   runs $n_U$ times less often.

Combining (1) and (2), the per-(cat, bin) cost drops from
$k_F$ rotations × $k_F$-wide pruning to $k_O$ collapsed passes (width
$k_{\mathrm{eff}} = k_O + 1$) plus $1$ uncollapsed pass (width $k_F$,
weighted ×$n_U$).

For fused ascertainment, claim (1) gives a collapsed
spike-on-O contribution but the constant-state pseudo-characters
indexed by $s \in U$ break the class-constancy hypothesis (H);
ascertainment for those pseudo-chars must stay uncollapsed (or use a
finer partition — see §6).

---

## Assumptions

1. **Substitution model.** The F81-Het per-(bin, rot) Q-matrix is
   $Q_{ij} = \mu (\pi_j - \mathbb{1}\{i = j\})$ with
   $\mu = 1/(1 - \sum_j \pi_j^2)$. Transition probabilities are
   $P_{ij}(t) = \pi_j (1 - d) + \delta_{ij} d$ with $d = e^{-\mu t}$.
   Encoded at `src/mcmc_likelihood.cpp:1877` ($\mu$),
   `src/mcmc_likelihood.cpp:1892–1894` ($d$, $1-d$).
2. **Stationary frequencies.** Pi is taken to be the per-component
   stationary distribution; rotations are over which state receives the
   "high" $\beta$ frequency. Encoded at
   `src/mcmc_likelihood.cpp:1862–1875` (build pi for one bin × rot).
   For $k \geq 3$: $\pi_{\mathrm{rot}} = \beta$,
   $\pi_s = r := (1-\beta)/(k_F-1)$ for $s \neq \mathrm{rot}$.
   For $k = 2$: $n_{\mathrm{Rot}} = 1$ (no rotation); the model
   degenerates and the collapse case is moot (no $U$ states with
   $n_U \geq 2$ unless $k_F \geq k_O + 2$, which is impossible at
   $k_F = 2$).
3. **Tree, edges, tip data, missing data.** Same as the JC proof
   (`f81-collapse.md` Assumptions 1, 2, 4): tip states remapped to
   $\{0,\ldots,k_O-1\}$, missing $= -1$, $t \geq 0$.
4. **Root distribution.** The pruning kernel uses the per-component
   stationary $\pi$ at the root (`src/mcmc_likelihood.cpp:1949–1952`).
   This is the F81 stationary, *not* the JC uniform.
5. **Mixture averaging.** Site likelihood per character is
   $L_c = \frac{1}{R} \sum_{\mathrm{cat},\,\mathrm{bin},\,\mathrm{rot}}
   L_c^{(\mathrm{cat},\,\mathrm{bin},\,\mathrm{rot})}$
   with $R = n_{\mathrm{Cat}} \cdot n_{\mathrm{BetaCat}} \cdot
   n_{\mathrm{Rot}}$. Encoded at
   `src/mcmc_likelihood.cpp:1971–1977` (`avg = siteLikPtr[c] *
   inv_comp`; `inv_comp = 1/totalComp`).
6. **Rotation count.** $n_{\mathrm{Rot}} = k_F$ for $k_F \geq 3$; the
   collapse design is for this regime. (`src/mcmc_likelihood.cpp:1774`.)
7. **Dispatch gate.** Collapse only ever fires when $k_F > k_O + 1$,
   equivalently $n_U \geq 2$.
8. **Profile contingency.** The whole design is **contingent on
   profile evidence** that `pruning_f81_het_acrv_flat` is a bottleneck
   in Het-enabled runs (design stub §C). The proof below is correct
   irrespective of profile evidence, but implementation should not
   proceed without it.

---

## Statement (formal)

Let $T$, $\rho$, $L$, $\{t_e\}$, $x : L \to O \cup \{-1\}$ be as in
the JC proof. Fix a single (cat, bin) pair; let $\beta \in (0, 1)$ be
the bin value, $r = (1-\beta)/(k_F - 1)$.

For each $\mathrm{rot} \in S = O \cup U$, write $\pi^{(\mathrm{rot})}$
for the spike-on-rot pi vector. The per-component site likelihood is

$$
L^{(\mathrm{rot})}(x) \;=\; \sum_{s \in S} \pi^{(\mathrm{rot})}_s \cdot
a^{(\mathrm{rot})}_\rho(s),
$$

with $a^{(\mathrm{rot})}_v$ the F81 conditional likelihood under
$\pi^{(\mathrm{rot})}$.

Define the **partition $\Pi$** and the *average-over-class* lumped CL

$$
\hat a^{(\mathrm{rot})}_v(I) := \frac{1}{|I|} \sum_{s \in I}
a^{(\mathrm{rot})}_v(s), \quad I \in \Pi.
$$

For $\mathrm{rot} \in O$, define the lumped pi
$\hat\pi^{(\mathrm{rot})}(I) = \sum_{s \in I} \pi^{(\mathrm{rot})}_s$.
Concretely: $\hat\pi^{(\mathrm{rot})}(\{j\}) = \beta$ if $j = \mathrm{rot}$,
$\hat\pi^{(\mathrm{rot})}(\{j\}) = r$ for the other $j \in O$, and
$\hat\pi^{(\mathrm{rot})}(U) = n_U \cdot r$.

**Claim A (spike-on-O lumped likelihood).** For every
$\mathrm{rot} \in O$, the F81 site likelihood equals the lumped-chain
site likelihood:

$$
L^{(\mathrm{rot})}(x) = \sum_{I \in \Pi} \hat\pi^{(\mathrm{rot})}(I)\,
\hat a^{(\mathrm{rot})}_\rho(I), \tag{A}
$$

with $\hat a^{(\mathrm{rot})}_\rho$ computable by post-order recursion
on a $k_{\mathrm{eff}}$-wide buffer using the lumped transition kernel
(7) derived below.

**Claim B (spike-on-U exchange).** For every $\mathrm{rot}, \mathrm{rot}'
\in U$,

$$
L^{(\mathrm{rot})}(x) = L^{(\mathrm{rot}')}(x). \tag{B}
$$

Consequently, $\sum_{\mathrm{rot} \in U} L^{(\mathrm{rot})}(x) =
n_U \cdot L^{(u^\star)}(x)$ for any fixed $u^\star \in U$ — i.e.
the $n_U$ U-rotations can be replaced by *one* uncollapsed pass × $n_U$.

**Claim C (overall savings).** Per (cat, bin), the total work is:

- $k_O$ collapsed pruning passes of width $k_{\mathrm{eff}} = k_O + 1$
  (claim A);
- 1 uncollapsed pruning pass of width $k_F$ (claim B representative).

Total flops per (cat, bin): $(k_O \cdot k_{\mathrm{eff}}^2 + k_F^2) /
(k_F \cdot k_F^2) = (k_O (k_O + 1)^2 + k_F^2) / k_F^3$ times the
uncollapsed cost. For $k_F = 12$, $k_O = 2$: $(2 \cdot 9 + 144) / 1728
= 162/1728 \approx 0.094$, i.e. ~10.6× per-(cat,bin) speedup on the
likelihood pass. Ascertainment savings are smaller — see §6.

---

## Proof

### Lemma F1 (spike-on-O lumpability: Kemeny–Snell).

Fix $\mathrm{rot} = r^\star \in O$. Then $\pi^{(\mathrm{rot})}_s$ is
constant ($= r$) for every $s \in U$.

For F81, $Q_{ss'} = \mu (\pi_{s'} - \mathbb{1}\{s = s'\})$. Compute the
row sums into class $J \in \Pi$:

$$
\sum_{s' \in J} Q_{s s'} = \mu\,\pi(J) - \mu \pi_s\,\mathbb{1}\{s \in J\}.
$$

For $s \in I$, we need this to be independent of the choice of $s \in I$:

- **$I$ is a singleton $\{j\} \subseteq O$.** Trivially constant in $s$.
- **$I = U$**, $J \neq U$. Then $s \in U$, $s \notin J$, so the indicator
  is 0 and the sum is just $\mu \pi(J)$ — independent of $s$. ✓
- **$I = U$**, $J = U$. The sum is $\mu \pi(U) - \mu \pi_s$, but
  $\pi_s = r$ is constant on $U$ (by hypothesis $r^\star \in O$).
  Independent of $s \in U$. ✓

By Kemeny & Snell (1960) Theorem 6.3.2, the F81 chain is **strongly
lumpable** on $\Pi$ when $\mathrm{rot} \in O$. The lumped generator
is

$$
\hat Q^{(\mathrm{rot})}_{IJ} \;=\; \begin{cases}
\mu\,\hat\pi(J), & J \neq I, \\
\mu\,\hat\pi(I) - \mu\,\pi_s, & J = I,\text{ any }s \in I.
\end{cases} \tag{6}
$$

Note: when $I = \{j\}$ is a singleton, $\hat\pi(I) = \pi_j$ and
$\hat Q^{(\mathrm{rot})}_{II} = -\mu(1 - \pi_j)$. When $I = U$,
$\pi_s = r$ for any $s \in U$, so
$\hat Q^{(\mathrm{rot})}_{UU} = \mu(n_U r - r) = \mu (n_U - 1) r$,
which is a valid CTMC generator on $\Pi$ (rows sum to zero by direct
check). The lumped transition kernel is

$$
\hat P^{(\mathrm{rot})}(t) = \exp(\hat Q^{(\mathrm{rot})} t).
$$

Because every off-diagonal $\hat Q^{(\mathrm{rot})}_{IJ}$ has the F81
form $\mu \hat\pi(J)$ (only the diagonal carries the $-\mu \pi_s$
correction), the lumped chain on $\Pi$ is itself an F81 chain on $k_O
+ 1$ states with frequencies $(\pi_0, \ldots, \pi_{k_O - 1},
\hat\pi(U)) = (r, \ldots, r, \beta, r, \ldots, r, n_U r)$ (with $\beta$
at position $r^\star$ if $r^\star \in O$).

Therefore the lumped transition probabilities admit the same closed
form as the full chain, with **the same $\mu$ and the same $d = e^{-\mu
t}$** (this is critical — $\mu$ is computed from
$\sum_s \pi_s^2 = \beta^2 + (k_F - 1)r^2$ on the full chain, *not*
from the lumped pi — see §6 caveat 1).

$$
\hat P^{(\mathrm{rot})}_{IJ}(t) = \hat\pi(J)\,(1 - d) +
\delta_{IJ}\,d. \tag{7}
$$

### Lemma F2 (spike-on-O pruning identity).

The proof is structurally identical to JC Lemmas 2–5
(`f81-collapse.md:194–289`). The class-constancy hypothesis
(H): $a_v(s)$ constant on each $I \in \Pi$ propagates from tips
(Lemma 4 of the JC proof; uses Assumption 3 = data-remap) up the tree
under any chain that is strongly lumpable on $\Pi$. F81 with
$\mathrm{rot} \in O$ qualifies by Lemma F1.

The only line that differs from JC is the root marginalisation:

$$
\sum_s \pi^{(\mathrm{rot})}_s\,a^{(\mathrm{rot})}_\rho(s)
= \sum_{I \in \Pi} \hat\pi^{(\mathrm{rot})}(I)\,
\hat a^{(\mathrm{rot})}_\rho(I).
$$

By (H): for $I = \{j\}$, $\hat a(I) = a(j)$ and the contribution is
$\pi_j a(j)$. For $I = U$, $\hat a(U) = (1/n_U)\sum_{s \in U} a(s)$ and
$\hat\pi(U) = n_U r$, so the contribution is $r \cdot \sum_{s \in U}
a(s) = \sum_{s \in U} \pi_s a(s)$ (since $\pi_s = r$ on $U$). $\square$

### Lemma F3 (spike-on-U exchange via label permutation).

Fix two rotations $u, u' \in U$. Let $\sigma : S \to S$ be the
transposition swapping $u$ and $u'$ (and fixing everything else). $\sigma$
acts on the F81 chain and on pi:

- $\pi^{(u')}_s = \pi^{(u)}_{\sigma(s)}$ — i.e. $\sigma$ pulls back pi.
- The F81 generator is symmetry-equivariant under any permutation of
  state labels: $Q^{(u')}_{s,s'} = Q^{(u)}_{\sigma(s),\sigma(s')}$, since
  $Q_{ss'} = \mu(\pi_{s'} - \mathbb{1}\{s = s'\})$ and $\sigma$ commutes
  with the indicator.
- The Felsenstein recursion is similarly equivariant: if all tip CLs
  are permuted by $\sigma$, every internal CL is permuted by $\sigma$.

Since $u, u' \in U$ and tip data lies in $O$, **$\sigma$ leaves every
tip CL invariant** (tip CL is $\mathbb{1}\{s = j\}$ for some $j \in O$,
which $\sigma$ — acting on $U$ — does not touch; or all-ones for
missing). Equivalently, $a^{(u')}_\ell(s) = a^{(u)}_\ell(s)$ for all
tips $\ell$ and all $s$.

By induction up the tree (equivariance), $a^{(u')}_v(s) =
a^{(u)}_v(\sigma(s))$. Wait — there is a subtlety. We have

$$
a^{(u')}_v(s) \;=\; \prod_c \sum_{s'} P^{(u')}_{ss'}(t_{vc})
a^{(u')}_c(s').
$$

By equivariance, $P^{(u')}_{ss'} = P^{(u)}_{\sigma(s)\sigma(s')}$ and by
the inductive hypothesis $a^{(u')}_c(s') = a^{(u)}_c(\sigma(s'))$ would
need to hold. The base case (tips): we need
$a^{(u')}_\ell(s') = a^{(u)}_\ell(\sigma(s'))$. Since tip CL is
$\mathbb{1}\{s' = j\}$ with $j \in O$ and $\sigma$ fixes $O$ pointwise,
$\mathbb{1}\{s' = j\} = \mathbb{1}\{\sigma(s') = j\}$ iff $s'$ and
$\sigma(s')$ are either both $j$ or both not $j$. Since $j \in O$ and
$\sigma$ only permutes elements of $U$:

- If $s' \in O$: $\sigma(s') = s'$, so equality holds trivially.
- If $s' \in U$: $\sigma(s') \in U$ (still in $U$, just possibly
  different). Both $s'$ and $\sigma(s')$ are in $U$, hence neither
  equals $j \in O$. So $\mathbb{1}\{s' = j\} = 0 =
  \mathbb{1}\{\sigma(s') = j\}$. ✓

Similarly for missing-data tips ($a_\ell \equiv 1$, trivially
invariant). So the base case holds: $a^{(u')}_\ell(s) = a^{(u)}_\ell(\sigma(s))$
for all tips and all $s$.

By induction (the recursion preserves equivariance under any
permutation):

$$
a^{(u')}_v(s) = a^{(u)}_v(\sigma(s)). \tag{8}
$$

At the root:

$$
L^{(u')}(x) = \sum_s \pi^{(u')}_s a^{(u')}_\rho(s)
= \sum_s \pi^{(u)}_{\sigma(s)} a^{(u)}_\rho(\sigma(s))
= \sum_{s} \pi^{(u)}_s a^{(u)}_\rho(s) = L^{(u)}(x),
$$

using $\sigma$-bijectivity for the reindexing. $\square$

### Combining F2 and F3.

For each (cat, bin), the sum over rotations decomposes as

$$
\sum_{\mathrm{rot} \in S} L^{(\mathrm{rot})}(x)
= \sum_{r^\star \in O} L^{(r^\star)}(x) + n_U \cdot L^{(u^\star)}(x)
$$

for any fixed $u^\star \in U$. The first sum has $k_O$ terms, each
computable by **collapsed** pruning on $\Pi$ at width $k_{\mathrm{eff}} =
k_O + 1$ (Lemma F2). The second is a single uncollapsed pass at width
$k_F$, with weight $n_U$ in the final mixture average.

Dividing by $n_{\mathrm{Cat}} \cdot n_{\mathrm{BetaCat}} \cdot
n_{\mathrm{Rot}} = n_{\mathrm{Cat}} \cdot n_{\mathrm{BetaCat}} \cdot k_F$
recovers the F81-Het mixture site likelihood
(`src/mcmc_likelihood.cpp:1971–1977`). $\blacksquare$

---

## §6. Fused-ascertainment caveat (the partial part)

The fused-asc machinery
(`src/mcmc_likelihood.cpp:1818–1848`, `1923–1942`, `1955–1964`)
propagates $k_F$ pseudo-characters in parallel, the $s$-th of which
has tip CL $e_s$ (all tips constant in state $s$). The output is
$P(\text{constant site}) = \sum_s \pi_s \cdot P(\text{all tips}=s)$,
the same weighted sum across the $k_F$ pseudo-chars.

**For pseudo-chars indexed by $s \in O$:** tip CL = $\mathbb{1}\{s' = s\}$,
$s \in O$. Class-constancy on $\Pi$ holds (JC proof Lemma 4: this is
exactly an observed-state tip). The proof of (A) applies: under
spike-on-$O$ rotation, the asc pseudo-char-$s$ pruning collapses on
$\Pi$.

**For pseudo-chars indexed by $s \in U$:** tip CL = $\mathbb{1}\{s' = s\}$
with $s \in U$. On the class $U$, $a_\ell(s') = 1$ if $s' = s$ and 0
otherwise — **not class-constant on $U$**. Hypothesis (H) fails at the
base case. The pruning identity (A) does **not** hold for these
pseudo-chars.

**Resolution options.**

1. **Mixed-collapse design.** Compute the $k_O$ "main CL +
   $s \in O$ asc CLs" together on the collapsed $\Pi$ buffer (spike-on-O
   rotations only); compute the $n_U$ "$s \in U$ asc CLs" on the
   uncollapsed buffer (still under spike-on-O rotations — these are not
   class-constant on $U$). By Lemma F3 (label symmetry), all $n_U$
   $s \in U$ pseudo-chars give the **same** $P(\text{all tips}=s)$
   value under any single spike-on-$O$ rotation, *and* under the single
   $u^\star$-rotation representative for spike-on-U; so just **one**
   uncollapsed asc pseudo-char per rotation, multiplied by $n_U$,
   suffices. This gives roughly half the asc savings (only spike-on-O
   half collapses).

2. **Two-state partition for asc.** Use a finer partition
   $\Pi' = \{\{0\}, \ldots, \{k_O-1\}, \{s^\star\}, U \setminus
   \{s^\star\}\}$ of size $k_O + 2$ for ascertainment only; spike-on-O
   rotations have pi constant on both $\{s^\star\}$ (trivially) and
   $U \setminus \{s^\star\}$ (still $= r$), so strong lumpability
   holds on $\Pi'$. The $s^\star$-asc pseudo-char is class-constant on
   $\Pi'$. For asc indexed by $s \in U \setminus \{s^\star\}$, use
   label symmetry to argue they all equal the $s^\star$ value. Width
   $k_O + 2$ vs $k_F$ still saves substantially for $n_U$ large.

**The proof above does not pick between (1) and (2).** This is the
"design choice" that must be made before kernel writing. Either is
correct; (2) is finer and cheaper but adds plumbing.

### §6.1 Singleton ascertainment.

`het_singleton_site_prob` at `src/mcmc_likelihood.cpp:2256–2276`
currently returns 0 unconditionally — informative coding is not yet
supported under Het. The collapse design therefore does not need to
treat singleton-asc, but **flag**: when M-052 (TODO at
`src/mcmc_likelihood.cpp:2270`) lands, the same mixed-collapse
machinery should apply.

> Caveat (6.1): exact F81 singleton-site probability is a TODO. Closing
> this requires a separate proof (analogous to the constant-site fused
> mechanism). Should be picked up by lane L? / role math-prover once
> M-052 is implemented.

---

## §7. Lumped transition-kernel formulae (kernel skeleton)

For implementation, the per-rotation closed-form lumped P(t) for
$\mathrm{rot} = r^\star \in O$ is:

Let pi-vector on the lumped chain be $\tilde \pi = (\pi_0, \ldots,
\pi_{k_O-1}, n_U r)$ at positions $0, \ldots, k_O - 1, k_O$ (the last
index = the lumped column). Note $\sum_I \tilde\pi(I) = (k_O - 1) r +
\beta + n_U r = \beta + (k_F - 1) r = \beta + (1 - \beta) = 1$. ✓

The mu used in the lumped P(t) is the **same mu as the full chain** —
$\mu = 1/(1 - \sum_{s \in S} \pi_s^2) = 1/(1 - \beta^2 - (k_F - 1)
r^2)$. **Critical:** do not recompute mu from $\tilde\pi$ — that
gives the wrong rate. (See §8 edge case.)

The lumped P(t) entries (from (7)):

- Diagonal singleton $I = \{j\}$, $j \in O$:
  $\hat P_{II}(t) = \pi_j (1 - d) + d$.
- Off-diagonal singleton $I = \{i\} \to J = \{j\}$, both in $O$:
  $\hat P_{IJ}(t) = \pi_j (1 - d)$.
- Singleton $I = \{j\}$ to lumped $J = U$:
  $\hat P_{IU}(t) = n_U r (1 - d) = \tilde\pi_U (1 - d)$.
- Lumped $I = U$ to singleton $J = \{j\}$, $j \in O$:
  $\hat P_{UJ}(t) = \pi_j (1 - d)$. (no diagonal contribution.)
- Lumped $I = U$ to lumped $J = U$:
  $\hat P_{UU}(t) = n_U r (1 - d) + d = \tilde\pi_U (1 - d) + d$.

These satisfy row-stochasticity:
$\sum_J \hat P_{IJ}(t) = (1 - d) \sum_J \tilde\pi_J + d = (1 - d) + d = 1$
for any starting class $I$ (regardless of whether $I$ is singleton or
lumped). ✓

In particular, **the lumped P has the F81 form** with pi replaced by
$\tilde\pi$ — so the same O(k) Felsenstein hoisting trick used at
`src/mcmc_likelihood.cpp:1899–1921` (precompute `sum_pi_cl =
sum_j pi[j] * clCh[j]`; broadcast `base = (1-d) * sum_pi_cl`; add
`d * clCh[i]` per column) **carries over verbatim** with $\tilde\pi$
replacing $\pi$ on a $k_{\mathrm{eff}}$-wide buffer. The recursion
needs a per-column "size weight" only at the root marginalisation
(matching JC-collapsed practice — `f81-collapse.md` §C).

### Pseudocode (NOT C++).

```
function pruning_f81_het_acrv_flat_collapsed(...):
    // Standard prologue: nEdge, nTip, nChar, nCat, kFull, kObs,
    // partition states, etc.
    kEff = kObs + 1
    nU = kFull - kObs

    // Init lumped tip CLs (data-remap invariant: tips in [0, kObs-1])
    for tip in tips:
        for c in chars:
            if state[tip,c] == -1:    // missing
                for i in 0..kEff-1: cl[tip, c, i] = 1.0
            else:
                cl[tip, c, state[tip,c]] = 1.0  // observed state
                // lumped col cl[tip, c, kObs] = 0.0 by default

    // (Optional) Asc tips: split into O-pseudo-chars (kObs of them) +
    // one representative U-pseudo-char (will be multiplied by nU).
    // O-pseudo-chars use the kEff-wide buffer.
    // U-pseudo-char uses a separate kFull-wide buffer (under spike-on-O
    // rotations) -- see §6 option 1. Alternatively, the §6 option 2
    // kEff+1-wide partition for asc.

    site_lik = zeros(nChar)
    const_prob_sum = 0

    for cat in 1..nCat:
        rateA = rate_multipliers[cat]
        for bi in 1..nBetaCat:
            beta = betaBins[bi]
            r = (1 - beta) / (kFull - 1)
            sum_pi_sq_full = beta**2 + (kFull-1)*r**2
            mu = 1 / (1 - sum_pi_sq_full)

            // ===== spike-on-O rotations: kObs collapsed passes =====
            for rstar in 0..kObs-1:
                // lumped pi: pi_tilde[j] = r for j in O\{rstar},
                //            pi_tilde[rstar] = beta,
                //            pi_tilde[kObs] = nU * r        // lumped col
                build pi_tilde of length kEff

                // (Identical inner loop to current F81 kernel,
                // operating on kEff-wide CL buffer)
                for n in internal_nodes: init_flag[n] = 0
                for e in edges (post-order):
                    t = edge_length[e] * rateA
                    one_minus_d = -expm1(-mu * t)
                    d = 1 - one_minus_d

                    for c in chars:
                        sum_pi_cl = sum_{i=0..kEff-1} pi_tilde[i] * clCh[c,i]
                        base = one_minus_d * sum_pi_cl
                        for i in 0..kEff-1:
                            clPar[c,i] = base + d * clCh[c,i]
                            // (or *= for non-first child)

                // root: site_lik[c] += sum_i pi_tilde[i] * clRoot[c,i]
                for c in chars:
                    site_lik[c] += sum_i pi_tilde[i] * clRoot[c,i]

                // O-asc pseudo-chars: kObs of them, propagate on the
                // same kEff buffer; root contributes pi_tilde[s_asc]
                // weight times sum_i pi_tilde[i] * ascRoot[s,i]
                // (analogous to JC fused-asc, §6 option 1).
                ... (asc accumulation; const_prob_sum += ...)

                // U-asc pseudo-char representative: ONE uncollapsed
                // kFull-wide pruning, multiplied by nU. Use existing
                // uncollapsed asc inner loop.
                ... (uncollapsed asc; const_prob_sum += nU * ...)

            // ===== spike-on-U: ONE uncollapsed representative × nU =====
            // Use existing uncollapsed F81 kernel for ONE rotation
            // (e.g., rot = kObs), then multiply contribution by nU.
            uncollapsed_pruning(rot=kObs, pi=spike-on-kObs, mu=mu, ...)
            for c in chars:
                site_lik[c] += nU * <root sum for this rotation>
            // U-asc under U-rotation: by Lemma F3 again, identical
            // for any pseudo-char index; one representative × nU × nU.
            // Or equivalently kFull-1 of the asc pseudo-chars collapse
            // to one expression. (Spell out in implementation.)
            ... (asc under spike-on-U; const_prob_sum += ...)

    // Average over (nCat × nBetaCat × kFull) components -- DENOMINATOR
    // is still kFull, not kEff, because we've expanded all rotations.
    inv_comp = 1 / (nCat * nBetaCat * kFull)
    for c in chars: logLik += log(site_lik[c] * inv_comp)
    if outConstProb: *outConstProb = const_prob_sum * inv_comp

    return logLik
```

**Key invariants enforced by the pseudocode.**

- `mu` is computed from the FULL pi-vector ($\beta^2 + (k_F-1) r^2$).
- The denominator in the mixture average is still $n_{\mathrm{Cat}}
  \cdot n_{\mathrm{BetaCat}} \cdot k_F$ — not $k_O + 1$. The
  spike-on-U pass is weighted ×$n_U$ to compensate for not running it
  $n_U$ times.
- The collapsed buffer is $k_{\mathrm{eff}}$-wide. Existing flat-buffer
  stride must accommodate; gating $k_F > k_O + 1$ ensures
  $k_{\mathrm{eff}} < k_F$, so no buffer-size growth.

---

## §8. Edge cases

1. **$k_O = k_F - 1$ ($n_U = 1$).** Then $\Pi$ has $k_F$ singletons —
   the trivial partition. Lemma F2 still holds, but the collapse gives
   no width savings ($k_{\mathrm{eff}} = k_F$). Spike-on-U becomes a
   single rotation (no factor savings either). Dispatch should gate
   this out with `kObsMax + 1 < kFull` (strict), as the JC dispatch
   does. **No bug, no savings, gate out.**

2. **$k_O = 1$ ($n_U = k_F - 1$).** Maximal collapse. Partition $\Pi
   = \{\{0\}, U\}$, $k_{\mathrm{eff}} = 2$. The single spike-on-O
   rotation runs at width 2; one spike-on-U representative runs at
   width $k_F$. Asymptotically the spike-on-U cost dominates ($k_F^2$
   vs $1 \cdot 4$); but it runs once per (cat, bin) instead of $n_U =
   k_F - 1$ times, so the overall savings remain ~$(k_F - 1)$×.
   Lemmas F1–F3 hold; no edge concerns.

3. **$k_F = 2$.** $n_{\mathrm{Rot}} = 1$
   (`src/mcmc_likelihood.cpp:1774`). Either $k_O = 1$ (gates collapse
   on with $k_{\mathrm{eff}} = 2 = k_F$ — trivial, no savings) or
   $k_O = 2$ (gate doesn't fire). The collapse design is moot for
   $k_F = 2$. The $k_F = 2$ pi-construction at
   `src/mcmc_likelihood.cpp:1862–1868` is asymmetric (uses
   `gain_base`/`loss_base` from `baseRL`) and doesn't fit the
   spike-on-rot template — but again, this is in the no-savings
   regime so no design impact.

4. **Branch length $t = 0$.** $d = 1$, $1 - d = 0$. $\hat P(0) = I$
   trivially. The pruning identity holds.

5. **Long branches $t \to \infty$.** $d \to 0$. $\hat P_{IJ} \to
   \tilde\pi(J)$ — class-independent. Lumped CL at any internal node
   becomes $\hat a_v(I) = \prod_c \tilde\pi^\top \hat a_c$ — constant
   in $I$. Equilibrium behaviour. ✓

6. **Beta bin at $\beta = 1/k_F$ (uniform pi).** Then $r = (1 -
   1/k_F)/(k_F - 1) = 1/k_F = \beta$. Pi-vector is uniform $1/k_F$;
   all rotations give the same pi. The chain is JC. Lumpability holds
   (Lemma F1 trivially: pi constant everywhere) and Lemma F3 is
   trivial. The collapsed F81 design reduces to the JC collapse already
   proven. ✓

7. **Beta bin at extreme value ($\beta \to 0$ or $\beta \to 1$).**
   `compute_het_bins` guards degenerate bins at
   `src/mcmc_likelihood.cpp:1705–1707`. As long as $\beta \in (0, 1)$
   strictly (which the Beta-quantile binning enforces — boundaries are
   never exact 0 or 1 in IEEE), mu is finite and the proof applies.
   At $\beta = 1$ exactly, $r = 0$ and the chain becomes absorbing —
   pi[rot] = 1, all others 0, and mu = 1/(1 - 1) → infinity. Not
   reachable in practice.

8. **All tips missing.** $a_\ell \equiv 1$ on $S$; class-constant on
   every class of every partition; site lik = 1. Same as JC §7 case 6.

9. **Numerical: `mu * t` very small.** Already guarded by `expm1` form
   (`src/mcmc_likelihood.cpp:1893`) — FAST-EXP-001 patch landed in
   `59276ad`. Collapsed kernel should re-use the same `expm1` form;
   noted in pseudocode.

10. **Spike-on-O rotation index $r^\star \notin O$ never happens by
    construction.** The dispatch always iterates `rot in 0..kFull-1`;
    the "spike-on-O vs spike-on-U" split is purely in the new kernel's
    control flow. Adding the split is a logical operation, not a data
    dependency.

---

## §9. Implementation cross-check (for the future kernel)

This is a forward-looking cross-check, since the kernel does not yet
exist. The dispatch sites that would change:

- **Gibbs sweep persite** at `src/mcmc.cpp:3743–3791` (current JC
  dispatch). Currently the `useCollapse` flag is only set when
  `!useHet`. The new dispatch would also set `useCollapse` when
  `useHet && kObsMax + 1 < kFull`, and dispatch to a new
  `pruning_f81_het_acrv_persite_collapsed`.
- **Partition-LL homogeneous** at `src/mcmc_likelihood.cpp:2458–2480`.
  Currently `useHet` branch (line 2461 onwards: `compute_het_bins(...)`
  then `pruning_f81_het_acrv_flat(...)`) is unconditional in the Het
  arm. Add a `useCollapse = useHet && (kObsMax + 1 < kFull)` gate;
  when set, dispatch the new collapsed Het kernel.
- **Partition-LL transformational allSame** at
  `src/mcmc_likelihood.cpp:2562–2602`. Same pattern.
- **Partition-LL transformational sub-group** at
  `src/mcmc_likelihood.cpp:2676–2724`. Same pattern.

**Workspace sizing.** The fused-asc U-pseudo-char path needs a
$k_F \times k_F$-strided buffer for the uncollapsed inner loop, in
addition to the $k_{\mathrm{eff}}$-wide main CL buffer. Net memory ≤
existing F81-Het kernel (which already allocates $k_F \times k_F$ asc
stride at `src/mcmc_likelihood.cpp:1822`).

**Singleton-site coding == 2.** Out of scope (see §6.1); kernel still
uses uncollapsed `singleton_site_prob_jc` / Het variant as today.

---

## §10. Profile contingency

Per design stub §C and project memory `project_f81_collapse_design`,
this design is **not** to be implemented without first profiling a
Het-enabled MCMC run and showing that `pruning_f81_het_acrv_flat`
accounts for more than ~30% of wall time. The proof above is valid
irrespective of the profile, but the development effort (estimated
~600 lines of new C++ mirroring the JC collapse, plus tests) is only
justified by the profile evidence.

> Caveat (10): profile evidence pending. Should be picked up by lane
> L? / role numerical-auditor (or whoever has run Het profiles) before
> kernel implementation.

---

## Verdict

**HOLDS-PARTIAL.**

The F81-Het collapse is mathematically sound by a *hybrid* argument:

- **Likelihood (main CL):** A pruning-algebraic + label-symmetry hybrid
  splits the $k_F$ rotations into $k_O$ collapsible spike-on-O passes
  (width $k_{\mathrm{eff}}$, by strong lumpability — Lemma F1) plus a
  single uncollapsed spike-on-U representative weighted by $n_U$ (by
  label symmetry — Lemma F3). Per-(cat, bin) flops drop from $k_F^3$
  to $k_O (k_O + 1)^2 + k_F^2$ — roughly $(k_F / k_{\mathrm{eff}})$×
  speedup at large $k_F$, small $k_O$.

- **Fused ascertainment (k pseudo-chars):** Only the $k_O$
  $s \in O$-pseudo-chars are class-constant on $\Pi$; the $n_U$
  $s \in U$-pseudo-chars are not, so they must use the uncollapsed
  inner loop. By label symmetry within $U$, all $n_U$
  $U$-pseudo-chars give the same value per rotation, so only **one**
  uncollapsed asc representative per rotation × $n_U$ weight is
  needed. Ascertainment savings are smaller — only the $k_O$ O-asc
  passes collapse, the U-asc reps stay uncollapsed.

- **Singleton ascertainment:** Out of scope (currently returns 0;
  M-052 TODO).

- **Profile gate:** Implementation contingent on evidence that the
  Het kernel is a bottleneck (design stub §C).

**No patch produced.** This is a design proof; the kernel is ~600
lines of new C++ + tests, outside the trivial-fix policy. The proof
delivers:

1. Mathematical green light for the collapse design (Lemmas F1–F3).
2. Closed-form lumped transition kernel (§7 formulae).
3. Pseudocode skeleton with all invariants identified (§7 pseudocode).
4. Design decision flagged but not made: §6 option 1 (mixed-collapse
   asc) vs §6 option 2 (finer partition for asc). Choose before
   implementation.
5. Edge-case enumeration (§8) with no correctness blockers.
6. Profile contingency (§10) — do not implement before profiling.

**Out-of-lane follow-ups (for the next session/wave):**

> Caveat (V1): the §6 design choice between mixed-collapse asc (cheap,
> half-savings) and finer-partition asc (richer, fuller savings) needs
> to be made by whoever implements the kernel. Should be picked up by
> lane L? / role math-prover + numerical-auditor jointly (numerical
> for relative cost estimate; math for proof closure of option 2).

> Caveat (V2): profile evidence that
> `pruning_f81_het_acrv_flat` is a bottleneck. Should be picked up by
> lane L? / role numerical-auditor.

> Caveat (V3): once M-052 (exact F81 singleton-site probability)
> lands, the singleton-asc collapse needs its own proof
> (Lemma F3-analogue for the singleton pseudo-character family).
> Should be picked up by lane L? / role math-prover.

Status of this branch: read-only consultation of
`wip/f81-het-collapse-design`; no commits proposed.
