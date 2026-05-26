# JC(k) state-collapse via lumpability — proof

**Lane:** L1 (math-prover)
**Scope:** PR #2 commits `b99b0ae` (persite Gibbs sweep), `8d4638e` (flat-CL
partition path with fused ascertainment), `a966d8b` (dead-code removal).
**Author:** math-prover agent, 2026-05-26.

> **Scope clarification.** The lane title in the orchestrator brief says
> "F81-Het collapse". The three commits cited in the brief are explicitly
> **JC(k)** (i.e. F81 with uniform stationary frequencies) and explicitly
> *defer* the F81-Het case ("Het is deferred — JC lumpability holds only
> under equal stationary frequencies within the lumped class" —
> `b99b0ae` commit body, mirrored in `8d4638e`). Project memory
> `project_f81_collapse_design` confirms: "JC lumpability collapse landed in
> PR #2 commits b99b0ae/8d4638e/a966d8b (~3.7× kernel speedup); F81-Het
> version deferred — design stub on branch `wip/f81-het-collapse-design`."
>
> This proof therefore covers exactly the JC(k) collapse that landed.
> §5 ("JC limit") collapses to a tautology and §6 records what would change
> for an F81-Het extension.

---

## Theorem (informal title)

Let $J = (J_t)_{t \geq 0}$ be a continuous-time Markov chain on
$S = \{0, 1, \ldots, k_F - 1\}$ with the JC(k_F) generator
$Q_{ij} = \beta/(k_F - 1)$ for $i \neq j$, $Q_{ii} = -\beta$, where
$\beta = k_F/(k_F-1)$ is the unit-time-substitution normalisation already
in use. Fix an observed state count $k_O$ with $1 \leq k_O < k_F$ and
assume the tip data uses only states in $O = \{0, \ldots, k_O - 1\}$. Let
$U = \{k_O, \ldots, k_F-1\}$ be the unseen states, $|U| = n_U = k_F - k_O$.
Consider the lumped chain on the partition
$\Pi = \{\{0\}, \{1\}, \ldots, \{k_O-1\}, U\}$
of size $k_{\mathrm{eff}} = k_O + 1$.

Then for any rooted tree with non-negative branch lengths and any character
whose tip states lie in $O \cup \{-1\}$ (missing), the Felsenstein
site-likelihood computed on the lumped chain — with the explicit
"average-over-class" CL convention used in
`src/likelihood.cpp::pruning_jc_collapsed:152–161` — equals the site
likelihood of the full chain exactly (no normalising constant, $C = 1$).

The same identity holds with ACRV gamma rate categories applied as a
per-category rate multiplier in front of the entire $Q$.

---

## Assumptions

1. **Substitution model.** The per-edge transition matrix is
   $P(t) = \exp(Q t)$ with $Q$ the JC(k_F) generator above; equivalently
   $P_{ii}(t) = 1/k_F + (1 - 1/k_F)\,e^{-k_F t/(k_F-1)}$ and
   $P_{ij}(t) = 1/k_F - (1/k_F)\,e^{-k_F t/(k_F-1)}$ for $i \neq j$.
   This formula is used at `src/likelihood.cpp:213–215`,
   `src/mcmc_likelihood.cpp:984–986`, `src/mcmc_likelihood.cpp:1150–1152`,
   `src/mcmc_likelihood.cpp:704–706`, etc.
2. **Branch lengths.** $t \geq 0$ for every edge; finite. The kernels do
   not guard $t < 0$ but the surrounding code clamps. $t = 0$ gives the
   identity transition, a degenerate but valid case (§7).
3. **Stationary / root distribution.** $\pi = (1/k_F, \ldots, 1/k_F)$,
   the uniform JC stationary. The kernels hard-code this via
   `NumericVector rootFreqs(kStates, 1.0/kStates)` at the dispatch sites
   (`src/mcmc_likelihood.cpp:2457`, `2564`, `2680`) and via `inv_k = 1/kFull`
   in the collapsed roots
   (`src/mcmc_likelihood.cpp:1037`, `1201`, `src/likelihood.cpp:248`).
4. **Tip data invariant.** Tip states are remapped to be contiguous in
   $\{0, \ldots, k_O - 1\}$ before reaching pruning, where $k_O$ is the
   number of distinct *non-missing* states observed for the character /
   batch in question. Missing data is encoded as $-1$. See
   `src/likelihood.cpp:200` ("state is guaranteed to be in 0..kObs-1 by the
   data-remap invariant") and the dispatch in
   `src/mcmc_likelihood.cpp:2459–2461`, `2566–2567`, `2682–2684`, where
   $k_O^{\max}$ is the maximum `part.kObsLocal[ci]` over the batch /
   sub-group.
5. **Batch dispatch.** When several characters share the kernel call, the
   dispatch uses $k_O^{\max}$ = `max(part.kObsLocal[…])`. Characters with
   smaller individual $k_O$ ride along, with their unused upper observed
   columns simply never populated. This is a sub-partition of $\Pi$ — see
   §6.
6. **Het is out of scope.** F81-Het has non-uniform stationary frequencies
   *within* what would be the lumped class, breaking strong lumpability of
   the JC partition. The dispatch never sends Het to the collapsed kernel
   (the `useHet` branches at `src/mcmc_likelihood.cpp:2429–2447`, `2536–2560`,
   `2651–2675` keep the uncollapsed F81-Het kernel).
7. **Singleton ascertainment.** Coding $== 2$ adds a *separate*
   $P(\text{singleton})$ correction via the uncollapsed
   `singleton_site_prob_jc` helper (`src/mcmc_likelihood.cpp:2495`, `2602`,
   `2724`). Singleton probability depends only on $k_F$ by full-chain JC
   symmetry, so this is left on the uncollapsed path by design.

---

## Statement (formal)

Let $T = (V, E)$ be a rooted tree with root $\rho$, leaf set $L$, edge
lengths $\{t_e\}_{e \in E}$. For a single character with tip assignment
$x : L \to O \cup \{-1\}$, Felsenstein's pruning algorithm gives

$$
\mathcal{L}(x) \;=\; \sum_{s \in S} \pi_s \cdot a_\rho(s),
$$

where $a_\rho : S \to \mathbb{R}_{\geq 0}$ is the conditional likelihood
(CL) vector at the root, computed by the post-order recursion

$$
a_v(s) \;=\; \prod_{c \in \mathrm{ch}(v)} \Big(\sum_{s' \in S} P_{s s'}(t_{vc})\, a_c(s')\Big), \qquad
a_\ell(s) = \mathbb{1}\{x_\ell = s\} \text{ or } 1 \text{ (missing)},\; \ell \in L.
$$

Define the lumped CL on $\Pi$ by the *average-over-class* convention

$$
\hat a_v(I) \;:=\; \frac{1}{|I|} \sum_{s \in I} a_v(s), \qquad I \in \Pi. \tag{1}
$$

Define the lumped transition kernel by

$$
\hat P_{I J}(t) \;=\; \sum_{s' \in J} P_{s s'}(t), \quad s \in I \text{ arbitrary,}
\qquad I, J \in \Pi, \tag{2}
$$

provided the right-hand side does not depend on the choice of $s \in I$
(strong lumpability — Lemma 1 below verifies this for JC).

Define $\hat\pi(I) = |I|/k_F$ (uniform stationary on $S$, summed over class).

**Claim.** Under Assumptions 1–7,

$$
\mathcal{L}(x) \;=\; \sum_{I \in \Pi} \hat\pi(I)\, \hat a_\rho(I) \tag{3}
$$

where $\hat a_\rho$ is computed by the post-order recursion on the lumped
chain

$$
\hat a_v(I) \;=\; \prod_{c \in \mathrm{ch}(v)} \Big(\sum_{J \in \Pi} \hat P_{IJ}(t_{vc})\, \hat a_c(J)\Big), \qquad
\hat a_\ell \text{ as in (1) applied to the full-chain tip CL.}
\tag{4}
$$

That is, **(3) = (4) for every character**, exactly (not up to a constant).

---

## Proof

### Lemma 1 (JC is strongly lumpable under any partition).

For JC(k_F), $Q_{ij} = \beta/(k_F-1)$ for $i \neq j$. For any partition
$\Pi$ of $S$, and any $I \neq J$ in $\Pi$,

$$
\sum_{s' \in J} Q_{s s'} = \frac{\beta\,|J|}{k_F - 1}, \qquad s \in I,
$$

which does not depend on $s$. For $I = J$,
$\sum_{s' \in I} Q_{s s'} = -\beta + \beta(|I|-1)/(k_F-1)$ — also
independent of $s \in I$. By Kemeny & Snell (1960), *Finite Markov Chains*,
Theorem 6.3.2 (strong-lumpability criterion: the row-sums $\sum_{s' \in J}
Q_{ss'}$ are constant in $s \in I$ for every $J \in \Pi$), $J$ restricted
to states is strongly lumpable on $\Pi$, and the lumped chain on $\Pi$ is
itself Markovian with generator

$$
\hat Q_{IJ} = \sum_{s' \in J} Q_{s s'} \;=\; \frac{\beta\,|J|}{k_F - 1} \quad (I \neq J),
\qquad
\hat Q_{II} = -\beta + \frac{\beta\,(|I|-1)}{k_F-1}. \tag{5}
$$

Exponentiating, $\hat P(t) = \exp(\hat Q t)$ satisfies

$$
\hat P_{IJ}(t) = \sum_{s' \in J} P_{s s'}(t) = |J| \cdot P_{ss'}(t), \quad s \in I, s' \in J, \tag{6}
$$

since all off-diagonal entries of $P(t)$ are equal to $P_{ij}(t) = p_{\mathrm{diff}}(t)$
and all diagonals equal $p_{\mathrm{same}}(t) = 1/k_F + (1 - 1/k_F)e^{-k_F t/(k_F-1)}$.
Concretely:

$$
\hat P_{IJ}(t) = \begin{cases}
p_{\mathrm{same}}(t) + (|I|-1)\,p_{\mathrm{diff}}(t), & J = I,\\
|J|\,p_{\mathrm{diff}}(t), & J \neq I.
\end{cases} \tag{7}
$$

Independence of representative $s \in I$ is the strong-lumpability
condition. $\square$

### Lemma 2 (Pruning identity at an internal node).

Let $v \in V$ with child $c$, edge length $t$. Suppose
$\hat a_c(I) = (1/|I|)\sum_{s' \in I} a_c(s')$ for all $I \in \Pi$. Then,
for $s \in I$,

$$
\sum_{s' \in S} P_{s s'}(t)\, a_c(s') \;=\; \sum_{J \in \Pi} \hat P_{IJ}(t)\, \hat a_c(J). \tag{8}
$$

*Proof.* The left side decomposes as $\sum_{J \in \Pi} \sum_{s' \in J}
P_{ss'}(t) a_c(s')$. Since $P_{ss'}(t)$ depends on $s' \in J$ only through
whether $s' = s$ — and we will pick $s \in I$, so the only "diagonal" case
is $J = I$:

- For $J \neq I$, $s' \in J$ means $s' \neq s$, so $P_{ss'}(t) =
  p_{\mathrm{diff}}(t)$, constant in $s'$. Thus
  $\sum_{s' \in J} P_{ss'}(t) a_c(s') = p_{\mathrm{diff}}(t) \cdot |J|
  \cdot \hat a_c(J) = \hat P_{IJ}(t)\, \hat a_c(J)$ by (7).
- For $J = I$, split the sum at $s' = s$:
  $\sum_{s' \in I} P_{ss'}(t) a_c(s') = p_{\mathrm{same}}(t) a_c(s) +
  p_{\mathrm{diff}}(t) \sum_{s' \in I,\, s' \neq s} a_c(s')$
  $= p_{\mathrm{same}}(t) a_c(s) + p_{\mathrm{diff}}(t) (|I| \hat a_c(I) - a_c(s))$.

But the **right side** of (8), summed over the diagonal contribution
$J = I$, gives $\hat P_{II}(t) \hat a_c(I) =
(p_{\mathrm{same}}(t) + (|I|-1) p_{\mathrm{diff}}(t)) \hat a_c(I)$.

These are not equal *for arbitrary $a_c$*. They are equal **iff**
$a_c(s) = \hat a_c(I)$ for all $s \in I$ — i.e. the full-chain CL is
constant within each class $I \in \Pi$.

So (8) is not unconditional; it requires the inductive hypothesis

$$
a_c(s) = \hat a_c(I) \text{ for all } s \in I \in \Pi. \tag{H}
$$

We now show (H) holds at every node *for the data in scope*. $\square$

### Lemma 3 (Class-constancy is preserved by pruning, given tip-class-constancy).

Suppose (H) holds at every child of $v$. Then (H) holds at $v$.

*Proof.* By Felsenstein recursion, for $s \in I$,

$$
a_v(s) = \prod_{c \in \mathrm{ch}(v)} \Big( \sum_{s' \in S} P_{s s'}(t_{vc}) a_c(s') \Big).
$$

By the Lemma-2 computation applied to each child, the bracketed term for
child $c$ equals (under (H) at $c$)

$$
B_c(s) := p_{\mathrm{same}}(t) a_c(s) + p_{\mathrm{diff}}(t) \big(\sum_{s' \neq s} a_c(s')\big)
$$

$= (p_{\mathrm{same}}(t) - p_{\mathrm{diff}}(t)) a_c(s) + p_{\mathrm{diff}}(t)
\sum_{s'} a_c(s')$.

The second summand is independent of $s$. The first summand, restricted
to $s \in I$, gives $(p_{\mathrm{same}} - p_{\mathrm{diff}}) \hat a_c(I)$
under (H) at $c$, also independent of choice of $s \in I$. So $B_c(s)$
is constant in $s \in I$. Products of class-constant functions are
class-constant, so $a_v(s)$ is constant in $s \in I$, i.e. (H) at $v$. $\square$

### Lemma 4 (Tips satisfy (H)).

For a tip $\ell$ with observation $x_\ell \in O \cup \{-1\}$ and the
partition $\Pi = \{\{0\}, \ldots, \{k_O-1\}, U\}$:

- If $x_\ell$ is missing, $a_\ell(s) = 1$ for all $s \in S$. Trivially
  class-constant on every $I \in \Pi$.
- If $x_\ell = j \in O$, $a_\ell(s) = \mathbb{1}\{s = j\}$.
  - On singleton class $\{i\}$ for $i \neq j$: $a_\ell(s) = 0$ for the one
    $s \in \{i\}$ — class-constant.
  - On singleton class $\{j\}$: $a_\ell(s) = 1$ — class-constant.
  - On lumped class $U$: $a_\ell(s) = 0$ for all $s \in U$ (since $j \notin
    U$) — class-constant. $\square$

This is exactly where the data-remap invariant (Assumption 4) is used: if
the data contained a state in $U$, $a_\ell$ would not be class-constant
on $U$ and the lumping would fail.

### Lemma 5 (Root identity).

Under (H) at $\rho$,

$$
\sum_{s \in S} \pi_s a_\rho(s) = \sum_{I \in \Pi} \hat\pi(I) \hat a_\rho(I).
$$

*Proof.* $\sum_{s \in S} \pi_s a_\rho(s) = \sum_{I \in \Pi} \sum_{s \in I}
(1/k_F) a_\rho(s) = \sum_I (1/k_F) |I| \hat a_\rho(I) = \sum_I (|I|/k_F)
\hat a_\rho(I) = \sum_I \hat\pi(I) \hat a_\rho(I)$. $\square$

### Putting it together.

By Lemma 4, tips satisfy (H). By Lemma 3 (induction up the tree), every
internal node satisfies (H). By Lemma 2 + (H), the lumped recursion (4)
produces $\hat a_v(I) = (1/|I|) \sum_{s \in I} a_v(s)$ at every node.
By Lemma 5 at the root, the lumped marginalisation (3) reproduces the
full-chain site likelihood exactly. **Constant $C = 1$.**

ACRV with $n_C$ rate categories applies a per-category multiplier $r_c$
to all branch lengths, so the entire derivation runs unchanged for each
category, and the site likelihood is the average over categories. The
ACRV kernel
(`src/mcmc_likelihood.cpp:1140–1211` in `pruning_jc_acrv_flat_collapsed`,
`src/mcmc_likelihood.cpp:686–795` in `persite_collapsed_impl`) does
exactly this — outer loop over `cat`, inner pruning identical to the
non-ACRV form, average via `siteLikPtr[c] * (1/nCat)`.

$\blacksquare$

---

## Implementation cross-check

The reference implementation that the proof models is
`src/likelihood.cpp::pruning_jc_collapsed:167–256`. The three hot-path
kernels added by the cited commits replicate this math with progressively
more performance machinery; we check each line by line.

### A. Tip initialisation (Lemma 4)

| Step                              | `pruning_jc_collapsed` (ref)     | `persite_collapsed_impl` (b99b0ae) | `pruning_jc_flat_collapsed` (8d4638e) | `pruning_jc_acrv_flat_collapsed` (8d4638e) |
|-----------------------------------|----------------------------------|------------------------------------|---------------------------------------|--------------------------------------------|
| Observed state $s \in O$: $CL[s] = 1$, rest 0 | `src/likelihood.cpp:201`         | `src/mcmc_likelihood.cpp:680`      | `src/mcmc_likelihood.cpp:947`         | `src/mcmc_likelihood.cpp:1110`             |
| Missing: $CL[s] = 1$ on all $k_{\mathrm{eff}}$ entries | `src/likelihood.cpp:197–198`    | `src/mcmc_likelihood.cpp:677–678`  | `src/mcmc_likelihood.cpp:943–944`     | `src/mcmc_likelihood.cpp:1107–1108`        |

All four use the *average-over-class* CL convention (lumped col $= 1$ for
missing, not $n_U$), matching (1) and Lemma 4 ("$a_\ell$ class-constant").

### B. Edge update (Lemma 2 + 3)

The implementation uses the standard JC O(k) trick
$\mathrm{new\_CL}[i] = p_{\mathrm{diff}} \sigma + (p_{\mathrm{same}} -
p_{\mathrm{diff}}) \mathrm{CL}[i]$, where $\sigma$ is the row-sum
$\sum_{s'} P_{ss'} \cdot CL[s']$ reduced via JC symmetry. For the lumped
chain the relevant sum is

$$
\sigma_{\mathrm{eff}} := \sum_{s' \in S} P_{s s'}(t) \cdot a_c(s')
\cdot \frac{1}{|J(s')|} \cdot |J(s')|
= \sum_{J} \hat P_{IJ}(t) \hat a_c(J),
$$

but in the **lumped recursion (4)** we just need the row sum *in the
lumped CL representation*. From (7),

$$
\sigma_{\mathrm{eff}} = p_{\mathrm{diff}}(t) \sum_{J} |J| \hat a_c(J) +
(p_{\mathrm{same}} - p_{\mathrm{diff}})(t) \cdot \hat a_c(I).
$$

The factor $\sum_J |J| \hat a_c(J) = \sum_{j < k_O} 1 \cdot \hat a_c(\{j\})
+ n_U \cdot \hat a_c(U)$ is computed at:

- `src/likelihood.cpp:221–223` (CL[ch][offset..]: `for (j<kObs) sum += …`
  then `sum += n_U * CL[ch][offset+kObs]`)
- `src/mcmc_likelihood.cpp:758–760` (`persite_collapsed_impl`)
- `src/mcmc_likelihood.cpp:994–996` (`pruning_jc_flat_collapsed`)
- `src/mcmc_likelihood.cpp:1160–1162` (`pruning_jc_acrv_flat_collapsed`)

Then `clPar[i] = p_diff * sum_eff + diff_coeff * clCh[i]` at:

- `src/likelihood.cpp:225`
- `src/mcmc_likelihood.cpp:761`, `773`
- `src/mcmc_likelihood.cpp:997–998`, `1007–1008`
- `src/mcmc_likelihood.cpp:1163–1164`, `1173–1174`

These match (4) under (H): for $i \in I$,
$\hat a_v(I) = \sigma_{\mathrm{eff}}$, which equals
$p_{\mathrm{diff}} \cdot \sum_J |J| \hat a_c(J) + (p_{\mathrm{same}} -
p_{\mathrm{diff}}) \cdot \hat a_c(I)$, and the code writes the same value
into every column $i$ of class $I$ — preserving (H) post-update. The
"singleton classes have only one column" representation collapses the
"all $i \in I$" loop into one column for $I \in \{\{0\}, \ldots,
\{k_O-1\}\}$ and one for $I = U$.

**Note on the tip-edge fast path** (persite_collapsed_impl only,
`src/mcmc_likelihood.cpp:716–751`): when the child is a tip, `sum_eff`
collapses to $1$ (observed) or $k_F$ (missing). The comment at
`src/mcmc_likelihood.cpp:721–724` derives:

- Missing child, $\Sigma_{\mathrm{eff}} = k_O + n_U = k_F$, so
  $\mathrm{new}[i] = p_{\mathrm{diff}} k_F + (p_{\mathrm{same}} -
  p_{\mathrm{diff}}) \cdot 1$. By definition of the JC formula,
  $p_{\mathrm{diff}} k_F = 1 - e^{-k_F t/(k_F-1)}$ and
  $p_{\mathrm{same}} - p_{\mathrm{diff}} = e^{-k_F t/(k_F-1)}$. Sum = 1.
  Code: `for i: clPar[i] = 1.0` at `src/mcmc_likelihood.cpp:725`.
- Observed child at state $s$, $\Sigma_{\mathrm{eff}} = 1$, so
  $\mathrm{new}[i] = p_{\mathrm{diff}} + (p_{\mathrm{same}} -
  p_{\mathrm{diff}}) \mathbb{1}\{i = s\} = p_{\mathrm{diff}}$ for $i \neq
  s$, $p_{\mathrm{same}}$ for $i = s$. Code:
  `src/mcmc_likelihood.cpp:730–731`. $\checkmark$

The lumped col here is always one of the "$i \neq s$" cases (since
observed $s$ is never the lumped class), so it gets $p_{\mathrm{diff}}$.
This is correct: lumped-class CL for an observed tip is
$\hat a_\ell(U) = 0$ pre-edge; post-edge, $\hat a_v(U) =
\sum_J \hat P_{UJ}(t) \hat a_\ell(J)$ = $|J(s)| \cdot p_{\mathrm{diff}}(t)
\cdot \hat a_\ell(\{s\})$ = $p_{\mathrm{diff}}(t) \cdot 1$ since $|\{s\}|
= 1$. $\checkmark$

### C. Root marginalisation (Lemma 5)

`pruning_jc_collapsed` (ref): `sum_eff = Σ_{s<kObs} CL[root][s] + n_U
CL[root][kObs]`; `site_lik = (1/kFull) * sum_eff`. Source:
`src/likelihood.cpp:243–252`.

In all three new kernels: identical lines

- `src/mcmc_likelihood.cpp:780–783` (persite_collapsed_impl per-cat)
- `src/mcmc_likelihood.cpp:1032–1037` (pruning_jc_flat_collapsed)
- `src/mcmc_likelihood.cpp:1196–1202` (pruning_jc_acrv_flat_collapsed
  per-cat)

This matches (3): $\sum_I \hat\pi(I) \hat a_\rho(I) = (1/k_F)\sum_{j < k_O}
\hat a_\rho(\{j\}) + (n_U/k_F) \hat a_\rho(U) = (1/k_F)\,
\Sigma_{\mathrm{eff}}$. $\checkmark$

### D. Fused ascertainment (commit 8d4638e only)

The fused asc machinery propagates a *second* CL — the "pseudo-character
all tips in state 0" — alongside the main CL on the same tree, using the
same edge updates. At the root, summing $\hat\pi(I) \hat a^{\mathrm{asc}}_\rho(I)$
gives $P(\text{all tips in state 0}) =: q_0$. By full-chain JC symmetry
across the $k_F$ states, $P(\text{constant site}) = k_F \cdot q_0$.

In the **uncollapsed** kernel `pruning_jc_flat`,
`src/mcmc_likelihood.cpp:268–272` writes `*outConstProb = sum_i ascRoot[i]`
with the comment "`sum_i cl[i] = k × P(const-in-state-0) = P(const)`" —
i.e. the kernel exploits that the uncollapsed root marginalisation with
uniform $\pi = 1/k$ summed times $k$ exactly cancels: $k \cdot \sum_i
(1/k) \cdot \mathrm{ascRoot}[i] = \sum_i \mathrm{ascRoot}[i]$.

In the **collapsed** kernel, the analogous root sum is
`sum_eff = Σ_{s<kObs} ascRoot[s] + n_U * ascRoot[lumpIdx]`. By the same
proof as (3),

$$
q_0 = \sum_I \hat\pi(I) \hat a^{\mathrm{asc}}_\rho(I) = (1/k_F) \cdot \Sigma_{\mathrm{eff}}.
$$

So $P(\text{const}) = k_F \cdot q_0 = \Sigma_{\mathrm{eff}}$. The code
writes `*outConstProb = sum_eff` at `src/mcmc_likelihood.cpp:1047–1048`
(non-ACRV) and accumulates per-cat sums divided by `nCat` at
`src/mcmc_likelihood.cpp:1209` and `1221–1222` (ACRV) — exactly matching
the uncollapsed convention (`outConstProb` is the average over rate
categories, directly usable in $\log(1 - P(\text{const}))$). $\checkmark$

The pseudo-asc tips are initialised with state 0 in column 0, lumped col
$= 0$: `src/mcmc_likelihood.cpp:972` and `1135` set `ascBufPtr[tip * kEff]
= 1.0`. This matches the "all tips in state 0" pseudo-character; the
ascertainment correction then proceeds with identical edge / root logic
to the main CL, and Lemma 4 / Lemmas 1–5 apply unchanged.

### E. Dispatch criteria

- `b99b0ae` Gibbs sweep dispatch at `src/mcmc.cpp:3743–3791` calls the
  collapsed persite kernel when `ko >= 2` (i.e. $k_F > k_O + 1$, so
  $k_{\mathrm{eff}} < k_F$). $\checkmark$
- `8d4638e` partition-LL dispatch sites at
  `src/mcmc_likelihood.cpp:2458–2473` (homogeneous JC),
  `2562–2580` (transformational allSame),
  `2676–2697` (transformational sub-group):
  all use `kObsMax{Local,Sub} + 1 < k{Full,p0,p}` gate.
  All take `kObsMax` over the batch/sub-group, which is correct: the
  partition $\Pi$ has class $U = \{k_O^{\max}, \ldots, k_F - 1\}$; chars
  with smaller individual $k_O^{(c)} < k_O^{\max}$ simply have empty
  upper observed columns (state never observed → tip CL has 0 there →
  carries through pruning), which is a finer partition than $\Pi$ — and
  Lemma 1 says JC is strongly lumpable under *any* partition, so the
  coarser $\Pi$-pruning still gives the correct full-chain site
  likelihood.

---

## §5. JC limit

**Claim.** When the full-chain stationary frequencies are uniform on all
$k_F$ states (the JC limit), the collapsed pruning reduces to the
pre-PR-#2 JC pruning, in the precise sense that the output log-likelihood
is identical.

This is the entire content of the proof above: the collapsed kernel
is *always* operating in the JC limit by construction (Assumption 3,
hard-coded $\pi_s = 1/k_F$). There is no "limit" to take — the F81
generalisation is the would-be Het extension, which is deferred and not
in PR #2.

For empirical verification, `tests/testthat/test-jc-collapse-hotpath.R`
runs a 384-case sweep covering $k_O \in \{1, \ldots, 4\}$, $k_F$ up to
14, both ACRV settings, missing data, and tree sizes 5/12/25 tips,
showing max LL difference $2.27 \times 10^{-13}$ and max constprob
difference $1.39 \times 10^{-16}$ — consistent with the proof's
prediction of zero up to floating-point reordering.

---

## §6. F81-Het extension (out of scope, recorded for completeness)

For F81-Het, the stationary frequencies $\pi_s$ are *not* uniform. The
strong-lumpability criterion (Lemma 1) requires the row sums
$\sum_{s' \in J} Q_{ss'}$ to be constant in $s \in I$. F81-Het has
$Q_{ss'} = \mu (\pi_{s'} - \mathbb{1}\{s = s'\})$ for some normalising
$\mu$, so

$$
\sum_{s' \in J} Q_{s s'} = \mu\, \pi(J) - \mu \pi_s \mathbb{1}\{s \in J\}.
$$

If $J = U$ and $\pi$ is *not* constant on $U$ — which is the typical
F81-Het case — the term $-\mu \pi_s \mathbb{1}\{s \in J\}$ varies in
$s \in U$, so $\sum_{s' \in U} Q_{s s'}$ varies in $s \in U$, and strong
lumpability fails on the partition $\Pi$.

Strong lumpability of F81-Het on $\Pi$ is recoverable only if $\pi_s$
is constant on $U$ (equivalently, $\pi$ depends only on the partition
class). The design stub on
`wip/f81-het-collapse-design` should constrain $\pi$ to take a single
shared value on the unseen-states class. This is a *modelling*
constraint, not a free numerical optimisation.

The deferral note in the `b99b0ae` and `8d4638e` commit bodies is
therefore exactly right: F81-Het collapse needs additional structure
that JC gets for free.

---

## §7. Edge cases

1. **Equilibrium frequencies with zeros (in the *full* chain).** Not
   possible in JC: $\pi = (1/k_F, \ldots, 1/k_F)$ has no zeros. The
   collapsed `\hat\pi(I) = |I|/k_F` is also strictly positive for every
   non-empty class. ✓

2. **Trivial collapse, $k_O = k_F - 1$, $n_U = 1$.** Partition $\Pi$ has
   $k_F$ singletons; the "lumped" class $U$ has $|U| = 1$. The collapsed
   kernel runs at $k_{\mathrm{eff}} = k_F$, providing zero savings. The
   dispatch gates this out with `kObsMax + 1 < kFull` (strict inequality)
   at the three partition sites and `ko >= 2` (equivalent to
   $k_{\mathrm{eff}} \leq k_F - 1$) in the persite Gibbs sweep — so the
   trivial case is **never dispatched**. If it were (e.g. via the
   `test_flat_jc_collapsed` harness), it would still be correct: the
   $n_U = 1$ collapsed kernel is bit-equivalent to the uncollapsed
   kernel up to a rename of the last column. The kernel comment at
   `src/mcmc_likelihood.cpp:804–807` notes this.

3. **Maximal collapse, $k_O = 1$, $n_U = k_F - 1$.** Partition $\Pi =
   \{\{0\}, \{1, \ldots, k_F-1\}\}$, $k_{\mathrm{eff}} = 2$. All tips
   are in state 0 (the only observed state) or missing. The kernel runs
   at width 2 regardless of $k_F$. Lemmas 1–5 hold; verified by the
   sweep test's $k_O = 1$ corner case (`test-jc-collapse-hotpath.R`).
   Note: this case is degenerate w.r.t. the *parsimony-aware* prior
   downstream, but the *likelihood* is correctly computed.

4. **Branch length $t = 0$.** $P(0) = I$, $p_{\mathrm{same}}(0) = 1$,
   $p_{\mathrm{diff}}(0) = 0$, $\hat P_{IJ}(0) = \mathbb{1}\{I = J\}$. The
   pruning identity holds trivially. The code uses
   $\exp(-k_F \cdot 0 / (k_F - 1)) = 1$, giving the correct degenerate
   transition.

5. **Very long branches, $t \to \infty$.** $p_{\mathrm{same}}(t) \to
   1/k_F$, $p_{\mathrm{diff}}(t) \to 1/k_F$, $\mathrm{diff\_coeff} \to 0$.
   Then $\mathrm{new}[i] = p_{\mathrm{diff}} \Sigma_{\mathrm{eff}}$ — all
   $k_{\mathrm{eff}}$ columns receive the same value, which is the
   correct asymptotic equilibrium ($a_v(s) \propto 1$). The lumped
   recursion gives the same: $\hat P_{IJ}(\infty) = |J|/k_F = \hat\pi(J)$,
   so $\hat a_v(I) = \sum_J \hat\pi(J) \hat a_c(J)$ — class-independent,
   as expected.

6. **All tips missing for a character.** Every $a_\ell$ is identically 1
   on $S$ → (H) holds trivially → root sum $= \prod_e 1 = 1$ → site
   likelihood $= \sum_s \pi_s \cdot 1 = 1$. The collapsed kernel produces
   the same: tips are uniform 1 across $k_{\mathrm{eff}}$ entries,
   internal nodes propagate as 1's (by the "missing tip = identity-on-1"
   computation above), root sum $= \Sigma_{\mathrm{eff}} = k_O + n_U =
   k_F$, site lik $= k_F \cdot (1/k_F) = 1$. $\checkmark$

7. **$k_O = 0$ (no non-missing tips).** Excluded by Assumption 4 and by
   the dispatch precondition `kObsMaxLocal + 1 < kFull` (with
   $k_O^{\max} = 0$ this is `0 + 1 < k_F`, which *does* fire when
   $k_F \geq 2$). The reference kernel guards this at
   `src/likelihood.cpp:173–175` with
   `if (kObs < 1) Rcpp::stop(...)`; the three new kernels do **not**
   re-check. *Is this a bug?* In practice the upstream
   `part.kObsLocal[ci]` is `>= 1` for any character that participates in
   pruning (a character with no non-missing tips is degenerate and
   filtered earlier). But: if upstream ever passed $k_O^{\max} = 0$, the
   collapsed kernel would run at $k_{\mathrm{eff}} = 1$ (single lumped
   column $U = \{0, \ldots, k_F-1\}$, $n_U = k_F$), which still satisfies
   Lemmas 1–5 (the trivial partition $\Pi = \{S\}$ is always strongly
   lumpable) — root sum $= n_U \cdot \hat a_\rho(S) = k_F \cdot 1 = k_F$,
   site lik $= 1$ — also correct. So the missing guard is **not**
   a correctness bug, just a defensive-programming gap.

8. **Singleton-only ascertainment for partitions with $n_U \geq 1$.**
   The collapsed kernel produces $P(\text{const})$ correctly (above);
   $P(\text{singleton})$ uses the *uncollapsed* `singleton_site_prob_jc`
   helper, which is a single per-partition call. This is correct
   because singleton-site probability is a sum of $k_F$ terms with full
   JC symmetry — the collapse provides no measurable savings, and
   re-deriving it under the partition would require a separate proof
   that the "single tip in observed state, all others in (any) other
   state" event has class-constant CL. (It does, by the same Lemma-4
   argument restricted to "all observed" tips. But the gain is not
   worth the additional code path.) ✓ no bug.

9. **`kFull == kObs + 1` (i.e. $n_U = 1$).** Already covered in case 2.
   The kernel runs correctly but provides no speedup; dispatch gates it
   out.

10. **`kEff = 1` (no observed states distinguishable).** Covered in
    case 7. Correct but defensively undefended.

---

## Verdict

**Watertight with caveats.**

The JC(k) lumpability collapse is mathematically correct, and the three
hot-path kernels (`persite_collapsed_impl`, `pruning_jc_flat_collapsed`,
`pruning_jc_acrv_flat_collapsed`, all in `src/mcmc_likelihood.cpp`) match
the derivation line-for-line against the reference
`src/likelihood.cpp::pruning_jc_collapsed`. The fused-ascertainment
extension in `8d4638e` is also correct: the output `outConstProb`
matches the uncollapsed convention (directly usable, no caller-side
$k_F$ multiplication or rate-category division).

**Caveats:**

1. **Lane title mismatch.** The lane is titled "F81-Het collapse" but
   PR #2 implements JC(k) collapse only. The F81-Het case (Assumption 6
   violated) is *not* covered by the lumpability argument here. The
   commit messages and project memory `project_f81_collapse_design`
   acknowledge this; the lane title should be updated to "JC(k) collapse
   / lumpability proof" before this proof is filed as covering the
   intended scope. (No patch — this is an orchestrator/labeling concern,
   not a code issue.) See §6 for what the F81-Het extension would
   require.

2. **Missing `kObs < 1` defensive guard in the three new kernels.**
   The reference `pruning_jc_collapsed` (`src/likelihood.cpp:173–175`)
   has an `Rcpp::stop` guard; the three new kernels do not. The math
   shows this is harmless (degenerate $k_O = 0$ would still produce the
   correct site likelihood of 1), so it is *not a correctness bug*,
   merely a defensive-programming gap. No patch.

3. **No `kFull == kObs` check.** The kernels would silently produce
   wrong dimensions (lumped column has $n_U = 0$ multiplicity, which is
   degenerate but not flagged). The dispatch gates this out via
   `kObsMax + 1 < kFull`, so it cannot occur in practice. Not a bug; no
   patch.

**No patch produced.** Worktree state: no edits made, only
`dev/red-team/proofs/f81-collapse.md` added.
