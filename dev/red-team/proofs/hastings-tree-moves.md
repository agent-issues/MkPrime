# Hastings ratios and detailed balance for tree-topology proposals

**Lane:** L3 (red-team)
**Author:** math-prover (agent run 2026-05-26)
**Scope:** MkPrime tree-topology proposals in `src/tree_moves.cpp`,
`src/proposals.cpp`, and `src/mcmc.cpp`.

## Inventory

Nine topology proposals are present in the codebase. The first four are
"plain" MH proposals (uniform candidate selection, Jacobian-only Hastings
correction); the last five are likelihood-biased Gibbs–MH proposals (the
sampling weights depend on the likelihood; the Hastings ratio inherits the
forward/reverse weight ratio).

| Move | Type | Implementation | Hastings claim |
|------|------|----------------|----------------|
| `ProposeNni`              | symmetric MH | `tree_moves.cpp:23–102`     | $\log H = 0$ |
| `ProposeSpr`              | MH+Jacobian  | `proposals.cpp:22–167`      | $\log H = \log \ell_{\text{reg}}/\ell_{\text{merge}}$ |
| `ProposeTbr`              | MH+Jacobian  | `tree_moves.cpp:259–500`    | $\log H = \log(\ell_{\text{reg}}/\ell_{\text{merge}}) + \log(\ell_{\text{sub}}/\ell_{\text{mergeSub}})$ |
| `swap_subtrees_impl`      | symmetric MH | `tree_moves.cpp:130–172`    | $\log H = 0$ |
| `gibbs_spr_impl`          | Gibbs–MH     | `mcmc.cpp:899–1170`         | exact reversibility (Gibbs) |
| `gibbs_subtree_swap_impl` | Gibbs–MH     | `mcmc.cpp:1654–1878`        | exact reversibility (Gibbs) |
| `weighted_spr_impl`       | Gibbs–MH + Jacobian | `mcmc.cpp:2551–2807` | branch-fraction component only |
| `weighted_subtree_swap_impl` | Gibbs–MH + Jacobian | `mcmc.cpp:2825–3009` | branch-fraction component only |
| `pspr_proposal_impl`      | parsimony-biased MH | `mcmc.cpp:3288–3465`        | Jacobian + parsimony bias correction |

Proofs are organised in the same order.

---

## Assumptions (shared)

A1. **Tree representation.** Trees are stored as a `phylo` edge matrix in
canonical preorder (TreeTools convention). Tips have node indices
$1 \dots n$; internal nodes have indices $n+1 \dots 2n-1$; the root is
node $n+1$. The root is a *trifurcation* (three children), matching the
output of `ape::unroot()` followed by `TreeTools::Preorder()`. Every
non-root internal node has exactly two children. (Verified empirically for
$n \in \{4,\dots,8\}$.)

A2. **Edge count.** For an unrooted binary tree on $n \geq 3$ tips, the
edge matrix has $E = 2n - 3$ rows ($n$ pendant edges + $n - 3$ internal
edges).

A3. **Branch lengths.** Stored as `rel_br_lengths`, a simplex summing to 1
that, multiplied by `tree_length`, gives absolute edge lengths. Total tree
length is held fixed by every topology move; topology moves redistribute
mass on the simplex.

A4. **Topology prior.** Flat over labelled unrooted binary topologies. This
matches the convention under which the chi-squared test in
`tests/testthat/test-tbr-detailed-balance.R` is well-posed.

A5. **Branch-length prior.** Treated as a measurable function of the
simplex point and tree length; all relevant priors are jointly absolutely
continuous in the simplex coordinates, so the Jacobian factors below are
unambiguous.

A6. **Selection of internal edges.** "Internal edge" in this code means an
edge row $i$ with `parent[i] > nTip && child[i] > nTip`. Under A1 the root
has parent = 0 in the edge matrix, so root–incident edges have
$\text{parent}=n+1>n$ and may qualify if their child is internal. The
internal-internal edge count thus equals $n - 3$ — the textbook count of
internal edges in an unrooted binary tree (Lemma N1 below).

---

# 1. NNI proposal (`ProposeNni`)

## Theorem (NNI symmetry).
The NNI proposal in `tree_moves.cpp:23–102` is a symmetric proposal on
unrooted labelled binary topologies. Equivalently, $\log H_{\text{NNI}} = 0$.

## Lemma N1 (internal-edge count).
Under A1–A2, the number of edge rows with both endpoints internal equals
$n - 3$.

*Proof.* The unrooted binary tree on $n$ tips has $n - 3$ internal edges
and $n$ pendant edges. The trifurcating root introduces no extra edges
(it sits *at* a node of the unrooted tree, not on an edge). Each pendant
edge has a tip endpoint (label $\leq n$). Each internal edge has two
internal endpoints (labels $> n$). $\square$

*Empirical confirmation.* $n=4,\dots,8$: edges = $\{5,7,9,11,13\}$ =
$\{2n-3\}$; internal-internal = $\{1,2,3,4,5\}$ = $\{n-3\}$.

## Lemma N2 (degree of selectable endpoints).
Let an internal edge $(u, v)$ be selected by `nni_proposal_impl`. Then:

- $v$ is non-root (since $\text{parent}=u$), so $v$ has exactly two
  children. `vChildRows` has size 2.
- $u$ has $|\text{children}(u)| - 1$ "siblings" of $v$:
    - If $u = \text{root}$: $u$ has 3 children, so `uSibRows` has size 2.
    - If $u \neq \text{root}$: $u$ has 2 children, so `uSibRows` has size 1.

## Proof of NNI symmetry.

Let $T$ be the current tree, $T'$ the proposed tree. The forward proposal
performs three independent uniform draws:

1. Pick an internal edge $(u, v)$ uniformly from $|N(T)| = n-3$ edges
   (`tree_moves.cpp:44–48`).
2. Pick $c \in \text{children}(v)$ uniformly (size 2)
   (`tree_moves.cpp:70–71`).
3. Pick $w \in \text{children}(u) \setminus \{v\}$ uniformly (size 1 if
   $u$ is non-root; size 2 if $u$ is root) (`tree_moves.cpp:72–73`).

Define
$$
m(u) \;=\; \begin{cases} 1 & \text{if } u \neq \text{root}, \\ 2 & \text{if } u = \text{root.} \end{cases}
$$
Then
$$
q(T' \mid T) \;=\; \sum_{(u,v,c,w) \text{ producing } T'} \frac{1}{n-3} \cdot \frac{1}{2} \cdot \frac{1}{m(u)+1}
$$
where the sum runs over selection tuples that produce $T'$.

**Two regimes for the unrooted NNI neighbour:**

*(a) Non-root–incident edge ($u \neq \text{root}$).* The edge $(u,v)$ in
the unrooted view connects two degree-3 vertices. Each has two "off-edge"
neighbours; swapping any of the 2×1 = 2 ordered pairs produces a distinct
unrooted neighbour (since here $m(u)=1$). Multiplicity 1 per neighbour. So
$$q(T' \mid T) = \frac{1}{(n-3) \cdot 2} \;=\; \frac{1}{2(n-3)}.$$

*(b) Root–incident edge ($u = \text{root}$).* Here $u$ has 3 children
($v, r_1, r_2$). The unrooted view sees $u$ as a degree-3 vertex with
off-edge neighbours $r_1, r_2$. There are 2 distinct unrooted NNI
neighbours of $T$ on edge $(u,v)$ — each obtainable by *two* of the
$2 \times 2 = 4$ raw swaps (the swap $(c_1, r_1)$ and $(c_2, r_2)$
produce the same unrooted tree, since they yield the same bipartition).
So
$$q(T' \mid T) \;=\; 2 \cdot \frac{1}{(n-3) \cdot 2 \cdot 2} \;=\; \frac{1}{2(n-3)}.$$

In both regimes, $q(T' \mid T) = 1/[2(n-3)]$.

**Reverse move.** NNI preserves: (i) the existence of edge $(u, v)$ — both
endpoints stay internal post-swap; (ii) the internal-internal count
$|N(T')| = n - 3$; (iii) the root-incidence status of $(u,v)$ in $T'$
($u$ remains the root, $v$ remains its child). Therefore $q(T \mid T')$ is
computed by the same formula, $1/[2(n-3)]$.

$$\log H_{\text{NNI}} \;=\; \log \frac{q(T \mid T')}{q(T' \mid T)} \;=\; 0. \qquad \square$$

## Implementation cross-check (NNI)

| Step | Code |
|------|------|
| Internal-edge enumeration | `tree_moves.cpp:31–35` |
| Uniform pick | `tree_moves.cpp:44–48` |
| `vChildRows` / `uSibRows` | `tree_moves.cpp:53–60` |
| Uniform child / sibling pick | `tree_moves.cpp:70–73` |
| Swap (reassign parents) | `tree_moves.cpp:81–82` |
| Canonical reorder | `tree_moves.cpp:88–96` |
| Hastings return = 0 | `tree_moves.cpp:101` (`logHastings = 0.0`) |

Implementation matches the derivation.

## Edge cases (NNI)

- **$n < 4$:** $n - 3 \leq 1$, so $|N(T)| \in \{0, 1\}$. For $n = 4$,
  there is exactly one internal edge and the move is forced — still
  symmetric (each direction has $q = 1/2$). For $n = 3$, no internal
  edges; the code returns `logHastings = R_NegInf` (rejection)
  (`tree_moves.cpp:37–42`).
- **Branch length 0:** NNI does not modify branch lengths; the absolute
  lengths flow through `preorder_weighted_impl` unchanged. No special
  case needed.

## Verdict (NNI)
**Watertight.** Implementation matches the derivation. The non-trivial
root-incident case is handled correctly by the implementation's uniform
sampling because the multiplicity-2 collapsing happens in *both*
directions identically.

---

# 2. SPR proposal (`ProposeSpr`)

## Theorem (SPR Hastings).
The SPR proposal in `proposals.cpp:22–167` produces moves whose Hastings
ratio is
$$\log H_{\text{SPR}} \;=\; \log \ell_{\text{reg}} - \log \ell_{\text{merge}},$$
where $\ell_{\text{reg}}$ is the absolute length of the chosen regraft
edge in $T$ and $\ell_{\text{merge}} = \ell(p \to u) + \ell(u \to w)$ is
the merged length of the two edges flanking the suppressed node in $T'$
(equivalently, the new regraft edge in the reverse move).

## Lemma S1 (candidate-set sizes).
Let $T$ be an unrooted binary tree on $n$ tips, root = $n+1$, edge count
$E = 2n-3$. Let edge $(u, v)$ be selected for pruning, and let
$d(v) = |\text{desc}(v) \cap \text{tips}|$ be the number of tips below
$v$.

- Eligible prune edges: $E - 2 = 2n - 5$ (exclude the 2 edges
  $\text{root} \to r_1$ and $\text{root} \to r_2$ that violate the
  $\text{parent} \neq \text{root}$ filter at `proposals.cpp:31–35`).
  Wait — under A1 the root is a *trifurcation* with 3 children. Three
  edges have $\text{parent} = \text{root}$, so eligible prune edges =
  $E - 3 = 2n - 6$. **Verify by running:** for $n=5$, $E = 7$, eligible
  = 4.

  Confirmed for $n = 5$: `nEdge = 7`, the 3 edges with parent = node 6
  (the root) are excluded. The Gibbs and TBR variants apply the identical
  filter.

  *General formula:* `eligiblePrune` $= 2n - 6 = E - 3$.

- Regraft candidates: edges with $\text{child}[i] \notin \text{desc}(v)$
  AND $\text{parent}[i] \neq u$ AND $\text{child}[i] \neq u$
  (`proposals.cpp:100–106`). The descendants of $v$ in the *rooted*
  representation consist of $v$ itself plus its subtree. Number of edges
  with child in $\text{desc}(v)$ = $2 d(v) - 2$ (the subtree has
  $d(v)$ tips, so $2 d(v) - 1$ edges; the edge $(u, v)$ has child = $v$
  itself, plus $2 d(v) - 2$ edges strictly inside the subtree, total
  $2 d(v) - 1$). The remaining filter excludes the two edges
  $(p, u)$ and $(u, w)$. Net candidates:
  $$|\text{cands}| \;=\; E - (2d(v) - 1) - 2 \;=\; 2n - 3 - 2d(v) + 1 - 2 \;=\; 2n - 4 - 2d(v).$$
- After SPR, the residual backbone tree (the original tree with $v$
  pruned and $u$ suppressed) has $E - 2$ edges. Regrafting on any
  candidate produces a binary tree $T'$ of the same topology class. In
  $T'$, the reverse move selects pruning edge $(u, v)$ (identical
  node identifiers post-suppression-and-reinsertion-by-the-same-name),
  and the analogous candidate count in $T'$ equals
  $2n - 4 - 2 d_{T'}(v)$. Since $v$'s subtree is moved as a unit,
  $d_{T'}(v) = d_T(v)$. **The candidate-set sizes match.**

## Lemma S2 (branch-length Jacobian).
Let $\ell_{\text{merge}} = \ell(p \to u) + \ell(u \to w)$ be the
combined length of the two edges to be merged when $u$ is suppressed; let
$\ell_{\text{reg}}$ be the length of the chosen regraft edge before
splitting. Let $\tau \sim \mathrm{Uniform}(0, 1)$ (`proposals.cpp:121`).
The forward proposal:

1. *Merge* $(p \to u, u \to w) \mapsto (p \to w)$ with length
   $\ell_{\text{merge}}$ (a deterministic 2-to-1 map).
2. *Split* $(a \to b) \mapsto (a \to u, u \to b)$ with lengths
   $(\tau \ell_{\text{reg}}, (1-\tau) \ell_{\text{reg}})$.

In the space of (sibling edge, parent-of-u edge, regraft edge) absolute
lengths $(\ell_w, \ell_p, \ell_r)$ mapped to (parent edge, u-to-b edge,
new merged edge) $(\ell'_p = \tau \ell_r, \ell'_{ub} = (1-\tau) \ell_r,
\ell'_{\text{merge}} = \ell_w + \ell_p)$ — the transformation has
Jacobian

$$
\left|\frac{\partial(\ell'_p, \ell'_{ub}, \ell'_{\text{merge}}, \tau)}{\partial(\ell_w, \ell_p, \ell_r, \tau)}\right|
= \left|\det \begin{pmatrix} 0 & 0 & \tau & \ell_r \\ 0 & 0 & 1-\tau & -\ell_r \\ 1 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \end{pmatrix}\right| = \ell_r.
$$

The reverse move reverses this transformation, generating a new $\tau'
\in (0, 1)$. Its Jacobian (analogous matrix with $\tau \to \tau'$ and
$\ell_r$ replaced by the *post-move* merged length, which equals
$\ell_{\text{merge}}$) is $\ell_{\text{merge}}$.

## Proof of SPR Hastings.

Combining Lemmas S1–S2, with $|N(T)| = (2n-6) \cdot (2n - 4 - 2 d(v))$
(uniform prune × uniform regraft) and $|N(T')| = |N(T)|$ (S1), we have

$$
q(T' \mid T) \;=\; \frac{1}{|N(T)|} \cdot 1\big|_{\tau \sim U(0,1)}, \qquad q(T \mid T') \;=\; \frac{1}{|N(T)|} \cdot 1\big|_{\tau' \sim U(0,1)}.
$$

The candidate-set ratio is 1. By the change-of-variables formula for the
Metropolis–Hastings acceptance ratio with deterministic component,

$$
\log H = \log \frac{q(T \mid T')}{q(T' \mid T)} + \log \left| \frac{\partial(\text{rev image})}{\partial(\text{fwd image})} \right|^{-1} = 0 + \log \frac{\ell_{\text{merge}}}{\ell_r}\cdot\text{(?)}.
$$

To be careful: Green (1995) §3.3 gives the dimension-matching MH
acceptance as $\alpha = \min\{1, R\}$ where
$$R \;=\; \frac{\pi(T') q(T \mid T')}{\pi(T) q(T' \mid T)} \cdot |J|$$
with $|J|$ the Jacobian of the deterministic part of the transformation
(here, the bijection between auxiliary variable $\tau$ in the forward
move and the new edge lengths). Working with the joint density on
$(\text{topology}, \ell_1, \dots, \ell_E, \tau)$: the forward move maps
$(\ell_w, \ell_p, \ell_r, \tau)$ to
$(\ell_{\text{merge}}, \tau \ell_r, (1-\tau)\ell_r, \tau')$ where
$\tau' \in (0,1)$ is the new auxiliary. The full transformation has
Jacobian magnitude $\ell_r$ (the merge component contributes 1; the
split contributes $\ell_r$ via $\det \begin{pmatrix} \tau & \ell_r \\ 1-\tau & -\ell_r \end{pmatrix} = -\ell_r$).

Hence
$$\log H = \log \ell_r - \log \ell_{\text{merge}} \;=\; \log \ell_{\text{reg}} - \log \ell_{\text{merge}}.$$

This matches Lakner et al. 2008 §2.2 eq. 7 (their "$L_a/L_d$" form). $\square$

## Implementation cross-check (SPR)

| Step | Code |
|------|------|
| Eligible-prune filter | `proposals.cpp:30–35` |
| Uniform prune pick | `proposals.cpp:44–47` |
| BFS descendants of $v$ | `proposals.cpp:79–97` |
| Candidate filter | `proposals.cpp:100–106` |
| Uniform regraft pick | `proposals.cpp:115–118` |
| $\tau \sim U(0,1)$ | `proposals.cpp:121` |
| $\ell_{\text{merge}}$, $\ell_{\text{reg}}$ | `proposals.cpp:129–130` |
| Suppress $u$, split regraft | `proposals.cpp:138–148` |
| Hastings | `proposals.cpp:161` (`logHastings = log(lRegraft) - log(lMerge)`) |

Implementation matches.

## Edge cases (SPR)

- **$n = 4$:** $|\text{eligiblePrune}| = 2n-6 = 2$; for any prune, the
  residual backbone has 3 edges, candidate count = $2n - 4 - 2 d(v) =
  4 - 2d(v)$, which is $\in \{0, 2\}$ depending on whether $v$ is
  internal ($d=2$, cands=0 → rejected by `proposals.cpp:108–113`) or a
  tip ($d=1$, cands=2). So for $n=4$ the only useful prune is of a tip.
- **$d(v) = n - 1$** (prune everything except one tip): cands = $2n - 4 -
  2(n-1) = -2 < 0$ — impossible. In fact $d(v) \leq n - 2$ since the
  prune edge has parent $\neq$ root and root has $\geq 1$ tip
  descendant per side.
- **Zero-length branches:** $\log H$ becomes $-\infty$ if
  $\ell_{\text{merge}} = 0$ or $\ell_{\text{reg}} = 0$. The MCMC harness
  is expected to reject such moves; A5 (priors absolutely continuous in
  simplex coordinates) ensures this is measure-zero.

## Verdict (SPR)
**Watertight.** Implementation matches; Hastings = $\log(\ell_r/\ell_m)$
as derived.

---

# 3. TBR proposal (`ProposeTbr`)

## Theorem (TBR Hastings).
The TBR proposal in `tree_moves.cpp:259–500` produces moves whose
Hastings ratio is
$$\log H_{\text{TBR}} \;=\; \log \frac{\ell_{\text{reg}}}{\ell_{\text{merge}}} \;+\; \log \frac{\ell_{\text{sub}}}{\ell_{\text{mergeSub}}},$$
where the first term is the SPR component (as in §2) and the second is
the subtree re-rooting Jacobian.

## Lemma T1 (re-rooting preserves $|\text{desc}(v)|$).
Re-rooting the pruned subtree at one of its internal edges produces a
new rooted representation with the *same set of tips* below $v$ and the
same total number of edges (i.e. $2 d(v) - 1$). $\square$ (Direct from
the bijection between rooted representations of an unrooted subtree.)

## Lemma T2 (candidate-count balance).
Both the forward and reverse moves have:

- Eligible prune edges = $2n - 6$ (Lemma S1; A1).
- Regraft candidates = $2n - 4 - 2 d(v)$ (Lemma S1 + T1).
- Subtree-edge count for re-rooting = $2 d(v) - 1$ (T1).

All three counts are state-independent given $v$ (which is invariant
across the move). Thus the candidate-count ratio is identically 1.

## Lemma T3 (subtree re-rooting Jacobian).
Re-rooting picks an internal edge $(x, y)$ of the subtree uniformly from
$2 d(v) - 1$ candidates (`tree_moves.cpp:335–337`), draws
$\sigma \sim U(0, 1)$ (`tree_moves.cpp:382`), and applies:

- $\ell(v \to a_1)$ is absorbed into $\ell(a_1 \to c_{\text{other}})$
  (`tree_moves.cpp:386–387`): new length $= \ell(v \to a_1) +
  \ell(\text{other})$.
- The selected edge $(x \to y)$ is split into $\ell(v \to x) =
  \sigma \ell_{\text{sub}}$ and $\ell(v \to y) = (1-\sigma)
  \ell_{\text{sub}}$ (`tree_moves.cpp:391–392, 417–418`).
- Path-reversal of intermediate edges (`tree_moves.cpp:399–413`)
  preserves edge lengths.

The combined Jacobian of this transformation in the space of subtree
edge lengths is $\ell_{\text{sub}}$ (the split contributes $\ell_{\text{sub}}$
exactly as in S2; the absorb step is a measure-preserving sum-merge with
Jacobian 1; path reversal is permutation-with-determinant-$\pm 1$).

In the reverse direction, the merged $a_1 \to c_{\text{other}}$ edge
of length $\ell_{\text{mergeSub}}$ is split back, giving Jacobian
$\ell_{\text{mergeSub}}$.

## Proof of TBR Hastings.
Composition of SPR (Lemma S1–S2) and the re-rooting transformation
(Lemma T1–T3), with the candidate-count ratio = 1 (Lemma T2):
$$\log H_{\text{TBR}} = (\log \ell_{\text{reg}} - \log \ell_{\text{merge}}) + (\log \ell_{\text{sub}} - \log \ell_{\text{mergeSub}}). \qquad \square$$

## Implementation cross-check (TBR)

| Step | Code |
|------|------|
| Eligible-prune filter | `tree_moves.cpp:275–283` |
| BFS descendants of $v$ + subEdgeRows | `tree_moves.cpp:296–313` |
| Pick subtree edge for re-rooting | `tree_moves.cpp:335–339` |
| Identify $a_1$ (path-top) | `tree_moves.cpp:355–372` |
| Draw $\sigma$ | `tree_moves.cpp:382` |
| Absorb $v\to a_1$ into $a_1\to c_{\text{other}}$ | `tree_moves.cpp:386–387` |
| Split $x \to y$ via $\sigma$ | `tree_moves.cpp:391–392, 416–418` |
| Path-reversal (length-preserving) | `tree_moves.cpp:400–413` |
| Regraft (SPR phase) | `tree_moves.cpp:455–474` |
| $\log H$ | `tree_moves.cpp:492–494` |

Implementation matches.

## Minor BFS redundancy (non-bug).
`tree_moves.cpp:306` reads `if (parent[i] == cur && isDesc[parent[i]])`.
Since `parent[i] == cur` and `cur` was popped from the BFS queue (which
only contains marked descendants), `isDesc[parent[i]] == isDesc[cur] ==
true` is tautologically satisfied. The clause is harmless — `subEdgeRows`
is still correctly populated — but the redundancy obscures the
invariant. The SPR analogue at `proposals.cpp:88–95` omits the redundant
clause and is clearer.

No patch required (the code is correct); flagged for stylistic cleanup.

## Edge cases (TBR)

- **$v$ is a tip** (`tree_moves.cpp:436–440`): re-rooting skipped,
  $\ell_{\text{sub}} = \ell_{\text{mergeSub}} = 1$, subtree-Jacobian
  contribution = 0; TBR reduces to SPR. Correct.
- **Subtree has no internal edges** ($d(v) = 1$, equivalent to $v$ tip):
  `nSubEdge = 0`, same path. Correct.
- **$x = v$** (`tree_moves.cpp:420–435`): the chosen subtree edge is a
  direct child of $v$, so no topology change is needed. Code draws
  $\sigma$ (to consume the same RNG count as the general path) and sets
  $\ell_{\text{mergeSub}} = \ell_{\text{sub}}$, giving subtree-Jacobian
  contribution = 0. *Caveat:* the code does **not** actually use
  $\sigma$ to redistribute lengths in this branch — the comment at
  `tree_moves.cpp:421–431` correctly identifies this as a no-op for
  topology and branch lengths. The Hastings ratio (subtree part = 0) is
  correct because no Jacobian factor is introduced. RNG-parity with
  the general case is maintained.

## Verdict (TBR)
**Watertight.** Implementation matches; the candidate-count cancellation
relies on Lemmas T1–T2, which hold under A1–A6. The detailed-balance
test in `tests/testthat/test-tbr-detailed-balance.R` (1M iters, 15-bin
$\chi^2$) is the right empirical check given the proven
$\log H_{\text{TBR}}$ above.

---

# 4. Plain subtree swap (`swap_subtrees_impl`)

## Theorem (subtree-swap symmetry).
The plain subtree swap in `tree_moves.cpp:130–172` is symmetric in
topology *and* branch lengths (the two edges' lengths swap with their
subtrees). $\log H = 0$.

## Proof.
The forward transformation is the involution
$(\text{parent}[\text{rowA}], \text{parent}[\text{rowB}]) \mapsto
(\text{parent}[\text{rowB}], \text{parent}[\text{rowA}])$,
$(\ell_A, \ell_B) \mapsto (\ell_B, \ell_A)$ (`tree_moves.cpp:148–151`).
Involutions have Jacobian $\pm 1$ (here, exact $-1$ on the length
component, $+1$ on the topology); $|J| = 1$.

The candidate set of valid swap partners is determined by
`get_valid_swap_partners_impl` (`tree_moves.cpp:176–226`), which excludes
descendants, ancestors, siblings, and the root partner. This set is
symmetric in nodeA / nodeB:
$$B \in \text{partners}(A) \iff A \in \text{partners}(B)$$
because all four exclusions (descendant, ancestor, sibling, root) are
invariant under the involutive swap (after swap, $A$'s descendants
include the new subtree below $A$, but the *node identifiers* in
$\text{partners}$ refer to the rooted-tree's relations, which are
swapped symmetrically).

If `swap_subtrees_impl` is called from a higher-level proposal that
selects (nodeA, nodeB), the proposal is symmetric iff that selection is
symmetric. `swap_subtrees_impl` itself is the implementation of the
deterministic swap; symmetry holds at the implementation level.

$\square$

## Implementation cross-check
| Step | Code |
|------|------|
| Find rows | `tree_moves.cpp:135–143` |
| Swap parent + length | `tree_moves.cpp:148–151` |
| Reorder | `tree_moves.cpp:158–166` |
| $\log H = 0$ | `tree_moves.cpp:171` |

Implementation matches.

## Verdict (subtree swap)
**Watertight** (subject to symmetric selection of `(nodeA, nodeB)` by the
caller, which holds for the dispatch through `gibbs_subtree_swap_impl` —
see §6).

---

# 5. GibbsSPR (`gibbs_spr_impl`)

## Theorem (GibbsSPR reversibility).
GibbsSPR samples the new topology from the conditional posterior over
the SPR neighbourhood plus the current state ("self") at inverse
temperature $\beta$. This is an exact Gibbs update on the (topology,
branch-fractions-at-regraft-edge) joint, *provided* the candidate-set
symmetry of Lemma S1 holds and the branch-fraction split is at the
fixed midpoint $\tau = 1/2$.

## Construction (`gibbs_spr_impl`).
1. Pick prune edge $(u, v)$ uniformly (`mcmc.cpp:916–918`).
2. Enumerate candidate regraft edges (Lemma S1; `mcmc.cpp:949–956`).
3. For each candidate, evaluate the likelihood at $\tau = 1/2$
   (`mcmc.cpp:1078, 1149–1152`).
4. Sample one of $\{$self, cand$_1$, ..., cand$_N\}$ proportional to
   $\exp(\beta \cdot \log L)$ (`mcmc.cpp:1120–1141`).

## Reversibility.
This is a *systematic-scan Gibbs* update on the joint conditional
$$\pi(T' \mid T_{\text{prune-residual}}, \tau = 1/2),$$
where $T_{\text{prune-residual}}$ is the unique unrooted backbone after
removing $v$'s subtree and suppressing $u$. By Lemma S1, this residual
backbone is identical in forward and reverse. Provided the prior on $T'$
*given* the prune residual is uniform over the candidate set ∪ {self},
the Gibbs draw is reversible by construction (Geman & Geman 1984; see
also Yang 2014 §8.3).

**Caveat (mild).** The fixed $\tau = 1/2$ means GibbsSPR is reversible
on the *topology marginal* conditional on a hard-coded mid-edge regraft
position. Since the *branch-fraction* component is deterministic
($\tau = 1/2$ both directions), the move targets a restricted invariant
distribution — one where regrafts happen at the mid-edge. The chain is
*not* irreducible on the full joint without additional branch-length
moves; combined with the simplex Dirichlet / branch-scale moves
(`proposals.cpp` Dirichlet, `mcmc.cpp:2264` weighted_branch_scale_impl),
this is unproblematic.

## Implementation cross-check (GibbsSPR)
| Step | Code |
|------|------|
| Prune-edge sample | `mcmc.cpp:916–918` |
| BFS desc | `mcmc.cpp:932–947` |
| Candidate filter | `mcmc.cpp:949–956` |
| Self-weight | `mcmc.cpp:1125, 1127` |
| Cand weights | `mcmc.cpp:1128–1131` |
| Sample | `mcmc.cpp:1134–1141` |
| Commit | `mcmc.cpp:1144–1166` |

## Verdict (GibbsSPR)
**Watertight with caveats.** Reversibility holds *conditional on* the
shared prune-residual backbone (Lemma S1) and a flat conditional prior on
the regraft edge. The hard-coded $\tau = 1/2$ restricts the joint
invariant distribution; this is compensated by separate branch-length
moves in the chain.

---

# 6. GibbsSubtreeSwap (`gibbs_subtree_swap_impl`)

## Theorem (GibbsSubtreeSwap reversibility).
For a randomly chosen node $A$ and the enumerated valid partners
$\text{partners}(A)$, GibbsSubtreeSwap samples uniformly with weights
$\exp(\beta \log L)$. Reversibility holds *provided*
$$\sum_{A: B \in \text{partners}(A)} \mathbf{1}[B = \text{partners}(A)_i] \;=\; \sum_{A': A \in \text{partners}(A')} \mathbf{1}[A = \text{partners}(A')_i].$$

## Lemma G1 (partner-set symmetry).
$B \in \text{partners}(A) \iff A \in \text{partners}(B)$
(`tree_moves.cpp:176–226`).

*Proof.* "Descendant of $A$ excluded" $\iff$ "ancestor of $A$ excluded
from $B$'s perspective" — the descendant/ancestor relation is the
inverse of itself across the pair. Sibling: $A, B$ siblings is symmetric
in $A, B$. Root: neither $A$ nor $B$ may be the root. All four
exclusions are symmetric. $\square$

## Detailed-balance.
By G1 the forward selection $(A, B)$ from the marginal $A$-distribution
(uniform over `nEdge` edge-child nodes) and the conditional uniform on
$\text{partners}(A)$, followed by the Gibbs weight, has the property
$$q(T' \mid T) = \frac{1}{\text{nEdge}} \cdot \frac{1}{|\text{partners}(A)|} \cdot \frac{e^{\beta L(T')}}{Z(A; T)},$$
where $Z(A; T) = e^{\beta L(T)} + \sum_{B' \in \text{partners}(A)}
e^{\beta L(T'(A, B'))}$.

For the reverse move on $(B, A)$ in $T'$, an analogous expression with
$|\text{partners}_{T'}(B)|$ and $Z(B; T')$ appears. By G1, $A \in
\text{partners}_{T'}(B) \iff B \in \text{partners}_T(A)$, but the
*sizes* $|\text{partners}_T(A)|$ and $|\text{partners}_{T'}(B)|$ are
not generally equal: a swap can change the ancestor/descendant structure
of nodes far from $A, B$.

**Therefore detailed balance is NOT trivially achieved.** This is a
known issue: GibbsSwap with asymmetric neighborhood needs a Hastings
correction $\log(|\text{partners}_T(A)| / |\text{partners}_{T'}(B)|) +
\log(\text{nEdge}_T / \text{nEdge}_{T'})$. The latter is 0
(`nEdge` invariant). The former is **not** in the code (`mcmc.cpp:1874`
sets `state->logLik = candLL[chosen]` and there's no Hastings factor).

**Resolution.** A closer look: the Gibbs ratio $w_{T'}/Z(A; T)$ in
forward and $w_T/Z(B; T')$ in reverse, combined with the
proposal-selection probabilities, gives the detailed-balance equation
$$\pi(T) q(T'|T) = \pi(T') q(T|T') \iff \pi(T) \frac{w_{T'}}{Z(A; T)} = \pi(T') \frac{w_T}{Z(B; T')}.$$
Under $\pi(T) \propto w_T = e^{\beta L(T)}$:
$$w_T \frac{w_{T'}}{Z(A;T)} \stackrel{?}{=} w_{T'} \frac{w_T}{Z(B; T')},$$
i.e. $Z(A;T) = Z(B; T')$. This requires the candidate-set $\{T\} \cup
\{T(A, B'') : B'' \in \text{partners}(A)\}$ in the forward and
$\{T'\} \cup \{T'(B, A'') : A'' \in \text{partners}_{T'}(B)\}$ in the
reverse to be *the same set of weighted topologies*. By G1 they're
"symmetric" but not generally equal — they differ in the included
partner sets when the swap changes ancestor/descendant relations.

**Therefore strictly the GibbsSubtreeSwap proposal is NOT detail-balanced
in the partner-set sense.** Whether the actual implementation produces a
correct invariant chain depends on the higher-level acceptance.
Inspection of `mcmc.cpp:1654–1878` shows there is **no** MH
accept/reject step — the chosen partner is *always committed*
(`mcmc.cpp:1874`). This is a true Gibbs draw, and reversibility requires
$Z(A; T) = Z(B; T')$.

**Open question.** I have not closed whether the partner-set sizes are
actually invariant under a single swap. Empirical: for $n = 5$ the
detailed-balance test on TBR passes; no analogous test for
GibbsSubtreeSwap exists in the slow-test suite (only
`tests/testthat/test-subtree-swap.R` which tests the plain
`swap_subtrees_impl`, and the determinism tests). A targeted simulation
would resolve this.

## Verdict (GibbsSubtreeSwap)
**Open question.** Cannot close detailed balance for the partner-set
sampling without either (a) proving partner-set count invariance under
swap, or (b) confirming that the existing tests cover this case
empirically.

*Recommendation:* extend `test-tbr-detailed-balance.R` to include a
`.flat_posterior_chain(MkPrime:::PlainSubtreeSwapDispatcher, ...)` and a
`MkPrime:::GibbsSubtreeSwap` dispatcher at $\beta = 0$ (flat
likelihood), and verify uniform stationarity over the 15 unrooted
topologies on 5 tips. Note this requires exposing a hook at $\beta = 0$
since the Gibbs sampler in MkPrime uses temperature $\beta = 1$ by
default.

---

# 7. WeightedSPR (`weighted_spr_impl`)

## Construction.
For each SPR-candidate regraft edge $i$, marginalise the likelihood
over $B$ branch-fraction bins (Beta-quantile midpoints) to produce a
marginal weight $M_i = \sum_{b=1}^{B} \exp(\beta L(T_i, \tau_b))$. Include
"self" as a separate point weight $M_{\text{self}} = \sum_b \exp(\beta
L(T, \tau_b))$ varying the *current* parent/sibling fraction.

Sample topology proportional to $M$; then sample bin within chosen
topology; then draw $\tau' \sim \mathrm{Beta}(\alpha_{\text{bin}},
\beta_{\text{bin}})$ centred on the bin midpoint.

## Theorem (WeightedSPR Hastings, claimed).
$$\log H = \log w_{\text{self}, b_{\text{old}}} + \log \mathrm{Beta}(f_{\text{old}}; \alpha_{\text{old}}, \beta_{\text{old}})
            - \log w_{\text{chosen}, b_{\text{new}}} - \log \mathrm{Beta}(f_{\text{new}}; \alpha_{\text{new}}, \beta_{\text{new}}).$$
The topology-selection components $M_{\text{chosen}} / \sum M$ in
forward and $M_{\text{self}_{T'}} / \sum M_{T'}$ in reverse are claimed
to cancel.

## Cancellation argument (claimed, examined).
The claim relies on:
(C1) The candidate set in forward and reverse comprises the same
underlying *unrooted topologies* (Lemma S1 — yes, holds).
(C2) The marginal weights $M_i$ are functions only of the resulting
topology and the residual backbone — i.e. $M_i$ in the forward direction
equals $M_{j(i)}$ in the reverse direction for the corresponding partner
index $j(i)$.
(C3) The "self" weight $M_{\text{self}}$ in forward, computed from the
*current* parent/sibling lengths varied across bins, equals the weight
of the chosen-candidate position in the reverse direction.

(C1) is Lemma S1. (C2) follows because $M_i$ integrates over
$\tau \in (0,1)$ at the regraft edge holding all other branch lengths
fixed (`mcmc.cpp:2666–2676`); the integration variable is local to that
edge, and the bin grid is the same in both directions.

(C3) is the **subtle** point. In forward, $M_{\text{self}}$ is computed by
varying $\tau \in \{$bin midpoints$\}$ at the (parent, sibling) edges of
$u$ (`mcmc.cpp:2624–2634`). In reverse, the "self" weight in $T'$ would
be computed similarly at $T'$'s (parent, sibling) edges of $u$ (after
the swap). These edges are *not* the same physical edges in the unrooted
tree — but the *marginalised likelihood* of the unrooted tree $T$ given
the residual backbone $T_{\text{res}}$ is the same in either direction.
The bin grid and Beta concentration are fixed parameters, so the
integration produces the same marginal weight.

**Conclusion.** The topology-selection cancellation is *plausible* but
depends on equating residual-backbone-conditional marginal weights
across the prune/regraft direction. I believe (C3) holds but have not
written a rigorous bijective argument; the empirical
`test-weighted-spr.R` provides partial coverage but does not directly
verify detailed balance.

## Implementation cross-check (WeightedSPR)
| Step | Code |
|------|------|
| Prune sample | `mcmc.cpp:2568–2571` |
| Self-marginal eval | `mcmc.cpp:2622–2634` |
| Candidate marginals (in-place + preorder) | `mcmc.cpp:2653–2686` |
| Mixed-sampling | `mcmc.cpp:2720–2727` |
| Bin sample | `mcmc.cpp:2730–2738` |
| Final fraction $f_{\text{new}}$ | `mcmc.cpp:2745–2747` |
| $\log H$ | `mcmc.cpp:2778–2781` |
| MH accept | `mcmc.cpp:2792–2805` |

The MH step at `mcmc.cpp:2794` provides a *real* accept/reject — so
this is a Gibbs–MH (not pure Gibbs). Therefore detailed balance does
not require exact topology-cancellation; any approximation error in the
proposal density is corrected by the MH ratio. As long as the *claimed*
$\log H$ matches the *actual* ratio of $q(T|T')/q(T'|T)$, the chain is
reversible.

## Verdict (WeightedSPR)
**Watertight with caveats.** The MH step at `mcmc.cpp:2794` ensures
correctness *modulo* the accuracy of `logHR`. Cancellation of the
topology component (claim C3) is plausible and the empirical test
suite does not flag detectable bias, but a rigorous bijection between
forward and reverse marginal weights is not written down here.
*This is an open sub-question worth a dedicated check.*

---

# 8. WeightedSubtreeSwap (`weighted_subtree_swap_impl`)

Analogous to WeightedSPR. The MH step at `mcmc.cpp:2996` corrects for
proposal-density mismatches.

## Verdict
**Watertight with caveats** (same reasoning as §7). The partner-set
symmetry issue from §6 (Open question) is sidestepped because the
WeightedSubtreeSwap uses an MH accept/reject; any bias in proposal
density is corrected by the acceptance ratio. **However**, this requires
the `logHR` formula at `mcmc.cpp:2980–2983` to actually equal
$\log q(T|T') - \log q(T'|T)$. The argument parallels §7 (C1–C3) but
also needs the partner-set size correction:

$$\log H_{\text{topo}} = \log \frac{|\text{partners}_{T'}(B)|}{|\text{partners}_T(A)|}$$

I do not see this term in the code (`mcmc.cpp:2980–2983` only contains
the branch-fraction component). **Potential bug:** if partner-set sizes
differ across a swap, the omitted term is non-zero. This may produce a
biased chain.

Concrete check: pick a node $A$ in a 5-tip tree, perform a swap with
partner $B$, recompute $|\text{partners}_{T'}(B)|$ vs $|\text{partners}_T(A)|$.

Empirical from the existing test suite: `test-weighted-subtree-swap.R`
likely covers basic functionality but not detailed balance at $\beta=0$.

## Verdict (WeightedSubtreeSwap)
**Open question / Potential implementation bug.** The Hastings ratio
appears to omit a $\log(|\text{partners}_{T'}(B)|/|\text{partners}_T(A)|)$
correction. Whether this matters numerically depends on how often
partner-set sizes differ across a swap. **No patch attached** — needs
empirical confirmation before claiming a bug.

---

# 9. pSPR (`pspr_proposal_impl`)

## Theorem (pSPR Hastings).
$$\log H_{\text{pSPR}} = \underbrace{\log \frac{\ell_{\text{reg}}}{\ell_{\text{merge}}}}_{\text{Jacobian (S2)}} + \underbrace{\log w_{\text{orig}} - \log w_{\text{chosen}} + \log Z_{\text{fwd}} - \log Z_{\text{rev}}}_{\text{parsimony bias correction}},$$
where $w_i = \exp(-\alpha (\text{Fitch}(T_i) - \text{Fitch}_{\min}))$,
$Z_{\text{fwd}} = \sum_{i \in \text{cands}} w_i$, and
$Z_{\text{rev}} = Z_{\text{fwd}} + w_{\text{orig}} - w_{\text{chosen}}$.

## Lemma P1 (residual-backbone invariance).
The pSPR forward and reverse moves share the same prune-residual
backbone tree (Lemma S1, T1). The set of candidate regraft edges in
$T_{\text{res}}$ is the same in both directions; thus the Fitch scores
$\text{Fitch}(T_i)$ enumerated over candidates are identical in the
forward and reverse direction — only the *excluded* element differs
(forward excludes the original-position virtual candidate; reverse
excludes the chosen-position candidate).

## Proof of pSPR Hastings.
By P1, the forward weights $\{w_i\}_{i \in \text{cands}}$ and the
reverse weights $\{w_j\}_{j \in \text{cands}'}$ are functions of the
same enumerated Fitch scores. The forward proposal probability for
choosing index $c$:
$$q(T' \mid T) = \frac{1}{|N_T|} \cdot \frac{w_c}{Z_{\text{fwd}}} \cdot 1\big|_{\tau \sim U(0,1)}.$$
The reverse proposal selects the same prune edge in $T'$ (preserved by
SPR — the edge $(u, v)$ exists in both), recovers the same residual
backbone, and chooses the original-position regraft (call it index $o$):
$$q(T \mid T') = \frac{1}{|N_{T'}|} \cdot \frac{w_o}{Z_{\text{rev}}} \cdot 1\big|_{\tau' \sim U(0,1)},$$
where $Z_{\text{rev}}$ now includes $w_o$ (which represents the
*original-position* candidate in the reverse direction) and excludes
$w_c$. So $Z_{\text{rev}} = (Z_{\text{fwd}} - w_c) + w_o = Z_{\text{fwd}}
+ w_o - w_c$.

Combining with the Jacobian (S2):
$$\log H = \log \frac{q(T|T')}{q(T'|T)} + \log \ell_{\text{reg}} - \log \ell_{\text{merge}}
       = \log w_o - \log w_c + \log Z_{\text{fwd}} - \log Z_{\text{rev}} + \log \ell_{\text{reg}} - \log \ell_{\text{merge}}. \square$$

## Implementation cross-check (pSPR)
| Step | Code |
|------|------|
| Fitch scoring + alpha | `mcmc.cpp:3286, 3376–3387` |
| $w_i = e^{-\alpha(s_i - s_{\min})}$ | `mcmc.cpp:3397–3401` |
| $\log H_{\text{brlen}}$ | `mcmc.cpp:3454` |
| $Z_{\text{rev}} = Z + w_{\text{orig}} - w_{\text{chosen}}$ | `mcmc.cpp:3456` |
| $\log H_{\text{pars}}$ | `mcmc.cpp:3457–3458` |

Implementation matches the derivation.

## Caveats (pSPR)
- The constant $\alpha = 0.1$ (`mcmc.cpp:3286`) controls bias strength;
  reversibility holds for any $\alpha$.
- `|N_T| = |N_{T'}|` by Lemmas S1, T1 (the prune-edge count is preserved,
  and the residual-backbone candidate count is preserved). So the
  $1/|N|$ factors cancel.
- The "original position" virtual candidate is *not enumerated* in the
  forward candidates (the forward only enumerates "different" regraft
  positions). But the reverse moves *would* enumerate the
  original-position as one of its candidates (the position where $u$
  was originally inserted in $T$). The code's treatment of `scoreOrig`
  and `wOrig` (`mcmc.cpp:3387, 3455`) handles this.

## Verdict (pSPR)
**Watertight.** Implementation matches.

---

# Composite verdict

| Move | Verdict |
|------|---------|
| NNI | Watertight |
| SPR | Watertight |
| TBR | Watertight |
| swap_subtrees_impl | Watertight |
| GibbsSPR | Watertight with caveats (fixed $\tau = 1/2$) |
| GibbsSubtreeSwap | **Open question** — partner-set asymmetry not closed |
| WeightedSPR | Watertight with caveats (topology-cancellation argument informal) |
| WeightedSubtreeSwap | **Open question / potential bug** — partner-set size correction may be missing |
| pSPR | Watertight |

## Headline finding

Two open questions deserve follow-up:

1. **GibbsSubtreeSwap (`mcmc.cpp:1654–1878`)** commits the chosen
   partner unconditionally without an MH step. Detailed balance
   requires $Z(A; T) = Z(B; T')$, which in turn requires
   $|\text{partners}_T(A)| = |\text{partners}_{T'}(B)|$ — not
   generally true after a swap. The pure-Gibbs form may be biased.

2. **WeightedSubtreeSwap (`mcmc.cpp:2825–3009`)** has an MH step, but
   its `logHR` (lines 2980–2983) omits the
   $\log(|\text{partners}_{T'}(B)|/|\text{partners}_T(A)|)$ correction.
   Whether this matters depends on the empirical distribution of
   partner-set size changes per swap.

A 5-tip flat-likelihood test analogous to `test-tbr-detailed-balance.R`
applied to `GibbsSubtreeSwap` (at $\beta = 0$) and
`WeightedSubtreeSwap` (at $\beta = 0$) would resolve both.

## Patch status

No patch attached. The pure SPR / NNI / TBR / pSPR moves are correctly
derived and implemented; the two flagged moves require analytic or
empirical follow-up before producing a bug-fix patch.

## References

- Lakner C, van der Mark P, Huelsenbeck JP, Larget B, Ronquist F. 2008.
  Efficiency of Markov chain Monte Carlo tree proposals in Bayesian
  phylogenetics. *Syst Biol* 57(1):86–103.
- Felsenstein J. 2004. *Inferring Phylogenies*. Sinauer. §27 (Bayesian
  inference).
- Green PJ. 1995. Reversible jump Markov chain Monte Carlo computation
  and Bayesian model determination. *Biometrika* 82(4):711–732.
- Aldous DJ. 2000. Mixing time for a Markov chain on cladograms.
  *Combinatorics, Probability and Computing* 9(3):191–204.
- Geman S, Geman D. 1984. Stochastic relaxation, Gibbs distributions, and
  the Bayesian restoration of images. *IEEE PAMI* 6(6):721–741.
- Yang Z. 2014. *Molecular Evolution: A Statistical Approach*. OUP. §8.3.
