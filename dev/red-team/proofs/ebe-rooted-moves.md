# Rooted-tree Hastings correctness of the ecology-mode topology moves (EBE)

**Lane:** math-prover. **Date:** 2026-05-29. **Scope:** the topology moves active in
ecology-aware (EBE) mode, asked as: *is each move's stated `logHastings` a correct
proposal-density ratio `log[q(R'→R)/q(R→R')]` — including branch-length
change-of-variables Jacobians — for a reversible MH sampler over the space of **rooted**
trees-with-branch-lengths under a **non-root-invariant** target?*

This question is load-bearing because EBE is non-stationary: re-rooting changes the
likelihood when `z≠0` (test contract T3). The relative proposal weight among different
rootings of one unrooted topology was **untestable** under the old reversible Mk model
(likelihood could not distinguish rootings), so these ratios have only ever been
exercised where their rooted component was invisible. A wrong rooted `q`-ratio biases
root inference invisibly to all legacy tests.

---

## Theorem (informal)

Each topology move active by default in ecology mode — **NNI**, **SPR**, **TBR**, **pSPR**
— supplies a `logHastings` value equal to `log[q(R'→R)/q(R→R')]` as a proposal over
**rooted** binary trees-with-branch-lengths, where `q` includes both the discrete
prune/reroot/regraft selection probabilities and the branch-length change-of-variables
Jacobian. Consequently, the standard MH acceptance probability
`min{1, [π(R')/π(R)]·exp(logHastings)}` defines a kernel reversible w.r.t. the EBE
posterior `π` over rooted trees, even though `π` is not root-invariant.

The two non-default moves **weighted_spr** and **weighted_subtree_swap** are *not*
registered in the default ecology move set; their rooted correctness is assessed
separately and tagged for follow-up.

---

## Assumptions

1. **State space.** The inferential object is a rooted binary tree on `nTip` labelled
   tips, with one distinguished **root node** of out-degree 2, plus a positive
   branch length on each of the `nEdge = 2·nTip − 2` edges. Branch lengths are stored
   as `relBrLengths` (a simplex) × a scalar `treeLength`; the moves here hold
   `treeLength` fixed and permute/transform `relBrLengths` (equivalently absolute
   lengths `absLen = treeLength · relBrLengths`).
2. **Edge representation.** A tree is a `(parent, child)` integer-vector pair. The
   root is the unique node that appears as a `parent` but never as a `child`
   (`TreeTools::renumber_tree.h:204`, `:273`: `if (!data.parent_of[i]) root_node = i`).
   Tips are nodes `1..nTip`; the root is node `nTip+1`; internal nodes are `> nTip`.
3. **Canonical reordering is a root-preserving relabelling.** Every move finishes by
   calling `TreeTools::preorder_weighted_impl(parent, child, absLen)`
   (`src/tree_moves.cpp:99,169,489`; `src/proposals.cpp:162`;
   `src/mcmc.cpp:3665,2901`). This relabels internal nodes into canonical preorder and
   carries branch lengths along **without changing which vertex is the root** (the root
   is recomputed by the same "parent-never-a-child" rule and the traversal is seeded
   from it). It is a bijection on rooted-trees-with-branch-lengths whose Jacobian is the
   identity (a pure relabelling + permutation of edge rows; lengths are carried, not
   transformed). *Verified empirically below: across thousands of moves the root node
   and its degree-2 status are invariant.*
4. **Prior is root-invariant** (established upstream; `dev/ecology/prior-reroot-check.R`,
   `ebe-spec.md §7`). The branch term is the constant `lfactorial(nEdge−1)` (flat
   Dirichlet) and the tree-length term is a Gamma on the preserved total. **This proof
   does not re-audit the prior.** It does mean the rooted target's topology factor is
   uniform on labelled rooted binary trees, so the entire directional pull comes from
   the EBE likelihood.
5. **Move Hastings ratios are likelihood-independent** (structural). The exported
   proposals return `logHastings` computed only from branch lengths and candidate
   counts; the EBE likelihood ratio is applied by the outer MH (`ebe-spec.md §7`).
   *(Exception: pSPR uses Fitch parsimony scores of the data in its proposal; addressed
   in its section. The two weighted_* moves fold the likelihood into their own internal
   accept/reject; addressed in their section.)*
6. **Auxiliary draws are `U(0,1)`** (`ChainRng::unif`), so they have density 1 on the
   unit interval and contribute no extra factor beyond the change-of-variables Jacobian.

If Assumption 3 failed (i.e. the reorder re-rooted the tree or rescaled lengths), every
verdict below would be void. It is checked explicitly and holds.

---

## Active ecology move set (scope confirmation)

Default MCMC control flags (`R/MkPrimeMCMC.R:343–352`):
`gibbsSpr=TRUE`, `gibbsSubtreeSwap=TRUE`, `tbr=TRUE`, `pSpr=TRUE`,
`weightedBranchScale=FALSE`, `weightedSpr=FALSE`, `weightedSubtreeSwap=FALSE`,
`blockGibbsBranch=FALSE`.

Move-weight registration (`R/RunMkPrime.R:3482–3539`): NNI and SPR are always added
(when `!fixTopology && nEdge≥5`); TBR/pSPR/weighted_*/gibbs_* are each behind an
`isTRUE(mcmc$<flag>)` guard. Eco self-disable: `gibbs_spr_impl` and
`gibbs_subtree_swap_impl` `return false` immediately when `data->ecologyAware`
(`src/mcmc.cpp:1024`, `:1791`), so although they carry nonzero weight they never fire in
eco mode.

**Therefore the default-active ecology topology moves are exactly:**

| name | R type | moveType (mcmc.cpp) | proposal impl |
|------|--------|---------------------|---------------|
| NNI  | `nni`  | 5  | `nni_proposal_impl` `src/tree_moves.cpp:33` |
| SPR  | `spr`  | 6  | `spr_proposal_impl` `src/proposals.cpp:32` (eco always uses full-eval branch: `src/mcmc.cpp:5049` requires `!ecologyAware` for the partial-CL path) |
| TBR  | `tbr`  | 17 | `tbr_proposal_impl` `src/tree_moves.cpp:270` |
| pSPR | `pspr` | 20 | `pspr_proposal_impl` `src/mcmc.cpp:3513` |

(The brief's "pSPR `mcmc.cpp:3513`" and "TBR `tree_moves.cpp:270`" name the *impl*
functions; the *dispatch* cases are 20 and 17 respectively — `src/mcmc.cpp:5162`,
`:5070`.)

**Non-default, opt-in, eco-capable:** `weighted_branch_scale` (12),
`weighted_spr` (13, `src/mcmc.cpp:2692`), `weighted_subtree_swap`
(14, `src/mcmc.cpp:2973`). Unlike gibbs_spr/gibbs_subtree_swap these have **no eco
self-disable**, because they re-evaluate the *full* eco likelihood via
`compute_full_loglik_at` and do their own internal MH (`return weighted_*_impl(...)` at
`src/mcmc.cpp:5135,5141,5144`). They simply receive no move weight by default.
Assessed in the final section.

---

## Statement (formal)

Let `S` be the set of rooted binary trees on `nTip` tips with positive edge lengths and
fixed total `treeLength`, equipped with the product of counting measure (topologies) and
Lebesgue measure (the branch-length simplex face). For a move `m` with proposal kernel
`q_m(R→·)`, write its returned value as `H_m(R→R')`. We prove, for
`m ∈ {NNI, SPR, TBR, pSPR}`:

$$ \exp\big(H_m(R\to R')\big) \;=\; \frac{q_m(R'\to R)}{q_m(R\to R')} \qquad \text{for all } R\xrightarrow{m}R' \text{ with } H_m>-\infty. $$

Each `q_m` factorises as `q_m = q_m^{\text{disc}} \cdot q_m^{\text{cont}}`, where
`q_m^{\text{disc}}` is the probability of the discrete selection (prune edge, any reroot
edge, regraft edge) and `q_m^{\text{cont}}` is the density of the new branch lengths
given the discrete choice and the `U(0,1)` auxiliaries. We show
`q_m^{\text{disc}}(R→R') = q_m^{\text{disc}}(R'→R)` (discrete symmetry) and that
`q_m^{\text{cont}}(R'→R)/q_m^{\text{cont}}(R→R')` equals `exp(H_m)` via an explicit
change-of-variables Jacobian.

---

## Proof

### Lemma A (root preservation; Assumption 3 made precise)

`preorder_weighted_impl` (`renumber_tree.h:246–314`) computes
`parent_of[child_i] = parent_i` for every edge, then identifies `root_node` as the
unique node with `parent_of == 0` (`:273`), seeds the preorder traversal at that node
(`:309`, `PreorderState(... root_node ...)`), and emits edges in visit order while
carrying `wt_above[child]` as the edge weight (`:158`). No edge is added or removed; no
length is altered; the root vertex is unchanged. Hence the map `T ↦ preorder(T)` is a
bijection on `S` that fixes the root and has identity Jacobian on branch lengths. ∎

*Numerical confirmation.* On a 6-tip tree, 2000 draws each of SPR/TBR/NNI all return
trees whose root is node 7 with degree 2 (the input's root), i.e. distinct root-nodes =
{7}, root-degrees = {2}. So all three moves are genuinely maps on **rooted** trees, never
silently re-rooting via the relabel step.

### Lemma B (the stochastic branch-split Jacobian is `lRegraft/lMerge`)

This is the standard "stochastic SPR/NNI" branch redistribution (e.g. Lakner et al. 2008,
*Syst. Biol.* 57:86; the BEAST/RevBayes subtree-slide-with-split family). Consider the
local operation that **destroys** three edge lengths and **creates** three:

- destroyed: source edges `a = l(p\!\to\!u)`, `c = l(u\!\to\!w)`, target edge `R = l(a_t\!\to\!b)`;
- created: merged `M = a + c`, and the split `x = \tau R`, `y = (1-\tau) R`, with
  `\tau \sim U(0,1)`.

Define the reverse free variable `\tau' = a/(a+c) = a/M`, so that the reverse draw
`\tau'` (uniform) reproduces `a = \tau' M`, `c = (1-\tau')M`, and the reverse merge of
`x,y` reproduces `R = x+y`. The map
$$\phi:(a,c,R,\tau)\longmapsto(M,x,y,\tau'),\quad M=a+c,\ x=\tau R,\ y=(1-\tau)R,\ \tau'=\tfrac{a}{a+c}$$
is a bijection on `(0,\infty)^3\times(0,1)`. Its Jacobian determinant is
$$\bigl|\det \partial\phi\bigr| = \frac{R}{a+c} = \frac{R}{M} = \frac{l_{\text{Regraft}}}{l_{\text{Merge}}}.$$
Since the forward and reverse auxiliaries are `U(0,1)` (density 1), the MH/RJ acceptance
contribution from the branch transformation is exactly `|\det\partial\phi| = R/M`, i.e.
`log(lRegraft) − log(lMerge)`. ∎

*Numerical confirmation.* Constructing `∂φ` explicitly and evaluating
`|det ∂φ|` against `R/(a+c)` over random `(a,c,R,τ)` gives ratio `1.0000` to 4 d.p. in
every trial (six trials shown in the run log). Independently, the forward+reverse
round-trip `[log(R)−log(M)] + [log(M)−log(R)] = 0` exactly.

This Jacobian is **rooting-agnostic**: it depends only on the three local lengths, not on
where the root is. So it is correct whether or not the move relocates a subtree across
the global root.

**`treeLength` is preserved, so the absolute-coordinate Jacobian equals the
relative-coordinate Hastings ratio the target requires.** The move conserves the sum of
the touched lengths: `a + c + R = M + x + y` (since `M=a+c` and `x+y=R`), and no other
edge changes; hence `treeLength = Σ absLen` is invariant and the moves operate within a
single tree-length slice. The state is `(treeLength, relBr)` with `absLen = treeLength ·
relBr`; because `treeLength` is held fixed, the map on `relBr` is the map on `absLen`
rescaled by the constant `treeLength` on both the destroyed and created coordinates, so
the rescaling cancels in the determinant: `|det ∂φ|` computed in absolute lengths equals
the ratio in relative coordinates. (Equivalently: parameterise the fixed-sum simplex face
and include the `τ'` row — the determinant is still `R/M`.) So Lemma B is the correct
Hastings factor for the actual `(treeLength, relBr)` state space, not merely for free
positive reals.

**Multiplicity of discrete paths cancels.** If more than one discrete selection
`(prune, regraft)` (or `(prune, reroot, regraft)` for TBR) maps `R→R'`, the same number
map `R'→R` by the prune-clade/reroot bijection used in the count tests, and each
contributes the *same* per-path density; the multiplicities therefore enter `q(R→R')`
and `q(R'→R)` identically and cancel in the ratio. (On binary trees with the canonical
reorder the realised forward move is uniquely invertible, so this is a formality, but it
pre-empts the question.)

### NNI (`src/tree_moves.cpp:33–113`) — `H = 0`

**Discrete proposal.** NNI selects an internal edge `(u→v)` uniformly among the
`nInternal` edges with both endpoints internal (`tree_moves.cpp:42–46,56`), then selects
one child of `v` (`vChildRows`) and one *other* child of `u` (`uSibRows`)
(`:64–87`) and swaps their parents (`:90–94`). The selection probability is
`1 / (nInternal · |vChildRows| · |uSibRows|)`.

- In a rooted binary tree every internal node has exactly 2 children. For an internal
  edge `(u→v)`: `v` has 2 children, and `u` has exactly 1 child other than `v` — this
  holds **including when `u` is the root**, because the root also has out-degree 2.
  Hence `|vChildRows|·|uSibRows| = 2·1 = 2` for *every* internal edge. *Numerically:
  over 200 random rooted trees (5–14 tips) this product is constant `= 2` across all
  internal edges, and `= 2` restricted to root-adjacent edges (`u==root`).*
- `nInternal` is invariant under NNI (a swap rearranges children but neither creates nor
  destroys an internal–internal edge). *Numerically: `nInternal = 4` before and after
  all 3000 NNI moves on the 6-tip tree.*

The reverse NNI uses the *same* internal edge `(u→v)` (still present) and swaps the same
pair back; its selection probability is again `1/(nInternal·2)`. Therefore
`q^{\text{disc}}(R→R') = q^{\text{disc}}(R'→R)`.

**Branch proposal.** NNI moves whole subtrees with their pendant edge lengths intact
(`:90–97`: only `newParent[cRow]`, `newParent[wRow]` change; `absLen` is carried
unchanged through `preorder_weighted_impl`). No length is split or merged ⇒ identity
Jacobian.

**Verdict cross-check.** `tree_moves.cpp:112` returns `logHastings = 0.0`. Both factors
are 1, so `H_{NNI} = 0` is the exact rooted ratio. *Numerically: `logHastings ∈ {0}`
over 5000 moves; root node 7 / degree 2 preserved; 5 distinct root bipartitions reached
(root-adjacent NNI does explore rootings, matching `dev/ecology/rootprobe.R`).*

> **Root-adjacency note.** Because the root node `nTip+1` satisfies `> nTip`, the edge
> `(root → internal-child)` qualifies as an "internal edge" and is eligible for NNI.
> Such a move swaps a grandchild of the root with the root's other child, **changing the
> root bipartition** — i.e. NNI legitimately samples over rootings. This is correct and
> desired in infer-root mode. (The fix-root option of `ebe-spec.md §7` would exclude
> `u==root` NNI; that exclusion is *not* applied in the default, and need not be, since
> the move is reversible as shown.)

### SPR (`src/proposals.cpp:32–178`) — `H = log(lRegraft/lMerge)`

**Discrete proposal.** SPR picks a prune edge uniformly among `eligiblePrune`
(`parent ≠ root`, `proposals.cpp:42–58`), then a regraft edge uniformly among
`candidates` (not in `v`'s subtree, not adjacent to `u`, `:110–129`). Selection
probability `1 / (nElig · nCand_fwd)`.

- `nElig = nEdge − (edges incident to root as parent) = nEdge − 2`, since the root always
  has out-degree 2. This is **constant** over all of `S`. *Numerically: `nElig = 8 =
  nEdge−2` before and after all 2000 SPR moves on the 6-tip tree.*
- `nCand_fwd`. The filter (`proposals.cpp:113–117`) excludes (a) every edge whose child
  lies in `desc(v)` and (b) the two edges incident to `u` other than the already-excluded
  `(u→v)`. The set in (a) is the `2m−2` *internal* subtree edges **plus the pendant edge
  `(u→v)`** = `2m−1` edges, where `m` is the number of tips in `v`'s subtree; (b) removes
  `(u→w)` and `(p→u)`. Hence
  $$n_{\text{Cand}} = nEdge - (2m-1) - 2 = nEdge - 2m - 1.$$
  *(Check against `/tmp/spr_struct.R`, `nEdge=10`: internal `v` with `m=2` gives
  `10−4−1 = 5` ✓ matching the printed `nCand=5`; a tip with `m=1` gives `10−2−1 = 7` ✓
  matching `nCand=7`.)* This depends only on `m` and `nEdge`, both preserved by the move.

**Discrete symmetry incl. across-root.** After regrafting, re-pruning the same subtree
`v` has the same `nElig = nEdge−2` and the same `nCand = nEdge − 2m − 1` (because the
subtree tip-count `m` is unchanged and `u` again has exactly the two adjacent residual
edges). *Numerically: on
the 6-tip tree, all 52 (prune, regraft) forward moves have `nCand_fwd == nCand_rev` when
the reverse prune is identified by clade content (label-invariant); 0 mismatches.* This
**includes moves that cross the global root**: 63% of SPR moves change the root
bipartition and 9 distinct root bipartitions are reached, and the 52/52 symmetric set
contains those. Hence `q^{\text{disc}}(R→R') = q^{\text{disc}}(R'→R)`.

**Branch proposal.** SPR merges `a = l(p→u)` and `c = l(u→w)` into `M = lMerge`
(`proposals.cpp:141,149–151`), and splits the regraft edge `R = lRegraft` into
`x = τR`, `y = (1−τ)R`, `τ ∼ U(0,1)` (`:132,154,159`). By **Lemma B** the
change-of-variables Jacobian is `R/M`, contributing `log(lRegraft) − log(lMerge)`.

**Verdict cross-check.** `proposals.cpp:172` sets exactly
`logHastings = log(lRegraft) − log(lMerge)`. The discrete factor is 1 (proven), so this
is the complete and correct rooted ratio. *Numerically the forward+reverse Jacobian
round-trip sums to 0.*

> **Across-root regraft (brief's explicit concern).** When the pruned subtree regrafts
> onto an edge on the *other* side of the global root, the rooted topology (and the root
> bipartition) changes. The merge `a+c→M` and split `R→(τR,(1−τ)R)` are unaffected by
> this — they are local. And the discrete count `nCand` is symmetric across-root (the
> 52/52 test covers these). So `log(lRegraft/lMerge)` **remains the full ratio** for
> across-root SPR. No correction is needed.

### TBR (`src/tree_moves.cpp:270–512`) — `H = log(lRegraft/lMerge) + log(lSubEdge/lMergeSub)` — **prime suspect, verified**

TBR = SPR with an extra **subtree re-rooting** between prune and regraft. Three discrete
choices now: prune edge (uniform over `nElig`), subtree re-root edge (uniform over
`nSubEdge` subtree edges, `:345–349`), regraft edge (uniform over `nCand`,
`:457–469`). Selection probability `1 / (nElig · nSubEdge · nCand_fwd)` when `v` is
internal.

**(i) `nElig` and `nCand` symmetry.** Identical argument to SPR: `nElig = nEdge−2`
(constant); `nCand = nEdge − 2m − 1` depends only on the subtree tip-count `m` and
`nEdge`, both preserved by the move. The TBR comment's claim
(`tree_moves.cpp:256,500–503`) that the candidate-count ratio is "identically 0" because
"re-rooting preserves `|desc(v)|` and `u` always has 3 adjacent edges" is **correct as a
rooted statement** — re-rooting the subtree permutes its *internal* edges but changes
neither its tip-set `m` nor the number of edges descending from `v` (`2m−1` including the
pendant `(u→v)`), and `u`'s three incident edges (parent, the two children) are untouched
by the reroot. (The comment's "`nEdge − |desc(v)| − 2`" uses `|desc(v)|` to mean the
descendant *node* count; what matters for the ratio is only that it is a function of `m`,
which it is.)

**(ii) `nSubEdge` symmetry.** `nSubEdge =` number of edges within `v`'s subtree `= 2m−2`
(`m` = subtree tips). Re-rooting preserves `m`, hence preserves `nSubEdge`.
*Numerically: across 40 random 12-tip trees, **804** (prune, subRow, regraft) general-case
moves (interior path-reversal lengths 1–5, so the `Transform 3` reversal loop
`tree_moves.cpp:412–425` is exercised) all satisfy `nSubEdge_fwd == nSubEdge_rev` **and**
`nCand_fwd == nCand_rev`; 0 mismatches.*

**(iii) Reverse reroot is reachable.** The reverse move must re-root the subtree back to
the original position by choosing the subtree edge that re-merges to the original split.
*Numerically: across 30 random 11-tip trees (510 general-case moves) that specific
reverse-reroot edge is always a member of the proposed tree's selectable subtree-edge set
(`subEdgeRows`); 510/510.* So the reverse selection has the same `1/nSubEdge`
probability and is reachable. Hence `q^{\text{disc}}(R→R') = q^{\text{disc}}(R'→R)`.

**(iv) Branch Jacobian (two factors).** The regraft contributes `R/M` (Lemma B). The
subtree-reroot has the **same local structure** as Lemma B: it merges `v`'s two child
edges `p,q` into `lMergeSub = p+q` (carried onto the path child `a₁`,
`tree_moves.cpp:398–399`) and splits the chosen subtree edge `S = lSubEdge` into
`σS` and `(1−σ)S`, `σ ∼ U(0,1)` (`:394,404,430`). By Lemma B applied to
`(p, q, S, σ)`, its Jacobian is `S/(p+q) = lSubEdge/lMergeSub`. The two transformations
act on disjoint edge sets, so the combined Jacobian is the product, hence the **sum** in
log-space. *Numerically: `[log(R)−log(M)] + [log(S)−log(p+q)]` forward and its reverse
sum to exactly 0.*

**(v) The `x == v` no-op branch (`tree_moves.cpp:432–447`).** When the chosen subtree
edge is directly below `v` (one of `v`'s two children), the reroot is the identity:
the code draws `σ` and **discards** it (`(void)rng.unif()`), makes no topology/length
change, and sets `lMergeSub = lSubEdge` so the sub-Jacobian term is `log(S)−log(S)=0`.
The discarded `σ` is `U(0,1)` in both directions and cancels. This branch reduces TBR to
plain SPR and is reversible by the SPR argument; the `nSubEdge` count still includes
these two edges symmetrically in forward and reverse.

**Verdict cross-check.** `tree_moves.cpp:504–506` returns exactly
`log(lRegraft) − log(lMerge) + log(lSubEdge) − log(lMergeSub)`. All discrete factors are
1 (proven) and both Jacobians are accounted for. **The TBR rooted `q`-ratio is correct.**
The brief's worry — that the explicit subtree re-rooting might make the candidate-count
ratio nonzero in the rooted space — is resolved: re-rooting permutes subtree-internal
edges but preserves `m`, hence preserves every count.

> **Why the old (root-invariant) tests could not catch a TBR error here.** Under
> root-invariant Mk, all rerootings of one unrooted topology have equal likelihood, so
> the chain's stationary distribution over the *unrooted* topology is insensitive to how
> proposal mass is split among rootings. A wrong `nSubEdge`/`nCand` ratio would have
> biased only the (invisible) distribution over rootings. Under EBE that distribution is
> the inferential target. The verification above is precisely the previously-untested
> part, and it passes.

### pSPR (`src/mcmc.cpp:3513–3691`) — `H = log(lRegraft/lMerge) + log(wOrig) − log(w_chosen) + log(sumW) − log(sumWRev)`

pSPR is SPR whose regraft target is drawn **non-uniformly**, weighted by a parsimony
bias `w_i = exp(−α·(score_i − score_min))`, `α = 0.1` (`mcmc.cpp:3511,3616–3635`). It
does **not** re-root any subtree (the subtree keeps its internal rooting); it is a guided
SPR.

**Discrete structure.** Prune edge uniform over `nElig` (`:3530–3539`). Regraft edge
drawn with probability `w_chosen / sumW` over the candidate set `F` (forward candidates
exclude the original position; `:3572–3635`). The forward discrete proposal is
`(1/nElig)·(w_chosen/sumW)`. The reverse: prune the same subtree (same residual tree),
draw the *original* position with probability `w_orig / sumW_rev`, where the reverse
candidate set `Rev` excludes the chosen position and so has normaliser
`sumW_rev = sumW + w_orig − w_chosen` (`:3679–3682`).

**Three facts make this a correct rooted ratio:**

1. **The residual tree is identical forward and reverse.** Pruning `v` yields the same
   residual tree regardless of where `v` was attached, so `F` and `Rev` differ only by
   the swapped excluded edge, and all candidate scores/weights are computed on the same
   residual tree.
2. **The SPR candidate set is symmetric over rooted trees** (proven in the SPR section:
   `nElig` constant, `nCand` symmetric incl. across-root). pSPR uses the *same* candidate
   filter (`mcmc.cpp:3572–3578` ≡ `proposals.cpp:113–117`), so the candidate *sets*
   correspond one-to-one forward/reverse.
3. **Fitch parsimony is root-invariant.** The weights derive from `fitch_score_all` /
   `fitch_score_candidates` (`mcmc.cpp:3604,3613`), standard unrooted Fitch (Fitch 1971;
   Felsenstein 2004, ch. 1): the minimum number of state changes on the *unrooted* tree;
   the root adds a length-0 fictitious split contributing no change. *Numerically
   confirmed via an independent oracle (`phangorn::parsimony`): score `= 28` for the
   original tree and for re-rootings on `t1`, `t3`, `t5` — invariant.* Hence
   `score_i`, `w_i`, `sumW`, `w_orig`, `sumW_rev` are all functions of the rooted tree
   that take the same value under any rooting; they are well-defined and match the
   forward/reverse bookkeeping irrespective of where the root sits, including when the
   move crosses the root.

**Branch Jacobian.** Same merge/split as SPR ⇒ `R/M` (Lemma B), term
`log(lRegraft) − log(lMerge)` (`mcmc.cpp:3680`).

**Assembling the ratio.**
$$
\frac{q(R'\to R)}{q(R\to R')}
= \underbrace{\frac{R}{M}}_{\text{branch (Lemma B)}}
\cdot
\underbrace{\frac{(1/nElig)\,(w_{\text{orig}}/sumW_{\text{rev}})}{(1/nElig)\,(w_{\text{chosen}}/sumW)}}_{\text{discrete}}
= \frac{R}{M}\cdot\frac{w_{\text{orig}}}{w_{\text{chosen}}}\cdot\frac{sumW}{sumW_{\text{rev}}},
$$
whose log is exactly `mcmc.cpp:3680–3685`. **The pSPR rooted `q`-ratio is correct**,
because every parsimony quantity it uses is root-invariant and the SPR substrate is
root-symmetric.

> **Caveat (1).** pSPR's correctness as written is *inherited* from (a) the SPR
> discrete-symmetry result above and (b) Fitch root-invariance. Both are established
> here. No residual gap. (This is the one move whose proposal reads the data, but only
> through a root-invariant statistic, so Assumption 5's spirit is preserved.)

---

## Implementation cross-check (per theorem step)

| Step | Claim | Source |
|------|-------|--------|
| Lemma A | reorder fixes root, identity Jacobian on lengths | `renumber_tree.h:204,273,309,158` |
| Active set | NNI/SPR always on; TBR/pSPR default-on; gibbs_* eco-skip; weighted_* default-off | `R/MkPrimeMCMC.R:343–352`, `R/RunMkPrime.R:3482–3539`, `mcmc.cpp:1024,1791` |
| SPR eco uses full-eval | partial-CL path gated `!ecologyAware` | `mcmc.cpp:5049` |
| NNI `H=0` | constant branching factor 2, `nInternal` invariant, lengths carried | `tree_moves.cpp:42–94,112` |
| SPR `H=log(R/M)` | `nElig=nEdge−2`, `nCand` symmetric, merge/split | `proposals.cpp:42–58,110–129,141,154,159,172` |
| TBR `H=log(R/M)+log(S/(p+q))` | + `nSubEdge=2m−2` symmetric, reverse reroot reachable, 2 Jacobians | `tree_moves.cpp:256,345–349,394–430,457–469,500–506` |
| pSPR `H=...+log(wOrig/wChosen)+log(sumW/sumWRev)` | residual tree identical; Fitch root-invariant; SPR substrate symmetric | `mcmc.cpp:3530–3539,3572–3635,3680–3685,3613` |
| Lemma B | `|det ∂φ| = R/M`; sub-version `S/(p+q)` | derived; numeric `ratio=1.0000` |

No divergence found between the derivations and the code for any of the four
default-active moves.

---

## Edge cases

1. **`nTip = 5` (smallest tree with NNI/SPR active, `nEdge=8`).** `nElig = 6`,
   `nInternal ≥ 1`. All counts well-defined; arguments hold. (Registration guard is
   `nEdge ≥ 5`, `RunMkPrime.R:3482`; satisfied for `nTip ≥ 4`.)
2. **Pruned subtree is a single tip (`v ≤ nTip`).** SPR: `|descEdges(v)|=0`, fine. TBR:
   the reroot is skipped (`lSubEdge=lMergeSub=1`, term 0, `tree_moves.cpp:448–452`), and
   TBR degenerates to SPR — reversible by the SPR argument. *(A tip↔tip TBR is a pure
   SPR; the reverse is also tip-pruned SPR.)*
3. **`u`'s parent `p` is the root.** SPR/TBR still well-defined: `parentRow` finds the
   `(p→u)` edge (`p=root` is fine, the guard is on `parent[pruneRow]≠root`, i.e. `u≠root`,
   not `p≠root`). Suppressing `u` makes `p→w`; `p` remains the (degree-2) root because it
   keeps two children (`w` and its other child). *Confirmed by Lemma A's numeric root
   check, which exercised many such configurations.*
4. **Across-root regraft / reroot.** Covered explicitly in SPR and TBR sections; 63%
   (SPR) of moves change the root bipartition and remain reversible.
5. **`x == v` reroot in TBR.** No-op branch, σ discarded symmetrically; covered above.
6. **Degenerate proposal (`candidates` or `internalRows` empty).** Each impl returns
   `logHastings = R_NegInf` (`proposals.cpp:50,86,121`; `tree_moves.cpp:51,75`;
   `tree_moves.cpp:291,306,464`), and the dispatcher rejects (`if(!R_FINITE(logHastings))
   return false`, `mcmc.cpp:5062,5075,5167`). A `−∞` proposal is never accepted, so it
   does not perturb reversibility. Such configurations do not arise for `nTip ≥ 4`
   binary trees in practice.
7. **`σ` or `τ` numerically at `{0,1}`.** `U(0,1)` from `ChainRng::unif` returns values
   in `[0,1)`; a `0` split would zero a branch length. Not a Hastings-correctness issue
   (the Jacobian formula is exact for any `τ∈(0,1)`); it is a numerical-robustness
   concern for the *likelihood* at a zero-length edge. Tagged below.

---

## Non-default eco-capable moves (assessed, not in critical path)

`weighted_branch_scale` (12), `weighted_spr` (13, `mcmc.cpp:2692–2955`),
`weighted_subtree_swap` (14, `mcmc.cpp:2973–...`) are **off by default**
(`MkPrimeMCMC.R:349–351`) and so do not bear on the infer-root validity *as shipped*.
For completeness, should a user enable them in eco mode:

- They are **internally self-contained MH moves**: each computes `newLogLik` via
  `compute_full_loglik_at` (which, in eco mode, is the full EBE likelihood), computes the
  proposed prior, forms `logAlpha = β·(newLogLik − logLik) + (newLogPrior − logPrior) +
  logHR`, and accepts/rejects internally (`mcmc.cpp:2937–2954`), then `return`s a bool
  from the dispatcher (`mcmc.cpp:5141,5144`). They therefore do **not** double-apply or
  bypass the EBE likelihood — unlike the eco-disabled gibbs_spr/gibbs_subtree_swap, which
  used a *blind* candidate evaluator. This is why they correctly lack an eco self-disable.
- Their `logHR` is the **branch-fraction component only** (`mcmc.cpp:2919–2922`); the
  topology-selection normaliser is claimed to cancel (`Z = Z'` "by symmetry of the
  candidate set", `mcmc.cpp:2685–2686`). The candidate **set** symmetry is the SPR result
  proven here, so the sets correspond. **However**, the claim `Z = Z'` (equal forward and
  reverse marginal-weight normalisers `sumM`) is a *stronger* statement than set symmetry:
  the per-candidate marginal weights are full eco likelihoods, and `sumM` (forward,
  centred on `R`) need not equal `sumM'` (reverse, centred on `R'`) in general for a
  guided move. The same potential gap exists in the blind `gibbs_spr` family. Because
  these moves are disabled in eco mode by default, I do **not** close this here.

> **Caveat (2).** The `Z = Z'` topology-cancellation assumption of `weighted_spr` /
> `weighted_subtree_swap` (and structurally `gibbs_spr` / `gibbs_subtree_swap`) is not
> proven for a likelihood-guided proposal over rooted trees. It is moot by default
> (these moves carry zero weight in eco mode), but **must be verified before any eco run
> enables `weightedSpr`/`weightedSubtreeSwap`**. Closing it requires either a detailed-
> balance proof that the marginal normaliser cancels for guided moves, or a missing
> `+ log(sumM) − log(sumM')` correction term. Should be picked up by lane L? / role
> **math-prover** (analytic) with a **mcmc-diagnostician** detailed-balance simulation as
> the empirical gate.

---

## Caveats / downstream tags

> **Caveat (3).** This proof establishes **per-move reversibility** (the spec's
> "per-move detailed balance" requirement) analytically + numerically. The spec also asks
> for an **exact-rooted-posterior** confirmation: run the eco-mode kernel on a tiny tree
> against an independently enumerable rooted target and check the sampled rooted-tree
> distribution matches. That is an empirical end-to-end check, not an analytic step, and
> is **out-of-lane** for math-prover. Should be picked up by role
> **mcmc-diagnostician** (e.g. a 4–5-tip exhaustive rooted-tree target with `z≠0,
> φ≠1`, comparing MCMC frequencies to the normalised `π`). My per-move results predict it
> will pass; it is the right belt-and-braces check before trusting published root
> posteriors.

> **Caveat (4).** Zero-length-edge robustness (edge case 7) under `τ,σ→0`: the Hastings
> Jacobian is exact, but a near-zero `lRegraft·τ` edge can produce `ex = exp(−2·tEff)→1`
> in `ebe_mkn_P` and a (finite) but extreme likelihood. This is a *numerical* concern,
> not a Hastings-correctness one. Should be picked up by role **numerical-auditor**.

---

## Verdict

Per the brief's required **per-move verdict** (rooted-`q`-ratio CORRECT, or BIASED with
correction):

| Move | Default-active in eco? | Stated `logHastings` | Rooted-`q`-ratio verdict |
|------|:----------------------:|----------------------|--------------------------|
| **NNI** | yes | `0` | **CORRECT** — symmetric discrete selection (branching factor 2 incl. root-adjacent), lengths carried, identity Jacobian. |
| **SPR** | yes | `log(lRegraft/lMerge)` | **CORRECT** — `nElig=nEdge−2` constant, `nCand` symmetric incl. across-root (52/52), Lemma B Jacobian `=R/M`. |
| **TBR** | yes | `log(lRegraft/lMerge)+log(lSubEdge/lMergeSub)` | **CORRECT** — the "candidate-count ratio identically 0" claim holds *as a rooted statement*; `nSubEdge=2m−2` and `nCand` both reroot-invariant (804/804), reverse reroot reachable (510/510), two independent Lemma-B Jacobians. The prime suspect is cleared. |
| **pSPR** | yes | `log(lRegraft/lMerge)+log(wOrig/wChosen)+log(sumW/sumWRev)` | **CORRECT** — guided SPR on a root-symmetric candidate set with **root-invariant Fitch** weights (verified vs `phangorn`). |

**Overall: Watertight (with caveats).** For the four moves active by default in ecology
mode, each stated `logHastings` is the exact `log[q(R'→R)/q(R→R')]` over rooted
trees-with-branch-lengths, including all branch Jacobians and discrete selection
probabilities. Combined with the (separately proven) root-invariant prior and the
EBE likelihood applied by the outer MH, the default ecology kernel is reversible w.r.t.
the rooted EBE posterior. **Infer-root is valid as-is for the default move set; no move
needs to be fixed or disabled.**

The caveats do not affect this conclusion: Caveat (2) concerns only the **non-default**
`weighted_*` moves (zero weight in eco mode — must be re-checked if ever enabled);
Caveats (3)–(4) are empirical/numerical belt-and-braces checks for other lanes that my
analysis predicts will pass.

**No patch attached** — no implementation bug found; nothing to fix in the default path.

### Assumptions not enforced by the code (for the record)
- Assumption 3 (root preservation by reorder) is a property of `TreeTools`, not asserted
  in `src/`; it is verified here numerically and by reading `renumber_tree.h`. If a future
  `TreeTools` upgrade changed the root-selection rule, every verdict would need re-checking.
- Assumption 6 (`U(0,1)` auxiliaries) holds for `ChainRng::unif`; if a move were switched
  to a non-uniform split proposal, the Lemma-B Jacobian would acquire the auxiliary's
  density ratio.

---

### References
- Fitch, W. M. (1971). Toward defining the course of evolution. *Syst. Zool.* 20:406–416. (Fitch parsimony, root-invariance.)
- Felsenstein, J. (2004). *Inferring Phylogenies*, ch. 1 (parsimony) and ch. 16 (MCMC). Sinauer.
- Lakner, C., van der Mark, P., Huelsenbeck, J. P., Larget, B., Ronquist, F. (2008). Efficiency of Markov chain Monte Carlo tree proposals in Bayesian phylogenetics. *Syst. Biol.* 57:86–99. (Stochastic NNI/SPR/TBR branch-split proposals and their Hastings/Jacobian factors.)
- Yang, Z. (2014). *Molecular Evolution: A Statistical Approach*, ch. 7–8. OUP.
- `dev/ecology/ebe-spec.md §7` (rooting decision); `dev/ecology/prior-reroot-check.R` (prior root-invariance, upstream).
