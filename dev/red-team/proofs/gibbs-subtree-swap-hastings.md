# GSWAP-001 — the Hastings ratio of `gibbs_subtree_swap`

Checked against `src/mcmc.cpp` at `66510e1` (`origin/main`, 2026-09-19).
Resolves issue #21.

**Verdict: the move as committed at `66510e1` is not π-invariant.** The
correction is an MH accept step with

$$\alpha \;=\; \min\!\Big(1,\ \frac{Z_x(A)}{Z_y(A)}\Big)$$

where $Z$ is the selection normaliser of the anchored swap neighbourhood.

## 1. The kernel

`gibbs_subtree_swap_impl` (and its `_het` / `_full` siblings, which differ
only in how likelihoods are evaluated):

1. Draw an anchor uniformly from the `nEdge` edge-children, i.e. uniformly
   over the non-root nodes (`src/mcmc.cpp`, "Pick a random node").
2. `get_valid_swap_partners_impl` (`src/tree_moves.cpp:178-226`) returns every
   edge-child $w$ that is not a descendant of $A$ (nor $A$ itself), not an
   ancestor of $A$, and not a sibling of $A$. The relation is symmetric.
3. Evaluate $\mathrm{LL}(S(x,A,B'))$ for every $B' \in P_x(A)$, where
   $S(x,A,B)$ transposes the parent assignments and the stem lengths of the
   two subtrees.
4. Draw from $\{\text{stay}\} \cup P_x(A)$ with weights
   $\exp(\beta\,\mathrm{LL})$, normaliser
   $Z_x(A) = e^{\beta \mathrm{LL}(x)} + \sum_{B'} e^{\beta \mathrm{LL}(S(x,A,B'))}$.
5. Commit the draw. **There is no accept/reject step.**

## 2. Detailed balance

Condition on the anchor. The swap leaves the edge set intact — it only
re-attaches two subtrees — so there is a canonical bijection between the
anchor choices at $x$ and at $y = S(x,A,B)$, fixing $A$ and $B$; and `nEdge`
is invariant, so the anchor probability is $1/nEdge$ in both directions. It
therefore suffices to balance matched anchors term by term. (The same $y$ is
also reachable by anchoring on $B$, since $S(x,B,A) = S(x,A,B)$; that branch
balances separately, with its own ratio $Z_x(B)/Z_y(B)$. A sum of reversible
kernels is reversible.)

$S(\cdot,A,B)$ is an involution, and $B$ is the unique element of $P_y(A)$
with $S(y,A,B) = x$ — up to the π-null set where two swapped subtrees are
isomorphic *with identical branch lengths*. So

$$q_A(x \to y) = \frac{e^{\beta \mathrm{LL}(y)}}{Z_x(A)}, \qquad
  q_A(y \to x) = \frac{e^{\beta \mathrm{LL}(x)}}{Z_y(A)} .$$

The swap transposes two entries of `relBrLengths` and leaves `treeLength`
alone: a permutation, so the Jacobian is exactly 1, and the flat
`Dirichlet(1, ..., 1)` branch prior (`src/mcmc.cpp`, `cpp_log_prior`) is
exchangeable over edges, so the prior ratio is 1 as well. With no topology
prior, $\pi \propto e^{\beta \mathrm{LL}}$ on the candidate set, and

$$\frac{\pi(y)\,q_A(y \to x)}{\pi(x)\,q_A(x \to y)}
 = \frac{e^{\beta \mathrm{LL}(y)} e^{\beta \mathrm{LL}(x)} / Z_y(A)}
        {e^{\beta \mathrm{LL}(x)} e^{\beta \mathrm{LL}(y)} / Z_x(A)}
 = \frac{Z_x(A)}{Z_y(A)} .$$

The implemented kernel uses $\alpha = 1$, which is correct only if
$Z_x(A) = Z_y(A)$.

> **Guard.** The prior ratio is 1 *because the branch prior is exchangeable*.
> Lengths stay with their slots, not their subtrees, so a length can move
> between an internal and a pendant edge. A compound-Dirichlet prior with
> distinct internal and pendant concentrations would give this move a
> prior-ratio term it does not currently carry.

## 3. Why the normalisers differ

`gibbs_spr` gets away without this term because both endpoints of an SPR
prune to the *same residual tree*, so a selection distribution that is a
function of that tree alone has the same normaliser in both directions. A
swap has no such shared object.

Write $d(A)$ for the number of edges from the storage root to $A$, and
$\mathrm{desc}^*(A)$ for $A$ together with its descendants. The tree is
stored unrooted, with node $nTip+1$ a **trifurcation**, so a child of the
root has *two* siblings and `nEdge = 2 nTip - 3`. Then

$$|P_x(A)| = nEdge - |\mathrm{desc}^*(A)| - g(d(A)), \qquad
  g(1) = 2, \quad g(d) = d \ \ (d \ge 2),$$

verified exhaustively against `get_valid_swap_partners_impl` over every
rooted representation of every unrooted binary topology for
$nTip = 4, \dots, 8$. Since $g(1) = g(2)$, sizes agree iff
$d(A) = d(B)$ or $\{d(A), d(B)\} = \{1, 2\}$. Sizes first differ at
$nTip = 5$; at $nTip = 6$ they differ for **57.4%** of valid ordered pairs.
Issue #21's "at six taxa the forward and reverse neighbourhoods are nearly
always the same size" is therefore wrong, and so is the opposite reading —
size equality does not rescue the move, because the *sets* differ too.

The set argument is the load-bearing one, and it bites at $nTip = 4$, where
every neighbourhood has the same size. With root $5\{t_1, t_2, m\}$,
$m\{t_3, t_4\}$, stem lengths $l_1 \dots l_5$, anchor $A = t_1$, partner
$B = t_3$:

- $P_x(t_1) = \{t_3, t_4\}$; the third member of the neighbourhood,
  $S(x, t_1, t_4)$, is the split $\{t_1,t_3\}|\{t_2,t_4\}$ with pendant
  lengths $t_1 = l_5,\ t_2 = l_2,\ t_3 = l_4,\ t_4 = l_1$.
- $P_y(t_1) = \{t_3, t_2\}$; its third member, $S(y, t_1, t_2)$, is the *same*
  split with pendant lengths $t_1 = l_2,\ t_2 = l_4,\ t_3 = l_1,\ t_4 = l_5$.

Same topology, different length assignment, hence different likelihood, hence
$Z_x(t_1) \neq Z_y(t_1)$. The neighbourhood is a *root-dependent* object
attached to a *root-independent* target.

## 4. Verification

**Exact, $nTip = 4$.** The finite kernel was built state by state from the
source (24 reachable states: rooted representation plus a fixed multiset of
five distinct stem lengths; `LL` an arbitrary root-invariant function).

| kernel | $TV(\pi P, \pi)$ | max $|\pi P - \pi| / \pi$ | max DB residual |
|---|---|---|---|
| as implemented | 1.86e-02 | 2.16e-01 | 2.17e-03 |
| with $\alpha = \min(1, Z_x/Z_y)$ | 3.90e-17 | 2.94e-16 | 5.20e-18 |

Individual states are up to **48%** off target.

**Empirical, on the compiled kernel, $nTip = 5$.** Equal stem lengths make
the reachable state space finite (45 states), so the target
$\pi(s) \propto e^{\mathrm{LL}(s)}$ is exact. The defect over-visits
high-probability states, i.e. gives a positive slope of
$\log(\text{empirical}/\text{target})$ on $\log(\text{target})$:

| iterations | | uncorrected | corrected |
|---|---|---|---|
| 60 000  | TV | 0.0293 | 0.0163 |
|         | slope (t) | +0.1006 (**+7.8**) | +0.0163 (+1.5) |
| 300 000 | TV | 0.0296 | 0.0065 |
|         | slope (t) | +0.0922 (**+10.5**) | +0.0014 (+0.35) |

The uncorrected TV does not shrink with the run length; the corrected one
falls as $1/\sqrt{N}$ and the slope is indistinguishable from zero. Gate:
`tests/testthat/test-gibbs-subtree-swap-invariance.R`.

## 4a. Why the earlier evidence missed this

`dev/red-team/findings-archive.md` carries `SWAP-001` as **REFUTED**, severity
downgraded HIGH → LOW, on `dev/red-team/heavy-tests/subtree-swap-db.R`
(Hamilton job 17304185, n = 6). That row is superseded. Three reasons the
harness could not have seen this, all visible in the harness itself:

1. **It never ran the compiled move.** `subtree-swap-db.R` is a pure R driver
   that reproduces the draw using `get_valid_swap_partners_cpp` and
   `swap_subtrees_cpp`, so it tests the analytic partner-count claim, not
   `gibbs_subtree_swap_impl`. Nothing in the commit path — including the
   branch-length convention defect in §6 — was reachable.
2. **It ran at β = 0 only**, where every candidate weight is 1 and `Z` reduces
   to `1 + |P|`. The set asymmetry of §3, which bites even where the sizes
   agree, is structurally invisible there. This is the trap
   `dev/red-team/proofs/hastings-tree-moves.md` §5 names.
3. **Its own verdict says "not detected", not "refuted"**:
   `subtree-swap-db.R` prints `SWAP-001/SWAP-002 not detected.` The archive's
   prose row says as much too — "the partner-set log-ratio is genuinely absent
   from the move" — before downgrading anyway.

`dev/red-team/findings-archive.md` is frozen, so the stale row stays; this
section is where a reader who starts there should end up.

## 4b. Deviation from the filed gate

`dev/plans/2026-08-12-gibbs-spr-fix-design.md` §7 and issue #21 both ask for an
arm in `subtree-swap-db.R` at "a taxon count large enough for the forward and
reverse neighbourhoods to differ". No arm was added. The premise behind that
request is false — the kernel already fails at n = 4, where every neighbourhood
has the same size, so taxon count is the wrong dial — and the harness's β = 0,
pure-R design cannot reach the defect at any n. The gate is instead the exact
n = 4 kernel and the finite-state n = 5 chain above, both of which drive the
compiled move at β = 1.

## 5. Cost

The reverse normaliser needs $|P_y(A)|$ extra likelihood evaluations on a
tree the forward downpass did not cache, so the move costs roughly twice what
it did. There is no cheaper correct variant: unlike `gibbs_spr`, no
reparameterisation makes $Z_x = Z_y$, because the neighbourhood is not an
equivalence class.

## 6. A second defect, found in passing

`gibbs_subtree_swap_impl_het` committed the *opposite* branch-length
convention from the one it scored. `evaluate_swap_impl`
(`src/gibbs_partial_cl.h:844+`) attaches `nodeB` at `pA` with `lenA` and
`nodeA` at `pB` with `lenB` — each subtree takes the stem of the slot it
moves *into*. The main and full paths commit that (transposing parents and
lengths together); the Q-heterogeneity path transposed `child` and gave each
subtree the stem it arrived *with*. `state->logLik` was then the likelihood
of a tree other than the one in `state`, drifting by up to **6.5 log-units**
per accepted swap and corrupting every downstream MH ratio. Unifying the
three paths behind one commit removes it by construction; gated by the second
test in the file above.
