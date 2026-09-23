# WSWAP-001 — the Hastings ratio of `weighted_subtree_swap`

Derived against `src/mcmc.cpp` at `c6b9fc9` (`fix/tree-move-hastings`, before
this fix), 2026-09-21. Fixed on `fix/weighted-move-hastings`. Part of issue
#158; confirms the archived `SWAP-002`, and corrects its fix sketch.

**Verdict: three terms were wrong, and none of them cancels.** The ratio needs
the bin mixture (as in `weighted_spr`), scored with the *reverse* topology's
weights in the reverse direction, and the ratio `Z_x(A) / Z_y(A)` of the
anchored neighbourhoods' selection normalisers (as in `gibbs_subtree_swap`).
The reverse neighbourhood has to be enumerated in full, which roughly doubles
the cost of the move.

The move is **default-off** (`weightedSubtreeSwap = FALSE` in
`MkPrimeMCMC()`).

## 1. The kernel

`weighted_subtree_swap_impl` (moveType 14), `BranchBins` as in
`weighted-spr-hastings.md` §1. Write `l_A`, `l_B` for the absolute stem lengths
of the subtrees rooted at `A` and `B`, and `S(x, A, B, f)` for `x` with the
parent assignments of `A` and `B` transposed, `l_A = f t` and
`l_B = (1 - f) t`, where `t = l_A + l_B`.

1. Anchor on `A`, uniform over the `nEdge` non-root nodes.
2. `P_x(A)`: every node that is not `A`, a descendant, an ancestor or a sibling
   of `A` (`get_valid_swap_partners_impl`).
3. Weights: `w^x_{B,b} = exp(beta LL(S(x, A, B, mids[b])))` for each
   `B ∈ P_x(A)` and bin `b`, and `w^x_0 = exp(beta LL(x))` for staying put;
   `Z_x(A) = w^x_0 + Σ_{B,b} w^x_{B,b}`.
4. Draw `(B, b')` with probability `w^x_{B,b'} / Z_x(A)` (staying put is a
   no-op), then `f' ~ Beta(mids[b'] conc + 1, (1 - mids[b']) conc + 1)`.
5. Propose `y = S(x, A, B, f')`; accept by MH, including the prior ratio.

## 2. The reverse move

Pair the anchors as in `gibbs-subtree-swap-hastings.md` §2: the swap keeps
the edge set, so anchoring on `A` at `y` is matched with anchoring on `A` at
`x`, with probability `1 / nEdge` both ways. (Anchoring on `B` reaches the
same topology through a different fraction, `1 - f'`, and balances
separately. A sum of reversible kernels is reversible.)

`B ∈ P_y(A)`: the swap rewrites only the parent edges of `A` and `B`, and
neither lies on the other's path to the root, so in `y` `B` is still neither a
descendant, an ancestor nor a sibling of `A`. And `S(y, A, B, f) = x` exactly
when `f = l_A / t`, where `t` is unchanged by the move. So

    q_A(x -> y) = (1/nEdge) · Σ_b w^x_{B,b} Beta(f' | b) / Z_x(A)
    q_A(y -> x) = (1/nEdge) · Σ_b w^y_{B,b} Beta(f  | b) / Z_y(A),   f = l_A / t

As in `weighted_branch_scale`, the bin is integrated out, and parameterising
the pair `(l_A, l_B)` by `(f, t)` makes the dimension-matching map
`(f, u = f') -> (f', u' = f)` a transposition with `|J| = 1`. Hence

    log HR = log Σ_b w^y_{B,b} Beta(f | b)  - log Z_y(A)
           - log Σ_b w^x_{B,b} Beta(f' | b) + log Z_x(A)

## 3. Why nothing cancels

**The normalisers.** `Z_x(A) ≠ Z_y(A)` for the reasons of
`gibbs-subtree-swap-hastings.md` §3: the anchored partner set depends on the
depth of `A` below the storage root, and the configurations its partners reach
depend on where `A` sits. That section's `nTip = 4` counter-example, where every
neighbourhood has the same size but different members, applies unchanged.

**The bin weights.** `w^x_{B,b}` scores `y`'s topology at each bin midpoint;
`w^y_{B,b}` scores `x`'s. They are different functions of `b`, so unlike in
`weighted_branch_scale` the two mixtures do not share weights.

**What the code had.** Before this change,

    log HR = log w^x_0 + log Beta(f_d | bin(f_d)) - log w^x_{B,b'} - log Beta(f' | b')

on the forward normaliser's offset, with `f_d = absLen[rowB] / t`. Three
things are wrong with that:

1. `w^x_0`, the forward *stay* weight, stands in for the reverse *swap*
   weight `w^y_{B,·} / Z_y(A)`; no reverse neighbourhood is ever evaluated.
2. A single component stands in for each mixture.
3. `f_d` was read *after* `absLen[rowB]` had been overwritten with
   `(1 - f') t`, so `f_d = 1 - f'`. The bins are symmetric about ½, so the
   "reverse" component is just the containing-bin component *at `f'`*. The
   old ratio never looked at the current fraction at all.

**The archived fix sketch is wrong too.** `SWAP-002` proposed adding
`log |P_y(B)| - log |P_x(A)|`. The matched reverse anchor is `A`, not `B`; and
the ratio is of normalisers, not partner counts. The two coincide only at
`beta = 0`, and even there only approximately: `Z = 1 + nBins |P|`.

## 4. An implementation trap: internal labels

`preorder_weighted_impl` renumbers internal nodes. If `A` or `B` is internal,
the node that carries `A`'s label in the *canonicalised* `y` is a different
node, and a reverse neighbourhood enumerated there answers the wrong question.
A first draft of this fix did exactly that. It passed every length statistic at
`beta = 0` but spent 18.2% of its time on three-cherry topologies against a
target of 14.3%. The reverse neighbourhood is therefore enumerated on the
un-reordered proposal; canonicalisation waits for the commit, as in
`gibbs_subtree_swap`. The replicate gate of §5 does not see this draft; a
single 40 000-iteration chain does (z = 5.9 on the three-cherry frequency), and
is gated.

## 5. Verification

**`beta = 0`, `pi K^n = pi`.** `tests/testthat/helper-prior-invariance.R`,
`nTip = 6`: 1500 replicates of 25 moves, each starting from an exact prior
draw. The storage root is part of what the move sees, so each draw puts the
trifurcation at a uniformly chosen internal node.

| statistic | uncorrected | corrected | target | KS p uncorrected | KS p corrected |
|---|---|---|---|---|---|
| smallest fraction | 0.00422 | 0.01253 | 0.01233 | **5.9e-212** | 0.42 |
| largest fraction  | 0.4026  | 0.3150  | 0.3143  | **1.0e-142** | 0.57 |

(The corrected run used 3000 replicates.) The uncorrected smallest fraction
lands 66% below its target mean. A 60 000-iteration chain of the corrected
move spends 14.1% of its time on three-cherry topologies (target 1/7) and has
`t1` on the storage root 25.3% of the time (target ¼).

**The normaliser term alone is almost invisible at `beta = 0`.** A build that
has the bin mixture but drops `Z_x / Z_y` passes the `beta = 0` gate. On the
number of valid swap pairs in the storage representation, the statistic that
term reweights, it gives `t = -0.6` at `nTip = 6` and `t = 2.0` at
`nTip = 8`. This is why `SWAP-002` was archived as REFUTED: its harness
(`subtree-swap-db.R`) ran at `beta = 0` and looked only at topology
frequencies. It never saw the bin-mixture error either, which is glaring in
the lengths (KS p = 6e-212 above).

**`beta = 1`, exact topology posterior.** Five tips give fifteen unrooted
topologies. Each one's posterior mass is integrated over the Dirichlet edge
fractions by plain Monte Carlo. The likelihood is invariant to the storage
root, so that mass is the chain's target however the root wanders. A
20 000-iteration chain of moveType 14 alone, with batch-means z-scores over
50 batches:

| build | TV to target | max \|z\| |
|---|---|---|
| uncorrected | 0.21 | 24.8 |
| bin mixture, no `Z_x / Z_y` | 0.13 | 17.6 |
| corrected | 0.023 | 2.3 |

At 40 000 iterations and 4000 Monte Carlo draws per topology the corrected
chain gives TV 0.012, max |z| 1.6. Gate:
`tests/testthat/test-weighted-subtree-swap-invariance.R`: the replicate and
long-chain `beta = 0` checks, and this one.

## 6. Cost

The reverse neighbourhood costs `|P_y(A)| · nBins` further likelihood
evaluations, roughly as many as the forward one, so the move now costs about
twice as much. As with `gibbs_subtree_swap`, no reparameterisation makes
`Z_x = Z_y`: the neighbourhood is not an equivalence class.
