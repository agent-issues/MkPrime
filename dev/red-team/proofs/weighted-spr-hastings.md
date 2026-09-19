# GSPR-003 — the Hastings ratio of `weighted_spr`

Checked against `src/mcmc.cpp` at `66510e1` (`origin/main`, 2026-09-19).
Resolves issue #20.

**Verdict, in two parts.** The term issue #20 suspected — a length-scaling
factor — is *present and correct*: the Jacobian `lReg / lMerge` is there, and
the selection-normaliser cancellation that the source comment called
"ASSERTED, NOT VERIFIED" is exact. What is missing is different: the
branch-fraction proposal density is a **mixture over bins**, and the code
used a single component of it.

## 1. The kernel

`weighted_spr_impl`, with `BranchBins` (`src/mcmc_state.h:35-55`): `nBins`
(default 10) quantiles of `Beta(0.25, 0.25)`, `mids[b]` the bin midpoint,
`concentration = 2 * nBins`.

1. Draw a prune edge uniformly from the edges whose parent is not the root.
   That leaves `nEdge - 3` of them — the root trifurcates — the same count in
   both directions.
2. Prune `(u -> v)`; `parentRow = (p -> u)` of length `l_p`, `sibRow` of
   length `l_s`, `lMerge = l_p + l_s`, `fOld = l_p / lMerge`.
3. For each of `{self}` and every candidate edge, form a marginal weight
   `m = sum_b exp(beta * LL(configuration at bin midpoint b))`.
4. Draw a candidate proportional to `m` (a self-draw is a no-op), then a bin
   `b_new` proportional to `candW[chosen][b]`, then
   `f_new ~ Beta(mids[b_new] * conc + 1, (1 - mids[b_new]) * conc + 1)`.
5. Accept by MH.

## 2. The normaliser cancels — exactly

Let `R` be the residual tree after deleting `subtree(v)` and suppressing `u`,
carrying the merged edge `(p -> sib)` of length `lMerge`.

**The candidate filter is `E(R)` minus the merged edge.** It drops every edge
whose child lies in `subtree(v)` (that is `pruneRow` plus `subtree(v)`'s
internal edges), then `sibRow` (`parent == u`) and `parentRow`
(`child == u`). Every surviving edge is an edge of `R` with an unchanged
length. And "self" — which re-splits `lMerge` across `parentRow`/`sibRow` —
*is* regrafting onto the merged edge. So the forward enumeration is exactly
`E(R) x bins`.

**The reverse produces the same `R`.** In `y`, `u` sits on the old edge `rr`,
split into `f_new * lReg` and `(1 - f_new) * lReg`, while `parentRow` carries
`lMerge`. Pruning `v` from `y` merges those two back into `lReg`, restoring
edge `rr`, and leaves the merged edge at `lMerge`. So `R_y = R_x` with
identical lengths, the reverse enumerates the same `E(R) x bins` over the
same trees, and

    candW_y[merged][b] = selfW_x[b],    sumM_y = sumM_x,

term for term — the shared `exp(-beta * globalMax)` offset included, since
the two log-likelihood multisets are identical.

Note what the reverse move *is*: reaching `x` from `y` is a **candidate**
draw, onto the merged edge, not `y`'s self-draw. That is why `selfW` legibly
carries the reverse weight.

So the `!! INCOMPLETE` comment removed by this change was wrong about its own
prime suspect. Its worry — "self is enumerated over bins of `lMerge` while
candidates are over bins of `lReg`" — is precisely *why* the cancellation
works: for the merged edge, the regraft edge length genuinely is `lMerge`.

## 3. What is actually missing

The bin is an auxiliary variable and is **not carried in the state**, so it
must be integrated out. The Beta components have full support on `(0, 1)`, so
every bin could have produced the fraction that was drawn. Writing
`elig = nEdge - 3`:

    q(x -> y) = (1/elig) * (mCand[c]/sumM) * sum_b (candW[c][b]/mCand[c]) * dbeta(f_new | b)
              = (1/(elig * sumM)) * sum_b candW[c][b] * dbeta(f_new | b)

    q(y -> x) = (1/(elig * sumM)) * sum_b selfW[b] * dbeta(f_old | b)

`mCand[c]` and `mSelf` cancel *within* each direction, between the
topology-selection factor and the bin-selection factor; only `sumM_y = sumM_x`
from §2 is needed across directions. Hence

    logHR = log( sum_b selfW[b] * dbeta(f_old | b) )
          - log( sum_b candW[c][b] * dbeta(f_new | b) )
          + log(lReg) - log(lMerge)

which is what `bin_mixture_logdensity` now computes.

**Orientation check.** The fraction always denotes the *above* part: the
forward candidate sets `absLen[rr] = mids[b] * lReg` on the edge `par(rr) -> u`,
the reverse self sets `trialAbs[parentRow] = mids[b] * lMerge` on `p -> u`,
and `fOld = absLen[parentRow] / lMerge` is read before any mutation. No
`1 - f` is needed anywhere.

**Jacobian.** On the slice `l_p + l_s + l_Reg = S` (`treeLength` is
untouched), with `a = f_new * (S - lMerge)`,

    |d(lMerge, a, f_old) / d(l_p, l_s, f_new)| = lReg / lMerge

matching the block-diagonal derivation already in the source. It is
scale-invariant, so identical in relative or absolute coordinates, and the
prior ratio is identically zero for this move (uniform simplex, fixed
`treeLength`), so nothing is double-counted.

## 4. How big is the omission?

With `nBins = 10`, `conc = 20`, the components are nowhere near confined to
their own bins:

| bin | interval | mass inside own bin |
|---|---|---|
| 1 | [0.0000, 0.0012] | 0.023 |
| 3 | [0.0187, 0.0905] | 0.501 |
| 5 | [0.2548, 0.5000] | 0.766 |

Even with equal weights, `log(mixture / single component)` has median 0.54
and reaches 1.89; with peaked weights (log-sd 3) the 95th percentile is 1.66
and the maximum 8.6 log-units. Because bin 1's component puts 97.7% of its
mass *outside* bin 1, the forward-sampled bin and the backward
containing-bin are frequently different bins — the two directions were not
even applying the same rule, which also rules out the "retain the bin as an
auxiliary variable" escape: that construction needs the reverse to draw its
auxiliary from the same conditional, not to take a point mass at the
containing bin.

## 5. Verification

`pi K^B = pi` at `beta = 0`, where the target is the prior —
Uniform(labelled unrooted topologies) x Dirichlet(1, ..., 1) on the edge
fractions — which is exactly i.i.d.-sampleable. Each of 1500 replicates
starts from an exact draw and takes 25 moves; the reference is 15000 i.i.d.
draws. `nTip = 6`.

| statistic | uncorrected mean | corrected mean | target mean | KS p (uncorrected) | KS p (corrected) |
|---|---|---|---|---|---|
| smallest fraction | 0.00919 | 0.01242 | 0.01224 | **3.7e-22** | 0.35 |
| largest fraction  | 0.33117 | 0.31633 | 0.31422 | **6.0e-11** | 0.15 |
| internal total    | 0.33231 | 0.33001 | 0.33336 | 0.020 | 0.47 |

The uncorrected kernel drives fractions toward the extremes: the smallest
comes out 25% below its target mean. Gate:
`tests/testthat/test-weighted-spr-invariance.R`.

## 6. Not fixed here

`weighted_branch_scale_impl`, `block_gibbs_branch` and
`weighted_subtree_swap_impl` build their Hastings ratios from the same
single-component pattern and are very likely to carry the same defect. Their
*normaliser* structure has not been derived, and a partial fix that looks
complete is worse than a known-open one, so they are left alone: issue #158.
