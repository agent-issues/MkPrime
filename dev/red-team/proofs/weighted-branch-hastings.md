# WBS-001 — the Hastings ratios of `weighted_branch_scale` and `block_gibbs_branch`

Derived against `src/mcmc.cpp` at `c6b9fc9` (`fix/tree-move-hastings`, before
this fix), 2026-09-21. Fixed on `fix/weighted-move-hastings`. Part of issue #158.

**Verdict.** The selection normaliser cancels exactly in both moves, and so
do the bin weights themselves. The only defect is the one
`weighted-spr-hastings.md` found in `weighted_spr`: the fraction's proposal
density is a mixture over bins, and the code used one component forward and
the *containing* bin's component backward. Replacing both with
`bin_mixture_logdensity()` is the whole fix.

Both moves are **default-off** (`weightedBranchScale = FALSE`,
`blockGibbsBranch = FALSE` in `MkPrimeMCMC()`).

## 1. The kernel

`weighted_branch_scale_impl` (moveType 12), with `BranchBins` as in
`weighted-spr-hastings.md` §1 (`nBins` quantiles of `Beta(0.25, 0.25)`,
`mids[b]`, `conc = 2 nBins`):

1. Draw an ordered pair of distinct edges `(i, j)`: `i` uniform over `nEdge`,
   `j` uniform over the other `nEdge - 1`.
2. `s = r_i + r_j` (relative lengths), `f = r_i / s`.
3. For each bin `b`, `LL_b = LL(x with r_i = mids[b] s, r_j = (1 - mids[b]) s)`;
   `w_b = exp(beta (LL_b - max_b LL_b))`, zero where `LL_b` is not finite.
4. Draw `b' ∝ w_b`, then `f' ~ Beta(mids[b'] conc + 1, (1 - mids[b']) conc + 1)`.
5. Set `r_i = f' s`, `r_j = (1 - f') s`; accept by the generic MH step, which
   adds the prior ratio.

`block_gibbs_branch_sweep_impl` (moveType 15) visits every edge `i` in a random
permutation, draws a partner `j` uniformly from the rest, and applies steps
2–5 to the pair with an inline MH step.

## 2. Nothing but the fraction changes

Every input to the weights — the topology, `treeLength`, every `r_k` with
`k ∉ {i, j}`, and the sum `s` — is left unchanged by the move. So from
`y = x[f -> f']` the reverse move, using the same ordered pair, computes the
**same** `LL_b` multiset, the same offset, the same `w_b` and the same
`W = Σ_b w_b`. The pair probability `1 / (nEdge (nEdge - 1))` is also common.

With the bin integrated out (it is not carried in the state, and every
component has full support on `(0, 1)`):

    q(x -> y) = 1/(nEdge (nEdge-1)) · Σ_b (w_b / W) · Beta(f' | b)
    q(y -> x) = 1/(nEdge (nEdge-1)) · Σ_b (w_b / W) · Beta(f  | b)

**Jacobian.** Parameterise the slice `r_i + r_j = s` by the fraction. The map
`(f, u = f') -> (f', u' = f)` is a transposition, `|J| = 1`; the factor `s`
that converts a density in `f` to a density in `(r_i, r_j)` is the same at both
endpoints. Hence

    log HR = log Σ_b w_b Beta(f | b) - log Σ_b w_b Beta(f' | b)

with `W` and the pair probability cancelling. The code before this change had

    log w[bin(f)] + log Beta(f | bin(f)) - log w[b'] - log Beta(f' | b')

where `bin(f)` is the bin *containing* `f` — a point mass the forward draw
never uses (bin 1's component puts about 98% of its mass outside bin 1).

**Clamp.** `f'` is clamped to `[1e-8, 1 - 1e-8]`, an atom the mixture does not
describe. Its mass is below `2e-7` per draw at `nBins = 10`, as in
`weighted_spr`, and is ignored.

## 3. `block_gibbs_branch`

Each pair update is the kernel of §2 with the weights computed from the
sweep's *current* state, so each is π-reversible on its own. A random-scan
sweep is a composition of π-invariant kernels, hence π-invariant (though not
reversible); drawing the partner at random makes each step a mixture of
reversible kernels.

> **Guard.** The sweep's inline accept step omits the prior ratio. That is
> exact only because the relative-length prior is `Dirichlet(1, ..., 1)`,
> constant on the simplex (`cpp_log_prior`, `lgamma(nEdge)`). A non-flat
> branch prior would need the ratio added here; `weighted_branch_scale`
> already gets it from the generic MH step.

## 4. Verification

`pi K^n = pi` at `beta = 0`, where the target is the prior and is exactly
i.i.d.-sampleable (`tests/testthat/helper-prior-invariance.R`). In this table
each of 1500 replicates at `nTip = 6` starts from an exact draw; the reference
is 15000 i.i.d. draws.

| move (moves per replicate) | smallest fraction: uncorrected | corrected | target | KS p uncorrected | KS p corrected |
|---|---|---|---|---|---|
| `weighted_branch_scale` (25) | 0.00768 | 0.01222 | 0.01238 | **8.8e-50** | 0.99 |
| `block_gibbs_branch` (25 sweeps) | 0.00755 | 0.01260 | 0.01223 | **8.7e-49** | 0.55 |

The uncorrected kernels push fractions to the extremes: the smallest lands
38% below its target mean and the largest 14–17% above. The corrected
`weighted_branch_scale` was additionally run at 16 000 replicates on two seeds
(every KS p ≥ 0.038; the one low value was not reproduced on the other seed).

Gate: `tests/testthat/test-weighted-branch-invariance.R`, with 1000
replicates of 25 moves (`weighted_branch_scale`) or 5 sweeps
(`block_gibbs_branch`). On the uncorrected build its smallest KS p-values are
5e-22 and 3e-34. At `beta = 0` every weight is 1,
so the gate cannot see how the weights enter the mixture; that is pinned
against an independent R-side sum in `test-weighted-spr-invariance.R`, and the
weights' direction-independence is structural (§2).
