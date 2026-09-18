# Provenance — read before reusing anything in this directory

Audited 2026-09-18 against the open correctness issues on `agent-issues/MkPrime`.

## Produced by a sampler that was not π-invariant (issue #19)

Driver: `data-raw/hamilton/run_one.R`. It calls `MkPrimeMCMC()` (line 164) with **no
`fixTopology`** and **no `gibbsSpr` override**, so every run in `summary/` was made under
free topology with `gibbsSpr = TRUE` — the exact configuration issue #19 covers.

`gibbs_spr` committed a deterministic `0.5 * lReg` edge split with no Metropolis–Hastings
step. Under the flat Dirichlet branch prior the set of trees with two exactly-equal
incident edges is π-null, and every accepted move landed in it, so the chain did not target
the posterior. The fix is PR #30.

**What that distorts:** adjacent edge-length fractions at regraft sites, pulled toward
equality; and the **topology posterior**, because candidates were scored at the τ = ½ point
value while the current state was scored at its adapted lengths — under-weighting
topologies whose best attachment fraction is far from mid-edge. `tree_length` is affected
only at second order (the sum is preserved). Magnitude unquantified.

## NOT affected by issues #25 / #26

Checked, not assumed: the `mkp_eg` streamed logs carry no `rate_neo` column, so these
datasets contain **no neomorphic characters**. `compute_partition_scales` therefore returns
`(1, 1)` and the missing transformational scale factor is identically zero here. The `mk`
arms additionally pin `k` through `knownStates`, so the case-25 Gibbs k′ sweep — the second
symptom of #25 — never runs for them either.

## The judgement this directory needs

The headline result drawn from these runs is **EG-003**: `u_post` is prior-dominated, with
`Spearman(u_post, u_true) ≈ −0.1` under both the geometric and empirical-geometric priors
(see `analysis/eg001_upost_compare.{R,rds,png}`). That is a statement about a
**per-character discrete marginal**, while #19 distorts **branch fractions and topology**.

The conclusion is therefore plausibly robust to the defect — but plausibly is not
established, and nobody has checked whether a distorted topology posterior feeds back into
the `u` marginal. **Do not cite EG-003 as settled without either re-running a subset under
the #30 fix, or arguing explicitly why the `u` marginal is insensitive to it.**

Everything else here should be treated as contaminated unless the same argument is made for
it specifically.
