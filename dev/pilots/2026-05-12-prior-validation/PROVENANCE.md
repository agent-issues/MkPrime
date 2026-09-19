# Provenance — read before reusing anything in this directory

Audited 2026-09-18 against the open correctness issues on `agent-issues/MkPrime`.
Amended 2026-09-19: both open questions below are now answered, and the answer
to the second one is that EG-003's headline correlation was wrong for a reason
nothing here had considered.

## EG-003's ρ ≈ −0.1 was an artefact of the analysis, not of the sampler (#54)

`run_one.R` builds its character matrix from `sort(list.files(…, "^chr[0-9]+\.nex$"))`,
which is **lexical** — `chr1, chr10, chr11, …, chr2` — so `kPrime_i` is the i-th
lexically sorted file. `ground_truth.csv` is in **numeric** order.
`analysis/eg001_upost_compare.R` paired the two positionally, so every published
per-character statistic compared one character's posterior against a different
character's truth.

The falsifier needs no modelling assumption: `k' ≥ kObs` holds by construction,
so `u_post < 0` is impossible for a correctly aligned pair.

| | published pairing | corrected pairing |
|---|---|---|
| geo: `u_post < 0` | **248** / 1300 char-tasks | **0** |
| EG: `u_post < 0` | **115** / 1300 | **0** |
| geo Spearman | −0.119 | **+0.308** |
| EG Spearman | −0.111 | **+0.310** |
| geo mean `u_post` | 0.136 | 0.136 (exactly invariant) |
| EG mean `u_post` | 1.313 | 1.313 (exactly invariant) |

**Survives:** "`u_post` is prior-dominated" *for the geometric arm, as a
statement about location* — geo shrinks `u` to ≈ 0.14 against a `u_true` mean
of 1.246. EG is well calibrated in location.

**False:** "`u_post` does not track `u_true`". Both arms track at ρ ≈ +0.31.
The published figure was a permutation null by construction: the analysis could
not have detected tracking had there been any.

The means are *exactly* invariant because the error permutes a multiset. That
is the internal check on the correction, and it is why every marginal claim
survives while every paired one does not. **This is a re-analysis, not a
re-run** — the stored samples are fine, because the sampler does not care what
order the characters are in.

`analysis/eg001_upost_compare.R` now reorders the ground truth into the
sampler's order and asserts the zero-violation property, so the pairing cannot
silently regress.

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

## The #19 question is now answered — for the `u` marginal only

This section used to say: *do not cite EG-003 as settled without either re-running a subset
under the #30 fix, or arguing explicitly why the `u` marginal is insensitive to it.* The
re-run has been done.

`dev/pilots/2026-09-18-gspr-eg003-sensitivity/` runs a matched `gibbsSpr` TRUE/FALSE grid —
8 datasets × 3 seeds × 2 arms, 200k iterations — and bounds the effect on the per-character
`u` marginal at **Δ mean `u_post` ≤ 0.011** and **Δρ ≤ 0.022** at 95% confidence, at or
below the Monte Carlo noise floor. So the `u` marginal is insensitive to #19, as
conjectured.

**Not cleared:** the **topology** posteriors — `cid_mk`, `cid_mkp`, `cid_mkp_eg` and
`thinned_trees` — which were not compared between arms, and which is where #19 does its
damage. Treat those as contaminated.

Note the two findings are independent: #19 does not touch EG-003, and the reason to doubt
the published EG-003 correlation is the character-ordering artefact above, which has nothing
to do with the sampler.

Everything else here should be treated as contaminated unless the same argument is made for
it specifically.

## One further reason these numbers are not current

The `u_post` values predate the Model A (`priorVariant = "unconditional"`) flip in `2b3c054`
/ `5f12d9c`. Anyone wanting values from current code needs a re-run for **that** reason,
independently of everything above.
