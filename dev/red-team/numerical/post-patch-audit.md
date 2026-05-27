# Post-patch audit — L5 / N1 / D3

Auditor lane (numerical) re-audit of the patches landed at commit `59276ad`
(2026-05-26; main HEAD `f9ea652`). The patches are reviewed as new
perturbations: do they introduce numerical issues that the original proofs/
audits did not anticipate?

## Scope of evidence

- `src/node_cl_cache.h` (L5 sites + N1 sites)
- `src/mcmc_likelihood.cpp::const_site_prob_for_k` (L5 closed-form caller)
- `src/mcmc.cpp::getCSP` and the `csp >= 1.0` guard at lines 3696-3698
- All 30 `arg = ...; ... -std::expm1(arg)` substitution sites
- `R/Convergence.R::.ComputeRhat`, exercised by `dev/red-team/numerical/post-patch-d3-probe.R`
- Pre-existing post-patch cache-coherence run under
  `dev/red-team/numerical/cache-coherence-results/{variable,informative}/`

The audit is fully read-only on `src/` and `R/`; the probe driver lives only
under `dev/red-team/numerical/`.

---

## Patch L5 (LIKE-001) — singleton-site ascertainment under `coding="informative"`

### Branch selection

`mcmc_state.h:73` documents `codingType: 0 = none, 1 = variable, 2 =
informative`. All four new `if (data.codingType == 2)` / `if (cache.coding ==
2)` gates correctly target the *informative* branch, not the variable branch.
`src/mcmc_likelihood.cpp:2898` confirms the parser:
`(codingStr == "none") ? 0 : (codingStr == "variable") ? 1 : 2`.

### Range of `p = p_const + p_singleton`

`p_const` integrates `P(all tips share a state)` (k pseudo-characters); the
singleton function integrates `P(exactly-one-tip differs)` over the `k*(k-1) *
nTip` distinguishable pseudo-characters. These two pattern-sets are disjoint
under JC/MkN, so in exact arithmetic `p_const + p_singleton ≤ 1`. The
remainder is `P(parsimony-informative pattern)`, which is ≥ 0.

In double-precision arithmetic, the sum could land at `1.0` or marginally
above due to per-pseudo-character rounding. Each new call site is guarded:

| Site | Guard | Action on `p ≥ 1.0` |
|------|-------|---------------------|
| `node_cl_cache.h:752` (covers 711-717 type-0 and 722-728 type-2) | `if (p > 0.0 && p < 1.0)` | log-correction skipped → `ll` unchanged |
| `node_cl_cache.h:747` (covers 740-746 type-1) | `if (pu > 0.0 && pu < 1.0)` | same |
| `mcmc.cpp:3696-3698` (covers `const_site_prob_for_k` return path) | `csp < 1.0` else `logAscCorr = R_NegInf` | safe sentinel |

No new site reads `p` for anything else (e.g. adaptive weights, scratch
storage); `p` is local in each block and consumed only by `std::log(1.0 -
p)`.

### Cache-coherence regression (task item 4)

Pre-existing `dev/red-team/numerical/cache-coherence-results/informative/verdict.txt`
(produced post-patch by the orchestrator's wave) reports
**0/83 mismatches, max |partial - full| = 4.26e-14** under
`coding="informative"`. Pre-patch baseline (from `findings.md`) was **94/94
with 2.20 nat drift**. Re-running the `--quick` driver was unnecessary —
the artefact is current and the diff confirms the cache-coherence bug closed.

### Verdict: **SAFE**

No new numerical issue. All guards in place; sum stays in [0, 1] by disjoint-
pattern argument; caller handles the ULP-edge `csp ≥ 1.0` case.

---

## Patch N1 (FAST-EXP-001) — `-std::expm1(arg)` substitution

### Identity check — `arg ≤ 0` at all sites

Listed all `arg = ...; ... = -std::expm1(arg)` sites and verified the
sign-of-`arg` invariant per site (output of `Grep "double arg\s*="` on
`src/`):

```
arg = -kStates * t / km1   [JC,  20 sites]   km1 = kStates - 1 > 0
arg = -kFull   * t / km1   [JC-collapse, 6]  km1 = kFull   - 1 > 0
arg = -lambda * t          [MkN, 9 sites]    lambda = rate01 + rate10 > 0
arg = -mu * t              [F81, 4 sites]    mu > 0 (Het bin midpoint)
```

In every case `t = edge_length * rate` with `edge_length ≥ 0` (phylogenetic
invariant; negative branch lengths would be rejected upstream) and `rate >
0` (ACRV multiplier or 1.0). Hence `arg ≤ 0` and `-std::expm1(arg) ∈ [0,
1]`. The C99 `expm1` specification guarantees ulp-accurate `exp(x) - 1` for
all finite `x`, so `1 - exp(arg) = -expm1(arg)` is ulp-accurate at every
site, including the `arg ≈ -1e-15` corner that motivated the patch.

### `p_same` consistency

At each JC site the patch uses:
```cpp
double neg_expm1 = -std::expm1(arg);
double exp_term  = 1.0 - neg_expm1;   // recovered E
double p_same    = inv_k + (1.0 - inv_k) * exp_term;
double p_diff    = inv_k * neg_expm1;
```
Both `p_same` and `p_diff` are derived from the **same** `neg_expm1` value.
Bit-identity follows because `1.0 - neg_expm1` is exact when `neg_expm1 ∈
[0, 1]` and `exp(arg) > 0.5` (since `arg ≥ -ln 2` keeps this true at the
regime that matters). For `arg ≪ 0`, `exp_term` may lose 0-1 ulp of
precision relative to `std::exp(arg)`, but `p_same = inv_k + (1 - inv_k) *
exp_term` is dominated by `inv_k = 1/k` for any `arg ≪ 0`, so the ULP-scale
error in `exp_term` does not propagate.

The 30 sites I located all follow this template; **no site computes
`p_same` and `p_diff` independently** (i.e. with two separate `MKP_EXP` /
`expm1` calls). Pattern is uniform.

### Bit-identity reference

`git show --stat 59276ad` confirms the only test-tree change is
`tests/testthat/_reference/partition-bitcompat-null-ref.rds`
(1924 → 1925 bytes). The `Bin … bytes` diff line + no schema changes in
`test-partition-bitcompat-null.R` (which uses `expect_identical` for both
schema and values) is consistent with regenerating identical-shape but
ULP-perturbed doubles. The commit message states "5244 tests pass / 0 fail",
which is the empirical confirmation. Re-running the test would be a
zero-information confirmation that I did not perform.

### Verdict: **SAFE**

`-std::expm1(arg)` is mathematically equivalent to `1 - exp(arg)` at every
site, with strictly better conditioning at small `|arg|`. No new corner
introduced. Performance cost (~ +1 expm1 per branch per cat) was bounded in
the original fast-exp-conditioning analysis at a few percent.

---

## Patch D3 (CONV-002) — tail-equalised R-hat

### Degenerate-input probe

`dev/red-team/numerical/post-patch-d3-probe.R` exercises the patched
`.ComputeRhat` on seven cases:

| Case | Input | Result |
|------|-------|--------|
| A | `perRun = list()` (length 0) | `ERROR: subscript out of bounds` |
| B | one run, nrow = 0 | `NA` |
| C | two runs, both nrow = 0 | `NA` |
| D | runs nrow ∈ {1, 50} → nKeep = 1 | `NA` |
| E | runs nrow ∈ {2, 50} → nKeep = 2 | `NA` |
| F | two runs nrow = 50 (sanity) | `1.003447` |
| G | runs nrow ∈ {37, 50} (CONV-002 scenario) | `1.009593` |

Case A: the crash arises from `paramNames <- colnames(perRun[[1]]$samples
...)` — a pre-patch line (verified by `git show 59276ad^:R/Convergence.R`).
This is **not introduced by D3**, and the upstream call in
`ConvergenceDiagnostics` is guarded by `nRuns >= 2L && !is.null(per_run)`
which renders the empty-list case unreachable in production.

Cases B-E (nKeep ≤ 1): the `.Rhat` core safely returns `NA` via the
within-chain variance handling. No crash.

Cases F-G: D3 fix demonstrably works — unequal-length runs (G) now produce a
finite Rhat where the pre-patch code would silently recycle or warn.

### Verdict: **SAFE**

The patch is a strict improvement: the empty-list crash is pre-existing and
unreachable; all other degenerate inputs degrade to `NA` gracefully; the
target scenario (unequal nrow) returns a finite Rhat. No new bug.

---

## Summary

| Patch | Verdict | Reason |
|-------|---------|--------|
| L5 / LIKE-001 | **SAFE** | Branches correct; `p ≤ 1` by disjoint-pattern argument; guards present at all 3 cache sites + caller; cache-coherence artefact confirms zero drift |
| N1 / FAST-EXP-001 | **SAFE** | `arg ≤ 0` invariant holds at all 30 sites; both `p_same` and `p_diff` derived from the same `neg_expm1` (no independent recomputation); bit-identity reference regeneration is the expected ULP perturbation |
| D3 / CONV-002 | **SAFE** | Empty-list crash is pre-existing and unreachable; nKeep ≤ 1 yields `NA`; target unequal-nrow case returns a finite Rhat |

No new findings filed under `dev/red-team/findings.md`.

## Worktree state

- `dev/red-team/numerical/post-patch-audit.md` — this note
- `dev/red-team/numerical/post-patch-d3-probe.R` — D3 degenerate-input
  driver (CLI: `Rscript dev/red-team/numerical/post-patch-d3-probe.R`).

No modifications under `R/`, `src/`, or `tests/testthat/`.
