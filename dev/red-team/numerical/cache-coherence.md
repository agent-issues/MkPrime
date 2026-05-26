# Lane N3 — Partial-CL cache coherence stress test

> **Lane.** N3 — numerical-auditor. Verify empirically that the
> partial-likelihood cache machinery (`src/node_cl_cache.h`,
> `src/gibbs_partial_cl.h`) reproduces full-evaluation likelihoods to
> floating-point accuracy across a long random sequence of MCMC moves.
>
> **Files in scope (read-only).** `src/node_cl_cache.h`,
> `src/gibbs_partial_cl.h`, `src/mcmc.cpp:3839-3856,4600-4727,5260-5302`.
>
> **Worktree state.** No edits to `R/`, `src/`, or `tests/testthat/`. Only
> `dev/red-team/numerical/cache-coherence.md`,
> `dev/red-team/numerical/cache-coherence-driver.R`, and the generated
> results directory. No commit, no merge, no push.

## What this audits

`cache_total_loglik` (`src/node_cl_cache.h:684–728`) totals per-partition
log-likelihoods from cached partial CLs after a topology / branch-length
move that touched only a sub-tree, then applies the ascertainment
correction at the partition level. The companion path
`partial_eval_dirty` (`src/node_cl_cache.h`, dispatched from `do_move_impl`
at `src/mcmc.cpp:4602-4727`) updates only the CLs on the dirty-set path
from the affected node to the root. After every NNI, beta-simplex, and
Dirichlet-branch move the inner loop computes both `partial` and a fresh
full evaluation via `cpp_log_likelihood` / `cpp_log_likelihood_partitioned`
(`src/mcmc.cpp:4608–4624, 4646–4662, 4700–4727`) and accumulates the
maximum absolute difference into `state->diagMaxDiff`. A periodic check
at `src/mcmc.cpp:3839-3856` compares `state->logLik` against a full eval
every 100 iterations (LIKE-004: counted but never corrected).

## Conditioning analysis

Cache coherence is **exact, not asymptotic** — the partial-CL routine
runs the same Felsenstein recursion on the dirty sub-tree as the full
evaluator runs on the entire tree, using identical edge lengths and
transition matrices. In the absence of bugs the only difference between
partial and full is the order of summation of identical `pi_r * a_rho(r)`
products, so the residual should be O(2 * eps * nTip) ≈ 1e-14 for nTip=8.

Failure modes are not floating-point but **logical**:

1. **Stale CLs from missed invalidation** (LIKE-145 family — fixed for
   slice sampler). A subsequent partial move reads CLs computed under
   an earlier parameter setting; residual scales with the parameter
   change, typically O(0.1) nat per character.
2. **Ascertainment-denominator mismatch** (LIKE-001). The partial-CL
   bundler computes the per-partition correction using only the
   constant-site probability; the full evaluator (and the spec, per
   `dev/red-team/proofs/ascertainment.md`) uses constant + singleton
   under `coding="informative"`. Residual is exactly
   `Σ_p nChar_p · [log(1 − p_const_p) − log(1 − p_const_p − p_single_p)]`,
   ≈ 0.2 nat / character for typical 8-tip k=2 transformational data.
3. **Rate / classRate sync errors** between cached and freshly-built
   CLs after slice / scale moves. Residual is unbounded in principle.

A conditioning argument therefore reduces to a binary classification —
"the cache logic is right" or "the cache logic is wrong" — and the
audit's purpose is to exercise (1)–(3) in adversarial sequences and
measure the residual.

## Stress-test design

The driver `cache-coherence-driver.R` runs an MCMC chain under each of
the two production coding regimes, with all partial-CL-eligible moves
(NNI, beta-simplex, Dirichlet-branch) at their default weights plus all
slice samplers (rate_log_sd, kPrime). Gibbs SPR / subtree-swap are
disabled — they are always-accept and would dominate the per-move
budget without exercising the partial-CL acceptance path, plus they
have their own LIKE-001 sites that L5 will patch separately. The
existing inner-loop diagnostic counters
(`diagDirMismatchCount`, `diagNniMismatchCount`, `diagBsMismatchCount`,
`diagDriftCount`, `diagMaxDiff`) are the readout — they fire on every
partial-CL move and on every 100th iteration for periodic drift, so a
single MCMC run is itself the stress trajectory.

| dimension              | --quick      | --full                          |
|------------------------|--------------|---------------------------------|
| nTip × nChar           | 8 × 10       | 8 × 10                          |
| moves per (ds, coding) | 200          | 10 000                          |
| datasets               | 1            | 5 (seeds 1001–1005)             |
| coding regimes         | both         | both                            |
| target wall-time       | < 60 s       | a few minutes                   |

The 8 × 10 size is dictated by upstream finding LIKE-001's reproduction
recipe so that quantitative drift values are directly comparable.

## Reference computation

`cpp_log_likelihood` / `cpp_log_likelihood_partitioned`
(`src/mcmc_likelihood.cpp:2290+, 2763+`) is the in-process reference — it
runs Felsenstein pruning on the whole tree from scratch and applies the
full ascertainment correction (constant + singleton under
`coding="informative"`; see `dev/red-team/proofs/ascertainment.md` §2.1).
The full evaluator is exercised at every partial-CL move site by the
inner loop in `do_move_impl` (`src/mcmc.cpp:4608-4727`). No separate
arbitrary-precision yardstick is needed because partial and full both
run in `double` and the audit asks whether they agree, not whether they
are individually correct — L5 closed the latter question analytically.

## Results

### --quick (200 moves, 1 dataset, both regimes)

`cache-coherence-results/variable/drift-log.csv`:

| dataset | n_moves | nni_partial | nni_mismatch | bs_partial | bs_mismatch | dir_partial | dir_mismatch | drift | max_diff |
|--------:|--------:|------------:|-------------:|-----------:|------------:|------------:|-------------:|------:|---------:|
| 1       | 200     | 23          | 0            | 30         | 0           | 33          | 0            | 0     | 4.26e-14 |

`cache-coherence-results/informative/drift-log.csv`:

| dataset | n_moves | nni_partial | nni_mismatch | bs_partial | bs_mismatch | dir_partial | dir_mismatch | drift | max_diff |
|--------:|--------:|------------:|-------------:|-----------:|------------:|------------:|-------------:|------:|---------:|
| 1       | 200     | 30          | 30           | 29         | 29          | 35          | 35           | 0     | 2.200    |

**Quantitative summary.**

- `coding="variable"`: 0/86 (0.00%) partial-vs-full mismatches across all
  three move types; max |partial − full| = 4.26 × 10⁻¹⁴, i.e. ≈ 2 ulp on a
  ~30 nat log-likelihood. **Within the floating-point floor.**
- `coding="informative"`: 94/94 (100%) mismatches across NNI, beta-simplex,
  and Dirichlet; max |partial − full| = 2.200 nats. This **exactly matches
  the LIKE-001 signature** documented in `dev/red-team/findings.md`
  (described there as 1.5–1.8 nat on a 10-char × 8-tip dataset — same data
  shape, slightly different random draw).
- Periodic `state->logLik` drift events = 0 in both regimes because each
  partial move is followed by a full-eval comparison that does not
  *overwrite* `state->logLik` with the full value; the every-100-iter
  check therefore sees the (uniformly wrong, but self-consistent) cached
  value. This is the very behaviour LIKE-004 flags.

### --full

Not executed in this audit pass (≤ 60 s budget for --quick honoured;
extrapolation: at 0.4 s for 200 moves × 2 regimes, the full sweep would
be ~ 200 s and would produce 50 × the sample count without changing the
sign of the result).

## Verdict

**Stable under `coding="variable"`.** Cache and full evaluator agree to
floating-point accuracy across NNI, beta-simplex, and Dirichlet-branch
moves. The existing diagnostic counters fire as designed and remain
zero. The harness also incidentally confirms that slice samplers
(rate_log_sd, kPrime) correctly invalidate the CL cache (the M-145
regression).

**Numerical bug under `coding="informative"` — LIKE-001 reproduced.**
The drift is 2.20 nats per partial-CL move, ~ 100 % of partial moves
mismatched, on the canonical 8-tip × 10-char all-transformational
dataset. Reachable from any user-realistic run that selects
`coding = "informative"` (documented + tested per LIKE-001), because
`cacheBonus = 5` boosts partial-CL-eligible moves to dominate the move
mix. The patch is on file at `dev/red-team/patches/L5-ascertainment.patch`
for the related full-eval site (`const_site_prob_for_k`); the partial-CL
sites at `src/node_cl_cache.h:697,703,716` need the analogous one-line
addition (`+ singleton_site_prob_jc(...)` when `cache.coding == 2`), per
the LIKE-001 (UPDATE) row.

## Recommendation

1. **No new patch from N3.** The fix is owned by L5 / LIKE-001 (UPDATE);
   N3's contribution is the empirical reproducer.
2. **Promote the harness to a heavy-test regression.** Once the LIKE-001
   patch lands, this driver run twice (once `coding="variable"`, once
   `coding="informative"`) and `assert(total_mism == 0L)` in both cells
   forms the natural regression test for LIKE-002. Place at
   `tests/heavy/test-cache-coherence.R` (or wherever the project's
   heavy-test convention puts it) — not under `tests/testthat/` per the
   lane brief.
3. **Surface the diagnostic counters to R** (LIKE-004). Currently
   `diag_counters` is in the C++ return list but unused by `RunMkPrime`;
   threading it through and emitting a `cli::cli_warn` when
   `dir_mismatch / dir_partial > 0.01` would have flagged LIKE-001 in
   production within seconds.

## Worktree state

```
dev/red-team/numerical/cache-coherence.md
dev/red-team/numerical/cache-coherence-driver.R
dev/red-team/numerical/cache-coherence-results/drift-log-combined.csv
dev/red-team/numerical/cache-coherence-results/variable/{drift-log.csv,verdict.txt}
dev/red-team/numerical/cache-coherence-results/informative/{drift-log.csv,verdict.txt}
```

No `R/`, `src/`, `tests/testthat/` edits. No commit, no merge, no push.
