# Review: fix/partition-rate-normalisation @ c739d0c

(base: af58922, the merge-base with main)

## Verdict

**Block.** The C++ `cpp_partition_log_likelihood`, the node-CL cache, and the R-side `.MkpLogLikelihood` are correctly patched. But the **partition-likelihood cache (`state->partLogLik`) is updated incorrectly for rate_neo (move type 3), neo_joint (move type 18), and the rate_neo slice sampler (paramIdx = 3)**. These paths still recompute only `data->neoPartIndices` after a rate_neo change, even though under the new RB-style formula rate_neo *also* shifts `transScale` and therefore changes the trans / known partitions' log-likelihoods. The MH acceptance ratio is computed against a stale `newLogLik`, biasing the chain. Cache-invalidation logic in `do_move_impl` and `slice_scalar_impl` was widened correctly; the analogous logic in the *partial-LL* recompute branch was missed.

## Spec recap

`dev/rb-equivalence/notes/partition-rate-and-acrv-audit.md` (§Issue 1) requires MkPrime to adopt RB's symmetric two-sided scaling `neoScale = r/(1+r) · n/n_neo`, `transScale = 1/(1+r) · n/n_trans`, applied consistently at every site where the old one-sided `edgeLen * rate_neo` (or unscaled `edgeLen` for trans) was used. The fix must extend cache invalidation: rate_neo now affects *all* units, not just neomorphic. The Bactrian-scale / slice proposals operate on the parameter, not derived per-edge quantities, so Hastings ratios stay unchanged — that's correct.

## Findings

### F1. partLogLik partial-cache update on rate_neo / neo_joint moves is stale on trans — severity: **blocker**

**Where:** the three partial-LL branches that share the predicate `(moveType == 1 || 3 || 18)` or `(paramIdx == 1 || 3)`. **Case 1 / paramIdx 1 (rate_loss) stays correct** — rate_loss only enters via MkN stationary frequencies and the Q-matrix, both neo-only. The bug is specifically that **cases 3 (rate_neo) and 18 (neo_joint) and paramIdx 3** were left inside that "neo only" branch when the new RB-style scaling also makes them affect trans partitions:
- `src/mcmc.cpp:4815-4831` — `do_move_impl` partial-LL branch: `case 1: case 3: case 18:` recomputes only `data->neoPartIndices` and assumes `newPC[pi]` for trans is unchanged. After the patch, transScale = 1/(1+r) · n/n_trans is r-dependent, so trans partition log-likelihoods *do* change when r changes. Cases 3 and 18 must move out of this branch; case 1 should remain.
- `src/mcmc.cpp:3083-3093` — `eval_slice_target` `(paramIdx == 1 || paramIdx == 3)` branch; same construction. paramIdx 3 must split off.
- `src/mcmc.cpp:3163-3175` — `slice_scalar_impl` accept-branch (same predicate); only neo partitions written into `state->partLogLik`.

**Claim:** When `state->partLogLik` is populated (which is the steady state — `fill_partition_cache` populates it at startup, `do_move_impl` default branch repopulates it after every full eval, see line 4848-4854), a rate_neo Bactrian scale move (case 3), the neo_joint move (case 18), or a rate_neo slice update reuses the stale trans contribution in `newLogLik`. The MH acceptance ratio is therefore wrong; the chain samples from a biased posterior.

**Why it matters:** This silently undoes the patch's intended correction on every mixed (neo+trans) dataset. Acceptance probabilities are biased; rate_neo and tree_length posteriors will *not* match the RB-equivalent target the patch is meant to deliver. The behaviour is most acute when partial-LL caching is on (default).

**Evidence:** Source-trace only — building the package on Windows requires a full Rtools toolchain run that exceeds the review budget. The diagnosis is unambiguous, however:

- `neoPartIndices` is populated only with type==0 partitions: `src/mcmc_likelihood.cpp:2958` (`if (pinfo.type == 0) d->neoPartIndices.push_back(pi);`).
- The cache-invalidation paths at `src/mcmc.cpp:3197-3201` and `src/mcmc.cpp:4887-4892` were correctly widened from `invalidate_neo_cls` → `invalidate_all_cls` for case 3 / 18 — but those widen the *node-CL* cache, not the *partition-LL* cache. The two caches are independent (`partLogLik` is a `std::vector<double>` of partition log-likelihoods; `nodeCL` is the per-internal-node CL tensor). The patch widened one but not the other.
- The non-patched analogous paths for `rateLogSd` (line 4893, "all units stale") and topology / `tree_length` moves (default branch at 4848-4854) correctly recompute every partition — that's the contract this patch needed to extend to case 3 / 18.

Reproduction sketch left under `dev/red-team/reviews/fix-partition-rate-normalisation/repro-01-partial-cache-bug.R` (script computes the partition scales at two rate_neo values to make the divergence concrete; full build-and-run not executed due to host toolchain).

**Suggested fix:** In all three sites, when `partLogLik` is populated and the changed parameter is rate_neo (case 3, case 18, or paramIdx == 3), recompute **every** partition log-likelihood, not just `neoPartIndices`. The simplest path mirrors what the default branch in `do_move_impl` (line 4848-4854) already does for tree_length and topology moves. The neo_joint move (case 18) also touches `rateLoss`, but `rateLoss` only affects neo partitions — so neo_joint needs the full sweep too, driven by the rate_neo component.

### F2. Tests are weak — half are tautological, the other half soft — severity: **major**

**Where:** `tests/testthat/test-partition-rate.R`.

**Claim:** Several tests would have passed on `main` (pre-fix); the remaining tests don't tightly bound what the C++ code is doing.

- Lines 25-65 (`partition rate weighted mean = 1`, `degenerates correctly`, `matches RB factorisation`): all three exercise the **R helper `.partition_scales`** defined at the top of the test file (lines 12-23). They do not call into `compute_partition_scales` in `src/mcmc_state.h`. Pure R-on-R identity — would pass with no C++ change at all. Regression-noise, not regression-detectors.
- Lines 107-144 (R↔C++ direct-eval parity): `cpp_log_likelihood_xptr` does route through `cpp_log_likelihood` → `cpp_partition_log_likelihood` and so exercises the new `compute_partition_scales` helper on the C++ side (`src/mcmc_likelihood.cpp:3049-3061`). The R side was also patched. So this test is a real regression detector against a *half-revert* (someone reverting just R or just C++) — that has value. But the test would still pass if both sides were wrong in the same way (e.g. if both implemented `neoScale = r · n/n_neo` and `transScale = n/n_trans`, mathematically wrong but symmetric). It guards consistency, not correctness of the math.
- Lines 153-179 (`cache (via short MCMC) and direct-eval agree at chain end`): the assertion is `expect_true(all(is.finite(res$samples[, "log_post"])))` — a finite-value check, not an equality check. Per the external-reviewer profile, "expect_true(is.finite(x))" is almost never the right contract for a numerical patch. This test would pass even if `partLogLik` were arbitrarily wrong, as long as values stayed finite. And the test is skipped on CI (`skip_if_not(... !nzchar(Sys.getenv("CI", "")))`) — so it doesn't run in the place it would catch the F1 bug.
- Lines 194-216 (`rate_neo affects likelihood via BOTH partitions`): asserts only that `ll(0.1) ≠ ll(1) ≠ ll(10)`. Pre-fix, rate_neo scaled neo edges, so changing rate_neo changes the neo partition's contribution — `ll_low ≠ ll_one ≠ ll_high` was already true on `main`. Tautology.
- Lines 222-262 (`extreme rate_neo: trans pruning sees small effective branches`): the only test that genuinely tests the trans-side effect. It uses 1 neo + 8 highly-informative trans chars; the claim `expect_lt(ll_huge, ll_one - 1)` likely fires post-fix because transScale collapses. Pre-fix, the single neo character also saturates at rate_neo=1e4, so the same inequality may *also* fire pre-fix — the test does not isolate the mechanism. To be a real regression test it should be run against a hidden ground truth (e.g. a transformational-only sub-eval with edges pre-scaled by the expected transScale).

**Coverage gaps:** None of the new tests exercises (a) `partLogLik` partial-cache updates on case 3 / 18, (b) the slice sampler accept path, (c) the node-CL cache invalidation on a rate_neo move under a long enough chain that the cache is in-state. Each of these is a patch site the author edited.

**Suggested fix:** add a partial-cache regression test: populate the cache (via `fill_partition_cache`), mutate `state->rateNeo`, recompute the *full* likelihood directly, and assert equality with the partial-cache-aware update path. If the F1 bug fix lands, this test becomes the regression detector for it. (Implementation requires exposing the partial update path to R, which the package already does via `cpp_log_likelihood_xptr` plus the `partLogLik` field — a 30-line testthat addition.)

### F3. cpp_log_likelihood_partitioned path: undocumented behaviour change for hasNeo == true — severity: minor

**Where:** `src/mcmc_likelihood.cpp:2867-2871`. The partitioned-API calls `cpp_partition_log_likelihood(..., /* rateNeo */ 1.0, ...)`.

**Claim:** With rateNeo = 1 and `hasNeo == true`, the new `compute_partition_scales(1.0, nNeo, nTrans)` returns `neoScale = 0.5 · n/n_neo`, `transScale = 0.5 · n/n_trans` — **not** (1.0, 1.0). So at `rateNeo = 1` plus `hasNeo == true`, the partitioned path now produces a different per-partition effective edge length than the pre-fix partitioned path did. The comment at line 2816-2820 says "the normalisation reduces to multiplication by 1.0" only in the `hasNeo == false` case, but reads as though it covers both.

**Why it matters:** The partitioned API (Layer 1) is the partition-aware code path used by the Casali driver. The commit message claims `hasNeo == false production workloads (Casali) are unaffected` — true — but says nothing about hasNeo == true users of the partitioned path. Currently they get a posterior shift identical in magnitude to the non-partitioned-path shift; the *intent* of the fix is that this shift is the *desired* RB-style behaviour, but the comment hasn't been updated to say so.

**Suggested fix:** Either edit the comment at `src/mcmc_likelihood.cpp:2812-2832` to state explicitly that hasNeo == true partitioned callers also see a deterministic posterior shift identical to non-partitioned hasNeo == true, or update the eta_neo-aware partitioned path (currently deferred — see the Layer 1 comment at 2828-2832) to thread `etaNeo` into `compute_partition_scales` so the partitioned semantics remain consistent with the §5.2 design intent.

### F4. Behaviour-change paper trail is partial — severity: minor

**Where:** Commit message + `dev/red-team/findings.md`.

**Claim:** The behaviour change is documented in the commit message body — good. But:
- There is no `NEWS.md` entry (the file does not exist in the repo, so this is debatable, but the memory note `feedback_no_oversample.md` precedent says behaviour changes warrant an entry somewhere durable).
- The R-side `MkpLogLikelihood` function comment is updated (`R/likelihood.R:90-95`) — good.
- The C++ `cpp_partition_log_likelihood` carries inline rationale at `src/mcmc_likelihood.cpp:2341-2348` — good.
- No memory note (`~/.claude/projects/.../memory/`) records that posterior samples on mixed datasets pre-c739d0c are not comparable to post-c739d0c samples. This will catch out a user comparing a smoke-run posterior from yesterday with one from tomorrow. The audit document already records the magnitude (~31% tree_length inflation pre-patch on the 635 by_nt_9v dataset), so the paper trail exists in `dev/rb-equivalence/notes/partition-rate-and-acrv-audit.md` — the gap is just that there is no inbound link from `findings.md` or the user's memory to that note for the "incomparable posteriors" claim.

**Suggested fix:** add a short stanza to `feedback_*.md` memory or a NEWS entry: `fix/partition-rate-normalisation changes the meaning of tree_length and rate_neo on mixed (neo+trans) datasets; pre- and post-c739d0c posterior samples are not directly comparable on those parameters. Trans-only and neo-only datasets are unchanged.`

## Coverage notes

**Patch sites exercised by new tests:**
- The R `.MkpLogLikelihood` direct-eval path (lines 107-144) — covered.
- The C++ `cpp_log_likelihood` / `cpp_partition_log_likelihood` direct-eval path (same test) — covered.
- Degenerate (trans-only) collapse to identity scales (lines 270-289) — covered.

**Patch sites NOT exercised by new tests:**
- `state->partLogLik` partial-cache update paths in `do_move_impl` (line 4815-4831) — **uncovered, where F1 lives.**
- `slice_scalar_impl` partial-cache update paths (line 3083-3093, 3163-3175) — **uncovered.**
- `populate_cache_full` / `populate_cache` symmetric rateScale setting (`src/node_cl_cache.h:583-592, 638-649`) — only indirectly via R↔C++ parity.
- `cache_total_loglik` neoAscEl/transAscEl ascertainment correction (`src/node_cl_cache.h:704-712`) — only indirectly.
- `gibbs_spr_impl` / `gibbs_subtree_swap_impl` CLGroup rateScale (`src/mcmc.cpp:1003, 1014, 1043, 1275, 1284, 1308, 1719, 1729, 1757, 1945, 1954, 1978, 5406, 5413, 5435`) — uncovered.
- `validate_swap_partial_cl` (`src/mcmc.cpp:5396, 5406, 5413, 5435`) — uncovered.

**Tautological tests:** lines 25-65 (R-on-R helper identities), lines 194-216 (rate_neo changes likelihood, true pre-fix), lines 270-289 (trans-only invariance is incidental pre-fix and explicit post-fix — would pass on `main`).

**Real regression detectors:** lines 107-144 (R↔C++ parity, but only against itself); lines 153-179 (cache parity, but skipped on CI and asserts only finite-ness); lines 222-262 (extreme rate_neo, but does not isolate trans mechanism).

## Things I checked and they were fine

- `compute_partition_scales` helper algebra: `nNeo · neoScale + nTrans · transScale = nNeo + nTrans` for any `r > 0`, verified at `src/mcmc_state.h:154-170` against the audit document's eq. (C). Degenerate branches return (1.0, 1.0) for nNeo == 0 or nTrans == 0 — correct.
- `cpp_partition_log_likelihood`: every `edgeLen` arg into a pruning / ascertainment / singleton site has been replaced with `scaledEdge` (12 pruning + 4 ascertainment + 4 singleton sites). The single use of bare `edgeLen` left in the function after the patch is the pre-scaling loop itself at `src/mcmc_likelihood.cpp:2353-2354`. Confirmed via grep over the file post-patch: 0 remaining bare-`edgeLen` uses inside the partition body.
- `populate_cache_full` and `populate_cache` both call `compute_partition_scales` and set `unit.rateScale` symmetrically for neo (`pScales.neo`) and trans (`pScales.trans`). Re-derive on re-populate (line 638) is correct — rateNeo may have changed since the previous populate.
- `cache_total_loglik` builds separate `neoAscEl` / `transAscEl` arrays and routes type==0 to neoAscEl, type==1/2 to transAscEl. Correct.
- Hastings ratios on `rate_neo`: Bactrian-scale (case 3) at line 4163 uses `logHastings = log(mult)` — the standard log-scale-MH Jacobian. The Bactrian proposal is on `log(rate_neo)`. Because `compute_partition_scales` is a deterministic function of `rate_neo`, no extra Jacobian is needed when the per-edge effective rates are derived from it (the parameter, not the derived quantity, is the MCMC state). Author's claim verified.
- LogNormal prior on rate_neo (`src/mcmc.cpp:250`, `R::dlnorm(rateNeo, ...)`) untouched — correct, the prior is on the parameter itself.
- Cache invalidation widening at `src/mcmc.cpp:3197-3201` (slice case 3) and `src/mcmc.cpp:4887-4892` (do_move case 3/18) is correct for the *node-CL* cache. The bug in F1 is that the *partition-LL* cache wasn't analogously widened — these are independent caches.
- The `gibbs_partial_cl.h` `CLGroup::rateScale` field is correctly threaded through the F81 / JC apply-transition kernels at lines 366, 474, 495, 543, 715, 845 — the trans rateScale now propagates symmetrically.
- Build / include order: `mcmc_state.h` defines `PartitionScales` and `compute_partition_scales` as a free `static inline` helper before the workspace structs; `node_cl_cache.h` and `mcmc.cpp` and `mcmc_likelihood.cpp` all include `mcmc_state.h` transitively (or directly). No forward-declaration issues observed; the author's claim that `devtools::test()` compiled cleanly is consistent with this.
