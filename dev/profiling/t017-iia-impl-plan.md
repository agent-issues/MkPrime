# T-017 Phase 2a — RNG threading (impl plan)

**Status:** READY-FOR-IMPL (2026-05-22). Plan produced by an Opus read-only Plan subagent on a worktree branched off `worktree-ecology-aware`. To be enacted by a Sonnet impl subagent.

**Base branch for impl:** `worktree-ecology-aware` (has T-017 Phase 1 + 1b + T-018 + T-019 landed).

**Impl branch:** `t017-iia-rng-threading`.

**Scope statement:** Thread a per-chain RNG (`ChainRng`, `std::mt19937_64`-backed) through every R-RNG site reachable from the per-chain dispatch loop body in `run_mcmc_batch_cpp` (mcmc.cpp:6193 onward — the `for (int ch = 0; ch < nChains; ++ch)` block). Code remains **serial** after Phase 2a; same-seed reproducibility within the new build is required; bit-identity vs the pre-Phase-2a serial build is NOT required.

---

## 0. Boundary statement (read this first)

Inside scope (must use `ChainRng`):
- Every R-RNG call inside the per-chain dispatch loop body in `run_mcmc_batch_cpp` (`src/mcmc.cpp:6193-6377` on the ecology branch — covers the move-select, `do_move_impl`, `slice_scalar_impl`, `slice_kprime_hyper_impl` calls plus the int_walk character pick).
- Every R-RNG call inside functions transitively called from that loop body:
  - `do_move_impl` (mcmc.cpp:4779)
  - `slice_scalar_impl` (mcmc.cpp:3249)
  - `slice_kprime_hyper_impl` (mcmc.cpp:3387)
  - `bactrian_perturbation` (mcmc.cpp:191) — 15 in-mcmc call sites + R-export `bactrian_draws` (see §6)
  - `bactrian_2d_perturbation` (mcmc.cpp:212) — 4 in-mcmc call sites + R-export `bactrian_2d_draws`
  - `gibbs_spr_impl{,_full,_het}` (mcmc.cpp:981 / 1273 / 1578)
  - `gibbs_subtree_swap_impl{,_full,_het}` (mcmc.cpp:1751 / 1988 / 2244)
  - `weighted_branch_scale_impl` (mcmc.cpp:2372)
  - `block_gibbs_branch_sweep_impl` (mcmc.cpp:2493)
  - `weighted_spr_impl` (mcmc.cpp:2661)
  - `weighted_subtree_swap_impl` (mcmc.cpp:2942)
  - `pspr_proposal_impl` (mcmc.cpp:3482)
  - `gibbs_kprime_sweep_impl{,_ecology}` (mcmc.cpp:3680 / 3677 / 4193)
  - `block_kprime_shift_impl` (mcmc.cpp:4464)
  - `gibbs_z_sweep_impl` (mcmc.cpp:4624)
  - The two C-style `unif_rand()`/`R::rbeta`/`R::rgamma` calls in `proposals.cpp` and `tree_moves.cpp` that are called from any of the above (full list in §1).
- The two `R::unif_rand()` sites inside the chain loop body itself: move-index pick (mcmc.cpp:6193) and `int_walk` char pick (mcmc.cpp:6213). Note: the swap-proposal sites at :6322/:6328 also live inside the outer iter loop — they run **between** chain-loop iterations (after the per-chain loop, before sample collection). For Phase 2a they could go either way; recommendation: keep them as `R::unif_rand()` since they execute serially and are outside the future parallel region. This stays clean.

Outside scope (keep R RNG):
- Warmup/tuning code that lives at `run_mcmc_batch_cpp` scope but outside the chain loop body (e.g. tuning updates, the initial state init via `init_mcmc_state`, adaptive scheduler updates).
- Any sampling/diagnostic code that runs in the outer-iter scope but after the chain loop closes.
- Pure-density Rmath: `R::dgamma`, `R::dbeta`, `R::dlnorm`, `R::dexp`, `R::pbeta`, `R::qbeta`. These are stateless. They live in `cpp_log_prior` and in `gibbs_partial_cl.h` (lines 1060-1068). See §3.
- `seed_chain_rngs` itself draws `base_seed` from R-RNG once, serially, inside `RNGScope`.

---

## 1. File-by-file change list with LoC budget

### `src/chain_rng.h` (NEW, ~70 LoC)

```cpp
#ifndef MKPRIME_CHAIN_RNG_H
#define MKPRIME_CHAIN_RNG_H

#include <cstdint>
#include <cmath>
#include <random>

// Per-chain RNG. One instance per MCMC chain, owned by run_mcmc_batch_cpp.
// All distribution parameter orderings match Rmath (R::r*) so call-site
// rewrites are mechanical: R::rgamma(shape, scale) -> rng.rgamma(shape, scale).
struct ChainRng {
  std::mt19937_64 gen;

  explicit ChainRng(uint64_t seed) : gen(seed) {}

  // U(0, 1) -- matches R::unif_rand() / unif_rand()
  double unif() {
    // std::generate_canonical for highest-quality double in [0, 1).
    return std::generate_canonical<double, std::numeric_limits<double>::digits>(gen);
  }

  // Normal(mu, sd) -- matches R::rnorm(mu, sd)
  double rnorm(double mu, double sd) {
    std::normal_distribution<double> d(mu, sd);
    return d(gen);
  }

  // Gamma(shape, scale) -- matches R::rgamma(shape, scale).
  // NOTE: std::gamma_distribution(alpha, beta) uses beta as SCALE
  // (same convention as Rmath), confirmed by libstdc++ docs.
  double rgamma(double shape, double scale) {
    std::gamma_distribution<double> d(shape, scale);
    return d(gen);
  }

  // Beta(a, b) -- matches R::rbeta(a, b). std lacks beta; build via
  // two gammas: X ~ Gamma(a, 1), Y ~ Gamma(b, 1), Beta = X / (X + Y).
  double rbeta(double a, double b) {
    double x = rgamma(a, 1.0);
    double y = rgamma(b, 1.0);
    double s = x + y;
    if (s <= 0.0) return 0.5;  // numerical safety; should not occur for valid a,b > 0
    return x / s;
  }
};

#endif // MKPRIME_CHAIN_RNG_H
```

**Decision: match Rmath parameter ordering.** Rationale: every existing call site already uses Rmath's `(shape, scale)`, `(a, b)`, `(mu, sd)` order. Matching it means call-site rewrites are pure textual substitution (`R::rgamma(s, c)` -> `rng.rgamma(s, c)`) with zero risk of arg-swap bugs. std::gamma_distribution's first arg is "alpha" (shape) and second is "beta" — in libstdc++ this "beta" is interpreted as scale, matching Rmath. (Cross-check: cppreference says "beta = scale".)

**API surface is complete.** Restricted-grep of the in-scope functions (§0 list) shows the ONLY R-RNG kinds invoked are: `unif_rand` / `R::unif_rand`, `R::norm_rand` (zero hits inside in-scope code — only `R::rnorm(mu, sd)` is used), `R::rnorm`, `R::rgamma`, `R::rbeta`. No `R::rexp`, `R::rbinom`, `R::rlnorm`, `R::rpois`, `R::rchisq` in scope. Confirm with the regex used in §1 below before adding to ChainRng.

### `src/mcmc.cpp` (~220 LoC modified, ~30 LoC added, net +30)

**Add at the top of file (after existing includes):**
```cpp
#include "chain_rng.h"
```

**`bactrian_perturbation` (mcmc.cpp:191) — keep + overload (see §2):**
```cpp
// Keep existing zero-arg version for the R-exported wrapper bactrian_draws.
static inline double bactrian_perturbation() { ... existing ... }

// New overload for ChainRng-threaded callers (15 in-mcmc sites).
static inline double bactrian_perturbation(ChainRng& rng) {
  double z = rng.rnorm(0.0, BACTRIAN_SD);
  double raw = (rng.unif() < 0.5) ? (BACTRIAN_M + z) : (-BACTRIAN_M + z);
  return raw * BACTRIAN_SCALE;
}
```
Same pattern for `bactrian_2d_perturbation` (mcmc.cpp:212).

**RNG sites in mcmc.cpp that go through ChainRng (line numbers on ecology branch HEAD; expect ±5 drift if rebased):**

The complete list from `git show worktree-ecology-aware:src/mcmc.cpp | grep -nE "R::(unif_rand|norm_rand|rgamma|rbeta|rnorm)" + R-RNG calls in in-scope functions only:

| Line | Site | Function | Replace with |
|---|---|---|---|
| 192 | `R::rnorm(0.0, BACTRIAN_SD)` | bactrian_perturbation(ChainRng&) | `rng.rnorm(0.0, BACTRIAN_SD)` |
| 193 | `R::unif_rand()` | " | `rng.unif()` |
| 215 | `R::rnorm(0.0, 1.0)` | bactrian_2d_perturbation(ChainRng&) | `rng.rnorm(0.0, 1.0)` |
| 216 | `R::rnorm(0.0, 1.0)` | " | `rng.rnorm(0.0, 1.0)` |
| 221 | `R::unif_rand()` | " | `rng.unif()` |
| 223 | `R::unif_rand()` | " | `rng.unif()` |
| 1007 | `R::unif_rand()` | gibbs_spr_impl helper | `rng.unif()` |
| 1225 | `R::unif_rand()` | gibbs_spr_impl_het | `rng.unif()` |
| 1286 | `R::unif_rand()` | gibbs_spr_impl_het | `rng.unif()` |
| 1537 | `R::unif_rand()` | gibbs_spr_impl_het | `rng.unif()` |
| 1589 | `R::unif_rand()` | gibbs_spr_impl_full | `rng.unif()` |
| 1683 | `R::unif_rand()` | gibbs_spr_impl_full | `rng.unif()` |
| 1766 | `R::unif_rand()` | gibbs_subtree_swap_impl | `rng.unif()` |
| 1940 | `R::unif_rand()` | gibbs_subtree_swap_impl_het | `rng.unif()` |
| 1993 | `R::unif_rand()` | gibbs_subtree_swap_impl_het | `rng.unif()` |
| 2203 | `R::unif_rand()` | gibbs_subtree_swap_impl_het | `rng.unif()` |
| 2249 | `R::unif_rand()` | gibbs_subtree_swap_impl_full | `rng.unif()` |
| 2312 | `R::unif_rand()` | gibbs_subtree_swap_impl_full | `rng.unif()` |
| 2384, 2386 | `R::unif_rand()` x2 | weighted_branch_scale_impl | `rng.unif()` |
| 2434 | `R::unif_rand()` | weighted_branch_scale_impl | `rng.unif()` |
| 2449 | `R::rbeta(alphaNew, betaNew)` | weighted_branch_scale_impl | `rng.rbeta(...)` |
| 2464 | `R::dbeta(...)` | weighted_branch_scale_impl | **KEEP** (density, stateless) |
| 2466 | `R::dbeta(...)` | " | **KEEP** |
| 2507, 2530 | `R::unif_rand()` x2 | block_gibbs_branch_sweep_impl | `rng.unif()` |
| 2570, 2584, 2598, 2600, 2614 | mix of unif/rbeta/dbeta | block_gibbs_branch_sweep_impl | unif/rbeta -> rng; dbeta KEEP |
| 2678, 2830, 2842, 2855, 2889, 2891, 2898, 2909 | mix | weighted_spr_impl | unif/rbeta -> rng; dbeta KEEP; cpp_log_prior call KEEP |
| 2951, 3040, 3053, 3066, 3098, 3100, 3107, 3118 | mix | weighted_subtree_swap_impl | unif/rbeta -> rng; dbeta KEEP; cpp_log_prior KEEP |
| 3166, 3261, 3265, 3285, 3291 | unif + cpp_log_prior | slice_scalar_impl (`logZ = logY0 + log(unif)`; bounds; cpp_log_prior calls) | unif -> rng; cpp_log_prior KEEP |
| 3399, 3408, 3422, 3438, 3446 | similar | slice_scalar_impl | unif -> rng; cpp_log_prior KEEP |
| 3505, 3598, 3608 | unif | pspr_proposal_impl | rng.unif() |
| 4059, 4079 | unif x2 | gibbs_kprime_sweep_impl | rng.unif() |
| 4157, 4367, 4387, 4443, 4451, 4489, 4569 | mix | gibbs_kprime_sweep_impl_ecology + block_kprime_shift_impl | unif -> rng; cpp_log_prior KEEP |
| 4707 | unif | gibbs_z_sweep_impl | rng.unif() |
| 4744 | cpp_log_prior | " | KEEP |
| 4779-5982 | **`do_move_impl` body** — must take `ChainRng& rng` arg | (see §2 for signature) | every bactrian_perturbation() call becomes `bactrian_perturbation(rng)`; every R::unif_rand becomes `rng.unif()` (sites at 4905, 4910, 4911, 4917, 4923, 4939, 4955, 4957, 5040, 5068, 5069, 5122, 5128, 5148, 5158, 5189, 5197, 5208, 5217, 5240, 5242, 5263, 5276, 5304, 5318, 5343, 5370) |
| 5830 | unif | (in cpp_log_likelihood diagnostic, mcmc.cpp body — verify if reachable from hot path; treat as out-of-scope if it's in a diagnostic-only function) | INVESTIGATE; likely KEEP |
| **6193** | `R::unif_rand() * tw` (move-index select) | **chain loop body** | `rngs[ch].unif() * tw` |
| **6213** | `R::unif_rand() * nTrans` (int_walk char pick) | **chain loop body** | `rngs[ch].unif() * nTrans` |
| 6322, 6328 | unif (swap proposals) | OUTSIDE chain loop body, serial | **KEEP** (serial; documented in §0) |

**LoC impact in mcmc.cpp:** ~110 line edits + 30 new lines for the `ChainRng& rng` arg threading through 17 function signatures. Internal-checkpoint candidate after the first 3 functions (ChainRng + bactrian + slice_scalar_impl).

### `src/proposals.cpp` (~12 LoC modified)

RNG sites:
- Line 44: `unif_rand()` (bare C form) — function = `mk_split_subtree_proposal` (called from tree-move impls in scope)
- Line 115: `unif_rand()` — same
- Line 121: `unif_rand()` — same
- Line 211: `unif_rand()` — likely `dirichlet_simplex_impl` (called from `do_move_impl`)
- Line 227: `R::rbeta(alpha, betaPar)` — same
- Line 300: `R::rgamma(alphaFwd[i], 1.0)` — same
- Line 355: `unif_rand()` — `local_dirichlet_impl` (?)
- Line 394: `unif_rand()` — `beta_simplex_impl` start-edge pick

Impl agent action: read each function's signature, add `ChainRng& rng`, replace `unif_rand()` → `rng.unif()`, `R::rbeta(a,b)` → `rng.rbeta(a,b)`, `R::rgamma(s,c)` → `rng.rgamma(s,c)`. Update declarations in mcmc.cpp lines 158-168 (the forward decls for `beta_simplex_impl`, `dirichlet_simplex_impl`, `local_dirichlet_impl`) and the 3-4 in-mcmc call sites to pass `rng`.

### `src/tree_moves.cpp` (~12 LoC modified)

RNG sites (all `unif_rand()` bare C form):
- Line 45: `pickInternal` (in `nni_proposal_impl`)
- Line 70, 72: `pickV`, `pickU` (in `nni_proposal_impl`)
- Line 281: `pickPrune` (in `spr_proposal_impl` or `tbr_proposal_impl`)
- Line 335: `pickSub`
- Line 382: `sigma = unif_rand()`
- Line 425: `(void)unif_rand()` (RNG order-preserving no-op draw)
- Line 455: `pickRegraft`
- Line 460: `tau = unif_rand()`

These functions are called from `do_move_impl`. Forward decls live at mcmc.cpp lines 141-155 (`nni_proposal_impl`, `spr_proposal_impl`, `swap_subtrees_impl`, `tbr_proposal_impl`). All four take `ChainRng& rng` after Phase 2a.

### `src/mcmc_state.h` (NO CHANGE in Phase 2a)

Phase 2b will add per-chain workspace fields. Phase 2a leaves the struct untouched.

### `src/RcppExports.cpp` (auto-regen)

If any [[Rcpp::export]] function signatures change. Per §6 we keep `bactrian_draws` and `bactrian_2d_draws` signatures unchanged, so RcppExports likely needs no regeneration. Run `Rcpp::compileAttributes()` and check the diff.

### Tests (see §4 for the full strategy)

Updates likely in 8-18 files. Pre-flagged list in §4.

### `dev/profiling/findings.md` (~5 LoC)

T-017-IIa row.

### `dev/profiling/log.md` (~40 LoC)

Round entry.

**Total Phase 2a LoC budget:** ~330 net (~370 added / ~40 removed). Internal checkpoint at 250 LoC = commit `ChainRng` header + `bactrian_*` overloads + `do_move_impl` signature + ~3 in-scope helper signatures.

---

## 2. Signature-refactor decision

For `do_move_impl`, `slice_scalar_impl`, `slice_kprime_hyper_impl`, and all the `gibbs_*_impl` / `weighted_*_impl` / `block_*_impl` / `pspr_proposal_impl` / `gibbs_z_sweep_impl` helpers:

**Recommendation: modify signatures in place (no overload).** Reasons:
- All are `static` (file-local linkage) — no external callers, no ABI concern.
- All have exactly ONE call site each (from `do_move_impl` body or `run_mcmc_batch_cpp`).
- An overload would require maintaining two divergent code paths forever; we'd never delete the old one.
- The diff is mechanical and lossless — `data, state, ...` becomes `data, state, rng, ...`.

For `bactrian_perturbation` and `bactrian_2d_perturbation`:

**Recommendation: KEEP the existing zero-arg version AND add a ChainRng-arg overload.** Reasons:
- `bactrian_draws` and `bactrian_2d_draws` are R-exported (mcmc.cpp:199, 235) and used by `tests/testthat/test-bactrian.R` and `test-joint-2d.R` to validate the kernel's statistical properties under R's seed.
- Forcing those tests to construct a ChainRng would either (a) require exporting ChainRng to R (unjustified scope creep) or (b) fail the tests for an unrelated reason.
- 15 + 4 = 19 in-mcmc call sites switch to the new overload `bactrian_perturbation(rng)`; the two R-export wrappers keep calling `bactrian_perturbation()`.

Cost: ~10 extra LoC (two short overloads). Saves test churn on 2 well-tested files.

---

## 3. `cpp_log_prior` Rmath handling

**Recommendation: Option A (leave `cpp_log_prior` untouched in Phase 2a).**

Rationale:
- Phase 2a is still serial. `R::dgamma`, `R::dbeta`, `R::dlnorm`, `R::dexp` are stateless under serial execution. They can `Rf_error` longjmp on out-of-domain, but that's already the case in the pre-2a code and no current test exercises it.
- The `Rf_error` UB-from-non-master-thread concern only materialises when chain dispatch goes parallel (Phase 2b). Phase 2b's plan will handle it then (likely by pre-validating inputs and returning `R_NegInf` directly).
- Pre-emptive input validation in Phase 2a is non-trivial: there are 12 distinct density calls across `cpp_log_prior` (mcmc.cpp:393-595), each with its own valid-domain check. Doing it now risks subtle prior-evaluation behaviour changes that are hard to verify without comparing against a known-good baseline.
- The user's explicit preference (per the brief): **minimal scope for 2a; defer Rmath-validation work to 2b unless trivial**. This is not trivial.

Defer to Phase 2b's plan, which should add a §"Rmath validation in parallel region" item.

---

## 4. Tests likely needing expected-value updates

### Strategy (operational rule for impl agent)

Don't try to triage 67 set.seed-using files up front. Instead:
1. After all RNG threading is in place, build with `devtools::load_all()`.
2. Run the full suite once: `devtools::test()`.
3. For each FAIL: open the test, look at the failing assertion.
   - If it's `expect_equal(x, <hard-coded-numeric>)` after a `set.seed()` + MCMC trace → **pin-the-seed**; record old and new value, update.
   - If it's a property test (Rhat < threshold, ESS > threshold, density-match within tolerance, dimension check, S3 class check) that's failing → **semantic regression**; STOP and diagnose.
4. Budget: ≤18 test updates expected. If >18 expected-value updates seem needed, the RNG threading probably has a bug — pause and re-grep for missed sites.

### Pre-flagged likely-affected files (PIN-THE-SEED candidates)

Files that run MCMC traces with pinned expected numerics — read these first if needed:

| File | Why it's likely affected |
|---|---|
| `tests/testthat/test-mcmc-engine.R` | Runs MCMC, checks `nrow(result$samples) == 60` (safe, dim-only) but also `acceptance` values per move — those may be pinned numerically; verify |
| `tests/testthat/test-mcmc-recovery.R` | 8 set.seed sites; recovery-after-checkpoint relies on RNG-stream stability |
| `tests/testthat/test-checkpoint.R` | 20+ set.seed sites; checkpoint round-trip may compare sample values |
| `tests/testthat/test-tempering.R` | PT swap acceptance counts may be pinned |
| `tests/testthat/test-three-phase.R` | Likely runs an MCMC with pinned summary stats |
| `tests/testthat/test-stepping-stone.R` | Marginal-likelihood estimates depend on RNG order |
| `tests/testthat/test-streaming.R` | May check stream content exactly |
| `tests/testthat/test-proposals.R` | Direct proposal tests; check whether they hit beta_simplex_impl / dirichlet_simplex_impl with pinned expected values |
| `tests/testthat/test-tree-moves.R` | Similar to proposals |
| `tests/testthat/test-gibbs-spr.R`, `test-gibbs-batched.R`, `test-gibbs-kprime.R`, `test-gibbs.R`, `test-gibbs-het.R` | Gibbs paths use RNG; trace tests likely pinned |
| `tests/testthat/test-weighted-spr.R`, `test-weighted-branch-scale.R`, `test-weighted-subtree-swap.R` | weighted_*_impl funcs are in scope |
| `tests/testthat/test-block-gibbs-branch.R` | block_gibbs_branch_sweep_impl is in scope |
| `tests/testthat/test-pspr.R` | pspr_proposal_impl in scope |
| `tests/testthat/test-mcmc-engine.R`, `test-burnin.R`, `test-adapt-thin.R`, `test-auto-thin.R` | Engine-level integration |
| `tests/testthat/test-cache-aware-scheduling.R`, `test-m092-adaptive-scheduler.R` | Scheduler decisions influenced by RNG-driven move selection |
| `tests/testthat/test-tbr.R`, `test-tbr-detailed-balance.R`, `test-subtree-swap.R`, `test-nni-preorder.R` | Tree-move impls hit RNG |
| `tests/testthat/test-bg-hyperparameter-moves.R`, `test-empirical-geometric-prior.R`, `test-beta-geometric-prior.R`, `test-logseries-prior.R` | kPrime hyper moves use slice_kprime_hyper_impl + gibbs_kprime |
| `tests/testthat/test-m125-dirichlet-branch.R`, `test-m127-partial-eval-dirichlet.R`, `test-m126-warmup-rho.R`, `test-m140-m141.R`, `test-m151-streaming-multirun.R`, `test-m161-selective-cache.R`, `test-m170-csp-symmetry.R`, `test-m177-warmup-eta.R`, `test-m090-move-wiring.R` | Specific-move regression tests, likely pinned to RNG order |

### Pre-flagged confirmed-safe (SEMANTIC tests — should still pass)

| File | Why it's safe |
|---|---|
| `tests/testthat/test-omp-determinism.R` | Tests `cpp_log_likelihood_ecology` only; no RNG path (see §7) |
| `tests/testthat/test-bactrian.R` | Tests `bactrian_draws()` (R-export, uses R's RNG via the zero-arg overload kept per §2) |
| `tests/testthat/test-joint-2d.R` | Same — `bactrian_2d_draws()` keeps R-RNG path |
| `tests/testthat/test-fast-exp.R`, `test-rate-matrix.R`, `test-acrv.R` | No MCMC, pure-math |
| `tests/testthat/test-likelihood.R`, `test-full-loglik.R`, `test-jc-collapse.R`, `test-jc-collapse-hotpath.R` | Likelihood-only, no RNG-driven trace |
| `tests/testthat/test-priors.R`, `test-ecology-prior.R`, `test-validation.R`, `test-ecology-likelihood.R` | Density-eval, no RNG path |
| `tests/testthat/test-ess-rhat.R`, `test-treeESS.R`, `test-adaptive-tree-ess.R` | Property tests; pass on any reasonable trace |
| `tests/testthat/test-convergence.R`, `test-convergence-kprime.R` | Convergence diagnostics on trace properties |
| `tests/testthat/test-MkPrimeData.R`, `test-partition.R`, `test-fused-ascertainment.R`, `test-ascertainment.R`, `test-corrections.R` | Data-prep / pure-math |
| `tests/testthat/test-ticker-helpers.R`, `test-progress.R`, `test-progress-table.R`, `test-smoke.R`, `test-edge-cases.R`, `test-runmkprime-dots.R`, `test-bayesian-module.R`, `test-stopping.R` | Plumbing / UX tests, no RNG-pin |
| `tests/testthat/test-relabel-ecology.R`, `test-ecology-scaffold.R`, `test-ecology-plumbing.R`, `test-ecology-logging.R`, `test-ecology-smoke.R` | Plumbing |
| `tests/testthat/test-swap-partial-cl.R`, `test-node-cl-cache.R` | Cache validity, no RNG-pin |
| `tests/testthat/test-tree-thin.R` | Subsetting logic |
| `tests/testthat/test-posterior-multi.R`, `test-parallel.R`, `test-parallel-resume.R`, `test-interrupt-recovery.R` | Mostly orchestration; some traces may be pinned — verify |
| `tests/testthat/test-t018-adaptive-moves.R` | Recently added; likely property-based |

**New test to add: `tests/testthat/test-chain-rng-determinism.R`** (~30 LoC, see §8 pass criteria).

---

## 5. Seed-stream design

```cpp
// Inside run_mcmc_batch_cpp, BEFORE the outer iter loop and BEFORE any
// future parallel region. RNGScope is already active here (the function
// is called via [[Rcpp::export]]).
uint64_t base_seed = static_cast<uint64_t>(
  R::unif_rand() * static_cast<double>(std::numeric_limits<uint32_t>::max())
);

std::vector<ChainRng> rngs;
rngs.reserve(nChains);
// Weyl-mix constant: 2^64 / phi, where phi is the golden ratio.
// Standard mixing constant (used by SplitMix64, etc.) — spreads chain
// seeds away from base_seed + ch which can be too close at the low bits.
constexpr uint64_t WEYL_MIX = 0x9E3779B97F4A7C15ULL;
for (int ch = 0; ch < nChains; ++ch) {
  uint64_t chain_seed = base_seed ^ (static_cast<uint64_t>(ch) * WEYL_MIX);
  rngs.emplace_back(chain_seed);
}
```

Rationale:
- One `R::unif_rand()` draw consumes from R's seed stream, so `set.seed()` in R correctly determines the whole MCMC trace.
- Casting `unif()` to `uint32_t` range gives 2^32 distinct base seeds — plenty for distinguishing different R-seed states; we don't need 2^64 because R's RNG state itself is 2^32-bounded for the Mersenne Twister.
- The Weyl mix (the golden-ratio constant) breaks low-bit correlation between adjacent chain seeds. For nChains=4 a naive `base_seed + ch` would have 30 of 32 bits identical across chains — pathological for some PRNGs and a bad smell even for `std::mt19937_64`. Weyl mixing makes the chain seeds maximally far apart.
- No XOR with `base_seed = 0` collision: even if `R::unif_rand() = 0` exactly (effectively never), `WEYL_MIX * ch` gives distinct nonzero seeds for ch >= 1; chain 0 would seed with 0 but mt19937_64 handles that fine.

Note for impl agent: declare `WEYL_MIX` as a `static constexpr` at file scope or inside `run_mcmc_batch_cpp`; either is fine.

---

## 6. Build matrix & smoke-test commands

The impl agent should run these BEFORE pushing:

```r
# 1. Compile-only check.
devtools::load_all(".")

# 2. Targeted smoke (fast, focuses on RNG-threaded paths).
devtools::test(filter = "mcmc-engine|gibbs|weighted|tree-moves|proposals|bactrian|joint-2d|chain-rng-determinism|omp-determinism")

# 3. Full suite (slow; ~5-10 min).
devtools::test()
```

**Bench (must run before push):**

```sh
# In a bash shell, from the worktree root:
MKP_TIMING_NCHAINS=4 OMP_NUM_THREADS=1 Rscript dev/profiling/drivers/rodent_aware_timing.R
```

Expected wall: **174 s ± 9 s** (±5% of the pre-2a baseline). If wall is >185 s, investigate (ChainRng overhead vs Rmath should be negligible, but std::gamma_distribution can be slower than Rmath::rgamma in some libstdc++ versions — measure).

5 timed runs, report median.

---

## 7. Bit-identity test impact prediction

`tests/testthat/test-omp-determinism.R` (3 test_that blocks, 5 assertions): **all pass under Phase 2a.**

Reasoning:
- The test exercises `MkPrime:::.CppLogLikelihoodEcology` — a pure C++ log-likelihood function with no RNG path.
- Phase 2a does NOT change `cpp_log_likelihood_ecology`, `cpp_partition_log_likelihood_ecology`, the inner pruners, or any code reachable from `.CppLogLikelihoodEcology`.
- The test's `set.seed(seedZ)` is used to generate a random z-matrix BEFORE the call — it does NOT seed a ChainRng (which doesn't exist for likelihood-eval paths).
- OMP=1 vs OMP=4 bit-identity is determined entirely by the partial-sum merge order, which Phase 2a doesn't touch.

If this test fails after Phase 2a, it indicates a code path inadvertently rewired through ChainRng that shouldn't have been — a sign of scope creep.

---

## 8. Pass criteria

Phase 2a is **DONE** when ALL of the following hold:

1. **Coverage:** Every R-RNG call inside any function listed in §0 "Inside scope" goes through `ChainRng`. Verified by:
   ```sh
   git diff worktree-ecology-aware -- src/mcmc.cpp src/proposals.cpp src/tree_moves.cpp \
     | grep -E "^-\s*(R::(unif_rand|rnorm|rgamma|rbeta)|\bunif_rand\(\))"
   ```
   should show ~85+ deletions, balanced by ~85+ `rng.*` additions.

   Negative check:
   ```sh
   git show HEAD:src/mcmc.cpp \
     | awk '/^static bool do_move_impl/,/^}/' \
     | grep -E "R::(unif_rand|rnorm|rgamma|rbeta)\(|\bunif_rand\(\)"
   ```
   should return **zero** lines (no in-scope R-RNG sites remaining in `do_move_impl`).

2. **Test suite:** ≥ 228 / 246 tests pass. Up to 18 test files may need expected-value updates (per §4 budget). 6 expected skips persist.

3. **Determinism:** A new test in `tests/testthat/test-chain-rng-determinism.R` proves same-seed reproducibility within the new build:
   ```r
   test_that("Same R seed produces identical MCMC samples (T-017-IIa)", {
     mkd <- MkPrimeData(...small fixture...)
     set.seed(20260522)
     r1 <- run_mcmc(mkd, nChains = 4, nBatch = 200, ...)
     set.seed(20260522)
     r2 <- run_mcmc(mkd, nChains = 4, nBatch = 200, ...)
     expect_identical(r1$samples, r2$samples)
   })
   ```

4. **Bench:** Rodent aware nChains=4 OMP=1: **≤ 183 s** wall (174 s ± 5%). 5 runs, median. No regression.

5. **No silent skips:** the `omp-determinism` test still passes (§7).

---

## 9. Anti-pattern reminders for impl agent

These are not new — they're the lessons from Round 13/14 silent-death incidents. **Read these and internalise.**

- **Commit early.** As soon as `chain_rng.h` compiles and `bactrian_perturbation(ChainRng&)` plus the first 2-3 in-scope helpers (`slice_scalar_impl` is a small, contained target — do it first) are wired and the suite still builds, **commit**. This is your salvage point. If the session is killed at 70% through `do_move_impl`, the dispatcher can salvage the working part and dispatch a follow-up for the rest.

- **Internal checkpoint at 250 LoC.** Track diff size with `git diff --stat`. At 250 LoC added, stop, commit (even if intermediate-state), and proceed.

- **Don't preserve bit-identity vs serial pre-2a.** The bar is **same-seed reproducibility within the new build**, not "same trace as the old build". Trying to match the old build's RNG order would require replicating Rmath's internal Marsaglia-polar / Wichmann-Hill specifics — not feasible.

- **Don't write a Write/Edit before reading.** When you find a function signature in mcmc.cpp, READ ~30 lines of context first to see what fields of `state` it reads, what local RNG calls it makes, what helper functions it calls. Then edit. The do_move_impl function is 1200+ LoC — read it once start to end before modifying.

- **Forward decls.** mcmc.cpp lines 141-168 declare proposal/simplex helpers. When you change those impls' signatures (in proposals.cpp, tree_moves.cpp), you MUST also change the forward decl in mcmc.cpp. Compile errors at that point will be clear; just don't forget.

- **R-export wrappers.** `bactrian_draws` and `bactrian_2d_draws` at mcmc.cpp:199 and :234 are [[Rcpp::export]] — they must keep the same R-visible signature. Per §2, they keep calling `bactrian_perturbation()` (no arg). DO NOT thread a ChainRng through them.

- **The cpp_log_prior temptation.** Don't refactor `cpp_log_prior` in Phase 2a (see §3). Even small "while I'm here" edits to it expand scope and risk subtle prior-eval changes. Phase 2b's plan will address it.

- **Watch for `unif_rand()` (no `R::`).** proposals.cpp and tree_moves.cpp use the bare C form. A regex of `R::unif_rand` will miss them.

- **Bench BEFORE committing your final state.** If the bench regresses >5%, debug before push.

- **30-min status check expectation:** if you've been working for 30 min without a commit on `t017-iia-rng-threading`, the dispatcher will SendMessage to check. Always respond with "still working, here's diff stat" or "stuck on X". Silent death is the failure mode the user is most worried about — explicit liveness signals prevent it.
