# T-017 Phase 2 Design — Threaded per-chain dispatch

**Status:** PENDING (2026-05-22). Self-contained brief for a Plan-then-action
chip / subagent sequence to enact post-compact. Reads cold from this doc
+ cross-references.

**Branch base:** `worktree-ecology-aware` (has T-014/T-015/T-017 Phase 1+1b/T-018).

---

## Headline

Phase 1 + Phase 1b landed OpenMP infrastructure at the **partition pruner**
level, but rodent's workload (2 partitions, 50/50 neo-vs-trans imbalance)
caps parallel speedup at 2× and OMP overhead per call cancels it. **The
only remaining code-only path to ≥2× wall on rodent is parallel CHAIN
dispatch** — parallelise the `for (int ch = 0; ch < nChains; ++ch)` loop
in `src/mcmc.cpp:6129` so that 4 chains run on 4 cores simultaneously,
each chain serial internally.

**Projected: 3-3.5× wall at nChains=4** (4-way chain parallel, bounded
by serial swap + sample collection ~8 ms / 370 ms outer iter, per T-015
chrono). Hits Hamilton viability target.

**Estimated LoC:** 700-1100 net. **RNG order will change**, breaking many
test expected values.

**This is the high-risk follow-up.** Plan-then-action sequence + advisor
gating + Sonnet impl with explicit early-status-check is the proven
pattern (Phase 1 used this — see Round 14 in `log.md`).

---

## Why Phase 2 is harder than Phase 1

Phase 1 parallelised inside `cpp_log_likelihood_ecology` — a pure-C++
region with std::vector workspaces. Safe.

Phase 2 parallelises the per-chain move dispatch loop, which calls
`do_move_impl` / `slice_scalar_impl` / `slice_kprime_hyper_impl`.
These call:

- **47 R-RNG sites** (`R::unif_rand`, `R::norm_rand`, `R::rgamma`,
  `R::rbeta`, etc.) — mutate `.Random.seed`, **UNSAFE in parallel**.
- **74 Rcpp alloc sites** (`NumericVector(n)`, `IntegerMatrix(r,c)`,
  `clone(...)`) — touch R's memory manager + PROTECT stack, **UNSAFE
  in parallel**.
- Various `R::dgamma`/`R::dbeta`/`R::dlnorm` density calls — pure math,
  probably safe, but call `Rf_error` on out-of-domain → longjmp UB from
  non-master thread.

Full inventory in the Phase 1 design doc (`t017-design.md`, Section 1).
That inventory is fresh and saves the Plan agent a re-scan.

---

## Implementation strategy (recommended split)

### Phase 2a: RNG threading only (~300 LoC, drop-in-compatible)

- Create `src/chain_rng.h` defining `ChainRng` struct:
  ```cpp
  struct ChainRng {
    std::mt19937_64 gen;
    explicit ChainRng(uint64_t seed) : gen(seed) {}
    double unif() { /* ... */ }
    double rnorm(double mu, double sd) { /* ... */ }
    double rgamma(double shape, double scale) { /* ... */ }
    double rbeta(double a, double b);  // X / (X+Y) via two gammas
    double rexp(double rate) { /* ... */ }
    int    rbinom(int n, double p) { /* ... */ }
  };
  ```
- Thread `ChainRng&` through `do_move_impl` and all its helpers (47
  RNG sites). Most sites are inside `bactrian_perturbation` /
  `bactrian_2d_perturbation` (defined at `mcmc.cpp:183, 204`) — those
  helpers themselves call `R::rnorm`. Thread `ChainRng&` through them.
- In `run_mcmc_batch_cpp`, draw `base_seed` once via `R::unif_rand`
  (under the existing RNGScope, before any parallel region), then
  construct `std::vector<ChainRng> rngs; for (int ch = 0; ch < nChains; ++ch)
  rngs.emplace_back(base_seed + ch);`.
- Replace `R::unif_rand()` etc. with `chainRng.unif()` etc. at every
  threaded site.

**Phase 2a alone does NOT enable parallelism** — chain dispatch stays
serial. But it lays the groundwork without breaking anything visible.

**Bit-identity is BROKEN** vs. pre-2a serial code. Many tests have
hard-coded expected values that depend on the specific RNG draw order.
Those tests will need expected-value updates.

**Test impact** (from Phase 1 design doc Section 5):
- Pin-the-seed tests in `test-tempering.R`, `test-mcmc-engine.R`,
  `test-checkpoint.R`, `test-mcmc-recovery.R`, `test-gibbs-kprime.R` etc.
  All need expected-value updates.
- Genuine semantic tests (`test-convergence.R`, `test-validation.R`,
  `test-likelihood.R`, etc.) should still pass.
- ~10-15 test files affected based on `set.seed` pattern grep.

### Phase 2b: Alloc hoisting + parallel chain dispatch (~400-700 LoC)

Once 2a's foundation is in place:

- Move per-call Rcpp allocations in `do_move_impl` into per-chain
  workspaces on `McmcState`. ~15 new fields: `NumericVector edgeLenBuf`,
  `IntegerVector workParBuf`, etc.
- Add `#pragma omp parallel for schedule(static, 1) num_threads(min(nChains, omp_get_max_threads()))`
  around the inner per-chain loop at `src/mcmc.cpp:6129`.
- Each thread reads `states[ch]` and `rngs[ch]` — no contention because
  per-chain.
- `state->ecoCL` and `state->nodeCL` caches are per-chain (no contention).
- Chain swap at :6243 stays SERIAL (post-parallel-region).
- Sample collection at :6250 stays SERIAL.

---

## Critical risks (from advisor's pre-flight analysis)

1. **R-API in supposedly-safe paths.** `cpp_log_prior` calls
   `R::dbeta`/`R::dlnorm`/`R::dgamma` — pure math but `Rf_error` on
   out-of-domain → longjmp UB from non-master thread. Mitigation:
   pre-validate inputs; return `R_NegInf` directly without calling
   Rmath on out-of-domain args. Or: keep `cpp_log_prior` serial,
   parallelise only the likelihood call.

2. **Rcpp PROTECT stack corruption.** Concurrent Rcpp object
   construction from parallel threads corrupts R's PROTECT stack.
   **Strict audit needed**: every line under the chain-dispatch
   `#pragma omp parallel for` must be Rcpp-allocator-free. Phase 2b's
   alloc hoisting is the mitigation — move all Rcpp construction
   outside the parallel region.

3. **Determinism break under thread schedule.** Per-chain `ChainRng`
   keyed by `base_seed + chain_index` is deterministic regardless of
   thread schedule. As long as no shared state is written from parallel
   threads, results are reproducible (not bit-identical to serial code,
   but same-seed→same-result).

4. **glibc malloc contention.** Same risk as Phase 1 — but Phase 2b
   eliminates per-call Rcpp alloc by hoisting workspaces to McmcState.
   Inner pruner's std::vector workspaces (in `cpp_log_likelihood_ecology`)
   are already per-thread via Phase 1b. Should be OK.

5. **OpenMP not on Hamilton.** Hamilton gcc 9+ has libgomp by default.
   Verify pre-deploy via `module list && which gcc && gcc -fopenmp test.c`.

6. **Tests that depend on RNG order.** Many — see list above. These
   need expected-value updates. NOT a sign of regression if the
   updated value still passes the same property check (Rhat, ESS, etc.).

7. **OMP overhead per outer iter.** Spawning a parallel region per
   outer iter (×500 iter × nBatches) — overhead is bounded (~10-100 µs
   per region, total ~5-50 ms across a 500-iter run). Negligible vs.
   the ~370 ms outer-iter cost.

---

## Dispatch strategy (Plan-then-action, proven pattern)

This is the proven pattern from Round 14. Two-stage chip:

### Stage 1: Plan subagent (Opus, worktree, read-only)

Brief: "Read `dev/profiling/t017-ii-design.md` + the Phase 1 design at
`dev/profiling/t017-design.md` (Section 1 has the R-API inventory).
Produce a concrete implementation plan for Phase 2a (RNG threading only),
including:
- Exact list of files to modify and LoC budget per file
- Exact list of test files needing expected-value updates (grep for
  `set.seed` + `expect_equal` patterns)
- Whether to refactor `bactrian_perturbation` helpers' signatures or
  add a new `bactrian_perturbation_rng(ChainRng&)` overload
- Whether to handle `cpp_log_prior` Rmath calls in the parallel region
  or keep serial
- A specific bench-target methodology
Write design to `dev/profiling/t017-ii-impl-plan.md`. Do not implement —
Plan agents are read-only."

### Advisor pass on plan

Open the plan + call `advisor()`. Look for:
- Hidden RNG sites the Plan agent missed (sub-helpers, internal funcs)
- Test-list completeness
- Phase 2a-vs-2b sequencing (2a should be drop-in-compatible; 2b is
  the parallel-dispatch landing)

### Stage 2: Impl subagent (Sonnet, worktree, foreground with 30-min status check)

Brief: "Implement Phase 2a per `dev/profiling/t017-ii-impl-plan.md`.
Branch off `worktree-ecology-aware` as `t017-ii-rng-threading`. Bit-
identity vs serial pre-T-017-II is NOT required (RNG order changes);
same-seed reproducibility within the new build IS required. Update
test expected values for the files identified in the plan. Push when
228+/246 tests pass (allow up to 10 tests to need expected-value
updates — track which). Bench: rodent_aware_timing.R nChains=4
OMP=1 (still serial in Phase 2a) — should match pre-2a baseline
within ±5%. If after 30 min there's no commit, ping with SendMessage."

### Stage 3 (after 2a lands): Phase 2b dispatch

Once 2a is verified, file as APPLIED, then dispatch a second chip for
Phase 2b (alloc hoisting + parallel chain dispatch + OMP wiring).
Phase 2b brief will be in `dev/profiling/t017-iib-impl-plan.md`
(written by an Opus Plan agent in Stage 3a, mirroring this structure).

---

## Bench target

`dev/profiling/drivers/rodent_aware_timing.R MKP_TIMING_NCHAINS=4 OMP_NUM_THREADS=4`.

- Baseline (current worktree-ecology-aware, Phase 1+1b landed):
  174 s wall (2.9 iter/s) at OMP=1, ~191 s at OMP=4 (Phase 1+1b made
  things worse on rodent — OMP overhead without gain).
- Phase 2a target: ≤ 180 s OMP=1 (no regression — chain dispatch
  still serial, only RNG path changed).
- Phase 2b target: **≤ 60 s wall at nChains=4 OMP=4** (≥ 3× improvement
  vs. Phase 1+1b OMP=1 baseline). This is the brief's ≥ 2× target,
  comfortably met.

Hamilton viability after Phase 2b: at 4× wall reduction, the 30 k iter
aware warmup that currently projects 12 h on Hamilton would drop to ~3 h,
clearly under the v2 8 h walltime limit. **Phase 2 is the key to
production aware on Hamilton.**

---

## File references

- `src/mcmc.cpp:6005-6249` — `run_mcmc_batch_cpp`, the PT loop, where
  `#pragma omp parallel for` over chains goes
- `src/mcmc.cpp:6129` — the `for (int ch = 0; ch < nChains; ++ch)` loop
- `src/mcmc.cpp:4751-5982` — `do_move_impl` (1232 LoC, target of RNG
  threading + alloc hoisting)
- `src/mcmc.cpp:3241-3450` — `slice_scalar_impl`, `slice_kprime_hyper_impl`
- `src/mcmc.cpp:183, 204` — `bactrian_perturbation`, `bactrian_2d_perturbation`
  (called from 23 sites in do_move_impl, internally call `R::rnorm`)
- `src/proposals.cpp` — 2 RNG sites + 12 alloc sites
- `src/tree_moves.cpp` — 0 RNG, 21 alloc sites
- `src/mcmc_state.h` — McmcState; Phase 2b adds per-chain workspaces here
- `dev/profiling/t017-design.md` — Phase 1 design (has full R-API
  inventory in Section 1)
- `dev/profiling/log.md` Round 14 — Phase 1 round notes including
  agent silent-death pattern and salvage lessons

---

## Don't-do list

- Don't attempt Phase 2 in one shot. Split into 2a (RNG only) then 2b
  (alloc hoisting + parallel). Lessons from Round 14: 950 LoC overshoot
  on a 430-LoC budget caused agent death.
- Don't try to preserve bit-identity vs. pre-2a serial code. Document
  the break and update expected values.
- Don't run RNG-stateful R API in parallel region — use ChainRng.
- Don't run Rcpp constructors in parallel region — hoist to McmcState
  (Phase 2b).
- Don't run `Rf_error` / `Rcpp::stop` paths in parallel region — replace
  with `throw std::runtime_error` caught at the parallel-region boundary,
  or pre-validate inputs.
- Don't skip the advisor gate between Plan and Impl stages.
- Don't dispatch as a single Opus agent without the Plan-then-action
  scaffolding — silent death is the known failure mode for big-refactor
  briefs without scaffolding.
