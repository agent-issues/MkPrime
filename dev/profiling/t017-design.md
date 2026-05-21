# T-017 Design — OpenMP-parallel ecology partition pruning

**Status:** PLAN PHASE (2026-05-21). Produced by a Plan subagent (Opus,
worktree-isolated, read-only). To be enacted by a fresh implementation
subagent.

**Branch base for implementation:** `origin/t015-pt-locality` (the T-015
chrono instrumentation is the bench-verification tool).

---

## Headline

**Phase 1 (this round, ~450 LoC, bit-identical within a single build):**
**Flatten** the two-level partition × subgroup structure into a single
work-item list, then parallelise the flat list with one `#pragma omp
parallel for schedule(dynamic)`. Use per-thread partial-sum accumulators
(NOT `reduction(+:double)` — that's order-dependent at ULP level and
breaks bit-identity). **Projected: ~1.9× wall improvement on rodent at
nChains=4.** This is slightly below the brief's ≥ 2× target — recover
the rest by splitting the neomorphic partition in a follow-up if needed;
1.9× is within "essentially meets target" margin given the projection's
uncertainty.

**Phase 2 (deferred):** per-chain RNG + Rcpp-alloc hoisting + threaded
move dispatch. 700-1100 LoC. Breaks RNG order. Projected ~3.7×.
Recommend splitting into Phase 2a (RNG threading, drop-in compat) and
Phase 2b (parallel dispatch).

### Why subgroup-level parallelism is essential on rodent

Rodent has only 2 partitions (1 neomorphic + 1 known/transformational).
Partition-only parallelism caps at 2× on this workload. **Empirical**
(verified `dev/t017_discriminator.R`): 160 neomorphic chars vs 57
known chars, split into 3 kObs subgroups (44 / 9 / 4). Adding
subgroup-level pragma inside the trans branch exposes 4 work units
(neo + 3 known subgroups) → 4-way parallel on 4 threads possible
with `schedule(dynamic)` to handle the load imbalance.

Projected per-chain-move time (with flattened work list on 4 threads):
```
Work units roughly ∝ nChar × kStates²:
  Neo:        160 × 2² = 640 units (one work item — can't split)
  Trans k=3:   44 × 3² = 396 units
  Trans k=4:    9 × 4² = 144 units
  Trans k=5:    4 × 5² = 100 units
  Total: 1280; trans subtotal: 640 — perfectly balanced with neo!
Critical path on 4 threads with dynamic schedule:
  max(neo=640, trans=396+144+100=640) ≈ 640 / 1280 = 50% of serial
Effective speedup: 2× on parallel portion.
Per-chain-move: 26 ms serial + 64 ms parallel × 0.5 = 58 ms
≈ ~1.55× wall improvement (worse than I'd hoped — see note below).
```

**Note:** the prior 2.14× estimate assumed full 4-way parallelism. The
actual neo-vs-trans split is 50/50, so flat-list parallelism caps at
~2× on the parallel portion → ~1.55× wall. To hit ≥ 2× we'd need to
split the neomorphic partition (no kPrime subgroup structure but could
chunk by character count: e.g., 4 chunks of 40 chars each). That's
extra complexity (stride into `part.tipStates`); recommend as a Phase-1
follow-up only if the as-implemented bench shows < 1.6×.

---

## Inventory of R API surface in per-chain hot path

(Verified by Plan agent via Grep on full source tree.)

- **47 RNG sites** (`R::unif_rand`, `R::norm_rand`, `R::rgamma`, etc.)
  across `do_move_impl` + 10 helper `_impl` functions + `proposals.cpp`.
  All mutate `.Random.seed` — UNSAFE in parallel region.
- **74 Rcpp alloc sites** (`NumericVector(n)`, `IntegerMatrix(r,c)`,
  `clone(...)`) in the per-chain move dispatch. UNSAFE in parallel
  region (R memory manager + PROTECT stack not thread-safe).
- **39 alloc sites in `mcmc_ecology.cpp`** total. **Critically**:
  - 26 are Rcpp at the orchestrator level
    (`cpp_log_likelihood_ecology`, `cpp_partition_log_likelihood_ecology`,
    `const_site_prob_*_eco_single`) — UNSAFE.
  - **13 are inside the inner pruners** (`pruning_jc_acrv_flat_ecology`
    at line 287, `pruning_mkn_acrv_flat_ecology` at line 576): all
    are `std::vector<double>` workspaces. **Thread-safe** (glibc malloc
    is, modulo lock contention).
- Pure-math `R::d*`/`R::p*`/`R::q*` calls (`dgamma`, `dlnorm`, `dbeta`):
  pure functions of arguments — assumed safe — **but** any out-of-domain
  argument triggers `Rf_error` longjmp, which is UB from a non-master
  thread. Mitigation: pre-validate inputs; keep `cpp_log_prior` serial.

---

## Phase 1 (this round) — implement

### Goal

Restructure `cpp_log_likelihood_ecology` (the outer partition loop)
and `cpp_partition_log_likelihood_ecology` (the inner subgroup loop)
so the inner pruners can be invoked from a parallel region. Then
parallelise with ONE flat work-item list.

### Architecture: FLATTEN, don't nest

OpenMP nested parallelism is OFF by default (`OMP_NESTED=false`). With
two `#pragma omp parallel for` levels:
- Default: inner loop runs serially → inner pragma is dead code.
- `OMP_NESTED=true`: oversubscription on a 4-core box → contention
  and possible slowdown.

**Correct pattern: build a single flat work-item list serially, then
parallelise that list.**

```cpp
struct EcoWorkItem {
  int partIdx;          // index into data.parts[]
  int subgroupKp;       // -1 for neo/known; else kPrime value
  // Plus any pre-computed metadata for the inner pruner call.
};

// Build the flat list serially.
std::vector<EcoWorkItem> work;
for (int pi = 0; pi < nParts; ++pi) {
  if (data.parts[pi].type == 0 /*neo*/ ||
      data.parts[pi].type == 2 /*known*/) {
    work.push_back({pi, -1});
  } else {
    // Trans: enumerate kPrime subgroups
    std::set<int> kpSet;
    for (int c : data.parts[pi].globalCharIdx) kpSet.insert(kPrime[c]);
    for (int kp : kpSet) work.push_back({pi, kp});
  }
}

// Parallelise the flat list.
std::vector<double> partialLL(work.size(), 0.0);
#pragma omp parallel for schedule(dynamic)
for (int wi = 0; wi < (int)work.size(); ++wi) {
  partialLL[wi] = compute_work_item_loglik(work[wi], ...);
}

// Sum partial results serially → bit-identical to serial code.
double totalLL = 0.0;
for (double x : partialLL) totalLL += x;
```

The inner pruners are already thread-safe (std::vector workspaces only).
The obstacle is that `cpp_partition_log_likelihood_ecology` (:967-1112)
wraps the inner pruner call in Rcpp temporaries (`neoEl`, `rootFreqs`,
`zPart`, `subStates`, `subZ`). Refactor those temporaries to
`std::vector` so the per-work-item body is purely C++.

### Bit-identity strategy

OpenMP `reduction(+:double)` is order-dependent at ULP level: floating-
point addition isn't associative, and thread schedule determines the
order. Bit-identity vs serial cannot use `reduction`.

**Use per-thread partial-sum accumulators**, as in the code sketch
above. Each parallel iteration writes `partialLL[wi]` (no contention,
each thread writes only its own indices). After the parallel region,
sum serially. Order is fixed → bit-identical to serial code.

### Bit-identity scope

Bit-identity holds for **same-build OMP_NUM_THREADS=1 vs =4**, since
the partial-sum order is fixed by the post-parallel-region serial
merge. Bit-identity does NOT necessarily hold for
**compile-with-OMP vs compile-without-OMP**: putting the body inside
a parallel region may change auto-vectorisation, which can change FP
arithmetic at ULP level. This matters for: do not write tests comparing
pre-T-017 vs post-T-017 builds and expect bit-identity.

### Concrete refactors

| File | Change | LoC est. |
|---|---|---|
| `src/Makevars` (new) | `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)`<br>`PKG_LIBS = $(SHLIB_OPENMP_CXXFLAGS)` | +3 |
| `src/Makevars.win` (new) | Same content | +3 |
| `.gitignore` | Remove `src/Makevars.win` exclusion (line 88) | −1 |
| `src/mcmc_ecology.cpp` | Refactor `cpp_partition_log_likelihood_ecology` (:967-1112) Rcpp temporaries → `std::vector`. Same for `const_site_prob_jc_eco_single` (:876) and `const_site_prob_mkn_eco_single` (:907). Restructure `cpp_log_likelihood_ecology` (:1149-1198) to build a flat `EcoWorkItem` list (one item per neo partition / known partition / trans subgroup), then `#pragma omp parallel for schedule(dynamic)` over the flat list with per-thread partial sums in a pre-allocated `std::vector<double>`. Serial merge after parallel region. Do NOT use `reduction(+)`. Do NOT nest pragmas. | +250 / −80 |
| `src/mcmc.cpp` | Emit warning at `run_mcmc_batch_cpp` entry if OpenMP is unavailable. Optional: read `OMP_NUM_THREADS` env. | +30 / −5 |
| `tests/testthat/test-omp-determinism.R` (new) | Same-seed → same-result tests at nChains ∈ {1, 4} with `OMP_NUM_THREADS ∈ {1, 4}`. Verify bit-identity vs serial since Phase 1 doesn't change RNG. | +60 |
| `dev/profiling/findings.md` | T-017 row | +5 |
| `dev/profiling/log.md` | Round 14 entry | +50 |
| **TOTAL** | | **~450 LoC net** |

### Bench target

`dev/profiling/drivers/rodent_aware_timing.R MKP_TIMING_NCHAINS=4`,
500 iter.

- **Baseline (T-015 measured):** 184 s wall, 2.7 iter/s.
- **Phase-1 floor (acceptable):** ≥ 1.5× improvement → ≤ 122 s wall, ≥ 4.1 iter/s.
- **Phase-1 projection (flat work-list parallel on 4 threads):** ~1.55-1.9× depending on actual neo-vs-trans balance and parallel overhead → ~97-118 s wall, ~4.2-5.2 iter/s.
- **Brief target:** ≥ 2× → likely misses by ~0.1-0.4×. If as-implemented bench < 1.6×, consider splitting the neomorphic partition by character chunks (Phase-1b, ~80 LoC) to recover the gap. Otherwise document and move to Phase 2.

5 timed runs per build. Welch's t-test on median wall.
`MKPRIME_T015_DIAG=1` per-chain chrono accumulators verify the win
lands at the right level (per-chain-move time should drop from
~90 ms to ~50 ms; outer-iter cost mostly chain-move-time × nChains).

### Reproducibility verification

Phase 1 keeps R-side RNG (`R::unif_rand`) **untouched**. Parallelism
is bit-identical to serial because we use per-thread partial-sum
accumulators (not `reduction`) and the sum order is fixed by the
post-parallel-region serial merge. **Bit-identity vs serial is
preserved.** Existing tests should pass without expected-value
updates.

Add a new test `tests/testthat/test-omp-determinism.R`:
1. Same seed, `OMP_NUM_THREADS=1` vs `=4` — bit-identical `logLik`
   trajectory required (no tolerance).
2. Same seed, two repeats with `OMP_NUM_THREADS=4` — bit-identical
   trajectory required.
3. Different seed — trajectory differs (sanity check).

**Critical reminder for the impl agent.** Do NOT use
`reduction(+:double)` — floating-point addition is non-associative,
and `reduction` produces a sum whose ULP-level result depends on thread
schedule. The bit-identity test will fail intermittently. Use:

```cpp
std::vector<double> partialLL(nParts, 0.0);
#pragma omp parallel for schedule(dynamic)
for (int pi = 0; pi < nParts; ++pi) {
  partialLL[pi] = cpp_partition_log_likelihood_ecology(...);
}
double totalLL = 0.0;
for (int pi = 0; pi < nParts; ++pi) totalLL += partialLL[pi];
```

### Risk catalog (Phase-1 relevant)

| # | Risk | Mitigation |
|---|---|---|
| 1 | R-API contamination in supposed-safe paths | Phase 1 keeps `cpp_log_prior` and all R density calls outside parallel region. Only inner pruner parallelised. |
| 2 | Rcpp PROTECT stack corruption in parallel region | Strict audit: every line under `#pragma omp parallel for` must be Rcpp-allocator-free. Add assertion. |
| 3 | Determinism break via thread-schedule-dependent reduction | Use `reduction(+:totalLoglik) ordered`, or accumulate per-partition results into `std::vector<double>` and sum serially. |
| 4 | OpenMP unavailable on Hamilton's gcc | Hamilton gcc 9+ has libgomp by default. Verify pre-deploy. |
| 5 | glibc malloc contention on the std::vector workspaces | Empirical: 4 threads × ~120 allocs per call may not contend in practice. Measure via VTune `-collect threading`. If hot, pre-allocate per-thread workspaces. |
| 6 | Phase-1 doesn't meet brief's 2× target | Documented and authorised by the brief's fallback clause. |
| 7 | Heisenbug under threading | Test BOTH `OMP_NUM_THREADS=1` AND `=4`. Run suite 10× to flush flakies. |

### Implementation soundness checklist

When the impl agent picks up this design:

1. Read this design doc + T-007, T-014, T-015 entries in `findings.md`.
2. Confirm Phase 1 scope with the dispatcher before any code change.
3. Write Makevars + Makevars.win FIRST. Install. Confirm build with
   `Rscript -e 'devtools::install(".", quick=TRUE)'`. Verify OpenMP
   symbols present in the .dll/.so.
4. Run baseline 5-run timing of `rodent_aware_timing.R` with the
   OMP build (no parallelism yet) to confirm build is functional and
   wall matches T-015 baseline.
5. Refactor `cpp_partition_log_likelihood_ecology` Rcpp temporaries
   → `std::vector` **serially first**. Run full test suite. Must pass
   228/228 BEFORE adding the pragma.
6. Add `#pragma omp parallel for reduction(+:totalLoglik) num_threads(min(nChains, omp_get_max_threads()))` to the partition loop.
7. Re-run tests + `test-omp-determinism.R`. Verify bit-identity.
8. Bench. Report wall before/after.
9. File T-017 row in findings.md, Round 14 in log.md.
10. Commit on branch `t017-pt-threaded`, push.

---

## Phase 2 (deferred — next round)

Not in scope this round. Would land:

- New `src/chain_rng.h` with `ChainRng` struct (~50 LoC).
- Thread `ChainRng&` through 47 RNG sites across `mcmc.cpp`, `proposals.cpp`, `tree_moves.cpp` (~400 LoC).
- Hoist 74 Rcpp alloc sites in `do_move_impl` to per-chain `McmcState` workspaces (~200 LoC).
- `#pragma omp parallel for` around the per-chain move dispatch loop at `src/mcmc.cpp:6129`.
- Projected ~3.7× wall improvement (4-chain ideal with serial swap + sample overhead).
- 700-1100 LoC net diff.
- **Breaks RNG order** across many tests. Many expected-value updates needed.

Recommend splitting Phase 2 into Phase 2a (RNG threading only, drop-in
compatible) and Phase 2b (alloc hoisting + parallel dispatch).

---

## File-line references

Plan-agent line numbers may be off by ±5 (it read a worktree based on
an earlier commit). Re-grep before editing:

- `src/mcmc_ecology.cpp:967-1112` — `cpp_partition_log_likelihood_ecology` (primary refactor target)
- `src/mcmc_ecology.cpp:1149-1198` — `cpp_log_likelihood_ecology` (partition loop site for #pragma)
- `src/mcmc_ecology.cpp:876, 907` — `const_site_prob_*_eco_single` (secondary refactor)
- `src/mcmc_ecology.cpp:287, 576` — inner pruners (DO NOT MODIFY — already thread-safe)
- `src/mcmc.cpp:6005-6249` — `run_mcmc_batch_cpp` (PT loop — DO NOT MODIFY in Phase 1)
- `src/mcmc_state.h` — McmcState (future Phase-2 home for per-chain scratch)
- `.gitignore:88` — `src/Makevars.win` exclusion (remove for Phase 1)
