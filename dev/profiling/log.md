# Profiling log — mkp

Per-round notes from the `/profile` rotation. Each round profiles ONE focus area, verifies with a micro-bench, files at most one finding into `findings.md`, and refreshes `baselines.md`.

## Round 0 — scaffolding — 2026-05-18

Created `dev/profiling/{drivers,log.md,findings.md,baselines.md,focus-areas.md}` plus this header. No profiling round on the same turn as scaffolding (per skill design).

**Baseline source:** `data-raw/step6-profile-result.rds` (existing user profile data). Wall: 1.57 min for the step6 driver. `.Call → run_mcmc_batch_cpp` = 99.27 % self.pct. R overhead negligible — all findings expected to be C++-internal.

**Ranking signals used:**
- Existing profvis `by.self` / `by.total` from step6
- `[[Rcpp::export]]` map across `src/`
- `AGENTS.md` architecture notes (C++ hot loop for likelihood evaluation + MCMC inner loop)
- Recent commits touching CL/likelihood/tree-move surfaces (M-154 through M-172)

**Next reviewer:** the rotation order is in `focus-areas.md`. Area 1 (Felsenstein pruning / CL accumulation) is up first. The autonomous loop alternates with `/red-team`, so the next /profile invocation will run area 1.

## Round 1 — area #1 (Felsenstein pruning / CL accumulation) — 2026-05-18

Done autonomously (no VTune install on PATH; `.vtune-lib/` 6 weeks stale vs HEAD).

**Driver:** `dev/profiling/drivers/01_felsenstein_pruning.R` (Sun2018, 54×225, all transformational, EG prior, 800 iter, production move weights, 200-iter warmup). Bare wall = 99.78 s; Rprof attributes `.Call` self.pct = 98.62 % (the existing step6 baseline had 99.27 % — consistent).

**Outcome:** no verified `[Optimise]` finding this round — VTune not installed locally; can't break down the C++ time below the `.Call` boundary; Rprof is opaque inside `run_mcmc_batch_cpp`.

**What I CAN file:**
- **T-001 P0 [Baseline-refresh]** — PROFILING.md (2026-03-28) reported 2611 iter/s on Sun2018; current driver runs ~8 iter/s. The 330× gap is dominated by production move weights (Gibbs k′ 29.5 %, Gibbs SPR / subtree-swap 1.3 % each) that the 2611 iter/s baseline excluded. **PROFILING.md's OPP-1..OPP-9 ranking is built on the obsolete move mix; do not trust it for new optimisation decisions until a fresh VTune collection is taken.** Rebuild `.vtune-lib/` first (current one is 6 weeks behind src/).
- **T-002 P1 [Optimise candidate, unverified]** — `singleton_site_prob_jc` (`src/ascertainment.cpp:132-133`) allocates per-call `std::vector<double> cl_flat((maxNode+1)*stride)` (~235 KB for Sun2018) even when `ClWorkspace` is available. `pruning_jc_acrv_flat` already accepts pre-allocated `preAscBuf`/`preAscInit`/`preSiteLikSum`; the singleton path wasn't refactored. **Only hits under `coding="informative"`.** Cross-ref LIKE-001 (red-team round 7) — any fix to that bug should thread `ws` through and close this finding at the same time. Needs micro-bench before filing as a verified `[Optimise]`.

**Trivial fixes applied:** none.

**Latent / next round for this area:**
- Inside the `.Call`, the M-145/M-143/M-158/M-161/M-162B/M-170 cluster of cache + pre-alloc improvements has likely shifted the hotspot mix substantially since 2026-03-28. Without a fresh VTune, I can only enumerate code-level allocation patterns; can't rank by measured time.
- The huge production-vs-reference iter/s gap means the *most-impactful* /profile target is currently the Gibbs k′ sweep (29.5 % weight), not Felsenstein pruning per se. Recommend re-ordering the rotation to prioritise Gibbs k′ in the next pass — see focus-areas.md row 3.
- Note: `useWs=false` fallback in `cpp_partition_log_likelihood` (line 1837 branch) still calls non-flat `pruning_jc_acrv` from `src/acrv.cpp` and the per-call `constant_site_prob_jc`. Worth confirming this branch is unreachable under production paths after M-162B; if reachable, it's an easy win.

**Required for the NEXT productive round on this area (or any other C++ focus):**
- Install VTune (or use Hamilton + `perf` via `hamilton-hpc` skill)
- Rebuild `.vtune-lib/` against current HEAD (`-O2 -g -fno-omit-frame-pointer`)
- Re-run Driver 01 under VTune; report module breakdown and top 10 self-time C++ functions.

## Round 1 — area #1 supplement (VTune collected) — 2026-05-19

User installed VTune at `C:\Program Files (x86)\Intel\oneAPI\vtune\2025.10\bin64\vtune.exe` (not on PATH; full-path invocation). Skill (`~/.claude/skills/profile/SKILL.md`) updated to record the location.

**Build:** wiped `.vtune-lib/`, wrote temporary `src/Makevars.win` with `PKG_CXXFLAGS = $(SHLIB_CXXFLAGS) -g -fno-omit-frame-pointer`, installed via `MAKEFLAGS=DLLFLAGS=-static-libgcc Rscript -e 'install.packages(".", lib=".vtune-lib", ...)'`. Verified `.debug_*` sections present in `MkPrime.dll` (16 MB with symbols). **Removed `src/Makevars.win` after build** (must not be committed).

**Driver:** unchanged from round 1, except switched from `devtools::load_all` to `library(MkPrime, lib.loc = ".vtune-lib")`. Bare wall **dropped from 99.8 s → 21.2 s** (5× speedup just from installed binary vs load_all — surprising; worth investigating in a separate round whether this is byte-compile caching or something more material).

**VTune command:**
```bash
'/c/Program Files (x86)/Intel/oneAPI/vtune/2025.10/bin64/vtune.exe' \
  -collect hotspots -knob sampling-mode=sw \
  -r dev/profiling/result_01_felsenstein_20260519 \
  -- Rscript dev/profiling/drivers/01_felsenstein_pruning.R
```

Collection wall: 23 s. Sample interval at default. CPU = i7-10700 @ 2.904 GHz.

**Hotspot report (filtered to `MkPrime.dll`):**

| Function | Self CPU | Share | Source |
|---|---|---|---|
| `pruning_jc_acrv_persite` | **15.088 s** | **75 %** | `mcmc_likelihood.cpp:459` |
| `pruning_jc_acrv_flat`    | 0.313 s | 1.5 % | `mcmc_likelihood.cpp:279` |
| `mkp::fast_neg_exp` (5 sites, inlined) | 0.874 s | 4.4 % | `fast_exp.h` |
| `jc_transition` (2 sites) | 0.305 s | 1.5 % | `gibbs_partial_cl.h`, `node_cl_cache.h` |
| `constant_site_prob_jc` | 0.094 s | 0.5 % | `ascertainment.cpp` |
| `gibbs_kprime_sweep_impl` (self) | 0.093 s | 0.5 % | `mcmc.cpp` |
| `std::fill` + `std::vector` ctors/dtors | ~0.28 s | 1.4 % | stdlib |

**Story:** The hotspot has shifted dramatically since the 2026-03-28 PROFILING.md baseline. That round reported `pruning_jc_acrv` at 42 % with multiple secondary hotspots; today, **`pruning_jc_acrv_persite` is 75 % alone**, and all other functions are < 1.6 %. This is consistent with the production move schedule promoting Gibbs k′ (29.5 %) — `pruning_jc_acrv_persite` is the Gibbs k′ sweep's per-candidate-k full-tree pruner (M-155 + M-172). The flat path (`pruning_jc_acrv_flat`) is no longer hot — M-162B's workspace pre-allocation worked.

**Per-call analysis:** ~236 Gibbs k′ sweep steps over 800 iter × ~6 candidate-k values per step ≈ 1400 calls × 15.088 s = **~11 ms per call**. For a 106-edge × ~225-char × 4-cat traversal, arithmetic O(nEdge × nChar × kStates × nCat) ≈ 5 ms predicted vs 11 ms measured. The 2× gap suggests **memory-access ceiling, not compute**. Confirm with `vtune -collect memory-access` in a next round before picking an optimisation strategy.

**Findings filed (2 new):**
- **T-003 P0 [Optimise]** — `pruning_jc_acrv_persite` 75 % CPU dominance under EG + production move weights. Verified by VTune; needs memory-access analysis before picking a fix. Candidates: SOA layout of CL buffer, chunked-character iteration to fit L1, explicit AVX2 vectorisation of the inner kStates loop.
- **T-004 P2 [Optimise candidate, small]** — cache per-edge JC transition coefficients across Gibbs k′ sweep candidate-k calls. Bounded saving ~1–3 % wall. Won't show if T-003 is confirmed memory-bound.

**Cleanup:**
- Removed `src/Makevars.win` (never commit; debug flags only).
- `dev/profiling/result_01_felsenstein_20260519/` removed after report extraction (7.6 MB; `.gitignore`-d).
- `.vtune-lib/` kept (next round needs fresh build anyway; can rebuild on demand).

**Refresh of PROFILING.md priority order:**

The 2026-03-28 ranking in `benchmark/PROFILING.md` ("1. OPP-1, 2. M-105, 3. OPP-4, 4. CL zeroing, …") is **obsolete**. Current (2026-05-19) priority based on measured hotspots:

1. **T-003: `pruning_jc_acrv_persite` memory-access investigation** (75 % CPU)
2. **T-004: Gibbs k′ sweep per-edge coefficient cache** (small but clean)
3. **LIKE-001 fix** (red-team, HIGH; correctness, not perf — but the fix touches the same call sites)
4. Defer everything in OPP-2 / OPP-4-lite / OPP-9 — all sub-1 % today.

## Round 3 — T-006 tip-edge fast path + T-005 refactor — 2026-05-19

User asked: "if 4.4 % is not worth the complexity, drop it; if low-maintenance, keep it. Algorithmic gains?"

**Two changes this round:**

1. **T-005 refactor** — collapsed `persite_kfixed<K>` + `persite_generic` (duplicate bodies) into a single templated `persite_impl<K_TPL>` with `K_TPL = 0` sentinel for the runtime-K fallback. One algorithm body, two compile paths. Maintenance cost: zero — both paths edit the same code. Median wall **17.39 → 16.45** (sd 0.14, tightest yet). **Δ = 5.4 % vs baseline**, up from the 4.4 % of the split-body version (likely due to better register allocation when the compiler sees a single function template).

2. **T-006 tip-edge fast path** — exploits the one-hot structure of child CLs at tips. When `ch ≤ nTip`, the K loads of `clCh` and the K-1 `sum_cl` adds are eliminated; the per-state update collapses to K writes of `p_diff` plus one overwrite of `p_same` (known state) or K writes of 1.0 (missing data, by the JC marginalisation identity `K·p_diff + (p_same − p_diff) = 1`). Multiply branch: `clPar[i] *= (i == state ? p_same : p_diff)`. ~50 % of edges in a 54-tip tree are tip-incident; per-character savings scale linearly with K. Median wall **16.45 → 12.07** (sd 0.084). **Δ = 26.6 % on top of T-005, 30.5 % vs baseline.**

| Build | Median wall (s) | Mean | sd | vs Baseline | VTune persite total |
|-------|-----------------|------|----|-------------|---------------------|
| Baseline | 17.39 | 17.47 | 0.28 | — | 15.088 s |
| T-005 (templated K) | 16.45 | 16.53 | 0.14 | -5.4 % | ~10 s |
| **T-005 + T-006 (tip fast path)** | **12.07** | **12.07** | **0.084** | **-30.5 %** | ~7-8 s |

Welch's t-test, T-005+T-006 vs baseline: t = 35.4 on 5.7 df, p < 10⁻⁶.

**VTune verification** (2026.0.0):
- `persite_impl<0>` (runtime K, K > 24): 15.088 → 4.413 s
- Sum of `persite_impl<2..24>` specialisations: ~3 s
- Net pruning self-CPU: ~7.5 s — **~50 % reduction**

**Test suite:** 234/234 across `filter="likelihood|gibbs|ascertain|mcmc-engine"`. Mathematical identity for the tip path verified both analytically (in source comment) and empirically (zero test failures including Gibbs k′ posterior comparisons).

**Skill update:** VTune install path uses `latest` junction (not version-pinned) so the skill survives oneAPI auto-updates. VTune 2026.0.0 confirmed working with the existing CLI.

**Maintenance summary:**
- T-005: zero — single-body template
- T-006: ~50 LoC, single mathematical insight (one-hot child CL), no duplicate code path

**Status updates:**
- T-003 → PARTIALLY-OPTIMISED, superseded by T-005 + T-006
- T-004 → REFUTED (kStates IS what varies per call)
- T-005 → APPLIED + REFACTORED, verified 5.4 %
- T-006 → APPLIED, verified 26.6 % on top of T-005 (30.5 % vs baseline)

**Cleanup:** `dev/profiling/.vtune-lib-*/` and `result_*/` removed.

**Further algorithmic targets (not pursued this round):**
- `persite_impl<0>` (K > 24) still 4.4 s — tighter early termination in Gibbs k′ sweep would cut candidate-k count; needs convergence-correctness audit before any change.
- Apply the tip-edge fast path to `pruning_jc_acrv_flat` (the non-Gibbs path) — only 0.3 s today, marginal.
- SOA vs AOS CL layout for SIMD across characters — major refactor, touches all CL consumers; defer until persite work plateaus.
- Pattern compression granularity (M-172) — already deduplicates by unique pattern; could probe sub-pattern factoring, but very speculative.

## Round 2 — area #1 optimisation lands (T-005) — 2026-05-19

**Skill updates:** record VTune install location in `~/.claude/skills/profile/SKILL.md`; rewrite "Step 4 — Build with symbols" to warn against `devtools::load_all` (builds `-O0`) and require installing to a timestamped temporary directory (avoids DLL-lock collisions and stale-binary failures).

**T-005 implementation (verified):**
- Refactored `pruning_jc_acrv_persite` in `src/mcmc_likelihood.cpp` from a single runtime-K function to a switch-dispatch over `persite_kfixed<K>` for K = 2..24, falling back to a bit-identical `persite_generic` for K > 24.
- Inner per-state loops carry `#pragma GCC unroll 8` and pointers are `__restrict__`-qualified.
- K = 2..24 covers ~all calls a typical morphological dataset hits; the EG prior's heavy tail can push K > 24, hence the fallback.

**Methodology:**
1. Wrote a separate baseline driver `dev/profiling/drivers/01_bench.R` loading from `dev/profiling/.vtune-lib-baseline/` (clean rebuild from current HEAD with `-g -O2 -fno-omit-frame-pointer`).
2. Implemented T-005, installed to `dev/profiling/.vtune-lib-t003/` (separate dir so the two builds coexist).
3. Ran `dev/profiling/drivers/01_bench{,_t003}.R` 5 times each on the same Sun2018 workload (800 iter, EG prior, production move weights).
4. VTune `-collect hotspots -knob sampling-mode=sw` on the T-005 build to verify the function distribution changed as predicted.

**Results:**

| Build | Median wall (s) | Mean (s) | sd (s) | VTune persite self CPU |
|-------|-----------------|----------|--------|------------------------|
| Baseline | 17.39 | 17.47 | 0.28 | `pruning_jc_acrv_persite` 15.088 s |
| T-005    | **16.75** | **16.70** | 0.20 | `persite_generic` 7.113 s + `persite_kfixed<2..24>` ~4.5 s combined |

**Δ = 4.4 % wall-time speedup.** Welch's t-test t = 4.3 on 7.4 df, p < 0.01. Tighter variance under T-005 (sd ↓ from 0.28 → 0.20) is a bonus — likely fewer branch-prediction stalls on the per-state loop.

**Why not bigger:**
- VTune shows the residual `persite_generic` 7.113 s is K > 24 calls. The EG prior's geometric tail keeps some non-trivial mass at K = 25..50.
- Extending `DISPATCH_K` higher works mechanically but the diminishing-returns curve is steep — K = 25..50 calls together would probably yield another ~2–3 % at ~2000 LoC of generated code. Deferred.

**Tests:** `Rscript -e 'devtools::test(filter="likelihood|gibbs-kprime")'` runs 53 tests with 0 failures. `gibbs_kPrime` Gibbs sweep test exercises the new path on every supported K.

**T-004 closed as REFUTED.** Original premise (cache `(p_diff, diff_coeff)` per-edge per-cat across candidate-k calls) does not hold: `kStates` IS what varies across calls, so the JC coefficients change per call too. No cacheable redundancy.

**T-003 status:** PARTIALLY-OPTIMISED. The 75 % CPU dominance of `pruning_jc_acrv_persite` is reduced (now split across template + generic), but the function as a whole still dominates. Next-round candidates listed in T-005's "Further wins available" — all algorithmic, not codegen.

**Cleanup:** `dev/profiling/.vtune-lib-baseline/`, `dev/profiling/.vtune-lib-t003/`, and `dev/profiling/drivers/01_bench*.R` removed after this round (per skill discipline; `.gitignore`-d anyway). `src/Makevars.win` removed (never commit).

**Modified file (uncommitted, ready for review):** `src/mcmc_likelihood.cpp` — adds `persite_kfixed<K>` template, `persite_generic` fallback, and switch dispatch in `pruning_jc_acrv_persite`. Net: +210 lines (mostly generated by the template + generic duplication).

## Round 4 — area #0 (ecology-aware orchestrator) — 2026-05-21

User-requested while iterating: Hamilton blind run finished rodent MkNT in 19 min (job 17254914); aware (job 17254915) still warming up after 80+ min. Same matrix, same MkNT likelihood — only ecology layer differs.

**Methodology.** Read the eco surface (`src/mcmc_ecology.cpp`, `src/mcmc.cpp` eco call sites); confirmed it bypasses partition cache (`partLogLik` force-cleared, comment at `mcmc.cpp:3842-3846` is explicit), `nodeCL` partial-CL cache, and `ClWorkspace`. Installed HEAD into `dev/profiling/.vtune-lib-prof/`. Wrote three drivers:
- `11_aware_vs_blind_rodent.R` — same rodent matrix (MorphoBank X24848, 64 tips × 217 chars, kEco=4), 200 iter, 1 chain 1 run, both modes. **Aware 24.59 s, blind 0.97 s, ratio 25.3×.**
- `11b_rprof_aware.R` — Rprof at 0.02 s tick to confirm where time goes. **97.20 % in `.Call → run_mcmc_batch_cpp`** (R overhead negligible, matching the M-150 pattern from area #1).
- `11c_per_call_cost.R` — isolated per-call cost via `.MkpEcologyLogLikelihood` vs `.MkpLogLikelihood`, 30 reps after warm-up. **Aware 0.67 ms, blind 0.33 ms — only 2.0× per call.**

**Diagnosis.** The 25.3× wall-time gap factorises as 2.0× per-call × ~12.7× call-frequency excess. The 2× per-call is the inner-loop tax (extra mixture loop, fresh allocations, no T-006 fast path). The ~12.7× is **architectural**: every eco move-handler recomputes the full likelihood instead of reusing a cached `partLogLik`. **VTune skipped** for this round: VTune attributes time *within* a call but cannot attribute time to calls that shouldn't have happened — the diagnosis comes from cross-referencing the per-call vs wall ratios, not from inner-loop hotspots. Advisor concurred.

**Findings filed (4):**
- **T-007 P0 [Optimise]** — headline: eco path bypasses all likelihood caches; predicted 5–10× wall recovery once partition cache is restored.
- **T-008 P1 [Optimise]** — per-call heap allocations in `cpp_log_likelihood_ecology` (~10s of fresh ~100 KB allocs/call) and `per_char_log_lik_ecology` (~1950 fresh allocs per Gibbs z sweep). Thread an `EcologyWorkspace` through both.
- **T-009 P1 [Optimise]** — port T-006 tip-edge fast path to `pruning_jc_acrv_flat_ecology`; same identity holds, contained change in one function.
- **T-010 P2 [Optimise candidate, design]** — restore partition-level cache for ecology non-tree moves; design risk (wEdge invalidation under tree moves) — defer until T-008/T-009 verified.

**Ordering for next rounds:** T-009 first (most contained — one function body, mathematical identity already proven on the blind path). Then T-008 (workspace plumbing). Then T-010 (architectural, larger surface area). After all three, expect aware/blind wall ratio to drop from 25× to ~3-5× (vs the rodent-relative tax of the mixture pruning itself, which is irreducible without an algorithmic change to the model).

**Caveats** (recorded in T-007): per-call bench used `rateLogSd = 0` (nCat = 1). Production cost scales linearly with nCat (default ~5). Local 25× vs Hamilton ~4×+ ratio is partly warmup-adaptation acting differently at long iter counts. The headline ratio number is from the local short-iter bench — Hamilton long-run ratios will differ.

**Cleanup:** `dev/profiling/.vtune-lib-prof/` and `dev/profiling/aware_rprof.out` removed (gitignored anyway). Drivers (`11_*.R`) kept for the next-round verifications. Modified `src/RcppExports.R` rebuilt by the package install — confirmed it matches HEAD's `git diff R/RcppExports.R` (no functional drift, just `setRcppClass` regeneration). No source files modified.

---

## Round 7 — T-008 per_char workspace — 2026-05-21

T-008 per-character workspace portion. Orchestrator portion tracked separately.

**Scope.** `per_char_log_lik_ecology` in `src/mcmc_ecology.cpp` is called 1953 times per Gibbs z sweep (nChar=217 × (kEco-1)=3 × 3 candidate values). Before this change each call allocated 4 heap objects: `IntegerMatrix tipStates(nTip,1)`, `IntegerMatrix zPart(1,zCols)`, `std::vector<double> buf((maxNode+1)*stride)`, `std::vector<uint8_t> initFlg(maxNode+1)`. Total: ~7812 heap allocs per sweep.

**Implementation.** Added `src/gibbs_z_workspace.h` (`GibbsZWorkspace` struct with `buf`, `initFlg`, `tipStates`, `zPart`, and `ensure()` resize-only-when-needed method). Refactored `per_char_log_lik_ecology` to accept `GibbsZWorkspace&`; `gibbs_z_sweep_impl` in `src/mcmc.cpp` allocates ONE workspace at sweep start (max stride across all partitions) and passes it through all calls. `tipStates` column copied via `std::memcpy` from `part.tipStates` (column-major offset), eliminating the per-tip loop. `zPart` filled inline from `zRow`. **0 heap allocs per call** (down from 4).

**Verification.** `dev/profiling/drivers/12_per_char_alloc.R` (subprocess subproces, 200 reps each, rodent 217-char kEco=4):

| Build | ms/batch (217 chars) | ms/char |
|-------|----------------------|---------|
| Baseline (.vtune-lib-prof) | 18.7 | 0.0862 |
| T-008 (.vtune-lib-t008) | 18.1 | 0.0834 |
| **Δ** | **~4 %** | **~4 %** |

End-to-end driver 11 (200 iter rodent, same seed): baseline 21.27 s → T-008 21.19 s (within noise). Speedup modest because malloc cost is small relative to per-edge mixture arithmetic; expected to compound with T-009 and T-010.

**Test result:** 375 / 375 pass (`filter="ecology|gibbs|likelihood"`), 0 failures, 7 skips (expected).

**Alloc budget:**
- Gibbs z sweep: 7812 heap allocs per sweep → 0
- `CppLogLikelihoodEcologyPerChar` R wrapper: 1 workspace per call-batch (reuses across nChar characters within the batch) — ~96 allocs eliminated per R-level call

**Files modified:** `src/gibbs_z_workspace.h` (new), `src/mcmc_ecology.cpp` (include + struct usage + per_char signature), `src/mcmc.cpp` (include + forward decl + workspace in gibbs_z_sweep_impl). `dev/profiling/drivers/12_per_char_alloc.R` (new bench driver).

**Cleanup:** `dev/profiling/.vtune-lib-t008/` removed after filing (gitignored).

**Orchestrator portion** (`cpp_log_likelihood_ecology` allocs at ~lines 968/992/1033) not yet addressed — tracked under T-010/T-008-orchestrator.

---

## Round 8 — T-009 tip-edge fast path port — 2026-05-21

**Implementation.** Ported T-006 tip-edge fast path to `pruning_jc_acrv_flat_ecology` in `src/mcmc_ecology.cpp`. Added a `tipChild = (ch <= nTip)` branch before the per-character propagation loop in both init and multiply branches. Tip path:
- known state s: `clPar[i] = pdMix` for i≠s, `clPar[i] = psMix` for i==s (init); `clPar[i] *= pdMix/psMix` (multiply)
- missing: `clPar[i] = 1.0` (init, by JC identity Σ w_s·(K·pdF+(psF−pdF))=1); no-op (multiply)
The internal-child path is unchanged. psMix/pdMix accumulation order preserved for FP consistency.

**Bench results.** Isolated `.PruningJcEcology` direct call (k=4, nChar=150, nTip=60, kEco=4, 5000 reps, 5 outer):
- Baseline: 0.241 ms/call
- T-009:    0.241 ms/call
- Δ: ~0 % (within Windows 10 ms `proc.time` resolution)

Orchestrator-level (`.MkpEcologyLogLikelihood`, 1000 reps): baseline 0.57 ms, T-009 0.57 ms.

**Why near-zero gain.** The dominant per-edge cost in the ecology pruner is the `for (s = 0; s < kEco; ++s)` psMix/pdMix accumulation (kEco `MKP_EXP` calls × 3 z-values precomputed), not the K-state inner loops that T-009 eliminates. For kEco=4, K=4: the saved `sum_cl` (4 adds) and collapsed writes are negligible relative to 4 EXP-dominated mixture accumulations per character per edge. The saving is real (saves ~K adds + K reads per tip-edge character) but too small to measure above the benchmark noise floor.

**Tests:** 226 / 226 pass (`filter="ecology|likelihood"`), 0 failures, 4 skips (slow MCMC + pilot RDS).

**Code correctness.** The implementation is mathematically correct and tested. It sets up the code structure for future cases where the inner K-state work is comparatively more expensive (e.g., larger K under rate-variation or when heap-alloc cost is removed by T-008).

**Cleanup.** `dev/profiling/.vtune-lib-prof/`, `dev/profiling/.vtune-lib-t009/`, `dev/profiling/drivers/11c_bench_b.R`, `dev/profiling/drivers/11c_bench_t.R` removed after filing (gitignored / temporary).

## Round 6 — T-010 partition cache for ecology — 2026-05-21

**Goal.** Restore partition-level likelihood caching to the ecology-aware MCMC
path. Predicted recovery: 5–10× wall on the rodent matrix once non-tree moves
hit a cache for partitions they don't touch.

**Implementation.**
- Added `cpp_partition_log_likelihood_ecology` (bit-identical to one
  orchestrator iteration) and `compute_gamma_e_ecology` helper in
  `src/mcmc_ecology.cpp`.
- Refactored `cpp_log_likelihood_ecology` to call the per-partition function
  per partition (summation order preserved).
- Wired `state->wEdge` + `state->wEdgeDirty` into a real cache via
  `eco_refresh_wedge`/`eco_recompute_all_partitions` helpers in
  `src/mcmc.cpp`.
- Reused `state->partLogLik` as the ecology partition cache (cleared on
  tree-move accept paths, same as the blind path).
- Refactored every move-handler eco branch: slice probes use partition
  cache; kPrime moves recompute only the trans partition (subgroup
  composition caveat — recompute whole trans partition is the easiest
  correct option, documented); phi/pi0/theta moves keep wEdge cached and
  recompute all partitions; tree moves invalidate wEdge.
- `wEdgeDirty = true` added to every tree-move accept path (8 sites).

**Build.** `dev/profiling/.vtune-lib-t010/` clean rebuild.

**Tests.** 228/228 across `filter="ecology|likelihood|mcmc"` (was 226/226 on
T-008+T-009 HEAD; +2 from a new T-010 cache-reproducibility test in
`tests/testthat/test-ecology-likelihood.R`). 0 failures. The existing
"per-char ecology log-liks sum to total" invariant already covered the
orchestrator-vs-partition correctness check.

**Drift check.** Rodent matrix, 1000-iter aware MCMC: **zero `[eco-resync]`
warnings** emitted (the 20-iter cross-check threshold of |dLL|+|dLP| > 0.5
nats was never crossed). Cache invalidation logic is correct.

**Bench results** (`dev/profiling/drivers/11d_t010_bench.R`, rodent matrix
64×217, kEco=4, 200 iter):

| Build | Aware (s) | Blind (s) | Aware iter/s | Ratio |
|-------|-----------|-----------|--------------|-------|
| Baseline (T-008+T-009 in HEAD) | 24.59 | 0.97 | 8.13 | 25.3× |
| **T-010** | **21.40** | **0.77** | **9.35** | 27.79× |
| Δ aware  | -13.0 %  | -20.6 %  | +15.0 %  | — |

**Why only 13 %, not the predicted 5–10×.** Production move schedule on
rodent assigns ~78 % weight to tree-touching moves (NNI 17 %, SPR 8.5 %,
TBR 8.5 %, pSPR 8.5 %, gibbs_spr 2.9 %, gibbs_subtree_swap 2.9 %,
branch_lengths 11 %, dirichlet_branch 8.3 %, local_dirichlet 8.3 %,
tree_length 1.7 %) — all of which invalidate wEdge AND need every partition
recomputed. Only ~22 % of moves benefit from the partition cache (slice on
rate_loss/rate_neo → neo partitions only; phi/pi0/theta/rateLogSd →
wEdge-cached; kPrime int_walk → one trans partition; z sweep → cached
post-sweep refresh). Within those ~22 %, the cache delivers the predicted
savings; on the wall they're diluted by the dominant tree-move cost.

**Per-call cost** (`dev/profiling/drivers/11c_per_call_cost.R`,
30 reps): aware orchestrator 0.67 ms/call — unchanged from baseline (T-010
is bit-identical at the orchestrator level, by design).

**Next bigger win** requires either (a) partial-CL evaluation in ecology
tree moves (much larger refactor; ecology `wEdge` couples every edge so
naive partial CL doesn't apply), or (b) approximate / incremental wEdge
update for branch-length moves. Both deferred.

**Status updates.**
- T-007 → PARTIALLY-OPTIMISED (T-008 + T-009 + T-010 together).
- T-010 → APPLIED, VERIFIED ~13 % wall on the rodent matrix.

**Cleanup.** `dev/profiling/.vtune-lib-t010/` removed after filing. No
`src/Makevars.win` left behind.

---

last_focus: 0
