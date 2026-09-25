# Profiling log — mkp

Per-round notes from the `/profile` rotation. Each round profiles ONE focus area,
verifies with a micro-bench, files at most one finding as a **GitHub issue**
(label `profiling`; the issue number is the id), and refreshes `baselines.md`.

**Findings are never written to a file.** `findings-archive.md` is frozen
anti-duplication memory, not a work list. Round records, `[AT-LIMIT]` verdicts and
same-round fixes stay here — they are records, not work.

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

## Round 9 — T-011 partial-CL for ecology tree moves — 2026-05-21

**Goal.** Predicted ~5× wall recovery on rodent via Lever 1 (node-CL cache
for the ecology pruner) per the T-011 spec.

**Outcome.** **Foundation landed; 0 % measured wall delta**; no regressions.
Reframed as "T-012-ready scaffolding" — see T-011 in findings.md.

**Methodology.**
1. Pre-implementation advisor consultation. Prediction: spec's "simplest"
   Lever 1 (invalidate-everything-on-wEdgeDirty) gives 0 % wall because
   T-010's partition cache already captures every cache-eligible non-tree
   move. Meaningful wall recovery requires either Lever 2 lite (per-edge
   wEdge dirty detection) or an M-158-style save/restore around tree-move
   evals — both are ~500 LoC of code I didn't have appetite to ship and
   stress-test under the round budget.
2. Implemented the cache structure + invalidation skeleton anyway, as the
   spec's required scope. Verified per-call cost unchanged, drift clean,
   tests green.
3. Post-implementation advisor consultation. Confirmed reading: don't ship
   dead code, document the scaffolding for what it is, commit honestly.
4. Removed dead `pruning_*_ecology_cached` helpers from `src/mcmc_ecology.cpp`
   per advisor; kept header + invalidation wiring + design notes.

**Verification.**
- `dev/profiling/drivers/11e_t011_bench.R` (rodent 64×217, 200 iter):
  baseline (T-010 HEAD) aware **20.38 s**, blind **0.76 s**, ratio **26.82×**
  → T-011 aware **20.47 s**, blind **0.79 s**, ratio **25.91×**. Within
  bench noise either direction.
- `dev/profiling/drivers/11g_t011_drift.R` (5000-iter aware MCMC):
  wall **473.98 s**, **0 `[eco-resync]` warnings**. Cache invalidation
  is correct at every wired site.
- `Rscript -e 'devtools::test(filter="ecology|likelihood|mcmc")'`:
  **228 / 228 pass**, 0 failures, 6 skips (matches T-010 baseline).
- `dev/profiling/drivers/11c_per_call_cost.R`: per-call aware
  orchestrator cost unchanged (~0.33–0.67 ms, within noise of T-010).

**Code landed.**
- `src/ecology_cl_cache.h` (~360 LoC): EcoCacheUnit, EcoCLCache, build /
  detect / snapshot / equality helpers. All `[[maybe_unused]]` with a
  STATUS comment block explaining the deferral. Future Lever 2 round
  consumes this directly.
- `src/mcmc.cpp`: include header, add `EcoCLCache ecoCL;` to McmcState,
  17 `state->ecoCL.invalidate_*` calls mirroring every existing
  `state->nodeCL.invalidate_*` call site.
- `src/mcmc_ecology.cpp`: include header. Dead cached pruner helpers
  removed per advisor before commit.

**Side-finding spawned.** `gibbs_kprime_sweep_impl` uses the blind pruner
under ecology mode (no wEdge / phi / z consultation). Possible
correctness bug — spawned as a separate task to audit.

**T-007 status:** **stays PARTIALLY-OPTIMISED.** The 25× wall gap is
essentially unchanged.

**Next round recommendation.** T-012: NNI + branch-simplex partial-CL
with save/restore via `EcoCLCache`. Expected gain ~8–12 % wall on
NNI alone (17 % move weight × ~50 % dirty-set reduction). Add SPR
later once the NNI path is proven; SPR/TBR/pSPR have wider dirty
sets and lower marginal return.

**Cleanup.** `dev/profiling/.vtune-lib-t011/` removed after filing.
No `src/Makevars.win` left behind.

---

## Round 10 — T-012 partial-salvage of dead background agent — 2026-05-21

**Goal.** Implement Lever 1 + Lever 2 lite for the ecology pruner (per the
T-011 finding's "next round recommendation": NNI partial-CL with
save/restore via `EcoCLCache`, expected ~8–12 % wall).

**Method.** Dispatched an Opus background subagent (`agentId
a80bae8b738eab40c`) in an isolated worktree with the full T-012 brief.
The agent worked for ~14 minutes (11:22 → 11:36) modifying
`src/ecology_cl_cache.h` (+513 LoC), `src/mcmc_ecology.cpp` (+285 LoC),
and creating `dev/profiling/drivers/12a_bitident.R`; then went silent.
By 11:55 there were no R/build/g++ processes in the tasklist, the
agent's JSONL transcript was 0 bytes (never flushed), and
`TaskOutput a80bae8b738eab40c` reported "no task found" — the harness
had reaped the task entry. Worktree contents survived.

**Salvage.** Inspected the diff from `worktree-ecology-aware` parent.
What landed:
- `EcoCacheUnit` + `EcoCLCache` allocation, per-(cat, node) CL storage
  with tip init, and a full postorder downpass.
- `populate_eco_cache_full` orchestrator (per-partition, per-subgroup
  units; bit-identical at the root to the legacy pruner).
- `eco_cache_total_loglik` (sums per-unit root logliks with relabel /
  ascertainment corrections, mirroring the legacy partition wrapper).
- `eco_cache_partial_eval_nni` + `find_dirty_nni_eco` + dirty-edge
  walker over `detect_dirty_child_nodes` (the T-011 helper, now
  consumed) + `save_dirty_eco_cls` / `restore_dirty_eco_cls`.
- R-callable `.CppLogLikelihoodEcologyCached` exporting the full
  `populate_eco_cache_full` + `eco_cache_total_loglik` path for
  bit-identity testing.

What did NOT land:
- Wiring into `do_move_impl` for NNI. `src/mcmc.cpp` is untouched.
- Bit-identity testing of the **partial-eval** path
  (`eco_cache_partial_eval_nni`) — only the full-eval path can be
  driven from R as of this commit.
- The 5000-iter drift stress and the production-mode wall bench: with
  no consultation wiring, both reduce to the T-011 baseline (which
  was verified clean).
- The agent's own `12a_bitident.R` was broken (referenced a
  non-existent `mkd$dataPtr`); rewritten by the salvage step using the
  `prepare_mcmc_data` + test-helper pattern from
  `tests/testthat/test-ecology-plumbing.R`.

**Verified.**
- `dev/profiling/drivers/12a_bitident.R`: legacy
  `.CppLogLikelihoodEcology` vs cached `.CppLogLikelihoodEcologyCached`
  on a 10-tip / 6-character / kEco=3 random fixture: **|diff| = 0
  exactly**. The full-eval cache mechanics are correct.
- `Rscript -e 'devtools::test(filter="ecology|likelihood|mcmc")'`:
  **228 / 228 pass**, 0 failures, 6 expected skips.
- `dev/profiling/drivers/11e_t011_bench.R` (200-iter rodent):
  aware **20.58 s**, blind **0.78 s**, ratio **26.38×** — within bench
  noise of the T-011 baseline (20.47 / 0.79 / 25.91×).

**Decision.** Commit the salvaged work as T-012 PARTIAL — bit-identity
gate is the right correctness anchor for the full-eval path; partial-eval
correctness and production wiring belong to a follow-up (T-013).
The dead-code anti-pattern advisor flagged on T-011 doesn't apply here
because the new functions ARE called via `CppLogLikelihoodEcologyCached`
and the bit-identity driver exercises them on every run.

**T-007 status:** stays **PARTIALLY-OPTIMISED**. The 25× aware-vs-blind
gap is essentially unchanged. The infrastructure to close it is now
substantially in place; T-013 is the remaining work.

**Cleanup.** `dev/profiling/.vtune-lib-t012/` removed after filing. No
`src/Makevars.win` left behind.

---

## Round 11 — T-013 NNI partial-CL wiring (in conversation) — 2026-05-21

**Goal.** Verify the T-012 partial-eval mechanics by exposing them to R,
then wire NNI partial-CL into production via `state->ecoCL`. Spec target:
aware wall ≤ 18 s (≥ 12 % improvement vs T-012 baseline 20.58 s); ≥ 5 %
floor for reporting.

**Method.** The first T-013 background agent died after ~15 min having
added the R-callable `.CppPartialEvalEcologyNNI` wrapper +
`12b_partial_bitident.R` driver + `12c_rodent_dirty_fraction.R`
analysis, but never wired production. Salvaged the agent's diff
(`src/mcmc_ecology.cpp` modifications + the two new drivers) and
completed the work in the foreground.

**Bit-identity gate (P0a).** 20 random in-place NNI swaps on a
10-tip / 6-char / kEco=3 fixture — every sample passed with exact
`|partial − fresh_new| = 0` and `|restored − fresh_old| = 0`. The
mechanics (populate, dirty-set walk, partial recompute, save/restore)
are mathematically correct.

**Empirical dirty fraction (P0b).** `12c_rodent_dirty_fraction.R` on
the rodent topology (60 tips, 59 internal nodes, 118 edges, kEco=4):
30 random NNI swaps. Dirty count mean 45.7, median 54, max 57
(fraction of nInternal: mean 0.77, max 0.97). wEdge-dirty edges
mean 89.1, median 109 (~92 % of edges). **Naive partial-eval saving
1 − dirty/nInternal = 7.6 % median.** This is the dealbreaker: even
with perfect partial-eval, the rodent move geometry leaves so much
of the tree dirty that the savings are tiny.

**Production wiring landed.** Lazy `populate_eco_cache_full` before
NNI proposals + partial-eval branch in the main MH eco accept block
(mirrors the M-158 NNI pattern with `update_topo_nni` on the cache's
TreeNav, partial-eval call, save scratch on the cache for rollback,
revert TreeNav on fallback). Reject path: restore CLs + reverse swap
symmetrically. Wiring is gated behind `MKPRIME_ECO_PARTIAL_CL=1`
env var (default OFF) — see Wall below for why.

**Wall (P0c).**
- T-012 baseline: aware 20.58 s, blind 0.78 s, ratio 26.4×.
- T-013 gate OFF (default): aware 20.76 s, blind 0.77 s, ratio 27.0×.
- T-013 gate ON: aware 20.69 s, blind 0.78 s, ratio 26.5×.

All within bench noise. With gate ON the partial-eval path engages
on every NNI; on rodent the dirty fraction is so high that the
dirty-walk + save/restore overhead approximately cancels the per-call
savings.

**5000-iter drift gate.** With wiring engaged (no env gate),
wall 506.22 s, **0 `[eco-resync]` warnings**. Correctness confirmed.

**Decision: gate behind env var.** Default OFF avoids the small
overhead on the rodent workload (where partial-eval is net ~0 %).
Setting `MKPRIME_ECO_PARTIAL_CL=1` engages the partial-eval path —
useful for benchmarking on larger trees where dirty fraction should
shrink (path-to-root grows ~log(nTip) while nInternal grows ~nTip).
This was the cleanest way to land the infrastructure without forcing
a possible regression on production runs.

**Tests.** 228 / 228 pass (`filter="ecology|likelihood|mcmc"`), 0
failures, 6 expected skips. Bit-identity driver `12b_partial_bitident.R`
PASSES (|diff| = 0 across all 20 samples + rollback).

**T-007 status.** Stays **PARTIALLY-OPTIMISED**. The 25× wall ratio
is essentially unchanged. The architectural foundation for partial-CL
is now complete; further wall recovery on rodent-sized trees would
require either (a) algorithmic redesign of the wEdge dependency
(approximate marginals that don't propagate globally) or (b)
per-character CL caching with a restructured pruner — both are major
multi-week refactors out of scope for the current focus rotation.

**Cleanup.** `dev/profiling/.vtune-lib-t013/` removed after filing.
No `src/Makevars.win` left behind.

---

last_focus: 14

## Round 11 (T-014) — 2026-05-21

**Target:** "aware PT amplification" — aware costs 4.2× more per outer iter
under `nChains=4` PT than `nChains=1` (10.6 → 2.5 iter/s); blind only 2.2×
(228 → 101.6 iter/s). User-requested while triaging stalled Hamilton job
17258771.

**Method.** No VTune this round — code inspection of the PT outer loop
(`src/mcmc.cpp:6116-6219`) plus targeted modulus toggle on the eco
drift-resync gate at `:6188`. Driver: `dev/profiling/drivers/rodent_aware_timing.R`
with `MKP_TIMING_NCHAINS={1,4}` (added a one-line env-var hook for chain
count).

**Found.** The drift-resync gate fires for every chain when
`iter % 20 == 0`, calling `cpp_log_likelihood_ecology` + `cpp_log_prior`
from scratch (uncached). Bench: switching `% 20` → `% 200` gives
**Δ +12 % wall** at nChains=4 (2.5 → 2.8 iter/s) and ~0 % at nChains=1.

This accounts for ~25 % of the excess PT amplification; the remaining
~75 % is per-chain ecology state (zMatrix / wEdge / gammaE) defeating
inter-chain cache locality — a much larger refactor (filed as T-015
candidate in the finding, not pursued this round).

**Filed:** T-014 in findings.md, status OPEN, kind [Optimise], priority P0.

**Hamilton implication.** Even with the 12 % fix, aware nChains=4 PT
projects 12 h on rodent at 30 k iter warmup, hitting the v2 walltime
limit. User now has a clean choice: nChains=1 + nRuns=4 (no PT, viable
today) or wait for shared-arena cache (T-015).

**Cleanup.** src/mcmc.cpp bench patch reverted to `% 20`. No
`.vtune-lib-*/` created (not needed for this round). No
src/Makevars.win touched.

---

## Round 12 (T-015) — 2026-05-21

**Target.** Investigate the "remaining ~75 % of aware-PT amplification"
that T-014 left on the table, with four candidate hypotheses (cache
eviction between chains, hidden per-chain rebuild, PT swap cost, false
sharing) and a brief that explicitly said *don't presume a fix exists*.

**Method.** Advisor pre-consult flagged a math-framing issue in the
T-014 write-up: the "1.89× excess amplification" was computed as
(aware 4.24× slowdown) / (blind 2.24× slowdown), treating blind's
super-linear scaling as the no-amplification baseline. Linear scaling
under PT is 4× per outer iter (you do 4× the chain-moves), not the 2.24×
that blind achieves. Blind's 2.24× **already contains a cache-sharing
bonus**; aware can't be expected to match it unless aware is also
memory-bound (which it isn't — aware does ~4× more EXP per edge per cat
than blind, by T-007's per-call cost analysis). Advisor recommended:
chrono-instrument the PT loop before VTune, since hotspots can't
discriminate between "compute-bound, scales linearly" and "memory-bound,
scales sub-linearly".

**Instrumentation.** Added env-gated (`MKPRIME_T015_DIAG=1`) per-chain
chrono accumulators inside the PT loop at `src/mcmc.cpp:6121-6249`:
one bracket per move dispatch (`do_move_impl` / `slice_scalar_impl` /
`slice_kprime_hyper_impl`), one per drift-resync gate firing, and a
single accumulator around the `std::swap(*states[iPair], *states[jPair])`
at `:6243`. Emitted to stderr at end of each batch via REprintf.
Inert when env unset; instrumentation cost adds <1 % under the gate.

**Bench.** `dev/profiling/drivers/rodent_aware_timing.R` (rodent
60×217, kEco=4, 500 iter) at `MKP_TIMING_NCHAINS ∈ {1, 2, 4}` with
diag enabled. Numbers are per-chain-move mean μs averaged across
chains and both batches (200 + 300 iter):

| nChains | Blind iter/s | Aware iter/s | Blind per-move μs | Aware per-move μs |
|---|---|---|---|---|
| 1 | 222 | 11.5 | 1475 | ~82700 |
| 2 | 160 | 5.3  | ~1563 | ~91300 |
| 4 | 102 | 2.7  | ~1659 | ~90600 |

**Per-chain-move growth from 1ch → 4ch:** blind +12.5 %, aware +9.5 %.
The small uniform ~10 % growth affects both modes symmetrically — **no
aware-specific locality cost**. **Caveat:** per-chain weight adaptation
produced different move mixes across the 1ch/2ch/4ch runs (e.g. blind
1ch nni:11.7 % vs 4ch nni:4.7 %; aware 1ch nni:15.9 % vs 4ch nni:1.1 %),
so per-chain-move μs averages over different work distributions. The
12.5 % vs 9.5 % gap is within this confound; what matters for T-015 is
the pattern — comparable growth between modes — not the exact %s.

The "super-linear" pattern blind shows in outer-iter throughput
(222 → 102 = 2.18× slowdown vs 4× linear) is explained entirely by
**fixed per-outer-iter R overhead** (sample save, streamed log write,
R callback at batch boundary). On 4 chains:
- Blind 4ch outer iter = 9.8 ms; 4 × per-move (1.66 ms) = 6.6 ms; **3.2 ms is fixed overhead (33 %)**.
- Aware 4ch outer iter = 370 ms; 4 × per-move (90.6 ms) = 362 ms; **only 8 ms is fixed overhead (2 %)**.

Multiplying chains scales the move work but the fixed overhead is
already paid. For blind that fixed overhead is comparable to its tiny
per-move cost, so adding chains looks sub-linear; for aware it's
negligible, so adding chains looks linear.

**Swap cost (hypothesis #3).** `swap_total_ns` was 0.000 s under both
modes across all chain counts — `std::swap(McmcState)` is move-based
(Rcpp SEXP refcount-swap + std::vector pointer-swap) and registered
below the chrono resolution. <0.001 % of wall.

**Drift-resync (hypothesis #2 residual).** At default `% 200` cadence
the gate fires 2-3 times per 500-iter batch. Aware 4ch: 232 ms total
out of 184 s wall = 0.13 %. T-014's fix already collapsed this; no
remaining headroom.

**Lazy populate (hypothesis #2 primary).** `populate_eco_cache_full`
at `src/mcmc.cpp:4868` is gated to NNI moves AND only fires when
`MKPRIME_ECO_PARTIAL_CL=1`. Default unset → ecoCL is dormant. Not
applicable to the default workload.

**False sharing (hypothesis #4).** Single-threaded PT loop; not
applicable.

**Verdict.** T-015 as phrased is **REFUTED**. Per-chain-move time is
flat-to-within-3% between modes when normalised by chain count. The
remaining "~75 % of excess amplification" in T-014's framing was a
math artefact of comparing blind's super-linear scaling (driven by
fixed R overhead) to aware's near-linear scaling (driven by move time
dominating). There is no aware-specific PT locality cost to recover.

**Implication for T-014's framing.** T-014's "shared arena" candidate
for the remaining 1.69× excess is misconceived; there is no recoverable
excess. T-014's drift-gate fix is real and stands (~12 % wall at 4ch),
but it accounts for *all* the recoverable PT-specific cost on rodent,
not 25 %.

**Implication for production.** Hamilton viability is unchanged from
T-014's conclusion: aware nChains=4 PT on rodent runs at ~linear in
nChains. For a 30 k iter warmup budget under the v2 8 h walltime,
nChains=1 + nRuns=4 remains the viable production option. PT amplifies
walltime by ~4× as expected from doing 4× the work, not because of any
fixable locality cost.

**Filed.** T-015 in findings.md, status **REFUTED-AS-PHRASED**,
kind `[Investigation]`, priority P1 (closes a candidate that would
otherwise have absorbed a major refactor).

**Tests.** 228 / 228 pass (`filter="ecology|likelihood|mcmc"`), 0
failures, 6 expected skips. The instrumentation is bit-identity-safe
(only writes to local accumulators when env-gated).

**What's NOT done (out of scope for T-015 but worth recording).**
The orchestrator-level heap allocations in
`cpp_partition_log_likelihood_ecology` (`src/mcmc_ecology.cpp:1005,
1029`) and in `const_site_prob_*_eco_single` (`:854, :893, :924`) still
allocate ~800 KB per call uniformly per chain. This is a *T-007*
follow-up (uniform per-chain cost), not T-015 (excess at nChains>1).
Adding T-008-style workspace pre-allocation here would shave a fraction
of aware wall regardless of nChains. Estimated 5-10 % wall, but the
ecology orchestrator's call frequency dropped substantially after
T-010's partition cache, so the gain may be smaller in practice.
Filing as candidate **T-016**? — defer to the next /profile rotation.

**Cleanup.** Diagnostic instrumentation kept in-tree behind the
`MKPRIME_T015_DIAG=1` gate (inert by default, costs ~5 LoC of variable
declarations + ~25 LoC of accumulators + emission). Re-usable for future
PT-related investigations.

---

## Round 13 (T-018 — adaptive move-weight decay; salvaged) — 2026-05-21

**Target.** After T-015 closed the PT-locality investigation, the
user picked options 1 (threaded PT) and 2 (adaptive moves) from the
"how do we actually improve aware performance" conversation. Two
background subagents launched in worktrees. **Both died silently**
(known Opus-on-big-refactor failure mode — same pattern as the T-012
agent death in Round 10). The Sonnet T-018 agent made ~85 % progress
before dying. Salvaged from `mkp/.claude/worktrees/agent-a0b72a05663949b39/`.

**Salvage strategy.** Unlocked dead worktrees with `git worktree
unlock`, copied `R/RunMkPrime.R` and `tests/testthat/test-t018-
adaptive-moves.R` to `/tmp/`, removed the worktrees with `git worktree
remove --force`, deleted the stale branch refs, created a fresh
`t018-adaptive-moves` branch off `origin/t015-pt-locality`, applied
the copied files. The dev/t018_*.R scratch scripts were NOT copied
(throwaway agent-only files).

**What landed (R-side only, no C++ changes).** `.DecayLowAcceptMoves()`
in `R/RunMkPrime.R:4209` — multiplicative decay (default ×0.7) on free
moves with batch-level cold-chain acceptance < 2 % and ≥ 30 proposals
in the batch. Floor at 0.1 × initial weight per move. Re-normalises
to preserve the unpinned budget. Wired at `:1098-1107` inside the
warmup loop, after `.AdaptMoveWeights()`. Frozen at warmup-to-sampling
transition by the existing phase guard. `MKPRIME_ADAPT_DIAG=1` enables
per-decay `message()` to stderr. 120 LoC R + 13 unit tests.

**Verified.**
- 241 / 241 tests pass (228 baseline + 13 new T-018), 7 skips
  (1 new gated behind `MKPRIME_SLOW_TESTS=true`).
- Rodent 500-iter bench: aware 45.76 s (10.9 iter/s) vs. T-015
  baseline 43.56 s (11.5 iter/s) — within ±5 % noise.
- The `[T-018 decay]` diagnostic did NOT fire on the rodent run.
  Confirmed function is exported and call site executes; the
  default thresholds (n_min=30 proposals, accept_floor=2 %) simply
  don't trigger on rodent with warmup=100-200 iter. Low-weight
  moves don't accumulate 30 proposals; high-weight moves accept
  above 2 %.

**Honest framing.** T-018 is a defensive correctness improvement
that lands the infrastructure for acceptance-aware move-weight
adaptation. On rodent 500-iter it's a no-op (bench null). The
production benefit is conditional: long warmups (≥ 2000 iter) on
harder datasets where pathological low-acceptance patterns can
develop in slice or branch-length moves. Filing with the bench
null openly disclosed; if a future workload triggers decay we
can revisit.

**Filed.** T-018 in findings.md, status APPLIED, kind [Optimise],
priority P2 (because bench impact on the production target is null;
defensive infrastructure rather than wall recovery).

**Sibling task T-017 (threaded PT).** Opus agent died with near-zero
progress (single chore commit, no threading work). Same death pattern
as T-012 Opus agent. Re-launching with the same brief is likely to
fail the same way. User opted to schedule a 90-min cron wakeup to
launch a Plan-then-action sequence — a Plan subagent first scopes
the R-API thread-safety landscape (read-only, lower death risk),
then a fresh implementation agent works against a concrete design.

**Cleanup.** Dead agent worktrees removed (`git worktree remove
--force`), stale branches deleted (`git branch -D t017-pt-threaded
t018-adaptive-moves`). Salvaged files in /tmp left for now (small).

last_focus: 18


---

## Recovered from `to-do.md` (2026-09-18)

`to-do.md` was retired: its specific tasks became GitHub issues, and its
standing-task rows were rotation records, not work. The accumulated round
history of the **Performance profiling** standing task is preserved verbatim below. It predates
the Discussions migration and is history, not a queue.

**Standing: Performance profiling.** Profile the compiled MCMC hot path using VTune (see `r-package-profiling` skill) or `bench::mark()` microbenchmarks. Identify the current top hotspot after OPP-1–6. Check whether `pruning_jc_flat` / `pruning_jc_acrv_flat` show further vectorisation opportunities, whether chain-swap overhead is visible at scale, or whether R↔C++ boundary crossings dominate for small datasets. File any actionable findings as new `M-nnn` tasks. When completed, record the focus and key finding in Notes and reset to OPEN. Priority: same dynamic rule as S-RED. | Last run: 2026-03-31 round 5 (B, M-166). Focus: **VTune re-profile post M-156/M-157/M-158/M-164 optimizations.** Sun2018 (54 taxa, 225 all-trans, nCat=6), 5000 iters, 337s CPU. SW sampling (no admin). Compared to round 4 baseline (`vtune-out/`). **Key findings:** (1) `exp()` → `fast_neg_exp` (M-156): `_expl_internal` dropped from 11.7% to 0.1%; `fast_neg_exp` at 1.8% total — net ~10% CPU savings, the dominant improvement. (2) `pruning_jc_acrv_persite` share rose 59.1%→75.7% (absorbs ex-exp time). (3) `constant_site_prob_jc` (ascertainment) rose 6.6%→15.2% — now the clear #2 bottleneck. (4) `pruning_jc_acrv_flat` dropped 3.9%→0.6% (M-157 fused ascertainment, M-158 SPR partial CL, M-159 cache-aware scheduling). (5) M-164 pre-filter effect not separately quantifiable — multiple changes between baseline and current; per-iteration sweep cost roughly unchanged (pre-filter overhead may offset savings at this dataset scale). (6) Combined Gibbs sweep (persite+ascertainment+exp-like) = 92.7% of CPU (up from 77.4%). **Next bottleneck: `constant_site_prob_jc` at 15.2%** — batched ascertainment would give biggest remaining improvement. Results saved in `vtune-out-m166/`. Prev: round 4 (B) — baseline profile; filed M-156 fast_exp opportunity.

## Measure-first pass: #22, #11, #18 — 2026-09-25

Cloud container (Intel Xeon @ 2.10 GHz, 2 cores used), package built from
`main` @ `1ad163a` plus the #145 change (no effect on these paths). Workload:
Sun2018 (54 taxa, 225 trans chars, 3 partitions at kObs 2/3/4), NJ start
tree, `nCat = 4`. Drivers in `dev/profiling/drivers/`, all runnable with
`TREESEARCH_SRC=<TreeSearch checkout>` when TreeSearch is not installed.

**#22 — k' sweep bound past K: already fixed (#66, `1638678`).**
`22_kprime_sweep_bound.R` counts candidates per character with
`kprime_sweep_candidates()` on the plain-geometric arm (the only arm
truncated at K) over K in {10, 30, 100, 200} x p in {0.01, 0.1, 0.5}: zero
characters are evaluated past `nEff = K - kObs + 1` anywhere on the grid. The
per-partition bound #22 proposes is `partKoMax` in
`compute_per_kprime_log_lik`. For scale, one sweep costs 0.044 s at K = 30
and 3.1-3.5 s at K = 256 (the pre-#66 enumeration). Nothing left to land.

**#11 — `cache_total_loglik` unit re-scan: [AT-LIMIT].**
`11_cache_units_rescan.R` times the `nParts x nUnits` scan against a
precomputed per-partition index at Sun2018's 3 partitions: 0.006 us saved per
call at 1 unit per partition, 0.055 us at a pessimistic 10 distinct k' per
partition. The same call does at least one `constant_site_prob_jc` pass
(12.5 us at k = 2, including R call overhead), so the scan is < 0.5 % of the
call and far less of an evaluation (the scan is skipped under
`coding = "none"`).

**#18 — `singleton_site_prob_jc` per-call buffer: [AT-LIMIT].**
`18_singleton_alloc.R` compiles the kernel twice (fresh zero-filled vector per
call vs caller-owned scratch with only tips re-zeroed), 21 paired
replicates of 2000 calls: 75.9 vs 75.6 us per call, paired delta median
-0.04 us (IQR -2.8 to 3.0 us). Inside noise. The allocation is ~93 KB and the
pruning pass is O(nEdge * nTip * k) over it, so the malloc and zero-fill are
lost in the arithmetic.

**Lead for the next ascertainment round (not filed).** The per-call buffer is
not the cost, the call count is: under `coding = "informative"` Sun2018 makes
~111 singleton calls per likelihood evaluation (one per binary
transformational partition plus one per distinct missing-data mask,
`asc_probs_masked`), each a full O(nTip) pseudo-character pruning pass. A full
evaluation costs 8.3 ms informative against 0.65 ms variable on the same tree,
so the masked uninformative-mass correction is ~90 % of an informative
evaluation. Batching masks into one pass is where a real saving would be.
