# Profiling baselines — mkp

Reference timings; refreshed each round. `/profile regress` rebuilds the same drivers and flags any > 10 % regression vs the value here.

| Area | Driver | Wall (s) | Top hotspot | Hotspot self % | As-of |
|------|--------|----------|-------------|----------------|-------|
| step6 (legacy) | `data-raw/step6_profile_eg.R` | 94.3 | `.Call → run_mcmc_batch_cpp` | 99.27 | 2026-05-18 |
| 1 (Felsenstein pruning) | `dev/profiling/drivers/01_felsenstein_pruning.R` | 17.39 (baseline) / 16.45 (T-005) / **12.07 (T-005 + T-006)** | `pruning_jc_acrv_persite` 75 % (VTune baseline); after T-005+T-006: `persite_impl<0>` 4.4 s + kfixed<2..24> ~3 s combined | 75 → ~40 (VTune) | 2026-05-19 |
| 0 (Eco orchestrator wall) | `dev/profiling/drivers/11_aware_vs_blind_rodent.R` | aware 24.59 (baseline) → 21.40 (T-010, -13.0 %) → 20.47 (T-011) → **20.58 (T-012 partial-salvage, cache not yet consulted in production)** / blind 0.97 → 0.77 → 0.79 → 0.78 (200 iter, rodent MkNT, kEco=4, nTip=64, nChar=217) | `.Call → run_mcmc_batch_cpp` | 97.20 (Rprof, baseline) | 2026-05-21 (T-012) |
| 0 (Eco orchestrator per-call) | `dev/profiling/drivers/11c_per_call_cost.R` | aware 0.57 ms / blind 0.20 ms (1000 reps, nCat=1, post T-008+T-009); T-009 adds no measurable delta — dominant cost is psMix/pdMix mixture loop not K-state inner work | — (isolated likelihood call) | 2.85× per call | 2026-05-21 (post T-008+T-009) |
| 0 (Eco pruning direct, T-009) | `.PruningJcEcology` direct call, synthetic data | 0.241 ms/call baseline, 0.241 ms/call T-009 (k=4, nChar=150, nTip=60, kEco=4, 5000 reps, 5 outer) | — | ~0 % delta | 2026-05-21 |
| T-008 (per_char batch) | `dev/profiling/drivers/12_per_char_alloc.R` | baseline 18.7 ms/batch (217 chars) → T-008 18.1 ms/batch (~4 % speedup); 0 heap allocs per sweep (was 7812) | — (isolated per-char helper) | ~4 % batch speedup; 0 heap allocs | 2026-05-21 |

Driver 1 yields ~8 iter/s under empirical_geometric + production move weights on Sun2018 (54 taxa, 225 chars, all transformational). PROFILING.md (2026-03-28) reported **2611 iter/s** on the same dataset under a reduced move schedule (NNI+SPR+BetaSimplex+kPrime+scalar only); the 330× slowdown is driven primarily by Gibbs k′ sweep (29.5 % weight) and Gibbs SPR / subtree-swap (1.3 % each) — these were explicitly excluded from the 2611 iter/s reference run.
