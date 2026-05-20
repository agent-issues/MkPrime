# Profiling baselines — mkp

Reference timings; refreshed each round. `/profile regress` rebuilds the same drivers and flags any > 10 % regression vs the value here.

| Area | Driver | Wall (s) | Top hotspot | Hotspot self % | As-of |
|------|--------|----------|-------------|----------------|-------|
| step6 (legacy) | `data-raw/step6_profile_eg.R` | 94.3 | `.Call → run_mcmc_batch_cpp` | 99.27 | 2026-05-18 |
| 1 (Felsenstein pruning) | `dev/profiling/drivers/01_felsenstein_pruning.R` | 17.39 (baseline) / 16.45 (T-005) / **12.07 (T-005 + T-006)** | `pruning_jc_acrv_persite` 75 % (VTune baseline); after T-005+T-006: `persite_impl<0>` 4.4 s + kfixed<2..24> ~3 s combined | 75 → ~40 (VTune) | 2026-05-19 |

Driver 1 yields ~8 iter/s under empirical_geometric + production move weights on Sun2018 (54 taxa, 225 chars, all transformational). PROFILING.md (2026-03-28) reported **2611 iter/s** on the same dataset under a reduced move schedule (NNI+SPR+BetaSimplex+kPrime+scalar only); the 330× slowdown is driven primarily by Gibbs k′ sweep (29.5 % weight) and Gibbs SPR / subtree-swap (1.3 % each) — these were explicitly excluded from the 2611 iter/s reference run.
