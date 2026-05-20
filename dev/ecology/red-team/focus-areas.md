# Red-team focus-area rotation — MkPrime ecology-aware extension

Project: `MkPrime` (worktree `worktree-ecology-aware`). Rotation tuned to
the ecology-aware NT model added on this branch. Recent bug pattern:
silent accumulator drift, boundary handling in priors, move dispatch
inconsistencies after refactors. Highest-leverage uncovered areas come
first in rotation.

| # | Area                                            | Files in scope                                                                                                                  | Key questions                                                                                                                                                                                                       |
|---|-------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 1 | Ecology likelihood + prior math                 | `src/mcmc_ecology.cpp`; `src/mcmc.cpp` (`cpp_log_prior`); `R/MkPrimeModel.R` (`LogPrior`); `inst/simulations/ecology/sim3-simulate.R`; `vignettes/ecology-details.qmd` | Do simulator, R prior, and C++ prior describe the same model? `gamma_e` consistency. wEdge marginal computation. Boundary cases (theta, pi0 ∈ {0, 1}). Reference-ecology contract. Per-edge mixture application. |
| 2 | MCMC proposals + MH ratios                      | `src/mcmc.cpp` (move dispatch + scale_phi/pi0/theta + gibbs_z); `R/RunMkPrime.R` (`.kMoveTypes`, `.BuildMoves`)                  | Jacobians for log/logit transforms. Bactrian symmetry. gibbs_z conditional sampling. Full-recompute requirements. Move weights / `gibbsZEvery`. Per-ecology mode `phi[refE]` handling.                              |
| 3 | State synchronisation (logLik / logPrior)       | `src/mcmc.cpp` (all `state->logLik =` / `state->logPrior =`); `src/mcmc_ecology.cpp` (`cpp_log_likelihood_ecology`); `compute_full_loglik_at`; partition-cache early returns | Which moves write non-eco logLik into state in eco mode? `wEdgeDirty` honoured? Drift diagnostic correct in eco? Periodic resync correctness. Slice samplers via `slice_scalar_impl` bypass dispatcher. |
| 4 | **State init + `initOverrides` + `.InitState`** | `R/RunMkPrime.R` (`.InitState`, `.InitRun`, `.InitMcmcChain`, `RunMkPrime`); `src/mcmc.cpp` (`init_mcmc_state`, `fill_partition_cache`)        | Does R-computed `log_lik` from `.MkpEcologyLogLikelihood` match C++ `cpp_log_likelihood_ecology` byte-for-byte? Does `initOverrides` propagate every field correctly? `fill_partition_cache` early-return semantics. Phase-1 chain state shape. |
| 5 | Checkpoint / resume                             | `R/RunMkPrime.R` (`ResumeMkPrime`, `.LoadCheckpoint`, `.SaveCheckpoint`); `src/mcmc.cpp` (`get_mcmc_state`, checkpoint serialization) | Does resume preserve eco state (phi, pi0, theta, z)? Does it re-init `state->logLik` correctly? Move-tuning state restored? Random-seed continuity?                                                                  |
| 6 | Numeric stability + edge cases                  | `src/mcmc_ecology.cpp` (pruning + gamma_e); `src/mcmc.cpp` (log-sum-exp, MH ratios)                                              | Behaviour at phi → 0+ or large phi (~30). pi0 / theta near boundaries. Very small `rate_neo`. `tree_length` collapse. Wide `rate_log_sd`. Log-sum-exp underflow / overflow. NaN propagation.                       |
| 7 | Adaptive move-scale tuning during warmup        | `R/RunMkPrime.R` (warmup logic, scale tuning, ETA, convergence criteria); `src/mcmc.cpp` (scale persistence across batches)      | Why did warmup historically run for 7000 iter when `maxWarmup=1000`? Are eco moves tuned correctly? Stabilisation criteria sane? Does tuning push the chain out of high-posterior basins?                          |
| 8 | Simulator vs model alignment (full)              | `inst/simulations/ecology/sim3-helpers.R`, `sim3-simulate.R`, `sim3-multirep.R`, `sim3-helpers.R` (tree construction)            | Tree topology + edge ecology mapping. Seed reproducibility. Unit consistency between simulator and model (gamma_e normalisation, baseRate). `kEco ≠ 2` extension hazards.                                          |
| 9 | Data setup + input validation                   | `R/MkPrimeData.R` (`.ExtractEcology`, ecology arg parsing); `R/MkPrimeModel.R` (validator block)                                 | Does `ecology` arg accept all documented forms (column index / named vector / integer vector)? Validation tight on kEco, refEcology, taxon alignment? Coercion warnings?                                            |
| 10| Tests coverage + correctness                    | `tests/testthat/test-ecology-*.R`; `tests/testthat/test-m092*.R`; `tests/testthat/test-m140*.R`                                  | Do tests exercise the eco code paths actually run on the cluster? Are seeds pinned? What's NOT covered? Are assertions strict enough to catch the bugs we just found?                                              |
| 11| HPC dispatch + Hamilton scripts                 | `inst/hamilton/sim3-multirep-v3/*`; `inst/hamilton/sim3-truth-init/*`                                                            | Error handling. Resource correctness. Log management. DOS line endings. Idempotency. Env var handling. SSH/scp gotchas.                                                                                              |
| 12| Vignette math vs implementation                 | `vignettes/ecology-details.qmd`                                                                                                  | Does the documented model match the code? Are claims substantiated by sim results referenced? Are the boundary discussions consistent with the actual prior bounds?                                                |

## Rotation rationale (post-round-1 state)

Rounds 1-3 (likelihood/prior, moves/MH, accumulator) ran 2026-05-14 and
surfaced ~10 bugs, most fixed. Next areas chosen by leverage:

- **Area 4 (state init)** prioritised because (a) `initOverrides` was
  written this session and is fresh code; (b) advisor flagged "if
  state->logLik at iter 0 differs from fresh recompute, the bug is
  initialisation — separate from accumulator drift"; (c) recent 113-nat
  gap at sample 1 might originate there.
- **Area 5 (checkpoint/resume)** next because resume is the production
  path for multirep extensions and has never been audited.
- **Area 6 (numeric stability)** next because chains have visited phi=33
  and TL=346 — pruning under those is untested.
- **Area 7 (warmup tuning)** because the "warmup ran for 7000 iter"
  anomaly is unresolved; the chain leaves truth-init during warmup tuning.

After 7-12, rotation wraps to 1.
