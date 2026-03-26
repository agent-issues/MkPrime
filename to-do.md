# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

---

## Phase 4: Tree search — COMPLETE

All Phase 4 tasks (M-022 through M-028) are done. See `completed-tasks.md`.

## Phase 5: Parallel tempering + convergence

### Dependencies

```
M-029 → M-030 → M-031
M-029 → M-032 → M-033 → M-034
M-029 → M-035
M-032 + M-033 → M-036
```

### Tasks

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-029 | P1 | ASSIGNED (A) | Parallel tempering core: multi-chain + temperature ladder |
| M-030 | P1 | OPEN | Chain swap proposals between adjacent temperatures |
| M-031 | P2 | OPEN | Adaptive temperature tuning (target swap acceptance) |
| M-032 | P1 | OPEN | Independent runs (`nRuns` outer loop) |
| M-033 | P1 | OPEN | Convergence monitoring (ESS + PSRF via coda) |
| M-034 | P2 | OPEN | Stopping rules (ESS/PSRF/time/iter thresholds) |
| M-035 | P2 | OPEN | Checkpointing (save/restore full MCMC state) |
| M-036 | P2 | OPEN | MkPosterior multi-run updates (combined summaries, diagnostics) |

### Task details

**M-029: Parallel tempering core**
- `nChains` parameter in MkPrimeMCMC (default 4)
- `heat` parameter (default 0.2) — temperature of hottest chain
- Temperature ladder: geometric spacing β_i = heat^((i-1)/(nChains-1))
  - Chain 1: β=1 (cold), Chain nChains: β=heat
- Heated MH acceptance: log_alpha = β*(logLik_new - logLik_old) + (logPrior_new - logPrior_old) + logHastings
- State stores unheated log_lik + log_prior; heated post computed on the fly
- Per-chain: state, acceptance tracking, tuning (independent adaptation)
- Only cold chain (β=1) stores samples
- All chains start from same initial state (temperatures cause divergence)
- Progress bar shows cold chain log_post
- Backward compatible: nChains=1 gives identical behavior to Phase 4

**M-030: Chain swap proposals**
- After each iteration cycle, propose swap between random adjacent chains
- Swap acceptance: min(1, exp((β_i - β_j)(logLik_j - logLik_i)))
- Swap = exchange full states between chains (cheap: pointer swap)
- Track swap acceptance rates per adjacent pair
- `swap_interval` parameter (default 1 = every iteration)

**M-031: Adaptive temperature tuning**
- During warmup, adjust `heat` to target 23–30% swap acceptance
- Only adjust if swap acceptance significantly off target
- Simple multiplicative: if swap rate too low, increase heat (hotter = easier swaps)
- Freeze temperatures after warmup

**M-032: Independent runs (nRuns)**
- `nRuns` parameter in MkPrimeMCMC (default 2)
- Each run = independent set of nChains tempered chains
- Different random starting states per run (random tree + perturbed params)
- Sequential execution (parallel later via future/furrr if needed)
- Each run contributes cold chain samples independently

**M-033: Convergence monitoring (ESS + PSRF)**
- Compute ESS per parameter from each cold chain (coda::effectiveSize)
- Compute PSRF across cold chains from different runs (coda::gelman.diag)
- Requires nRuns >= 2 for PSRF
- Monitoring interval: check every N iterations (default 1000)
- Console reporting of min ESS and max PSRF

**M-034: Stopping rules**
- Parameters: `max_iter`, `max_time` (seconds), `min_ESS`, `max_PSRF`
- Check at monitoring intervals (from M-033)
- Stop when: (min_ESS >= threshold AND max_PSRF <= threshold) OR max_iter OR max_time
- Report reason for stopping
- Default: no early stopping (only max_iter applies)

**M-035: Checkpointing**
- `checkpoint_file` parameter in MkPrimeMCMC (default NULL)
- Save full state (all chains, all runs, tuning, samples so far) to RDS
- Checkpoint at monitoring intervals
- `resume_from` parameter: restore from checkpoint and continue
- Handle version/config mismatches gracefully

**M-036: MkPosterior multi-run updates**
- Store per-run results (samples, trees, acceptance)
- Combined posterior: concatenate cold chain samples across runs
- summary() includes: ESS per parameter, PSRF per parameter, min ESS, max PSRF
- print() shows per-run and combined stats
- plot() shows traces per run (different colors)
- Convergence warning if PSRF > 1.05 or ESS < 200
