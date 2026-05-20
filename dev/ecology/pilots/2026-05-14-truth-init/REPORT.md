# Truth-init diagnostic — verdict

**Date:** 2026-05-14
**SLURM jobs:** 17157462 (corrected; 10k iter, 1.8 min)
**Config:** Sim 3 v3 rep 1 (seed 20260602), aware model, chain init at
(phi=4, pi0=0.75, theta=0.999, z = truth, tree = truth, TL=13.5)
**Result file:** `inst/ecology/simulations/truth-init-results/`

## TL;DR

Truth doesn't hold. The chain leaves truth during warmup and reports
log_post values **300+ nats higher than reality** (stale log_post
accumulator). Once corrected, truth is genuinely the posterior peak —
the chain wandered into a *worse* region but doesn't know it.

Per the advisor decision tree: this is the **"mixing + sigmaPhi
tightening"** branch, not the "tighten pi0" branch. The model spec is
fine; the chain machinery is the problem.

## Numbers

### Logged values (chain.log)

| Iter  | log_post | log_lik | TL  | phi  | pi0  | θ₁    | rate_neo |
| ----- | -------- | ------- | --- | ---- | ---- | ----- | -------- |
| 7020  | **−4950**| −4727   | 169 | 25.9 | 0.23 | 0.995 | 1.89     |
| 10000 | −5363    | −5086   | 413 | 24.2 | 0.19 | 0.987 | 0.05     |
| median (last 50%) | **−5117** | −4859 | 280 | 33.1 | 0.15 | 0.986 | 0.10 |

### Recomputed (manual call to `.MkpEcologyLogLikelihood` + `LogPrior`)

| State                          | log_lik | log_prior | log_post |
| ------------------------------ | ------- | --------- | -------- |
| **Truth (sim params, true tree)** | −4617   | −457      | **−5074** |
| Sample 1 (using saved z, params, tree)         | −4707   | −549      | **−5257** |

**Drift: sample 1 is 183 nats below truth.** The chain reports it 125
nats *above* truth. Stale-log discrepancy: **+307 nats** at sample 1.

## What this means

1. **Model spec is fine.** Truth's log_post is the highest in the
   computed pair. Tightening priors to bias the chain *away from
   drifted regions* is not necessary — drifted regions are already lower
   posterior. (Beta(7, 3) pi0 prior is fine; data overwhelms it but
   posterior mode is correctly near truth.)

2. **Chain leaves truth and doesn't return.** With minWarmup=200,
   maxWarmup=1000 the chain warmed up for 7000 iter (max_iter override?
   needs investigation) and by first sample was already at phi=26.

3. **log_post accumulator is stale by 300 nats.** This is a serious
   diagnostic bug — chain.log readers can't tell convergence from
   chaos. Likely cause: incremental delta updates to log_post (`new =
   old + Δ`) accumulate floating-point error over many iter without
   periodic from-scratch recompute.

## Next moves (in priority order)

1. **Diagnose stale log_post.** Cheap: run a chain, recompute log_lik
   from scratch every K iter, log alongside the in-memory value. If
   drift grows monotonically, an accumulator periodic re-sync fixes
   it. If drift is one-shot, find the move that breaks the delta.
2. **Tighten sigmaPhi 1.0 → 0.5.** With phi=33 at 7σ rather than 3.5σ,
   the drifted phi mode becomes ~15 nats harder to reach — combined
   with mixing fix (#3), should close the gap.
3. **Add joint (log phi, log tree_length) Bactrian move.** Currently
   phi moves independently of TL; high-phi excursions stretch TL to
   compensate (sample 1: phi=26, TL=169). A correlated move would let
   the chain slide back to truth along the ridge.
4. **Re-run truth-init** under fixes — verify truth holds for at least
   5000 post-warmup iter.
5. **Re-run multirep** with the fixes — expect closing the CID_to_truth
   gap that the first multirep showed (blind 0.476 → aware 0.673).

## Files

- `inst/ecology/simulations/truth-init-results/` — raw outputs
- `inst/ecology/simulations/truth-init-diag.R` — chain.log dump
- `inst/ecology/simulations/truth-init-deep.R` — LogPrior decomposition
- `inst/ecology/simulations/truth-init-deeper.R` — warmup_trace + truth recompute
- `inst/ecology/simulations/truth-init-verify.R` — sample 1 recompute (the smoking gun)
- `inst/ecology/hamilton/sim3-truth-init/` — HPC dispatch scaffolding
