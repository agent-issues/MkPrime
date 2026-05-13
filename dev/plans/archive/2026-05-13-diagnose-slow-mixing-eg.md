# Plan: diagnose slow mixing under `empirical_geometric` (Opus, 2026-05-12)

Symptom: 26 Hamilton tasks (array 17140607, heat=0.1) at 2h: median minESS=27,
max 73, 0/26 reach 100. logP swings 30-50. Per-sample ratio plateaus ~5e-5.
mh_logit_p is working (38% accept, 55% scheduler weight); the bottleneck is
downstream.

## Order: cheap first, fix last

### Step 1 — Decompose logP into logLik vs logPrior  (30 min)

**Highest payoff diagnostic.** Tells us *which* term is swinging.

- Add per-sample dump of `state->logLik` and `state->logPrior` separately
  to the cold-chain sample writer (R/RunMkPrime.R ~L2013).
- Plot traces for one stuck run from array 17140607.
- Decision:
  - **logLik swings**, prior quiet → topology / k' moving (go to Step 4)
  - **logPrior swings**, logLik quiet → k' moving along a flat-likelihood
    ridge driven by the prior alone (go to Step 5)
  - **Both swing, anti-correlated** → joint (topology, k') multimodality
    (go to Step 6, joint move design)

### Step 2 — Per-move acceptance + proposal counts  (15 min)

Already tracked in `chain_acceptance` (R/RunMkPrime.R ~L2052).

- `mh_logit_p` healthy at 20-40%? <5% → scaleTuning too large; >60% →
  proposals too small.
- `gibbs_kPrime` acceptance ~1.0 *but k' doesn't actually change* → full
  conditional extremely peaked, "k' locked" failure mode.
- `block_kPrime` low acceptance → k' values jointly constrained by topology
  (corroborates joint-multimodality hypothesis).

### Step 3 — Swap acceptance on tempering ladder  (10 min)

Pull `swap_rates` from result. Target 20–40% per adjacent pair.

- <5% or >80% → PT broken, fix ladder before anything else (PT failure
  invalidates downstream reasoning).

### Step 4 — Topology ESS vs scalar ESS  (1 h)

`ConvergenceDiagnostics(posterior, trees = TRUE)` gives `treeEss`.

- If topology ESS / sample << minESS / sample: tree-mixing problem hidden
  by scalar-focused minESS. Raise SPR/pSPR weights under `empirical_geometric`
  — `gibbs_kPrime`'s 55% may crowd them out.
- If topology ESS healthy but scalar minESS low: bottleneck is a specific
  scalar. Per-parameter ESS in diagnostics reveals which.

### Step 5 — Is the ESS estimator the problem?  (1 h)

`R/ess.R` uses Geyer + Vehtari; pathological on integer slow-drifting traces.

- Take one cold trace of `p`, thin by 1/10/100, recompute `.Ess()`, plot.
- Estimate stable across thinnings → real ESS roughly as reported.
- ESS falls with thinning → autocov tail truncated too eagerly, real ESS
  *higher* than reported (we're pessimistic).
- ESS rises with thinning → long-range autocorrelation missed, *worse* than
  reported.
- Also: which scalar dominates minESS? May be `betaScale` or one rate.

### Step 6 — Wall-time profile of gibbs_kPrime  (1-2 h)

Scheduler allocates 20% by *proposals*, not wall time. If gibbs_kprime_sweep
costs >50% of cycles, it's starving everything.

- VTune: `amplxe-cl -collect hotspots -duration 30 -- Rscript run_one.R`
- If gibbs_kprime_sweep_impl >50% of cycles: reduce its weight from
  `nTrans` to `nTrans / 4` at R/RunMkPrime.R L2788.

### Step 7 — Prior-data mismatch sanity check  (30 min)

- Plot `empiricalNObs$body` against the 26 stuck datasets' kObs distributions.
- If 3+ datasets have kObs median in the tail of the empirical body, the
  prior fundamentally fights the data; no sampler fix will help. Need to
  re-derive `empiricalNObs` from a matching corpus.

## Quick wins to try in parallel

1. **Lock mh_logit_p weight floor at 5%** — confirm `scalarTypes` floor
   actually applies (R/RunMkPrime.R L2899).
2. **Raise starting scaleTuning for mh_logit_p** if Step 2 shows low accept.
3. **Longer warmup** under empirical_geometric — NJ start may be far from
   the (different) posterior region.
4. **Check partition cache thrash**: how many NNI moves use `nniInPlace`
   partial eval vs full recompute? gibbs_kPrime invalidates `nodeCL` every
   call.

## Design sketches if real fix needed

**A. Joint (p, k') update** (medium effort)
- New move case 31: propose p on logit, then resample k' from conditional
  given p_new (call `gibbs_kprime_sweep_impl` under proposed p), accept/reject
  the package using marginal likelihood after integrating over k'.
- Moves along the (p, k') ridge, breaking the rate-limiting correlation.

**B. Cached egLogPriorByK by p-bin** (cheap)
- `gibbs_kprime_sweep_impl` recomputes `egLogPriorByK[m]` every call
  (mcmc.cpp L3386-3425). Cache keyed on discretised p (~0.01 bins),
  invalidate on bin crossing. ~5-10% throughput win.

**C. Joint (k', topology) SPR** (large; only if Step 1 anti-correlated)
- Weighted SPR that rescales k' on the moved subtree's affected chars.

## Critical files

- `src/mcmc.cpp` (cpp_log_prior emp_geom: L269-330; case 30: L4108-4225;
  gibbs_kprime_sweep: L3340-3600)
- `R/RunMkPrime.R` (move scheduling, .BuildMoves, .kMoveTypes)
- `R/ess.R`, `R/Convergence.R` (ESS estimator)
- `R/MkPrimeMCMC.R` (warmup / tuning defaults)
