# MCMC moves / MH ratio red-team audit (2026-05-14)

**Summary.** The Bactrian kernel, log/logit transforms, and Jacobians on
`scale_phi` / `scale_pi0` / `scale_theta` are correct. The gibbs_z
sampler is mathematically sound (log-sum-exp shift, theta=1 edge case
safe). Three real bugs and one suspect: (B1) `scale_phi` can perturb
`phi[refEcology]` in `per_ecology` mode, breaking the identifiability
constraint; (B2) the phi prior in `cpp_log_prior` sums `dlnorm` over
**all** kEco entries including the reference index, compounding B1 and
inflating the prior by a constant in global-mode-disguised-as-ref runs;
(B3) `gibbsZEvery` is stored in `McmcData` but **never read** anywhere
in `mcmc.cpp` — the documented "Gibbs every N gens" gate does not
exist; gibbs_z fires whenever the weighted schedule draws it. (S1) The
2% weight floor list omits `scale_theta` and `gibbs_z`, so
`scale_theta` can be starved when total weights are large. The 7020-iter
"warmup" anomaly is explained by the Tuning phase (`thin × 100 × 4`
extra iter post-stabilisation) — not a bug, but easy to miss.

---

## Bugs

### B1. `scale_phi` mutates `phi[refEcology]` in per_ecology mode
`src/mcmc.cpp:4544-4554`.

```cpp
phiOldIdx = (nPhi == 1) ? 0 : static_cast<int>(R::unif_rand() * nPhi);
```

In `per_ecology` mode `state->phi` has length `kEco` and `phi[refE]`
is meant to be inert (always 1) — confirmed by `R/RunMkPrime.R:2782`
("phi[refEcology] is inert") and `src/mcmc_ecology.cpp:712` (`if (s ==
refEcology) gammaE[s] = 1.0`). The proposal picks any of the kEco
entries with equal probability, so 1/kEco of `scale_phi` proposals
perturb the reference and break the identifiability constraint
relative to the likelihood model. Symptom: the chain drifts because
the over-parameterised phi can compensate for true phi via the ref
index.

**Fix:** skip `refE` in index draw; or sample from `{0, ..., nPhi-1} \
{refE}` directly. In global mode (nPhi=1) this is a no-op.

### B2. Phi prior includes the reference index
`src/mcmc.cpp:410-412`.

```cpp
for (int i = 0; i < phi.size(); ++i) {
  lp += R::dlnorm(phi[i], 0.0, data.sigmaPhi, 1);
}
```

Loop iterates over all phi entries. In per_ecology mode this counts
`phi[refE]` (nominally =1, log-density `-log(sigmaPhi*sqrt(2π))`),
which is a constant offset in normal use — but interacts with B1: when
scale_phi flips the ref index, the prior moves and the MH ratio sees
phantom prior pressure. Even with B1 fixed, the loop should exclude
refE to match the model's identifiability story.

**Fix:** `for (int i = 0; i < phi.size(); ++i) if (i != refE) lp += ...`
(global mode keeps nPhi=1 and refE=-1 / ignored).

### B3. `gibbsZEvery` is plumbed but unused
`R/MkPrimeModel.R:101,189`; `src/mcmc_state.h:156`;
`src/mcmc_likelihood.cpp:2122,2223`; **zero matches** for `gibbsZEvery`
in `src/mcmc.cpp`.

The field is set from R, stored on `McmcData`, and never consulted by
the move dispatch in `run_mcmc_batch_cpp` or the case-36 branch. The
doc comment ("Number of MCMC generations between Gibbs sweeps")
describes behaviour that does not exist. In practice gibbs_z fires
purely on weighted draw — currently weight=10 in `.BuildMoves`
(`R/RunMkPrime.R:3074`). The header histogram `gibbs_z:1.0%` is the
draw rate, not 1/50.

**Fix:** either (a) honour the period by skipping case-36 proposals
that arrive before `iter % gibbsZEvery == 0`, returning a no-op; or
(b) delete `gibbsZEvery` from the model and doc the behaviour as
weight-based.

## Suspect

### S1. `scale_theta` and `gibbs_z` excluded from the 2% weight floor
`R/RunMkPrime.R:3099-3101`.

```r
scalarTypes <- c("scale", "int_walk", "gibbs_p", "scale_p", "logit_scale_p",
                  "slice", "kprime_alpha", "kprime_beta",
                  "scale_phi", "scale_pi0")
```

`scale_theta` is a 1D scale move identical in spirit to `scale_phi` /
`scale_pi0` but not listed, so it does **not** receive the weight floor.
When the kPrime moves dominate raw weight (large nTrans), scale_theta's
share collapses below the 2% backstop. The reported ~0.9% sampling rate
for scale_theta is consistent with floor non-application. With nIter=10k
this is ~90 theta proposals, plausibly under-sampled.

`gibbs_z` doesn't need a 1D scalar floor (it's a sweep), but absent
from the list is fine; flag here only because the same omission may be
accidental.

**Fix:** add `"scale_theta"` to `scalarTypes`.

---

## Looks correct

- **Bactrian symmetry / Jacobians.** Kernel at `src/mcmc.cpp:115-119`
  is symmetric about 0 (fair ±M coin), so q-ratio = 1 in transformed
  space. `scale_phi` uses `logHastings = log(mult)` (Jacobian for log
  → phi) — correct given the LogNormal prior on phi-space. `scale_pi0`
  / `scale_theta` use `log(p_new*(1-p_new)) - log(p_old*(1-p_old))` —
  correct logit Jacobian for target on pi0/theta-space (where the
  Beta prior lives).
- **Boundary handling.** Logit transform on pi0/theta keeps proposals
  in (0,1) by construction; explicit `!(pi0New > 0 && pi0New < 1)`
  guards at 4573, 4623 catch overflow ties.
- **Full likelihood recompute.** scale_phi falls through to the
  bottom-of-function path at `src/mcmc.cpp:4879-4883`, routing to
  `cpp_log_likelihood_ecology`. scale_pi0 / scale_theta call
  `cpp_log_likelihood_ecology` directly inside their case (4582-4585,
  4631-4634). No partial-CL shortcut for ecology, as expected.
- **gibbs_z categorical.** `src/mcmc.cpp:4023-4045` uses max-shift
  log-sum-exp; `theta=1` exactly gives `log_disc=-Inf` →
  `pCum[2]=exp(-Inf)=0` cleanly; only-all-`-Inf` is guarded at 4029.
  After-sweep full logLik/logPrior re-evaluation (4052-4061) prevents
  cell-by-cell drift.
- **Move id renumbering.** `case 30` (mh_logit_p) and `case 34-37`
  (eco) consistent between `src/mcmc.cpp` switch, `.kMoveTypes`
  (`R/RunMkPrime.R:3127-3156`), and comment at 4081-4085.
- **Rollback paths.** `phiOldIdx >= 0` test at 4672, 4703, 5016
  restores phi correctly on logHastings / prior / MH rejection.

## Could not verify in time budget

- **Adaptive scale tuning during warmup.** Whether scale_phi's
  `scaleTuning` grows too aggressively in early warmup, sending the
  chain out of the truth basin before mixing settles. Worth a quick
  trace of `chainScaleTunings[, scale_phi_idx]` across batches in the
  drifted run.
- **Warmup "early stop at 7020" claim.** Most likely explanation:
  `mcmc$warmup` (maxWarmup) was overridden by `MkPrimeMCMC` defaults
  (`R/MkPrimeMCMC.R:355-358`: NULL → `nIter/2`) or the
  Tuning-phase budget `mcmc$thin * 100 * 4` (≈6000 iter at thin=15)
  pushed the first sample past iter 7000. Not a moves/MH bug.
