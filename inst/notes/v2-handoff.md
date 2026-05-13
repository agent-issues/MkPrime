# v2 ecology-aware implementation — handoff

State at end of this session. The v2 redesign (γ normalisation, reference
ecology, asymmetric slab) is approximately 80% implemented. All math
machinery compiles cleanly and likelihood evaluation works correctly on
realistic data; the MCMC inner loop hangs and needs targeted debugging.

## What works (verified)

* **Compilation:** clean build with `devtools::load_all(recompile = TRUE)`.
  Warnings only (unused legacy helpers).
* **MkPrimeData:** `mkd$refEcology` (tip-frequency-most-common) and
  `mkd$ecoTipFraction` populated correctly. Verified `kEcology = 2`,
  `refEcology = 0` on the Sim 3 convergent fixture.
* **LogPrior (R-side):** v2 asymmetric slab with per-ecology `theta_e`,
  Beta hyperprior on `theta_e`, prior on `pi0` and `phi`. zMatrix is
  `nChar × (kEco − 1)`.
* **MkPrimeModel:** new args `thetaAlpha`, `thetaBeta` (defaults 2, 2)
  validated, plumbed into returned list, printed in `print.MkPrimeModel`.
* **C++ rate computations:**
  - `gamma_e_compute(pi0, theta_e, phi) → π₀ + (1−π₀)(θ φ + (1−θ)/φ)`
  - `trans_rate_factor(z, e, refEcology, phi, gamma_e)`: 1 on ref;
    `μ(z, φ) / γ_e` on non-ref
  - `mkn_rates_for_state`: same; both gain and loss scale by 1/γ_e
* **Pruning kernels:** `pruning_jc_acrv_flat_ecology` and
  `pruning_mkn_acrv_flat_ecology` take `refEcology` and `gammaE`,
  iterate ecology columns with the (s == refEcology ? 1 : μ/γ) logic.
  zMatrix column lookup via `(s < refEcology) ? s : (s − 1)`.
* **Ascertainment:** `const_site_prob_*_eco_single` updated with same.
* **Likelihood orchestrator:** `cpp_log_likelihood_ecology` computes
  `gammaE[]` once per call, passes refEcology and gammaE through to all
  kernel and ascertainment calls.
* **R wrappers** `.PruningJcEcology`, `.PruningMknEcology`,
  `.ConstSiteProb*Ecology` take refEcology/theta/pi0 with sensible
  defaults.
* **Single-call likelihood:** `.MkpEcologyLogLikelihood(...)` returns
  finite log_lik (verified `-103.5` on a 10-char Sim 3 fixture).
* **Blind path (`ecologyAware = FALSE`)**: end-to-end smoke passes,
  unchanged.

## What's broken (one specific issue)

`RunMkPrime(..., ecologyAware = TRUE)` **hangs in the MCMC inner loop**.
Single-call likelihood is finite; chain init reaches `init_mcmc_state`
successfully; first iteration of the C++ MCMC loop never returns.

Diagnostic so far:

* Blind path runs normally → not a generic MCMC bug.
* Likelihood is finite on the initial state → not an `R_NegInf` death
  spiral on retry.
* The hang happens before any output ticker, suggesting it's in
  the very first move's processing.
* Likely culprits (in priority order):
  1. **`scale_pi0` (case 31) or `scale_theta` (case 33) likelihood
     recompute** entering a state where γ_e blows up (e.g. pi0 → 0 with
     phi >> 1, or theta numerically pegged to 0/1) and never returning
     finite log_lik → infinite retry loop in the move-adapter.
  2. **Move-weight scorer running with stale move tuning** for the new
     `scale_theta` move ID 33; if the move count for ID 33 is 0 or 1
     and the weight goes to zero/divergent, the scheduler may loop.
  3. **`gibbs_z_sweep_impl` race** when zCols = 1 (the typical
     binary-ecology case) — the per-cell loop may have an off-by-one
     in the `zCol` mapping.
* `alwaysAcceptTypes` for `scale_pi0` was set to TRUE in v1 (always
  accept; no likelihood factor). v2 makes pi0 affect γ_e, so this
  must be FALSE under v2. Search in `R/RunMkPrime.R` for the
  `alwaysAcceptTypes` definition — confirm scale_pi0 is removed.
  Same check for scale_phi.
* Move table entry for `scale_theta` (ID 33) may not be registered
  in the R move table at `R/RunMkPrime.R:~2960`. If the C++ code
  emits move 33 but R didn't allocate a name/weight slot, the
  scheduler could loop trying to find a runnable move.

## What's not done

* **`scale_theta` move table registration in R** (Step F partial)
* **`gibbs_z_sweep_impl` C++ update to use ecologyToZCol mapping**
  — currently iterates 0..kEco-1 which is wrong for v2; should
  iterate 0..kEco-2 (zMatrix columns) with the (s < refEcology ?
  s : s+1) reverse mapping back to ecology states. Almost certainly
  the source of the hang in combination with the gibbs_z move weight.
* **Tests** — `tests/testthat/test-ecology-*.R` still expect v1 z
  dimensionality and will fail. The smoke-test fixture also uses
  random matrices that produce -Inf likelihoods (orthogonal bug,
  predates v2).
* **Re-run Sim 1, Sim 2, Sim 3, Sim 3b** under v2 once MCMC works.
* **Rodent empirical run.**

## Suggested next-session debugging plan

The hang is somewhere in the C++ MCMC inner loop. Fixed so far in this
session:
- `R/MkPrimeMCMC.R`: ecology move names registered for validation
- `R/RunMkPrime.R`: `scale_pi0` removed from `alwaysAcceptTypes`
  (v1-isms; under v2 it's a full MH move)

Still to do — DEBUGGING:

1. **First, add Rcpp::Rcout printfs** to each ecology move case
   (30, 31, 32, 33) and the move dispatcher in `src/mcmc.cpp`.
   Rebuild. Re-run the minimum reproducer. The last printf before
   silence pinpoints the move that hangs.
2. **Likely culprits** (priority order):
   - `gibbs_z_sweep_impl` calls `per_char_log_lik_ecology` with the
     z vector but does NOT recompute wEdge between cells (correct —
     wEdge doesn't depend on z). However, per_char_log_lik_ecology
     might be calling pruning_*_ecology with an inconsistent
     state. Verify the kp argument is correct for neomorphic chars.
   - `scale_phi` (case 30) recomputes likelihood via the generic
     flow at ~line 4670 onwards. Verify that flow uses the
     ecology-aware code path: `cpp_log_likelihood_ecology` not the
     blind `cpp_log_likelihood`. There's a `bool eco = data->ecologyAware`
     check that should route correctly — verify by inspection.
   - `scale_pi0` (case 31) and `scale_theta` (case 33) — these are
     full MH with explicit `cpp_log_likelihood_ecology` calls.
     The internal state (theta, pi0) is rolled back on rejection;
     verify rollback logic for both.
3. **Register scale_theta in R move table** (Step F partial):
   `R/RunMkPrime.R` lines ~2953-2975 has the move list. Add:
   ```r
   list(name = "scale_theta", type = "scale_theta", target = "theta",
        weight = 2)
   ```
   And in the move ID map at ~3052-3054: `scale_theta = 33L`.
4. Once the chain runs, immediately re-run Sim 1 (12 min) and check
   that under v2 the aware logL ≥ blind logL on the same data —
   this is the key v2-vs-v1 sanity check.

Other partial work that may interact with the hang:
- `gibbs_z_sweep_impl` re-saves logLik/logPrior after the sweep
  (lines ~4032-4040); this recompute is full ecology likelihood.
  Should be fine but verify it's not entering a NaN loop.
- The advisor noted γ_e *cancels* in the Gibbs categorical
  ratio for a single cell — the current code computes `ll[v]`
  three times per cell via `per_char_log_lik_ecology` which DOES
  include γ_e effects. This is over-computing but should still
  give correct sampling — not a hang cause.

## Minimum reproducer

```r
setwd("C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware")
devtools::load_all(".", quiet = TRUE)
library("TreeTools")
source("inst/simulations/ecology/sim3-helpers.R")
source("inst/simulations/ecology/sim3-simulate.R")
set.seed(1)
tr <- .BuildConvergentTree()
eco <- .ConvergentEcology(tr)
edgeEc <- .AssignEdgeEcology(tr, eco)
z <- matrix(0L, 30L, 2L); z[1:10, 2] <- 1L
dat <- .SimulateMkPrimeEcology(tr, edgeEc, z, phi = 4,
  type = c(rep("neomorphic", 10), rep("transformational", 20)),
  baseRate = 0.5)
pd <- MatrixToPhyDat(dat)
mkd <- MkPrimeData(pd, ecology = eco[tr$tip.label], neomorphic = 1:10)
# This returns finite log_lik:
.MkpEcologyLogLikelihood(tr, mkd,
  kPrime = rep(2L, mkd$nChar),
  rate_loss = 1, rate_log_sd = 0, nCat = 1, rate_neo = 1, relabel = FALSE,
  phi = 1.0, zMat = matrix(0L, mkd$nChar, 1), magnitudeMode = "global",
  coding = "none", refEcology = mkd$refEcology, theta = 0.5, pi0 = 0.7)
#> [1] -103.4996

# This hangs:
model <- MkPrimeModel(ecologyAware = TRUE, kPrimePrior = "geometric",
                      coding = "variable")
mcmc <- MkPrimeMCMC(nIter = 20L, nChains = 1L, nRuns = 1L,
                    thin = 5L, treeThin = 5L,
                    minWarmup = 0L, maxWarmup = 0L, logFile = NULL,
                    checkpointFile = NULL)
tr2 <- Preorder(ape::rtree(mkd$nTip, tip.label = mkd$taxon_names))
tr2$edge.length <- rep(0.1, nrow(tr2$edge))
res <- RunMkPrime(mkd, tree = tr2, model = model, mcmc = mcmc)
```
