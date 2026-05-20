# v2 verification & follow-up runbook

To execute in order once the debug agent reports the MCMC hang is fixed.

## 0. Verify the fix

```r
setwd("C:/Users/pjjg18/GitHub/mkp/.claude/worktrees/ecology-aware")
devtools::load_all(".", quiet = TRUE)
library("TreeTools")
source("inst/ecology/simulations/sim3-helpers.R")
source("inst/ecology/simulations/sim3-simulate.R")
set.seed(1); tr <- .BuildConvergentTree(); eco <- .ConvergentEcology(tr)
edgeEc <- .AssignEdgeEcology(tr, eco)
z <- matrix(0L, 30L, 2L); z[1:10, 2] <- 1L
dat <- .SimulateMkPrimeEcology(tr, edgeEc, z, phi = 4,
  type = c(rep("neomorphic", 10), rep("transformational", 20)),
  baseRate = 0.5)
pd <- MatrixToPhyDat(dat)
mkd <- MkPrimeData(pd, ecology = eco[tr$tip.label], neomorphic = 1:10)
res <- RunMkPrime(mkd, tree = tr,
  model = MkPrimeModel(ecologyAware = TRUE, kPrimePrior = "geometric",
                       coding = "variable"),
  mcmc = MkPrimeMCMC(nIter = 500L, nChains = 1L, nRuns = 1L,
                     thin = 100L, treeThin = 100L,
                     minWarmup = 200L, maxWarmup = 400L,
                     logFile = NULL, checkpointFile = NULL))
# Should complete in ~20s and produce a finite log_likelihood trace.
```

## 1. Sim 1 (null) re-run, ~12 min

Validates that under H0 the v2 hyperparameters stay near prior.

```bash
Rscript inst/ecology/simulations/sim1-null.R
```

Pass criteria:
* `phi` posterior CI includes 1, mean ≈ 1
* `pi0` posterior CI overlaps 0.7 (prior mean), mean ≈ 0.7
* `theta` posterior near 0.5 (Beta(2,2) prior with no signal data)
* No NaN/Inf in any trace

## 2. Sim 2 (recovery), ~12 min

```bash
Rscript inst/ecology/simulations/sim2-recovery.R
```

Pass criteria (v2 should improve over v1 here):
* `pi0` posterior well below 0.7 (recovers that ~half of cells have z≠0)
* `phi` posterior CI excludes 1, mean ≈ 3 (v2 should fix the skewed-edge
  degeneracy that pinned v1 phi at 1)
* `theta` posterior near 1 (all effect cells are z=enc in this sim)

## 3. Sim 3 v1/v2 comparison head-to-head, ~70 min

```bash
Rscript inst/ecology/simulations/sim3-v1v2-compare.R 200000
```

This is the headline v2 sanity check. Critical pass criterion:
* **aware logL median ≥ blind logL median**. In v1, aware was ~76 nats
  worse than blind (rate-time ridge). v2 should close that gap because
  the ridge is gone.
* `tree_length` distribution should be tight (IQR similar to blind), not
  blow up to 1000+ as in v1
* `pi0` posterior < 0.7 (detects the ecology signal)
* `phi` posterior CI excludes 1

## 4. Sim 3b (negative control), ~70 min

```bash
Rscript inst/ecology/simulations/sim3b-mcmc.R
```

Pass criteria:
* `pi0` posterior near prior (no spurious effect)
* `phi` posterior CI includes 1
* Bipartition support: blind ≈ aware (ecology layer doesn't bias)

## 5. Multi-rep Sim 3 (paper-quality), local pilot

```bash
Rscript inst/ecology/simulations/sim3-multirep.R 4 20260601 100000
```

4 replicates × 100k iter ≈ 4 hours total. Look for consistent across-rep
P(true) > P(wrong) and aware ≥ blind logL gap.

## 6. Multi-rep Sim 3, HPC

20+ replicates × 200k iter. Dispatch via hamilton-hpc skill. Headline
figure: paired blind/aware P(true)/P(wrong) across replicates.

## 7. Rodent empirical run

```bash
Rscript inst/ecology/scripts/rodent-ecology-v2.R ~/downloads/mbank_X24848_2026-5-9-1135.nex 500000
```

~3 hr locally or ~1 hr HPC. Outputs:
* `inst/ecology/scripts/rodent-ecology-v2-result.rds`
* `inst/ecology/scripts/rodent-ecology-v2-zPost.csv`
* `inst/ecology/scripts/rodent-ecology-v2-heatmap.pdf`

## 8. Vignette finalisation

* Fill in actual sim numbers in `vignettes/ecology-details.qmd`
  simulation evidence section (sim1-3 + 3b results)
* Fill in heatmap and top-20 characters in
  `vignettes/rodent-ecology.qmd`
* Render: `quarto::quarto_render("vignettes/ecology-details.qmd")` and
  `quarto::quarto_render("vignettes/rodent-ecology.qmd")`

## 9. Tests

Update v2-incompatible tests in `tests/testthat/test-ecology-*.R`
(zMatrix dimensions, theta presence, refEcology, γ normalisation
reference). Add the "γ → 1 reduces to blind when π₀ = 1" regression
test (`v2-redteam.md` Concern 2 — verifies ascertainment correction
is γ-normalised).

## 10. Performance pass (low priority)

See `inst/ecology/notes/v2-gibbs-z-perf.md` for the gibbs_z γ-cancellation
optimization. ~3× speedup on the Gibbs sweep, second-order win
compared to fixing the mixing.

## Paper figures pipeline

Once steps 1-7 complete, generate paper-ready figures:

```bash
Rscript inst/ecology/simulations/sim-figures.R
```

Produces:
* `sim3-headline.pdf` — single-rep blind vs aware
* `sim3b-negctrl.pdf` — matched null
* `sim3-multirep-n20.pdf` — multi-rep paired comparison

The rodent heatmap (panel 3 in `rodent-ecology.qmd`) is the empirical
headline figure.
