# Ecology-aware MkPrime simulation suite

Validation simulations for the ecology-aware extension of MkPrime.
Each script is self-contained and writes a `*-result.rds` artefact and
one or more `*.log` chain files into the same directory. Run from the
package root with `Rscript inst/ecology/simulations/<script>.R`.

## Scripts

| Script | Purpose |
|---|---|
| `sim3-helpers.R` | Tree / ecology builders for the convergent-clades design (sourced by the others). |
| `sim3-simulate.R` | Forward Mk' simulator with ecology-driven per-edge rates (sourced by the others). |
| `sim3-pilot.R` | Parsimony pilot grid: identifies the Goldilocks `(nEco, nBase, phi, stem, root)` config where the convergence signal overwhelms blind Mk while leaving enough phylogenetic signal for the aware model to recover. |
| `sim1-null.R` | Null sim: no ecology effect; expect `pi0 ≈ 0.7`, `phi ≈ 1`, true tree recovered. |
| `sim2-recovery.R` | Recovery sim: half of transformational characters driven by `z = 1` at `phi = 3`; expect `|log phi|` away from 0, `pi0 << 0.7`. |
| `sim3-mcmc.R` | Headline single-replicate Sim 3: blind vs aware on convergent-ecology data at the Goldilocks config. |
| `sim3b-mcmc.R` | Negative control for Sim 3: same tree and tip ecology, but `z = 0` everywhere. Expect blind and aware to behave similarly. |
| `sim3-multirep.R` | Multi-replicate Sim 3 for paper-quality statistics. CLI args: `nReps seedBase nIter`. Local pilot at `nReps = 4`; HPC run at `nReps >= 20`. |
| `sim3-analyse.R` | Post-hoc analyser: takes a result file (default `sim3-mcmc-result.rds`), reports per-chain CID, bipartition support, and "closer to truth vs. closer to wrong" tally. |
| `sim-figures.R` | Headline PDF panels for the paper: single-replicate bar/box, multi-replicate paired comparison. |

## The convergent-clades design

Tree: four 4-tip clades $A, B, C, D$, true topology $((A, C), (B, D))$.
Tip ecology: $A_*, B_* \to 1$; $C_*, D_* \to 0$. Branch lengths tuned
(`stemBranch = 0.10`, `rootBranch = 0.15`, `tipBranch = 0.5`) so the
baseline character signal alone supports the truth under parsimony but
the convergent characters drive blind Mk to $((A, B), (C, D))$.

Goldilocks config: `nEco = 60` neomorphic + `nBase = 180`
transformational characters, `phi = 4`.

## Reading the results

Each `*-result.rds` contains:

- `tree`, `eco`, `z`, `charType`, `phi`, `data`, `edgeEco` — the
  ground-truth simulator inputs.
- `resBlind` and / or `resAware` — `RunMkPrime` outputs.

`sim3-analyse.R inst/ecology/simulations/sim3-mcmc-result.rds` prints
the headline statistics. `sim-figures.R` generates the corresponding
PDFs.
