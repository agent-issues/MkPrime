# MkPrime — validation datasets

## Small annotated datasets

In `../neotrans/inst/matrices/`, with Excel metadata files that classify
characters as Neomorphic / Transformational:

| Project | Taxa | Chars | Neo | Trans | Taxon |
|---------|------|-------|-----|-------|-------|
| 3832 | 10 | 27 | 12 | 15 | Canthyloscledidae |
| 950  | 12 |  9 |  ? |  ? | Hexacorallia |
| 4789 | 13 | 12 |  ? |  ? | Agelacrinitinae |
| 1271 | 25 | 33 |  ? |  ? | Amaltheidae |
| 3408 | 19 | 30 |  ? |  ? | Galericini |

Excel filename pattern: `Project{N}_{author}.xlsx`. Column
`"Character Pattern"` contains `"Neomorphic"` / `"Transformational"`.

## Primary benchmark

**Sun2018** — 54 taxa, 225 chars (all transformational), via the
`TreeSearch` package's bundled datasets. Used for:

- VTune profiling (15k iterations, ~90s CPU)
- M-155 Gibbs kPrime speedup benchmarks
- Warmup stabilisation (M-131) validation

## Topology-multimodal targets

**None available.** `dev/benchmarks/bimodal-target.R` builds and verifies a
conflicting-signal simulation (18 taxa, 300 characters; a focal clade grafted
into the left or right arm of a symmetric backbone, half the characters
simulated on each), intended as the ESJD proposal-scheduling benchmark.

**Its `valley` stage FAILS.** The pooled matrix is fitted best by neither
generating tree but by the compromise topology with the clade on the central
edge: measured barriers -14.7 nats (one-spine variant, MCMC over 25
placements), -10.0 / -23.9 nats (two-arm variant, MCMC), and negative in 20 of
21 proxy cells spanning branch length, state count and clade size. Committing
to the supported island gains only +3.0 nats from its own 150 characters while
costing +28.1 on the conflicting 150.

So do **not** use it as a mixing benchmark, and do not expect any simulation
under the inference model to be topology-multimodal: balance and bimodality pull
against each other. For a peaky target, prefer an empirical matrix already known
to be multimodal (Lakner et al. 2008's DS1-DS11 proposal-benchmark set; Whidden
& Matsen 2015 identify DS1 as two-peaked), run as `knownStates = 4` Mk(4).

## Hyoliths vignette

`vignettes/hyoliths.qmd` runs a full Sun2018 example end-to-end. Run output
files (`hyoliths_*.log`, `hyoliths_*.ckp`, `hyoliths_*_trees.nwk`) are
gitignored.
