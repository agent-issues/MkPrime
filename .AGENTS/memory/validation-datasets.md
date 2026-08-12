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

**Its `valley` stage FAILS: 7.2 nats against a 10-nat bar.** The islands are
real but shallow — on the balanced default, peak A -4541.3, peak B -4535.0, best
crossing topology -4548.4 (MCMC over all 25 clade placements). The barrier comes
almost entirely from the constructor's rejection-sampling balance step: drawn
without it, the pooled matrix is fitted *best* by the compromise topology
(-14.7 nats one-spine MCMC; -10.0 / -23.9 two-arm MCMC; negative in 20 of 21
proxy cells spanning branch length, state count, clade size and a capped focal
stem). Committing to the supported island then gains only +3.0 nats from its own
150 characters while costing +28.1 on the conflicting 150.

So do **not** treat it as a confirmed mixing benchmark, and do not expect a
simulation under the inference model to be reliably topology-multimodal: a 50/50
mixture of two trees is usually best explained by one intermediate tree. For a
peaky target, prefer an empirical matrix: Whidden &
Matsen (2015, Syst. Biol. 64:472) show topological peaks *do* occur in posteriors
from real data, and measure MCMCMC's ability to cross the valleys, using the
standard DS1-DS8 proposal-benchmark set of Lakner et al. (2008). Those are
nucleotide matrices of 27-67 taxa, so they would enter MkPrime as 4-state
characters with `knownStates = 4` — verify the peaks survive Mk(4) before
relying on them.

## Hyoliths vignette

`vignettes/hyoliths.qmd` runs a full Sun2018 example end-to-end. Run output
files (`hyoliths_*.log`, `hyoliths_*.ckp`, `hyoliths_*_trees.nwk`) are
gitignored.
