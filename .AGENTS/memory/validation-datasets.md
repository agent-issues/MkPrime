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

## Hyoliths vignette

`vignettes/hyoliths.qmd` runs a full Sun2018 example end-to-end. Run output
files (`hyoliths_*.log`, `hyoliths_*.ckp`, `hyoliths_*_trees.nwk`) are
gitignored.
