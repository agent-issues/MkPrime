# MkPrime — architecture

R orchestration + Rcpp C++ hot loop, modelled on StratoBayes.

## R-side responsibilities

- Setup: parse phyDat input, classify characters (neomorphic vs.
  transformational), build partitions, configure priors and proposals.
- Adaptation: scalar tuning (sigma, scale), per-move acceptance tracking,
  warmup phases (slice-width adaptation, scalar mixing).
- Convergence: native rank-normalized R-hat (Vehtari et al. 2021), FFT-based
  ESS, tree-topology pseudo-ESS, stopping rules.
- Tempering: heated chains, swap proposals between chains.
- I/O: Tracer-compatible TSV streaming, Newick tree writes, checkpoint files.
- GUI: Shiny module under `R/BayesianModule.R` for TreeSearch integration.

## C++-side responsibilities

- Felsenstein pruning under JC(k') with analytical `P(t) = exp(Qt)`.
- Per-partition likelihood, conditional likelihood caching with
  dirty-flag invalidation (`node_cl_cache.h`).
- MCMC inner loop: propose → evaluate → accept/reject, batch runner.
- Tree moves (SPR, NNI, TBR, pSPR) and Gibbs topology moves.
- Scalar/branch proposals (scale, Dirichlet, beta-simplex, slice).
- Among-character rate variation (lognormal, discrete categories).
- Ascertainment correction (constant-site / singleton-site).
- Relabelling correction for Mk' (avoids model-graph cycles).

## Key design decisions

1. **Unrooted trees.** Stationary root frequencies make all models —
   including asymmetric MkN — time-reversible; root placement is irrelevant.
2. **Model k'_i not u_i.** True states (k') avoid model-graph cycles where
   u depends on kObs.
3. **No kMax.** JC eigendecomposition is O(1) analytical; prior + relabelling
   correction naturally constrain k'.
4. **Per-character k'_i from shared hyperprior.** Pools information across
   characters of the same type.
5. **`coding = "variable"` first.** `"informative"` deferred to Phase 7.

## Three character models, one analysis

| Type | Model | Key feature |
|------|-------|-------------|
| Neomorphic       | MkN (asymmetric binary)  | Shared `rate_loss` across characters |
| Transformational | Mk' (infer true k)       | Per-character k'_i from shared hyperprior |
| Known            | Mk(k)                    | User-specified k, no k' inference |

Partitions exist per character type and per k' (or k) bin.

## Public headers

`inst/include/MkPrime/` holds the header-only kernels shared with downstream
packages: `transition.h` (JCRates, MkNRates, F81), `models.h`, `pruning.h`
(`Pruner<Model>`: dirty-flag CL cache, checkpoint/accept/reject), `acrv.h`,
`corrections.h`, `ascertainment.h`. The MCMC hot path (`mcmc_likelihood.cpp`,
`ascertainment.cpp`, `gibbs_partial_cl.h`, `node_cl_cache.h`) uses the
transition and ACRV kernels; the R-level `pruning_*()` functions are thin
`Pruner` wrappers. `NodeCLCache` (MCMC) is separate from `Pruner`: it is bound
to `McmcData`, k'-grouping and unrooted `TreeNav`.

## Shared file convention

`src/init.cpp` and `RcppExports.cpp` are auto-generated. After running
`Rcpp::compileAttributes()`, verify the package compiles and exports register
cleanly.

`src/Makevars.win` is **never committed** — it is reserved for transient
profiling/debug flag overrides.

## R file layout (high level)

Functional separation: data → model → mcmc-config → run → posterior →
diagnostics. See `mkp/AGENTS.md` for the full `R/` listing.
