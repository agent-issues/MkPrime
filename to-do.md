# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

---

## Phase 8: MCMC performance — C++ inner loop

The MCMC inner loop currently runs entirely in R. Every iteration incurs
R function-call overhead, list allocation/copy, GC pressure, and redundant
computation. This phase ports the hot path to C++, using TreeTools C++
headers where available.

### 8a: Quick R-side wins (no C++ changes)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-055 | P1 | ASSIGNED (B) | **Eliminate redundant `ape::reorder.phylo` in `MkpLogLikelihood()`.** The tree is *always* in postorder (invariant maintained by `.InitState()` and all topology proposals). The `ape::reorder.phylo(tree, "postorder")` call at likelihood.R:63 is pure waste — called on every single likelihood evaluation. Remove it; add a documented invariant comment. Also remove the `inherits()` and `match.arg()` validation from the hot-path code by creating an internal `.MkpLogLikelihood()` that skips checks, called from `.DoMove()`; the exported `MkpLogLikelihood()` retains validation for user-facing use. |
| M-056 | P1 | ASSIGNED (B) | **Pre-compute move weight vector.** `vapply(moves, [[, numeric(1), "weight")` is recomputed inside the innermost loop (per-chain, per-run, per-iteration) at RunMkPrime.R:135. Weights are constant except during adaptation (every 200 iterations). Compute once before the main loop and update only inside the adaptation block. |
| M-057 | P2 | OPEN | **Replace remaining `ape::reorder.phylo` calls with `TreeTools::Postorder`.** Affects proposals.R:171 (NNI), proposals.R:280 (SPR), RunMkPrime.R:65 (init), RunMkPrime.R:483 (resume), MkPrimeModel.R:104 (Fitch). `TreeTools::Preorder()` benchmarks at ~44us vs ape's ~70us. `TreeTools::Postorder()` is ~94us (still builds postorder from preorder), so use `Preorder()` where possible and adapt the C++ pruning to accept preorder edge ordering if beneficial. Add `TreeTools` to `Imports` in DESCRIPTION (already in `Suggests`). |

### 8b: Port proposals to C++

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-058 | P2 | OPEN | **C++ NNI proposal.** Port `ProposeNni()` (proposals.R:124-177) to C++. Current R implementation does: find internal edges via `which()`, pick random internal edge, find children/siblings via `which()`, swap two edge endpoints, call `ape::reorder.phylo`. In C++: operate directly on integer parent/child vectors, use TreeTools `preorder_edges_and_nodes` header for reordering. Return new edge matrix + relative branch lengths + logHastings. Avoids R list allocation and `ape::reorder.phylo` overhead per NNI proposal. |
| M-059 | P2 | OPEN | **C++ SPR proposal.** Port `ProposeSpr()` (proposals.R:202-289) to C++. The R version's `.Descendants()` BFS (lines 294-307) uses `c()` concatenation in a while loop — classic R antipattern. Use TreeTools `descendant_edges` C++ header or write a direct BFS on integer arrays. The rest of the SPR logic (prune node, suppress, regraft, reorder) is straightforward edge-index manipulation. Return new edge matrix + relative branch lengths + logHastings. |
| M-060 | P2 | OPEN | **C++ BetaSimplex proposal.** Port `ProposeBetaSimplex()` (proposals.R:24-80) to C++. Currently allocates a new vector per call and computes two `dbeta()` evaluations. For branch-length simplex proposals on large trees (nEdge can be hundreds), this vector allocation is significant. C++ version modifies in-place and returns logHastings. |

### 8c: C++ MCMC inner loop

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-061 | P2 | OPEN | **C++ `.DoMove()` — propose/evaluate/accept cycle.** This is the single highest-impact change. Currently `.DoMove()` (RunMkPrime.R:754-829) is called millions of times and each call: (1) copies the full state list (`proposed <- state`), (2) dispatches through R `switch()`, (3) calls the proposal function in R, (4) calls `LogPrior()` in R, (5) constructs a temporary tree and calls `MkpLogLikelihood()`, (6) returns a new R list. On rejection (~70-80% of calls), all allocations are garbage. A C++ `do_move()` would hold state in a C++ struct, call C++ proposals directly, call the existing C++ pruning functions, and return only accept/reject + updated state — eliminating all per-iteration R overhead. Depends on M-058, M-059, M-060. Also requires porting `LogPrior()` to C++. |
| M-062 | P2 | OPEN | **C++ outer iteration loop.** Once `.DoMove()` is in C++, port the per-chain/per-run/per-iteration triple loop (RunMkPrime.R:127-198) to C++. The C++ loop handles: move selection (weighted sampling), per-chain state updates, chain swap proposals, acceptance tracking. Call back to R only for: progress display (every N iterations), sample storage (every thin-th iteration post-warmup), adaptation (every 200 iterations during warmup), and checkpointing. This eliminates R interpreter overhead from the hot loop entirely. Depends on M-061. |
| M-063 | P3 | OPEN | **C++ state struct with pre-allocated workspace.** Design a `McmcState` C++ struct holding: edge matrix (parent/child int vectors), branch lengths, all scalar parameters, kPrime vector, cached log-likelihood, cached log-prior, and pre-allocated conditional-likelihood workspace (currently `std::vector<std::vector<double>>` reallocated per likelihood call in likelihood.cpp). The CL workspace should be allocated once at MCMC init and reused across iterations, avoiding millions of heap allocations. |

### 8d: Likelihood computation efficiency

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-064 | P3 | OPEN | **Partial likelihood recalculation for parameter-only moves.** When only a scalar parameter changes (tree_length, rate_loss, rate_log_sd, rate_neo, p), the tree topology and tip states are unchanged. For `tree_length` and `rate_neo` moves, only branch lengths change — the transition probability matrices change but the tree structure doesn't. For `rate_log_sd` moves, only the ACRV rate categories change. For `kPrime` moves, only one partition sub-group changes. Currently `MkpLogLikelihood()` recomputes everything from scratch. Implement move-type-specific partial recomputation: cache partition likelihoods and only recompute the affected partition(s). |
| M-065 | P3 | OPEN | **Avoid re-extracting parent/child/edgeLength from tree every likelihood call.** `MkpLogLikelihood()` extracts `tree$edge[,1]`, `tree$edge[,2]`, `tree$edge.length` as R vectors, then passes them to C++. In the C++ inner loop (M-061+), these would be held directly in the C++ state struct, eliminating this per-call overhead. This task is partly subsumed by M-063 but is listed separately for tracking. |

---

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing insufficient). |

---

## User issues (triaged)

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-066 | P2 | OPEN | **Post-hoc burnin tuning on `MkPosterior`.** Add a `SetBurnin(posterior, burnin)` function that trims samples/trees to post-burnin rows and stores the burnin fraction in the object. Add an `AutoBurnin(posterior)` function that scans possible burnin fractions and selects the one that maximises ESS while minimising PSRF (or a weighted combination). Update `print.MkPosterior()`, `ConvergenceDiagnostics()`, and the `plot` method to respect the stored burnin. By default the posterior should expose only post-burnin samples from all cold chains/runs. |
| M-067 | P2 | OPEN | **Vignette references via `inst/REFERENCES.bib` + Rd macros.** Move all vignette bibliography entries to `inst/REFERENCES.bib`. Add `Rdpack` to `Imports` and wire up `\insertRef{KEY}{mkp}` macros in Rd/roxygen where relevant. Check `../TreeSearch` for the Rdpack template and DESCRIPTION/NAMESPACE boilerplate. |
| M-068 | P2 | OPEN | **Rogue taxon suppression in hyoliths vignette.** In the summary section of `vignettes/hyoliths.qmd`, call `Rogue::QuickRogue()` on the posterior tree sample to identify and prune rogue taxa before computing and displaying the consensus tree. Show the rogue taxon list and compare consensus stability before/after pruning. |

---

## Phase 6b: TreeSearch integration (deferred)

*Not yet broken into tasks.*

Planned work:
- TreeSearch GUI integration hook ("Bayesian (Mk')" mode in EasyTrees)

## Posterior & analysis improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-066 | P1 | OPEN | **Automatic burnin selection for MkPosterior.** The posterior object should support tunable burnin that maximizes ESS while minimizing PSRF. Add a `burnin` parameter (number of post-warmup samples to discard) with an automatic selection method that optimizes the ESS/PSRF trade-off. The `MkPosterior` accessors (samples, trees, summary, etc.) should default to returning only post-burnin samples from all runs combined. Provide a method to adjust burnin after the fact (e.g. `set_burnin(posterior, n)` or a `burnin` argument to accessor functions). Currently, only warmup is discarded; additional burnin is common practice when chains take time to find the typical set even after warmup ends. |
| M-067 | P2 | OPEN | **Vignette and package documentation references via inst/REFERENCES.bib + Rdpack.** Create `inst/REFERENCES.bib` with all cited references (Sun et al. 2018, Lewis 2001, Xie et al. 2011, etc.). Set up Rdpack for `\insertRef{}` in roxygen docs. Update `hyoliths.qmd` to use `bibliography: ../inst/REFERENCES.bib` (or copy to vignettes/). Check `../TreeTools/` for a working template of this setup. Add `Rdpack` to Imports in DESCRIPTION. |
| M-068 | P2 | OPEN | **Rogue taxon suppression in hyoliths.qmd consensus tree.** In the tree summary section of `vignettes/hyoliths.qmd`, use `Rogue::QuickRogue()` to identify rogue taxa, then exclude them from the consensus tree (e.g. `Consensus(trees[, -rogues])` or `ape::drop.tip(consensus, rogues)`). Display which taxa were identified as rogues and show the cleaned consensus. Add `Rogue` to Suggests in DESCRIPTION. |
