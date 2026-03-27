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

### 8b: Port proposals to C++

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-058 | P2 | DONE (C) | **C++ NNI proposal.** Port `ProposeNni()` (proposals.R:124-177) to C++. Current R implementation does: find internal edges via `which()`, pick random internal edge, find children/siblings via `which()`, swap two edge endpoints, call `ape::reorder.phylo`. In C++: operate directly on integer parent/child vectors, use TreeTools `preorder_edges_and_nodes` header for reordering. Return new edge matrix + relative branch lengths + logHastings. Avoids R list allocation and `ape::reorder.phylo` overhead per NNI proposal. |
| M-059 | P2 | DONE (C) | **C++ SPR proposal.** Port `ProposeSpr()` (proposals.R:202-289) to C++. The R version's `.Descendants()` BFS (lines 294-307) uses `c()` concatenation in a while loop — classic R antipattern. Use TreeTools `descendant_edges` C++ header or write a direct BFS on integer arrays. The rest of the SPR logic (prune node, suppress, regraft, reorder) is straightforward edge-index manipulation. Return new edge matrix + relative branch lengths + logHastings. |
| M-060 | P2 | DONE (C) | **C++ BetaSimplex proposal.** Port `ProposeBetaSimplex()` (proposals.R:24-80) to C++. Currently allocates a new vector per call and computes two `dbeta()` evaluations. For branch-length simplex proposals on large trees (nEdge can be hundreds), this vector allocation is significant. C++ version modifies in-place and returns logHastings. |

### 8c: C++ MCMC inner loop

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-061 | P2 | DONE (C) | **C++ `.DoMove()` — propose/evaluate/accept cycle.** This is the single highest-impact change. Currently `.DoMove()` (RunMkPrime.R:754-829) is called millions of times and each call: (1) copies the full state list (`proposed <- state`), (2) dispatches through R `switch()`, (3) calls the proposal function in R, (4) calls `LogPrior()` in R, (5) constructs a temporary tree and calls `MkpLogLikelihood()`, (6) returns a new R list. On rejection (~70-80% of calls), all allocations are garbage. A C++ `do_move()` would hold state in a C++ struct, call C++ proposals directly, call the existing C++ pruning functions, and return only accept/reject + updated state — eliminating all per-iteration R overhead. Depends on M-058, M-059, M-060. Also requires porting `LogPrior()` to C++. |
| M-062 | P2 | DONE (B) | **C++ outer iteration loop.** Once `.DoMove()` is in C++, port the per-chain/per-run/per-iteration triple loop (RunMkPrime.R:127-198) to C++. The C++ loop handles: move selection (weighted sampling), per-chain state updates, chain swap proposals, acceptance tracking. Call back to R only for: progress display (every N iterations), sample storage (every thin-th iteration post-warmup), adaptation (every 200 iterations during warmup), and checkpointing. This eliminates R interpreter overhead from the hot loop entirely. Depends on M-061. |
| M-063 | P3 | DONE (B) | **C++ state struct with pre-allocated workspace.** |

### 8d: Likelihood computation efficiency

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-064 | P3 | DONE (C) | **Partial likelihood recalculation for parameter-only moves.** When only a scalar parameter changes (tree_length, rate_loss, rate_log_sd, rate_neo, p), the tree topology and tip states are unchanged. For `tree_length` and `rate_neo` moves, only branch lengths change — the transition probability matrices change but the tree structure doesn't. For `rate_log_sd` moves, only the ACRV rate categories change. For `kPrime` moves, only one partition sub-group changes. Currently `MkpLogLikelihood()` recomputes everything from scratch. Implement move-type-specific partial recomputation: cache partition likelihoods and only recompute the affected partition(s). |
| M-065 | P3 | DONE (C) | **Avoid re-extracting parent/child/edgeLength from tree every likelihood call.** `MkpLogLikelihood()` extracts `tree$edge[,1]`, `tree$edge[,2]`, `tree$edge.length` as R vectors, then passes them to C++. In the C++ inner loop (M-061+), these would be held directly in the C++ state struct, eliminating this per-call overhead. This task is partly subsumed by M-063 but is listed separately for tracking. |

---

## Progress display & UX

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-076 | P2 | DONE (C) | **`MkpWatchLog()` — live trace plot from disk log files.** Polling watcher that reads one or more log files (from `logFile` in `MkPrimeMCMC`) and redraws trace plots each cycle. Ctrl+C stops cleanly and returns last data invisibly. `MkLogPaths()` helper expands base log name to per-run paths. Auto-selects key scalar params (excludes br_, log_likelihood). |
| M-075 | P2 | DONE (C) | **Streaming log output (`logFile` / `bufferSize`).** Implement Tracer-compatible TSV output during MCMC. `logFile` writes scalar parameter samples to disk in batches (`bufferSize`), flushed at each buffer boundary. `ReadMkLog()` reads log files back into R. `nRuns > 1` creates per-run log files (`run_1.log`, `run_2.log`, …). Checkpoint includes streaming state for seamless resume. |
| M-074 | P2 | DONE (C) | **Per-parameter progress table during MCMC.** At each `checkEvery` interval (post-warmup), print a convergence table showing ESS and PSRF (when `nRuns >= 2`) for each key parameter — mirroring the format of `print.MkpDiagnostics`. Extend `.CheckConvergence()` to return full per-parameter `ess` and `psrf` vectors (and work for `nRuns = 1`, ESS-only). Add `.PrintProgressTable()` helper in `Convergence.R`. Update `doCheck` blocks in both `RunMkPrime()` and `ResumeMkPrime()`. |
| M-071 | P1 | DONE (C) | **Progress bar redesign.** Replace `{cli::pb_bar} current/total \| accept% ` with `iter \| ESS: N \| logP: N \| PSRF: N`. A bar is misleading because the run works toward a *convergence* condition, not a fixed iteration count. Make `nIter` default to `Inf` in `MkPrimeMCMC()` (rely on `minEss`/`maxPsrf`/`maxTime` stopping criteria). Update `ResumeMkPrime` progress display to match. Update vignette MCMC configs and docs. |
| M-072 | P2 | DONE (C) | **hyoliths.qmd tree-summary-demo chunk fixes.** (a) Add `par(mar = rep(0, 4))` before `plot(consensus, ...)` — do this via a chunk `fig.par` option or an explicit `par()` call. (b) The plot title says "rogues excluded" even when no rogue detection has run; make the title conditional on whether rogues were actually identified and removed. |

## Convergence diagnostics

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-073 | P1 | DONE (C) | **Tree ESS in ConvergenceDiagnostics.** Add `trees = TRUE` argument. If `treess` and `TreeDist` installed and `!isFALSE(trees)`, compute topology ESS using `treess::treess(perRunTrees, TreeDist::RobinsonFoulds, methods = treess::getESSMethods(TRUE))`. Subsample each run to ≤1000 trees; if `interactive()` and total trees > threshold, emit a progress message. Extract `frechetCorrelationESS` and `medianPseudoESS`, sum across runs. `print.MkpDiagnostics` always shows a topology row (NA when not computed). Add `treess` and `TreeDist` to Suggests. |

---

## Posterior & analysis improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-066 | P1 | DONE (C) | **Automatic burnin selection for MkPosterior.** The posterior object should support tunable burnin that maximizes ESS while minimizing PSRF. Add a `burnin` parameter (number of post-warmup samples to discard) with an automatic selection method that optimizes the ESS/PSRF trade-off. The `MkPosterior` accessors (samples, trees, summary, etc.) should default to returning only post-burnin samples from all runs combined. Provide a method to adjust burnin after the fact (e.g. `SetBurnin(posterior, n)` or a `burnin` argument to accessor functions). Currently, only warmup is discarded; additional burnin is common practice when chains take time to find the typical set even after warmup ends. |
| M-067 | P2 | DONE (C) | **Vignette and package documentation references via inst/REFERENCES.bib + Rdpack.** Create `inst/REFERENCES.bib` with all cited references (Sun et al. 2018, Lewis 2001, Xie et al. 2011, etc.). Set up Rdpack for `\insertRef{}` in roxygen docs. Update `hyoliths.qmd` to use `bibliography: ../inst/REFERENCES.bib` (or copy to vignettes/). Check `../TreeTools/` for a working template of this setup. Add `Rdpack` to Imports in DESCRIPTION. |
| M-068 | P2 | DONE (C) | **Rogue taxon suppression in hyoliths.qmd consensus tree.** In the tree summary section of `vignettes/hyoliths.qmd`, use `Rogue::QuickRogue()` to identify rogue taxa, then exclude them from the consensus tree. Display which taxa were identified as rogues and show the cleaned consensus. Add `Rogue` to Suggests in DESCRIPTION. |

---

## Phase 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing insufficient). |

---

## Phase 6b: TreeSearch integration (deferred)

*Not yet broken into tasks.*

Planned work:
- TreeSearch GUI integration hook ("Bayesian (Mk')" mode in EasyTrees)
