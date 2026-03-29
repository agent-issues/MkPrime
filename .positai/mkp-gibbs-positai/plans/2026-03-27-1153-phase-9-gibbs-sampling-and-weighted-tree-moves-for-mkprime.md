# Sun 2018 stress test + red-team review + optimization roadmap

## Overview

Three activities:
1. **Sun 2018 mixing comparison** — repeat M-091 on the Sun 2018 hyolith dataset (54 taxa, 225 chars) to see whether Gibbs/Weighted moves show ESS/s gains on a larger tree.
2. **Red-team correctness review** — systematic audit of Phase 9 C++ move implementations for correctness bugs.
3. **Optimization roadmap** — profile-informed priorities for the next round of performance work.

---

## Part 1: Sun 2018 stress test

### Setup
- Load Sun 2018 from TreeSearch: `TreeSearch::ReadCharacters(system.file("datasets/Sun2018.nex", package = "TreeSearch"))`
- Build NJ starting tree, ensure tip labels match
- Use `MkPrimeModel()` defaults

### Configurations (same structure as M-091)

| Config | Moves enabled | Key flag |
|--------|--------------|----------|
| (a) Baseline | NNI + SPR + standard | defaults |
| (b) +Gibbs | + GibbsSPR + GibbsSubtreeSwap | `gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE` |
| (c) +Weighted | + WeightedSPR | `weightedSpr = TRUE` |
| (d) +BlockGibbs | (b) + blockGibbsBranch | `blockGibbsBranch = TRUE` |

### Parameters
- `nIter = 20000`, `warmup = 10000`, `thin = 5` → 2000 post-warmup samples
- Same seed across configs
- Single chain, single run (no tempering), debug build
- Record wall time, ESS(log_posterior), ESS(tree_length), acceptance rates per move

### Expected outcome
On 54 taxa, standard SPR acceptance should be much lower (~0.1–0.5%), making Gibbs SPR's higher acceptance rate more valuable. The critical question: does the ESS/s crossover point favour Gibbs moves on this larger tree?

---

## Part 2: Red-team correctness review

Systematic code audit of all Phase 9 move implementations. Focus areas:

### 2a. Gibbs SPR (lines 324–475)

**Candidate set completeness:** Current state included via `wOrig` — proper discrete Gibbs. Need to verify no valid SPR rearrangements are silently dropped by the descendant/incidence filter.

**Candidate set symmetry:** For pure discrete Gibbs with acceptance = 1, the candidate set must be the same size in both directions (forward prune from current position, reverse prune from proposed position). For unrooted binary trees, the number of valid SPR positions after pruning is constant for a given n. Verify.

**Branch-length handling:** tau = 0.5 midpoint insertion means the regraft edge is split equally. The reverse move re-merges via lMerge = absLen[parentRow] + absLen[sibRow]. This is the same total length. No Jacobian needed.

### 2b. Weighted Branch Scale (lines 591–699)

**Hastings ratio:** `logHastings = log(w_oldBin) + logBeta(f_old|...) - log(w_chosenBin) - logBeta(f_new|...)`. The bin weights do NOT cancel — they are the reverse/forward bin selection probabilities. This is the correct MH-within-Gibbs ratio. ✓

### 2c. Block Gibbs Branch Sweep (lines 703–861)

**Sweep semantics:** Each pair accepted/rejected independently with `currentLL` updated after each accepted pair. This is a valid composition of MH kernels (systematic scan Gibbs). ✓

**Global static `s_branchBins`:** Thread safety concern if OpenMP is ever added. Currently safe (separate processes via `future`). Flag for future.

### 2d. Weighted SPR (lines 864–1127)

**Hastings ratio:** Topology selection cancels because the same candidate set appears in both directions. Ratio reduces to branch-fraction component. Verify candidate set symmetry (same concern as 2a).

**`fOld` computation:** `fOld = absLen[parentRow] / lMerge` — the fraction assigned to the parent edge at the current regraft point. Must match what the reverse move would compute.

### 2e. Memory allocation in candidate loops

GibbsSPR and WeightedSPR `clone()` Rcpp vectors for each candidate (lines 407–409, 980–981). For 54 taxa with ~104 edges, this is ~104 clones per Gibbs SPR call. Profile whether this dominates.

---

## Part 3: Optimization roadmap

### 3a. Partial likelihood reuse in Gibbs moves (high impact)
When pruning subtree v and evaluating regraft at different positions, CLs for nodes not on the path between old/new positions are unchanged. A "lazy pruning" approach could reduce per-candidate cost from O(N × C) to O(D × C) where D is path length.

### 3b. Reduce clone() overhead
Replace per-candidate `clone()` with a single working buffer that gets mutated and restored. For GibbsSPR this would eliminate ~N IntegerVector + NumericVector allocations per call.

### 3c. Vectorize bin evaluations in weighted moves
For B bin midpoints differing only in two edge lengths, compute CL "prefix" once and only vary the "suffix" per bin.

### 3d. Production build comparison
Re-run at -O2 before investing in algorithmic work. Relative costs may shift significantly.

### Priority order
1. Production build comparison (low effort, high information)
2. Profile the hot path with bench/VTune (medium effort, informs all decisions)
3. Reduce clone() overhead (low–medium effort, moderate payoff)
4. Partial likelihood reuse in Gibbs moves (high effort, large payoff on big trees)
5. Bin evaluation vectorization (medium effort, moderate payoff)

---

## Implementation steps

1. Load Sun 2018 dataset, set up configs (a)–(d)
2. Run all four configs, collect wall time + ESS + acceptance rates
3. Present comparison table
4. Walk through red-team findings in prose
5. Summarize optimization priorities based on results + code review
