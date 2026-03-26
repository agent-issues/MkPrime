# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

---

## Phase 5: Parallel tempering + convergence — COMPLETE

All Phase 5 tasks (M-029 through M-036) are done. See `completed-tasks.md`.

## Phase 6: Progress display + GUI hooks — COMPLETE (core)

M-037 through M-039 done. Core progress display infrastructure in place:
- Callback architecture (plot_every, progress_fn)
- Console trace plots (mkp_trace_plot)
- PNG output for Shiny polling (mkp_png_progress)

## Phase 6b: TreeSearch integration (deferred)

*Not yet broken into tasks.*

Planned work:
- TreeSearch GUI integration hook ("Bayesian (Mk')" mode in EasyTrees)

## Phase 7: Extensions

### 7a: `coding = "informative"` ascertainment correction

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-044 | P1 | OPEN | C++ `singleton_site_prob_jc()` — compute P(autapomorphy) by running n singleton pseudo-characters through pruning. With ACRV averaging. |
| M-045 | P1 | OPEN | C++ `singleton_site_prob_mkn()` — asymmetric version: 2n pseudo-characters (each tip as sole 0 among 1s, and sole 1 among 0s). |
| M-046 | P1 | OPEN | R-level `coding = "informative"` support: `MkPrimeModel()` accepts `"informative"`; `MkpLogLikelihood()` computes P(uninformative) = P(constant) + P(singleton) and applies correction. |
| M-047 | P1 | OPEN | Tests: validate `coding = "informative"` against RevBayes reference likelihoods on hyoliths dataset. Unit tests for singleton prob functions. |

### 7b: Partition-specific rate scalars

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-048 | P2 | OPEN | Add `rate_neo` parameter (rate multiplier for neomorphic partition, default 1.0). Multiply branch lengths by rate scalar in likelihood. LogNormal(0, 2) prior. Scale proposal. Transformational rate fixed at 1.0 for identifiability. |
| M-049 | P2 | OPEN | Tests: partition rate scalar recovery on simulated data; verify identifiability constraint. |

### 7c: Stepping-stone marginal likelihood

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-050 | P2 | OPEN | Implement `mkp_stepping_stone()`: run MCMC at K power posteriors (β from 0→1), compute marginal likelihood via stepping-stone estimator (Xie et al. 2011). Return log-marginal-likelihood with SE. |
| M-051 | P2 | OPEN | Tests: verify marginal likelihood on simple tree where analytical solution is tractable. |

### 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing is insufficient). |
