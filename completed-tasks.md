# MkPrime — Completed Tasks

| ID | Description | Agent | Date | Notes |
|----|-------------|-------|------|-------|
| M-001 | Package skeleton (DESCRIPTION, Rcpp stub, testthat harness) | A | 2026-03-26 | `693c6ee` |
| M-002 | `MkPrimeData()` — phyDat input and character classification | A | 2026-03-26 | `5afd39e`. 19 tests. |
| M-003 | Character partitioning by type and kObs | A | 2026-03-26 | `efc1ba0`. 10 new tests. |
| M-004 | Data validation and edge cases | A | 2026-03-26 | `46386b6`. Invariant drop, index remapping. |
| M-005 | Comprehensive unit tests for data layer | A | 2026-03-26 | `46386b6`. 70 total tests. |
| M-006 | JC(k') rate matrix and analytical P(t) | A | 2026-03-26 | `da4492f`. O(1) eigendecomposition. |
| M-007 | MkN asymmetric binary rate matrix and P(t) | A | 2026-03-26 | `da4492f`. rate_loss parameterization. |
| M-008 | Felsenstein pruning (JC + MkN) | A | 2026-03-26 | `8e0fd20`. Post-order traversal, hand-verified. |
| M-009 | ACRV discretized lognormal rate categories | A | 2026-03-26 | `f2d5ee3`. 6-category midpoint quantiles. |
| M-010 | Ascertainment bias correction (variable coding) | A | 2026-03-26 | `ffb5769`. Constant site probability. |
| M-011 | Mk' relabelling correction | A | 2026-03-26 | `e001c74`. lgamma-based, scalar + batch. |
| M-012 | Likelihood validation + mkp_loglikelihood() | A | 2026-03-26 | `25e08ad`. Matches phangorn to machine epsilon. 194 total tests. |
| M-013 | Fix per-character k' in likelihood | A | 2026-03-26 | `41105eb`. kPrime > kObs now uses correct JC(k') matrix. |
| M-014 | Prior specification & log-prior (MkPrimeModel) | A | 2026-03-26 | `41105eb`. Gamma/LogNormal/Geometric priors with boundary checks. |
| M-015 | Proposal functions (Scale, BetaSimplex, BoundedIntegerWalk) | A | 2026-03-26 | `41105eb`. Correct Hastings ratios, reversibility tested. |
| M-016 | MCMC state object & initialization | A | 2026-03-26 | `41105eb`. tree_length + rel_br_lengths parameterization. |
| M-017 | Single-chain MH engine (RunMkPrime) | A | 2026-03-26 | `41105eb`. R-side loop, weighted move schedule. |
| M-018 | Sample accumulation & progress bar | A | 2026-03-26 | `41105eb`. In-memory matrix, cli progress. |
| M-019 | Proposal adaptation | A | 2026-03-26 | `41105eb`. Multiplicative tuning during warmup. |
| M-020 | MkPosterior result object | A | 2026-03-26 | `41105eb`. print/summary/plot methods, ESS via coda. |
| M-021 | Integration validation | A | 2026-03-26 | `41105eb`. Posterior recovery on simulated data. 466 total tests. |
