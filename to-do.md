# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.
- Tasks awaiting GHA results: `PARKED (<Letter>, GHA <run_id>)`.
- Tasks with an open PR awaiting human merge: `PR #N (<Letter>)`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

---

## Phase 1: Foundation (package skeleton + data layer)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| M-003 | P1 | ASSIGNED (A) | M-006 | **`MkPrimeData()` — phyDat input and character classification.** Accept a `phyDat` object + optional character annotations. Auto-detect neomorphic (binary with state 0 = "absent") vs transformational (multi-state). Allow user override. Return a structured list with partitioned characters. | Input: `phyDat` from TreeTools/ape. |
| M-003 | P1 | ASSIGNED (A) | M-006 | **Character partitioning by kObs.** Within each type (neomorphic/transformational/known), group characters by their observed state count (kObs). Build the internal data structures needed by the likelihood engine: tip-state matrices per partition. | This is the data layout the C++ code will consume. |
| M-004 | P1 | OPEN | — | **Data validation and edge cases.** Handle: invariant characters (warn + drop), missing data (NA / `?` states), characters with only 1 observed state, empty taxa. Test each case. | |
| M-005 | P1 | OPEN | — | **Unit tests for data layer.** Test character classification, partitioning, edge cases. Use small hand-crafted `phyDat` objects. Validate round-trip: `phyDat` → `MkPrimeData()` → internal structures. | |

## Phase 2: Likelihood engine (C++)

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| M-006 | P1 | OPEN | M-008, M-009, M-010 | **JC(k') rate matrix and P(t).** Implement analytical eigendecomposition for JC(k'): `P_ij(t) = (1/k') + ((k'-1)/k') * exp(-k't/(k'-1))` for i=j, `(1/k') - (1/k') * exp(-k't/(k'-1))` for i≠j. Rcpp-exported function that returns P(t) matrix given k' and t. | O(1) computation — no matrix exponentiation needed. |
| M-007 | P1 | OPEN | M-008, M-009 | **MkN rate matrix and P(t).** Implement asymmetric 2-state Q-matrix with `rate_loss` parameter. Analytical P(t) for 2×2. Rcpp-exported. | Simpler than JC(k') but different parameterization. |
| M-008 | P1 | OPEN | M-010, M-011 | **Felsenstein pruning (C++).** Post-order traversal computing conditional likelihoods at each node. Input: tree topology (parent-child arrays), branch lengths, tip-state data, P(t) matrices. Output: log-likelihood at root. | Core hot-path code. Design the data layout carefully (column-major tip states, flat arrays). |
| M-009 | P2 | OPEN | M-011 | **ACRV: discretized lognormal rate categories.** Implement 6-category discretized lognormal (Wagner 2012). Given `rate_log_sd`, compute category rates and weights. Integrate into likelihood: sum over rate categories. | `rate_log_sd ~ Gamma(1, 1)` prior. |
| M-010 | P1 | OPEN | M-011 | **Ascertainment bias correction.** Implement `coding = "variable"` correction: subtract log-probability of constant-site patterns from the likelihood. Requires computing likelihood of all-state-j patterns for each rate category. | See Lewis 2001, Allman et al. 2008. |
| M-011 | P1 | OPEN | — | **Mk' relabelling correction.** Per-character correction: `log C(k', kObs) = logFact[kObs] - logFact[k'] + logFact[k' - kObs] + kObs * log(k')`. Add to log-likelihood for each transformational character. | Simple arithmetic — but critical for correctness. |
| M-012 | P1 | OPEN | — | **Likelihood validation.** Compare MkPrime likelihoods against RevBayes on a set of reference trees and datasets. Test: (a) Mk with known k, (b) MkN, (c) Mk' with fixed k', (d) with ACRV, (e) with ascertainment correction. Must match to within numerical tolerance (~1e-6). | Use `../mkprime/` reference data. May need to generate RevBayes reference values. |
