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

## Phase 2: Likelihood engine (C++) — nearly complete

| ID | Pri | Status | Blocks | Description | Notes |
|----|-----|--------|--------|-------------|-------|
| M-012 | P1 | ASSIGNED (A) | — | **Likelihood validation.** Compare MkPrime likelihoods against RevBayes on a set of reference trees and datasets. Test: (a) Mk with known k, (b) MkN, (c) Mk' with fixed k', (d) with ACRV, (e) with ascertainment correction. Must match to within numerical tolerance (~1e-6). | Use `../mkprime/` reference data. May need to generate RevBayes reference values. |
