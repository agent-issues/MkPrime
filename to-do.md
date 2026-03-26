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

## Phase 4: Tree search

| ID | Pri | Description | Status | Notes |
|----|-----|-------------|--------|-------|
| M-024 | P1 | SPR proposal on unrooted binary trees: prune-regraft with correct Hastings ratio, branch length redistribution | ASSIGNED (A) | Depends on M-022 |
| M-025 | P1 | Tree move integration: add NNI/SPR to move schedule, adaptation tuning, postorder reordering | OPEN | Depends on M-023, M-024 |
| M-026 | P2 | exp_steps default from parsimony score: Goloboff-Wagner heuristic for Gamma rate | OPEN | Uses TreeSearch if available |
| M-027 | P2 | Tree logging to Newick file during MCMC: append-mode write, optional path in MkPrimeMCMC | OPEN | |
| M-028 | P1 | Integration validation: simulate data on known tree, recover topology from posterior | OPEN | Depends on M-025 |
