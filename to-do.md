# MkPrime Task Queue

## How this works

- Tasks are sorted by priority (highest first within each status group).
- An agent claims a task by changing its status to `ASSIGNED (X)`.
- On completion, **delete** the row from this file and append a summary row
  to `completed-tasks.md`.

Task IDs use `M-nnn` prefix (MkPrime) to avoid collision with TreeSearch
`T-nnn` IDs.

---

## Phase 4: Tree search — COMPLETE

All Phase 4 tasks (M-022 through M-028) are done. See `completed-tasks.md`.

## Phase 5: Parallel tempering + convergence

*Not yet broken into tasks.*

Planned work:
- Temperature ladder (geometric spacing)
- Adaptive temperature tuning
- Multiple independent runs
- ESS monitoring (coda)
- PSRF monitoring (Gelman-Rubin)
- Stopping rules (ESS, PSRF, max iter, max time)
- Checkpointing
- C++ hot loop migration (move MCMC inner loop from R to C++)
