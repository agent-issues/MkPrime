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

## Phase 6: Progress display + GUI hooks

*Not yet broken into tasks.*

Planned work:
- PlotDuringMCMC: log-posterior traces, parameter traces
- PNG-based progress for Shiny integration
- TreeSearch GUI module (Suggests: MkPrime)

## Phase 7: Extensions (future)

*Not yet broken into tasks.*

Planned work:
- `coding = "informative"` ascertainment correction
- Partition-specific rate scalars
- Beta-distributed Q-matrix heterogeneity (siteMatrices)
- Stepping-stone / path sampling for model comparison
- TBR moves (if needed for mixing)
- HMC for branch lengths (if needed for mixing)
