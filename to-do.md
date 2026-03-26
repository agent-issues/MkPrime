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

## Phase 3: Basic MCMC (fixed topology)

*Not yet broken into tasks. Phase 2 just completed.*

Planned work:
- MH proposals for continuous parameters (Scale, BetaSimplex)
- BoundedIntegerWalk for k'_i
- Single chain, no tempering
- Sample accumulation to file
- Console progress bar
- MkPosterior result object
