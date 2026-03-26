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

*Not yet broken into tasks. Phase 3 just completed.*

Planned work:
- SPR and NNI proposals (C++)
- Compound Dirichlet branch length prior
- Uniform topology prior
- exp_steps default from Goloboff-Wagner parsimony score
- Tree logging (Newick to file)
- Validate topology moves preserve tree invariants
