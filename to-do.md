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

## Infrastructure

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-041 | P3 | OPEN | Update parent AGENTS.md: replace `issues.md` block protocol with `u.nnn` file scanning. Agents scan for `u.*` files at session start, triage → `to-do.md`, delete file. Apply to MkPrime, TreeSearch, and parent `GitHub/`. |
| M-042 | P3 | DONE | Create `inst/_pkgdown.yml` for MkPrime (model on TreeSearch). Add pkgdown GHA workflow. |

## EasyMkPrime app improvements

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-043 | P2 | DONE | Auto-detect neomorphic characters: pre-populate the "Neomorphic characters" input with indices where state `0` occurs. Characters with `0` and states other than `1` should be elevated to transformational. |

## Phase 6b: TreeSearch integration (deferred)

*Not yet broken into tasks.*

Planned work:
- TreeSearch GUI integration hook ("Bayesian (Mk')" mode in EasyTrees)

## Phase 7: Extensions (future)

*Not yet broken into tasks.*

Planned work:
- `coding = "informative"` ascertainment correction
- Partition-specific rate scalars
- Beta-distributed Q-matrix heterogeneity (siteMatrices)
- Stepping-stone / path sampling for model comparison
- TBR moves (if needed for mixing)
- HMC for branch lengths (if needed for mixing)
