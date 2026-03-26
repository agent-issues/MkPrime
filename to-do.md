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

### 7a–c: COMPLETE

M-044 through M-051 done. See `completed-tasks.md`.
- `coding = "informative"` ascertainment correction (singleton probs + R-level support)
- Partition rate scalar `rate_neo` (neomorphic-specific rate multiplier)
- Stepping-stone marginal likelihood (`mkp_stepping_stone()`)

### 7d: Deferred extensions

| ID | Priority | Status | Description |
|----|----------|--------|-------------|
| M-052 | P3 | OPEN | Beta-distributed Q-matrix heterogeneity (siteMatrices). |
| M-053 | P3 | OPEN | TBR moves (if mixing diagnostics show SPR is insufficient). |
| M-054 | P3 | OPEN | HMC for branch lengths (if MH mixing is insufficient). |
