# dev/profiling

Working directory for the `/profile` skill (one hot-path focus area per
invocation). The skill definition holds the doctrine; this README records what
lives here.

**The live work list is
[GitHub Issues on `agent-issues/MkPrime`](https://github.com/agent-issues/MkPrime/issues)**
labelled `profiling`. **The issue number is the finding id** — do not mint new
`T-nnn` ids. Legacy `T-nnn` ids stay resolvable as search keys in
`findings-archive.md` and `../red-team/migration-map.tsv`; never filter a
duplicate search on that prefix, because new issues carry none.

## Round records stay here, not in Discussions

A deliberate divergence from `/red-team`. Profiling has one serial rotation with
no per-area record cadence to key off, `log.md` doubles as the context
`baselines.md` regressions are read against, and Discussion categories are capped
at 25 per repo. Round records go in `log.md`.

## Layout

```
dev/profiling/
  README.md             # this file
  focus-areas.md        # ranked hot paths; statuses NEW / PROFILED / OPTIMISED / AT-LIMIT / SKIPPED
  log.md                # per-round records + last_focus
  baselines.md          # current timings, refreshed each round; the regression reference
  findings-archive.md   # FROZEN 2026-09-18 — anti-duplication memory only
  drivers/              # one runnable driver per focus area
```

## The verification bar

File a finding only after an **isolated micro-benchmark reproduces the predicted
delta**. A drop in the VTune hotspot share is not a verified win. An unverified
lead may be filed so it is tracked, but its body must say plainly that it is a
candidate and what measurement is owed first — see #18 for the shape.

**No file-then-close.** A finding fixed in the round that found it gets no issue:
record it in `log.md` and the fixing PR, and stop.

## Platform

Windows/VTune and Linux/gcc are different binaries (MinGW ucrt vs glibc) and both
are production targets. Never quote a whole-program percentage from one platform
alone: cross-check, or label the finding platform-specific.

`result_*/` directories and `.vtune-lib/` are gitignored and **deleted at the end
of each round** — they are hundreds of MB, and reusing them across code changes
produces misleading hotspots.

MkPrime specifics (profiling flags in `src/Makevars.win`, which must never be
committed; VTune paths; hot-path data) are in `../../.AGENTS/memory/performance.md`
and `../../AGENTS.md`.
