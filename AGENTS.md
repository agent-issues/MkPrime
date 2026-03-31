> **Before starting work in this directory, read [`../AGENTS.md`](../AGENTS.md)**
> for multi-agent coordination rules, build/test infrastructure, GHA workflows,
> and worktree discipline. That file is the authoritative reference for all
> cross-package agent operations.

# MkPrime — Agent Development Notes

## Current phase: post-core, optimization & GUI integration

All core phases (1–10) complete. 155+ tasks delivered. Current work:

- **Phase 6b** (TreeSearch GUI integration): M-080 assigned to Agent C.
- **M-155** (P1): Gibbs kPrime sweep batch-by-k' optimization — **DONE**
  (`224bda8`). Batched partition-level pruning + progressive early
  termination. 3.6× speedup on full Sun2018 (see "Performance notes").
- **M-131** (P2): Warmup stabilisation validation study (Hamilton HPC).
  BLOCKED on M-155.
- Standing tasks (S-RED, S-PROF, S-COORD) at P1.

Agents should:

- Monitor `to-do.md` for task selection.
- Check `coordination.md` for strategic context.
- Keep code clean and well-tested.

> **Worktree note:** If you are working in a feature worktree (e.g.
> `mkp-parallel`, `mkp-gibbs`), **always read and write coordination
> files from `../mkp/`**, not from your own worktree directory.
> Each feature branch contains a stale copy of these files from when
> the branch was cut. The `../mkp/` directory (the `main` worktree)
> is the single source of truth for:
> - `to-do.md` (task queue)
> - `coordination.md` (phase tracking)
> - `completed-tasks.md` (archive)
> - `u.nnn` issue files
> - `remote-jobs.md` (pending async jobs)
> - `.positai/plans/` (plan files)
>
> Feature-branch code changes go in your own worktree as usual.
> Coordination file changes go in `../mkp/` and are committed to `main`.

---

## Package overview

**MkPrime** implements Bayesian phylogenetic inference under the Mk'
(Mk-prime) model for discrete morphological characters. It supports
three character models in a single analysis:

| Type | Model | Key feature |
|------|-------|-------------|
| Neomorphic | MkN (asymmetric binary) | Shared `rate_loss` across characters |
| Transformational | Mk' (infer true k per character) | Per-character k'_i from shared hyperprior |
| Known state space | Mk(k) | User-specified k, no k' inference |

The architecture follows StratoBayes: C++ hot loop (Rcpp) for likelihood
evaluation and MCMC inner loop, R orchestration for setup, adaptation,
tempering, convergence monitoring, and I/O.

### Key design decisions

1. **Unrooted trees.** Stationary root frequencies make all models
   (including asymmetric MkN) time-reversible; root placement irrelevant.
2. **Model k'_i not u_i.** True states (k') avoids model-graph cycles
   where u depends on kObs.
3. **No kMax.** JC eigendecomposition is O(1) analytical; prior +
   relabelling correction naturally constrain k'.
4. **Per-character k'_i from shared hyperprior.** Pools information
   across characters.
5. **`coding = "variable"` first.** `"informative"` deferred to Phase 7.

Full design rationale: `../.positai/plans/mkprime-r.md`.

---

## Dependencies

| Type | Packages |
|------|----------|
| **Imports** | Rcpp, ape, TreeTools, cli |
| **Suggests** | treess (tree-space ESS), TreeSearch (dataset + GUI), testthat |
| **LinkingTo** | Rcpp |

Note: `coda` dependency was dropped — native rank-normalized R-hat
(Vehtari et al. 2021) and FFT-based ESS replace all `coda` calls.

---

## File layout

```
mkp/
├── DESCRIPTION
├── NAMESPACE
├── R/
│   ├── MkPrimeData.R        # phyDat input, character classification, partitioning
│   ├── MkPrimeModel.R       # Model specification: priors, ACRV options
│   ├── MkPrimeMCMC.R        # MCMC configuration object (nIter, thin, tuning)
│   ├── RunMkPrime.R         # Main entry: MCMC engine, move dispatch, adaptation
│   ├── MkPosterior.R        # Results object (samples, print/summary/plot)
│   ├── Convergence.R        # ESS, R-hat monitoring, stopping rules
│   ├── proposals.R          # Scale, BetaSimplex, BoundedIntegerWalk proposals
│   ├── SteppingStone.R      # Stepping-stone marginal likelihood estimation
│   ├── BayesianModule.R     # Shiny module for TreeSearch GUI integration
│   ├── EasyMkPrime.R        # Simplified entry point for GUI
│   ├── PlotDuringMCMC.R     # Live progress plots during MCMC
│   ├── streaming.R          # Streaming log output (Tracer-compatible TSV)
│   ├── burnin.R             # Burnin detection and removal
│   ├── ess.R                # Native ESS computation (FFT-based)
│   ├── treeESS.R            # Tree-topology ESS via pseudo-ESS
│   ├── likelihood.R         # R-side likelihood wrappers
│   ├── partition.R          # Partition management
│   ├── acrv.R               # R-side ACRV helpers
│   └── RcppExports.R        # Auto-generated
├── src/
│   ├── mcmc.cpp             # Hot loop: propose → evaluate → accept/reject, batch runner
│   ├── mcmc_likelihood.cpp  # Felsenstein pruning, per-partition likelihood, Gibbs helpers
│   ├── mcmc_state.h         # McmcData / McmcState structs, function declarations
│   ├── node_cl_cache.h      # Partial conditional likelihood cache (dirty-flag invalidation)
│   ├── gibbs_partial_cl.h   # Gibbs topology move CL helpers
│   ├── likelihood.cpp       # R-callable likelihood functions
│   ├── rate_matrix.cpp      # JC(k') construction, analytical P(t) = exp(Qt)
│   ├── tree_moves.cpp       # SPR, NNI, TBR, pSPR proposals in C++
│   ├── proposals.cpp        # Scalar/branch proposals (scale, Dirichlet, etc.)
│   ├── acrv.cpp             # Among-character rate variation (lognormal)
│   ├── ascertainment.cpp    # Constant-site / singleton-site corrections
│   ├── corrections.cpp      # Relabelling correction for Mk'
│   ├── tree_ess.cpp         # Tree-distance ESS computation
│   ├── init.cpp             # Package init (registration)
│   ├── fitch.h              # Fitch parsimony (for starting tree score)
│   └── RcppExports.cpp      # Auto-generated
├── inst/
│   └── REFERENCES.bib
├── tests/testthat/          # ~4000 tests
├── man/
└── vignettes/
    └── hyoliths.qmd         # Full worked example (Sun2018, 54 taxa)
```

---

## Stale process cleanup — REQUIRED PROTOCOL

On Windows, when a bash tool timeout fires or a conversation is
interrupted, child `Rscript.exe` processes are **not killed** — they
become orphans consuming CPU. These accumulate silently.

### When to kill stale processes

**At conversation start**, before doing any work:

```bash
taskkill //F //IM Rscript.exe 2>/dev/null
```

**Before every subprocess launch** (already in subprocess template below):

```bash
taskkill //F //IM Rscript.exe 2>/dev/null; sleep 1
```

This is safe because agent work never depends on a previously-launched
Rscript surviving across tool calls. If the user has their own Rscript
processes running, the `2>/dev/null` suppresses "not found" errors, and
the user can restart them. In practice the user's R work runs inside
the RStudio session (rsession.exe), not via Rscript.

### S-COORD check

S-COORD rounds should check for orphan Rscript processes as a standard
step:

```bash
tasklist //FI "IMAGENAME eq Rscript.exe" 2>/dev/null
```

If any are found, kill them and note it in the S-COORD log.

---

## Build workflows

### Single-agent phase (current)

Standard devtools workflows in a subprocess:

```bash
# Compile + load (always in subprocess, never RStudio session)
cd mkp
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all()"

# Run targeted tests
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all(); testthat::test_file('tests/testthat/test-foo.R')"

# Regenerate docs (use load_installed to avoid debug .o contamination)
Rscript -e "roxygen2::roxygenise(load_code = roxygen2::load_installed)"
```

### Multi-agent phase

When multiple agents are active, switch to the renamed-package build
system (`build-agent.sh` / `test-agent.sh`) following the parent
`AGENTS.md` protocol.

### Validation

Use **GitHub Actions** for full test suites and R CMD check. Local builds
are for targeted iteration only (build + run 1–2 specific test files).

---

## Running `RunMkPrime()` in subprocesses — REQUIRED PROTOCOL

`RunMkPrime()` with convergence criteria can run for hundreds of thousands
of iterations. On Windows, when the bash tool timeout fires, the child
Rscript process is **not killed** — it becomes an orphan consuming CPU.
Repeated retries accumulate orphan processes that must be killed manually
via Task Manager.

**Triple-guard every subprocess `RunMkPrime()` call:**

1. **`maxTime`** — always set in the `RunMkPrime()` call (e.g. `maxTime = 120`).
   This is the application-level stop; the MCMC will finish the current
   batch and exit cleanly.
2. **`setTimeLimit()`** — wrap the R code in `setTimeLimit(elapsed = 150)` as a
   backup in case `maxTime` isn't checked frequently enough.
3. **Tool timeout** — set the bash tool `timeout` parameter to ~180s (above
   the R-level limits so they fire first).

**Template:**

```bash
taskkill //F //IM Rscript.exe 2>/dev/null; sleep 1
Rscript -e '
  setTimeLimit(elapsed = 150)
  pkgbuild::compile_dll(debug = FALSE); devtools::load_all()
  # ... setup ...
  posterior <- RunMkPrime(..., maxTime = 120, ...)
  # ... diagnostics ...
' 2>&1
echo "EXIT: $?"
```

**Never** run `RunMkPrime()` with convergence criteria (`minEss`, `minTreeEss`,
`maxRhat`) without also setting `maxTime`. The `nIter` cap alone is not
sufficient — warmup can consume most of the iteration budget.

---

## Build failure recovery

| Symptom | Cause | Fix |
|---------|-------|-----|
| Debug `.o` contamination | `roxygen2::roxygenise()` default uses `debug=TRUE` | `rm -f src/*.o src/*.dll`, rebuild |
| "Access is denied" | Another R process has DLL loaded | Kill the process or wait |
| Namespace errors after rename | `RcppExports.cpp` stale | `Rscript -e "Rcpp::compileAttributes()"` |
| Orphan Rscript processes | Subprocess timeout didn't kill child | `taskkill //F //IM Rscript.exe` |

---

## VTune profiling

VTune is installed at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/`.
This PC is Intel i7-10700 (10th gen); hardware sampling works.

**When profiling, use the `r-package-profiling` skill** to locate VTune
and follow the full workflow (build with symbols, driver script,
collection, report). Key points:

- Override `DLLFLAGS` via `MAKEFLAGS` env var (not `src/Makevars.win`)
- Add `-g -fno-omit-frame-pointer` to `PKG_CXXFLAGS` in `src/Makevars.win`
- Remove profiling flags after collection

---

## Performance notes

### Gibbs kPrime sweep (S-PROF round 3, 2026-03-31)

**Pre-optimization** (per-character individual traversals):

| Move | Cost | Coverage |
|------|------|----------|
| Gibbs kPrime sweep (moveType 25) | 340 ms | All 99 trans chars |
| 99 × int_walk kPrime (moveType 7) | 33 ms | All 99 trans chars |
| Block kPrime shift (moveType 26) | 0.66 ms | All 99 trans chars (collective mode) |

Root cause of 10× overhead: `single_char_loglik_jc()` did per-character,
per-candidate-k' individual tree traversals with per-call heap allocation.

**Post-optimization (M-155, `224bda8`):** Batched partition-level pruning
with progressive early termination.

| Config | Chars | nCat | Old (est.) | New | Speedup |
|--------|-------|------|------------|-----|---------|
| all-trans | 225 | 1 | ~70 ms | 35.7 ms | ~2× |
| all-trans | 225 | 6 | ~770 ms | 214 ms | ~3.6× |
| with-neo | 31 | 1 | ~70 ms | 13.0 ms | ~5.4× |
| with-neo | 31 | 6 | ~107 ms | 71 ms | ~1.5× |

nCat scaling remains roughly linear. Speedup is larger with more
characters (better amortization of batched traversals).

### VTune hotspot profile (S-PROF round 4, 2026-03-31)

Sun2018 (54 taxa, 225 chars, all transformational, nCat=6 ACRV).
15k iterations, ~90s CPU time. VTune 2025.10, user-mode sampling.

| Function | Source | CPU Time | % |
|----------|--------|----------|---|
| `pruning_jc_acrv_persite` | mcmc_likelihood.cpp | 53.0s | 59.1% |
| `_expl_internal` (exp) | compiler runtime | 10.5s | 11.7% |
| `constant_site_prob_jc` | ascertainment.cpp | 5.9s | 6.6% |
| `pruning_jc_acrv_flat` | mcmc_likelihood.cpp | 3.5s | 3.9% |
| `jc_transition` | gibbs_partial_cl.h | 1.3s | 1.4% |
| `std::vector` copies | stl_vector.h | ~1.3s | ~1.5% |
| Rcpp bounds checks | traits.h | ~1.0s | ~1.1% |
| Everything else | | ~13s | ~14.5% |

Key findings:
1. **Gibbs kPrime sweep** (`persite` + ascertainment) dominates at ~65-70%
   of MkPrime CPU. Regular MH proposals (`flat`) are only 3.9%.
2. **`exp()` calls at 11.7%** — one per edge × rate category per traversal
   (line 370). Already amortized across characters; cannot factor across
   candidate k' values (eigenvalue depends on k). Fast approx exp could
   save ~5-8% of total.
3. **Ascertainment at 6.6%** — `constant_site_prob_jc` does separate tree
   traversals. Could batch like the main pruning.
4. **Vector copies at 1.5%** — heap allocation in hot path.
5. **Rcpp bounds checks at 1.1%** — `check_index` calls on operator[].

### Overall MCMC bottleneck (S-PROF round 2, 2026-03-28)

C++ Felsenstein pruning is ~90% of wall time. OPP-1–6 optimizations achieved
1.80× cumulative speedup. Diminishing returns on further pruning optimization.

---

## Known low-priority issues

### ETA estimation with tree-ESS-only convergence

When only `minTreeEss` is set (no `minEss`), the ETA estimation in the
progress ticker doesn't work. The ETA code defaults `etaTarget` to
`mcmc$minEss` which is NULL, causing `.EstimateEta()` to return NULL.
Low priority — unusual configuration.

### Tree ESS not enforced in log-based convergence

`.CheckConvergenceFromLogs()` (used by parallel runs and serial Phase 2
cross-run R-hat loop) cannot enforce `minTreeEss` because trees are not
stored in log files. Tree ESS is enforced per-run by `.CheckConvergence()`
in Phase 1. This is by design but creates an asymmetry.

---

## Validation datasets

Small annotated datasets in `../neotrans/inst/matrices/` with Excel
metadata files that classify characters as Neomorphic/Transformational:

| Project | Taxa | Chars | Neo | Trans | Taxon |
|---------|------|-------|-----|-------|-------|
| 3832 | 10 | 27 | 12 | 15 | Canthyloscledidae |
| 950 | 12 | 9 | ? | ? | Hexacorallia |
| 4789 | 13 | 12 | ? | ? | Agelacrinitinae |
| 1271 | 25 | 33 | ? | ? | Amaltheidae |
| 3408 | 19 | 30 | ? | ? | Galericini |

Excel files: `Project{N}_{author}.xlsx`, column `"Character Pattern"`
contains `"Neomorphic"` / `"Transformational"`.

Primary benchmark: **Sun2018** (54 taxa, 225 chars) via TreeSearch package.

---

## Reference material

| Resource | Location | Content |
|----------|----------|---------|
| Design plan | `../.positai/plans/mkprime-r.md` | Full model spec, architecture, phases |
| Mk' model math | `../revbayes-ms/.positai/expertise/mkprime-model.md` | Mathematical specification |
| RevBayes scripts | `../mkprime/` | Reference Mk' inference (RevBayes) |
| RevBayes C++ impl | `../revbayes-ms/` | Native `dnMkPrime` distribution |
| Best-practice phylo | `../neotrans/` | Compound Dirichlet, ACRV, MkN, etc. |
| StratoBayes arch | `../StratoBayes/` | Architectural template (MCMC engine, tempering) |

---

## Naming conventions

| Scope | Convention | Examples |
|-------|------------|---------|
| Exported functions | **PascalCase** | `RunMkPrime`, `MkPrimeData`, `ConvergenceDiagnostics` |
| Internal functions (dot-prefixed) | **`.PascalCase`** | `.InitState`, `.BuildMoves`, `.FinalizeModel` |
| Non-exported helpers (no dot) | **PascalCase** | `ProposeScale`, `LogPrior`, `DiscreteLognormalRates` |
| Function parameters | **camelCase** | `knownStates`, `fixTopology`, `treeLengthShape`, `checkEvery` |
| Local variables | **camelCase** | `nEdge`, `transIdx`, `charMatrix`, `startTree` |
| Model parameter names (MCMC state, column names) | **snake_case** | `tree_length`, `rate_loss`, `rate_log_sd`, `log_posterior` |
| S3 class names | **PascalCase** | `MkPrimeData`, `MkPrimeModel`, `MkPosterior` |

**Model parameter names** (`rate_loss`, `tree_length`, `rate_log_sd`, `kPrime`, `rel_br_lengths`, etc.)
are intentionally kept in snake_case because they are domain terminology that appears
in output column names, documentation, and RevBayes cross-references.

See `CONTRIBUTING.md` for the full rationale and worked examples.

---

## Key coordination files

All coordination files live in the **`main` worktree (`../mkp/`)**.
Feature worktrees contain stale copies — do not use them.

| File | Worktree to use | Purpose |
|------|-----------------|---------|
| `u.nnn` | `../mkp/u.nnn` | User issue files (agents triage → `to-do.md`, then delete) |
| `to-do.md` | `../mkp/to-do.md` | Task queue (active/open tasks only) |
| `remote-jobs.md` | `../mkp/remote-jobs.md` | Pending async jobs (Hamilton SLURM, long GHA) — check at `/assign` |
| `completed-tasks.md` | `../mkp/completed-tasks.md` | Archive of completed tasks |
| `coordination.md` | `../mkp/coordination.md` | Strategic plan and phase tracking |
| `agent-<letter>.md` | `./agent-<letter>.md` | Agent progress log (local-only, gitignored — one per worktree) |
| `AGENTS.md` | `../mkp/AGENTS.md` | This file — package conventions and architecture |
| `.positai/plans/` | `../mkp/.positai/plans/` | Plan files |
