# MkPrime — Agent Development Notes

You MUST read the `r-conventions` skill before writing any code.

## Dispatcher

Use the global `/dispatch` skill (`Skill(skill: "dispatch", ...)`). The skill
auto-detects per-repo overrides in `dev/dispatch/agent-brief.md` and
`dev/dispatch/ranker.txt` and falls back to the bundled defaults otherwise.
This repo does **not** ship its own `dispatch.sh` / `todo-lock.sh` — those
live in `~/.claude/skills/dispatch/`.

---

## Memory files (load on demand)

Technical reference material lives in `.AGENTS/memory/`. Load the relevant file
before starting work in that area:

| Memory file | Load when... |
|-------------|--------------|
| `architecture.md`         | Editing `src/*.cpp`/`.h`, adding Rcpp exports, reviewing R-level API or key design decisions |
| `performance.md`          | Running benchmarks, doing VTune profiling, interpreting hot-path data |
| `naming-conventions.md`   | Adding R files, writing roxygen, naming new functions or parameters |
| `validation-datasets.md`  | Picking a dataset for an MCMC test or validation study |
| `testing.md`              | Adding or modifying `tests/testthat/`, choosing test budgets |
| `subprocess-runmkprime.md`| Running `RunMkPrime()` from a bash subprocess (triple-guard for runaway) |

---

## Worktree note

If you are working in a feature worktree (e.g. `mkp-parallel`, `mkp-gibbs`),
**always read and write coordination files from `../mkp/`**, not from your
own worktree directory. Each feature branch contains a stale copy of these
files from when the branch was cut. The `../mkp/` directory (the `main`
worktree) is the single source of truth for:

- `to-do.md`, `coordination.md`, `completed-tasks.md`
- `u.nnn` issue files
- `remote-jobs.md`
- `dev/plans/`

Feature-branch code changes go in your own worktree as usual. Coordination
file changes go in `../mkp/` and are committed to `main`.

See `../AGENTS.md` → **Worktree discipline** for the full set of rules.

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

Architecture: C++ hot loop (Rcpp) for likelihood evaluation and the MCMC inner loop;
R for setup, adaptation, tempering, convergence monitoring, and I/O.

### Key design decisions

1. **Unrooted trees.** Stationary root frequencies make all models
   (including asymmetric MkN) time-reversible; root placement irrelevant.
2. **Model k'_i not u_i.** True states (k') avoids model-graph cycles
   where u depends on kObs.
3. **No kMax.** JC eigendecomposition is O(1) analytical; prior +
   relabelling correction naturally constrain k'.
4. **Per-character k'_i from shared hyperprior.** Pools information
   across characters.
5. **`coding = "variable"` first.** `"informative"` deferred.


---

## Dependencies

| Type | Packages |
|------|----------|
| **Imports** | Rcpp, ape, TreeTools, cli |
| **Suggests** | treess, TreeSearch (dataset + GUI), testthat |
| **LinkingTo** | Rcpp |

`coda` is **not** a dependency — native rank-normalized R-hat
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
│   ├── mcmc.cpp             # Hot loop: propose → evaluate → accept/reject
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
├── dev/
│   ├── plans/               # Plan files (active + archive/)
│   └── dispatch/            # Per-repo overrides for global /dispatch skill (agent-brief.md, ranker.txt)
├── .AGENTS/memory/          # Domain memory files (load on demand)
└── vignettes/
    └── hyoliths.qmd         # Full worked example (Sun2018, 54 taxa)
```

---

## Branch / worktree mapping (MkPrime)

See `../AGENTS.md` → **Current worktree mapping (MkPrime)** for the live
table. As of this writing:

| Directory | Branch |
|-----------|--------|
| `mkp`           | `main` |
| `mkp-gibbs`     | `feature/gibbs-weighted-moves` |
| `mkp-logseries` | `logseries-kprime-prior` |
| `mkp-parallel`  | `feature/parallel-runs` |
| `mkp-het`       | `feature/het-dirichlet-marginal` |
| `mkp-tbr`       | `feature/tbr-moves` |

---

## Build workflow (MkPrime-specific notes)

The parent `AGENTS.md` covers the canonical workflow (`build-agent.sh` /
`test-agent.sh`, renamed-DLL pattern, GHA-first validation). Two MkPrime
particulars:

1. **`RunMkPrime()` runs must be triple-guarded.** See
   `.AGENTS/memory/subprocess-runmkprime.md` for the exact template
   (`maxTime` + `setTimeLimit()` + bash `timeout`). Convergence-criterion
   runs (`minEss`, `minTreeEss`, `maxRhat`) **must** also set `maxTime`.

2. **Build failure recovery — mkp specifics:**

   | Symptom | Cause | Fix |
   |---------|-------|-----|
   | Namespace errors after rename | `RcppExports.cpp` stale | `Rscript -e "Rcpp::compileAttributes()"` |
   | Debug `.o` contamination | `roxygen2::roxygenise()` default uses `debug=TRUE` | `rm -f src/*.o src/*.dll`, rebuild |

For everything else (DLL locking, `roxygen` with `load_installed`, GHA dispatch,
parking, PR workflow), see parent `AGENTS.md`.

---

## VTune profiling

VTune is at `C:/Program Files (x86)/Intel/oneAPI/vtune/latest/bin64/`. This PC
is Intel i7-10700 (10th gen); hardware sampling works. Use the
`r-package-profiling` skill (`Skill r-package-profiling`) to locate VTune and
follow the full workflow. Key points specific to MkPrime:

- Override `DLLFLAGS` via `MAKEFLAGS` env var (not `src/Makevars.win`).
- Add `-g -fno-omit-frame-pointer` to `PKG_CXXFLAGS` in `src/Makevars.win`.
- **Remove profiling flags after collection.** `src/Makevars.win` must never
  be committed (see `.gitignore`).

Hot-path data lives in `.AGENTS/memory/performance.md`.

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
in Phase 1. By design — but an asymmetry.

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
| Active plans | `dev/plans/` | Per-task plan files (archive under `dev/plans/archive/`) |

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

**Model parameter names** (`rate_loss`, `tree_length`, `rate_log_sd`, `kPrime`,
`rel_br_lengths`, etc.) are intentionally kept in snake_case because they are
domain terminology that appears in output column names, documentation, and
RevBayes cross-references.

See `CONTRIBUTING.md` for the full rationale and worked examples.

---

## Key coordination files

All coordination files live in the **`main` worktree (`../mkp/`)**. Feature
worktrees contain stale copies — do not use them.

| File | Worktree to use | Purpose |
|------|-----------------|---------|
| `u.nnn`              | `../mkp/u.nnn`              | User issue files (agents triage → `to-do.md`, then delete) |
| `to-do.md`           | `../mkp/to-do.md`           | Task queue |
| `remote-jobs.md`     | `../mkp/remote-jobs.md`     | Pending async jobs |
| `completed-tasks.md` | `../mkp/completed-tasks.md` | Archive of completed tasks |
| `coordination.md`    | `../mkp/coordination.md`    | Strategic plan and phase tracking |
| `agent-<letter>.md`  | `./agent-<letter>.md`       | Agent progress log (local-only, gitignored) |
| `AGENTS.md`          | `../mkp/AGENTS.md`          | This file |
| `dev/plans/`         | `../mkp/dev/plans/`         | Plan files |
