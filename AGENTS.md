> **Before starting work in this directory, read [`../AGENTS.md`](../AGENTS.md)**
> for multi-agent coordination rules, build/test infrastructure, GHA workflows,
> and worktree discipline. That file is the authoritative reference for all
> cross-package agent operations.

# MkPrime — Agent Development Notes

## Current phase: greenfield development (Phase 1–3)

The package is being built from scratch. Agents should:

- Follow the phased plan in `coordination.md`. Phases are sequential;
  do not jump ahead unless the current phase is complete.
- Monitor `to-do.md` for task selection.
- Keep code clean and well-tested from the start — this is the foundation
  everything else builds on.
- Validate against known RevBayes likelihoods where possible (see
  `../mkprime/` for reference scripts and results).

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
| **Suggests** | coda (ESS/PSRF), treess (tree-space ESS), TreeSearch (dataset + GUI), testthat |
| **LinkingTo** | Rcpp |

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
│   ├── RunMkPrime.R         # Main entry point: single-chain MH engine
│   ├── proposals.R          # Scale, BetaSimplex, BoundedIntegerWalk proposals
│   ├── MkPosterior.R        # Results object (samples, print/summary/plot)
│   ├── Tempering.R          # (Phase 5) Parallel tempering
│   ├── Convergence.R        # (Phase 5) ESS, PSRF monitoring, stopping rules
│   ├── PlotDuringMCMC.R     # (Phase 6) Live progress plots
│   └── utils.R              # Helpers
├── src/
│   ├── likelihood.cpp       # Felsenstein pruning, per-partition likelihood
│   ├── rate_matrix.cpp      # JC(k') construction, analytical P(t) = exp(Qt)
│   ├── tree_moves.cpp       # SPR, NNI proposals in C++
│   ├── branch_moves.cpp     # Branch length proposals
│   ├── mcmc_engine.cpp      # Hot loop: propose → evaluate → accept/reject
│   ├── acrv.cpp             # Among-character rate variation
│   └── RcppExports.cpp      # Auto-generated
├── inst/
│   └── REFERENCES.bib
├── tests/testthat/
├── man/
└── vignettes/
```

---

## Build workflows

### Single-agent phase (Phases 1–3)

During early development with one agent, use standard devtools workflows
in a subprocess:

```bash
# Compile + load (always in subprocess, never RStudio session)
cd mkp
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all()"

# Run targeted tests
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all(); testthat::test_file('tests/testthat/test-foo.R')"

# Regenerate docs
Rscript -e "roxygen2::roxygenise()"
```

Once C++ code exists, switch to:

```bash
Rscript -e "roxygen2::roxygenise(load_code = roxygen2::load_installed)"
```

to avoid debug `.o` contamination.

### Multi-agent phase (Phase 4+)

When multiple agents are active, switch to the renamed-package build
system (`build-agent.sh` / `test-agent.sh`) following the parent
`AGENTS.md` protocol. Add MkPrime support to those scripts at that time.

### Validation

Use **GitHub Actions** for full test suites and R CMD check. Local builds
are for targeted iteration only (build + run 1–2 specific test files).

---

## Build failure recovery

| Symptom | Cause | Fix |
|---------|-------|-----|
| Debug `.o` contamination | `roxygen2::roxygenise()` default uses `debug=TRUE` | `rm -f src/*.o src/*.dll`, rebuild |
| "Access is denied" | Another R process has DLL loaded | Kill the process or wait |
| Namespace errors after rename | `RcppExports.cpp` stale | `Rscript -e "Rcpp::compileAttributes()"` |

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
