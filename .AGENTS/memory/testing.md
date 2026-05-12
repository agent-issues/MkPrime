# MkPrime — testing

Tests live in `tests/testthat/`. ~4000 R-level tests.

## Test budgets (Windows runaway-prevention)

Per parent `AGENTS.md` → **Subprocess timeout discipline**, every Rscript
subprocess must be triple-guarded. For tests that run MCMC:

| Context | Max nIter | Rationale |
|---------|-----------|-----------|
| Single-run unit test | ≤ 1 000 | Fast iteration during development |
| Multi-run (nRuns ≥ 2) test | ≤ 500 | Multiplied by nRuns |
| Stopping-rule tests | Use `maxTime = 0.5` | Self-limiting wall-clock |

**Expensive optional computations** (tree-ESS via distance matrices, etc.)
should be gated behind an env var (e.g. `MKP_TREE_ESS_TESTS=1`) and skipped
by default.

## Targeted test runs

During development, run only the tests you need. For the renamed-DLL
multi-agent setup (parent `AGENTS.md`):

```bash
bash test-agent.sh mkp <Letter> test-foo
```

For single-agent direct iteration:

```bash
Rscript -e "
  setTimeLimit(elapsed = 90, transient = FALSE)
  pkgbuild::compile_dll(debug = FALSE); devtools::load_all()
  testthat::test_file('tests/testthat/test-foo.R')
"
```

**Never** run the full test suite locally — use GHA (`agent-check.yml`).

## Convergence-criterion tests

Any test that uses `minEss`, `minTreeEss`, or `maxRhat` **must** also pass
`maxTime`. Without `maxTime`, warmup can consume the entire `nIter` budget
and the test hangs.
