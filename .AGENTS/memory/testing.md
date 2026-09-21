# MkPrime — testing

Tests live in `tests/testthat/`. ~4000 R-level tests.

## Keeping the log quiet

`test_check("MkPrime")` should print test outcomes and nothing else, so
`tests/testthat/setup.R` pins `options(MkPrime.verbosity = 0)` and attaches
`ape`, `TreeTools` and `shiny` once, quietly. Two rules follow:

 - **Do not call `library()` in a test file.** Everything the suite needs is
   already attached; a repeat `library("TreeTools")` re-emits `ape`'s startup
   banner, and `library("MkPrime")` errors outright ("namespace is already
   attached") under the renamed-DLL harness.
 - **Never work around output with a bare `suppressMessages()`.** Either the
   output is the thing under test -- assert on it -- or it should not be
   produced.

Helpers in `helper-limits.R`:

| Helper | Use when |
|--------|----------|
| `local_mkp_verbosity(level = 1)` | The test asserts on console output. Raises the level for that block only. |
| `allow_warning(expr, regexp)` | A warning is an expected artefact of the fixture (a 50-iteration warmup never stabilises; `Lobo.phy` has invariant characters) rather than a result. Matching warnings are muffled; everything else still surfaces. |
| `expect_prints(expr)` | Checking that a `print` method works. Asserts something reached stdout *or* the message stream, and swallows it -- stronger than `expect_no_error(print(x))`, which a silently-broken method would pass. |

Prefer `expect_warning()` when the warning *is* the result being tested;
`allow_warning()` exists for the cases where asserting it would make the test
fail on the rarer run where the warning does not fire.

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
