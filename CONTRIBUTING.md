# Contributing to MkPrime

## Naming conventions

MkPrime uses **PascalCase for functions** and **camelCase for variables**.
This applies everywhere in the package except for model parameter names,
which are kept in snake_case as domain terminology.

### Functions

| Context | Convention | Examples |
|---------|------------|---------|
| Exported functions | PascalCase | `RunMkPrime`, `MkPrimeData`, `MkPrimeModel`, `MkPrimeMCMC`, `ConvergenceDiagnostics`, `AutoDetectNeomorphic` |
| Internal helpers (no dot) | PascalCase | `ProposeScale`, `ProposeBetaSimplex`, `LogPrior`, `DiscreteLognormalRates` |
| Private helpers (dot-prefixed) | `.PascalCase` | `.InitState`, `.BuildMoves`, `.FinalizeModel`, `.FitchScore`, `.KeyParamCols` |
| S3 methods | `generic.ClassName` | `print.MkPosterior`, `summary.MkPrimeData` |

### Variables and parameters

| Context | Convention | Examples |
|---------|------------|---------|
| Function parameters | camelCase | `knownStates`, `fixTopology`, `treeLengthShape`, `checkEvery`, `checkpointFile` |
| Local variables | camelCase | `charMatrix`, `transIdx`, `hasNeo`, `startTree`, `nSavedPerRun` |
| Loop counters | camelCase or single-letter OK | `iter`, `run`, `ch`, `i`, `j` |

### Model parameter names (exception)

The MCMC state list, sample matrix column names, and related domain terms
use **snake_case** because they are scientific identifiers, not R variable names:

```
tree_length    rate_loss     rate_log_sd    rel_br_lengths
log_posterior  log_likelihood  kPrime        rate_neo
```

These names are stable across the codebase, appear in user-facing output
(`summary()`, `plot()`), and correspond to RevBayes parameter names.
Do **not** rename them to camelCase.

### Worked examples

```r
# GOOD
.InitState <- function(tree, mkd, model) {
  treeLength <- sum(tree$edge.length)
  relBr      <- tree$edge.length / treeLength
  hasNeo     <- any(mkd$type == "neomorphic")
  # state list uses domain terms (snake_case)
  state <- list(tree_length = treeLength, rel_br_lengths = relBr)
  state
}

MkPrimeModel <- function(treeLengthShape = 2, expSteps = NULL) { ... }

RunMkPrime <- function(data, tree, fixTopology = FALSE, knownStates = integer(0)) { ... }

# BAD — do not use these patterns
.init_state <- function(...) { ... }          # snake_case function
run_mk_prime <- function(...) { ... }         # snake_case exported function
MkPrimeModel <- function(tree_length_shape)   # snake_case parameter
charMatrix <- ...                             # PascalCase variable (should be camelCase)
```

## Code style

- Use base R pipe `|>`, not `%>%`.
- Prefer `cli::cli_abort()` / `cli::cli_warn()` over `stop()` / `warning()`.
- Keep comments focused on *why*, not *what*. Avoid restating the code.
- Write tests for all new functionality. Target file: `tests/testthat/test-<topic>.R`.

## Building and testing

See `AGENTS.md` for the build workflow. The short version for local iteration:

```bash
cd mkp
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all(); testthat::test_file('tests/testthat/test-foo.R')"
```

Never use `devtools::load_all()` in the active RStudio session — it locks the
DLL and prevents recompilation. Always build in a subprocess.
