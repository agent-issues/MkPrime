# MkPrime — naming conventions

| Scope | Convention | Examples |
|-------|------------|---------|
| Exported functions | **PascalCase** | `RunMkPrime`, `MkPrimeData`, `ConvergenceDiagnostics` |
| Internal functions (dot-prefixed) | **`.PascalCase`** | `.InitState`, `.BuildMoves`, `.FinalizeModel` |
| Non-exported helpers (no dot) | **PascalCase** | `ProposeScale`, `LogPrior`, `DiscreteLognormalRates` |
| Function parameters | **camelCase** | `knownStates`, `fixTopology`, `treeLengthShape`, `checkEvery` |
| Local variables | **camelCase** | `nEdge`, `transIdx`, `charMatrix`, `startTree` |
| Model parameter names (MCMC state, column names) | **snake_case** | `tree_length`, `rate_loss`, `rate_log_sd`, `log_posterior` |
| S3 class names | **PascalCase** | `MkPrimeData`, `MkPrimeModel`, `MkPosterior` |

## Why model parameters are snake_case

`rate_loss`, `tree_length`, `rate_log_sd`, `kPrime`, `rel_br_lengths`, etc.
are domain terminology that appears in:

- Output column names (Tracer-compatible TSV logs)
- User-facing documentation
- RevBayes cross-references (RevBayes uses snake_case for state vars)

Keeping these as snake_case in the R/C++ code avoids translation at I/O
boundaries and matches user expectations from related tooling.

## Full rationale and worked examples

See `CONTRIBUTING.md` in the package root.
