# Running `RunMkPrime()` from a bash subprocess

`RunMkPrime()` with convergence criteria can run for hundreds of thousands
of iterations. On Windows, when the bash tool timeout fires, the child
Rscript process is **not killed** — it becomes an orphan consuming CPU.
Repeated retries accumulate orphans that must be killed manually.

## Triple-guard every subprocess `RunMkPrime()` call

1. **`maxTime`** — always set in the `RunMkPrime()` call (e.g. `maxTime = 120`).
   Application-level stop: MCMC finishes the current batch and exits cleanly.
2. **`setTimeLimit()`** — wrap the R code in `setTimeLimit(elapsed = 150)` as
   a backup in case `maxTime` isn't checked frequently enough.
3. **Shell `timeout`** — wrap the Rscript with `timeout Ns Rscript -e "..."`.
   Rtools45 ships GNU `timeout`. This is the outermost guard.
4. **Tool timeout** — set the bash tool `timeout` parameter above the shell
   `timeout` so the shell layer fires first (e.g. shell 60s, tool 65 000 ms).

## Template

```bash
taskkill //F //IM Rscript.exe 2>/dev/null; sleep 1
timeout 180 Rscript -e '
  setTimeLimit(elapsed = 150, transient = FALSE)
  pkgbuild::compile_dll(debug = FALSE); devtools::load_all()
  # ... setup ...
  posterior <- RunMkPrime(..., maxTime = 120, ...)
  # ... diagnostics ...
' 2>&1
echo "EXIT: $?"
```

## Rule

**Never** run `RunMkPrime()` with convergence criteria (`minEss`, `minTreeEss`,
`maxRhat`) without also setting `maxTime`. The `nIter` cap alone is not
sufficient — warmup can consume the entire iteration budget.

## See also

Parent `AGENTS.md` → **Subprocess timeout discipline** for the general rule
(which this is a concrete instance of).
