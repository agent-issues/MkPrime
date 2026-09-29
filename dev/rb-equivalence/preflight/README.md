# Fixed-state likelihood pre-flight

This check asks: do MkPrime and RevBayes compute the same likelihood for the same state? Run it before any posterior comparison (`compare.R`). If they don't agree, the two samplers target different posteriors, and no amount of MCMC will fix that.

**How it works.** It takes K states (tree, branch lengths, `rate_loss`, `rate_neo`, `rate_log_sd`) from a finished RB cell and re-evaluates them in RevBayes under each model toggle:
- neo Q rescaled, or not;
- `coding = "variable"` or `"all"`;
- polymorphisms as partial ambiguity, or as `?`.

It then evaluates the same states in MkPrime.

```bash
# Locally, from this directory (in Git Bash prefix with MSYS_NO_PATHCONV=1 and
# give local paths as C:/..., or /nobackup/... is rewritten)
Rscript make_probe.R 3408 by_nt_9v --rb-dir=<cell dir> --out-dir=probe
scp probe/* hamilton8.dur.ac.uk:/nobackup/$USER/rt11/probe/
# On Hamilton, in that directory (seconds per cell)
bash run_probe.sh
# Locally, after copying out_*.txt back
Rscript check_probe.R 3408 by_nt_9v --probe-dir=probe
```

`check_probe.R` fails when MkPrime and RevBayes differ by more than `--tol` (1e-4 nats), with polymorphisms coded identically on both sides.

## Columns
- **`asRun`** compares against RevBayes as the harness actually runs it. It is non-zero wherever a matrix has polymorphic cells: MkPrime reads them as missing, RevBayes as partial ambiguity.
- **`depoly` and `noAsc`** compare with polymorphisms coded identically on both sides; these are the columns the pass/fail check uses.

## Recorded results (2026-09-24, RB `build-pr816`, 8 states per cell)

| cell | agent-issues main b1eb0ad | PR #223 (eecee51) |
|---|---|---|
| 950 by_nt_9v | FAIL, 1.6 nats | PASS, 7e-6 |
| 3408 by_nt_9v | FAIL, 4.3 | PASS, 5e-5 |
| 3408 by_nt_kv | FAIL, 16 | PASS, 8e-6 |

On main the failures are #212 (the neo rate scale) and #213 (the ascertainment mask). On 3408 by_nt_kv, the harness's k (#214) makes it worse still. After #223, `asRun` still differs by 1–3 nats on 3408, which is polymorphism handling.

Always align tips before calling `MkpLogLikelihood()`: it pairs tips with data rows by position (#224). `check_probe.R` does this with `RenumberTips()`.
