# ess-rhat-stress — Lane D3 heavy test

## What this tests

Empirically validates the correctness of MkPrime's effective sample size (ESS)
and rank-normalised split-chain R-hat implementations against adversarial MCMC
chains whose true integrated autocorrelation time (IAT) is *analytically*
known. The targets:

- `R/ess.R::.Ess`, `.EssVector`, `.EssMatrix`, `.EssMultiChain` — Geyer (1992)
  initial-positive / initial-monotone sequence + Vehtari et al. (2021)
  truncated estimator on FFT-based autocovariances.
- `R/Convergence.R::.Rhat`, `.ComputeRhat` — rank-normalised, split-chain,
  fold-and-take-max R-hat per Vehtari et al. (2021).
- `R/Convergence.R::.ComputeRhat` — regression for **CONV-002** (cbind-recycle
  on unequal-length per-run matrices).
- All of the above — regression for the **CONV-001 family** (constant chain
  must not yield `Rhat <= 1.01`).

This complements the existing `data-raw/step5_ess_estimator_check.R` AR(1)
sanity check by widening to AR(2), bimodal mixtures, near-non-stationary
random-walks, and the diagnostic-regression cases that the orchestrator flagged.

## Theoretical basis

| Process | True diagnostic | Source |
|---|---|---|
| AR(1)  `x_t = rho x_{t-1} + eps_t` | IAT = `(1+rho)/(1-rho)`, true ESS = `N(1-rho)/(1+rho)` | standard textbook (e.g. Geyer 1992 §3.3) |
| AR(2)  `x_t = phi1 x_{t-1} + phi2 x_{t-2} + eps_t` | IAT = `1 + 2 sum_{k>=1} rho_k`, ACF from Yule-Walker recursion | Brockwell & Davis (1991) §3 |
| Bimodal N(+/-mu,1), 2 chains each pinned to one mode | per-mode ESS = `N` (iid), cross-chain R-hat >> 1 (Vehtari §4.4) | Vehtari et al. (2021) |
| Near-non-stationary RW Metropolis (acceptance ~1%) | ESS/N << 0.1; cross-chain Rhat > 1.01 | qualitative |
| Constant chain (`x_t = c`) | between-chain var = 0, within-chain var = 0 -> Rhat is undefined; ESS undefined | Vehtari et al. (2021) §3 |

For the AR(1) sticky regime (rho >= 0.99), the truncated estimator is known to
be biased *low* because the initial-monotone-sequence truncation cuts off the
slowly-decaying autocorrelation tail before it integrates fully. This is a
documented property of Geyer's truncation; we therefore widen the per-rho
tolerance accordingly.

## Pass criterion

Per-scenario, relative-error or threshold-based:

| Scenario | Quantity | Tolerance | Status rule |
|---|---|---|---|
| AR(1) rho in {0.1, 0.5} | mean(measured ESS) vs `N(1-rho)/(1+rho)` | 10% | PASS if rel-err <= 0.10 |
| AR(1) rho = 0.9 | same | 10% | PASS if rel-err <= 0.10 |
| AR(1) rho = 0.99 | same | 30% | PASS if rel-err <= 0.30 |
| AR(1) rho = 0.999 | same | 60% | PASS if rel-err <= 0.60 (estimator strongly biased low) |
| AR(2) | mean(measured ESS) vs `N / IAT(phi1,phi2)` | 15% | PASS if rel-err <= 0.15 |
| Bimodal per-mode ESS | measured vs `N` | 15% | PASS if rel-err <= 0.15 |
| Bimodal cross-chain Rhat | measured | n/a | PASS if Rhat >= 1.1 |
| Near-non-stationary ESS | measured / `N` | n/a | PASS if measured/N < 0.10 |
| Near-non-stationary cross-chain Rhat | measured | 1.01 | PASS if Rhat > 1.01 |
| Sticky AR(1) (rho=0.99) | measured ESS vs truth | 30% | PASS if rel-err <= 0.30 |
| Constant `.EssVector` | output | n/a | PASS iff `NA` |
| Constant `.Rhat` (two const chains, identical) | output | n/a | PASS iff `NA` / `Inf` (i.e. not finite & <= 1.05) |
| Constant `.Rhat` (two const chains, different consts) | output | n/a | PASS iff `NA` / `Inf` / > 1.01 |
| `.ComputeRhat` w/ const per_run | each element | n/a | PASS iff all `NA` / `Inf` / > 1.01 |
| Unequal-length per_run (CONV-002) | warning/error emission | n/a | PASS iff no warning, no error, finite Rhat returned (post-fix expectation). **Pre-fix this row is expected to FAIL** with `warning = "number of rows of result is not a multiple of vector length"` — the test is a regression detector that turns green once `.ComputeRhat` performs tail-equalisation. |
| ESS-cap probe (white noise, small N) | measured/N | 1.5 | PASS if measured <= 1.5 N; WARN otherwise (potential ESS-CAP-001 finding) |

Overall verdict = `FAIL` if any row fails; else `WARN` if any row warns;
else `PASS`. Written as the last line of `verdict.txt` for grep-ability.

### Tolerance rationale

- **AR(1) low/mid rho (10%)**: Geyer's truncated estimator is approximately
  unbiased and consistent for stationary AR(1). 20 replicates of N=10^4 give
  measurement sd around 7%; 10% is a comfortable threshold without inviting
  false positives.
- **AR(1) high rho (30% / 60%)**: at rho=0.99, the true IAT is 199, so true
  ESS for N=10^5 is around 500. The truncated estimator cuts the geometric
  tail short — bias is documented to be 10-30%. For rho=0.999 (true IAT
  around 2000), only N=10^5 has enough samples to estimate IAT meaningfully,
  and the truncated estimator's bias can exceed 50%. These tolerances surface
  *new* problems (e.g. the cap firing inappropriately) without falsely
  flagging known bias.
- **Bimodal R-hat threshold 1.1**: Vehtari et al. recommend Rhat <= 1.01 for
  reliable inference. Two chains stuck in modes separated by 5sigma yield
  *huge* Rhat (often > 5); 1.1 is a generous lower bound for "diagnoses
  non-mixing".

## How to run

- **Quick** (<= 90s laptop, used by orchestrator to verify execution):
  ```
  Rscript dev/red-team/heavy-tests/ess-rhat-stress.R --quick
  ```
- **Full grid** (a few minutes on a laptop; no HPC needed):
  ```
  Rscript dev/red-team/heavy-tests/ess-rhat-stress.R
  ```

  Runs all rho x N combinations with 20 replicates each, all AR(2) cases, and
  all diagnostic-regression scenarios.

## Output interpretation

Outputs land in `dev/red-team/heavy-tests/ess-rhat-stress-results/`:

- `verdict.txt` — one row per scenario, last line is `OVERALL: PASS/FAIL/WARN`.
- `results.csv` — same content as a CSV (scenario, parameter, n, measured,
  truth, tolerance, rel_err, status, note).
- `summary.rds` — the same data.frame for programmatic consumption.

Quick mode covers a subset (one or two AR(1) rho, one AR(2) case, smaller N
for bimodal & sticky); the full grid covers the full matrix in the table
above.

## What a failure would mean

- **AR(1) low/mid-rho FAIL**: most likely cause is a regression in the
  FFT-based autocovariance (`.Autocovariance` in `R/ess.R`) or in the
  initial-positive-sequence loop (rho_even/rho_odd accounting). Less likely:
  an off-by-one in the truncation index `maxT`.
- **AR(1) high-rho FAIL (way beyond widened tolerance)**: indicates the
  safety cap (`tauHat < 1/log10(n) -> ess = n*log10(n)`) is firing
  inappropriately for sticky chains; or the truncation is firing far too
  early. Either is an estimator-tuning bug.
- **Bimodal cross-chain Rhat FAIL** (Rhat falsely <= 1.1): rank-normalisation
  bug or `.RhatClassical` denominator bug (varWithin computed across the
  pooled draws instead of per-chain).
- **Constant-chain FAIL** (Rhat <= 1.05 instead of NA/Inf): regression of
  the CONV-001 family. `.IsConstant` not catching the case, or
  `.RhatClassical`'s `varWithin == 0` guard removed.
- **CONV-002 FAIL** (warning emitted or error): `.ComputeRhat` still uses
  `do.call(cbind, ...)` on per-run vectors of differing length. Required fix
  is tail-equalisation (truncate to `min(nrow)` across runs).
- **ESS-cap probe WARN**: `.Ess`'s safety cap returning `n * log10(n)` for
  white-noise chains. Not a hard FAIL because the cap is intentional, but
  worth surfacing as a potential `ESS-CAP-001` finding (cap is too aggressive
  on short, well-mixed chains).

## Constraints honoured

- No commit, no push, no merge.
- Outputs land under `dev/red-team/heavy-tests/ess-rhat-stress-results/`
  (gitignored at the path level; not committed).
- Does not modify `R/`, `src/`, `tests/testthat/`.
- Uses `pkgload::load_all()` to access internal `.Ess`, `.Rhat`, `.ComputeRhat`
  (intentionally not exported).
- No HPC required — full grid runs in a few minutes on a laptop. No long-burst
  MCMC, no large /nobackup writes.
