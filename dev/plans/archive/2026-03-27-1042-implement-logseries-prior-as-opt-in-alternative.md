# Plan: Logseries k' prior as opt-in alternative

**Date:** 2026-03-27  
**Branch:** `logseries-kprime-prior`  
**Worktree:** `../mkp-logseries`

---

## Goal

Add a Logseries distribution as an opt-in alternative to the hierarchical
geometric for the k'_i prior on transformational characters.  This matches
the fixed prior used in RevBayes (`Logseries(c=0.7)`), enabling a clean
prior-matched comparison between the two implementations.

---

## Mathematical background

The logarithmic series (log-series) PMF:

```
P(k; c) = c^k / (k · (−log(1−c))),   k = 1, 2, 3, …,   c ∈ (0, 1)
```

Log form (used in `cpp_log_prior` / `LogPrior`):

```
log P(k; c) = k·log(c) − log(k) − log(−log(1−c))
```

The hard constraint `k'_i ≥ kObs_i` is enforced by returning `−Inf` for
invalid states (identical to the geometric treatment).  Because the MH
acceptance ratio divides two evaluations at the same `kObs_i`, the
truncation normalising constant cancels — **no explicit truncation
correction is needed in the MH ratio**.  The `−log(−log(1−c))` term is
included anyway so that the absolute `log_posterior` value reported in the
samples is correct.

`c` is a **fixed model parameter** (not sampled), stored in `MkPrimeModel`,
defaulting to 0.7 to match RevBayes.  There is no hyperparameter `p` when
`kPrimePrior = "logseries"`.

---

## Worktree setup

```bash
cd C:/Users/pjjg18/GitHub/mkp
git worktree add ../mkp-logseries logseries-kprime-prior
# All implementation work is done inside ../mkp-logseries
```

Merge back to `main` only after all tests pass.

---

## Scope of changes

### 1. `R/MkPrimeModel.R`

**`MkPrimeModel()` constructor:**
- Add parameter `kPrimePrior = "geometric"` (`match.arg` with `"logseries"`).
- Add parameter `kprimeLogseriesC = 0.7` (only relevant when
  `kPrimePrior = "logseries"`; ignored for geometric).
- Store both in the returned list.
- The `kprimeHyperA / kprimeHyperB` parameters remain but are only used for
  the geometric prior.
- Add a `cli::cli_warn()` if `kprimeLogseriesC` is supplied alongside
  `kPrimePrior = "geometric"` (unused parameter).

**`LogPrior()`:**
- After the existing boundary checks, add a new boundary:
  `if (model$kPrimePrior == "logseries" && (model$kprimeLogseriesC <= 0 || model$kprimeLogseriesC >= 1)) return(-Inf)`
- For the k' block (currently geometric), branch on `model$kPrimePrior`:
  - `"geometric"`: unchanged (uses `state$p`, Beta hyperprior).
  - `"logseries"`: compute `sum(k'_i * log(c) - log(k'_i)) - nTrans * log(-log(1-c))`.
    No `state$p` access; no Beta term.
- When `kPrimePrior = "logseries"`, skip the `state$p` boundary check entirely.

**`print.MkPrimeModel()`:**
- Update the k' bullet to show which prior is active:
  - `"k' prior: Geometric (Beta hyperprior: {a}, {b})"` or
  - `"k' prior: Logseries (c = {c})"`.

---

### 2. `src/mcmc_state.h` — `McmcData` struct

Add two fields:

```cpp
bool   kPriorLogseries;   // true = log-series; false = hierarchical geometric
double kprimeLogseriesC;  // c parameter (only used when kPriorLogseries = true)
```

---

### 3. `src/mcmc.cpp` — `cpp_log_prior()`

In the `hasTrans` block, replace the current unconditional geometric code
with a branch:

```cpp
if (hasTrans) {
  if (data.kPriorLogseries) {
    // Logseries: log P(k'_i; c) = k'_i*log(c) - log(k'_i) - log(-log(1-c))
    double c = data.kprimeLogseriesC;
    double logC   = std::log(c);
    double logNorm = std::log(-std::log1p(-c));  // log(-log(1-c))
    for (int i = 0; i < nTrans; ++i) {
      int kp = kPrime[data.transIdxGlobal[i]];
      lp += kp * logC - std::log(static_cast<double>(kp)) - logNorm;
    }
    // No p / Beta term
  } else {
    // Hierarchical geometric (existing code)
    if (p <= 0.0 || p >= 1.0) return R_NegInf;
    double sumU = 0.0;
    for (int i = 0; i < nTrans; ++i) {
      int gi = data.transIdxGlobal[i];
      sumU += (kPrime[gi] - data.kObs[gi]);
    }
    lp += nTrans * std::log(p) + sumU * std::log1p(-p);
    lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
  }
}
```

The `p` boundary check that currently appears at the top of `cpp_log_prior`
must also be guarded: skip it when `data.kPriorLogseries == true`.

---

### 4. `src/mcmc_likelihood.cpp` — `prepare_mcmc_data()`

Add two parameters to the signature:

```cpp
bool   kPriorLogseries,
double kprimeLogseriesC
```

Store them:

```cpp
d->kPriorLogseries   = kPriorLogseries;
d->kprimeLogseriesC  = kprimeLogseriesC;
```

Then run `Rcpp::compileAttributes()` to regenerate `RcppExports.cpp`.

---

### 5. `R/RunMkPrime.R`

**`.PrepareData()` (the block that calls `prepare_mcmc_data`):**
- Pass two new arguments at the end:
  ```r
  model$kPrimePrior == "logseries",
  model$kprimeLogseriesC %||% 0.7
  ```

**`.InitState()`:**
- Only set `state$p <- 0.5` when `model$kPrimePrior != "logseries"`.
- When logseries, omit `p` from the state list entirely.

**`.InitMcmcChain()`:**
- `init_mcmc_state(…)` still receives a `p` argument (C++ struct field exists
  for both paths); pass `state$p %||% 0.0` — the value is irrelevant when
  `kPriorLogseries = true`.

**`.BuildMoves()`:**
- Pass `kPrimePrior` as a new argument (or look it up from `model`).
- Only append the `"p"` move when `kPrimePrior != "logseries"`.

**`.ParamNames()`:**
- Pass `kPrimePrior` as a new argument.
- Only include `"p"` in the column-name vector when geometric.

---

### 6. `tests/testthat/test-logseries-prior.R` (new file)

Tests to write:

| # | Test name | What it checks |
|---|-----------|----------------|
| 1 | `MkPrimeModel stores logseries prior settings` | Constructor stores `kPrimePrior = "logseries"` and `kprimeLogseriesC` |
| 2 | `LogPrior logseries matches manual density` | `LogPrior()` equals hand-computed `sum(k*log(c) - log(k)) - nTrans*log(-log(1-c))` + other prior terms |
| 3 | `LogPrior logseries: higher k' has lower prior density` | Prior favours smaller k' for c < 1 |
| 4 | `LogPrior logseries: c out of bounds returns -Inf` | `c = 0` and `c = 1` both return `-Inf` |
| 5 | `LogPrior logseries: k' < kObs returns -Inf` | Boundary constraint enforced |
| 6 | `LogPrior logseries: no p in state` | `LogPrior()` runs correctly when `state$p` is absent |
| 7 | `RunMkPrime smoke test with logseries prior` | Short MCMC run (50 iter, small data) completes without error; samples contain no `p` column |
| 8 | `print.MkPrimeModel shows logseries` | Output contains "Logseries" |

---

## Build sequence

All commands run in `../mkp-logseries`:

```bash
# 1. After C++ changes:
Rscript -e "Rcpp::compileAttributes()"

# 2. Build + load:
Rscript -e "pkgbuild::compile_dll(debug = FALSE); devtools::load_all()"

# 3. Run targeted tests:
Rscript -e "pkgbuild::compile_dll(debug=FALSE); devtools::load_all(); \
  testthat::test_file('tests/testthat/test-logseries-prior.R')"

# 4. Full test suite:
Rscript -e "pkgbuild::compile_dll(debug=FALSE); devtools::load_all(); \
  devtools::test()"
```

---

## Merge & cleanup

```bash
cd C:/Users/pjjg18/GitHub/mkp

# Merge (fast-forward or no-ff as appropriate):
git merge logseries-kprime-prior --no-ff -m "M-NEW: Logseries k' prior as opt-in alternative"

# Remove worktree:
git worktree remove ../mkp-logseries
git branch -d logseries-kprime-prior
```

---

## Files changed summary

| File | Change type |
|------|-------------|
| `R/MkPrimeModel.R` | Modify |
| `R/RunMkPrime.R` | Modify (5 locations) |
| `src/mcmc_state.h` | Modify |
| `src/mcmc.cpp` | Modify |
| `src/mcmc_likelihood.cpp` | Modify |
| `src/RcppExports.cpp` | Auto-regenerated |
| `tests/testthat/test-logseries-prior.R` | New |

---

## Out of scope

- Free (sampled) `c` parameter — the RevBayes comparison needs fixed `c`.
- `"logseries"` option in `EasyMkPrime` — can be added in a follow-up.
- Vignette updates — deferred until this feature is merged.
