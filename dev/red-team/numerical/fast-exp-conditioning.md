# Lane N1 — `fast_neg_exp` + closed-form `P(t)` conditioning audit

## What this audits

- `src/fast_exp.h::fast_neg_exp` — a Taylor-based fast `exp(x)` for `x ≤ 0`,
  claimed to match `std::exp` within ~1e-15 relative error.
- The downstream **closed-form** transition-probability formulas it feeds:
  - JC(k): `src/rate_matrix.cpp:32`, `src/likelihood.cpp:84,213`,
    `src/gibbs_partial_cl.h:208`, `src/mcmc_likelihood.cpp:211,378,511,704,984,1150`,
    `src/node_cl_cache.h:197`, `src/acrv.cpp:87,205`,
    `src/ascertainment.cpp:71,160,254,358`.
  - MkN (2-state asymmetric): `src/rate_matrix.cpp:78`, `src/likelihood.cpp:312`,
    `src/gibbs_partial_cl.h:235`, `src/mcmc_likelihood.cpp:1399,1574`,
    `src/node_cl_cache.h:220`, `src/acrv.cpp:315`,
    `src/ascertainment.cpp:462,560`.
  - F81 (k-state asymmetric): `src/gibbs_partial_cl.h:260`,
    `src/mcmc_likelihood.cpp:1867,2051,2180`.

**No eigendecomposition path exists** in this codebase — the lane's "eigenvalue
near-degeneracy" axis turns into an analytic-formula-cancellation axis.

## Conditioning analysis

For all three models, P(t) reduces to a single scalar `E = exp(−α t)` and
linear combinations of `E` and `1 − E`. The two regimes of concern:

**(R1) `αt → 0`:** Computing `1 − E` directly catastrophically cancels because
`E` rounds to 1 in double precision once `αt ≲ 2.2e-16`. For `αt` in
`[1e-16, 1e-8]`, `1 − E` loses roughly `−log₁₀(αt)` significant digits. The
JC `p_diff = (1/k)·(1 − E)` and F81 `1 − d` inherit this loss. Standard fix:
`1 − exp(x) = −expm1(x)`, which is exact to ulp for any `x ≤ 0`.

**(R2) `αt → ∞`:** `E → 0` (underflow). `fast_neg_exp` gates at `x < −708`
and returns 0; `std::exp` returns subnormal for `x ∈ [−744, −708]`. The
difference is bounded by `2.2e-308`, with zero downstream impact: the
formulas reduce to `P_ij → π_j` cleanly in both cases.

**(R3) `fast_neg_exp` interior accuracy:** degree-11 Taylor on a reduced
interval `|r| ≤ ln(2)/2 ≈ 0.347` has truncation error `≤ 0.347¹²/12! < 1e-17`.
Argument reduction with split `LN2 = LN2_HI + LN2_LO` controls cancellation
in `r = x − n·ln(2)` to about half-ulp. The IEEE-754 exponent injection at
the end is exact. Expected `ε_rel ≲ 5e-15`.

**(R4) `fast_neg_exp` argument-reduction quirk:** the code uses
`int n = static_cast<int>(nd − 0.5)` to round `nd = x·log₂(e)` to the nearest
integer for `x ≤ 0`. This is **truncate-after-shift** rather than true
round-to-nearest. For `nd ∈ {…, −0.5, 0.5, …}` the two differ by 1 but
either choice keeps `|r| < ln(2)`, so the polynomial remains in its
accuracy window. Verified numerically (Test 1).

## Reference computation

`mpmath` at `dps=50` (≈ 166 binary digits). `Rmpfr` was not available in this
worktree's R install; the audit's R driver shells out to a Python driver that
re-implements `fast_neg_exp` **bit-exactly** in Python via `struct.pack`
IEEE-754 reinterpretation, so the reported errors are the actual errors of
the C function. Verified bit-equivalence against a hand-computed table of
known points.

## Stress-test design

Five tests on grids spanning the four conditioning regions:

| Test | Routine | Grid |
|------|--------|------|
| 1 | scalar `fast_neg_exp` vs `std::exp` | 50k log-uniform on `[1e-12, 708]` plus the boundary points `{−708, −709}` |
| 2 | JC(k) `(p_same, p_diff)` | `k ∈ {2,4,6,10}` × `rt ∈ {1e-12, …, 1e4}` |
| 3 | MkN(2) `P` | `rate_loss ∈ {1, 0.1, 0.01, 100, 1e6}` × `t ∈ {1e-12, …, 1e4}` |
| 4 | F81(k) `P` with ACRV | `k ∈ {2,4,6,10}` × `π` asymmetry `{1, 3, 100, 1e6}` × `σ ∈ {0.1, 1, 3, 10}` × `K=4` cats × `t ∈ {1e-10, …, 100}` |
| 5 | cancellation of `1 − E` | `k ∈ {2,4,6,10}` × `rt ∈ {1e-16, …, 1}` against `−expm1(arg)` fix |

Grid spans more than three points per axis (lane spec minimum).

## Results

CSVs under `dev/red-team/numerical/fast-exp-results/`; summary
in `summary.txt`. Headline numbers from the full grid:

| Test | Quantity | Worst `ε_rel` | Where | Verdict |
|------|---------|---------------|-------|--------|
| 1 | `fast_neg_exp` vs `mpmath`, x ≥ −708 | **8.6e-15** | `x ≈ −0.347` (Taylor edge) | Matches `<1e-14`, slightly above the in-file `<1e-15` claim (header comment is optimistic) |
| 1 | `std::exp` vs `mpmath` | 1.6e-16 | — | baseline ulp |
| 1 | underflow boundary `x ∈ [−744, −708]` | `fast` returns 0, `std::exp` returns subnormal `≤ 2.2e-308` | by design | benign |
| 2 | JC `p_diff` (current, full grid) | **5.0e-5** | `k=10, rt=1e-12` | algebraic cancellation in `1−E` |
| 3 | MkN `P` entries (current) | 4.4e-5 | `rate_loss=0.1, t=1e-12` | same cancellation |
| 4 | F81 `P` (current, σ=10 corner) | **1.0** (total loss) | k=2, t·μ ≈ 8e-20 | same cancellation |
| 5 | JC `p_diff` current, `rt ∈ [1e-16, 1]` | **3.9e-1** | `rt ≤ 1e-15` | cancellation, full digit loss |
| 5 | JC `p_diff` with `−expm1(arg)` | **2.4e-16** | — | restored to ulp |
| 5 | F81 `1−d` current | 1.7e-1 | small `μt` | cancellation |
| 5 | F81 `1−d` with `−expm1(arg)` | 1.8e-16 | — | restored to ulp |

### `fast_neg_exp` itself is fine

In the regime it claims to cover (`x ∈ [−708, 0]`), `fast_neg_exp` agrees
with `std::exp` to within `8.6e-15` relative error over 50k random samples
(50 evenly spaced + boundary points + log-uniform random). This is below
the lane's 1e-10 acceptance threshold by ~5 orders of magnitude. The in-file
header claim of `<1e-15` is mildly optimistic (true bound is closer to
1e-14); for safety I would document the bound as `≤ 1e-14`, but no
correctness impact follows.

### The cancellation in `1 − E` is **NOT** a `fast_neg_exp` bug

In Test 2/3, **`std::exp` and `fast_neg_exp` produce identical relative
errors on `p_diff`** because both return `E` to ulp accuracy and the
subtraction `1/k − (1/k)·E` is what loses digits. Recompiling with
`-DMKP_NO_FAST_EXP` (the validation flag in `fast_exp.h:81`) would not
mitigate any of the cancellation findings.

### Production reachability

Per `R/MkPrimeModel.R:113`, the default `rateLogSd` prior is `Gamma(1, 1)`
(exponential mean 1). σ ≥ 7 is in the tail but **not bounded out**. With
K=6 ACRV cats, the slowest category midpoint at σ=10 has rate
`≈ exp(−σ²/2 − σ·1.38) ≈ 1.3e-16`. Combined with a typical branch
`t ≈ 1e-2`, the product `rt ≈ 1.3e-18` — solidly in the regime where
`1 − E = 0` and `p_diff = 0`, causing site contributions on
state-discordant tip pairs to evaluate to `log(0) = −∞` and the whole
chain's likelihood to collapse to `−Inf`.

At more moderate σ ≤ 5 (the bulk of the posterior under the default
prior), the smallest `rt` is ~1e-9 and `p_diff` relative error stays at
`~1e-7` — below the per-site logL noise floor of a 30-taxon dataset but
not below the 1e-10 lane threshold.

## Verdict

**Stable with caveats.**

- `fast_neg_exp` itself meets the audit threshold (`ε ≤ 1e-14` versus the
  1e-10 acceptance) across its claimed operating range. The `<1e-15`
  header comment is mildly aspirational; `≤ 1e-14` is the demonstrated
  bound.
- The underflow cutoff at `x < −708` is benign (downstream P entries are
  the same in either case).
- The **algebraic cancellation** in the JC/MkN/F81 closed forms
  (`p_diff = inv_k − inv_k·E`, `1 − d` in F81) exceeds the 1e-10 threshold
  for `rt ≲ 1e-10` and reaches full digit loss for `rt ≲ 1e-15`. This is
  reachable under default priors when `rateLogSd ≳ 7` combines with short
  branches — a tail of the posterior, but not bounded out.
- This is **not a `fast_neg_exp` bug**: the same cancellation occurs with
  `std::exp`.

## Recommendation

Trivial patch under `dev/red-team/patches/N1-fast-exp.patch` rewrites the
five primary call sites to use `std::expm1`:

```cpp
// Before (loses up to all digits at small αt):
double exp_term = MKP_EXP(-α * t);
double p_diff   = inv_k - inv_k * exp_term;

// After (ulp-accurate everywhere):
double neg_expm1 = -std::expm1(-α * t);   // = 1 - exp(-α t), no cancellation
double exp_term  = 1.0 - neg_expm1;        // recovered, still ulp-accurate for E
double p_diff    = inv_k * neg_expm1;
double p_same    = inv_k + (1.0 - inv_k) * exp_term;
```

For F81 / MkN the same pattern applies (replace `1.0 - exp_term` with
`-std::expm1(arg)`). Performance: `std::expm1` is ~1.5x slower than
`std::exp` and ~3-4x slower than `fast_neg_exp`. **Use it only at the
single site per call where `1 − E` is needed**, keeping `MKP_EXP` for the
`E`-as-such uses. Net hot-path cost is one extra `expm1` per branch per
character class — small relative to the per-character pruning loop.

The patch is *not* committed; the auditor lane doesn't write to `src/`.

## Worktree state

Audit worktree is clean of `src/` modifications. Output files:

- `dev/red-team/numerical/fast-exp-conditioning.md` — this note
- `dev/red-team/numerical/fast-exp-driver.R` — R wrapper (Rscript entry)
- `dev/red-team/numerical/fast-exp-driver.py` — Python driver (does the
  work; called by the R wrapper)
- `dev/red-team/numerical/fast-exp-results/` — CSVs + `summary.txt`
- `dev/red-team/patches/N1-fast-exp.patch` — proposed mechanical fix
