# Plan: ESS-over-time panel + exclude invariant parameters

## Summary

Two related changes to live MCMC progress display:

**A.** Add a final ESS-over-time panel to `MkpTracePlot()`, with
per-parameter color coordination between trace panels and the ESS panel.

**B.** Exclude `rate_loss` from sample matrix, logs, and monitors when
there are no neomorphic characters (it's invariant in that case).

---

## Part A: ESS-over-time panel

### Current state

`MkpTracePlot()` shows one trace panel per key parameter
(`log_posterior`, `tree_length`, `rate_loss`, `rate_log_sd`, `p`).
All panels use the same color ("steelblue" for 1 run, per-run palette
for multiple runs). No ESS information is displayed in the plot.

### Design

1. **Stateful ESS history** via a package-level environment
   `.tracePlotEnv` in `PlotDuringMCMC.R`:
   - `essHistory`: list of named numeric vectors (one per callback call),
     each containing ESS for key params at that iteration
   - `iterHistory`: integer vector of iter values (x-axis)
   - `lastIter`: track for reset detection

   **Reset** when `info$inWarmup == TRUE` or `info$iter <= lastIter`
   (indicates new MCMC run). Requires `coda` — if missing, skip the
   ESS panel silently.

2. **Per-parameter color palette**: Assign each key parameter a
   distinct color (e.g., from `hcl.colors(n, "Dark 3")`).

   - **nRuns == 1**: Use the parameter's color as the line color in
     its trace panel AND the ESS panel. This directly fulfils the
     user's "same colour" requirement.
   - **nRuns > 1**: Trace panels continue using per-run colors (users
     need to distinguish runs). ESS panel uses per-parameter colors
     with a legend.

3. **New `.PlotEssPanel()` helper**: Draws ESS trajectories from
   `.tracePlotEnv$essHistory`. One line per parameter, colored to
   match. Compact legend (param names) in corner. Optional horizontal
   reference line at `minEss` target if available.

4. **Layout change**: `nPanels` incremented by 1 (the ESS panel) when
   post-warmup samples exist and `coda` is available.

5. **ESS computation**: Use `coda::effectiveSize()` on the combined
   `runSamples` for the key scalar parameters (excluding `br_*` and
   `kPrime_*`). This mirrors what `.CheckConvergence()` computes but
   is self-contained within the plot callback.

### Changes

| File | What |
|------|------|
| `R/PlotDuringMCMC.R` | Add `.tracePlotEnv`, `.PlotEssPanel()`, modify `MkpTracePlot()` to compute/store/plot ESS, assign per-parameter colors |

---

## Part B: Exclude invariant `rate_loss`

### Current state

`rate_loss` is always present in sample output:
- C++ (`mcmc.cpp` line 1791): `row[col++] = s0->rateLoss;` — always written
- C++ (`mcmc.cpp` line 1723): `nScalarCols = 5 + ...` — the `5` includes
  rate_loss
- R (`.ParamNames()`): always includes `"rate_loss"`
- R (`.StateToRow()`): always includes `state$rateLoss`
- R (`brColStart`): computed as `5L + pCols + hasNeo + qHet + nTrans + 1L`

When there are no neomorphic characters, `rate_loss` is never proposed
(no move is added by `.BuildMoves()`), so it flatlines at the initial
value (1.0). It adds a useless column to output and a confusing trace
panel.

### Design

Make `rate_loss` conditional on `hasNeo`, following the same pattern
used for `rate_neo`, `p`, and `beta_scale`.

### Changes

| File | What |
|------|------|
| `src/mcmc.cpp` | Make `rateLoss` in sample row conditional on `hasNeo`. Change base `nScalarCols` from `5` to `4`, add `+ (hasNeo ? 1 : 0)`. |
| `R/RunMkPrime.R` (`.ParamNames()`) | Only include `"rate_loss"` when `any(mkd$type == "neomorphic")`. |
| `R/RunMkPrime.R` (`.StateToRow()`) | Conditionally include `state$rateLoss` (like `rateNeoVal`). |
| `R/RunMkPrime.R` (`brColStart`) | Change from `5L + ...` to `4L + hasNeo + ...` (rate_loss now counts as part of the hasNeo block). |
| `R/RunMkPrime.R` (`.BuildProgressInfo`) | `rate_loss` in `currentState` only when hasNeo. |
| `R/PlotDuringMCMC.R` | Remove `"rate_loss"` from hardcoded `keyParams` list — it simply won't be in `colnames()` when absent. |
| `tests/` | Update any tests that assume rate_loss is always present. |

**Not changed:**
- C++ `McmcState` struct: `rateLoss` field remains (always initialized
  internally). Only the output serialization changes.
- Checkpoint format: `rate_loss` remains in checkpoint data for
  backward compatibility. On restore, if present it's loaded; if
  absent, defaults to 1.0.
- `.InitState()`: still sets `rate_loss = 1.0` internally.

---

## Task order

1. **Part B first** (invariant exclusion) — this reduces the set of
   parameters that Part A needs to handle, and avoids writing ESS
   code that then needs to be modified.

2. **Part A** (ESS panel) — builds on the clean parameter set.

3. **Test** — rebuild, run test suite, manual interactive check.

4. **Commit** — single commit with clear message.
