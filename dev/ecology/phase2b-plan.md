# Phase 2b wiring plan — replace the rate model (R-level), EBE v1

The KERNEL is already EBE-only (1A; GATE 1 green). The ecology config params
(`phi`/`pi0`/`theta`/`z`, `magnitudeMode`, `sigmaPhi`, `rho0*`, `theta*`, `gibbsZEvery`)
and the `LogPrior` ecology block are REUSED unchanged — there is **no rate-vs-equilibrium
toggle to remove** (`ecologyAware=TRUE` now simply *is* EBE). So Phase 2b is bounded:

## 1. Neomorphic-only `z` (the substantive piece — correctness, not cosmetics)
Under EBE, transformational chars get NO ecology effect (1A neutralized the JC tilt, R8).
But `z` is still initialised over ALL chars (`RunMkPrime.R:3409`,
`matrix(0L, nrow = mkd$nChar, ...)`), and the prior counts every row. Consequence: trans
z rows are inert in the LIKELIHOOD but each contributes `nTrans·log(pi0)` per ecology
column to the PRIOR → biases the `pi0` (sparsity) posterior toward high sparsity (~59
phantom z=0 rows on the rodent 160neo/59trans set) and pollutes the "fraction flagged"
interpretation; gibbs_z would also sample trans z from the prior (spurious flags).
**Fix:** restrict the z PRIOR and the gibbs_z SWEEP to neomorphic character rows.
- `R/MkPrimeModel.R` LogPrior z-loop (L749-768): count only neomorphic rows. Needs the
  neomorphic index set available in `mkd` (mkd$type == "neomorphic").
- `src/mcmc.cpp` `cpp_log_prior` ecology z block (~L407-430, 579+): MIRROR the R change —
  keep R↔C++ bit-identical. (This is the R5/EG-001-style sync discipline.)
- `src/mcmc.cpp` `gibbs_z_sweep_impl` (~L4659): sweep only neomorphic chars.
- Likelihood already ignores trans z (post-1A); the full nChar×(kEco-1) zMat can stay as
  the storage shape (trans rows pinned 0), OR shrink to neomorphic rows — decide by which
  keeps kernel z-indexing simplest. `.PruningMknEcology` already takes the neomorphic
  partition's z rows; `.MkpEcologyLogLikelihood` takes the full matrix and tilts only neo.
- **Verify with T6 prior-invariance + a pi0-posterior check** (trans-count must not move pi0).

**VERIFIED 2026-05-29 (read-only against code):**
- Issue is REAL. R z-prior loop `MkPrimeModel.R:749-768` counts `sum(z[,j]==0L)` etc. over
  ALL nChar rows; C++ `cpp_log_prior` `mcmc.cpp:602-617` iterates `c in 0..nCharZ` (=
  `zMat.nrow()`) for each ecology col. Both pool trans rows into the shared pi0/theta posterior.
- Neomorphic identity lives per-PARTITION: `PartInfo.type` (`mcmc_state.h:14`, 0=neo/1=trans/
  2=known) + `PartInfo.globalCharIdx` (local→global kPrime index). `zMatrix` is GLOBAL
  nChar×(kEco-1). So the fix needs a **neomorphic-global-char mask** built from the partitions
  (`isNeo[globalCharIdx[c]] = (part.type==0)`), then skip `!isNeo[c]` rows in the prior loop
  (R + C++) and the gibbs_z sweep. **KEEP the full nChar zMatrix shape (mask, do NOT shrink)** —
  preserves the kernel's global-char z-indexing (likelihood looks up z[globalChar, eco]).
- RESOLVED: (a) `gibbs_z_sweep_impl` (`mcmc.cpp:4715`) iterates ALL nChar rows
  (`for c in 0..nChar`). For trans/known chars the likelihood is z-invariant (R8), so
  `ll[0]==ll[1]==ll[2]` and the cell is sampled PURELY from the {log_pi0, log_enc, log_disc}
  prior (L4724-4726) → spurious z≠0 flags on non-ecology chars that also dilute the shared
  pi0/theta. (b) R `mkd$type` is a length-nChar char vector (`MkPrimeData.R:48,223`); the
  mask is just `which(mkd$type == "neomorphic")`. C++ has `data->charToPartition[c]` +
  `data->parts[pi].type`.

**EXACT EDITS (3, R↔C++ mirrored, all guarded on neomorphic == type 0):**
1. `cpp_log_prior` count loop (`mcmc.cpp:604`): build `isNeo` once from
   `data.parts[*].type==0` via `globalCharIdx` (or use `data.charToPartition[c]`), then
   `if (!isNeo[c]) continue;` inside `for (c=0; c<nCharZ; ++c)`.
2. `gibbs_z_sweep_impl` (`mcmc.cpp:4715`): at top of the `for (c...)` body, after
   `int pi = data->charToPartition[c];`, add `if (pi < 0 || data->parts[pi].type != 0) continue;`
   — leaves trans/known z rows inert at init 0 (no likelihood change; removes spurious
   flags + 3×per_char_log_lik calls).
3. R `LogPrior` z-loop (`MkPrimeModel.R:752-768`): `neoRows <- which(mkd$type=="neomorphic")`
   (guard length>0), then `zCol <- z[neoRows, j]`.
Net: trans/known z rows fully inert (likelihood R8 + prior masked + gibbs skipped). Keep
the full nChar zMatrix shape.

## 2. Docs / interpretation (no behaviour change)
`R/MkPrimeModel.R` roxygen (L83-158) and `print.MkPrimeModel` (L803-809) describe
`ecologyAware` as a "per-edge, per-character **rate** modifier". Reword to Ecology-Biased
Equilibrium: ecology tilts the F81 **equilibrium** of neomorphic binary chars
(`logit π1 += s_z·log φ`, decay λ fixed at 2), NOT the rate. `magnitudeMode` doc (L89):
"rate-modifier magnitude" → "equilibrium-tilt magnitude (log-odds shift log φ)".

## 3. Retire the 6 rate-model tests
`tests/testthat/test-ecology-likelihood.R` has 6 tests asserting EBE == old rate-model R
oracle (now correctly failing). Delete/relocate them so GATE 2's signal is clean
(advisor). Keep the EcologyNodeMarginals/EcologyEdgeWeights plumbing tests (still valid).

## 4. C++ dead-code cleanup (optional)
1A marked `mkn_rates_for_state`/`trans_rate_factor` `[[maybe_unused]]`; can delete.

## GATE 2 (after wiring)
Full unit suite + T5 (partial==full reconciliation after an accepted phi move and after
NNI) + T6 (prior-invariance under flat likelihood: φ-leg only at z=0 since φ now enters
the kernel; gibbs_z-leg at φ=1; pi0/theta legs prior-only, invariant). Retired tests gone.

**neo-z test must be RED-before-GREEN (advisor 2026-05-29).** The pi0/z-prior check only
verifies the fix if the test matrix CONTAINS transformational characters AND asserts the
z-prior / `pi0` posterior (or `cpp_log_prior` value at fixed state) is INVARIANT to the
trans-char z rows. A neo-only matrix passes before AND after → proves nothing. Concretely:
build a mixed matrix, set some trans z rows to ≠0, and assert `cpp_log_prior` (and R
`LogPrior`) is unchanged vs the all-trans-z=0 case. **Confirm this assertion FAILS on
current HEAD (trans rows counted) and PASSES after the mask.** Also a gibbs_z test: after a
sweep, trans/known z rows remain 0.

**Pre-recompile hygiene (advisor):** verify no live local R/Rscript holds the DLL before
`compile_dll()` (run via Rscript subprocess, never an active RStudio session). Sequence the
`--full` Hamilton submit (remote, no local lock) and any diagnostician pSPR-wiring AFTER the
local Phase 2b recompile to avoid DLL contention.

## NOTE — independent of this wiring
The root-inference move verification (math-prover + mcmc-diagnostician, dev/red-team/)
may return a verdict that a topology move's rooted Hastings ratio is biased → that's a
MOVE-code fix (tree_moves.cpp/proposals.cpp/dispatch), separate from this model wiring.
Reconcile those verdicts before GATE 2.
