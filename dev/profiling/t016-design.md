# T-016 Design — Orchestrator heap-alloc hoisting (mirror T-008)

**Status:** PENDING (2026-05-22). Self-contained brief for a chip / subagent
to enact post-compact. Reads cold from this doc + cross-references.

**Branch base:** `worktree-ecology-aware` (already has T-014/T-015/T-017/T-018).

---

## Headline

The aware ecology orchestrators (`cpp_partition_log_likelihood_ecology`,
`const_site_prob_jc_eco_single`, `const_site_prob_mkn_eco_single`) still
allocate ~800 KB of `std::vector` workspace per call uniformly per chain.
Mirror the M-162B / T-008 ClWorkspace pattern: pre-allocate the workspaces
in `McmcState`, thread them through, eliminate per-call malloc/free.

**Estimated impact:** 5-10 % wall reduction on aware nChains=1
(uniform per-chain cost, independent of nChains). Smaller on nChains=4
than Phase 1b would have given. **Not enough on its own** to make Hamilton
viable, but combines additively with future work.

**Estimated LoC:** ~100 net. Bounded refactor, established pattern.

---

## Why this is bounded and low-risk

This is exactly the T-008 work, ported to a different set of call sites.
T-008 successfully refactored `per_char_log_lik_ecology` (the Gibbs z
sweep per-character helper) the same way, eliminating ~7800 heap
allocations per Gibbs z sweep. See T-008 row in `findings.md` for the
GibbsZWorkspace pattern.

Phase 1b of T-017 already added `EcoThreadScratch` for the partition-loop
workers; the same struct can be extended (or a sibling struct created)
for the orchestrator entry-point allocations.

---

## Call sites to refactor

In `src/mcmc_ecology.cpp`:

### `cpp_partition_log_likelihood_ecology` (~:967-1112)

Per-call allocations:
- Line ~991: `IntegerMatrix zPart(nCharPart, zCols)` — Rcpp alloc
- Line ~1001-1006 (neo branch): `NumericVector neoEl(nEdge)`,
  `NumericVector rootFreqs(2)`, `std::vector<double> buf((maxNode+1)*stride, 0.0)`,
  `std::vector<uint8_t> initFlg(maxNode+1, 0u)`
- Line ~1015-1018 (neo ascertainment): `IntegerVector zVec(zCols)` per char in inner loop (!)
- Line ~1027-1030 (known branch): `NumericVector rootFreqs(kStates, 1/kStates)`,
  `std::vector<double> buf(...)`, `std::vector<uint8_t> initFlg(...)`
- Line ~1038-1042 (known ascertainment): `IntegerVector zVec(zCols)` per char
- Line ~1062-1080 (trans branch, per kPrime subgroup): `IntegerMatrix subStates`,
  `IntegerMatrix subZ`, `NumericVector rootFreqs`, `std::vector<double> buf`,
  `std::vector<uint8_t> initFlg`, and inner per-char `IntegerVector zVec`

The trans branch is called per `kp` subgroup — for rodent, 3 subgroups
per call → 3× the allocs from the trans branch alone.

### `const_site_prob_jc_eco_single` (:876)

- `IntegerMatrix tipStates(nTip, 1)` (all zeros for ascertainment)
- `IntegerMatrix zMat(1, zVec.size())`
- `NumericVector rootFreqs(kStates, 1/kStates)`
- `std::vector<double> buf((maxNode+1)*stride, 0.0)`
- `std::vector<uint8_t> initFlg(maxNode+1, 0u)`

### `const_site_prob_mkn_eco_single` (:907)

- Same as above plus `IntegerMatrix tipStates1` for the second
  pseudo-character pass

---

## Refactor sketch

Add an `EcoOrchestratorScratch` struct similar to `EcoThreadScratch`
but covering the orchestrator-entry-point allocations:

```cpp
struct EcoOrchestratorScratch {
  std::vector<int>     zPart;        // (maxNCharPart * zCols)
  std::vector<double>  neoEl;        // nEdge
  std::vector<double>  rootFreqs;    // maxKStates
  std::vector<double>  buf;          // (maxNode+1) * maxStride
  std::vector<uint8_t> initFlg;      // maxNode+1
  std::vector<int>     subStates;    // nTip * maxNSub (trans branch)
  std::vector<int>     subZ;         // maxNSub * zCols
  std::vector<int>     zVec;         // zCols (ascertainment per-char)
  std::vector<int>     tipStatesZero; // nTip * 1 (for const_site_prob)
  std::vector<int>     tipStatesOne;  // nTip * 1 (for mkn second pass)

  void allocate(int maxNode, int maxStride, int maxNCharPart, int maxNSub,
                int zCols, int nTip, int nEdge, int maxKStates);
};
```

Owned where? Three options:

(a) **Function-static** (like Phase 1b's scratch): simplest, but means
serial dispatch only (which is fine for now — chain loop is serial).
RECOMMENDED for Phase 1.

(b) **Per-chain on McmcState**: more invasive (signature change cascades),
but needed for Phase 2 (threaded chain dispatch). Defer until then.

(c) **ClWorkspace extension**: re-use existing `clWs` on McmcState. May
have size constraints (ClWorkspace is sized for the blind path). Skip.

Pick (a) for this round.

---

## Concrete steps

1. Branch off `worktree-ecology-aware` as `t016-orchestrator-scratch`.
2. Add `EcoOrchestratorScratch` struct in `src/mcmc_ecology.cpp` near the
   existing `EcoThreadScratch`.
3. Refactor `cpp_partition_log_likelihood_ecology` to use a function-static
   `EcoOrchestratorScratch` (grow-only-resize, same pattern as Phase 1b).
4. Refactor `const_site_prob_jc_eco_single` and `const_site_prob_mkn_eco_single`
   similarly (their own static scratch, or a shared one).
5. Verify 246/246 tests pass (bit-identity preserved).
6. Bench: rodent_aware_timing.R MKP_TIMING_NCHAINS=1 + nChains=4. Compare
   to current baseline (~46 s nChains=1, ~174 s nChains=4 with OMP=1).
   Expected: 5-10% improvement.
7. If bench shows ≥5% improvement, file as APPLIED in findings.md.
   Otherwise file as REFUTED with measured numbers.
8. Commit, push, stop.

---

## Bench expectation (honest)

T-008 itself only delivered ~4% on the per-char batch (because pruner
arithmetic dominated malloc). The orchestrator allocs are bigger
absolute (~800 KB per call vs ~4 KB for per_char), so the projected
gain is bigger — but T-010's partition cache already cuts the
orchestrator call frequency by ~3× (only tree-moves now invalidate
the cache). So the practical impact is bounded by:

- ~22% of moves use the partition cache (no orchestrator call) — no impact.
- ~78% of moves are tree-touching → orchestrator called fresh → uniform
  saving of ~5-10% per orchestrator call due to no-malloc.

Net rodent projection: 4-7% wall reduction. Don't oversell.

---

## What's out of scope

- Don't extend to the inner pruners (`pruning_*_flat_ecology`): they
  already use std::vector workspaces passed by pointer.
- Don't touch the partition cache logic (`state->partLogLik` machinery):
  this is purely about replacing allocator calls.
- Don't change any RNG calls or move handlers.
- Don't add new tests beyond ensuring existing 246 pass.

---

## File references

- `src/mcmc_ecology.cpp:967-1112` — `cpp_partition_log_likelihood_ecology`
- `src/mcmc_ecology.cpp:876` — `const_site_prob_jc_eco_single`
- `src/mcmc_ecology.cpp:907` — `const_site_prob_mkn_eco_single`
- `src/mcmc_ecology.cpp:~1290` (post-Phase-1b) — existing `EcoThreadScratch`
  pattern to mirror
- `src/gibbs_z_workspace.h` — T-008's `GibbsZWorkspace` pattern (for the
  per_char_log_lik_ecology case)
