# M-160: Delayed Rejection for SPR Topology Moves

## Summary

When an SPR proposal is rejected, attempt a cheaper NNI as a fallback
using the Delayed Rejection (DR) framework (Green & Mira 2001). The
DR acceptance ratio recycles the SPR likelihood (already computed) and
evaluates the NNI via the node CL cache (O(depth) cost). This gives
a second chance at tree movement on iterations where SPR fails, at
negligible additional cost.

## Background: DR acceptance ratio

Standard two-stage DR with symmetric, independent second proposal:

```
α₂(x, y, z) = min(1, exp(logα₂))

logα₂ = β·(logLik_z - logLik_x)
       + log(1 - α₁(z,y)) - log(1 - α₁(x,y))
```

Where:
- x = current tree
- y = rejected SPR proposal
- z = NNI proposal from x
- α₁(x,y) = min(1, exp(β·(logLik_y - logLik_x) + logHR_spr))
  — the standard SPR MH ratio (just rejected)
- α₁(z,y) = min(1, exp(β·(logLik_y - logLik_z) + logHR_spr))
  — the hypothetical SPR MH ratio from z to y

### Why this simplifies

1. **NNI is symmetric** (logHR_nni = 0) → q_NNI terms cancel.
2. **Prior is constant** for topology moves (flat Dirichlet on branch
   proportions, scalars unchanged) → prior terms cancel.
3. **Non-overlapping NNI** (NNI nodes disjoint from SPR prune/regraft
   nodes) → q_SPR(z,y) = q_SPR(x,y) and logHR_spr(z→y) = logHR_spr(x→y)
   → proposal density ratio = 1.
4. **logLik_y is recycled** from stage 1 (no recomputation).
5. **logLik_z uses node CL cache** (O(depth) partial eval).

### Numerical stability for log(1-α)

Use `log1mexp`: for log(1 - exp(a)) where a < 0:
- if a < -log(2): `log1p(-exp(a))`
- if a >= -log(2): `log(-expm1(a))`

When a >= 0 (α₁ = 1): log(1-α) = -∞ → DR rejects. This is correct:
if the SPR from z would be accepted to reach y, then z→x is the wrong
direction.

## Eligibility conditions

DR-NNI is attempted only when ALL of:
1. moveType == 6 (SPR)
2. SPR was evaluated and rejected (logAlpha was finite)
3. `state->nodeCL.valid` (node CL cache available for cheap NNI eval)
4. `!data->qHeterogeneity` (partial CL not supported for Q-het)
5. NNI proposal succeeds (valid internal edge, in-place safe)
6. **Non-overlap**: none of the 4 NNI nodes (u_nni, v_nni, c_nni, w_nni)
   is any of the SPR key nodes (u, v, a, b, parent_of_u)

When any condition fails, fall through to `return false` (standard SPR
rejection).

## Implementation plan

### Step 1: Extend `spr_proposal_impl` to return metadata

Add fields to the returned List:
- `"pruneParent"` (int u), `"pruneChild"` (int v)
- `"regraftParent"` (int a), `"regraftChild"` (int b)
- `"parentOfPruneParent"` (int p = parent of u — the node above prune point)

This is a few extra lines; the values are already available as local
variables in `spr_proposal_impl`.

### Step 2: Save SPR metadata in case 6

In the moveType 6 block, after extracting the proposal, save the
metadata into local variables (`sprU`, `sprV`, `sprA`, `sprB`, `sprP`).
These are only used if the SPR is later rejected.

### Step 3: DR-NNI block after rejection

Insert between the rollback block and `return false;` in `do_move_impl`:

```
if (moveType == 6 && drEligible) {
    // 1. Propose random in-place NNI (reuse case-5 logic)
    // 2. Check non-overlap against SPR nodes
    // 3. Update node CL cache for NNI
    // 4. Compute logLik_z via partial_eval_dirty
    // 5. Compute α₁(z,y) using recycled logLik_y
    // 6. Compute logα₂
    // 7. Accept → commit NNI topology + invalidate cache; return true
    // 8. Reject → rollback NNI + restore cache; return false
}
```

The NNI proposal code is extracted from case 5 into a helper to avoid
duplication. Only the in-place path is used (the non-in-place reorder
path is not cache-friendly and would be too expensive for a fallback).

### Step 4: Helper function for log(1-α)

Add `log1mexp(double a)` utility (stable computation of log(1-exp(a))
for a < 0). Used by the DR ratio.

### Step 5: Diagnostic counters

Add to McmcState:
- `diagDrAttempts`: number of DR-NNI attempts (SPR rejected + eligible)
- `diagDrAccepts`: number of DR-NNI acceptances
- `diagDrOverlap`: number skipped due to node overlap

Expose via `get_mcmc_state_cpp`.

### Step 6: R-side changes

Minimal:
- No new moveType needed (DR is part of SPR internally)
- No changes to move weights or adaptation
- DR accepts count as SPR accepts in the batch runner's acceptance
  matrix (the batch runner sees `do_move_impl` return true)
- Expose new diagnostic counters in `.PrintDiagnostics()` or similar

### Step 7: Testing

1. **Unit test: DR ratio arithmetic** — construct known logLik values,
   verify the computed logα₂ matches hand calculation.
2. **Unit test: non-overlap check** — construct SPR+NNI node sets with
   and without overlap, verify correct detection.
3. **Integration test** — run short MCMC with DR enabled, verify
   diagDrAttempts > 0 and diagDrAccepts > 0.
4. **Reversibility check** — on a small tree (5-8 tips), run many DR
   iterations and verify empirical detailed balance.
5. **Regression** — existing SPR-only tests must pass unchanged.

## Expected benefit

- SPR acceptance rate is typically 5-15%. DR gives ~80% of rejected SPR
  iterations a second chance via NNI (20% excluded by overlap).
- NNI evaluation via node CL cache costs ~5% of a full SPR evaluation.
- Net: significant increase in topology-changing moves per unit time,
  especially for cold chains where SPR acceptance is low.
- No change to the stationary distribution (DR preserves detailed
  balance).

## Files modified

| File | Changes |
|------|---------|
| `src/proposals.cpp` | Extend `spr_proposal_impl` return value |
| `src/mcmc.cpp` | DR-NNI block in `do_move_impl`; NNI helper extraction; `log1mexp` utility; diagnostic counters |
| `src/mcmc_state.h` | (if McmcState is there) DR diagnostic fields |
| `tests/testthat/test-delayed-rejection.R` | New test file |
