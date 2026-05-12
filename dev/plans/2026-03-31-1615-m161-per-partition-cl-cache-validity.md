# M-161: Per-Partition CL Cache Validity

**Task:** Avoid invalidating ALL CacheUnits when only a subset of partitions
is affected by a parameter change. Currently `nodeCL.valid = false`
invalidates everything, triggering a full `populate_cache_full()` on the next
partial-CL-eligible move.

**Primary benefit:** After `rate_loss` or `rate_neo` acceptance, only
neomorphic CacheUnits need recomputation. Transformational/known units'
CLs are preserved, skipping their expensive full downpass.

---

## Current Behaviour

Single `NodeCLCache.valid` flag. Any parameter change that affects likelihood
sets `valid = false`. Before the next NNI/beta_simplex/Dirichlet/SPR, the
entire cache is rebuilt from scratch: TreeNav, unit structure, allocation,
tip init, full downpass on ALL units.

## Proposed Design: Two-Level Validity

### Level 1: Global Validity

| Flag | Meaning | Triggers false |
|------|---------|---------------|
| `topoValid` | TreeNav (topology + nodeEdgeLen) matches current state | tree_length accepted; non-partial-CL topology moves (Gibbs SPR/swap, pSPR, TBR, weighted moves) |
| `structureValid` | Unit structure (partition→unit mapping, kStates per unit) matches current kPrime | kPrime int_walk accepted; Gibbs kPrime sweep; block kPrime shift |

When either is false → **full rebuild** (same as current `populate_cache_full`).

### Level 2: Per-Unit CL Validity

| Field | Location | Meaning |
|-------|----------|---------|
| `CacheUnit::clValid` | per-unit | Internal-node CLs in this unit are current |

When `topoValid` and `structureValid` are both true but some units have
`clValid = false` → **selective repopulation**: only recompute those units'
downpass. TreeNav, allocation, tip CLs all preserved.

### Ready Check

```cpp
bool NodeCLCache::ready() const {
  if (!topoValid || !structureValid) return false;
  for (const auto& u : units)
    if (!u.clValid) return false;
  return true;
}
```

This replaces all reads of the old `valid` flag.

---

## Invalidation Map

### Parameter → what's affected

| Event | topoValid | structureValid | Per-unit clValid |
|-------|-----------|----------------|------------------|
| tree_length accepted | ❌ | — | — (all stale anyway) |
| rate_loss accepted | ✓ | ✓ | neo units → ❌ |
| rateLogSd accepted | ✓ | ✓ | ALL units → ❌ |
| rate_neo accepted | ✓ | ✓ | neo units → ❌ |
| beta_scale accepted | ✓ | ✓ | ALL units → ❌ (Q-het; partial CL disabled anyway) |
| kPrime int_walk | ✓ | ❌ | — (structure stale) |
| kPrime Gibbs/block | ✓ | ❌ | — (structure stale) |
| neo_joint_scale (case 18) | ✓ | ✓ | neo units → ❌ |
| Non-partial topology moves | ❌ | — | — (topo stale) |
| Partial-CL moves (NNI, BS, Dir, SPR) | ✓ (updated in-place) | ✓ | ✓ (dirty nodes updated in-place) |

### Source-level invalidation sites (13 in mcmc.cpp)

| Line(s) | Function | Currently | New |
|----------|----------|-----------|-----|
| 893, 1202, 1346, 1601, 1857, 1966, 2532, 2738 | Gibbs/weighted/pSPR topology moves | `valid = false` | `topoValid = false` |
| 2250 | weighted_branch_scale | `valid = false` | `topoValid = false` |
| 2893 | slice_scalar_impl | `valid = false` | paramIdx-dependent (see below) |
| 3550 | gibbs_kprime_sweep | `valid = false` | `structureValid = false` |
| 3645 | block_kprime_shift | `valid = false` | `structureValid = false` |
| 4393 | generic MH acceptance | `valid = false` | moveType-dependent (see below) |

### slice_scalar_impl (line 2893) dispatch by paramIdx

```cpp
switch (paramIdx) {
  case 1: case 3:  // rate_loss, rate_neo
    invalidate_neo_cls(state->nodeCL);
    break;
  case 2:  // rateLogSd — ACRV rates change, all units stale
    invalidate_all_cls(state->nodeCL);
    break;
  default:  // tree_length, beta_scale, etc.
    state->nodeCL.topoValid = false;
    break;
}
```

### Generic MH acceptance (line 4393) dispatch by moveType

```cpp
if (!usedPartialCL && likChanges) {
  switch (moveType) {
    case 1: case 3: case 18:  // rate_loss, rate_neo, neo_joint
      invalidate_neo_cls(state->nodeCL);
      break;
    case 2:  // rateLogSd
      invalidate_all_cls(state->nodeCL);
      break;
    case 7:  // kPrime int_walk
      state->nodeCL.structureValid = false;
      break;
    case 16:  // beta_scale (Q-het; partial CL disabled)
      invalidate_all_cls(state->nodeCL);
      break;
    default:  // tree_length (0), topology, etc.
      state->nodeCL.topoValid = false;
      break;
  }
}
```

---

## Selective Repopulation

Replace `populate_cache_full()` with `populate_cache()`:

```
populate_cache(cache, data, parent, child, absEdgeLen, kPrime,
               rateLoss, rateLogSd, rateNeo):

  if (!topoValid || !structureValid):
    → full rebuild (same as current populate_cache_full)
    → set topoValid = structureValid = true, all clValid = true
    return

  # Topo + structure valid, but some units need CL recomputation
  # Check if ACRV rates are stale
  if (cachedRateLogSd != rateLogSd):
    recompute ACRV rates → cache.rates
    cache.cachedRateLogSd = rateLogSd

  for each unit in cache.units:
    if unit.clValid: continue

    # Update unit-specific parameters that may have changed
    if unit.isMkN:
      unit.rootFreqs = { rl/(1+rl), 1/(1+rl) }
      unit.rateScale = rateNeo

    # Rerun downpass (tip CLs are still valid, only internal nodes recomputed)
    full_downpass(unit, cache.topo, cache.rates, rateLoss)
    unit.clValid = true
```

**Key insight:** `full_downpass()` doesn't modify tip CLs (they're marked as
already initialized via `initFlg`). Tips only depend on character data, which
never changes. So we skip allocation and tip init for valid units entirely.

---

## Call-site Changes

### Pre-proposal cache population (mcmc.cpp ~line 3758)

Replace:
```cpp
if (...eligible... && !state->nodeCL.valid)
  populate_cache_full(...);
```

With:
```cpp
if (...eligible... && !state->nodeCL.ready())
  populate_cache(...);
```

### partial_eval_dirty (node_cl_cache.h ~line 4162 etc.)

Replace `state->nodeCL.valid` reads with `state->nodeCL.ready()`.
There are ~4 reads: lines 4162, 4190, 4222, 4243 (NNI, BS, Dir, SPR paths).

---

## Data Structure Changes

### CacheUnit (add 1 field)

```cpp
struct CacheUnit {
  // ... existing fields ...
  bool clValid = false;  // M-161: per-unit CL validity
};
```

### NodeCLCache (replace valid, add helpers)

```cpp
struct NodeCLCache {
  // ... existing fields ...

  // M-161: two-level validity
  bool topoValid = false;       // replaces old `valid`
  bool structureValid = false;  // unit structure matches current kPrime

  bool ready() const {
    if (!topoValid || !structureValid) return false;
    for (const auto& u : units)
      if (!u.clValid) return false;
    return true;
  }

  void invalidate_all() {
    topoValid = false;
  }

  void invalidate_structure() {
    structureValid = false;
  }

  void invalidate_neo_cls() {
    for (auto& u : units)
      if (u.isMkN) u.clValid = false;
  }

  void invalidate_all_cls() {
    for (auto& u : units)
      u.clValid = false;
  }
};
```

Remove old `bool valid = false;`.

---

## Expected Performance Impact

### Mixed dataset (e.g., 12 neo + 15 trans chars)
- After rate_loss/rate_neo acceptance: skip ~55% of downpass work (trans units)
- After rateLogSd acceptance: no improvement (all units affected)
- After tree_length: no improvement (full rebuild)

### All-transformational (Sun2018: 225 chars)
- rate_loss/rate_neo don't exist → no benefit from neo selectivity
- kPrime changes: still trigger full rebuild (structureValid = false)
- Marginal benefit from skipping TreeNav + allocation on rateLogSd changes

### Benefit scales with partition diversity
Most impactful for datasets with both neomorphic and transformational
characters. The 4789/3408/1271 validation datasets fall in this category.

---

## Testing Plan

1. **Unit test: selective repopulate correctness.**
   Build state with mixed neo + trans data. Run a few NNI moves (cache
   becomes valid). Accept a rate_loss change (neo units invalidated).
   Run `populate_cache()` (selective path). Compare `cache_total_loglik()`
   against `cpp_log_likelihood()` — must match to <1e-10.

2. **Unit test: ACRV selective repopulate.**
   Same as above but change rateLogSd (all units invalidated, but topo
   + structure valid). Verify selective repopulate matches full eval.

3. **Regression: existing test suites.**
   Run test-gibbs.R, test-gibbs-spr.R, test-gibbs-het.R,
   test-gibbs-kprime.R, test-gibbs-batched.R, test-likelihood.R,
   test-fused-ascertainment.R, test-m092-adaptive-scheduler.R.
   Zero regressions expected.

4. **Diagnostic counter.**
   Add `diagSelectivePopCount` to track how often the selective path
   fires (vs full rebuild). Report in the existing diagnostic output.

---

## Risk Assessment

**Low risk.** The selective path is purely additive — when either global flag
is false, we fall back to the existing full rebuild code. The selective path
only fires when topology and unit structure are known-valid, and only
recomputes units that need it.

The main correctness risk is forgetting to invalidate a unit when its
parameters change. Mitigated by:
- The diagnostic comparison (partial-CL vs full-eval) that already runs
  on every NNI and beta_simplex
- The regression test suite (~4000 tests)
- Conservative fallback: any unhandled moveType goes to `topoValid = false`

---

## Implementation Steps

1. Add `clValid` to `CacheUnit`, add `topoValid`/`structureValid`/helpers
   to `NodeCLCache`. Remove old `valid` field.
2. Add `populate_cache()` with selective path in `node_cl_cache.h`.
3. Update all 13 invalidation sites in `mcmc.cpp` to use granular methods.
4. Update 4 `cache.valid` reads to `cache.ready()` in mcmc.cpp.
5. Update `populate_cache_full()` to set new flags correctly.
6. Add diagnostic counter for selective repopulation.
7. Write new test (test-m161-selective-cache.R).
8. Run full test suite.
