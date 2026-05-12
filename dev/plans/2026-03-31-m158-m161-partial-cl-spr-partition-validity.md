# M-158 + M-161: Partial CL for SPR + Per-partition Cache Validity

## Overview

Two related optimizations to the node CL cache system:

- **M-158**: Enable partial CL evaluation for SPR moves (currently SPR always
  invalidates the entire cache and uses full O(nEdge) evaluation)
- **M-161**: Track per-CacheUnit validity so parameter-only changes (e.g.
  `rate_loss`) only invalidate affected units, preserving cache for unaffected
  partitions

## Current State

### Node CL Cache Architecture
- `NodeCLCache` in `node_cl_cache.h` holds per-node CLs in `CacheUnit` objects
- Each `CacheUnit` represents a partition (or kPrime-subgroup within a partition)
- `TreeNav` maintains topology as node-indexed parent/child/edge structure
- **Single `bool valid` flag** controls all-or-nothing cache validity
- Partial CL works for: NNI (moveType 5), beta_simplex (4), Dirichlet (23, 24)
- SPR (6), TBR (17), and all other topology/parameter moves invalidate fully

### SPR Flow (Case 6)
1. `spr_proposal_impl()` generates fully reordered `proposedParent/Child/RelBr`
2. Full likelihood evaluated on proposed topology
3. On acceptance: topology committed via `std::move()`, cache invalidated
4. On rejection: proposed vectors discarded, no rollback needed

### Key Constraint
`spr_proposal_impl()` calls `preorder_weighted_impl()` which fully reorders edge
arrays. After reordering, edge indices no longer correspond to TreeNav's
`edgeToPar[]` mapping. This is the central obstacle for partial CL with SPR.

---

## M-158: Partial CL Evaluation for SPR Moves

### Design: TreeNav-native SPR Proposal

Instead of trying to map `spr_proposal_impl()`'s reordered output back to the
TreeNav, implement the SPR proposal directly on the TreeNav structure.

#### New Functions in `node_cl_cache.h`

**1. `propose_spr_treenav()`**
Performs an SPR move proposal using only the TreeNav structure:
- Select random prune edge (node v with parent u, where u ≠ root)
- Identify sibling w of v under u, and p = parent(u)
- BFS on TreeNav to enumerate v's descendants (for exclusion)
- Select random eligible regraft edge (node b with parent a, not in v's subtree,
  not adjacent to u)
- Draw τ ~ Uniform(0,1)
- Compute Hastings ratio: `log(edgeLen[a→b]) - log(edgeLen[p→u] + edgeLen[u→w])`
- Return struct `SprProposal { u, v, p, w, a, b, tau, logHastings, nCandidates }`

**2. `update_topo_spr()`**
Apply SPR topology change to TreeNav:
```
Before: p→u→{v,w},  a→{b,...}
After:  p→{w,...},   a→u→{v,b}
```
Steps:
- Remove u from p's children; add w as p's child
- Remove b from a's children; add u as a's child
- Remove w from u's children; add b as u's child
- Update parentNode: `parentNode[w]=p`, `parentNode[u]=a`, `parentNode[b]=u`
- Update edge lengths:
  - `nodeEdgeLen[w] = old(nodeEdgeLen[u]) + old(nodeEdgeLen[w])` (merged)
  - `nodeEdgeLen[u] = τ × old(nodeEdgeLen[b])` (split: a→u portion)
  - `nodeEdgeLen[b] = (1-τ) × old(nodeEdgeLen[b])` (split: u→b portion)
- Save old values for reversal

**3. `reverse_topo_spr()`**
Undo TreeNav update using saved metadata. Exact reverse of `update_topo_spr()`.

**4. `find_dirty_spr()`**
Identify dirty set after SPR:
- Union of path from p→root and path from a→root, plus node u
- All three have changed children/subtree composition
- Sort in postorder (deepest first)
- Typical size: ~15–25 nodes for 54-tip tree (vs 8 for NNI)

**5. `treenav_to_preorder()`**
Reconstruct canonical `parent/child/relBrLengths` arrays from TreeNav.
Used on acceptance to commit topology to MCMC state.
- DFS traversal from root, emit edges in preorder
- O(nEdge) — only called on acceptance (~5–15% of SPR proposals)

#### New Field: `nodeEdgeLen[]` in TreeNav

Add `std::vector<double> nodeEdgeLen` to TreeNav: absolute edge length from
node to parent, indexed by node (not by edge array position).

- Built in `TreeNav::build()`: `nodeEdgeLen[c] = absEdgeLen[e]` for each edge
- Used by `recompute_dirty_nodes()` instead of `absEdgeLen[edgeToPar[ch]]`
- This decouples CL computation from edge array ordering
- Backward-compatible: existing NNI/Dirichlet partial CL still works, just
  using `nodeEdgeLen` instead of `edgeLen[edgeToPar[]]`

#### Modified: `recompute_dirty_nodes()`

Change edge length access from:
```cpp
double edgeLen = absEdgeLen[cache.topo.edgeToPar[ch]];
```
to:
```cpp
double edgeLen = cache.topo.nodeEdgeLen[ch];
```

This is the key change that enables TreeNav-native proposals without relying on
the edge array ordering. All existing partial CL moves (NNI, beta_simplex,
Dirichlet) must also update `nodeEdgeLen` when they change edge lengths.

For NNI: edge lengths don't change (only topology), so no `nodeEdgeLen` update
needed.

For beta_simplex and Dirichlet: they modify `state->relBrLengths[idx]`. Need
to also update `cache.topo.nodeEdgeLen[child[idx]]` at the same time.

#### Modified: `do_move_impl()` Case 6

```
if (SPR && cache.valid && !qHeterogeneity):
  1. Propose SPR on TreeNav → get {u,v,p,w,a,b,τ,logHR}
  2. Save old TreeNav edge lengths + topology at affected nodes
  3. update_topo_spr(treeNav, proposal)
  4. find_dirty_spr() → dirty set
  5. If dirty > 80% of tree: fallback to full eval (rare)
  6. save_dirty_cls()
  7. recompute_dirty_nodes() using nodeEdgeLen
  8. cache_total_loglik() → newLogLik
  9. usedPartialCL = true
  10. MH accept/reject:
      Accept: treenav_to_preorder() → commit parent/child/relBr
      Reject: reverse_topo_spr(), restore_dirty_cls()
else:
  Existing full-eval path (spr_proposal_impl + cpp_log_likelihood)
```

#### Diagnostics

Add `diagSprPartialCount`, `diagSprMismatchCount` to McmcState for
partial-vs-full comparison (same pattern as NNI/beta_simplex diagnostics).
Initially run both paths and compare; remove full-eval check once validated.

#### Edge Length Sync

When `nodeEdgeLen` is introduced, all places that modify edge lengths must
keep it in sync:

| Move | What changes | nodeEdgeLen update |
|------|-------------|-------------------|
| tree_length (0) | All edges scale | Invalidate cache (already does) |
| beta_simplex (4) | 2 edges | Update 2 entries in nodeEdgeLen |
| NNI (5) | Topology only | No edge len change (already handled) |
| SPR (6, partial) | 3 edges + topology | Handled by update_topo_spr |
| Dirichlet (23,24) | K edges | Update K entries in nodeEdgeLen |
| full_downpass | Rebuilds from scratch | nodeEdgeLen rebuilt in populate_cache_full |

### Expected Performance

For 54-tip tree (nEdge=105, depth≈7):
- SPR dirty set: ~15–25 nodes (two root paths + u)
- vs full evaluation: 105 edges
- **~5–7× speedup per SPR proposal evaluation**
- Plus: cache stays valid after SPR acceptance → subsequent NNI/Dirichlet moves
  don't need full repopulation

### Cache-Boost Integration (M-159)

After M-158, SPR (moveType 6) becomes a partial-CL-eligible move. Update the
cache-boost logic in `run_mcmc_batch_cpp` to include moveType 6 in the
boosted set (currently: 4, 5, 23, 24).

---

## M-161: Per-Partition CL Cache Validity

### Design: Per-Unit Validity Flags

Replace `NodeCLCache::valid` (single bool) with per-CacheUnit flags.

#### New Fields

```cpp
struct NodeCLCache {
  // Replace: bool valid = false;
  // With:
  std::vector<bool> unitValid;  // per-CacheUnit validity

  // Helper methods:
  bool allValid() const;   // all units valid
  bool anyValid() const;   // at least one unit valid
  void invalidateAll();    // set all to false
  void invalidateUnits(const std::vector<int>& indices);
};
```

#### Selective Invalidation

| Parameter change | Units to invalidate |
|-----------------|-------------------|
| Topology (SPR, TBR, NNI-reorder) | All units |
| Branch lengths (tree_length, beta_simplex, Dirichlet) | All units |
| `rate_loss` | Only neomorphic units (isMkN == true) |
| `rateNeo` | Only neomorphic units |
| `rateLogSd` | All units (ACRV rates shared) |
| `betaScale` | All units (but Q-het disables partial CL anyway) |

#### Selective Repopulation

`populate_cache_full()` → `populate_cache()` that skips already-valid units:

```cpp
for (int ui = 0; ui < units.size(); ++ui) {
  if (unitValid[ui]) continue;  // already valid, skip
  full_downpass(units[ui], topo, rates, rateLoss);
  unitValid[ui] = true;
}
```

This saves O(nEdge × nTransChars) when only neomorphic parameters changed.

#### Partial Eval with Mixed Validity

When some units are valid and some are invalid:
- **Option A (simpler)**: Repopulate invalid units before partial eval, then
  proceed as normal. Cost: full downpass for invalid units only.
- **Option B (complex)**: Mixed-mode eval — partial CL for valid units, full
  eval (via partition likelihood) for invalid units. More complex, more savings.

**Choose Option A** for simplicity. The repopulation cost for just the
neomorphic units is small (typically 10–15 characters vs 200+ transformational).

#### `anyValid()` for Scheduling

The cache-boost scheduler (M-159) needs to know if partial CL is possible:
- Currently checks `state->nodeCL.valid` (all-or-nothing)
- Change to `state->nodeCL.anyValid()` — partial CL possible if any unit is
  valid (will repopulate the rest on demand)

#### Modified Invalidation Points

All current `state->nodeCL.valid = false` become:
- For topology/branch changes: `state->nodeCL.invalidateAll()`
- For `rate_loss`/`rateNeo` only:
  `state->nodeCL.invalidateNeomorphic()` (new helper that only invalidates
  units where `isMkN == true`)
- For `rateLogSd`: `state->nodeCL.invalidateAll()` (ACRV affects all)

### Expected Benefit

Moderate, dataset-dependent:
- **Mixed datasets** (e.g., 12 neo + 15 trans): After `rate_loss` or `rateNeo`
  slice move, transformational cache preserved → next NNI only repopulates
  neomorphic units (12 chars) instead of all (27 chars). ~55% savings on
  repopulation.
- **All-transformational** (Sun2018): No benefit (no neomorphic partitions to
  selectively invalidate).
- **Scaling**: Benefit grows with partition diversity and frequency of
  neomorphic-only parameter moves.

---

## Implementation Order

1. **M-158 Phase 1**: Add `nodeEdgeLen[]` to TreeNav, update
   `recompute_dirty_nodes()` and beta_simplex/Dirichlet to use it.
   Run existing tests to verify backward compatibility.

2. **M-158 Phase 2**: Implement `propose_spr_treenav()`,
   `update_topo_spr()`, `reverse_topo_spr()`, `find_dirty_spr()`.
   Unit tests with manual SPR + verify dirty set.

3. **M-158 Phase 3**: Implement `treenav_to_preorder()`.
   Test: build TreeNav → SPR on TreeNav → reconstruct → compare with
   `spr_proposal_impl()` output.

4. **M-158 Phase 4**: Wire into `do_move_impl()` case 6 with diagnostic
   full-eval comparison. Run full test suite.

5. **M-158 Phase 5**: Add SPR to cache-boost set (M-159 integration).
   Remove diagnostic full-eval comparison after validation.

6. **M-161**: Replace `bool valid` with per-unit flags. Update all
   invalidation points. Modify `populate_cache_full()` for selective
   repopulation. Tests with mixed neo+trans datasets.

## Testing Strategy

### M-158 Tests
1. `propose_spr_treenav()` produces valid SPR moves (check node identities,
   Hastings ratio matches `spr_proposal_impl()`)
2. `update_topo_spr()` + `reverse_topo_spr()` roundtrip = identity
3. `find_dirty_spr()` contains expected nodes (prune parent, regraft parent,
   path to root)
4. `treenav_to_preorder()` output matches `preorder_weighted_impl()` applied
   to same topology
5. Partial CL loglik matches full eval for SPR (diagnostic comparison, like
   existing NNI/Dirichlet diagnostics)
6. Full MCMC run with SPR partial CL produces valid posteriors
7. Cache stays valid after SPR acceptance → verified by checking
   `state->nodeCL.valid` and immediate NNI partial eval

### M-161 Tests
1. After `rate_loss` change: neomorphic units invalidated, transformational
   units remain valid
2. After topology change: all units invalidated
3. Selective repopulation: only invalid units recomputed (verify via
   `diagCachePopCount` or timing)
4. Mixed dataset MCMC produces same posteriors as before

## Risk Assessment

**M-158** is high engineering effort with multiple interacting components. Key
risks:
- TreeNav SPR proposal must produce identical move distributions to
  `spr_proposal_impl()` (same eligible set, same Hastings ratio)
- `treenav_to_preorder()` must produce valid canonical preorder
- Edge length sync between `nodeEdgeLen` and `relBrLengths` must be watertight

**M-161** is moderate effort, lower risk. Main concern is ensuring all
invalidation paths are covered (same bug class as M-143/M-145).
