# M-121: Node-level CL dirty flags for standard NNI/SPR

## Problem

Standard NNI and SPR proposals do a full O(nEdge) Felsenstein downpass
for every proposal, even though most moves only invalidate CLs at a
small number of nodes. NNI swaps two children between adjacent nodes —
only those two nodes and their ancestors to the root need recomputation
(O(depth)). beta_simplex changes two branch lengths — same story.

The Gibbs partial CL framework (M-105/M-111) already solves this for
Gibbs moves via `caching_downpass` + `evaluate_candidate`, but that
design is per-invocation (cache built fresh each call, used for ~100
candidates, then discarded). Standard MH moves need a **persistent**
cache that survives across iterations, with selective dirty-node updates
and rollback on rejection.

## Scope

### In scope (this task)
- Persistent per-partition CL cache in `McmcState`
- Partial CL evaluation for **NNI** moves (in-place, node IDs stable)
- Partial CL evaluation for **beta_simplex** moves (2 edges change)
- Rollback of dirty CLs on rejection
- Full-evaluation fallback when cache is invalid
- ACRV support (per-category CL storage)

### Out of scope
- SPR/TBR partial CL (node reordering breaks cache; deferred)
- Q-heterogeneity (fall back to full eval; already handled by M-114 for Gibbs)
- Weighted/Gibbs/block moves (have their own CL management)

### Rationale for NNI + beta_simplex focus

NNI and beta_simplex together account for ~45% of moves in a typical
run. SPR/TBR change the canonical node numbering on acceptance (via
`preorder_weighted_impl`), invalidating the cache — so partial CL
for SPR would require either mapping old→new node IDs (O(nEdge), no
savings) or abandoning canonical reordering (breaks OPP-6b invariants).
NNI's in-place modification (OPP-6b) preserves node IDs, making it
the natural first target.

## Expected speedup

For a 54-tip tree (depth ≈ 7, nEdge = 105):
- Partial NNI/beta_simplex: O(depth) ≈ 7 nodes vs O(nEdge) ≈ 105 → ~15× per-move
- ~45% of moves benefit; cache valid ~82% of the time → ~37% of moves use partial
- Weighted average pruning speedup: ~1.5× overall
- Since pruning is ~90% of runtime: **~40% overall speedup**

## Design

### Data structures

```
struct NodeCLCache {
  // Per-partition flat CL buffers (persist across iterations)
  struct PartitionBuf {
    std::vector<double> cl;       // [nCat * (maxNode+1) * stride]
    int stride;                   // nChar * kStates for this partition
    int nCat;                     // 1 or nAcrvCat
  };
  std::vector<PartitionBuf> parts;

  // Tree navigation (persistent, updated incrementally for NNI)
  TreeNav topo;

  // Validity tracking
  bool valid = false;             // false → next eval must be full
  int maxNode = 0;

  // Rollback scratch (reused across iterations)
  std::vector<int> dirtyNodes;    // nodes to recompute (postorder)
  std::vector<double> savedCL;    // CL values before overwrite
};
```

Add `NodeCLCache nodeCL;` to `McmcState`.

### Cache lifecycle

| Event | Action |
|-------|--------|
| First full eval | Populate cache from pruning results; set `valid = true` |
| NNI proposed | Identify dirty nodes → save → partial eval → MH test |
| NNI accepted | Cache already updated; update `topo` in-place |
| NNI rejected | Restore saved CLs + parent values (existing O(1) rollback) |
| beta_simplex proposed | Same dirty-path pattern (2 edges → parents → root) |
| SPR/TBR accepted | Rebuild `topo` from new vectors; set `valid = false` |
| SPR/TBR rejected | Cache still valid (state unchanged) |
| Global param change (tree_length, rate_log_sd, beta_scale) | `valid = false` |
| rate_loss/rate_neo change | Invalidate MkN partitions only |
| kPrime change | Invalidate affected partition only |
| Gibbs/weighted moves | Set `valid = false` after accept (topology changed) |

### Partial evaluation algorithm (NNI)

```
partial_eval_nni(data, state, cache):
  // 1. Identify dirty set
  //    NNI swapped children of u and v. Dirty nodes:
  //    {v, u} ∪ ancestors(u → root)
  //    (v is child of u, so path is v → u → ... → root)
  dirty = collect_path(v, root)  // includes v, u, ancestors

  // 2. Save CLs at dirty nodes (all partitions, all categories)
  for each partition p:
    for each node n in dirty:
      for each cat c:
        copy cache.parts[p].cl[cat,node] → savedCL buffer

  // 3. Recompute dirty nodes in postorder
  for each node n in dirty (postorder = children before parents):
    for each partition p:
      for each cat c:
        for each child ch of n:
          apply transition(cache.cl[cat,ch], edgeLen[ch]) → contrib
        cache.cl[cat,n] = product of contributions

  // 4. Compute log-likelihood at root
  logLik = 0
  for each partition p:
    logLik += root_log_lik(cache.parts[p], root, rootFreqs, rates)
    // includes ACRV averaging and ascertainment correction

  return logLik
```

### Dirty node identification

**NNI:** After swapping `parent[cRow] = u` and `parent[wRow] = v`:
- Nodes v and u have modified child sets → both dirty
- All ancestors of u up to root are dirty (CL propagates)
- Since v is a child of u, the dirty set is just the path from v to root:
  `[v, u, g, ..., root]` — O(depth) nodes

**beta_simplex:** Two edges `idx1` and `idx2` have modified lengths:
- Parent nodes of those edges are dirty
- If both edges share a parent (common for adjacent edges): one path to root
- Otherwise: two paths that merge at their LCA, then continue to root
- Total: O(depth) nodes

### Integration into do_move_impl

```cpp
// In the likelihood evaluation section (line ~2810):

if (moveType == 5 && state->nodeCL.valid) {
  // NNI with valid cache → partial evaluation
  newLogLik = partial_eval_nni(data, state, v, u, ...);
  // newPC (partition cache) updated from node CL cache
} else if (moveType == 4 && state->nodeCL.valid) {
  // beta_simplex with valid cache → partial evaluation
  newLogLik = partial_eval_beta_simplex(data, state, bsIdx1, bsIdx2, ...);
} else {
  // Existing full evaluation path (unchanged)
  // After full eval, populate nodeCL cache from results
}

// On acceptance:
if (accepted) {
  // ... existing code ...
  if (moveType == 5) {
    // NNI: update topo in-place (swap children)
    state->nodeCL.topo.swapChildren(v, u, c, w);
  }
  if (topologyChanged) {
    // SPR/TBR: rebuild topo, invalidate cache
    state->nodeCL.topo.build(state->parent, state->child, data->nTip);
    state->nodeCL.valid = false;
  }
}

// On rejection:
if (rejected && partialEvalUsed) {
  // Restore saved CLs at dirty nodes
  restore_dirty_cls(state->nodeCL, savedCL, dirtyNodes);
}
```

### Populating cache from full evaluation

After a full evaluation (cache invalid or global param change), copy
the CL values from the pruning workspace into the persistent cache.
This requires the flat pruning functions to write into the cache's
partition buffers instead of (or in addition to) the ClWorkspace.

**Approach:** Pass each partition's cache buffer directly as the `buf`
pointer to `pruning_*_flat()`. The existing functions already accept
external buffer pointers — just point them at `cache.parts[pi].cl`
instead of `state->clWs.buf`. After the full eval, set `valid = true`.

This means the cache buffers need to be sized per-partition (each
with its own stride), unlike ClWorkspace which is a single shared
buffer sized to the maximum stride. The per-partition approach is
cleaner anyway since it avoids the stride-mismatch issue.

### ACRV handling

With ACRV, each partition stores `nCat` copies of CLs (one per rate
category). The cache layout is:
```
cl[(cat * (maxNode+1) + node) * stride + c * kStates + s]
```
Same as CLGroup. The partial eval iterates over categories when
recomputing dirty nodes, and averages across categories at the root.

### Ascertainment correction

The ascertainment correction (variable coding) requires computing the
probability of constant sites. This depends on the root CLs for
pseudo-characters (all-0, all-1, etc.), which must also be cached
and partially updated. Store a separate `constCL` buffer per partition
for constant-site pseudo-characters, with the same dirty-node update
logic.

## Implementation steps

### Step 1: NodeCLCache data structure
- Define `NodeCLCache` in `mcmc_state.h`
- Add `NodeCLCache nodeCL` member to `McmcState`
- Implement allocation function `allocate_node_cl_cache()`
- Call from `allocate_cl_workspace()` (or alongside it)

### Step 2: Populate cache from full evaluation
- Modify `cpp_partition_log_likelihood()` to optionally write CLs
  into an external per-partition buffer (NodeCLCache)
- After full eval in `do_move_impl`, mark `nodeCL.valid = true`
- Persistent TreeNav built from current topology

### Step 3: Dirty node identification helpers
- `find_dirty_nni(topo, v, u) → vector<int>` (path v → root)
- `find_dirty_beta_simplex(topo, edge1, edge2) → vector<int>`
- Both return nodes in postorder (children before parents)

### Step 4: Partial evaluation for NNI
- `partial_eval_nni()`: save dirty CLs, recompute path, return logLik
- Integrate into `do_move_impl` case 5
- Rollback: restore saved CLs on rejection
- Update TreeNav on acceptance

### Step 5: Partial evaluation for beta_simplex
- `partial_eval_beta_simplex()`: same pattern, 2 dirty paths
- Integrate into `do_move_impl` case 4

### Step 6: Cache invalidation logic
- After SPR/TBR accept: `nodeCL.valid = false`, rebuild TreeNav
- After global param change: `nodeCL.valid = false`
- After rate_loss/rate_neo: invalidate specific partitions
- After Gibbs/weighted accept: `nodeCL.valid = false`

### Step 7: Tests
- Unit: partial NNI eval matches full eval (random trees, multiple sizes)
- Unit: partial beta_simplex eval matches full eval
- Unit: cache invalidation triggers full eval correctly
- Unit: rollback produces identical state on rejection
- Integration: MCMC with partial CL produces valid posteriors
- Property: cache stays in sync after random move sequences

### Step 8: Benchmark
- Sun2018 (54 tips): measure iter/s improvement
- Compare NNI-only and mixed-move configurations
- Profile to confirm pruning time reduction

## Risks and mitigations

| Risk | Mitigation |
|------|-----------|
| Cache gets out of sync | Assertion: periodic full-eval comparison in debug mode |
| ACRV memory overhead | ~2-8 MB per chain; acceptable |
| Q-het complexity | Fall back to full eval (already scoped out) |
| SPR acceptance invalidates cache | Expected; full eval on next iteration rebuilds it |
| kPrime change resizes partition | Invalidate that partition's cache; handled by existing headroom logic |

## Files to modify

| File | Changes |
|------|---------|
| `src/mcmc_state.h` | Add `NodeCLCache` struct and member |
| `src/mcmc.cpp` | Allocation, do_move_impl routing, invalidation |
| `src/mcmc_likelihood.cpp` | Optional write-through to external cache buffer |
| `src/node_cl_cache.h` (new) | Partial eval functions, dirty node helpers |
| `tests/testthat/test-node-cl-cache.R` (new) | Validation tests |
