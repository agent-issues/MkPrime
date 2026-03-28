# M-053: TBR topology proposal for unrooted binary trees

**Date:** 2026-03-28
**Status:** Draft

---

## Setup: create worktree

```bash
cd mkp
git worktree add ../mkp-tbr feature/tbr-move
```

Update worktree mapping in `../AGENTS.md`. Merge from main before
starting code changes.

---

## 1. What TBR does (conceptual)

TBR = Tree Bisection and Reconnection. Given an unrooted binary tree:

1. Remove one internal edge, splitting the tree into two components.
2. In each component, the endpoint of the removed edge becomes a
   degree-2 node; suppress it (merge its two remaining edges).
3. In each component, pick a random edge and re-insert the suppressed
   node at a random point along that edge.
4. Reconnect the two components via the re-inserted nodes, using the
   original removed edge.

SPR is the special case where one component is left unchanged (the
re-insertion goes back to the same position). TBR moves both sides.

---

## 2. Worked example on a 6-tip tree

### 2a. Starting tree (ape format, rooted at 7)

```
         7 (root)
        / \
       8   1
      / \
     9   2
    / \
  10   3
  / \
 4    5
```

Edges (parent → child):
```
row 0: 7→8   row 1: 7→1
row 2: 8→9   row 3: 8→2
row 4: 9→10  row 5: 9→3
row 6: 10→4  row 7: 10→5
```
Internal edges (both endpoints > 6): rows 0, 2, 4.

### 2b. Bisect at row 2: u=8, v=9

Remove edge 8→9. Two components:
- T_A (u=8 side): nodes {7, 8, 1, 2}. u=8 has remaining edges:
  row 0 (7→8) and row 3 (8→2). Suppress 8: merge into 7→2 with
  length l_0 + l_3. Free row: row 3.
- T_B (v=9 side): nodes {9, 10, 3, 4, 5}. v=9 has remaining edges:
  row 4 (9→10) and row 5 (9→3). Suppress 9: merge into 10→3 with
  length l_4 + l_5. Free row: row 5.

### 2c. Regraft

T_A after suppression has edges: {row 0 (now 7→2), row 1 (7→1)}.
T_B after suppression has edges: {row 4 (now 10→3), row 6 (10→4), row 7 (10→5)}.

**Regraft u=8 into T_A:** say we pick row 1 (7→1), τ_A = 0.4.
Split: row 1 becomes 7→8 (len = 0.4 × l_1), row 3 becomes 8→1
(len = 0.6 × l_1). Row 0 keeps merged edge 7→2 (len = l_0 + l_3).

Wait — that gives node 8 children {1} and the u–v edge. Where's the
edge to 2? Node 8 needs degree 3 (two children + one parent, or
three children if root). Let me re-examine.

**Problem:** After suppressing u=8, row 0 becomes 7→2 (the merged
edge). When we regraft u=8 by splitting row 1 (7→1), we get:
- row 1: 7→8 (parent half of split)
- row 3: 8→1 (child half, using the free row)
- row 0: 7→2 (merged edge, unchanged)
- bisectRow 2: 8→9 (reconnect)

Node 8 now has: parent 7 (via row 1), children 1 (row 3) and 9
(row 2). Degree = 3. ✓
Node 7 now has: children 8 (row 1) and 2 (row 0). Degree = 2 if
root, but root must have degree 3 for unrooted tree. ✗

**This is the root problem I hit before.** Node 7 = root = nTip+1
originally had 3 children (8 and 1). After suppressing 8, we merged
row 0 into 7→2 and row 1 stayed as 7→1, giving root 7 two children
{2, 1}. Then we split row 1 to insert 8 between 7 and 1, giving
root 7 children {2, 8}. That's degree 2, which is wrong.

### 2d. Why SPR avoids this

SPR restricts: `parent[pruneRow] != root`. The pruned node u always
has a parent edge AND a child edge; suppression merges those two,
keeping u's parent (or root) at the same degree. TBR bisects an
internal edge where u could BE the root's child — and suppressing u
changes the root's degree.

### 2e. Root-safe TBR design

**Key insight:** In ape's directed format, the root (nTip+1) has
exactly 3 children. Suppressing a child of the root drops it to 2
children, which is invalid. This happens when u = a direct child of
root.

Two options:
**(A)** Exclude bisection edges adjacent to the root (like SPR does).
**(B)** Handle the root case specially.

**Choice: (A) — exclude root-adjacent edges.** This means we only
bisect edges where NEITHER endpoint is the root. For an internal edge
(u,v) with u=parent, v=child, this means `parent[bisectRow] != root`.
Since `child[bisectRow]` is always a non-root internal node (only root
has no parent edge), we only need `parent[bisectRow] != root`.

**Impact on candidate count:** An unrooted binary tree with n tips has
n-3 internal edges. The root has 2 or 3 internal-node children,
contributing 2-3 root-adjacent internal edges. So we lose 2-3
candidates. For n≥8 this is mild; for n=5 (2 internal edges) it could
leave 0 candidates. Solution: return -Inf when no eligible edges exist
(same as NNI/SPR).

---

## 3. Algorithm (detailed)

**Precondition:** ape-format edge vectors (parent, child, relBrLengths),
where `parent[i]` is always an internal node, root = nTip+1 never
appears as a child, and edges are in some valid order.

```
tbr_proposal_impl(parent, child, nTip, treeLength, relBrLengths):

  root = nTip + 1

  1. Find eligible bisection edges:
     for each i: parent[i] > nTip AND child[i] > nTip AND parent[i] != root
     (Both endpoints internal, parent not root.)
     If empty → return -Inf.

  2. Pick random eligible edge. Let bisectRow = chosen index.
     u = parent[bisectRow], v = child[bisectRow].

  3. Find u's edges (excluding bisection edge):
     Since u != root, u has exactly:
       - 1 "parent edge" (upRow): child[upRow] == u
       - 1 "other child edge" (ucRow): parent[ucRow] == u && child[ucRow] != v
     Let p_u = parent[upRow], c_u = child[ucRow].
     l_mergeA = absLen[upRow] + absLen[ucRow].

  4. Find v's edges (excluding bisection edge):
     v is never root (v = child of bisection, never in parent-column-only).
     v has exactly:
       - 1 "parent edge": this IS the bisection edge (parent[bisectRow] == u, child == v).
         But we already excluded it. v's other connections:
       - 2 "child edges": parent[i] == v for exactly 2 rows.
     Let vcRow0, vcRow1 be those rows.
     c_v0 = child[vcRow0], c_v1 = child[vcRow1].
     l_mergeB = absLen[vcRow0] + absLen[vcRow1].

  5. BFS to partition nodes into T_A (u's side) and T_B (v's side),
     treating the tree as undirected but excluding the bisection edge.
     T_A: start from p_u, c_u, u. T_B: start from c_v0, c_v1, v.

  6. Build candidate regraft edges for each subtree:
     T_A candidates = edges fully within T_A, excluding upRow and ucRow,
     PLUS the virtual merged edge (p_u ↔ c_u) with length l_mergeA.
     T_B candidates = edges fully within T_B, excluding vcRow0 and vcRow1,
     PLUS the virtual merged edge (c_v0 ↔ c_v1) with length l_mergeB.

  7. Pick random regraft edge in each subtree. Sample τ_A, τ_B ~ U(0,1).

  8. Apply edge rewrites (see §4 below).

  9. Reorder via TreeTools::preorder_weighted_impl().

  10. logHastings = log(l_regraftA) + log(l_regraftB) - log(l_mergeA) - log(l_mergeB).
```

---

## 4. Edge rewrite cases

Six edge rows are affected: bisectRow, upRow, ucRow, vcRow0, vcRow1,
and (if regraft is not the merged edge) the regraft-target row. The
key invariant: after rewriting, every internal node must appear as
parent in exactly the right number of edges (3 for root, 2 for others)
and as child in exactly 1 edge (0 for root).

### 4a. Suppress u and regraft u

**Suppress u:** Merge upRow and ucRow into one edge p_u → c_u.
- Use upRow for the merged edge: `child[upRow] = c_u`, `len[upRow] = l_mergeA`.
- ucRow becomes a "free" row.

**Regraft u into T_A:**

*Case A1: regraft onto merged edge (p_u ↔ c_u).*
Re-insert u between p_u and c_u at fraction τ_A.
- upRow: p_u → u, len = τ_A × l_mergeA.
- ucRow: u → c_u, len = (1 - τ_A) × l_mergeA.
This is the identity case for T_A's topology (u goes back between its
original neighbors, but at a different position).

*Case A2: regraft onto a different T_A edge (regRow, with p_r → c_r).*
The merged edge stays: upRow keeps p_u → c_u, len = l_mergeA.
Split regRow: regRow becomes p_r → u, len = τ_A × l_regraftA.
ucRow (free): u → c_r, len = (1 - τ_A) × l_regraftA.

**Check degrees (Case A2):**
- u: parent p_r (via regRow), children c_r (ucRow) and v (bisectRow). Degree 3. ✓
- p_u: was parent of u (via upRow), now parent of c_u. Degree unchanged. ✓
- p_r: was parent of c_r, now parent of u. Degree unchanged. ✓

### 4b. Suppress v and regraft v

**Suppress v:** Merge vcRow0 and vcRow1 into one edge.
Need to choose direction. Since both c_v0 and c_v1 could be tips or
internal, and we need parent > nTip: if one is a tip, the other must
be parent; if both are internal, either works (preorder_weighted_impl
will fix it). Safe choice: pick max(c_v0, c_v1) as parent.
- vcRow0: max(c_v0,c_v1) → min(c_v0,c_v1), len = l_mergeB.
- vcRow1 becomes free.

**Regraft v into T_B:**

*Case B1: regraft onto merged edge.*
- vcRow0: v → min(c_v0,c_v1), len = τ_B × l_mergeB.
  Wait — v needs to be a child of something. But v's only parent edge
  is the bisection edge (u → v). So v appears as child only in
  bisectRow. In T_B, v is the "local root" of the subtree; its edges
  are all parent edges (v → ...).
  
  Re-insert v between c_v0 and c_v1:
- vcRow0: v → c_v0, len = τ_B × l_mergeB.
- vcRow1: v → c_v1, len = (1 - τ_B) × l_mergeB.
  Node v: parent u (bisectRow), children c_v0 and c_v1. Degree 3. ✓

*Case B2: regraft onto a different T_B edge (regRowB, with p_s → c_s).*
The merged edge stays: vcRow0 keeps merged direction, len = l_mergeB.
Split regRowB: regRowB becomes p_s → v, len = τ_B × l_regraftB.
vcRow1 (free): v → c_s, len = (1 - τ_B) × l_regraftB.

**Check v's degree (Case B2):**
- v: parent p_s (via regRowB)... wait, v already has parent u via
  bisectRow. Two parent edges for v — that's illegal!

**This is wrong.** v is a child of u (via bisectRow). If we also make
v a child of p_s (via split regRowB), v has two parents. The fix: when
splitting regRowB, the direction should be v → c_s and p_s → v, but v
being child of p_s conflicts with v being child of u.

The correct approach: split so that v becomes the parent of c_s, and
the upper half of the split still points to v's old parent p_s.
Concretely:
- regRowB: p_s → v, len = τ_B × l_regraftB.
  This DOES make v a child of p_s, which conflicts with bisectRow.

Hmm. Let me reconsider. In the final tree, v has exactly one parent
edge (bisectRow: u → v) and exactly two child edges. So when regrafting
v into T_B, v should appear as PARENT in two new edges, not as a child.

The issue is that in T_B (as a connected component), v was the root of
that subtree (its only parent was u via the bisection edge, which was
removed). So all paths in T_B flow from v downward. When we suppress v
and then regraft v, we're re-inserting v as an internal node in T_B
where v has two downward edges.

Correct regraft for Case B2:
- After suppressing v, the merged edge replaces the path through v.
  Let's say vcRow0 holds the merged edge.
- To regraft v onto edge (p_s → c_s):
  regRowB: p_s → v (v becomes child of p_s for this edge).
  But v is ALSO child of u via bisectRow. Two parents. ✗

**The real problem:** In the original tree, v's parent is u (via
bisectRow). After TBR, v's parent might be a different node if v moved
within T_B. But bisectRow always says u → v. So v can't change its
parent in T_B — its parent is always u.

**Resolution:** After TBR, bisectRow should connect u to v with
u → v. In T_B, v is the local root. All of v's connections within T_B
must be v → (something), never (something) → v. So when we split an
edge (p_s → c_s) to insert v, the result must be (p_s → c_s still, but
going through v), which means: p_s → v and v → c_s. But p_s → v means
v is a child of p_s, conflicting with u → v.

**This means v cannot be regrafted as an intermediary on a directed
path that flows toward it.** In T_B, all edges flow away from v (v is
T_B's root). Suppressing v merges v's two child edges. Regrafting v
back means inserting v on some edge that flows away from v — but that
creates a cycle in the directed sense.

**Correct reformulation:** When we regraft v into T_B, we don't insert
v as an intermediary on an existing directed edge. Instead, we:
1. Detach one subtree from the merged edge and re-parent it under v.
2. v remains the root of T_B (child of u via bisectRow).

Concretely for Case B2 (regraft onto edge p_s → c_s):
- The merged edge stays (vcRow0).
- "Cut" edge p_s → c_s: detach c_s from p_s.
- Re-parent c_s under v: vcRow1 becomes v → c_s.
- Maintain p_s → c_s row... but c_s is detached. Instead:
  regRowB becomes p_s → v... no, that gives v two parents again.

I think the issue is that for the v-side, because v is T_B's root,
the regraft operation is fundamentally different from the u-side.
For u, we can use the same SPR-like "insert on an edge" technique
because u is an interior node. For v, we need a different operation
that re-roots a portion of T_B under v.

**Simpler correct approach:** Decompose TBR as two sequential SPR-like
operations, respecting the directed structure:
1. SPR-prune u from T_A (suppress u, merge upRow+ucRow), then
   SPR-regraft u somewhere in T_A.
2. For T_B, the operation is: remove one of v's children (say the
   subtree at c_v0), reattach it elsewhere in T_B, creating a new
   internal node that becomes a child of v. This is just an SPR within
   T_B rooted at v.

Actually — that's not TBR. That's two independent SPRs. TBR is more
specific: it bisects and reconnects. Let me go back to the definition.

**Correct TBR re-formulation for directed trees:**

After bisecting at (u, v), we have two subtrees. The key insight:
we don't need to suppress and re-insert u and v. Instead:

1. Disconnect u from its neighbors (within T_A), then reconnect u at a
   new position in T_A. This is exactly an SPR-prune-and-regraft of u's
   subtree within T_A.
2. Similarly, disconnect v from its children and reconnect v's children
   at a new arrangement within T_B.

But that changes v's children, not v's position. Hmm.

**Let me re-read the standard TBR definition for unrooted trees,
working entirely in the undirected setting, then convert to directed
at the end.**

### TBR on an undirected tree:

1. Remove edge {u, v}. Two components T_A (containing u) and T_B
   (containing v).
2. In T_A, node u has degree 2 (was 3, lost the u-v edge). Suppress u:
   merge its two remaining edges into one.
3. In T_B, node v has degree 2. Suppress v similarly.
4. Now T_A and T_B are proper trees without u and v.
5. Pick any edge e_A in T_A, insert u on e_A (split at fraction τ_A).
   u gets degree 2 from the split.
6. Pick any edge e_B in T_B, insert v on e_B (split at τ_B). v gets
   degree 2.
7. Reconnect u-v. Both u and v go from degree 2 to degree 3.

This is clear in the undirected setting. The challenge is step 5-7 in
ape's directed format.

**Key realization:** After step 7, the undirected tree is well-defined.
We then need to ROOT it at node nTip+1 (the ape root) and assign
directions (parent → child). The function
`TreeTools::preorder_weighted_impl(parent, child, lengths)` does
exactly this — it takes an edge list (which may have inconsistent
directions) and produces a properly rooted, preorder-ordered tree.

So the strategy is:
1. Work with cloned parent/child/absLen vectors.
2. Perform the undirected graph manipulations (suppress, regraft)
   without worrying about edge directions.
3. Call `preorder_weighted_impl` to fix directions, renumber, and
   reorder.

**Does `preorder_weighted_impl` tolerate arbitrary edge directions?**
Looking at the NNI implementation (lines 78-96 of tree_moves.cpp):
it modifies `newParent` (swapping which node is the parent of which)
and then calls `preorder_weighted_impl(newParent, newChild, absLen)`.
The NNI swap can put a tip in the parent column temporarily — and
preorder_weighted_impl handles it.

This is confirmed by the subtree-swap code (M-084, lines 146-166): it
swaps parent assignments freely and relies on preorder_weighted_impl
to sort out the result.

**So yes — we can do arbitrary undirected manipulations on the edge
vectors and let preorder_weighted_impl fix everything.**

---

## 5. Revised algorithm (undirected, preorder-fixed)

```
tbr_proposal_impl(parent, child, nTip, treeLength, relBrLengths):

  root = nTip + 1

  1. Find eligible bisection edges:
     parent[i] > nTip AND child[i] > nTip AND parent[i] != root
     If empty → return -Inf.

  2. Pick random bisection edge. bisectRow, u = parent, v = child.

  3. Absolute lengths. Clone parent, child, absLen.

  4. Find u's other two edges (excluding bisectRow):
     (a) upRow: where child[i] == u. Neighbor: a1 = parent[upRow].
     (b) ucRow: where parent[i] == u && child[i] != v. Neighbor: a2 = child[ucRow].
     l_mergeA = absLen[upRow] + absLen[ucRow].

  5. Find v's two child edges (excluding bisectRow):
     vcRow0, vcRow1 where parent[i] == v.
     Neighbors: b1 = child[vcRow0], b2 = child[vcRow1].
     l_mergeB = absLen[vcRow0] + absLen[vcRow1].

  6. Suppress u: merge upRow + ucRow into one edge.
     newChild[upRow] = a2. newAbsLen[upRow] = l_mergeA.
     (ucRow is now free — its contents will be overwritten.)

  7. Suppress v: merge vcRow0 + vcRow1 into one edge.
     newParent[vcRow0] = b1. newChild[vcRow0] = b2.
     (Or b2→b1 — direction doesn't matter, preorder fixes it.)
     newAbsLen[vcRow0] = l_mergeB.
     (vcRow1 is now free.)

  8. BFS to find T_A edges and T_B edges (on the suppressed graph,
     excluding bisectRow, ucRow, vcRow1).

  9. Build candidate lists:
     candA = {edges fully in T_A, including the merged edge at upRow}
     candB = {edges fully in T_B, including the merged edge at vcRow0}

  10. Pick regraft edge in T_A (index regRowA, length l_regA).
      Pick regraft edge in T_B (index regRowB, length l_regB).
      Sample τ_A, τ_B ~ U(0,1).

  11. Regraft u into T_A:
      If regRowA == upRow (regraft onto merged edge):
        newChild[upRow] = u.        // a1 → u (upper half; a1 = parent[upRow])
        newAbsLen[upRow] = τ_A × l_mergeA.
        newParent[ucRow] = u.       // u → a2 (lower half)
        newChild[ucRow] = a2.
        newAbsLen[ucRow] = (1-τ_A) × l_mergeA.
      Else:
        c_r = newChild[regRowA].
        newChild[regRowA] = u.      // p_r → u
        newAbsLen[regRowA] = τ_A × l_regA.
        newParent[ucRow] = u.       // u → c_r
        newChild[ucRow] = c_r.
        newAbsLen[ucRow] = (1-τ_A) × l_regA.
        // upRow still holds the merged edge a1 → a2.

  12. Regraft v into T_B:
      If regRowB == vcRow0 (regraft onto merged edge):
        newParent[vcRow0] = v.      // v → b2 (or b1, doesn't matter)
        newChild[vcRow0] = b2.
        newAbsLen[vcRow0] = τ_B × l_mergeB.
        newParent[vcRow1] = v.      // v → b1
        newChild[vcRow1] = b1.
        newAbsLen[vcRow1] = (1-τ_B) × l_mergeB.
      Else:
        c_s = newChild[regRowB].
        newChild[regRowB] = v.      // p_s → v
        newAbsLen[regRowB] = τ_B × l_regB.
        newParent[vcRow1] = v.      // v → c_s
        newChild[vcRow1] = c_s.
        newAbsLen[vcRow1] = (1-τ_B) × l_regB.
        // vcRow0 still holds merged edge.

  13. bisectRow: newParent[bisectRow] = u, newChild[bisectRow] = v.
      (May already be correct; ensure it.)

  14. Reorder:
      auto po = TreeTools::preorder_weighted_impl(
          newParent, newChild, newAbsLen);
      Extract ordParent, ordChild, ordRelBr.

  15. logHastings = log(l_regA) + log(l_regB) - log(l_mergeA) - log(l_mergeB).

  16. Return.
```

### Why preorder_weighted_impl fixes the direction issue

In step 11 (Case A2), after `newChild[regRowA] = u`, u appears as a
child of p_r. But u is also the parent of v (bisectRow) and parent of
c_r (ucRow). In the directed sense, u has parent p_r and children
{c_r, v}. That's correct (degree 3, 1 parent + 2 children).

In step 12 (Case B2), after `newChild[regRowB] = v`, v appears as a
child of p_s AND as a child of u (bisectRow). Two parent edges! But
`preorder_weighted_impl` will figure out the correct rooting starting
from node nTip+1, ignoring the input directions. The edge between p_s
and v will be oriented correctly based on which side of the root each
node falls.

**Wait — does preorder_weighted_impl actually ignore input
directions?** Based on the NNI code, after swapping `newParent[cRow]=u`
and `newParent[wRow]=v`, a node can temporarily have the wrong parent.
The function traverses from the root and re-assigns directions.

If it doesn't handle arbitrary directions, we'd need to ensure v has
only one parent edge before calling it. In that case, step 12 Case B2
needs: instead of making v a child of p_s, we need to reverse the edge
so v is parent of p_s... but then p_s has two parents. The fundamental
issue is that any edge re-insertion creates a temporary directed
inconsistency that must be resolved by re-rooting.

**Verification needed:** Read the TreeTools C++ source for
`preorder_weighted_impl` to confirm it handles arbitrary/inconsistent
edge directions by treating the graph as undirected from the root.

---

## 6. Hastings ratio derivation

**Forward proposal density:**
q(T → T') = (1/n_elig) × (1/n_candA) × (1/l_regA) ×
             (1/n_candB) × (1/l_regB)

where the 1/l factors come from the Jacobian of τ ~ U(0,1) → split
edge of length l into (τl, (1-τ)l).

**Reverse proposal density:**
q(T' → T) = (1/n_elig') × (1/n_candA') × (1/l_mergeA) ×
             (1/n_candB') × (1/l_mergeB)

**Cancellations:**
- n_elig = n_elig': excluding root-adjacent edges, the count depends
  only on the tree size, not topology.
  
  Actually — n_elig depends on which edges are internal. A TBR move
  can change which edges are internal (an edge between two internal
  nodes might become an edge between an internal node and a tip after
  rearrangement). Wait no — TBR doesn't change node identities, only
  connections. All internal nodes remain internal. So n_elig = n_elig'. ✓

- n_candA = n_candA': T_A has the same tips after bisecting the same
  edge in reverse, so the same number of edges. ✓
- n_candB = n_candB': Same argument. ✓

**Result:**
logHR = log(l_regA) + log(l_regB) - log(l_mergeA) - log(l_mergeB)

---

## 7. Files changed

| File | Change |
|------|--------|
| `src/tree_moves.cpp` | Add `tbr_proposal_impl()` + `tbr_proposal()` Rcpp export |
| `src/mcmc.cpp` | Add case 16 to `do_move_impl()`, forward declaration |
| `R/RunMkPrime.R` | Add `tbr = 16L` to `.kMoveTypes`, add TBR to `.BuildMoves()` |
| `R/proposals.R` | Add `ProposeTbr()` R wrapper |
| `tests/testthat/test-tbr.R` | New test file |

Move type = 16 (next free after block_gibbs_branch = 15).

---

## 8. Open question: preorder_weighted_impl

Before implementing, verify that `preorder_weighted_impl` treats the
input as an undirected graph rooted at nTip+1. If it requires valid
directed input, the v-side regraft (Case B2) needs a different
approach — probably re-rooting T_B at v before regraft.

**Check:** read TreeTools `inst/include/TreeTools/renumber_tree.h` for
the `preorder_weighted_impl` signature and traversal logic.

---

## 9. Tests

| Test | Correctness criterion |
|------|----------------------|
| Valid tree | Output has correct nEdge, nTip tips, postorder-compatible |
| Total length preserved | sum(absLen) before = sum(absLen) after |
| Hastings formula | Known 6-tip example: hand-compute l_merge and l_reg, verify logH |
| 50 random proposals on 10-tip tree | All produce valid trees |
| 4-tip tree | Reject (only 1 internal edge, which is root-adjacent → 0 eligible) |
| 5-tip tree | At most 1 eligible edge; verify valid output or -Inf |
| MCMC smoke test | RunMkPrime with TBR, 200 iters, returns valid posterior |

---

## 10. Implementation order

1. Verify `preorder_weighted_impl` behavior (read-only).
2. C++ `tbr_proposal_impl` + `tbr_proposal` in `src/tree_moves.cpp`.
3. Unit tests for proposal validity.
4. Dispatch case 16 in `src/mcmc.cpp`.
5. R-side registration in `.kMoveTypes` and `.BuildMoves()`.
6. `ProposeTbr()` wrapper in `R/proposals.R`.
7. MCMC smoke test.
