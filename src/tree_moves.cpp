// C++ implementations of MCMC tree topology proposals (NNI, SPR).
//
// NNI: for an internal edge (u, v), swap one child of v with one child of u.
// All operations on integer parent/child vectors; canonical preorder
// reordering via TreeTools::preorder_weighted_impl() (children sorted by
// smallest descendant, nodes renumbered in visit order).
//
// M-065: Added _impl versions that accept parent/child vectors directly
// (called from do_move_impl). The Rcpp-exported versions remain as thin
// wrappers that decompose the edge matrix and call _impl.
// T-017-IIa: impl functions take ChainRng& rng; R-exported wrappers seed from R.

#include <Rcpp.h>
#include "chain_rng.h"
#include <TreeTools/renumber_tree.h>
#include <vector>

using namespace Rcpp;

// Helper: construct a one-shot ChainRng from R's RNG (for R-exported wrappers).
// Draws exactly one R-RNG value.
static inline ChainRng rng_from_r() {
  uint64_t seed = static_cast<uint64_t>(
    unif_rand() * static_cast<double>(std::numeric_limits<uint32_t>::max()));
  return ChainRng(seed);
}


// ---------------------------------------------------------------------------
// NNI proposal — vector-based implementation
// ---------------------------------------------------------------------------

List nni_proposal_impl(ChainRng& rng,
                       IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths) {
  int nEdge = parent.size();

  // Find internal edges: both endpoints are internal nodes (> nTip)
  std::vector<int> internalRows;
  internalRows.reserve(nEdge / 2);
  for (int i = 0; i < nEdge; i++) {
    if (parent[i] > nTip && child[i] > nTip) {
      internalRows.push_back(i);
    }
  }

  if (internalRows.empty()) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // Pick a random internal edge
  int pickInternal = (int)(rng.unif() * (double)internalRows.size());
  if (pickInternal >= (int)internalRows.size())
    pickInternal = internalRows.size() - 1;
  int edgeRow = internalRows[pickInternal];
  int u = parent[edgeRow];
  int v = child[edgeRow];

  // Find v's children and u's other children (not v)
  std::vector<int> vChildRows, uSibRows;
  for (int i = 0; i < nEdge; i++) {
    if (parent[i] == v) {
      vChildRows.push_back(i);
    } else if (parent[i] == u && child[i] != v) {
      uSibRows.push_back(i);
    }
  }

  if (vChildRows.empty() || uSibRows.empty()) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // Pick one child of v and one sibling of v to swap
  int pickV = (int)(rng.unif() * (double)vChildRows.size());
  if (pickV >= (int)vChildRows.size()) pickV = vChildRows.size() - 1;
  int pickU = (int)(rng.unif() * (double)uSibRows.size());
  if (pickU >= (int)uSibRows.size()) pickU = uSibRows.size() - 1;

  int cRow = vChildRows[pickV];
  int wRow = uSibRows[pickU];

  // Swap: move cRow's subtree to u, wRow's subtree to v
  IntegerVector newParent = clone(parent);
  IntegerVector newChild = clone(child);
  newParent[cRow] = u;
  newParent[wRow] = v;

  // Absolute branch lengths; preorder_weighted_impl reorders both simultaneously.
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; i++) absLen[i] = treeLength * relBrLengths[i];

  auto po = TreeTools::preorder_weighted_impl(newParent, newChild, absLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector ordParent(nEdge), ordChild(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    ordParent[i] = ordEdge(i, 0);
    ordChild[i]  = ordEdge(i, 1);
  }
  NumericVector orderedRelBr = ordAbs / treeLength;

  return List::create(_["parent"] = ordParent,
                      _["child"] = ordChild,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = 0.0);
}


// ---------------------------------------------------------------------------
// Subtree swap (M-084)
//
// swap_subtrees_impl: exchange the positions of two subtrees.
//   nodeA and nodeB are child nodes (each has a parent edge).
//   Their parent connections are swapped; branch lengths swap with them so
//   total tree length is unchanged and the Jacobian is 1 (logHastings = 0).
//   Result is reordered to canonical preorder.
//   Caller must ensure nodeA/nodeB are non-nested, non-sibling, non-root
//   (use get_valid_swap_partners_impl).
//
// get_valid_swap_partners_impl: return all nodes w that are valid swap
//   partners for pruneNode — i.e. non-descendant, non-ancestor, non-sibling
//   of pruneNode (and pruneNode must not be the root).
// ---------------------------------------------------------------------------

// Returns the index of the edge row where child[i] == node, or -1.
static int find_child_row(const IntegerVector& child, int node) {
  for (int i = 0; i < child.size(); ++i)
    if (child[i] == node) return i;
  return -1;
}


// Internal impl: used by GibbsSubtreeSwap (M-086) and WeightedSubtreeSwap (M-089).
List swap_subtrees_impl(IntegerVector parent, IntegerVector child,
                        int nTip, double treeLength,
                        NumericVector relBrLengths,
                        int nodeA, int nodeB) {
  int nEdge = parent.size();
  int rowA  = find_child_row(child, nodeA);
  int rowB  = find_child_row(child, nodeB);

  if (rowA < 0 || rowB < 0 || rowA == rowB) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // Swap parent assignments and branch lengths
  IntegerVector newParent  = clone(parent);
  NumericVector newRelBr   = clone(relBrLengths);
  newParent[rowA] = parent[rowB];
  newParent[rowB] = parent[rowA];
  newRelBr[rowA]  = relBrLengths[rowB];
  newRelBr[rowB]  = relBrLengths[rowA];

  // Reorder to canonical preorder; preorder_weighted_impl handles both
  // topology and branch lengths in one pass.
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i) absLen[i] = treeLength * newRelBr[i];

  auto po = TreeTools::preorder_weighted_impl(newParent, child, absLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector ordParent(nEdge), ordChild(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    ordParent[i] = ordEdge(i, 0);
    ordChild[i]  = ordEdge(i, 1);
  }
  NumericVector ordRelBr = ordAbs / treeLength;

  return List::create(_["parent"]        = ordParent,
                      _["child"]         = ordChild,
                      _["rel_br_lengths"] = ordRelBr,
                      _["logHastings"]   = 0.0);
}


// Internal impl: enumerate valid swap partners for pruneNode.
std::vector<int> get_valid_swap_partners_impl(
    const IntegerVector& parent, const IntegerVector& child,
    int nTip, int pruneNode) {

  int nEdge  = parent.size();
  int pruneRow = find_child_row(child, pruneNode);
  if (pruneRow < 0) return {};   // pruneNode is root — no valid swaps

  int pruneParent = parent[pruneRow];
  int maxIdx = 2 * nTip + 2;    // safe upper bound on node indices

  // 1. Descendants of pruneNode (BFS)
  std::vector<bool> isDesc(maxIdx, false);
  isDesc[pruneNode] = true;
  if (pruneNode > nTip) {
    std::vector<int> q = {pruneNode};
    while (!q.empty()) {
      int cur = q.back(); q.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (parent[i] == cur && !isDesc[child[i]]) {
          isDesc[child[i]] = true;
          if (child[i] > nTip) q.push_back(child[i]);
        }
      }
    }
  }

  // 2. Ancestors of pruneNode (walk parent chain to root)
  std::vector<bool> isAnc(maxIdx, false);
  {
    int cur = pruneParent;
    while (cur > 0 && cur < maxIdx) {
      isAnc[cur] = true;
      int row = find_child_row(child, cur);
      if (row < 0) break;      // cur is root (no parent edge)
      cur = parent[row];
    }
  }

  // 3. Collect valid partners
  std::vector<int> partners;
  partners.reserve(nEdge / 2);
  for (int i = 0; i < nEdge; ++i) {
    int w = child[i];
    if (isDesc[w]) continue;               // descendant (incl. pruneNode)
    if (isAnc[w]) continue;                // ancestor
    if (parent[i] == pruneParent) continue; // sibling
    partners.push_back(w);
  }
  return partners;
}


// ---------------------------------------------------------------------------
// TBR proposal (M-053)
//
// Tree Bisection and Reconnection: SPR + subtree re-rooting.
//
//   1. Prune edge (u → v), suppress u (same as SPR).
//   2. If v is internal: re-root v's subtree at a random internal edge
//      (x → y) by moving v from its current position to between x and y.
//      This changes the rooted representation but preserves the unrooted
//      topology of the pruned subtree.
//   3. Regraft u on a random edge in the remaining tree (same as SPR).
//
// Hastings ratio:
//   logHR = log(lRegraft/lMerge)          [SPR Jacobian]
//         + log(lSubEdge/lMergeSub)       [subtree re-root Jacobian]
//
// The candidate-count ratio log(nRegraftReverse/nRegraftForward) is
// identically 0: re-rooting preserves |desc(v)| and u always has 3
// adjacent edges, so both counts equal nEdge - |desc(v)| - 2.
//
// When v is a tip the re-rooting is skipped and TBR degenerates to SPR.
// ---------------------------------------------------------------------------

// Helper: find edge row where child[i] == node.
static int tbr_find_child_row(const IntegerVector& child, int node) {
  for (int i = 0; i < child.size(); ++i)
    if (child[i] == node) return i;
  return -1;
}

List tbr_proposal_impl(ChainRng& rng,
                       IntegerVector parent, IntegerVector child,
                        int nTip, double treeLength,
                        NumericVector relBrLengths) {
  const int nEdge = parent.size();
  const int root  = nTip + 1;

  auto fail = [&]() {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  };

  // --- Phase A: Prune (same as SPR) ---

  // Eligible prune edges: parent != root
  std::vector<int> eligiblePrune;
  eligiblePrune.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (parent[i] != root) eligiblePrune.push_back(i);
  if (eligiblePrune.empty()) return fail();

  int pickPrune = (int)(rng.unif() * (double)eligiblePrune.size());
  if (pickPrune >= (int)eligiblePrune.size())
    pickPrune = (int)eligiblePrune.size() - 1;
  const int pruneRow = eligiblePrune[pickPrune];
  const int u = parent[pruneRow];
  const int v = child[pruneRow];

  // Find parentRow (p → u) and sibRow (u → w)
  int parentRow = -1, sibRow = -1, w = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (child[i] == u) parentRow = i;
    if (parent[i] == u && child[i] != v) { sibRow = i; w = child[i]; }
  }
  if (parentRow < 0 || sibRow < 0) return fail();

  // BFS: mark all descendants of v, collect subtree edge rows
  const int maxNode = 2 * nTip + 2;
  std::vector<bool> isDesc(maxNode, false);
  isDesc[v] = true;
  std::vector<int> subEdgeRows;   // edges within v's subtree
  if (v > nTip) {
    std::vector<int> queue = {v};
    while (!queue.empty()) {
      int cur = queue.back(); queue.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (parent[i] == cur && isDesc[parent[i]]) {
          isDesc[child[i]] = true;
          subEdgeRows.push_back(i);
          if (child[i] > nTip) queue.push_back(child[i]);
        }
      }
    }
  }

  // Absolute branch lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = treeLength * relBrLengths[i];

  const double lMerge = absLen[parentRow] + absLen[sibRow];

  // --- Phase B: Re-root pruned subtree ---

  double lSubEdge  = 0.0;
  double lMergeSub = 0.0;
  const int nSubEdge = (int)subEdgeRows.size();

  // Clone vectors for modification
  IntegerVector newParent = clone(parent);
  IntegerVector newChild  = clone(child);
  NumericVector newAbsLen = clone(absLen);

  if (v > nTip && nSubEdge > 0) {
    // Pick a random subtree edge
    int pickSub = (int)(rng.unif() * (double)nSubEdge);
    if (pickSub >= nSubEdge) pickSub = nSubEdge - 1;
    const int subRow = subEdgeRows[pickSub];
    const int x = parent[subRow];
    lSubEdge = absLen[subRow];

    // Identify v's two children
    int vChildRow1 = -1, vChildRow2 = -1;
    for (int i = 0; i < nEdge; ++i) {
      if (parent[i] == v) {
        if (vChildRow1 < 0) vChildRow1 = i;
        else vChildRow2 = i;
      }
    }

    lMergeSub = absLen[vChildRow1] + absLen[vChildRow2];

    if (x != v) {
      // General case: re-root subtree at edge (x → y).
      // Find path from v to x by walking upward from x.
      std::vector<int> pathNodes;  // x, aₘ, ..., a₁ (excludes v)
      {
        int cur = x;
        while (cur != v) {
          pathNodes.push_back(cur);
          int row = tbr_find_child_row(child, cur);
          if (row < 0) break;
          cur = parent[row];
        }
      }
      // pathNodes is [x, aₘ, ..., a₁] (bottom to top, excluding v)
      // The top of the path (last element) is a₁, which is v's direct child.

      if (pathNodes.empty()) {
        // x == v after all (shouldn't happen, but guard)
        lMergeSub = lSubEdge;
      } else {
        const int a1 = pathNodes.back();  // v's child on the path

        // Determine which vChildRow corresponds to a₁
        int vPathRow = -1, vOtherRow = -1;
        if (child[vChildRow1] == a1) {
          vPathRow = vChildRow1; vOtherRow = vChildRow2;
        } else {
          vPathRow = vChildRow2; vOtherRow = vChildRow1;
        }

        const double sigma = rng.unif();

        // Transform 1: v → c_other becomes a₁ → c_other
        //   parent changes v → a₁, length += L(v → a₁)
        newParent[vOtherRow] = a1;
        newAbsLen[vOtherRow] = absLen[vPathRow] + absLen[vOtherRow];

        // Transform 2: v → a₁ becomes v → x
        //   child changes a₁ → x, length = sigma × L(x → y)
        newChild[vPathRow] = x;
        newAbsLen[vPathRow] = sigma * lSubEdge;

        // Transform 3: reverse each interior path edge aᵢ → aᵢ₊₁
        // pathNodes = [x, aₘ, ..., a₂, a₁] (bottom to top)
        // Interior edges connect consecutive pairs from a₁ down to aₘ→x.
        // We reverse edges between pathNodes[k+1] → pathNodes[k] for
        // k = 0..(len-2).
        // pathNodes[len-1] = a₁, pathNodes[len-2] = a₂, ..., pathNodes[0] = x
        for (int k = 0; k < (int)pathNodes.size() - 1; ++k) {
          // Edge was pathNodes[k+1] → pathNodes[k], reverse it
          int fromNode = pathNodes[k + 1];
          int toNode   = pathNodes[k];
          // Find the edge row: parent == fromNode, child == toNode
          for (int i = 0; i < nEdge; ++i) {
            if (newParent[i] == fromNode && newChild[i] == toNode) {
              newParent[i] = toNode;
              newChild[i]  = fromNode;
              // length unchanged
              break;
            }
          }
        }

        // Transform 5: x → y becomes v → y
        //   parent changes x → v, length = (1 − sigma) × L(x → y)
        newParent[subRow] = v;
        newAbsLen[subRow] = (1.0 - sigma) * lSubEdge;
      }
    } else {
      // x == v: chosen edge is directly below v.
      // No topology change; just redistribute branch lengths with sigma.
      // Draw sigma to consume a RNG call (matching the forward/reverse symmetry)
      // but discard it — no topology or branch-length change when x == v.
      (void)rng.unif();
      // subRow is the chosen edge (v → y). The "other" child of v is the
      // other edge. We redistribute lSubEdge between the chosen edge and
      // combine with the other to define lMergeSub.
      // Since x == v, no path reversal. The only change is:
      // v → y gets length sigma × lSubEdge (already in subRow).
      // This is actually a no-op for topology. But we still need the
      // Jacobian to account for the sigma draw.
      // Set lMergeSub = lSubEdge so the sub-Jacobian = 0 (cancels).
      lMergeSub = lSubEdge;
    }
  } else {
    // v is a tip or subtree has no edges: TBR degenerates to SPR.
    lSubEdge  = 1.0;
    lMergeSub = 1.0;  // ratio = 1, log = 0
  }

  // --- Phase C: Regraft (same as SPR) ---

  // Candidate regraft edges: not in v's subtree, not adjacent to u
  std::vector<int> candidates;
  candidates.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[newChild[i]]) continue;
    if (newParent[i] == u || newChild[i] == u) continue;
    candidates.push_back(i);
  }
  if (candidates.empty()) return fail();
  const int nCand = (int)candidates.size();

  int pickRegraft = (int)(rng.unif() * (double)nCand);
  if (pickRegraft >= nCand) pickRegraft = nCand - 1;
  const int regraftRow = candidates[pickRegraft];
  const int b = newChild[regraftRow];

  const double tau = rng.unif();
  const double lRegraft = newAbsLen[regraftRow];

  // 1. Suppress u: (p -> u) becomes (p -> w)
  newChild[parentRow]  = w;
  newAbsLen[parentRow] = lMerge;

  // 2. Insert u on regraft edge: (a -> b) becomes (a -> u)
  newChild[regraftRow]  = u;
  newAbsLen[regraftRow] = tau * lRegraft;

  // 3. Reuse sibRow for (u -> b)
  newParent[sibRow] = u;
  newChild[sibRow]  = b;
  newAbsLen[sibRow]  = (1.0 - tau) * lRegraft;

  // Canonical preorder reordering
  auto po = TreeTools::preorder_weighted_impl(newParent, newChild, newAbsLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector ordParent(nEdge), ordChild(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    ordParent[i] = ordEdge(i, 0);
    ordChild[i]  = ordEdge(i, 1);
  }
  NumericVector ordRelBr = ordAbs / treeLength;

  // --- Hastings ratio ---
  // nRegraftReverse == nRegraftForward always: re-rooting preserves |desc(v)|,
  // u always has 3 adjacent edges in a binary tree, and the overlap (u→v) is
  // always 1.  So candidates = nEdge − |desc(v)| − 2 in both directions and
  // the candidate-count ratio term is identically 0.
  const double logHastings =
      std::log(lRegraft) - std::log(lMerge)
    + std::log(lSubEdge) - std::log(lMergeSub);

  return List::create(_["parent"]         = ordParent,
                      _["child"]          = ordChild,
                      _["rel_br_lengths"] = ordRelBr,
                      _["logHastings"]    = logHastings);
}


// Rcpp-exported wrapper for TBR (R testing)
// [[Rcpp::export]]
List tbr_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths) {
  int nEdge = edge.nrow();
  IntegerVector par(nEdge), ch(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    par[i] = edge(i, 0);
    ch[i]  = edge(i, 1);
  }
  ChainRng rng = rng_from_r();
  List result = tbr_proposal_impl(rng, par, ch, nTip, treeLength, relBrLengths);
  IntegerVector rp = result["parent"];
  IntegerVector rc = result["child"];
  IntegerMatrix outEdge(nEdge, 2);
  for (int i = 0; i < nEdge; ++i) {
    outEdge(i, 0) = rp[i];
    outEdge(i, 1) = rc[i];
  }
  return List::create(_["edge"]           = outEdge,
                      _["rel_br_lengths"] = result["rel_br_lengths"],
                      _["logHastings"]    = result["logHastings"]);
}


// Rcpp-exported wrapper: swap by node IDs (R testing / M-086/M-089 dispatch).
// [[Rcpp::export]]
List swap_subtrees_cpp(IntegerMatrix edge, int nTip, double treeLength,
                       NumericVector relBrLengths, int nodeA, int nodeB) {
  int nEdge = edge.nrow();
  IntegerVector par(nEdge), ch(nEdge);
  for (int i = 0; i < nEdge; ++i) { par[i] = edge(i,0); ch[i] = edge(i,1); }

  List res = swap_subtrees_impl(par, ch, nTip, treeLength, relBrLengths,
                                 nodeA, nodeB);
  // Reconstruct edge matrix for R compatibility
  IntegerVector rp = res["parent"], rc = res["child"];
  IntegerMatrix outEdge(nEdge, 2);
  for (int i = 0; i < nEdge; ++i) { outEdge(i,0) = rp[i]; outEdge(i,1) = rc[i]; }
  return List::create(_["edge"]           = outEdge,
                      _["rel_br_lengths"] = res["rel_br_lengths"],
                      _["logHastings"]    = res["logHastings"]);
}


// Rcpp-exported wrapper: return valid swap partners for pruneNode.
// [[Rcpp::export]]
IntegerVector get_valid_swap_partners_cpp(IntegerMatrix edge, int nTip,
                                           int pruneNode) {
  int nEdge = edge.nrow();
  IntegerVector par(nEdge), ch(nEdge);
  for (int i = 0; i < nEdge; ++i) { par[i] = edge(i,0); ch[i] = edge(i,1); }
  std::vector<int> v = get_valid_swap_partners_impl(par, ch, nTip, pruneNode);
  return IntegerVector(v.begin(), v.end());
}


// Rcpp-exported wrapper (for R callers via ProposeNni)
// [[Rcpp::export]]
List nni_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths) {
  int nEdge = edge.nrow();
  IntegerVector parent(nEdge), child(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    parent[i] = edge(i, 0);
    child[i] = edge(i, 1);
  }
  ChainRng rng = rng_from_r();
  List result = nni_proposal_impl(rng, parent, child, nTip, treeLength,
                                  relBrLengths);

  // Reconstruct edge matrix for R-side compatibility
  IntegerVector rp = result["parent"];
  IntegerVector rc = result["child"];
  IntegerMatrix outEdge(nEdge, 2);
  for (int i = 0; i < nEdge; ++i) {
    outEdge(i, 0) = rp[i];
    outEdge(i, 1) = rc[i];
  }
  return List::create(_["edge"] = outEdge,
                      _["rel_br_lengths"] = result["rel_br_lengths"],
                      _["logHastings"] = result["logHastings"]);
}
