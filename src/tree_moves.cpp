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

#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
#include <vector>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// NNI proposal — vector-based implementation
// ---------------------------------------------------------------------------

List nni_proposal_impl(IntegerVector parent, IntegerVector child,
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
  int pickInternal = (int)(unif_rand() * (double)internalRows.size());
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
  int pickV = (int)(unif_rand() * (double)vChildRows.size());
  if (pickV >= (int)vChildRows.size()) pickV = vChildRows.size() - 1;
  int pickU = (int)(unif_rand() * (double)uSibRows.size());
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
  List result = nni_proposal_impl(parent, child, nTip, treeLength,
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
