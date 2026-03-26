// C++ implementations of MCMC tree topology proposals (NNI, SPR).
//
// NNI: for an internal edge (u, v), swap one child of v with one child of u.
// All operations on integer parent/child vectors; postorder reordering via
// TreeTools::postorder_order() rather than ape::reorder.phylo().
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

  // Absolute branch lengths for reorder-safe reconstruction
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; i++) {
    absLen[i] = treeLength * relBrLengths[i];
  }

  // Build temporary edge matrix for postorder_order()
  IntegerMatrix tmpEdge(nEdge, 2);
  for (int i = 0; i < nEdge; i++) {
    tmpEdge(i, 0) = newParent[i];
    tmpEdge(i, 1) = newChild[i];
  }
  IntegerVector order = TreeTools::postorder_order(tmpEdge);

  IntegerVector ordParent(nEdge), ordChild(nEdge);
  NumericVector orderedRelBr(nEdge);
  for (int i = 0; i < nEdge; i++) {
    int j = order[i] - 1;
    ordParent[i] = newParent[j];
    ordChild[i] = newChild[j];
    orderedRelBr[i] = absLen[j] / treeLength;
  }

  return List::create(_["parent"] = ordParent,
                      _["child"] = ordChild,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = 0.0);
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
