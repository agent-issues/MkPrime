// C++ implementations of MCMC tree topology proposals (NNI, SPR).
//
// NNI: for an internal edge (u, v), swap one child of v with one child of u.
// All operations on integer parent/child vectors; postorder reordering via
// TreeTools::postorder_order() rather than ape::reorder.phylo().

#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
#include <vector>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// NNI proposal (Nearest Neighbour Interchange)
// ---------------------------------------------------------------------------
//
// Inputs (all in ape edge-matrix convention, 1-based node indices):
//   edge         nEdge x 2 integer matrix (parent | child)
//   nTip         number of tips
//   treeLength   current total tree length
//   relBrLengths nEdge-length vector of relative branch lengths (sum = 1)
//
// Returns a list with:
//   edge         new edge matrix in postorder
//   rel_br_lengths  new relative branch lengths (consistent with new row order)
//   logHastings  0 (symmetric proposal) or -Inf if NNI is not applicable

// [[Rcpp::export]]
List nni_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths) {
  int nEdge = edge.nrow();

  // --- Find internal edges: both endpoints are internal nodes (> nTip) ---
  std::vector<int> internalRows;
  internalRows.reserve(nEdge / 2);
  for (int i = 0; i < nEdge; i++) {
    if (edge(i, 0) > nTip && edge(i, 1) > nTip) {
      internalRows.push_back(i);
    }
  }

  if (internalRows.empty()) {
    // Too few tips for NNI (n <= 3); return unchanged with -Inf Hastings
    return List::create(_["edge"] = edge,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // --- Pick a random internal edge ---
  int pickInternal = (int)(unif_rand() * (double)internalRows.size());
  if (pickInternal >= (int)internalRows.size()) pickInternal = internalRows.size() - 1;
  int edgeRow = internalRows[pickInternal];
  int u = edge(edgeRow, 0);
  int v = edge(edgeRow, 1);

  // --- Find v's children and u's other children (not v) ---
  std::vector<int> vChildRows, uSibRows;
  for (int i = 0; i < nEdge; i++) {
    if (edge(i, 0) == v) {
      vChildRows.push_back(i);
    } else if (edge(i, 0) == u && edge(i, 1) != v) {
      uSibRows.push_back(i);
    }
  }

  if (vChildRows.empty() || uSibRows.empty()) {
    // Degenerate tree — should not happen on a binary unrooted tree
    return List::create(_["edge"] = edge,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // --- Pick one child of v and one sibling of v (child of u) to swap ---
  int pickV = (int)(unif_rand() * (double)vChildRows.size());
  if (pickV >= (int)vChildRows.size()) pickV = vChildRows.size() - 1;
  int pickU = (int)(unif_rand() * (double)uSibRows.size());
  if (pickU >= (int)uSibRows.size()) pickU = uSibRows.size() - 1;

  int cRow = vChildRows[pickV];  // row of v's selected child
  int wRow = uSibRows[pickU];    // row of u's selected other child

  // --- Swap: move cRow's subtree to u, wRow's subtree to v ---
  IntegerMatrix newEdge = clone(edge);
  newEdge(cRow, 0) = u;
  newEdge(wRow, 0) = v;

  // --- Compute absolute branch lengths for reorder-safe reconstruction ---
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; i++) {
    absLen[i] = treeLength * relBrLengths[i];
  }

  // --- Reorder to postorder using TreeTools::postorder_order() ---
  IntegerVector order = TreeTools::postorder_order(newEdge);

  IntegerMatrix orderedEdge(nEdge, 2);
  NumericVector orderedRelBr(nEdge);
  for (int i = 0; i < nEdge; i++) {
    int j = order[i] - 1;  // 1-based → 0-based
    orderedEdge(i, 0) = newEdge(j, 0);
    orderedEdge(i, 1) = newEdge(j, 1);
    orderedRelBr[i] = absLen[j] / treeLength;
  }

  return List::create(_["edge"] = orderedEdge,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = 0.0);
}
