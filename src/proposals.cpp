// C++ implementations of MCMC proposals (SPR, BetaSimplex).
//
// NNI is in tree_moves.cpp (nni_proposal).
// SPR: prune-and-regraft with Jacobian Hastings ratio.
// BetaSimplex: redistribute mass between two simplex elements.
//
// M-065: Added spr_proposal_impl that accepts parent/child vectors directly.

#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
#include <vector>
#include <cmath>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// SPR proposal — vector-based implementation
// ---------------------------------------------------------------------------

List spr_proposal_impl(IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths) {
  const int nEdge = parent.size();
  const int root = nTip + 1;

  // Eligible prune edges: parent != root
  std::vector<int> eligiblePrune;
  eligiblePrune.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (parent[i] != root) {
      eligiblePrune.push_back(i);
    }
  }

  if (eligiblePrune.empty()) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  int pickPrune = (int)(unif_rand() * (double)eligiblePrune.size());
  if (pickPrune >= (int)eligiblePrune.size())
    pickPrune = eligiblePrune.size() - 1;
  const int pruneRow = eligiblePrune[pickPrune];
  const int u = parent[pruneRow];
  const int v = child[pruneRow];

  // u's parent: row where child == u
  int parentRow = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (child[i] == u) {
      parentRow = i;
      break;
    }
  }

  // v's sibling: row where parent == u, child != v
  int sibRow = -1;
  int w = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (parent[i] == u && child[i] != v) {
      sibRow = i;
      w = child[i];
      break;
    }
  }

  if (parentRow < 0 || sibRow < 0) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // BFS to find all descendants of v
  const int maxNode = nEdge + nTip + 2;
  std::vector<bool> isDesc(maxNode, false);
  isDesc[v] = true;
  if (v > nTip) {
    std::vector<int> queue;
    queue.push_back(v);
    while (!queue.empty()) {
      const int cur = queue.back();
      queue.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (parent[i] == cur) {
          isDesc[child[i]] = true;
          if (child[i] > nTip) {
            queue.push_back(child[i]);
          }
        }
      }
    }
  }

  // Candidate regraft edges: not in v's subtree, not adjacent to u
  std::vector<int> candidates;
  candidates.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[child[i]]) continue;
    if (parent[i] == u || child[i] == u) continue;
    candidates.push_back(i);
  }

  if (candidates.empty()) {
    return List::create(_["parent"] = parent,
                        _["child"] = child,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  int pickRegraft = (int)(unif_rand() * (double)candidates.size());
  if (pickRegraft >= (int)candidates.size())
    pickRegraft = candidates.size() - 1;
  const int regraftRow = candidates[pickRegraft];
  const int b = child[regraftRow];

  const double tau = unif_rand();

  // Absolute branch lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    absLen[i] = treeLength * relBrLengths[i];
  }

  const double lRegraft = absLen[regraftRow];
  const double lMerge = absLen[parentRow] + absLen[sibRow];

  // Perform SPR on cloned vectors
  IntegerVector newParent = clone(parent);
  IntegerVector newChild = clone(child);
  NumericVector newAbsLen = clone(absLen);

  // 1. Suppress u: (p -> u) becomes (p -> w)
  newChild[parentRow] = w;
  newAbsLen[parentRow] = lMerge;

  // 2. Insert u on regraft edge: (a -> b) becomes (a -> u)
  newChild[regraftRow] = u;
  newAbsLen[regraftRow] = tau * lRegraft;

  // 3. Reuse sibRow for (u -> b)
  newParent[sibRow] = u;
  newChild[sibRow] = b;
  newAbsLen[sibRow] = (1.0 - tau) * lRegraft;

  // Canonical preorder reordering (topology + branch lengths in one pass)
  auto po = TreeTools::preorder_weighted_impl(newParent, newChild, newAbsLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector ordParent(nEdge), ordChild(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    ordParent[i] = ordEdge(i, 0);
    ordChild[i]  = ordEdge(i, 1);
  }
  NumericVector orderedRelBr = ordAbs / treeLength;

  const double logHastings = std::log(lRegraft) - std::log(lMerge);

  return List::create(_["parent"] = ordParent,
                      _["child"] = ordChild,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = logHastings);
}


// Rcpp-exported wrapper (for R callers via ProposeSpr)
// [[Rcpp::export]]
List spr_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths) {
  int nEdge = edge.nrow();
  IntegerVector parent(nEdge), child(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    parent[i] = edge(i, 0);
    child[i] = edge(i, 1);
  }
  List result = spr_proposal_impl(parent, child, nTip, treeLength,
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


// ---------------------------------------------------------------------------
// BetaSimplex proposal
// ---------------------------------------------------------------------------

// Internal helper — OPP-5: modifies x in-place, avoids Rcpp::List allocation.
// Returns false if the proposal is degenerate (total <= 0); sets logHastings.
// Output params for O(1) rollback: outOther = index of second modified element,
// outOldIdx/outOldOther = original values before modification.
bool beta_simplex_impl(NumericVector& x, int index, double tuning,
                       double& logHastings, int& outOther,
                       double& outOldIdx, double& outOldOther) {
  const int n = x.size();
  if (n < 2) { logHastings = 0.0; outOther = index; return true; }

  int other = (int)(unif_rand() * (double)(n - 1));
  if (other >= index) ++other;
  if (other >= n) other = n - 1;
  if (other == index) other = (index + 1) % n;
  outOther = other;

  const double oldA = x[index];
  const double oldB = x[other];
  outOldIdx   = oldA;
  outOldOther = oldB;
  const double total = oldA + oldB;
  if (total <= 0.0) { logHastings = 0.0; return false; }

  const double oldF    = oldA / total;
  const double alpha   = oldF * tuning + 1.0;
  const double betaPar = (1.0 - oldF) * tuning + 1.0;
  const double newF    = R::rbeta(alpha, betaPar);

  x[index] = newF * total;
  x[other] = (1.0 - newF) * total;

  const double logFwd  = R::dbeta(newF, alpha, betaPar, 1);
  const double revAlpha = newF * tuning + 1.0;
  const double revBeta  = (1.0 - newF) * tuning + 1.0;
  logHastings = R::dbeta(oldF, revAlpha, revBeta, 1) - logFwd;
  return true;
}


// Rcpp-exported wrapper — keeps R-callable interface returning a List.
// [[Rcpp::export]]
List beta_simplex_proposal(NumericVector x, int index, double tuning) {
  const int n = x.size();
  if (n < 2) return List::create(_["value"] = x, _["logHastings"] = 0.0);

  NumericVector xNew = clone(x);
  double logHastings;
  int dummy; double d1, d2;
  if (!beta_simplex_impl(xNew, index, tuning, logHastings, dummy, d1, d2))
    return List::create(_["value"] = x, _["logHastings"] = 0.0);
  return List::create(_["value"] = xNew, _["logHastings"] = logHastings);
}
