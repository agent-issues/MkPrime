// C++ implementations of MCMC proposals (SPR, BetaSimplex).
//
// NNI is in tree_moves.cpp (nni_proposal).
// SPR: prune-and-regraft with Jacobian Hastings ratio.
// BetaSimplex: redistribute mass between two simplex elements.

#include <Rcpp.h>
#include <TreeTools/renumber_tree.h>
#include <vector>
#include <cmath>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// SPR proposal (Subtree Pruning and Regrafting)
// ---------------------------------------------------------------------------
//
// Inputs (all in ape edge-matrix convention, 1-based node indices):
//   edge         nEdge x 2 integer matrix (parent | child)
//   nTip         number of tips
//   treeLength   current total tree length
//   relBrLengths nEdge-length vector of relative branch lengths (sum = 1)
//
// Returns a list with:
//   edge            new edge matrix in postorder
//   rel_br_lengths  new relative branch lengths
//   logHastings     log(lRegraft) - log(lMerge), or -Inf if inapplicable

// [[Rcpp::export]]
List spr_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths) {
  const int nEdge = edge.nrow();
  const int root = nTip + 1;

  // Eligible prune edges: parent != root
  std::vector<int> eligiblePrune;
  eligiblePrune.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (edge(i, 0) != root) {
      eligiblePrune.push_back(i);
    }
  }

  if (eligiblePrune.empty()) {
    return List::create(_["edge"] = edge,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  int pickPrune = (int)(unif_rand() * (double)eligiblePrune.size());
  if (pickPrune >= (int)eligiblePrune.size()) pickPrune = eligiblePrune.size() - 1;
  const int pruneRow = eligiblePrune[pickPrune];
  const int u = edge(pruneRow, 0);  // node being detached (parent of pruned)
  const int v = edge(pruneRow, 1);  // pruned subtree root

  // u's parent: row where child == u
  int parentRow = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (edge(i, 1) == u) {
      parentRow = i;
      break;
    }
  }

  // v's sibling: row where parent == u, child != v
  int sibRow = -1;
  int w = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (edge(i, 0) == u && edge(i, 1) != v) {
      sibRow = i;
      w = edge(i, 1);
      break;
    }
  }

  if (parentRow < 0 || sibRow < 0) {
    return List::create(_["edge"] = edge,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  // BFS to find all descendants of v (for exclusion)
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
        if (edge(i, 0) == cur) {
          isDesc[edge(i, 1)] = true;
          if (edge(i, 1) > nTip) {
            queue.push_back(edge(i, 1));
          }
        }
      }
    }
  }

  // Candidate regraft edges: not in v's subtree, not adjacent to u
  std::vector<int> candidates;
  candidates.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[edge(i, 1)]) continue;
    if (edge(i, 0) == u || edge(i, 1) == u) continue;
    candidates.push_back(i);
  }

  if (candidates.empty()) {
    return List::create(_["edge"] = edge,
                        _["rel_br_lengths"] = relBrLengths,
                        _["logHastings"] = R_NegInf);
  }

  int pickRegraft = (int)(unif_rand() * (double)candidates.size());
  if (pickRegraft >= (int)candidates.size()) pickRegraft = candidates.size() - 1;
  const int regraftRow = candidates[pickRegraft];
  const int b = edge(regraftRow, 1);  // save original child of regraft edge

  const double tau = unif_rand();

  // Absolute branch lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    absLen[i] = treeLength * relBrLengths[i];
  }

  const double lRegraft = absLen[regraftRow];
  const double lMerge = absLen[parentRow] + absLen[sibRow];

  // Perform SPR on cloned edge matrix and branch lengths
  IntegerMatrix newEdge = clone(edge);
  NumericVector newAbsLen = clone(absLen);

  // 1. Suppress u: (p -> u) becomes (p -> w)
  newEdge(parentRow, 1) = w;
  newAbsLen[parentRow] = lMerge;

  // 2. Insert u on regraft edge: (a -> b) becomes (a -> u)
  newEdge(regraftRow, 1) = u;
  newAbsLen[regraftRow] = tau * lRegraft;

  // 3. Reuse sibRow for (u -> b)
  newEdge(sibRow, 0) = u;
  newEdge(sibRow, 1) = b;
  newAbsLen[sibRow] = (1.0 - tau) * lRegraft;

  // Reorder to postorder (no node renumbering)
  IntegerVector order = TreeTools::postorder_order(newEdge);

  IntegerMatrix orderedEdge(nEdge, 2);
  NumericVector orderedRelBr(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    int j = order[i] - 1;  // 1-based -> 0-based
    orderedEdge(i, 0) = newEdge(j, 0);
    orderedEdge(i, 1) = newEdge(j, 1);
    orderedRelBr[i] = newAbsLen[j] / treeLength;
  }

  const double logHastings = std::log(lRegraft) - std::log(lMerge);

  return List::create(_["edge"] = orderedEdge,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = logHastings);
}


// ---------------------------------------------------------------------------
// BetaSimplex proposal
// ---------------------------------------------------------------------------
//
// Picks element `index` (0-based) and a random other element, redistributes
// mass between them via a Beta draw. Other elements unchanged.
//
// Returns list(value, logHastings).

// [[Rcpp::export]]
List beta_simplex_proposal(NumericVector x, int index, double tuning) {
  const int n = x.size();
  if (n < 2) {
    return List::create(_["value"] = x, _["logHastings"] = 0.0);
  }

  // Pick random other element (index is 0-based)
  int other = (int)(unif_rand() * (double)(n - 1));
  if (other >= index) ++other;
  // Clamp (unif_rand can return exactly 1.0 on some platforms)
  if (other >= n) other = n - 1;
  if (other == index) other = (index + 1) % n;

  const double oldA = x[index];
  const double oldB = x[other];
  const double total = oldA + oldB;

  if (total <= 0.0) {
    return List::create(_["value"] = x, _["logHastings"] = 0.0);
  }

  const double oldF = oldA / total;
  const double alpha = oldF * tuning + 1.0;
  const double betaPar = (1.0 - oldF) * tuning + 1.0;
  const double newF = R::rbeta(alpha, betaPar);

  NumericVector xNew = clone(x);
  xNew[index] = newF * total;
  xNew[other] = (1.0 - newF) * total;

  // Hastings ratio: log q(old|new) - log q(new|old)
  const double logFwd = R::dbeta(newF, alpha, betaPar, 1);
  const double revAlpha = newF * tuning + 1.0;
  const double revBeta = (1.0 - newF) * tuning + 1.0;
  const double logRev = R::dbeta(oldF, revAlpha, revBeta, 1);

  return List::create(_["value"] = xNew,
                      _["logHastings"] = logRev - logFwd);
}
