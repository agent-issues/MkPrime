#ifndef MKPRIME_FITCH_H
#define MKPRIME_FITCH_H

// Fitch parsimony scoring for Mk characters (M-119).
//
// Uses bitmask state sets: state i → bit (1 << i).
// Missing data (-1) → all bits set.  Observed states < 32 (uint32_t).
//
// The tree is stored as (parent, child) vectors in any row order: the postorder
// is derived from the topology, not from the array, because
// fitch_score_candidates rewrites rows in place and the result is no longer in
// preorder.  Row order is immaterial up to degree 3; the sequential fold below
// is order-sensitive at a polytomy of degree 4 or more, which MkPrime's binary
// and trifurcating-root encodings never present.

#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <cstdint>

// Edge indices ordered so that every edge out of a node precedes the edge
// into it — a postorder, whatever order the rows are stored in.
inline void fitch_postorder(
    const int* parent, const int* child, int nEdge, int nNode,
    std::vector<int>& order)
{
  std::vector<int> firstEdge(nNode + 1, -1), nextEdge(nEdge, -1);
  std::vector<char> hasParent(nNode + 1, 0);
  // Descending, so each node's outgoing edges stay in ascending row order.
  for (int e = nEdge - 1; e >= 0; --e) {
    const int p = parent[e], ch = child[e];
    if (p < 1 || p > nNode || ch < 1 || ch > nNode)
      Rcpp::stop("fitch_postorder: node ID out of range");
    if (hasParent[ch])
      Rcpp::stop("fitch_postorder: node %d has two parents", ch);
    nextEdge[e] = firstEdge[p];
    firstEdge[p] = e;
    hasParent[ch] = 1;
  }

  int root = -1;
  for (int e = 0; e < nEdge; ++e)
    if (!hasParent[parent[e]]) { root = parent[e]; break; }
  if (root < 0)
    Rcpp::stop("fitch_postorder: edge array has no root");

  order.clear();
  order.reserve(nEdge);
  std::vector<int> stack;
  stack.reserve(nEdge);
  stack.push_back(root);
  while (!stack.empty()) {
    const int nd = stack.back();
    stack.pop_back();
    for (int e = firstEdge[nd]; e >= 0; e = nextEdge[e]) {
      order.push_back(e);
      stack.push_back(child[e]);
    }
  }
  // One parent per node means no cycle is reachable from the root, so a short
  // traversal is a forest or an unreachable cycle, never a hang.
  if ((int)order.size() != nEdge)
    Rcpp::stop("fitch_postorder: edge array is not connected");

  std::reverse(order.begin(), order.end());
}

// Compute total Fitch parsimony score across all partitions.
// parent/child: 1-indexed node IDs, length nEdge, any row order.
// parts: vector of {tipStates (nTip × nChar), kStates} pairs.
// nTip: number of tips (nodes 1..nTip).
//
// Returns sum of character step counts across all partitions.
inline int fitch_score_all(
    const int* parent, const int* child, int nEdge, int nTip,
    const std::vector<std::pair<Rcpp::IntegerMatrix, int>>& parts)
{
  // Every buffer below is indexed by node ID, and nTip is caller-supplied, so
  // the bound comes from the IDs actually present.
  int nNode = 2 * nTip;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > nNode) nNode = parent[e];
    if (child[e] > nNode) nNode = child[e];
  }
  int totalScore = 0;

  std::vector<int> order;
  fitch_postorder(parent, child, nEdge, nNode, order);

  for (const auto& pp : parts) {
    const Rcpp::IntegerMatrix& tips = pp.first;
    const int nChar = tips.ncol();
    // Missing data gets every bit, not just the first k: k may exceed 32
    // (known-k partitions pass k verbatim), and bits no tip carries cannot
    // change the score, because every state set is then either All or a
    // subset of the observed states.
    const uint32_t allBits = ~0u;

    // State sets: indexed [nodeID * nChar + charIdx]
    // Node IDs are 1-based, so allocate nNode+1 slots.
    std::vector<uint32_t> ss((nNode + 1) * nChar, 0u);

    // Initialize tips (1-indexed)
    for (int t = 0; t < nTip; ++t) {
      const int nodeID = t + 1;
      for (int c = 0; c < nChar; ++c) {
        int s = tips(t, c);  // 0-based state, or -1 for missing
        if (s >= 32) Rcpp::stop("Fitch scoring supports at most 32 states");
        ss[nodeID * nChar + c] = (s < 0) ? allBits : (1u << s);
      }
    }

    // Track how many children we've merged into each internal node
    std::vector<int> nSeen(nNode + 1, 0);

    int score = 0;

    for (int oi = 0; oi < nEdge; ++oi) {
      const int e = order[oi];
      const int p = parent[e];
      const int ch = child[e];

      if (nSeen[p] == 0) {
        // First child: copy state sets
        for (int c = 0; c < nChar; ++c)
          ss[p * nChar + c] = ss[ch * nChar + c];
      } else {
        // Fitch intersection/union
        for (int c = 0; c < nChar; ++c) {
          uint32_t inter = ss[p * nChar + c] & ss[ch * nChar + c];
          if (inter) {
            ss[p * nChar + c] = inter;
          } else {
            ss[p * nChar + c] |= ss[ch * nChar + c];
            ++score;
          }
        }
      }
      ++nSeen[p];
    }

    totalScore += score;
  }

  return totalScore;
}


// Compute Fitch scores for every candidate regraft position in a single pass.
//
// After pruning subtree v (detaching edge u→v), the residual tree R has
// nEdge-2 effective edges (the pruneRow and sibRow are repurposed).
// For each candidate regraft edge, we:
//   1. Apply the SPR in-place (suppress u, insert u on candidate edge)
//   2. Score via fitch_score_all
//   3. Restore the original topology
//
// This is O(nCand × nEdge × nChar).
//
// parent/child are MODIFIED IN PLACE and restored on normal return; a throw
// from the scorer leaves them mid-regraft, so callers pass a copy.
inline void fitch_score_candidates(
    int* parent, int* child, int nEdge, int nTip,
    const std::vector<std::pair<Rcpp::IntegerMatrix, int>>& parts,
    int pruneRow, int parentRow, int sibRow,
    int u, int v, int sibNode,
    const std::vector<int>& candidates,
    std::vector<int>& scores)
{
  const int nCand = (int)candidates.size();
  scores.resize(nCand);

  // Save original topology at the 3 modified rows
  const int origParent_pr = parent[pruneRow], origChild_pr = child[pruneRow];
  const int origParent_pa = parent[parentRow], origChild_pa = child[parentRow];
  const int origParent_si = parent[sibRow], origChild_si = child[sibRow];

  // Step 1: suppress u  (parentRow: p→u becomes p→sibNode)
  child[parentRow] = sibNode;
  // pruneRow (u→v) and sibRow (u→w) get repurposed below per candidate

  for (int ci = 0; ci < nCand; ++ci) {
    const int rr = candidates[ci];
    const int b  = child[rr];

    // Step 2: insert u on candidate edge (rr: a→b becomes a→u)
    child[rr]      = u;
    // sibRow becomes u→b
    parent[sibRow]  = u;
    child[sibRow]   = b;
    // pruneRow stays as u→v (unchanged from original)
    parent[pruneRow] = u;
    child[pruneRow]  = v;

    scores[ci] = fitch_score_all(parent, child, nEdge, nTip, parts);

    // Restore candidate edge
    child[rr] = b;
  }

  // Restore all 3 rows to original
  parent[pruneRow] = origParent_pr; child[pruneRow] = origChild_pr;
  parent[parentRow] = origParent_pa; child[parentRow] = origChild_pa;
  parent[sibRow] = origParent_si; child[sibRow] = origChild_si;
}

#endif  // MKPRIME_FITCH_H
