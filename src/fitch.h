#ifndef MKPRIME_FITCH_H
#define MKPRIME_FITCH_H

// Fitch parsimony scoring for Mk characters (M-119).
//
// Uses bitmask state sets: state i → bit (1 << i).
// Missing data (-1) → all bits set.  k ≤ 26 (fits in uint32_t).
//
// The tree is stored as (parent, child) vectors in preorder.
// Processing edges in reverse gives a valid postorder traversal.

#include <Rcpp.h>
#include <vector>
#include <cstdint>

// Compute total Fitch parsimony score across all partitions.
// parent/child: 1-indexed node IDs in preorder, length nEdge.
// parts: vector of {tipStates (nTip × nChar), kStates} pairs.
// nTip: number of tips (nodes 1..nTip).
//
// Returns sum of character step counts across all partitions.
inline int fitch_score_all(
    const int* parent, const int* child, int nEdge, int nTip,
    const std::vector<std::pair<Rcpp::IntegerMatrix, int>>& parts)
{
  const int nNode = 2 * nTip;  // nodes numbered 1..nNode for unrooted binary
  int totalScore = 0;

  for (const auto& pp : parts) {
    const Rcpp::IntegerMatrix& tips = pp.first;
    const int kStates = pp.second;
    const int nChar = tips.ncol();
    const uint32_t allBits = (1u << kStates) - 1;

    // State sets: indexed [nodeID * nChar + charIdx]
    // Node IDs are 1-based, so allocate nNode+1 slots.
    std::vector<uint32_t> ss((nNode + 1) * nChar, 0u);

    // Initialize tips (1-indexed)
    for (int t = 0; t < nTip; ++t) {
      const int nodeID = t + 1;
      for (int c = 0; c < nChar; ++c) {
        int s = tips(t, c);  // 0-based state, or -1 for missing
        ss[nodeID * nChar + c] = (s < 0) ? allBits : (1u << s);
      }
    }

    // Track how many children we've merged into each internal node
    std::vector<int> nSeen(nNode + 1, 0);

    int score = 0;

    // Reverse-preorder = postorder: children before parents
    for (int e = nEdge - 1; e >= 0; --e) {
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
// This is O(nCand × nEdge × nChar) — plenty fast for morphological data
// (≤500 chars, ≤300 tips, ≤100 candidates → ~15M bitwise ops, <5 ms).
//
// parent/child are MODIFIED IN PLACE and restored; caller's arrays unchanged
// on return.
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
