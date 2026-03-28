#ifndef MKPRIME_GIBBS_PARTIAL_CL_H
#define MKPRIME_GIBBS_PARTIAL_CL_H

// M-105: Partial CL (conditional likelihood) reuse for Gibbs SPR
//
// Instead of running a full-tree Felsenstein pruning for each of the ~100
// candidate SPR regraft positions, we:
//   1. Run ONE full downpass, caching per-node CLs and per-child contributions.
//   2. Compute the "residual tree" CLs (tree with the pruned subtree removed).
//   3. For each candidate regraft, propagate updates only along the
//      regraft-to-root path (~6 nodes for a 54-tip tree instead of ~107).
//
// Complexity: O(nEdge) setup + O(depth × nCand) evaluation
// vs. current: O(nEdge × nCand)
// For nCand=100, depth=6, nEdge=105: 20× reduction.

#include "mcmc_state.h"
#include <cstring>
#include <array>
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// TreeNav: parent–child navigation built from canonical preorder edge table
// ---------------------------------------------------------------------------
struct TreeNav {
  int nTip, nEdge, root, maxNode;

  // Per-node (1-indexed; index 0 unused):
  // ch0, ch1, ch2: children (-1 for tips / missing).
  // Unrooted tree root has 3 children (trifurcation); all others have 2.
  std::vector<int> ch0, ch1, ch2;
  std::vector<int> parentNode;    // parent node (-1 for root)
  std::vector<int> edgeToPar;     // edge index connecting node to its parent
  std::vector<double> edgeLen;    // absolute edge length for each edge

  void build(const IntegerVector& parent, const IntegerVector& child,
             const NumericVector& absEdgeLen, int nTip_ = 0) {
    nEdge = parent.size();
    maxNode = 0;
    for (int i = 0; i < nEdge; ++i) {
      if (parent[i] > maxNode) maxNode = parent[i];
      if (child[i]  > maxNode) maxNode = child[i];
    }
    // nTip must be passed explicitly for unrooted trees (trifurcating root)
    // where nEdge = 2*nTip-3 rather than 2*nTip-2.
    nTip = (nTip_ > 0) ? nTip_ : (nEdge + 2) / 2;
    root = nTip + 1;

    ch0.assign(maxNode + 1, -1);
    ch1.assign(maxNode + 1, -1);
    ch2.assign(maxNode + 1, -1);
    parentNode.assign(maxNode + 1, -1);
    edgeToPar.assign(maxNode + 1, -1);
    edgeLen.assign(nEdge, 0.0);

    for (int e = 0; e < nEdge; ++e) edgeLen[e] = absEdgeLen[e];

    // Fill parent/child maps.  Postorder (reverse) determines which child
    // is ch0 (first encountered), ch1 (second), ch2 (third = root only).
    for (int e = nEdge - 1; e >= 0; --e) {
      int p = parent[e], c = child[e];
      parentNode[c] = p;
      edgeToPar[c]  = e;
      if      (ch0[p] < 0) ch0[p] = c;
      else if (ch1[p] < 0) ch1[p] = c;
      else                  ch2[p] = c;  // root's 3rd child
    }
  }

  // Which slot (0, 1, or 2) does child `c` occupy at parent `p`?
  int childSlot(int p, int c) const {
    if (ch0[p] == c) return 0;
    if (ch1[p] == c) return 1;
    return 2;
  }
};


// ---------------------------------------------------------------------------
// CLGroup: cached conditional likelihoods for one (partition, kStates) group
// ---------------------------------------------------------------------------
struct CLGroup {
  int kStates, nChar, nCat, maxNode;
  int stride;       // nChar * kStates
  bool isMkN;       // true → asymmetric binary model (uses rateLoss)
  double rateLoss;  // only for isMkN
  double rateScale; // for neo: rateNeo; for others: 1.0

  IntegerMatrix tipData;  // nTip × nChar, 0-indexed states, -1 missing

  // Flat storage indexed as [(cat * (maxNode+1) + node) * stride + c*kStates+s]
  std::vector<double> inside;  // I[n]  = CL at node (product of children)
  std::vector<double> from0;   // F0[n] = contribution from ch0 at node n
  std::vector<double> from1;   // F1[n] = contribution from ch1 at node n
  std::vector<double> from2;   // F2[n] = contribution from ch2 (root only)

  double* I (int cat, int node) { return inside.data() + ((size_t)cat * (maxNode+1) + node) * stride; }
  double* F0(int cat, int node) { return from0.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }
  double* F1(int cat, int node) { return from1.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }
  double* F2(int cat, int node) { return from2.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }
  const double* I (int cat, int node) const { return inside.data() + ((size_t)cat * (maxNode+1) + node) * stride; }
  const double* F0(int cat, int node) const { return from0.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }
  const double* F1(int cat, int node) const { return from1.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }
  const double* F2(int cat, int node) const { return from2.data()  + ((size_t)cat * (maxNode+1) + node) * stride; }

  // Get F for a given slot (0, 1, or 2)
  double* Fslot(int cat, int node, int slot) {
    return (slot == 0) ? F0(cat, node) : (slot == 1) ? F1(cat, node) : F2(cat, node);
  }
  const double* Fslot(int cat, int node, int slot) const {
    return (slot == 0) ? F0(cat, node) : (slot == 1) ? F1(cat, node) : F2(cat, node);
  }

  // Product of all F-slots EXCEPT the one for `excludeSlot` at node
  // result[i] = product of F_j[i] for j != excludeSlot
  void productExcluding(int cat, int node, int excludeSlot,
                        double* result, int nSlots) const {
    std::fill(result, result + stride, 1.0);
    for (int s = 0; s < nSlots; ++s) {
      if (s == excludeSlot) continue;
      const double* f = Fslot(cat, node, s);
      for (int i = 0; i < stride; ++i) result[i] *= f[i];
    }
  }

  void allocate(int maxNode_, int nCat_, int nChar_, int kStates_) {
    maxNode = maxNode_; nCat = nCat_; nChar = nChar_; kStates = kStates_;
    stride = nChar * kStates;
    size_t total = (size_t)nCat * (maxNode + 1) * stride;
    inside.assign(total, 0.0);
    from0.assign(total, 0.0);
    from1.assign(total, 0.0);
    from2.assign(total, 0.0);
  }
};


// ---------------------------------------------------------------------------
// GibbsCLCache: the full cache for one Gibbs SPR/SubtreeSwap evaluation
// ---------------------------------------------------------------------------
struct GibbsCLCache {
  TreeNav topo;
  std::vector<CLGroup> groups;
  NumericVector rates;  // ACRV rate multipliers
  bool useAcrv;
  int coding;  // 0=none, 1=variable, 2=informative

  // Constant-site probability product per state, per group (for ascertainment)
  // constProd[grp][cat * kStates + s] = π_s × Π_edges P_ss(t_e * rates[cat])
  std::vector<std::vector<double>> constProd;

  // Temporary per-node buffers used during path propagation (avoids alloc)
  // Sized to max stride across all groups × nCat
  std::vector<double> tmpNode;

  bool ready = false;
};


// ---------------------------------------------------------------------------
// ACRV rate computation (matches cpp_acrv_rates in mcmc_likelihood.cpp)
// ---------------------------------------------------------------------------
inline NumericVector gibbs_acrv_rates(double rateLogSd, int nCat,
                                       const std::vector<double>& acrvZ) {
  if (rateLogSd <= 0.0) return NumericVector(nCat, 1.0);
  double mu = -rateLogSd * rateLogSd / 2.0;
  NumericVector rates(nCat);
  double total = 0.0;
  for (int i = 0; i < nCat; ++i) {
    rates[i] = std::exp(mu + rateLogSd * acrvZ[i]);
    total += rates[i];
  }
  for (int i = 0; i < nCat; ++i) rates[i] *= nCat / total;
  return rates;
}


// ---------------------------------------------------------------------------
// JC transition helpers
// ---------------------------------------------------------------------------

// Compute p_diff and diff_coeff for JC(k) with branch length t
inline void jc_trans_params(int k, double t,
                            double& p_diff, double& diff_coeff) {
  double inv_k = 1.0 / k;
  double exp_term = std::exp(-k * t / (k - 1.0));
  p_diff     = inv_k - inv_k * exp_term;     // = (1/k)(1 - exp)
  diff_coeff = exp_term;                       // = p_same - p_diff
}

// Apply JC transition to child CL → result buf.
// result[c*k+i] = p_diff * sum(cl[c*k+:]) + diff_coeff * cl[c*k+i]
inline void jc_transition(const double* cl, double* result,
                           int nChar, int kStates,
                           double p_diff, double diff_coeff) {
  for (int c = 0; c < nChar; ++c) {
    int off = c * kStates;
    double s = 0.0;
    for (int j = 0; j < kStates; ++j) s += cl[off + j];
    for (int i = 0; i < kStates; ++i)
      result[off + i] = p_diff * s + diff_coeff * cl[off + i];
  }
}

// MkN (asymmetric binary k=2) transition helpers
inline void mkn_trans_params(double rateLoss, double t,
                              double& P00, double& P01,
                              double& P10, double& P11) {
  double sum_rl = 1.0 + rateLoss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rateLoss / sum_rl;
  double lambda = rate01 + rate10;
  double exp_t  = std::exp(-lambda * t);
  double i01    = rate01 / lambda;
  double i10    = rate10 / lambda;
  P00 = i10 + i01 * exp_t;
  P01 = i01 - i01 * exp_t;
  P10 = i10 - i10 * exp_t;
  P11 = i01 + i10 * exp_t;
}

inline void mkn_transition(const double* cl, double* result,
                            int nChar,
                            double P00, double P01, double P10, double P11) {
  for (int c = 0; c < nChar; ++c) {
    int off = c * 2;
    double c0 = cl[off], c1 = cl[off + 1];
    result[off]     = P00 * c0 + P01 * c1;
    result[off + 1] = P10 * c0 + P11 * c1;
  }
}


// ---------------------------------------------------------------------------
// Caching downpass: full postorder traversal, storing I, F0, F1 per node
// ---------------------------------------------------------------------------
static void caching_downpass(CLGroup& grp, const TreeNav& topo,
                              const IntegerVector& parent,
                              const IntegerVector& child,
                              const NumericVector& rates) {
  int nEdge   = parent.size();
  int nTip    = topo.nTip;
  int kStates = grp.kStates;
  int nChar   = grp.nChar;
  int stride  = grp.stride;
  int maxNode = grp.maxNode;
  int nCat    = grp.nCat;

  // Init tip CLs (same for all categories)
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = grp.I(0, tip);
    std::fill(cl, cl + stride, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state = grp.tipData(tip - 1, c);
      int off   = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[off + s] = 1.0;
      } else {
        cl[off + state] = 1.0;
      }
    }
    for (int cat = 1; cat < nCat; ++cat)
      std::memcpy(grp.I(cat, tip), cl, stride * sizeof(double));
  }

  // Temporary buffer for one child's contribution (stack for small stride)
  std::vector<double> contrib(stride);

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Track how many children seen at each node (root may have 3)
    std::vector<uint8_t> nChildSeen(maxNode + 1, 0);

    // Postorder: reverse edge order
    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch_ = child[e];
      double t = topo.edgeLen[e] * rate * grp.rateScale;

      // Compute transition: contrib = P(t) × I[ch_]
      if (grp.isMkN) {
        double P00, P01, P10, P11;
        mkn_trans_params(grp.rateLoss, t, P00, P01, P10, P11);
        mkn_transition(grp.I(cat, ch_), contrib.data(), nChar, P00, P01, P10, P11);
      } else {
        double p_diff, diff_coeff;
        jc_trans_params(kStates, t, p_diff, diff_coeff);
        jc_transition(grp.I(cat, ch_), contrib.data(), nChar, kStates, p_diff, diff_coeff);
      }

      // Track how many children have been seen (0, 1, 2 for root's 3rd)
      if (nChildSeen[par] == 0) {
        std::memcpy(grp.F0(cat, par), contrib.data(), stride * sizeof(double));
        std::memcpy(grp.I(cat, par),  contrib.data(), stride * sizeof(double));
      } else if (nChildSeen[par] == 1) {
        std::memcpy(grp.F1(cat, par), contrib.data(), stride * sizeof(double));
        double* clPar = grp.I(cat, par);
        for (int i = 0; i < stride; ++i) clPar[i] *= contrib[i];
      } else {
        // 3rd child (root only)
        std::memcpy(grp.F2(cat, par), contrib.data(), stride * sizeof(double));
        double* clPar = grp.I(cat, par);
        for (int i = 0; i < stride; ++i) clPar[i] *= contrib[i];
      }
      nChildSeen[par]++;
    }
  }
}


// ---------------------------------------------------------------------------
// Compute residual-tree CLs along the path from grandparent(u) to root.
//
// After pruning subtree v from parent u:
//   - u is collapsed; sibNode connects directly to g = parent(u) via lMerge
//   - CLs change along the path g → root
//
// Stores updated I values in residual[cat][pathNode * stride ...].
// The path is returned in `residPath` (g → ... → root).
// ---------------------------------------------------------------------------
struct ResidualCL {
  std::vector<int> path;             // nodes from g (inclusive) to root (inclusive)
  std::vector<double> data;          // flat: [cat * nPathNodes * stride + pathIdx * stride + ...]
  int nPathNodes;
  int stride;
  int nCat;

  // At each path node, one child slot's contribution (F) was replaced by the
  // residual computation.  Store the slot index and the new F values so that
  // evaluate_candidate can use them instead of the stale original F-slots.
  //
  // At path[0] (= g): the pruned child u's slot is replaced by sibNode's
  //   contribution via lMerge.
  // At path[pi] (pi>0): the path-child's slot is replaced by the recomputed
  //   contribution from the updated child below.
  std::vector<int> replacedSlot;         // [nPathNodes]
  std::vector<double> replacedF;         // flat: [nCat * nPathNodes * stride]

  double* I(int cat, int pathIdx) {
    return data.data() + ((size_t)cat * nPathNodes + pathIdx) * stride;
  }
  const double* I(int cat, int pathIdx) const {
    return data.data() + ((size_t)cat * nPathNodes + pathIdx) * stride;
  }

  double* F(int cat, int pathIdx) {
    return replacedF.data() + ((size_t)cat * nPathNodes + pathIdx) * stride;
  }
  const double* F(int cat, int pathIdx) const {
    return replacedF.data() + ((size_t)cat * nPathNodes + pathIdx) * stride;
  }

  // Look up path index for a node; -1 if not on path
  int pathIndex(int node) const {
    for (int i = 0; i < nPathNodes; ++i)
      if (path[i] == node) return i;
    return -1;
  }
};


static void compute_residual_cl(
    ResidualCL& res,
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int u, int sibNode, double lMerge)
{
  int g = topo.parentNode[u];   // grandparent of pruned subtree
  int kStates = grp.kStates;
  int nChar   = grp.nChar;
  int stride  = grp.stride;
  int nCat    = grp.nCat;

  // Build path g → root
  res.path.clear();
  for (int cur = g; cur >= 1; cur = topo.parentNode[cur])
    res.path.push_back(cur);
  res.nPathNodes = (int)res.path.size();
  res.stride     = stride;
  res.nCat       = nCat;
  res.data.resize((size_t)nCat * res.nPathNodes * stride);
  res.replacedSlot.resize(res.nPathNodes);
  res.replacedF.resize((size_t)nCat * res.nPathNodes * stride);

  std::vector<double> contrib(stride);

  // Record the replaced slot at g: u's slot is replaced by sibNode via lMerge
  int uSlot = topo.childSlot(g, u);
  res.replacedSlot[0] = uSlot;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // At g: replace u's contribution with sibNode's contribution via lMerge
    if (grp.isMkN) {
      double P00, P01, P10, P11;
      mkn_trans_params(grp.rateLoss, lMerge * rate * grp.rateScale, P00, P01, P10, P11);
      mkn_transition(grp.I(cat, sibNode), contrib.data(), nChar, P00, P01, P10, P11);
    } else {
      double p_diff, diff_coeff;
      jc_trans_params(kStates, lMerge * rate * grp.rateScale, p_diff, diff_coeff);
      jc_transition(grp.I(cat, sibNode), contrib.data(), nChar, kStates, p_diff, diff_coeff);
    }

    // Store the replacement F at g
    std::memcpy(res.F(cat, 0), contrib.data(), stride * sizeof(double));

    // I_R[g] = new_sibNode_contrib × product of all other children's contributions
    double* Ig = res.I(cat, 0);
    int nChildG = (topo.ch2[g] >= 0) ? 3 : 2;
    grp.productExcluding(cat, g, uSlot, Ig, nChildG);
    for (int i = 0; i < stride; ++i) Ig[i] *= contrib[i];

    // Propagate up: for each node on the path from g to root
    for (int pi = 1; pi < res.nPathNodes; ++pi) {
      int node   = res.path[pi];
      int prev   = res.path[pi - 1];  // child on the path
      int prevSlot = topo.childSlot(node, prev);

      // Record which slot is replaced at this path node
      if (cat == 0) res.replacedSlot[pi] = prevSlot;

      // Contribution from prev (updated) through edge prev→node
      double t = topo.edgeLen[topo.edgeToPar[prev]] * rate * grp.rateScale;
      const double* prevI = res.I(cat, pi - 1);  // updated CL of prev
      if (grp.isMkN) {
        double P00, P01, P10, P11;
        mkn_trans_params(grp.rateLoss, t, P00, P01, P10, P11);
        mkn_transition(prevI, contrib.data(), nChar, P00, P01, P10, P11);
      } else {
        double p_diff, diff_coeff;
        jc_trans_params(kStates, t, p_diff, diff_coeff);
        jc_transition(prevI, contrib.data(), nChar, kStates, p_diff, diff_coeff);
      }

      // Store the replacement F at this path node
      std::memcpy(res.F(cat, pi), contrib.data(), stride * sizeof(double));

      // Combine with product of all unchanged sibling contributions
      double* Inode = res.I(cat, pi);
      int nChild = (topo.ch2[node] >= 0) ? 3 : 2;
      grp.productExcluding(cat, node, prevSlot, Inode, nChild);
      for (int i = 0; i < stride; ++i) Inode[i] *= contrib[i];
    }
  }
}


// ---------------------------------------------------------------------------
// Evaluate one candidate regraft.
//
// Given residual tree (v detached from u) and candidate regraft at edge (a→b)
// with length lReg split at tau=0.5, compute the per-category per-site
// likelihoods.  Returns the log-likelihood contribution for this CLGroup.
//
// The key complexity: the original tree's parent-child navigation differs from
// the SPR-modified tree.  In particular, the path from a to root in the
// original tree may pass through u (the pruned node), but in the new tree u
// is at the regraft point and the path skips from sibNode to g via lMerge.
// We handle this by detecting parNode==u during propagation and jumping to g.
// ---------------------------------------------------------------------------
static double evaluate_candidate(
    const CLGroup& grp, const TreeNav& topo,
    const ResidualCL& res,
    const NumericVector& rates,
    int v, int u,           // pruned subtree root and its parent
    int sibNode,            // u's other child (sibling of v)
    double lMerge,          // collapsed edge length sibNode→g
    int a, int b,           // regraft parent and child nodes
    double lHalfReg,        // lReg / 2 (half the regraft edge length)
    double lPrune)          // prune edge length (u→v)
{
  int kStates = grp.kStates;
  int nChar   = grp.nChar;
  int stride  = grp.stride;
  int nCat    = grp.nCat;
  double sc   = grp.rateScale;
  int g       = topo.parentNode[u];
  int uSlot   = topo.childSlot(g, u);

  // Root frequency: 1/k for JC, stationary for MkN
  double rootF0 = 0.5, rootF1 = 0.5;  // MkN
  double inv_k = 1.0 / kStates;       // JC
  if (grp.isMkN) {
    rootF0 = 1.0 / (1.0 + grp.rateLoss);
    rootF1 = grp.rateLoss / (1.0 + grp.rateLoss);
  }

  // Per-site likelihood sum across ACRV categories
  std::vector<double> siteLikSum(nChar, 0.0);

  // Buffers for per-character computations
  std::vector<double> from_v(stride), from_b(stride), Iu(stride), contrib(stride);
  std::vector<double> curI(stride);

  // Helper: apply transition matrix to src, write to dst
  auto applyTransition = [&](const double* src, double* dst, double t) {
    if (grp.isMkN) {
      double P00, P01, P10, P11;
      mkn_trans_params(grp.rateLoss, t, P00, P01, P10, P11);
      mkn_transition(src, dst, nChar, P00, P01, P10, P11);
    } else {
      double pd, dc;
      jc_trans_params(kStates, t, pd, dc);
      jc_transition(src, dst, nChar, kStates, pd, dc);
    }
  };

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Contribution from v (pruned subtree) at u
    applyTransition(grp.I(cat, v), from_v.data(), lPrune * rate * sc);

    // Contribution from b at u (via half regraft edge)
    int bResIdx = res.pathIndex(b);
    const double* Ib = (bResIdx >= 0) ? res.I(cat, bResIdx) : grp.I(cat, b);
    applyTransition(Ib, from_b.data(), lHalfReg * rate * sc);

    // I at u = from_v × from_b
    for (int i = 0; i < stride; ++i) Iu[i] = from_v[i] * from_b[i];

    // from_u at a: transition through half regraft edge from u to a
    applyTransition(Iu.data(), contrib.data(), lHalfReg * rate * sc);

    // Helper lambda: compute the product of all sibling contributions at node
    // 'par' EXCEPT the child in slot 'excludeSlot'.
    //
    // Three cases for each sibling slot:
    //   1. par is on the residual path and slot == replacedSlot → use res.F
    //   2. Sibling node is on the residual path → recompute P(t) × res.I
    //   3. Otherwise → use original cached F
    auto siblingProduct = [&](int par, int excludeSlot, double* out) {
      int nChild = (topo.ch2[par] >= 0) ? 3 : 2;
      std::fill(out, out + stride, 1.0);
      int parResIdx = res.pathIndex(par);
      for (int sl = 0; sl < nChild; ++sl) {
        if (sl == excludeSlot) continue;

        // Case 1: this slot was replaced during residual computation
        if (parResIdx >= 0 && sl == res.replacedSlot[parResIdx]) {
          const double* f = res.F(cat, parResIdx);
          for (int i = 0; i < stride; ++i) out[i] *= f[i];
          continue;
        }

        // Case 2: sibling node is on the residual path
        int cn = (sl == 0) ? topo.ch0[par] :
                 (sl == 1) ? topo.ch1[par] : topo.ch2[par];
        int cnResIdx = res.pathIndex(cn);
        if (cnResIdx >= 0) {
          double t_sib = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
          const double* sibI = res.I(cat, cnResIdx);
          // Use from_v as scratch buffer for the transition
          applyTransition(sibI, from_v.data(), t_sib);
          for (int i = 0; i < stride; ++i) out[i] *= from_v[i];
          continue;
        }

        // Case 3: use original cached F
        const double* f = grp.Fslot(cat, par, sl);
        for (int i = 0; i < stride; ++i) out[i] *= f[i];
      }
    };

    // At a: I_cand[a] = from_u_at_a × product_of_siblings(a, exclude b)
    int bSlot = topo.childSlot(a, b);
    siblingProduct(a, bSlot, curI.data());
    for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];

    // Propagate up from a to root.
    // Key: if we encounter u as a parent, skip it — u doesn't exist on
    // this path in the new topology.  Instead connect to g via lMerge.
    int curNode = a;
    int parNode = topo.parentNode[a];
    while (parNode >= 1) {
      if (parNode == u) {
        // Skip u: in the new topology, curNode connects to g via lMerge
        // (curNode must be sibNode or an ancestor of sibNode that leads to u)
        applyTransition(curI.data(), contrib.data(), lMerge * rate * sc);

        // At g: exclude u's slot (which doesn't exist in the new tree)
        siblingProduct(g, uSlot, curI.data());
        for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];

        curNode = g;
        parNode = topo.parentNode[g];
        continue;
      }

      int curSlot = topo.childSlot(parNode, curNode);

      // Transition from curNode to parNode
      double t = topo.edgeLen[topo.edgeToPar[curNode]] * rate * sc;
      applyTransition(curI.data(), contrib.data(), t);

      // I_cand[parNode] = from_cur × product_of_siblings(parNode, exclude cur)
      siblingProduct(parNode, curSlot, curI.data());
      for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];

      curNode = parNode;
      parNode = topo.parentNode[parNode];
    }

    // curI now holds I_cand[root].  Accumulate per-site likelihood.
    for (int c = 0; c < nChar; ++c) {
      double sl;
      if (grp.isMkN) {
        int off = c * 2;
        sl = rootF0 * curI[off] + rootF1 * curI[off + 1];
      } else {
        int off = c * kStates;
        sl = 0.0;
        for (int s = 0; s < kStates; ++s) sl += inv_k * curI[off + s];
      }
      siteLikSum[c] += sl;
    }
  }

  // Log-likelihood: sum_c log(siteLikSum[c] / nCat)
  double logLik = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = siteLikSum[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// ---------------------------------------------------------------------------
// Ascertainment correction via partial CL.
//
// P(constant site) = (1/nCat) × sum_cat sum_s π_s × CL_root(s, cat)
//
// where CL_root(s, cat) is the Felsenstein conditional likelihood for a
// constant-site pattern (all tips in state s).  This requires full pruning
// on k pseudo-characters, not a product of diagonal P_ss elements.
//
// Uses the partial CL framework: one caching downpass on the pseudo-chars,
// one residual computation, then O(depth × k²) per candidate evaluation.
// ---------------------------------------------------------------------------

// Create a pseudo-character CLGroup for constant-site patterns.
// Returns a CLGroup with kStates pseudo-characters, each with all tips in
// state s.  Same model parameters as the source group.
static CLGroup create_const_pseudo_group(const CLGroup& src, int nTip, int maxNode, int nCat) {
  CLGroup pg;
  pg.isMkN      = src.isMkN;
  pg.rateLoss   = src.rateLoss;
  pg.rateScale  = src.rateScale;
  int k = src.kStates;

  // Tip data: k pseudo-characters, each with all tips in state s
  IntegerMatrix tips(nTip, k);
  for (int c = 0; c < k; ++c)
    for (int t = 0; t < nTip; ++t)
      tips(t, c) = c;
  pg.tipData = tips;
  pg.allocate(maxNode, nCat, k, k);
  return pg;
}

// Evaluate P_const for a candidate regraft using partial CLs on pseudo-chars.
// Same algorithm as evaluate_candidate but returns P(constant site) instead
// of log-likelihood.
static double evaluate_const_prob(
    const CLGroup& pg,          // pseudo-character group
    const TreeNav& topo,
    const ResidualCL& res,
    const NumericVector& rates,
    int v, int u, int sibNode, double lMerge,
    int a, int b, double lHalfReg, double lPrune)
{
  int kStates = pg.kStates;
  int nChar   = pg.nChar;  // = kStates
  int stride  = pg.stride;
  int nCat    = pg.nCat;
  double sc   = pg.rateScale;
  int g       = topo.parentNode[u];
  int uSlot   = topo.childSlot(g, u);

  double inv_k = 1.0 / kStates;
  double rootF0 = 0.5, rootF1 = 0.5;
  if (pg.isMkN) {
    rootF0 = 1.0 / (1.0 + pg.rateLoss);
    rootF1 = pg.rateLoss / (1.0 + pg.rateLoss);
  }

  std::vector<double> from_v(stride), from_b(stride), Iu(stride), contrib(stride);
  std::vector<double> curI(stride);

  auto applyTransition = [&](const double* src, double* dst, double t) {
    if (pg.isMkN) {
      double P00, P01, P10, P11;
      mkn_trans_params(pg.rateLoss, t, P00, P01, P10, P11);
      mkn_transition(src, dst, nChar, P00, P01, P10, P11);
    } else {
      double pd, dc;
      jc_trans_params(kStates, t, pd, dc);
      jc_transition(src, dst, nChar, kStates, pd, dc);
    }
  };

  double totalConstProb = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    applyTransition(pg.I(cat, v), from_v.data(), lPrune * rate * sc);

    int bResIdx = res.pathIndex(b);
    const double* Ib = (bResIdx >= 0) ? res.I(cat, bResIdx) : pg.I(cat, b);
    applyTransition(Ib, from_b.data(), lHalfReg * rate * sc);

    for (int i = 0; i < stride; ++i) Iu[i] = from_v[i] * from_b[i];

    applyTransition(Iu.data(), contrib.data(), lHalfReg * rate * sc);

    // siblingProduct lambda (same as evaluate_candidate)
    auto siblingProduct = [&](int par, int excludeSlot, double* out) {
      int nChild = (topo.ch2[par] >= 0) ? 3 : 2;
      std::fill(out, out + stride, 1.0);
      int parResIdx = res.pathIndex(par);
      for (int sl = 0; sl < nChild; ++sl) {
        if (sl == excludeSlot) continue;
        if (parResIdx >= 0 && sl == res.replacedSlot[parResIdx]) {
          const double* f = res.F(cat, parResIdx);
          for (int i = 0; i < stride; ++i) out[i] *= f[i];
          continue;
        }
        int cn = (sl == 0) ? topo.ch0[par] :
                 (sl == 1) ? topo.ch1[par] : topo.ch2[par];
        int cnResIdx = res.pathIndex(cn);
        if (cnResIdx >= 0) {
          double t_sib = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
          const double* sibI = res.I(cat, cnResIdx);
          applyTransition(sibI, from_v.data(), t_sib);
          for (int i = 0; i < stride; ++i) out[i] *= from_v[i];
          continue;
        }
        const double* f = pg.Fslot(cat, par, sl);
        for (int i = 0; i < stride; ++i) out[i] *= f[i];
      }
    };

    int bSlot = topo.childSlot(a, b);
    siblingProduct(a, bSlot, curI.data());
    for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];

    int curNode = a;
    int parNode = topo.parentNode[a];
    while (parNode >= 1) {
      if (parNode == u) {
        applyTransition(curI.data(), contrib.data(), lMerge * rate * sc);
        siblingProduct(g, uSlot, curI.data());
        for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];
        curNode = g;
        parNode = topo.parentNode[g];
        continue;
      }
      int curSlot = topo.childSlot(parNode, curNode);
      double t = topo.edgeLen[topo.edgeToPar[curNode]] * rate * sc;
      applyTransition(curI.data(), contrib.data(), t);
      siblingProduct(parNode, curSlot, curI.data());
      for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];
      curNode = parNode;
      parNode = topo.parentNode[parNode];
    }

    // curI = CL at root for each pseudo-character (constant-state pattern)
    // P_const_cat = sum_s π_s × CL_root(s, cat)
    for (int c = 0; c < nChar; ++c) {
      double sl;
      if (pg.isMkN) {
        int off = c * 2;
        sl = rootF0 * curI[off] + rootF1 * curI[off + 1];
      } else {
        int off = c * kStates;
        sl = 0.0;
        for (int s = 0; s < kStates; ++s) sl += inv_k * curI[off + s];
      }
      totalConstProb += sl;
    }
  }

  return totalConstProb / nCat;
}

// (Old precompute_const_prod / adjusted_const_site_prob removed — replaced
// by create_const_pseudo_group + evaluate_const_prob above.)


#endif // MKPRIME_GIBBS_PARTIAL_CL_H
