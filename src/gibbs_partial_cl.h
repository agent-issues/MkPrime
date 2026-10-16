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
#include <algorithm>
#include <cstring>
#include <array>
#include <cmath>
#include <limits>
#include "mkn_rates.h"
#include "f81.h"
#include "../inst/include/MkPrime/acrv.h"
#include "../inst/include/MkPrime/scaling.h"

using namespace Rcpp;

// ---------------------------------------------------------------------------
// CL rescaling (scaling.h). Every CL vector here -- cached I, residual I, and
// the working vectors of a candidate evaluation -- has per-character shifts:
// its block for character c is the true CL times 2^shift[c]. A vector with
// no shifts recorded has all of them 0. A child's F-slot contribution
// carries the shifts of that child's I.
// ---------------------------------------------------------------------------

// Shifts of one working CL vector.
struct CLShifts {
  int nChar;
  std::vector<int> v;
  bool on = false;

  explicit CLShifts(int nChar_) : nChar(nChar_) {}
  void reset() { on = false; }
  int* ensure() {
    if (!on) {
      v.assign(nChar, 0);
      on = true;
    }
    return v.data();
  }
  void add(const int* s) {
    if (!s) return;
    int* p = ensure();
    for (int c = 0; c < nChar; ++c) p[c] += s[c];
  }
  const int* get() const { return on ? v.data() : nullptr; }
};

// Shifts of an array of CL vectors, one slot per vector.
struct ShiftStore {
  int nChar = 0;
  std::vector<int> v;  // allocated at the first ensure()
  std::vector<uint8_t> on;

  void allocate(size_t nSlot, int nChar_) {
    nChar = nChar_;
    v.clear();
    on.assign(nSlot, 0);
  }
  void clear_all() { std::fill(on.begin(), on.end(), 0); }
  bool any() const { return std::find(on.begin(), on.end(), 1) != on.end(); }
  const int* get(size_t i) const {
    return on[i] ? v.data() + i * nChar : nullptr;
  }
  void clear(size_t i) { on[i] = 0; }
  int* ensure(size_t i) {
    if (v.empty()) v.assign(on.size() * nChar, 0);
    int* p = v.data() + i * nChar;
    if (!on[i]) {
      std::fill(p, p + nChar, 0);
      on[i] = 1;
    }
    return p;
  }
  void add(size_t i, const int* s) {
    if (!s) return;
    int* p = ensure(i);
    for (int c = 0; c < nChar; ++c) p[c] += s[c];
  }
};

// Rescales each character block of cl (nChar blocks of k) that has fallen
// below kRescaleBelow; ensure() gives the shifts to add to.
template <class Ensure>
inline void rescale_small(double* cl, int nChar, int k, Ensure ensure) {
  for (int c = 0; c < nChar; ++c) {
    double* block = cl + c * k;
    double sum = 0.0;
    for (int s = 0; s < k; ++s) sum += block[s];
    if (sum < MkPrime::kRescaleBelow) {
      const int sh = MkPrime::rescale_cl(block, k);
      if (sh) ensure()[c] += sh;
    }
  }
}

// cur *= f, then, if `check`, rescale cur as rescale_small().
inline void mul_rescale(double* cur, const double* f, int nChar, int k,
                        CLShifts& sh, bool check) {
  if (!check) {
    for (int i = 0; i < nChar * k; ++i) cur[i] *= f[i];
    return;
  }
  for (int c = 0; c < nChar; ++c) {
    double* block = cur + c * k;
    const double* fb = f + c * k;
    double sum = 0.0;
    for (int s = 0; s < k; ++s) {
      block[s] *= fb[s];
      sum += block[s];
    }
    if (sum < MkPrime::kRescaleBelow) {
      const int s = MkPrime::rescale_cl(block, k);
      if (s) sh.ensure()[c] += s;
    }
  }
}

// Candidate paths are first multiplied out without rescaling, and again with
// it if any character's root site likelihood is below kUncheckedFloor
// (scaling.h). Cached CLs are rescaled or sum to at least kRescaleBelow, so
// the checked pass cannot underflow either.
inline bool roots_too_small(const double* rootSl, int nChar) {
  for (int c = 0; c < nChar; ++c) {
    if (!(rootSl[c] >= MkPrime::kUncheckedFloor)) return true;
  }
  return false;
}

inline int shift_at(const int* s, int c) { return s ? s[c] : 0; }

// x * 2^-shift, the true value of a scaled site likelihood.
inline double unscale(double x, int shift) {
  return shift ? std::ldexp(x, -shift) : x;
}


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
  // M-158: node-indexed edge lengths (decoupled from edge array ordering).
  // nodeEdgeLen[c] = absolute edge length from node c to parentNode[c].
  // Enables partial CL recomputation without requiring edge array indices.
  std::vector<double> nodeEdgeLen;

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
    nodeEdgeLen.assign(maxNode + 1, 0.0);

    for (int e = 0; e < nEdge; ++e) edgeLen[e] = absEdgeLen[e];

    // Fill parent/child maps.  Postorder (reverse) determines which child
    // is ch0 (first encountered), ch1 (second), ch2 (third = root only).
    for (int e = nEdge - 1; e >= 0; --e) {
      int p = parent[e], c = child[e];
      parentNode[c] = p;
      edgeToPar[c]  = e;
      nodeEdgeLen[c] = absEdgeLen[e];  // M-158
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
  // Whether the group carries its own ascertainment term; a transformational
  // group's is its partition's joint term (joint_asc.h).
  bool ownAsc = false;
  double rateLoss;  // only for isMkN
  double rateScale; // partition-rate scale (audit Issue 1):
                    //   neoScale = r/(1+r) * (n_neo+n_trans)/n_neo  for neo
                    //   transScale = 1/(1+r) * (n_neo+n_trans)/n_trans  for others
                    // Reduces to 1.0 when n_neo == 0 or n_trans == 0.

  // M-114: F81 transitions for Q-heterogeneity
  bool useF81 = false;
  std::vector<double> f81Pi;  // frequency vector, sized to kStates
  double f81Mu = 0.0;         // 1 / (1 - Σπ²)

  IntegerMatrix tipData;  // nTip × nChar, 0-indexed states, -1 missing

  // Characters by missing-data mask; for a pseudo-character group, the masks
  // its pseudo-characters stand for (k per mask, in this order).
  MaskTally masks;

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
  // Shifts of I(cat, node), in slot ISlot(cat, node).
  ShiftStore Ish;
  size_t ISlot(int cat, int node) const {
    return (size_t)cat * (maxNode + 1) + node;
  }
  const int* IS(int cat, int node) const { return Ish.get(ISlot(cat, node)); }

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
    Ish.allocate((size_t)nCat * (maxNode + 1), nChar);
    f81Pi.assign(kStates, 0.0);
  }

  // Apply transition matrix: dst = P(t) × src for all characters.
  // Dispatches to F81, MkN, or JC based on model flags.
  void transition(const double* src, double* dst, double t) const;

  // Root-frequency-weighted site likelihood: Σ_s π_s × rootCL[c*k+s]
  double rootSiteLik(const double* rootCL, int c) const;
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
  NumericVector rates(nCat);
  MkPrime::lognormal_rates(rateLogSd, nCat, acrvZ.data(), rates.begin());
  return rates;
}


// ---------------------------------------------------------------------------
// JC / MkN transition helpers (kernels in inst/include/MkPrime/transition.h)
// ---------------------------------------------------------------------------

// p_diff and diff_coeff (= exp_term) for JC(k) with branch length t
inline void jc_trans_params(int k, double t,
                            double& p_diff, double& diff_coeff) {
  const MkPrime::JCProbs p = MkPrime::jc_probs(k, t);
  p_diff     = p.pDiff;
  diff_coeff = p.expTerm;
}

using MkPrime::jc_transition;

// MkN (asymmetric binary k=2) transition helpers
inline void mkn_trans_params(double rateLoss, double t,
                              double& P00, double& P01,
                              double& P10, double& P11) {
  const MkPrime::MkNProbs p = MkPrime::MkNRates(rateLoss).probs(t);
  P00 = p.p00;
  P01 = p.p01;
  P10 = p.p10;
  P11 = p.p11;
}

using MkPrime::mkn_transition;

// --- CLGroup::transition / rootSiteLik implementations ---

inline void CLGroup::transition(const double* src, double* dst, double t) const {
  if (useF81) {
    f81_transition(src, dst, nChar, kStates, f81Pi.data(), f81Mu, t);
  } else if (isMkN) {
    double P00, P01, P10, P11;
    mkn_trans_params(rateLoss, t, P00, P01, P10, P11);
    mkn_transition(src, dst, nChar, P00, P01, P10, P11);
  } else {
    double p_diff, diff_coeff;
    jc_trans_params(kStates, t, p_diff, diff_coeff);
    jc_transition(src, dst, nChar, kStates, p_diff, diff_coeff);
  }
}

inline double CLGroup::rootSiteLik(const double* rootCL, int c) const {
  if (useF81) {
    int off = c * kStates;
    double sl = 0.0;
    for (int s = 0; s < kStates; ++s) sl += f81Pi[s] * rootCL[off + s];
    return sl;
  } else if (isMkN) {
    double rf0 = rateLoss / (1.0 + rateLoss);  // π₀ = rl/(1+rl)
    double rf1 = 1.0 / (1.0 + rateLoss);       // π₁ = 1/(1+rl)
    int off = c * 2;
    return rf0 * rootCL[off] + rf1 * rootCL[off + 1];
  } else {
    double inv_k = 1.0 / kStates;
    int off = c * kStates;
    double sl = 0.0;
    for (int s = 0; s < kStates; ++s) sl += inv_k * rootCL[off + s];
    return sl;
  }
}


// ---------------------------------------------------------------------------
// Caching downpass: full postorder traversal, storing I, F0, F1 per node
// ---------------------------------------------------------------------------
// One pass of caching_downpass(). Unchecked, it rescales nothing and returns
// false if any character's root site likelihood is below kRescaleBelow.
// Otherwise no node can have fallen below it either (the largest CL entry
// never increases towards the root), so the result is exactly what a checked
// pass gives.
static bool caching_downpass_pass(CLGroup& grp, const TreeNav& topo,
                                  const IntegerVector& parent,
                                  const IntegerVector& child,
                                  const NumericVector& rates, bool check) {
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

      // ch_ is complete once its own edge is reached.
      if (check && ch_ > nTip) {
        const size_t slot = grp.ISlot(cat, ch_);
        rescale_small(grp.I(cat, ch_), nChar, kStates,
                      [&]() { return grp.Ish.ensure(slot); });
      }

      // Compute transition: contrib = P(t) × I[ch_]
      grp.transition(grp.I(cat, ch_), contrib.data(), t);

      if (check) {
        if (nChildSeen[par] == 0) grp.Ish.clear(grp.ISlot(cat, par));
        grp.Ish.add(grp.ISlot(cat, par), grp.IS(cat, ch_));
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
    if (!check) {
      const double* rootCl = grp.I(cat, topo.root);
      for (int c = 0; c < nChar; ++c) {
        if (!(grp.rootSiteLik(rootCl, c) >= MkPrime::kRescaleBelow)) {
          return false;
        }
      }
    }
  }
  if (!check) grp.Ish.clear_all();
  return true;
}

static void caching_downpass(CLGroup& grp, const TreeNav& topo,
                              const IntegerVector& parent,
                              const IntegerVector& child,
                              const NumericVector& rates) {
  if (!caching_downpass_pass(grp, topo, parent, child, rates, false)) {
    caching_downpass_pass(grp, topo, parent, child, rates, true);
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

  // Shifts of I(cat, pathIdx), in slot cat * nPathNodes + pathIdx. F at
  // path[0] carries the shifts of I(cat, sibNode) in the cached group; F at
  // path[pi] those of I(cat, pi - 1).
  ShiftStore Ish;
  int sibNode;
  const int* IS(int cat, int pathIdx) const {
    return Ish.get((size_t)cat * nPathNodes + pathIdx);
  }
  const int* FS(const CLGroup& grp, int cat, int pathIdx) const {
    return pathIdx ? IS(cat, pathIdx - 1) : grp.IS(cat, sibNode);
  }

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


// One pass of compute_residual_cl(); unchecked as caching_downpass_pass(),
// which is exact only when the cached group carries no shifts.
static bool compute_residual_cl_pass(
    ResidualCL& res,
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int u, int sibNode, double lMerge, bool check)
{
  int g = topo.parentNode[u];   // grandparent of pruned subtree
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
  res.Ish.allocate((size_t)nCat * res.nPathNodes, grp.nChar);
  res.sibNode = sibNode;

  std::vector<double> contrib(stride);
  const int nChar = grp.nChar, k = grp.kStates;

  // Shifts of slot `pi` of the residual path: the cached children's at
  // `node`, except `exclude`, plus `extra`; then rescale I(cat, pi).
  auto finishPathNode = [&](int cat, int pi, int node, int exclude,
                            const int* extra) {
    const size_t slot = (size_t)cat * res.nPathNodes + pi;
    res.Ish.clear(slot);
    const int kids[3] = { topo.ch0[node], topo.ch1[node], topo.ch2[node] };
    for (int sl = 0; sl < 3; ++sl) {
      if (sl != exclude && kids[sl] >= 0) {
        res.Ish.add(slot, grp.IS(cat, kids[sl]));
      }
    }
    res.Ish.add(slot, extra);
    rescale_small(res.I(cat, pi), nChar, k,
                  [&]() { return res.Ish.ensure(slot); });
  };

  // Record the replaced slot at g: u's slot is replaced by sibNode via lMerge
  int uSlot = topo.childSlot(g, u);
  res.replacedSlot[0] = uSlot;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // At g: replace u's contribution with sibNode's contribution via lMerge
    grp.transition(grp.I(cat, sibNode), contrib.data(),
                   lMerge * rate * grp.rateScale);

    // Store the replacement F at g
    std::memcpy(res.F(cat, 0), contrib.data(), stride * sizeof(double));

    // I_R[g] = new_sibNode_contrib × product of all other children's contributions
    double* Ig = res.I(cat, 0);
    int nChildG = (topo.ch2[g] >= 0) ? 3 : 2;
    grp.productExcluding(cat, g, uSlot, Ig, nChildG);
    for (int i = 0; i < stride; ++i) Ig[i] *= contrib[i];
    if (check) finishPathNode(cat, 0, g, uSlot, grp.IS(cat, sibNode));

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
      grp.transition(prevI, contrib.data(), t);

      // Store the replacement F at this path node
      std::memcpy(res.F(cat, pi), contrib.data(), stride * sizeof(double));

      // Combine with product of all unchanged sibling contributions
      double* Inode = res.I(cat, pi);
      int nChild = (topo.ch2[node] >= 0) ? 3 : 2;
      grp.productExcluding(cat, node, prevSlot, Inode, nChild);
      for (int i = 0; i < stride; ++i) Inode[i] *= contrib[i];
      if (check) finishPathNode(cat, pi, node, prevSlot, res.IS(cat, pi - 1));
    }
    if (!check) {
      const double* rootCl = res.I(cat, res.nPathNodes - 1);
      for (int c = 0; c < nChar; ++c) {
        if (!(grp.rootSiteLik(rootCl, c) >= MkPrime::kRescaleBelow)) {
          return false;
        }
      }
    }
  }
  return true;
}

static void compute_residual_cl(
    ResidualCL& res,
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int u, int sibNode, double lMerge)
{
  if (grp.Ish.any() ||
      !compute_residual_cl_pass(res, grp, topo, rates, u, sibNode, lMerge,
                                false)) {
    compute_residual_cl_pass(res, grp, topo, rates, u, sibNode, lMerge, true);
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
//
// mergedEdge (GSPR-004): evaluate the regraft onto the MERGED edge of the
// residual tree — the (g, sibNode) pair created by detaching v — which is the
// subtree's own position.  The caller passes a = g, b = sibNode, and
// lHalfReg = lMerge / 2.  That edge does not exist as a row of the original
// tree, so the only navigational difference is the slot excluded at a: u's
// slot (where the reattached u sits in the new tree) instead of
// topo.childSlot(a, b), which is undefined for the pair.  Everything else —
// from_b via grp.I(cat, sibNode) (off the residual path), the sibling
// product at g, and the upward propagation — is the generic code path.
// ---------------------------------------------------------------------------
// If siteLikAccum is non-null (M-114 Q-het streaming mode): accumulate
// per-site raw likelihoods into it (no log, no /nCat) and return 0.0.
// If null (default): return the log-likelihood as before.
static double evaluate_candidate(
    const CLGroup& grp, const TreeNav& topo,
    const ResidualCL& res,
    const NumericVector& rates,
    int v, int u,           // pruned subtree root and its parent
    int sibNode,            // u's other child (sibling of v)
    double lMerge,          // collapsed edge length sibNode→g
    int a, int b,           // regraft parent and child nodes
    double lHalfReg,        // lReg / 2 (half the regraft edge length)
    double lPrune,          // prune edge length (u→v)
    double* siteLikAccum = nullptr,
    bool mergedEdge = false)
{
  int nChar   = grp.nChar;
  int stride  = grp.stride;
  int nCat    = grp.nCat;
  double sc   = grp.rateScale;
  int g       = topo.parentNode[u];
  int uSlot   = topo.childSlot(g, u);

  // Per-site likelihood accumulator: external (Q-het streaming) or local.
  // The external one holds true likelihoods, so a character below ~1e-308
  // underflows there; the local one carries shifts.
  std::vector<double> localSiteLik;
  std::vector<int> localSiteShift;
  double* accum;
  if (siteLikAccum) {
    accum = siteLikAccum;  // caller manages zeroing
  } else {
    localSiteLik.assign(nChar, 0.0);
    localSiteShift.assign(nChar, 0);
    accum = localSiteLik.data();
  }

  // Buffers for per-character computations
  std::vector<double> from_v(stride), from_b(stride), Iu(stride), contrib(stride);
  std::vector<double> curI(stride);
  CLShifts shU(nChar), shCur(nChar), shBelow(nChar);
  std::vector<double> rootSl(nChar);
  const int k = grp.kStates;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Unchecked first; again with rescaling if the root is small.
    for (int pass = 0; pass < 2; ++pass) {
      const bool check = pass > 0;
      // Contribution from v (pruned subtree) at u
      grp.transition(grp.I(cat, v), from_v.data(), lPrune * rate * sc);

      // Contribution from b at u (via half regraft edge)
      int bResIdx = res.pathIndex(b);
      const double* Ib = (bResIdx >= 0) ? res.I(cat, bResIdx) : grp.I(cat, b);
      grp.transition(Ib, from_b.data(), lHalfReg * rate * sc);

      // I at u = from_v × from_b
      for (int i = 0; i < stride; ++i) Iu[i] = from_v[i] * from_b[i];
      shU.reset();
      shU.add(grp.IS(cat, v));
      shU.add((bResIdx >= 0) ? res.IS(cat, bResIdx) : grp.IS(cat, b));
      if (check) {
        rescale_small(Iu.data(), nChar, k, [&]() { return shU.ensure(); });
      }

      // from_u at a: transition through half regraft edge from u to a
      grp.transition(Iu.data(), contrib.data(), lHalfReg * rate * sc);

      // Helper lambda: compute the product of all sibling contributions at node
      // 'par' EXCEPT the child in slot 'excludeSlot', with their shifts.
      //
      // Three cases for each sibling slot:
      //   1. par is on the residual path and slot == replacedSlot → use res.F
      //   2. Sibling node is on the residual path → recompute P(t) × res.I
      //   3. Otherwise → use original cached F
      auto siblingProduct = [&](int par, int excludeSlot, double* out,
                                CLShifts& sh) {
        int nChild = (topo.ch2[par] >= 0) ? 3 : 2;
        std::fill(out, out + stride, 1.0);
        sh.reset();
        int parResIdx = res.pathIndex(par);
        for (int sl = 0; sl < nChild; ++sl) {
          if (sl == excludeSlot) continue;

          // Case 1: this slot was replaced during residual computation
          if (parResIdx >= 0 && sl == res.replacedSlot[parResIdx]) {
            const double* f = res.F(cat, parResIdx);
            for (int i = 0; i < stride; ++i) out[i] *= f[i];
            sh.add(res.FS(grp, cat, parResIdx));
            continue;
          }

          // Case 2: sibling node is on the residual path
          int cn = (sl == 0) ? topo.ch0[par] :
                   (sl == 1) ? topo.ch1[par] : topo.ch2[par];
          int cnResIdx = res.pathIndex(cn);
          if (cnResIdx >= 0) {
            double t_sib = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            const double* sibI = res.I(cat, cnResIdx);
            grp.transition(sibI, from_v.data(), t_sib);
            for (int i = 0; i < stride; ++i) out[i] *= from_v[i];
            sh.add(res.IS(cat, cnResIdx));
            continue;
          }

          // Case 3: use original cached F
          const double* f = grp.Fslot(cat, par, sl);
          for (int i = 0; i < stride; ++i) out[i] *= f[i];
          sh.add(grp.IS(cat, cn));
        }
      };

      // At a: I_cand[a] = from_u_at_a × product_of_siblings(a, exclude b).
      // For the merged edge (a = g, b = sibNode), b occupies u's original slot.
      int bSlot = mergedEdge ? uSlot : topo.childSlot(a, b);
      siblingProduct(a, bSlot, curI.data(), shCur);
      shCur.add(shU.get());
      mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);

      // Propagate up from a to root.
      // Key: if we encounter u as a parent, skip it — u doesn't exist on
      // this path in the new topology.  Instead connect to g via lMerge.
      int curNode = a;
      int parNode = topo.parentNode[a];
      while (parNode >= 1) {
        if (parNode == u) {
          grp.transition(curI.data(), contrib.data(), lMerge * rate * sc);
          std::swap(shBelow, shCur);
          siblingProduct(g, uSlot, curI.data(), shCur);
          shCur.add(shBelow.get());
          mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);
          curNode = g;
          parNode = topo.parentNode[g];
          continue;
        }

        int curSlot = topo.childSlot(parNode, curNode);
        double t = topo.edgeLen[topo.edgeToPar[curNode]] * rate * sc;
        grp.transition(curI.data(), contrib.data(), t);
        std::swap(shBelow, shCur);
        siblingProduct(parNode, curSlot, curI.data(), shCur);
        shCur.add(shBelow.get());
        mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);
        curNode = parNode;
        parNode = topo.parentNode[parNode];
      }

      for (int c = 0; c < nChar; ++c) {
        rootSl[c] = grp.rootSiteLik(curI.data(), c);
      }
      if (check || !roots_too_small(rootSl.data(), nChar)) break;
    }

    // curI now holds I_cand[root].  Accumulate per-site likelihood.
    const int* shRoot = shCur.get();
    for (int c = 0; c < nChar; ++c) {
      const double sl = rootSl[c];
      if (siteLikAccum) {
        accum[c] += unscale(sl, shift_at(shRoot, c));
      } else {
        MkPrime::add_scaled(accum[c], localSiteShift[c], sl,
                            shift_at(shRoot, c));
      }
    }
  }

  if (siteLikAccum) return 0.0;  // caller takes log after all components

  // Standalone mode: log-likelihood = sum_c log(siteLikSum[c] / nCat)
  double logLik = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = accum[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += MkPrime::log_scaled(avg, localSiteShift[c]);
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
static CLGroup create_const_pseudo_group(
    const CLGroup& src, int nTip, int maxNode, int nCat,
    const std::vector<std::vector<uint8_t>>& missingMasks) {
  CLGroup pg;
  pg.isMkN      = src.isMkN;
  pg.rateLoss   = src.rateLoss;
  pg.rateScale  = src.rateScale;
  pg.masks      = src.masks;
  int k = src.kStates;
  const bool complete = src.masks.ids.empty() || src.masks.complete();
  const int nMask = complete ? 1 : (int)src.masks.ids.size();

  // Tip data: k pseudo-characters per mask, each with all observed tips in
  // state s and the mask's missing tips marginalised
  IntegerMatrix tips(nTip, k * nMask);
  for (int mi = 0; mi < nMask; ++mi) {
    const uint8_t* missing =
      complete ? nullptr : missingMasks[src.masks.ids[mi]].data();
    for (int c = 0; c < k; ++c)
      for (int t = 0; t < nTip; ++t)
        tips(t, mi * k + c) = (missing && missing[t]) ? -1 : c;
  }
  pg.tipData = tips;
  pg.allocate(maxNode, nCat, k * nMask, k);
  return pg;
}

// Sum over the pseudo-group's masks of n_m * log(1 - P_m), where P_m is the
// sum of the mask's k per-pseudo-character accumulators over denom.
// Terms with P_m >= 1 are dropped.
static double pseudo_asc_log1m(const CLGroup& pg, const double* accum,
                               double denom) {
  const int k = pg.kStates;
  double sum = 0.0;
  for (size_t mi = 0; mi < pg.masks.ids.size(); ++mi) {
    double constP = 0.0;
    for (int c = 0; c < k; ++c) constP += accum[mi * k + c];
    constP /= denom;
    if (constP < 1.0) sum += pg.masks.counts[mi] * std::log(1.0 - constP);
  }
  return sum;
}

// Evaluate P_const for a candidate regraft using partial CLs on pseudo-chars.
// Same algorithm as evaluate_candidate but returns P(constant site) instead
// of log-likelihood.  mergedEdge as in evaluate_candidate.
// If constProbAccum is non-null (M-114 Q-het streaming), accumulate the raw
// per-pseudo-character constant-site likelihoods and return 0.0.
static double evaluate_const_prob(
    const CLGroup& pg,          // pseudo-character group
    const TreeNav& topo,
    const ResidualCL& res,
    const NumericVector& rates,
    int v, int u, int sibNode, double lMerge,
    int a, int b, double lHalfReg, double lPrune,
    double* constProbAccum = nullptr,
    bool mergedEdge = false)
{
  int nChar   = pg.nChar;  // = kStates
  int stride  = pg.stride;
  int nCat    = pg.nCat;
  double sc   = pg.rateScale;
  int g       = topo.parentNode[u];
  int uSlot   = topo.childSlot(g, u);

  std::vector<double> from_v(stride), from_b(stride), Iu(stride), contrib(stride);
  std::vector<double> curI(stride);
  CLShifts shU(nChar), shCur(nChar), shBelow(nChar);
  std::vector<double> rootSl(nChar);
  const int k = pg.kStates;

  double totalConstProb = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Unchecked first; again with rescaling if the root is small.
    for (int pass = 0; pass < 2; ++pass) {
      const bool check = pass > 0;
      pg.transition(pg.I(cat, v), from_v.data(), lPrune * rate * sc);

      int bResIdx = res.pathIndex(b);
      const double* Ib = (bResIdx >= 0) ? res.I(cat, bResIdx) : pg.I(cat, b);
      pg.transition(Ib, from_b.data(), lHalfReg * rate * sc);

      for (int i = 0; i < stride; ++i) Iu[i] = from_v[i] * from_b[i];
      shU.reset();
      shU.add(pg.IS(cat, v));
      shU.add((bResIdx >= 0) ? res.IS(cat, bResIdx) : pg.IS(cat, b));
      if (check) {
        rescale_small(Iu.data(), nChar, k, [&]() { return shU.ensure(); });
      }

      pg.transition(Iu.data(), contrib.data(), lHalfReg * rate * sc);

      // siblingProduct lambda (same as evaluate_candidate)
      auto siblingProduct = [&](int par, int excludeSlot, double* out,
                                CLShifts& sh) {
        int nChild = (topo.ch2[par] >= 0) ? 3 : 2;
        std::fill(out, out + stride, 1.0);
        sh.reset();
        int parResIdx = res.pathIndex(par);
        for (int sl = 0; sl < nChild; ++sl) {
          if (sl == excludeSlot) continue;
          if (parResIdx >= 0 && sl == res.replacedSlot[parResIdx]) {
            const double* f = res.F(cat, parResIdx);
            for (int i = 0; i < stride; ++i) out[i] *= f[i];
            sh.add(res.FS(pg, cat, parResIdx));
            continue;
          }
          int cn = (sl == 0) ? topo.ch0[par] :
                   (sl == 1) ? topo.ch1[par] : topo.ch2[par];
          int cnResIdx = res.pathIndex(cn);
          if (cnResIdx >= 0) {
            double t_sib = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            const double* sibI = res.I(cat, cnResIdx);
            pg.transition(sibI, from_v.data(), t_sib);
            for (int i = 0; i < stride; ++i) out[i] *= from_v[i];
            sh.add(res.IS(cat, cnResIdx));
            continue;
          }
          const double* f = pg.Fslot(cat, par, sl);
          for (int i = 0; i < stride; ++i) out[i] *= f[i];
          sh.add(pg.IS(cat, cn));
        }
      };

      int bSlot = mergedEdge ? uSlot : topo.childSlot(a, b);
      siblingProduct(a, bSlot, curI.data(), shCur);
      shCur.add(shU.get());
      mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);

      int curNode = a;
      int parNode = topo.parentNode[a];
      while (parNode >= 1) {
        if (parNode == u) {
          pg.transition(curI.data(), contrib.data(), lMerge * rate * sc);
          std::swap(shBelow, shCur);
          siblingProduct(g, uSlot, curI.data(), shCur);
          shCur.add(shBelow.get());
          mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);
          curNode = g;
          parNode = topo.parentNode[g];
          continue;
        }
        int curSlot = topo.childSlot(parNode, curNode);
        double t = topo.edgeLen[topo.edgeToPar[curNode]] * rate * sc;
        pg.transition(curI.data(), contrib.data(), t);
        std::swap(shBelow, shCur);
        siblingProduct(parNode, curSlot, curI.data(), shCur);
        shCur.add(shBelow.get());
        mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);
        curNode = parNode;
        parNode = topo.parentNode[parNode];
      }

      for (int c = 0; c < nChar; ++c) {
        rootSl[c] = pg.rootSiteLik(curI.data(), c);
      }
      if (check || !roots_too_small(rootSl.data(), nChar)) break;
    }

    const int* shRoot = shCur.get();
    for (int c = 0; c < nChar; ++c) {
      double sl = unscale(rootSl[c], shift_at(shRoot, c));
      if (constProbAccum)
        constProbAccum[c] += sl;
      else
        totalConstProb += sl;
    }
  }

  if (constProbAccum) return 0.0;
  return totalConstProb / nCat;
}

// (Old precompute_const_prod / adjusted_const_site_prob removed — replaced
// by create_const_pseudo_group + evaluate_const_prob above.)


// ---------------------------------------------------------------------------
// M-111: Partial CL evaluation for subtree swap
//
// A subtree swap (nodeA ↔ nodeB) changes CLs at two parent nodes
// (pA, pB) and propagates changes up two paths that merge at their LCA.
// Unlike SPR, no node is removed — both subtrees stay in the tree.
// Branch lengths stay at their parent positions (swap only moves children).
//
// Algorithm:
//   1. Cache all CLs via caching_downpass (shared across candidates)
//   2. For each candidate nodeB:
//      a. Find LCA of pA and pB via precomputed pathA flags
//      b. Process three segments bottom-up:
//         - segA (pA → LCA): swap nodeA→nodeB at pA, propagate
//         - segB (pB → LCA): swap nodeB→nodeA at pB, propagate
//         - shared (LCA → root): merge both changes, propagate
//      c. Accumulate root CL into site likelihoods
//
// constProbMode: false → return log-likelihood
//                true  → return P(constant site) for ascertainment
// ---------------------------------------------------------------------------
// M-114: added siteLikAccum / constProbAccum for Q-het streaming.
// When siteLikAccum is non-null and !constProbMode, per-site raw likelihoods
// are accumulated into siteLikAccum and return 0.0.
// When constProbAccum is non-null and constProbMode, per-pseudo-char raw
// const-prob likelihoods are accumulated and return 0.0.
static double evaluate_swap_impl(
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int nodeA, int nodeB,
    int pA, int slotA, double lenA,
    const std::vector<int>& pathA,
    const std::vector<int>& pathAIdx,
    bool constProbMode,
    double* siteLikAccum = nullptr,
    double* constProbAccum = nullptr)
{
  int nChar   = grp.nChar;
  int stride  = grp.stride;
  int nCat    = grp.nCat;
  double sc   = grp.rateScale;

  int pB    = topo.parentNode[nodeB];
  int slotB = topo.childSlot(pB, nodeB);
  double lenB = topo.edgeLen[topo.edgeToPar[nodeB]];

  // --- Find LCA: walk from pB upward until hitting a node on pathA ---
  std::vector<int> pathBbelow;
  pathBbelow.reserve(16);
  int lcaIdxA = -1;
  for (int n = pB; n >= 1; n = topo.parentNode[n]) {
    int ai = pathAIdx[n];
    if (ai >= 0) {
      lcaIdxA = ai;
      break;
    }
    pathBbelow.push_back(n);
  }

  int nSegA = lcaIdxA;
  int nSegB = (int)pathBbelow.size();

  // --- Buffers ---
  std::vector<double> curI(stride);
  std::vector<double> contrib(stride);
  std::vector<double> lastA(stride);
  std::vector<double> lastB(stride);
  std::vector<double> prevShared(stride);
  CLShifts shCur(nChar), shLastA(nChar), shLastB(nChar), shPrev(nChar);
  std::vector<double> rootSl(nChar);
  const int k = grp.kStates;

  // Accumulators; only the local one carries shifts (see evaluate_candidate)
  std::vector<double> localSiteLik;
  std::vector<int> localSiteShift;
  double totalConstProb = 0.0;
  double* accum = nullptr;
  if (constProbMode) {
    // constProbAccum: accumulate raw const-prob per pseudo-char if streaming
    // (otherwise accumulate into totalConstProb as before)
  } else if (siteLikAccum) {
    accum = siteLikAccum;
  } else {
    localSiteLik.assign(nChar, 0.0);
    localSiteShift.assign(nChar, 0);
    accum = localSiteLik.data();
  }

  // curI *= contrib; the last factor of a node also rescales when checking.
  bool check = false;
  auto combine = [&](bool last) {
    if (last) {
      mul_rescale(curI.data(), contrib.data(), nChar, k, shCur, check);
    } else {
      for (int i = 0; i < stride; ++i) curI[i] *= contrib[i];
    }
  };

  auto childAt = [&](int n, int sl) -> int {
    return (sl == 0) ? topo.ch0[n] : (sl == 1) ? topo.ch1[n] : topo.ch2[n];
  };

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Unchecked first; again with rescaling if the root is small.
    for (int pass = 0; pass < 2; ++pass) {
      check = pass > 0;
      // ============ Segment A: pathA[0 .. lcaIdxA-1] ============
      for (int ai = 0; ai < nSegA; ++ai) {
        int n = pathA[ai];
        int nChild = (topo.ch2[n] >= 0) ? 3 : 2;
        std::fill(curI.begin(), curI.end(), 1.0);
        shCur.reset();

        for (int sl = 0; sl < nChild; ++sl) {
          int cn = childAt(n, sl);

          if (n == pA && sl == slotA) {
            // SWAP: nodeB replaces nodeA at pA
            grp.transition(grp.I(cat, nodeB), contrib.data(),
                            lenA * rate * sc);
            shCur.add(grp.IS(cat, nodeB));
          } else if (ai > 0 && cn == pathA[ai - 1]) {
            // Propagated from previous segA node
            double t = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            grp.transition(lastA.data(), contrib.data(), t);
            shCur.add(shLastA.get());
          } else {
            std::memcpy(contrib.data(), grp.Fslot(cat, n, sl),
                        stride * sizeof(double));
            shCur.add(grp.IS(cat, cn));
          }
          combine(sl + 1 == nChild);
        }
        std::memcpy(lastA.data(), curI.data(), stride * sizeof(double));
        std::swap(shLastA, shCur);
      }

      // ============ Segment B: pathBbelow[0 .. nSegB-1] ============
      for (int bi = 0; bi < nSegB; ++bi) {
        int n = pathBbelow[bi];
        int nChild = (topo.ch2[n] >= 0) ? 3 : 2;
        std::fill(curI.begin(), curI.end(), 1.0);
        shCur.reset();

        for (int sl = 0; sl < nChild; ++sl) {
          int cn = childAt(n, sl);

          if (n == pB && sl == slotB) {
            // SWAP: nodeA replaces nodeB at pB
            grp.transition(grp.I(cat, nodeA), contrib.data(),
                            lenB * rate * sc);
            shCur.add(grp.IS(cat, nodeA));
          } else if (bi > 0 && cn == pathBbelow[bi - 1]) {
            // Propagated from previous segB node
            double t = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            grp.transition(lastB.data(), contrib.data(), t);
            shCur.add(shLastB.get());
          } else {
            std::memcpy(contrib.data(), grp.Fslot(cat, n, sl),
                        stride * sizeof(double));
            shCur.add(grp.IS(cat, cn));
          }
          combine(sl + 1 == nChild);
        }
        std::memcpy(lastB.data(), curI.data(), stride * sizeof(double));
        std::swap(shLastB, shCur);
      }

      // ============ Shared segment: pathA[lcaIdxA .. end] ============
      int nPathA = (int)pathA.size();
      for (int ai = lcaIdxA; ai < nPathA; ++ai) {
        int n = pathA[ai];
        int nChild = (topo.ch2[n] >= 0) ? 3 : 2;
        std::fill(curI.begin(), curI.end(), 1.0);
        shCur.reset();

        for (int sl = 0; sl < nChild; ++sl) {
          int cn = childAt(n, sl);
          bool handled = false;

          // Swap conditions (fire at LCA when pA or pB IS the LCA)
          if (n == pA && sl == slotA) {
            grp.transition(grp.I(cat, nodeB), contrib.data(),
                            lenA * rate * sc);
            shCur.add(grp.IS(cat, nodeB));
            handled = true;
          }
          if (!handled && n == pB && sl == slotB) {
            grp.transition(grp.I(cat, nodeA), contrib.data(),
                            lenB * rate * sc);
            shCur.add(grp.IS(cat, nodeA));
            handled = true;
          }

          // Dirty children feeding into LCA from segments below
          if (!handled && nSegA > 0 && ai == lcaIdxA &&
              cn == pathA[lcaIdxA - 1]) {
            double t = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            grp.transition(lastA.data(), contrib.data(), t);
            shCur.add(shLastA.get());
            handled = true;
          }
          if (!handled && nSegB > 0 && ai == lcaIdxA &&
              cn == pathBbelow.back()) {
            double t = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            grp.transition(lastB.data(), contrib.data(), t);
            shCur.add(shLastB.get());
            handled = true;
          }

          // Previous shared-segment node
          if (!handled && ai > lcaIdxA && cn == pathA[ai - 1]) {
            double t = topo.edgeLen[topo.edgeToPar[cn]] * rate * sc;
            grp.transition(prevShared.data(), contrib.data(), t);
            shCur.add(shPrev.get());
            handled = true;
          }

          if (!handled) {
            std::memcpy(contrib.data(), grp.Fslot(cat, n, sl),
                        stride * sizeof(double));
            shCur.add(grp.IS(cat, cn));
          }

          combine(sl + 1 == nChild);
        }

        std::memcpy(prevShared.data(), curI.data(), stride * sizeof(double));
        std::swap(shPrev, shCur);
      }

      for (int c = 0; c < nChar; ++c) {
        rootSl[c] = grp.rootSiteLik(curI.data(), c);
      }
      if (check || !roots_too_small(rootSl.data(), nChar)) break;
    }

    // ============ Root CL → accumulate site likelihoods ============
    // The root's shifts are in shPrev after the last swap.
    const int* shRoot = shPrev.get();
    for (int c = 0; c < nChar; ++c) {
      double sl = rootSl[c];
      if (constProbMode) {
        sl = unscale(sl, shift_at(shRoot, c));
        if (constProbAccum)
          constProbAccum[c] += sl;
        else
          totalConstProb += sl;
      } else if (siteLikAccum) {
        accum[c] += unscale(sl, shift_at(shRoot, c));
      } else {
        MkPrime::add_scaled(accum[c], localSiteShift[c], sl,
                            shift_at(shRoot, c));
      }
    }
  }  // end category loop

  // ============ Return ============
  if (constProbMode) {
    if (constProbAccum) return 0.0;
    return totalConstProb / nCat;
  }
  if (siteLikAccum) return 0.0;  // caller takes log after all components

  double logLik = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = accum[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += MkPrime::log_scaled(avg, localSiteShift[c]);
  }
  return logLik;
}

// Convenience wrappers
static inline double evaluate_swap_candidate(
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int nodeA, int nodeB,
    int pA, int slotA, double lenA,
    const std::vector<int>& pathA,
    const std::vector<int>& pathAIdx)
{
  return evaluate_swap_impl(grp, topo, rates, nodeA, nodeB,
                            pA, slotA, lenA, pathA, pathAIdx, false);
}

static inline double evaluate_swap_const_prob(
    const CLGroup& grp, const TreeNav& topo,
    const NumericVector& rates,
    int nodeA, int nodeB,
    int pA, int slotA, double lenA,
    const std::vector<int>& pathA,
    const std::vector<int>& pathAIdx)
{
  return evaluate_swap_impl(grp, topo, rates, nodeA, nodeB,
                            pA, slotA, lenA, pathA, pathAIdx, true);
}


// ---------------------------------------------------------------------------
// M-114: Beta discretization for Q-heterogeneity
//
// Duplicated from mcmc_likelihood.cpp::compute_het_bins (which is file-local).
// Discretizes Beta(α, (k-1)α) into nBins equal-probability bins, each
// represented by its conditional mean.
// ---------------------------------------------------------------------------
static void gibbs_compute_het_bins(double alpha, int k, int nBins,
                                    double* bins) {
  double a = alpha;
  double b = (k - 1.0) * alpha;
  for (int i = 0; i < nBins; ++i) {
    double lo = R::qbeta((double)i / nBins, a, b, 1, 0);
    double hi = R::qbeta((double)(i + 1) / nBins, a, b, 1, 0);
    if (hi - lo < 1e-15) {
      bins[i] = 0.5 * (lo + hi);
      continue;
    }
    double p_lo = R::pbeta(lo, a + 1.0, b, 1, 0);
    double p_hi = R::pbeta(hi, a + 1.0, b, 1, 0);
    double denom = R::pbeta(hi, a, b, 1, 0) - R::pbeta(lo, a, b, 1, 0);
    if (denom < 1e-300) {
      bins[i] = 0.5 * (lo + hi);
    } else {
      // A mean can round above 1, whose rate would overflow.
      bins[i] = std::min(1.0, (a / (a + b)) * (p_hi - p_lo) / denom);
    }
  }
}


// ---------------------------------------------------------------------------
// M-114: F81 component parameter construction for Q-heterogeneity
//
// Mirrors the π construction in pruning_f81_het_acrv_flat (mcmc_likelihood.cpp).
// Sets grp.f81Pi[], grp.f81Mu, grp.useF81 = true.
// ---------------------------------------------------------------------------
static void set_f81_component(CLGroup& grp, double betaVal, int rot,
                               double baseRL) {
  int k = grp.kStates;
  double sumPiSq = 0.0;

  if (k == 2) {
    // Neomorphic: asymmetric binary using rateLoss
    double sum_rl = 1.0 + baseRL;
    double gain_base = 1.0 / sum_rl;
    double loss_base = baseRL / sum_rl;
    double gain_b = gain_base * 2.0 * betaVal;
    double loss_b = loss_base * 2.0 * (1.0 - betaVal);
    double total_rate = gain_b + loss_b;
    grp.f81Pi[1] = gain_b / total_rate;
    grp.f81Pi[0] = 1.0 - grp.f81Pi[1];
    sumPiSq = grp.f81Pi[0] * grp.f81Pi[0] + grp.f81Pi[1] * grp.f81Pi[1];
  } else {
    // Symmetric k≥3: one elevated frequency at position 'rot'
    double r = (1.0 - betaVal) / (k - 1.0);
    for (int s = 0; s < k; ++s) grp.f81Pi[s] = r;
    grp.f81Pi[rot] = betaVal;
    sumPiSq = betaVal * betaVal + (k - 1.0) * r * r;
  }

  grp.f81Mu  = 1.0 / (1.0 - sumPiSq);
  grp.useF81 = true;
}

// Convert per-site raw likelihood accumulators to log-likelihood.
// siteLikSum[c] was accumulated across totalComp components; divide by
// totalComp and take log. The accumulators hold true likelihoods, without
// shifts, so a character below the normal range (~1e-308) cannot be scored:
// its candidate would silently weigh nothing, so stop the run instead.
static double siteLikAccum_to_logLik(const double* siteLikSum, int nChar,
                                      int totalComp) {
  double logLik = 0.0;
  double inv_comp = 1.0 / totalComp;
  for (int c = 0; c < nChar; ++c) {
    double avg = siteLikSum[c] * inv_comp;
    if (avg < std::numeric_limits<double>::min()) {
      Rcpp::stop(
        "A character's likelihood is below 1e-308, which the "
        "Q-heterogeneity Gibbs topology moves cannot represent (#346). Set "
        "`gibbs_spr` and `gibbs_subtree_swap` to 0 in `moveWeights`, or drop "
        "near-random characters with many states.");
    }
    logLik += std::log(avg);
  }
  return logLik;
}


#endif // MKPRIME_GIBBS_PARTIAL_CL_H
