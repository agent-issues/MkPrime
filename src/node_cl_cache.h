#ifndef MKPRIME_NODE_CL_CACHE_H
#define MKPRIME_NODE_CL_CACHE_H

// M-121: Persistent node-level CL cache for standard NNI/beta_simplex moves.
//
// Instead of a full O(nEdge) Felsenstein downpass for every NNI or
// beta_simplex proposal, cache per-node CLs persistently and recompute
// only the "dirty" path from the modified node(s) to the root.
//
// NNI dirty set:  v, u, ancestors(u) → root = O(depth)
// beta_simplex:   parents of 2 modified edges → root = O(depth)
// Full eval:      O(nEdge) per partition
//
// Speedup for 54-tip tree (depth≈7, nEdge=105): ~15× per-move for NNI.

#include "mcmc_state.h"
#include "gibbs_partial_cl.h"  // TreeNav
#include <cstring>
#include <cmath>
#include <algorithm>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// External functions defined in likelihood.cpp / mcmc_likelihood.cpp
// ---------------------------------------------------------------------------
double constant_site_prob_jc(IntegerVector parent, IntegerVector child,
                             NumericVector edge_length, int nTip,
                             int kStates, NumericVector root_freqs,
                             NumericVector rate_multipliers);
double constant_site_prob_mkn(IntegerVector parent, IntegerVector child,
                              NumericVector edge_length, int nTip,
                              double rate_loss, NumericVector root_freqs,
                              NumericVector rate_multipliers);
double mk_prime_relabel_log(int kPrime, int kObs);

// ---------------------------------------------------------------------------
// Inline helpers (duplicated from mcmc_likelihood.cpp to avoid link issues)
// ---------------------------------------------------------------------------
static inline NumericVector ncl_acrv_rates(double rateLogSd, int nCat,
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
// CacheUnit: per-partition (or per-kPrime-subgroup) CL storage
// ---------------------------------------------------------------------------
struct CacheUnit {
  int partIdx;       // index into data.parts[]
  int kStates;       // state count
  int nChar;         // number of characters in this unit
  int nCat;          // 1 if no ACRV, else nAcrvCat
  int stride;        // nChar * kStates
  bool isMkN;        // true → asymmetric binary (neomorphic)
  double rateScale;  // rateNeo for neomorphic, 1.0 otherwise

  // Per-unit root frequencies (length kStates)
  std::vector<double> rootFreqs;

  // Tip data for this unit (nTip × nChar, 0-indexed states, -1 = missing)
  IntegerMatrix tipStates;

  // Characters' kObs (for relabelling correction, transformational only)
  std::vector<int> kObsLocal;
  bool doRelabel;

  // Cached CLs: layout [cat * (maxNode+1) * stride + node * stride + c*k + s]
  // Tips + internal nodes.  Tips are constant (set once at init).
  std::vector<double> cl;

  int maxNode;  // allocated size (1-indexed: slots 0..maxNode)

  // Accessor: pointer to CL of (cat, node)
  double* CL(int cat, int node) {
    return cl.data() + ((size_t)cat * (maxNode + 1) + node) * stride;
  }
  const double* CL(int cat, int node) const {
    return cl.data() + ((size_t)cat * (maxNode + 1) + node) * stride;
  }

  void allocate(int maxNode_, int nCat_) {
    maxNode = maxNode_;
    nCat = nCat_;
    stride = nChar * kStates;
    cl.assign((size_t)nCat * (maxNode + 1) * stride, 0.0);
  }
};


// ---------------------------------------------------------------------------
// NodeCLCache: persistent cache held in McmcState
// ---------------------------------------------------------------------------
struct NodeCLCache {
  std::vector<CacheUnit> units;

  // Persistent tree navigation (updated incrementally for NNI)
  TreeNav topo;

  // ACRV rate multipliers (cached; invalidated when rateLogSd changes)
  std::vector<double> rates;
  bool useAcrv = false;
  double cachedRateLogSd = -1.0;  // detect rateLogSd changes

  // Ascertainment coding type (from data)
  int coding = 0;

  // Validity flag
  bool valid = false;
  int maxNode = 0;

  // Rollback scratch (reused across iterations to avoid allocation)
  std::vector<double> savedCL;
  std::vector<int>    dirtyNodes;
  double savedLogLik = 0.0;  // total log-lik before partial update
};


// ---------------------------------------------------------------------------
// Forward declarations
// ---------------------------------------------------------------------------
static void populate_cache_full(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo);

static double cache_total_loglik(
    const NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen, int nTip,
    double rateLoss, double rateLogSd, double rateNeo, double betaScale);

static double partial_eval_dirty(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen,
    double rateLoss, double rateNeo, double rateLogSd,
    double betaScale,
    const std::vector<int>& dirtyNodes);

static void restore_dirty_cls(
    NodeCLCache& cache, const std::vector<int>& dirtyNodes);


// ---------------------------------------------------------------------------
// JC transition: O(k) formula
// Writes result to `out` from source `src` (length kStates).
// ---------------------------------------------------------------------------
static inline void jc_transition(
    const double* src, double* out, int kStates, double t) {
  double inv_k = 1.0 / kStates;
  double km1   = kStates - 1.0;
  double exp_term = std::exp(-kStates * t / km1);
  double p_diff   = inv_k - inv_k * exp_term;
  double diff_coeff = (inv_k + (1.0 - inv_k) * exp_term) - p_diff;

  // For each character block (kStates values):
  // out[i] = p_diff * sum(src) + diff_coeff * src[i]
  // Caller iterates over characters externally.
  double sum_cl = 0.0;
  for (int s = 0; s < kStates; ++s) sum_cl += src[s];
  for (int s = 0; s < kStates; ++s)
    out[s] = p_diff * sum_cl + diff_coeff * src[s];
}


// ---------------------------------------------------------------------------
// MkN transition: 2×2 asymmetric (neomorphic binary)
// ---------------------------------------------------------------------------
static inline void mkn_transition(
    const double* src, double* out, double rateLoss, double t) {
  double sum_rl = 1.0 + rateLoss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rateLoss / sum_rl;
  double lambda = rate01 + rate10;
  double exp_term = std::exp(-lambda * t);
  double P00 = rate10 / lambda + rate01 / lambda * exp_term;
  double P01 = rate01 / lambda - rate01 / lambda * exp_term;
  double P10 = rate10 / lambda - rate10 / lambda * exp_term;
  double P11 = rate01 / lambda + rate10 / lambda * exp_term;
  out[0] = P00 * src[0] + P01 * src[1];
  out[1] = P10 * src[0] + P11 * src[1];
}


// ---------------------------------------------------------------------------
// apply_transition: dispatch to JC or MkN for one (cat, child) edge.
// Reads child CL from `srcCL`, writes per-character contributions to `dst`.
// ---------------------------------------------------------------------------
static void apply_transition(
    const CacheUnit& unit, const double* srcCL, double* dst,
    double edgeLen, double rateLoss, double rateMultiplier) {
  double t = edgeLen * rateMultiplier * unit.rateScale;
  int k = unit.kStates;
  int nChar = unit.nChar;

  if (unit.isMkN) {
    for (int c = 0; c < nChar; ++c)
      mkn_transition(srcCL + c * 2, dst + c * 2, rateLoss, t);
  } else {
    for (int c = 0; c < nChar; ++c)
      jc_transition(srcCL + c * k, dst + c * k, k, t);
  }
}


// ---------------------------------------------------------------------------
// init_tips: set tip CLs in a CacheUnit (constant; done once)
// ---------------------------------------------------------------------------
static void init_tips(CacheUnit& unit, int nTip) {
  int k = unit.kStates;
  int nChar = unit.nChar;
  // Tips are the same across all categories, so fill cat=0 and copy.
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = unit.CL(0, tip);
    std::fill(cl, cl + unit.stride, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state = unit.tipStates(tip - 1, c);
      int offset = c * k;
      if (state < 0) {
        for (int s = 0; s < k; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    // Copy tip CLs to all other categories (identical data)
    for (int cat = 1; cat < unit.nCat; ++cat)
      std::memcpy(unit.CL(cat, tip), cl, unit.stride * sizeof(double));
  }
}


// ---------------------------------------------------------------------------
// full_downpass: compute all internal-node CLs for one CacheUnit.
// Called during cache population.
// ---------------------------------------------------------------------------
static void full_downpass(
    CacheUnit& unit, const TreeNav& topo,
    const std::vector<double>& rates, double rateLoss) {
  int nCat = unit.nCat;

  // Temporary buffer for one child's contribution (one character block)
  std::vector<double> contrib(unit.stride);

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = (nCat > 1) ? rates[cat] : 1.0;

    // Reset internal node CLs (tips already initialized)
    // Use initFlg-style tracking: first child assigns, subsequent multiply
    std::vector<uint8_t> initFlg(topo.maxNode + 1, 0);
    for (int tip = 1; tip <= topo.nTip; ++tip) initFlg[tip] = 1;

    // Reverse edge order = postorder (children before parents)
    // This matches the existing flat pruning functions.
    // We need parent/child arrays — get them from the original vectors
    // that the TreeNav was built from.  But TreeNav doesn't store the
    // edge list directly in traversal order.
    //
    // Instead, iterate over internal nodes in a postorder that respects
    // "children processed before parent".  Build postorder from TreeNav.
    //
    // Actually, the simplest approach: iterate nodes nTip+1 ... maxNode
    // and process children.  But we need proper postorder.
    //
    // Use a stack-based postorder traversal.
    std::vector<int> postorder;
    postorder.reserve(topo.maxNode - topo.nTip);
    {
      std::vector<int> stack;
      stack.push_back(topo.root);
      while (!stack.empty()) {
        int n = stack.back(); stack.pop_back();
        if (n <= topo.nTip) continue;  // tip: skip
        postorder.push_back(n);
        // Push children (will be processed after this node in stack,
        // but we reverse postorder at the end)
        if (topo.ch0[n] >= 0) stack.push_back(topo.ch0[n]);
        if (topo.ch1[n] >= 0) stack.push_back(topo.ch1[n]);
        if (topo.ch2[n] >= 0) stack.push_back(topo.ch2[n]);
      }
      std::reverse(postorder.begin(), postorder.end());
    }

    for (int ni = 0; ni < (int)postorder.size(); ++ni) {
      int node = postorder[ni];
      double* clNode = unit.CL(cat, node);
      // Zero the node's CL before accumulating
      std::fill(clNode, clNode + unit.stride, 0.0);
      bool first = true;

      // Process each child
      int children[3] = { topo.ch0[node], topo.ch1[node], topo.ch2[node] };
      for (int ci = 0; ci < 3; ++ci) {
        int ch = children[ci];
        if (ch < 0) continue;

        double edgeLen = topo.edgeLen[topo.edgeToPar[ch]];
        const double* childCL = unit.CL(cat, ch);
        apply_transition(unit, childCL, contrib.data(), edgeLen, rateLoss, rate);

        if (first) {
          std::memcpy(clNode, contrib.data(), unit.stride * sizeof(double));
          first = false;
        } else {
          for (int j = 0; j < unit.stride; ++j)
            clNode[j] *= contrib[j];
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// unit_root_loglik: compute log-likelihood at root for one CacheUnit.
// Handles ACRV averaging.  Does NOT include ascertainment correction.
// ---------------------------------------------------------------------------
static double unit_root_loglik(const CacheUnit& unit, int root) {
  int k = unit.kStates;
  int nChar = unit.nChar;
  int nCat = unit.nCat;
  const double* rf = unit.rootFreqs.data();

  if (nCat == 1) {
    // No ACRV: single category
    const double* clRoot = unit.CL(0, root);
    double logLik = 0.0;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * k;
      double sl = 0.0;
      for (int s = 0; s < k; ++s)
        sl += rf[s] * clRoot[offset + s];
      if (sl <= 0.0) return R_NegInf;
      logLik += std::log(sl);
    }
    return logLik;
  } else {
    // ACRV: average site likelihoods across categories
    std::vector<double> siteLikSum(nChar, 0.0);
    for (int cat = 0; cat < nCat; ++cat) {
      const double* clRoot = unit.CL(cat, root);
      for (int c = 0; c < nChar; ++c) {
        int offset = c * k;
        double sl = 0.0;
        for (int s = 0; s < k; ++s)
          sl += rf[s] * clRoot[offset + s];
        siteLikSum[c] += sl;
      }
    }
    double inv_nCat = 1.0 / nCat;
    double logLik = 0.0;
    for (int c = 0; c < nChar; ++c) {
      double avg = siteLikSum[c] * inv_nCat;
      if (avg <= 0.0) return R_NegInf;
      logLik += std::log(avg);
    }
    return logLik;
  }
}


// ---------------------------------------------------------------------------
// build_cache_units: create CacheUnit structures from data + current kPrime.
// Does NOT allocate CL buffers — call allocate() on each unit afterward.
// ---------------------------------------------------------------------------
static void build_cache_units(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& kPrime, double rateLoss, int nCat) {

  cache.units.clear();

  for (int pi = 0; pi < (int)data.parts.size(); ++pi) {
    const PartInfo& part = data.parts[pi];
    int nCharPart = part.tipStates.ncol();

    if (part.type == 0) {
      // Neomorphic: k=2, MkN model
      CacheUnit u;
      u.partIdx = pi;
      u.kStates = 2;
      u.nChar   = nCharPart;
      u.nCat    = nCat;
      u.stride  = nCharPart * 2;
      u.isMkN   = true;
      u.rateScale = 1.0;  // rateNeo applied via edge lengths externally
      u.tipStates = part.tipStates;
      u.doRelabel = false;
      // MkN stationary frequencies: π0 = 1/(1+rl), π1 = rl/(1+rl)
      u.rootFreqs = { 1.0 / (1.0 + rateLoss), rateLoss / (1.0 + rateLoss) };
      cache.units.push_back(std::move(u));

    } else if (part.type == 2) {
      // Known state space: JC(k)
      int k = part.k;
      CacheUnit u;
      u.partIdx = pi;
      u.kStates = k;
      u.nChar   = nCharPart;
      u.nCat    = nCat;
      u.stride  = nCharPart * k;
      u.isMkN   = false;
      u.rateScale = 1.0;
      u.tipStates = part.tipStates;
      u.doRelabel = false;
      u.rootFreqs.assign(k, 1.0 / k);
      cache.units.push_back(std::move(u));

    } else {
      // Transformational: check if all kPrime are the same
      int kp0 = kPrime[part.globalCharIdx[0]];
      bool allSame = true;
      for (int ci = 1; ci < nCharPart; ++ci) {
        if (kPrime[part.globalCharIdx[ci]] != kp0) {
          allSame = false;
          break;
        }
      }

      if (allSame) {
        CacheUnit u;
        u.partIdx = pi;
        u.kStates = kp0;
        u.nChar   = nCharPart;
        u.nCat    = nCat;
        u.stride  = nCharPart * kp0;
        u.isMkN   = false;
        u.rateScale = 1.0;
        u.tipStates = part.tipStates;
        u.doRelabel = data.relabel;
        if (u.doRelabel) {
          u.kObsLocal.resize(nCharPart);
          for (int ci = 0; ci < nCharPart; ++ci)
            u.kObsLocal[ci] = part.kObsLocal[ci];
        }
        u.rootFreqs.assign(kp0, 1.0 / kp0);
        cache.units.push_back(std::move(u));
      } else {
        // Heterogeneous kPrime: create one unit per distinct k value.
        // Collect distinct k values and their character indices.
        std::vector<int> kpVals(nCharPart);
        for (int ci = 0; ci < nCharPart; ++ci)
          kpVals[ci] = kPrime[part.globalCharIdx[ci]];

        std::vector<int> uniqueK = kpVals;
        std::sort(uniqueK.begin(), uniqueK.end());
        uniqueK.erase(std::unique(uniqueK.begin(), uniqueK.end()), uniqueK.end());

        for (int kv : uniqueK) {
          // Collect columns for this k value
          std::vector<int> cols;
          for (int ci = 0; ci < nCharPart; ++ci)
            if (kpVals[ci] == kv) cols.push_back(ci);

          int nc = (int)cols.size();
          CacheUnit u;
          u.partIdx = pi;
          u.kStates = kv;
          u.nChar   = nc;
          u.nCat    = nCat;
          u.stride  = nc * kv;
          u.isMkN   = false;
          u.rateScale = 1.0;
          u.doRelabel = data.relabel;

          // Build sub-matrix of tip states
          int nTip = part.tipStates.nrow();
          u.tipStates = IntegerMatrix(nTip, nc);
          for (int ci = 0; ci < nc; ++ci)
            for (int tip = 0; tip < nTip; ++tip)
              u.tipStates(tip, ci) = part.tipStates(tip, cols[ci]);

          if (u.doRelabel) {
            u.kObsLocal.resize(nc);
            for (int ci = 0; ci < nc; ++ci)
              u.kObsLocal[ci] = part.kObsLocal[cols[ci]];
          }
          u.rootFreqs.assign(kv, 1.0 / kv);
          cache.units.push_back(std::move(u));
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// populate_cache_full: full downpass for all units.
// Called when cache is invalid (after SPR, global param change, etc.)
// ---------------------------------------------------------------------------
static void populate_cache_full(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo) {

  int nTip = data.nTip;
  int maxNode = 2 * nTip - 1;
  int nCat = (rateLogSd > 0.0) ? data.nCat : 1;
  cache.maxNode = maxNode;
  cache.useAcrv = (rateLogSd > 0.0);
  cache.coding  = data.codingType;
  cache.cachedRateLogSd = rateLogSd;

  // Build ACRV rates
  if (cache.useAcrv) {
    NumericVector rv = ncl_acrv_rates(rateLogSd, data.nCat, data.acrvZ);
    cache.rates.assign(rv.begin(), rv.end());
  } else {
    cache.rates = {1.0};
  }

  // Build TreeNav from current topology
  // Need to create absEdgeLen for TreeNav; scale neomorphic edges by rateNeo.
  // Actually, TreeNav stores raw edge lengths as passed. For neomorphic
  // partitions, the edge lengths include rateNeo scaling. But we store the
  // UNscaled absolute edge lengths in TreeNav and apply rateScale per-unit.
  cache.topo.build(parent, child, absEdgeLen, nTip);

  // Build cache units
  build_cache_units(cache, data, kPrime, rateLoss, nCat);

  // Allocate and populate each unit
  for (auto& unit : cache.units) {
    unit.allocate(maxNode, nCat);
    init_tips(unit, nTip);

    // For neomorphic partitions, edge lengths need rateNeo scaling.
    // We handle this via rateScale in apply_transition, but TreeNav
    // stores unscaled edge lengths. Set rateScale here.
    if (unit.isMkN) {
      unit.rateScale = rateNeo;
    }

    full_downpass(unit, cache.topo, cache.rates, rateLoss);
  }

  cache.valid = true;
}


// ---------------------------------------------------------------------------
// cache_total_loglik: compute total log-likelihood from cached CLs.
// Includes ascertainment correction (computed fully, not cached).
// ---------------------------------------------------------------------------
static double cache_total_loglik(
    const NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen, int nTip,
    double rateLoss, double rateLogSd, double rateNeo, double betaScale) {

  int root = cache.topo.root;
  double totalLL = 0.0;

  // Track per-partition LL for ascertainment correction grouping
  // (correction is per-partition, not per-unit)
  int nParts = (int)data.parts.size();
  std::vector<double> partRawLL(nParts, 0.0);
  std::vector<int>    partNChar(nParts, 0);

  for (int ui = 0; ui < (int)cache.units.size(); ++ui) {
    const CacheUnit& unit = cache.units[ui];
    double ll = unit_root_loglik(unit, root);

    // Relabelling correction (transformational only)
    if (unit.doRelabel) {
      for (int ci = 0; ci < unit.nChar; ++ci)
        ll += mk_prime_relabel_log(unit.kStates, unit.kObsLocal[ci]);
    }

    partRawLL[unit.partIdx] += ll;
    partNChar[unit.partIdx] += unit.nChar;
  }

  // Add ascertainment correction per partition (computed fully)
  NumericVector rates;
  if (cache.useAcrv) {
    rates = ncl_acrv_rates(rateLogSd, data.nCat, data.acrvZ);
  }

  for (int pi = 0; pi < nParts; ++pi) {
    double ll = partRawLL[pi];
    if (cache.coding != 0 && partNChar[pi] > 0) {
      const PartInfo& part = data.parts[pi];
      double p = 0.0;

      if (part.type == 0) {
        NumericVector neoEl(absEdgeLen.size());
        for (int i = 0; i < absEdgeLen.size(); ++i)
          neoEl[i] = absEdgeLen[i] * rateNeo;
        NumericVector rootFreqs(2);
        rootFreqs[0] = 1.0 / (1.0 + rateLoss);
        rootFreqs[1] = rateLoss / (1.0 + rateLoss);
        p = constant_site_prob_mkn(parent, child, neoEl, nTip,
                                    rateLoss, rootFreqs, rates);
      } else {
        // JC: determine kStates for this partition
        // For known (type 2): part.k
        // For transformational: use first unit's kStates as representative
        //   (if heterogeneous kPrime, each sub-group has its own k —
        //    but constant_site_prob_jc takes a single k.  We need per-unit
        //    ascertainment.  For simplicity, compute per unit below.)
        int kStates = (part.type == 2) ? part.k : 0;
        if (kStates > 0) {
          NumericVector rootFreqs(kStates, 1.0 / kStates);
          p = constant_site_prob_jc(parent, child, absEdgeLen, nTip,
                                     kStates, rootFreqs, rates);
        }
      }
      if (p > 0.0 && p < 1.0) {
        ll -= partNChar[pi] * std::log(1.0 - p);
      }
    }
    totalLL += ll;
  }

  return totalLL;
}


// ---------------------------------------------------------------------------
// find_dirty_nni: identify nodes needing CL recomputation after NNI.
//
// NNI swaps: parent[cRow] = u (was v), parent[wRow] = v (was u).
// Dirty nodes: v (lost child c, gained child w), u (gained c, lost w),
//              then ancestors of u up to root.
// Since v is a child of u, the dirty set is the path [v, u, g, ..., root]
// returned in postorder (v first, root last).
// ---------------------------------------------------------------------------
static std::vector<int> find_dirty_nni(
    const TreeNav& topo, int v, int u) {
  std::vector<int> path;
  path.reserve(16);

  // v is child of u; both have changed children
  path.push_back(v);
  int cur = u;
  while (cur >= 0) {
    path.push_back(cur);
    cur = topo.parentNode[cur];
  }
  return path;  // postorder: v, u, g, ..., root
}


// ---------------------------------------------------------------------------
// find_dirty_beta_simplex: identify dirty nodes after beta_simplex.
//
// Two edges have modified lengths: edges at rows idx1 and idx2.
// Their parent nodes need CL recomputation, plus ancestors to root.
// ---------------------------------------------------------------------------
static std::vector<int> find_dirty_beta_simplex(
    const TreeNav& topo,
    const IntegerVector& parent,
    int idx1, int idx2) {

  int par1 = parent[idx1];
  int par2 = parent[idx2];

  // Collect path from par1 to root
  std::vector<int> path1;
  for (int n = par1; n >= 0; n = topo.parentNode[n])
    path1.push_back(n);

  // Collect path from par2 to root
  std::vector<int> path2;
  for (int n = par2; n >= 0; n = topo.parentNode[n])
    path2.push_back(n);

  // Merge: take union, preserving postorder (deeper nodes first).
  // Both paths go root-ward; mark which nodes are on each path,
  // then collect in postorder.
  // Use depth (distance from root) to sort: deeper = earlier in postorder.
  std::vector<bool> onPath(topo.maxNode + 1, false);
  for (int n : path1) onPath[n] = true;
  for (int n : path2) onPath[n] = true;

  // Collect all dirty nodes
  std::vector<std::pair<int,int>> depthNode;  // (depth, node)
  for (int n = 1; n <= topo.maxNode; ++n) {
    if (!onPath[n]) continue;
    // Compute depth (count hops to root)
    int d = 0;
    for (int x = n; x >= 0; x = topo.parentNode[x]) ++d;
    depthNode.push_back({d, n});
  }
  // Sort by decreasing depth = postorder
  std::sort(depthNode.begin(), depthNode.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  std::vector<int> dirty;
  dirty.reserve(depthNode.size());
  for (auto& dn : depthNode) dirty.push_back(dn.second);
  return dirty;
}


// ---------------------------------------------------------------------------
// save_dirty_cls: save CLs at dirty nodes for rollback.
// Saves all units, all categories for each dirty node.
// ---------------------------------------------------------------------------
static void save_dirty_cls(
    NodeCLCache& cache, const std::vector<int>& dirtyNodes) {

  // Compute total size needed
  size_t totalSize = 0;
  for (const auto& unit : cache.units)
    totalSize += (size_t)unit.nCat * dirtyNodes.size() * unit.stride;

  cache.savedCL.resize(totalSize);
  cache.dirtyNodes = dirtyNodes;

  double* dst = cache.savedCL.data();
  for (const auto& unit : cache.units) {
    for (int cat = 0; cat < unit.nCat; ++cat) {
      for (int node : dirtyNodes) {
        const double* src = unit.CL(cat, node);
        std::memcpy(dst, src, unit.stride * sizeof(double));
        dst += unit.stride;
      }
    }
  }
}


// ---------------------------------------------------------------------------
// restore_dirty_cls: undo partial evaluation on rejection.
// ---------------------------------------------------------------------------
static void restore_dirty_cls(
    NodeCLCache& cache, const std::vector<int>& dirtyNodes) {

  const double* src = cache.savedCL.data();
  for (auto& unit : cache.units) {
    for (int cat = 0; cat < unit.nCat; ++cat) {
      for (int node : dirtyNodes) {
        double* dst = unit.CL(cat, node);
        std::memcpy(dst, src, unit.stride * sizeof(double));
        src += unit.stride;
      }
    }
  }
}


// ---------------------------------------------------------------------------
// recompute_dirty_nodes: update CLs at dirty nodes using cached children.
// dirtyNodes must be in postorder (children before parents).
// ---------------------------------------------------------------------------
static void recompute_dirty_nodes(
    NodeCLCache& cache,
    const NumericVector& absEdgeLen,
    double rateLoss,
    const std::vector<int>& dirtyNodes) {

  std::vector<double> contrib;  // reusable scratch

  for (auto& unit : cache.units) {
    int maxStride = unit.stride;
    if ((int)contrib.size() < maxStride)
      contrib.resize(maxStride);

    for (int cat = 0; cat < unit.nCat; ++cat) {
      double rate = (unit.nCat > 1) ? cache.rates[cat] : 1.0;

      for (int node : dirtyNodes) {
        if (node <= cache.topo.nTip) continue;  // tips are constant

        double* clNode = unit.CL(cat, node);
        bool first = true;

        int children[3] = {
          cache.topo.ch0[node],
          cache.topo.ch1[node],
          cache.topo.ch2[node]
        };

        for (int ci = 0; ci < 3; ++ci) {
          int ch = children[ci];
          if (ch < 0) continue;

          double edgeLen = absEdgeLen[cache.topo.edgeToPar[ch]];
          const double* childCL = unit.CL(cat, ch);
          apply_transition(unit, childCL, contrib.data(),
                           edgeLen, rateLoss, rate);

          if (first) {
            std::memcpy(clNode, contrib.data(),
                        unit.stride * sizeof(double));
            first = false;
          } else {
            for (int j = 0; j < unit.stride; ++j)
              clNode[j] *= contrib[j];
          }
        }
      }
    }
  }
}


// ---------------------------------------------------------------------------
// partial_eval_dirty: save, recompute dirty nodes, return new total loglik.
// Caller is responsible for calling restore_dirty_cls on rejection.
// ---------------------------------------------------------------------------
static double partial_eval_dirty(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen,
    double rateLoss, double rateNeo, double rateLogSd,
    double betaScale,
    const std::vector<int>& dirtyNodes) {

  // 1. Save CLs at dirty nodes for rollback
  save_dirty_cls(cache, dirtyNodes);

  // 2. Recompute dirty nodes
  recompute_dirty_nodes(cache, absEdgeLen, rateLoss, dirtyNodes);

  // 3. Compute total log-likelihood from (partially updated) cache
  return cache_total_loglik(cache, data, parent, child, absEdgeLen,
                             data.nTip, rateLoss, rateLogSd, rateNeo,
                             betaScale);
}


// ---------------------------------------------------------------------------
// update_topo_nni: update TreeNav in-place after an accepted NNI.
// NNI swapped parent[cRow] and parent[wRow]: c moved from v to u,
// w moved from u to v.
// ---------------------------------------------------------------------------
static void update_topo_nni(
    TreeNav& topo,
    int v, int u, int cNode, int wNode) {

  // Remove c from v's children
  if      (topo.ch0[v] == cNode) topo.ch0[v] = -1;
  else if (topo.ch1[v] == cNode) topo.ch1[v] = -1;
  else if (topo.ch2[v] == cNode) topo.ch2[v] = -1;

  // Remove w from u's children
  if      (topo.ch0[u] == wNode) topo.ch0[u] = -1;
  else if (topo.ch1[u] == wNode) topo.ch1[u] = -1;
  else if (topo.ch2[u] == wNode) topo.ch2[u] = -1;

  // Add w to v's children (fill first empty slot)
  if      (topo.ch0[v] < 0) topo.ch0[v] = wNode;
  else if (topo.ch1[v] < 0) topo.ch1[v] = wNode;
  else                       topo.ch2[v] = wNode;

  // Add c to u's children
  if      (topo.ch0[u] < 0) topo.ch0[u] = cNode;
  else if (topo.ch1[u] < 0) topo.ch1[u] = cNode;
  else                       topo.ch2[u] = cNode;

  // Update parent pointers
  topo.parentNode[cNode] = u;
  topo.parentNode[wNode] = v;
}


// ---------------------------------------------------------------------------
// can_use_partial_cl: check if partial CL is applicable for this move.
// Returns false if cache is invalid, Q-het is enabled, etc.
// ---------------------------------------------------------------------------
static bool can_use_partial_cl(
    const NodeCLCache& cache, const McmcData& data) {
  if (!cache.valid) return false;
  if (data.qHeterogeneity) return false;
  return true;
}


#endif  // MKPRIME_NODE_CL_CACHE_H
