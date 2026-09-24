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
#include "fast_exp.h"
#include <cstring>
#include <map>
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
double singleton_site_prob_jc(IntegerVector parent, IntegerVector child,
                              NumericVector edge_length, int nTip,
                              int kStates, NumericVector root_freqs,
                              NumericVector rate_multipliers);
double singleton_site_prob_mkn(IntegerVector parent, IntegerVector child,
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
  bool clValid = false;  // M-161: per-unit CL validity

  // Per-unit root frequencies (length kStates)
  std::vector<double> rootFreqs;

  // Tip data for this unit (nTip × nChar, 0-indexed states, -1 = missing)
  IntegerMatrix tipStates;

  // Characters by missing-data mask, for the ascertainment correction
  MaskTally masks;

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

  // M-161: two-level validity
  // Level 1: global topology + unit structure
  bool topoValid = false;       // TreeNav matches current topology + edge lengths
  bool structureValid = false;  // unit structure matches current kPrime grouping
  // Level 2: per-unit CLs (CacheUnit::clValid)

  int maxNode = 0;

  // Rollback scratch (reused across iterations to avoid allocation)
  std::vector<double> savedCL;
  std::vector<int>    dirtyNodes;
  double savedLogLik = 0.0;  // total log-lik before partial update

  // M-161: diagnostic counter for selective repopulation
  int diagSelectivePopCount = 0;

  // M-161: helpers
  bool ready() const {
    if (!topoValid || !structureValid) return false;
    for (const auto& u : units)
      if (!u.clValid) return false;
    return true;
  }

  void invalidate_all() {
    topoValid = false;
  }

  void invalidate_structure() {
    structureValid = false;
  }

  void invalidate_neo_cls() {
    for (auto& u : units)
      if (u.isMkN) u.clValid = false;
  }

  void invalidate_all_cls() {
    for (auto& u : units)
      u.clValid = false;
  }
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
// apply_transition: dispatch to JC or MkN for one (cat, child) edge.
// Reads child CL from `srcCL`, writes per-character contributions to `dst`.
// The transition probabilities depend only on (k, t), so they are computed
// once per edge rather than once per character.
// ---------------------------------------------------------------------------
static void apply_transition(
    const CacheUnit& unit, const double* srcCL, double* dst,
    double edgeLen, double rateLoss, double rateMultiplier) {
  double t = edgeLen * rateMultiplier * unit.rateScale;
  if (unit.isMkN) {
    double P00, P01, P10, P11;
    mkn_trans_params(rateLoss, t, P00, P01, P10, P11);
    mkn_transition(srcCL, dst, unit.nChar, P00, P01, P10, P11);
  } else {
    double pDiff, diffCoeff;
    jc_trans_params(unit.kStates, t, pDiff, diffCoeff);
    jc_transition(srcCL, dst, unit.nChar, unit.kStates, pDiff, diffCoeff);
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

        double el = topo.nodeEdgeLen[ch];  // M-158: node-indexed
        const double* childCL = unit.CL(cat, ch);
        apply_transition(unit, childCL, contrib.data(), el, rateLoss, rate);

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
      // MkN stationary frequencies: π0 = rl/(1+rl), π1 = 1/(1+rl)
      u.rootFreqs = { rateLoss / (1.0 + rateLoss), 1.0 / (1.0 + rateLoss) };
      u.masks = tally_masks(data, part.globalCharIdx, nCharPart);
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
      u.masks = tally_masks(data, part.globalCharIdx, nCharPart);
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
        u.masks = tally_masks(data, part.globalCharIdx, nCharPart);
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
          std::vector<int> unitGlobal(nc);
          for (int ci = 0; ci < nc; ++ci)
            unitGlobal[ci] = part.globalCharIdx[cols[ci]];
          u.masks = tally_masks(data, unitGlobal, nc);
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

  // Audit Issue 1: RB-style partition-rate normalisation. Both neo (MkN)
  // and trans/known units are scaled symmetrically so the nChar-weighted
  // mean partition rate equals 1.
  const PartitionScales pScales =
      compute_partition_scales(rateNeo, data.nNeo, data.nTrans);

  // Allocate and populate each unit
  for (auto& unit : cache.units) {
    unit.allocate(maxNode, nCat);
    init_tips(unit, nTip);

    // Per-unit rateScale (applied to edge lengths in apply_transition).
    unit.rateScale = unit.isMkN ? pScales.neo : pScales.trans;

    full_downpass(unit, cache.topo, cache.rates, rateLoss);
    unit.clValid = true;
  }

  cache.topoValid = true;
  cache.structureValid = true;
}


// ---------------------------------------------------------------------------
// populate_cache: smart wrapper that uses selective repopulation when
// topology and unit structure are still valid.  Falls back to full rebuild
// when either global flag is stale.  (M-161)
// ---------------------------------------------------------------------------
static void populate_cache(
    NodeCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absEdgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo) {

  // Fast path: full rebuild when topology or unit structure is stale
  if (!cache.topoValid || !cache.structureValid) {
    populate_cache_full(cache, data, parent, child, absEdgeLen,
                        kPrime, rateLoss, rateLogSd, rateNeo);
    return;
  }

  // --- Selective repopulation: topo + structure valid, some units stale ---
  cache.diagSelectivePopCount++;

  // Update ACRV rates if rateLogSd changed
  if (cache.cachedRateLogSd != rateLogSd) {
    cache.useAcrv = (rateLogSd > 0.0);
    if (cache.useAcrv) {
      NumericVector rv = ncl_acrv_rates(rateLogSd, data.nCat, data.acrvZ);
      cache.rates.assign(rv.begin(), rv.end());
    } else {
      cache.rates = {1.0};
    }
    cache.cachedRateLogSd = rateLogSd;
  }

  // Audit Issue 1: re-derive partition scales (rateNeo may have changed).
  const PartitionScales pScales =
      compute_partition_scales(rateNeo, data.nNeo, data.nTrans);

  // Repopulate only invalid units
  for (auto& unit : cache.units) {
    if (unit.clValid) continue;

    // Update unit-specific parameters that may have changed
    if (unit.isMkN) {
      unit.rootFreqs = { rateLoss / (1.0 + rateLoss),
                         1.0 / (1.0 + rateLoss) };
    }
    unit.rateScale = unit.isMkN ? pScales.neo : pScales.trans;

    // Rerun full downpass — tip CLs are still valid (character data never
    // changes), so only internal-node CLs are recomputed.
    full_downpass(unit, cache.topo, cache.rates, rateLoss);
    unit.clValid = true;
  }
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
  // Bug fix: always provide at least one rate category (1.0) so that
  // constant_site_prob_* functions iterate their pruning loop.
  // Previously, when ACRV was off, `rates` was empty (size 0), causing
  // nCat=0 inside constant_site_prob → loop skipped → correction = 0.
  NumericVector rates;
  if (cache.useAcrv) {
    rates = ncl_acrv_rates(rateLogSd, data.nCat, data.acrvZ);
  } else {
    rates = NumericVector(1, 1.0);
  }

  // Audit Issue 1: partition-rate normalisation also applies to
  // ascertainment correction — the constant-site probability must use the
  // same effective edge lengths as the pruning, so pre-scale by neoScale
  // (type 0) or transScale (type 1 / type 2).
  const PartitionScales pScales =
      compute_partition_scales(rateNeo, data.nNeo, data.nTrans);
  NumericVector neoAscEl(absEdgeLen.size());
  NumericVector transAscEl(absEdgeLen.size());
  for (int i = 0; i < absEdgeLen.size(); ++i) {
    neoAscEl[i]   = absEdgeLen[i] * pScales.neo;
    transAscEl[i] = absEdgeLen[i] * pScales.trans;
  }

  // The JC correction depends only on (k, missing-data mask), so units that
  // share one (across or within partitions) share one tree traversal.
  std::vector<double> jcAscByK;
  auto jcAscProb = [&](int k) {
    if (k >= (int)jcAscByK.size()) jcAscByK.resize(k + 1, -1.0);
    if (jcAscByK[k] < 0.0) {
      NumericVector rootFreqs(k, 1.0 / k);
      double pk = constant_site_prob_jc(parent, child, transAscEl, nTip,
                                        k, rootFreqs, rates);
      // LIKE-001 fix: informative coding adds the JC singleton probability.
      if (cache.coding == 2) {
        pk += singleton_site_prob_jc(parent, child, transAscEl, nTip,
                                     k, rootFreqs, rates);
      }
      jcAscByK[k] = pk;
    }
    return jcAscByK[k];
  };
  std::map<std::pair<int, int>, double> jcAscByKMask;
  auto jcAscProbsMasked = [&](int k, const std::vector<int>& ids) {
    std::vector<int> todo;
    for (int m : ids) {
      if (m != 0 && !jcAscByKMask.count({k, m})) todo.push_back(m);
    }
    if (!todo.empty()) {
      const std::vector<double> pNew = asc_probs_masked(
        data, parent, child, transAscEl, 1, k, 1.0, betaScale, rates, todo);
      for (size_t i = 0; i < todo.size(); ++i) {
        jcAscByKMask[{k, todo[i]}] = pNew[i];
      }
    }
    std::vector<double> p(ids.size(), 0.0);
    for (size_t i = 0; i < ids.size(); ++i) {
      if (ids[i] != 0) p[i] = jcAscByKMask[{k, ids[i]}];
    }
    return p;
  };

  for (int pi = 0; pi < nParts; ++pi) {
    double ll = partRawLL[pi];
    if (cache.coding != 0 && partNChar[pi] > 0) {
      for (int ui = 0; ui < (int)cache.units.size(); ++ui) {
        const CacheUnit& unit = cache.units[ui];
        if (unit.partIdx != pi) continue;
        int k = unit.kStates;
        if (k <= 0) continue;
        if (unit.isMkN) {
          NumericVector rootFreqs(2);
          rootFreqs[0] = rateLoss / (1.0 + rateLoss);
          rootFreqs[1] = 1.0 / (1.0 + rateLoss);
          double p = constant_site_prob_mkn(parent, child, neoAscEl, nTip,
                                            rateLoss, rootFreqs, rates);
          // LIKE-001 fix: informative coding adds the MkN singleton term.
          if (cache.coding == 2) {
            p += singleton_site_prob_mkn(parent, child, neoAscEl, nTip,
                                          rateLoss, rootFreqs, rates);
          }
          ll -= masked_asc_log1m(unit.masks, p, asc_probs_masked(
            data, parent, child, neoAscEl, 0, 2, rateLoss, betaScale, rates,
            unit.masks.ids), true);
        } else {
          ll -= masked_asc_log1m(unit.masks, jcAscProb(k),
                                 jcAscProbsMasked(k, unit.masks.ids), true);
        }
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
// find_dirty_dirichlet: identify dirty nodes after a K-element Dirichlet
// simplex proposal that modified edge lengths at arbitrary edge indices.
// Generalizes find_dirty_beta_simplex from 2 paths to K paths.
// ---------------------------------------------------------------------------
static std::vector<int> find_dirty_dirichlet(
    const TreeNav& topo,
    const IntegerVector& parent,
    const std::vector<int>& edgeIndices) {

  // Mark all nodes that lie on any path from a modified edge's parent to root
  std::vector<bool> onPath(topo.maxNode + 1, false);
  for (int idx : edgeIndices) {
    for (int n = parent[idx]; n >= 0; n = topo.parentNode[n])
      onPath[n] = true;
  }

  // Collect marked nodes with their depth for postorder sorting
  std::vector<std::pair<int,int>> depthNode;
  for (int n = 1; n <= topo.maxNode; ++n) {
    if (!onPath[n]) continue;
    int d = 0;
    for (int x = n; x >= 0; x = topo.parentNode[x]) ++d;
    depthNode.push_back({d, n});
  }

  // Sort by decreasing depth = postorder (children before parents)
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
// M-158: uses nodeEdgeLen (decoupled from edge array ordering).
// ---------------------------------------------------------------------------
static void recompute_dirty_nodes(
    NodeCLCache& cache,
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

          double el = cache.topo.nodeEdgeLen[ch];  // M-158
          const double* childCL = unit.CL(cat, ch);
          apply_transition(unit, childCL, contrib.data(),
                           el, rateLoss, rate);

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

  // 2. Recompute dirty nodes (uses nodeEdgeLen, M-158)
  recompute_dirty_nodes(cache, rateLoss, dirtyNodes);

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
  if (!cache.ready()) return false;
  if (data.qHeterogeneity) return false;
  return true;
}


// ===========================================================================
// M-158: Partial CL evaluation for SPR moves
// ===========================================================================

// ---------------------------------------------------------------------------
// SprMeta: metadata from an SPR proposal on TreeNav.
// ---------------------------------------------------------------------------
struct SprMeta {
  int u, v, p, w, a, b;   // node identities (see diagram below)
  double tau;              // split fraction for regraft edge
  double logHastings;      // log(lRegraft) - log(lMerge)
  bool valid;              // false if proposal failed (no eligible edges)

  // Saved edge lengths for reversal
  double oldLen_u;  // old nodeEdgeLen[u] = edge from u to old parent p
  double oldLen_w;  // old nodeEdgeLen[w] = edge from w to old parent u
  double oldLen_b;  // old nodeEdgeLen[b] = edge from b to old parent a

  // Saved child slots for reversal
  int oldSlot_p_u;   // which slot of p held u (-1 = ch0, etc.)
  int oldSlot_a_b;   // which slot of a held b
  int oldSlot_u_w;   // which slot of u held w
};

// ---------------------------------------------------------------------------
// Helper: replace child `old` with `rep` in node p's child slots.
// ---------------------------------------------------------------------------
static void replace_child(TreeNav& topo, int p, int oldChild, int newChild) {
  if      (topo.ch0[p] == oldChild) topo.ch0[p] = newChild;
  else if (topo.ch1[p] == oldChild) topo.ch1[p] = newChild;
  else if (topo.ch2[p] == oldChild) topo.ch2[p] = newChild;
}

// ---------------------------------------------------------------------------
// Helper: check if `target` is a descendant of `anc` in TreeNav.
// Uses parent-pointer upward walk (O(depth)).
// ---------------------------------------------------------------------------
static bool is_descendant_of(const TreeNav& topo, int target, int anc) {
  int cur = target;
  while (cur >= 0) {
    if (cur == anc) return true;
    cur = topo.parentNode[cur];
  }
  return false;
}

// ---------------------------------------------------------------------------
// propose_spr_treenav: select a random SPR move using TreeNav navigation.
//
// Diagram:
//   Before: ..→p→u→{v, w}   ..→a→{b, ..}
//   After:  ..→p→{w, ..}    ..→a→u→{v, b}
//
// u = node being re-grafted (parent of subtree v)
// v = root of pruned subtree (stays as child of u)
// p = old parent of u (gains w, loses u)
// w = sibling of v under u (promoted to p's child)
// a = parent end of regraft edge
// b = child end of regraft edge (becomes u's new child)
// ---------------------------------------------------------------------------
static SprMeta propose_spr_treenav(const TreeNav& topo) {
  SprMeta m;
  m.valid = false;
  m.logHastings = R_NegInf;

  int nTip = topo.nTip;
  int root = topo.root;

  // 1. Select random prune node v: any node whose parent u ≠ root
  //    (i.e. u has a parent p)
  std::vector<int> eligible;
  eligible.reserve(topo.nEdge);
  for (int node = 1; node <= topo.maxNode; ++node) {
    int par = topo.parentNode[node];
    if (par < 0) continue;       // root or unused
    if (par == root) continue;    // parent is root → can't prune
    eligible.push_back(node);
  }
  if (eligible.empty()) return m;

  int pick = (int)(R::unif_rand() * (double)eligible.size());
  if (pick >= (int)eligible.size()) pick = eligible.size() - 1;
  m.v = eligible[pick];
  m.u = topo.parentNode[m.v];
  m.p = topo.parentNode[m.u];

  // Find sibling w of v under u
  int ch[3] = { topo.ch0[m.u], topo.ch1[m.u], topo.ch2[m.u] };
  m.w = -1;
  for (int ci = 0; ci < 3; ++ci) {
    if (ch[ci] >= 0 && ch[ci] != m.v) { m.w = ch[ci]; break; }
  }
  if (m.w < 0) return m;  // shouldn't happen for valid binary tree

  // 2. Enumerate candidate regraft edges.
  //    Valid: node b where parentNode[b] = a, b is not in v's subtree,
  //    and b is not u and a is not u (not adjacent to u).
  std::vector<int> candidates;
  candidates.reserve(topo.nEdge);
  for (int node = 1; node <= topo.maxNode; ++node) {
    if (topo.parentNode[node] < 0) continue;  // root
    if (node == m.u) continue;                  // adjacent to u
    if (topo.parentNode[node] == m.u) continue; // adjacent to u (sibling of v or w)
    if (is_descendant_of(topo, node, m.v)) continue; // in v's subtree
    candidates.push_back(node);
  }
  if (candidates.empty()) return m;

  int pickR = (int)(R::unif_rand() * (double)candidates.size());
  if (pickR >= (int)candidates.size()) pickR = candidates.size() - 1;
  m.b = candidates[pickR];
  m.a = topo.parentNode[m.b];

  m.tau = R::unif_rand();

  // Hastings ratio: log(lRegraft) - log(lMerge)
  double lRegraft = topo.nodeEdgeLen[m.b];  // length of edge a→b
  double lMerge   = topo.nodeEdgeLen[m.u] + topo.nodeEdgeLen[m.w]; // p→u + u→w
  if (lRegraft <= 0.0 || lMerge <= 0.0) return m;
  m.logHastings = std::log(lRegraft) - std::log(lMerge);

  m.valid = true;
  return m;
}


// ---------------------------------------------------------------------------
// update_topo_spr: apply SPR topology change to TreeNav.
// Saves old state in SprMeta for reversal.
// ---------------------------------------------------------------------------
static void update_topo_spr(TreeNav& topo, SprMeta& m) {
  // Save old edge lengths
  m.oldLen_u = topo.nodeEdgeLen[m.u];
  m.oldLen_w = topo.nodeEdgeLen[m.w];
  m.oldLen_b = topo.nodeEdgeLen[m.b];

  double lRegraft = m.oldLen_b;

  // 1. Detach u from p: replace u with w in p's children
  replace_child(topo, m.p, m.u, m.w);
  topo.parentNode[m.w] = m.p;
  topo.nodeEdgeLen[m.w] = m.oldLen_u + m.oldLen_w;  // merged edge

  // 2. Insert u on regraft edge: replace b with u in a's children
  replace_child(topo, m.a, m.b, m.u);
  topo.parentNode[m.u] = m.a;
  topo.nodeEdgeLen[m.u] = m.tau * lRegraft;  // a→u

  // 3. b becomes child of u (replacing w)
  replace_child(topo, m.u, m.w, m.b);
  topo.parentNode[m.b] = m.u;
  topo.nodeEdgeLen[m.b] = (1.0 - m.tau) * lRegraft;  // u→b
}


// ---------------------------------------------------------------------------
// reverse_topo_spr: undo SPR topology change on rejection.
// ---------------------------------------------------------------------------
static void reverse_topo_spr(TreeNav& topo, const SprMeta& m) {
  // Reverse step 3: replace b with w in u's children
  replace_child(topo, m.u, m.b, m.w);
  topo.parentNode[m.b] = m.a;   // b's parent back to a
  topo.nodeEdgeLen[m.b] = m.oldLen_b;

  // Reverse step 2: replace u with b in a's children
  replace_child(topo, m.a, m.u, m.b);
  topo.parentNode[m.u] = m.p;   // u's parent back to p
  topo.nodeEdgeLen[m.u] = m.oldLen_u;

  // Reverse step 1: replace w with u in p's children
  replace_child(topo, m.p, m.w, m.u);
  topo.parentNode[m.w] = m.u;   // w's parent back to u
  topo.nodeEdgeLen[m.w] = m.oldLen_w;
}


// ---------------------------------------------------------------------------
// find_dirty_spr: identify dirty nodes after SPR.
// Dirty set = union of:
//   - node u (children changed)
//   - path from p (old parent) to root
//   - path from a (new parent) to root
// Returned in postorder (deepest first).
// ---------------------------------------------------------------------------
static std::vector<int> find_dirty_spr(
    const TreeNav& topo, int u, int p, int a) {

  // Mark all nodes on both paths to root + u itself
  std::vector<bool> onPath(topo.maxNode + 1, false);
  onPath[u] = true;
  for (int n = p; n >= 0; n = topo.parentNode[n]) onPath[n] = true;
  for (int n = a; n >= 0; n = topo.parentNode[n]) onPath[n] = true;

  // Collect with depth for postorder sorting
  std::vector<std::pair<int,int>> depthNode;
  for (int n = 1; n <= topo.maxNode; ++n) {
    if (!onPath[n]) continue;
    int d = 0;
    for (int x = n; x >= 0; x = topo.parentNode[x]) ++d;
    depthNode.push_back({d, n});
  }

  // Sort by decreasing depth = postorder (children before parents)
  std::sort(depthNode.begin(), depthNode.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });

  std::vector<int> dirty;
  dirty.reserve(depthNode.size());
  for (auto& dn : depthNode) dirty.push_back(dn.second);
  return dirty;
}


// ---------------------------------------------------------------------------
// treenav_to_preorder: reconstruct canonical parent/child/relBr arrays from
// TreeNav. Uses DFS preorder traversal.
//
// Returns false if reconstruction fails (shouldn't happen with valid TreeNav).
// Updates topo.edgeToPar and topo.edgeLen to match the new array ordering.
// ---------------------------------------------------------------------------
static bool treenav_to_preorder(
    TreeNav& topo, double treeLength,
    IntegerVector& outParent, IntegerVector& outChild,
    NumericVector& outRelBr) {

  int nEdge = topo.nEdge;
  outParent = IntegerVector(nEdge);
  outChild  = IntegerVector(nEdge);
  outRelBr  = NumericVector(nEdge);

  // DFS preorder: push root's children, visit in preorder
  std::vector<std::pair<int,int>> stack;  // (parent, child)
  // Push root's children in reverse order for correct preorder
  int rootCh[3] = { topo.ch2[topo.root], topo.ch1[topo.root], topo.ch0[topo.root] };
  for (int ci = 0; ci < 3; ++ci) {
    if (rootCh[ci] >= 0) stack.push_back({topo.root, rootCh[ci]});
  }

  int idx = 0;
  while (!stack.empty()) {
    auto [par, ch] = stack.back();
    stack.pop_back();

    if (idx >= nEdge) return false;  // shouldn't happen

    outParent[idx] = par;
    outChild[idx]  = ch;
    outRelBr[idx]  = topo.nodeEdgeLen[ch] / treeLength;

    // Update TreeNav edge mapping to match new array ordering
    topo.edgeToPar[ch] = idx;
    topo.edgeLen[idx]  = topo.nodeEdgeLen[ch];

    // Push children in reverse order (so first child is visited first)
    int childCh[3] = { topo.ch2[ch], topo.ch1[ch], topo.ch0[ch] };
    for (int ci = 0; ci < 3; ++ci) {
      if (childCh[ci] >= 0) stack.push_back({ch, childCh[ci]});
    }
    ++idx;
  }

  return (idx == nEdge);
}


#endif  // MKPRIME_NODE_CL_CACHE_H
