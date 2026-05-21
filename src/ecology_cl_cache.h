#ifndef MKPRIME_ECOLOGY_CL_CACHE_H
#define MKPRIME_ECOLOGY_CL_CACHE_H

// T-011: persistent per-node CL cache for the ecology-aware pruner.
//
// **STATUS (T-011 round, 2026-05-21): DESIGN SKETCH, NOT WIRED.**
//
// This header defines the cache layout, invalidation API, and helpers that a
// future "Lever 1 + Lever 2 lite" implementation would consume.  In the
// T-011 round only the structure and the invalidation calls land — the cache
// is never populated and the cached pruner variants are deferred.  The
// pre-implementation advisor consultation predicted (and the bench confirmed)
// that simple-invalidation Lever 1 delivers 0 % wall improvement on the
// rodent workload because the T-010 partition cache already captures every
// non-tree-move saving, and tree moves invalidate the whole cache via
// `invalidate_all()`.  Wall recovery requires per-edge wEdge dirty detection
// + an M-158-style save / restore around tree-move evaluations.  See
// `dev/profiling/findings.md` T-011 for the full rationale.
//
// Patterned after node_cl_cache.h (M-121 / M-158 / M-161), adapted for the
// ecology-mixture pruner.  The cache stores per-(partition-subgroup, cat,
// node) "mixed" CLs (each edge contribution is the kEco-weighted mixture of
// per-ecology-state transition probabilities, so the CL itself already folds
// in wEdge and z).  Validity tracking is per-node (one flag covers all cats).
//
// The cache exists per chain inside McmcState (parallel to nodeCL for the
// blind path).  Invalidation rules:
//
//   invalidate_all()       — topology/structure unknown; full rebuild needed
//   invalidate_structure() — kPrime composition shifted; rebuild subgroups
//   invalidate_all_cls()   — every node's CL stale (e.g. phi/pi0/theta/rateLogSd
//                            change, or wEdge dirty without per-edge tracking)
//   invalidate_neo_cls()   — only neomorphic unit CLs stale (rate_loss/rate_neo)
//
// On a future eval, the cache will store the wEdge / edgeLen snapshot it was
// populated against; later calls can use these to detect per-edge changes
// (Lever 2 lite) and invalidate only the affected ancestor nodes.  This is
// the only avenue for partial-CL savings on tree moves.

#include "mcmc_state.h"
#include "gibbs_partial_cl.h"  // TreeNav
#include "fast_exp.h"
#include <cstring>
#include <cmath>
#include <algorithm>
#include <map>

using namespace Rcpp;

// Forward decl from mcmc_likelihood.cpp.
double mk_prime_relabel_log(int kPrime, int kObs);

// ---------------------------------------------------------------------------
// EcoCacheUnit: per-(partition × kPrime-subgroup) per-node CL storage
// ---------------------------------------------------------------------------
struct EcoCacheUnit {
  int partIdx = -1;       // index into data.parts[]
  int kStates = 0;        // post-mixing CL state count
  int nChar   = 0;        // chars in this unit
  int nCat    = 1;        // ACRV categories
  int stride  = 0;        // nChar × kStates
  bool isMkN  = false;    // true → neomorphic (2-state asymmetric)
  double rateScale = 1.0; // rateNeo for neo, 1.0 otherwise
  bool doRelabel = false;

  // Tip data for this unit (nTip × nChar, 0-indexed, -1 = missing)
  IntegerMatrix tipStates;
  // Per-character kObs (transformational only, used for Mk' relabel correction)
  std::vector<int> kObsLocal;
  // 0-based global char indices (for ascertainment / per-character corrections)
  std::vector<int> globalCharIdx;
  // Root frequencies (length kStates)
  std::vector<double> rootFreqs;

  // CL storage: cl[(cat * (maxNode+1) + node) * stride]
  std::vector<double> cl;
  int maxNode = 0;

  double* CL(int cat, int node) {
    return cl.data() + ((size_t)cat * (maxNode + 1) + node) * stride;
  }
  const double* CL(int cat, int node) const {
    return cl.data() + ((size_t)cat * (maxNode + 1) + node) * stride;
  }

  void allocate(int maxNode_, int nCat_) {
    maxNode = maxNode_;
    nCat    = nCat_;
    stride  = nChar * kStates;
    cl.assign((size_t)nCat * (maxNode + 1) * stride, 0.0);
  }

  // Initialise tip CLs (constant across cats and across iterations).
  // Layout: clTip[c * kStates + state] = 1, rest 0 (or all 1 if missing).
  void init_tips(int nTip) {
    if (cl.empty()) return;
    for (int tip = 1; tip <= nTip; ++tip) {
      double* clT = CL(0, tip);
      std::fill(clT, clT + stride, 0.0);
      for (int c = 0; c < nChar; ++c) {
        int state = tipStates(tip - 1, c);
        int offset = c * kStates;
        if (state < 0) {
          for (int s = 0; s < kStates; ++s) clT[offset + s] = 1.0;
        } else {
          clT[offset + state] = 1.0;
        }
      }
      for (int cat = 1; cat < nCat; ++cat) {
        std::memcpy(CL(cat, tip), clT, stride * sizeof(double));
      }
    }
  }
};

// ---------------------------------------------------------------------------
// EcoCLCache: persistent cache held in McmcState
// ---------------------------------------------------------------------------
struct EcoCLCache {
  std::vector<EcoCacheUnit> units;
  // Per-partition-index lookup → indices into units[] (one trans partition can
  // span multiple subgroup units if kPrime is heterogeneous).
  std::vector<std::vector<int>> unitsByPart;

  // Persistent tree navigation
  TreeNav topo;
  int maxNode = 0;
  int nTip = 0;
  bool useAcrv = false;
  std::vector<double> rates;
  double cachedRateLogSd = -1.0;

  // wEdge / edgeLen snapshot used at the last full populate.  Future calls
  // compare against current to detect per-edge changes (Lever 2 lite).
  NumericMatrix cachedWEdge;
  NumericVector cachedEdgeLen;

  // Scalar param snapshot — change in any of these invalidates every CL.
  double cachedPi0 = -1.0;
  std::vector<double> cachedPhi;
  std::vector<double> cachedTheta;
  std::vector<double> cachedGammaE;

  // Per-node validity (single bool per node spans all cats).  nodeValid[n] == 1
  // means cl[*, n, *] is current with respect to the cached wEdge / params.
  std::vector<uint8_t> nodeValid;

  // Two-level validity gates.
  bool topoValid      = false;  // TreeNav matches current topology + edge lengths
  bool structureValid = false;  // unit composition matches kPrime grouping

  // Diagnostics
  int diagCachedSubtreeCount = 0;  // dirty nodes count, summed over evaluations
  int diagFullEvalCount = 0;
  int diagPartialEvalCount = 0;

  // ----- predicate helpers -----
  bool ready() const {
    return topoValid && structureValid;
  }

  // All units have a fully-valid set of node CLs.
  bool all_nodes_valid() const {
    if (!ready()) return false;
    for (int n = nTip + 1; n <= maxNode; ++n)
      if (!nodeValid[n]) return false;
    return true;
  }

  // ----- coarse invalidators -----
  void invalidate_all() {
    topoValid = false;
    structureValid = false;
    std::fill(nodeValid.begin(), nodeValid.end(), 0u);
  }

  void invalidate_structure() {
    structureValid = false;
    std::fill(nodeValid.begin(), nodeValid.end(), 0u);
  }

  // All node CLs stale (e.g. global parameter changed affecting factor table).
  // Structure intact.
  void invalidate_all_cls() {
    if (nodeValid.empty()) return;
    // Tips remain valid.
    for (int n = nTip + 1; n <= maxNode; ++n) nodeValid[n] = 0u;
  }

  // Only neomorphic units affected (rate_loss / rate_neo) — but our nodeValid
  // is global, so we invalidate it.  Future refinement could be per-unit.
  void invalidate_neo_cls() {
    invalidate_all_cls();
  }

  // ----- fine-grained invalidator (Lever 2 lite) -----
  //
  // Walk from each dirty edge's child upward, marking every ancestor invalid.
  // After this, `nodeValid` is correct for the new wEdge/edgeLen configuration.
  void invalidate_ancestors_of_dirty_edges(
      const std::vector<int>& dirtyChildNodes) {
    if (nodeValid.empty()) return;
    for (int ch : dirtyChildNodes) {
      // Mark the child itself stale (its outgoing edge has dirty factors —
      // but its CL is computed from ITS children, which are clean, so the
      // child is only stale if its EDGE TO PARENT changed.  Actually it's
      // the parent's CL that depends on the edge; the child's own CL is
      // computed from its own descendants).  Walk from parent upward.
      int p = topo.parentNode[ch];
      while (p >= 0) {
        if (!nodeValid[p]) break;  // already invalid → ancestors too
        nodeValid[p] = 0u;
        p = topo.parentNode[p];
      }
    }
  }

  // Compute a postorder list of currently-invalid internal nodes.  Used by
  // the partial-eval pruner.
  std::vector<int> dirty_internal_nodes_postorder() const {
    std::vector<std::pair<int,int>> depthNode;
    for (int n = nTip + 1; n <= maxNode; ++n) {
      if (nodeValid[n]) continue;
      int d = 0;
      for (int x = n; x >= 0; x = topo.parentNode[x]) ++d;
      depthNode.emplace_back(d, n);
    }
    std::sort(depthNode.begin(), depthNode.end(),
              [](const auto& a, const auto& b) { return a.first > b.first; });
    std::vector<int> out;
    out.reserve(depthNode.size());
    for (auto& dn : depthNode) out.push_back(dn.second);
    return out;
  }
};


// ---------------------------------------------------------------------------
// build_eco_cache_units: enumerate units for current kPrime grouping.
// Mirrors build_cache_units in node_cl_cache.h.  Does NOT allocate CL buffers
// (caller invokes allocate() per unit).
// ---------------------------------------------------------------------------
[[maybe_unused]] static void build_eco_cache_units(
    EcoCLCache& cache, const McmcData& data,
    const IntegerVector& kPrime) {

  cache.units.clear();
  int nParts = (int)data.parts.size();
  cache.unitsByPart.assign(nParts, {});

  for (int pi = 0; pi < nParts; ++pi) {
    const PartInfo& part = data.parts[pi];
    int nCharPart = part.tipStates.ncol();
    int nTipPart  = part.tipStates.nrow();

    if (part.type == 0) {
      // Neomorphic
      EcoCacheUnit u;
      u.partIdx = pi;
      u.kStates = 2;
      u.nChar   = nCharPart;
      u.isMkN   = true;
      u.tipStates = part.tipStates;
      u.doRelabel = false;
      u.globalCharIdx.assign(part.globalCharIdx.begin(),
                              part.globalCharIdx.end());
      // MkN stationary frequencies (set per evaluation via update_unit_params)
      u.rootFreqs = { 0.5, 0.5 };
      cache.unitsByPart[pi].push_back((int)cache.units.size());
      cache.units.push_back(std::move(u));

    } else if (part.type == 2) {
      // Known state space (single k)
      EcoCacheUnit u;
      u.partIdx = pi;
      u.kStates = part.k;
      u.nChar   = nCharPart;
      u.isMkN   = false;
      u.tipStates = part.tipStates;
      u.doRelabel = false;
      u.globalCharIdx.assign(part.globalCharIdx.begin(),
                              part.globalCharIdx.end());
      u.rootFreqs.assign(part.k, 1.0 / part.k);
      cache.unitsByPart[pi].push_back((int)cache.units.size());
      cache.units.push_back(std::move(u));

    } else {
      // Transformational — sub-group by kPrime
      std::map<int, std::vector<int>> byKp;
      for (int ci = 0; ci < nCharPart; ++ci) {
        int gi = part.globalCharIdx[ci];
        byKp[kPrime[gi]].push_back(ci);
      }
      for (auto& kv : byKp) {
        int kp = kv.first;
        const std::vector<int>& cols = kv.second;
        int nSub = (int)cols.size();

        EcoCacheUnit u;
        u.partIdx = pi;
        u.kStates = kp;
        u.nChar   = nSub;
        u.isMkN   = false;
        u.doRelabel = data.relabel;
        u.tipStates = IntegerMatrix(nTipPart, nSub);
        u.globalCharIdx.resize(nSub);
        u.kObsLocal.resize(nSub);
        for (int ci = 0; ci < nSub; ++ci) {
          for (int t = 0; t < nTipPart; ++t)
            u.tipStates(t, ci) = part.tipStates(t, cols[ci]);
          u.globalCharIdx[ci] = part.globalCharIdx[cols[ci]];
          u.kObsLocal[ci] = part.kObsLocal[cols[ci]];
        }
        u.rootFreqs.assign(kp, 1.0 / kp);
        cache.unitsByPart[pi].push_back((int)cache.units.size());
        cache.units.push_back(std::move(u));
      }
    }
  }
}


// ---------------------------------------------------------------------------
// detect_dirty_edges: scan cached vs current wEdge / edgeLen and return the
// list of "dirty" child nodes — children whose incoming edge has changed in
// either wEdge[e, :] or edgeLen[e].  Each such edge invalidates its parent's
// CL (and all ancestors of that parent).
//
// Comparison uses exact equality with a small absolute tolerance.  Tree moves
// regenerate wEdge from scratch so even tiny differences from FP rounding can
// appear; we use a conservative threshold (1e-12) to avoid spurious invalidations
// while still catching every real change.
// ---------------------------------------------------------------------------
[[maybe_unused]] static std::vector<int> detect_dirty_child_nodes(
    const EcoCLCache& cache,
    const IntegerVector& child,
    const NumericMatrix& newWEdge,
    const NumericVector& newEdgeLen) {

  std::vector<int> dirty;
  int nEdge = newEdgeLen.size();
  if (cache.cachedWEdge.nrow() != nEdge ||
      cache.cachedEdgeLen.size() != nEdge) {
    // No snapshot — every edge counts as dirty.
    dirty.reserve(nEdge);
    for (int e = 0; e < nEdge; ++e) dirty.push_back(child[e]);
    return dirty;
  }
  int kEco = newWEdge.ncol();
  const double TOL = 1e-12;
  for (int e = 0; e < nEdge; ++e) {
    bool dirtyEdge = false;
    if (std::abs(newEdgeLen[e] - cache.cachedEdgeLen[e]) > TOL) {
      dirtyEdge = true;
    } else {
      for (int s = 0; s < kEco; ++s) {
        if (std::abs(newWEdge(e, s) - cache.cachedWEdge(e, s)) > TOL) {
          dirtyEdge = true; break;
        }
      }
    }
    if (dirtyEdge) dirty.push_back(child[e]);
  }
  return dirty;
}


// ---------------------------------------------------------------------------
// snapshot_eco_inputs: copy current wEdge / edgeLen / params into the cache.
// Called after a successful full populate so subsequent change detection has
// a reference.
// ---------------------------------------------------------------------------
[[maybe_unused]] static void snapshot_eco_inputs(
    EcoCLCache& cache,
    const NumericMatrix& wEdge, const NumericVector& edgeLen,
    const NumericVector& phi, double pi0,
    const NumericVector& theta,
    const std::vector<double>& gammaE) {

  cache.cachedWEdge   = NumericMatrix(wEdge.nrow(), wEdge.ncol());
  std::memcpy(REAL(cache.cachedWEdge), REAL(wEdge),
              sizeof(double) * wEdge.nrow() * wEdge.ncol());
  cache.cachedEdgeLen = NumericVector(edgeLen.size());
  std::memcpy(REAL(cache.cachedEdgeLen), REAL(edgeLen),
              sizeof(double) * edgeLen.size());
  cache.cachedPhi.assign(phi.begin(), phi.end());
  cache.cachedPi0 = pi0;
  cache.cachedTheta.assign(theta.begin(), theta.end());
  cache.cachedGammaE = gammaE;
}


// ---------------------------------------------------------------------------
// eco_scalar_params_unchanged: cheap check that phi / pi0 / theta / gammaE
// match the snapshot.  If any differ, the per-edge factor table changes,
// so every cached node CL is stale.
// ---------------------------------------------------------------------------
[[maybe_unused]] static bool eco_scalar_params_unchanged(
    const EcoCLCache& cache,
    const NumericVector& phi, double pi0,
    const NumericVector& theta,
    const std::vector<double>& gammaE) {

  const double TOL = 1e-12;
  if (cache.cachedPhi.size() != (size_t)phi.size()) return false;
  for (int i = 0; i < phi.size(); ++i)
    if (std::abs(phi[i] - cache.cachedPhi[i]) > TOL) return false;
  if (std::abs(pi0 - cache.cachedPi0) > TOL) return false;
  if (cache.cachedTheta.size() != (size_t)theta.size()) return false;
  for (int i = 0; i < theta.size(); ++i)
    if (std::abs(theta[i] - cache.cachedTheta[i]) > TOL) return false;
  if (cache.cachedGammaE.size() != gammaE.size()) return false;
  for (size_t i = 0; i < gammaE.size(); ++i)
    if (std::abs(gammaE[i] - cache.cachedGammaE[i]) > TOL) return false;
  return true;
}


#endif  // MKPRIME_ECOLOGY_CL_CACHE_H
