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

  // T-013: rollback scratch for partial-eval NNI.  Holds saved CLs at dirty
  // nodes + the NNI swap parameters so the reject path can both restore
  // node CLs and reverse the TreeNav update.  Reused across iterations to
  // avoid per-call alloc.  Defined inline (EcoDirtyScratch struct appears
  // later in this header) via a forward-declared pointer-like wrapper:
  // we keep the buffers separate here and the struct provides typed access.
  std::vector<double> rollbackSavedCL;
  std::vector<int>    rollbackDirtyNodes;
  int rollbackNniV = -1, rollbackNniU = -1;
  int rollbackNniC = -1, rollbackNniW = -1;
  bool rollbackTopoUpdated = false;

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


// ===========================================================================
// T-012: cached ecology pruning helpers
// ===========================================================================
//
// The cached pruner stores per-(unit, cat, node) post-mixture CLs in
// `unit.cl`.  Each edge contribution under the ecology mixture is a
// kEco-weighted blend of per-ecology JC or MkN transition probabilities;
// the blend depends on the per-character z-row, the per-edge wEdge weights,
// and the gammaE / phi scalars.  Because the mixture is per-character, the
// transition is computed per-character at each edge — there is no "single
// matrix per edge" we can precompute.  The cache value is therefore the
// per-character CL after multiplying in all child contributions.
//
// Bit-identity vs the legacy flat pruner is preserved by:
//   * iterating ecology states `s` in ascending order (matches legacy)
//   * iterating rate categories `cat` outside the node loop (matches legacy)
//   * accumulating each unit's site-likelihoods over cats then averaging
//     with `1 / nCat` (matches legacy nCat normalisation)
//   * processing children in topo.ch0 / ch1 / ch2 order
//
// The cache is single-precision-equivalent: we operate in `double` throughout.

#include "gibbs_z_workspace.h"  // GibbsZWorkspace — not needed but mirrors blind path

// Forward declaration for the legacy const-site helper from mcmc_ecology.cpp.
// We re-use it for the per-character ascertainment correction (per-character
// pseudo-pruning).  Bit-identical to legacy.  Defined as a static inside
// mcmc_ecology.cpp; declared here as the inline helpers we add in the same
// translation unit can call it directly.


// ---------------------------------------------------------------------------
// eco_apply_mixture_jc: compute one (edge, cat) contribution under the JC-K
// ecology mixture and overlay it into `dst` (per-character buffer).
//
// For each character c, the per-state contribution from child CL `srcCL` to
// the parent is:
//   contrib[i] = pdMix(c) * sum_clCh(c) + (psMix(c) - pdMix(c)) * srcCL[c*k+i]
// where psMix/pdMix are the kEco-weighted sums of (ps/pd)Factor[z, s] under
// the character's z-row and the edge's wEdge row.
//
// If `firstChild` is true, the contribution is written; otherwise it is
// multiplied in (Felsenstein accumulation).
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_apply_mixture_jc(
    const EcoCacheUnit& unit,
    const double* srcCL, double* dst,
    int e, int nEdge, int kEco,
    int refEcology, int mode,
    const double* psFactor, const double* pdFactor,
    const double* wPtr, const int* zPtr,
    bool firstChild
) {
  int k = unit.kStates;
  int nChar = unit.nChar;
  // Special-case tip child detection isn't applied here; the caller knows
  // when src is a tip CL and can use the tip-edge fast path if desired.

  if (firstChild) {
    for (int c = 0; c < nChar; ++c) {
      double psMix = 0.0, pdMix = 0.0;
      for (int s = 0; s < kEco; ++s) {
        int z;
        if (s == refEcology) z = 0;
        else {
          int zCol = (s < refEcology) ? s : (s - 1);
          z = zPtr[c + zCol * nChar];
        }
        double w = wPtr[e + s * nEdge];
        psMix += w * psFactor[z * kEco + s];
        pdMix += w * pdFactor[z * kEco + s];
      }
      double diff_coeff = psMix - pdMix;
      int offset = c * k;
      double sum_cl = 0.0;
      for (int j = 0; j < k; ++j) sum_cl += srcCL[offset + j];
      for (int i = 0; i < k; ++i)
        dst[offset + i] = pdMix * sum_cl + diff_coeff * srcCL[offset + i];
    }
  } else {
    for (int c = 0; c < nChar; ++c) {
      double psMix = 0.0, pdMix = 0.0;
      for (int s = 0; s < kEco; ++s) {
        int z;
        if (s == refEcology) z = 0;
        else {
          int zCol = (s < refEcology) ? s : (s - 1);
          z = zPtr[c + zCol * nChar];
        }
        double w = wPtr[e + s * nEdge];
        psMix += w * psFactor[z * kEco + s];
        pdMix += w * pdFactor[z * kEco + s];
      }
      double diff_coeff = psMix - pdMix;
      int offset = c * k;
      double sum_cl = 0.0;
      for (int j = 0; j < k; ++j) sum_cl += srcCL[offset + j];
      for (int i = 0; i < k; ++i)
        dst[offset + i] *= pdMix * sum_cl + diff_coeff * srcCL[offset + i];
    }
  }
}

// ---------------------------------------------------------------------------
// eco_apply_mixture_mkn: same as eco_apply_mixture_jc but for the
// neomorphic 2-state mixture (full 2x2 matvec per character).
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_apply_mixture_mkn(
    const EcoCacheUnit& unit,
    const double* srcCL, double* dst,
    int e, int nEdge, int kEco,
    int refEcology, int /*mode*/,
    const double* Pfactor,  // size 3 * kEco * 4
    const double* wPtr, const int* zPtr,
    bool firstChild
) {
  (void)unit;
  int nChar = unit.nChar;
  if (firstChild) {
    for (int c = 0; c < nChar; ++c) {
      double P00m = 0, P01m = 0, P10m = 0, P11m = 0;
      for (int s = 0; s < kEco; ++s) {
        int z;
        if (s == refEcology) z = 0;
        else {
          int zCol = (s < refEcology) ? s : (s - 1);
          z = zPtr[c + zCol * nChar];
        }
        double w = wPtr[e + s * nEdge];
        const double* P = Pfactor + (z * kEco + s) * 4;
        P00m += w * P[0]; P01m += w * P[1];
        P10m += w * P[2]; P11m += w * P[3];
      }
      int off = c * 2;
      double cl0 = srcCL[off], cl1 = srcCL[off + 1];
      dst[off]     = P00m * cl0 + P01m * cl1;
      dst[off + 1] = P10m * cl0 + P11m * cl1;
    }
  } else {
    for (int c = 0; c < nChar; ++c) {
      double P00m = 0, P01m = 0, P10m = 0, P11m = 0;
      for (int s = 0; s < kEco; ++s) {
        int z;
        if (s == refEcology) z = 0;
        else {
          int zCol = (s < refEcology) ? s : (s - 1);
          z = zPtr[c + zCol * nChar];
        }
        double w = wPtr[e + s * nEdge];
        const double* P = Pfactor + (z * kEco + s) * 4;
        P00m += w * P[0]; P01m += w * P[1];
        P10m += w * P[2]; P11m += w * P[3];
      }
      int off = c * 2;
      double cl0 = srcCL[off], cl1 = srcCL[off + 1];
      dst[off]     *= P00m * cl0 + P01m * cl1;
      dst[off + 1] *= P10m * cl0 + P11m * cl1;
    }
  }
}

// ---------------------------------------------------------------------------
// eco_compute_jc_factors: fill psFactor / pdFactor (length 3*kEco each) for
// a single (cat, edge) cell.  Mirrors the corresponding block in
// pruning_jc_acrv_flat_ecology, so the values are bit-identical.
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_compute_jc_factors(
    int kStates, int kEco, int refEcology, int mode,
    const double* phiPtr,
    const std::vector<double>& gammaE,
    double tBase,
    double* psFactor, double* pdFactor
) {
  double inv_k = 1.0 / kStates;
  double km1 = (double)kStates - 1.0;
  for (int s = 0; s < kEco; ++s) {
    double gE = gammaE[s];
    for (int z = 0; z < 3; ++z) {
      double factor;
      if (s == refEcology) factor = 1.0;
      else {
        double p = (mode == 0) ? phiPtr[0] : phiPtr[s];
        double mu = (z == 0) ? 1.0 : (z == 1) ? p : 1.0 / p;
        factor = mu / gE;
      }
      double tEff = tBase * factor;
      double exV = MKP_EXP(-kStates * tEff / km1);
      psFactor[z * kEco + s] = inv_k + (1.0 - inv_k) * exV;
      pdFactor[z * kEco + s] = inv_k - inv_k * exV;
    }
  }
}

// ---------------------------------------------------------------------------
// eco_compute_mkn_factors: fill Pfactor (length 3*kEco*4) for one (cat, edge)
// cell.  Mirrors the corresponding block in pruning_mkn_acrv_flat_ecology.
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_compute_mkn_factors(
    int kEco, int refEcology, int mode,
    const double* phiPtr,
    const std::vector<double>& gammaE,
    double rateLoss,
    double tBase,
    double* Pfactor
) {
  double sum_rl = 1.0 + rateLoss;
  double r01_base = 2.0 / sum_rl;
  double r10_base = 2.0 * rateLoss / sum_rl;
  for (int s = 0; s < kEco; ++s) {
    double gE = gammaE[s];
    for (int z = 0; z < 3; ++z) {
      double mu01, mu10;
      if (s == refEcology) {
        mu01 = 1.0; mu10 = 1.0;
      } else {
        double p = (mode == 0) ? phiPtr[0] : phiPtr[s];
        if (z == 0) { mu01 = 1.0; mu10 = 1.0; }
        else if (z == 1) { mu01 = p; mu10 = 1.0 / p; }
        else { mu01 = 1.0 / p; mu10 = p; }
      }
      double r01, r10;
      if (s == refEcology) { r01 = r01_base; r10 = r10_base; }
      else { r01 = r01_base * mu01 / gE; r10 = r10_base * mu10 / gE; }
      double lam = r01 + r10;
      double pi0_ = r10 / lam;
      double pi1_ = r01 / lam;
      double ex = MKP_EXP(-lam * tBase);
      double P00 = pi0_ + pi1_ * ex;
      double P01 = pi1_ - pi1_ * ex;
      double P10 = pi0_ - pi0_ * ex;
      double P11 = pi1_ + pi0_ * ex;
      double* P = Pfactor + (z * kEco + s) * 4;
      P[0] = P00; P[1] = P01; P[2] = P10; P[3] = P11;
    }
  }
}


// ---------------------------------------------------------------------------
// eco_full_downpass_unit: write CLs for one EcoCacheUnit by running a full
// postorder downpass under the current wEdge / params.  Mirrors the legacy
// flat pruner edge-iteration order so the root CLs match bit-for-bit.
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_full_downpass_unit(
    EcoCacheUnit& unit,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen, int nTip,
    int kEco, int refEcology, int mode,
    const NumericVector& phi,
    const NumericMatrix& wEdge,
    const IntegerMatrix& zMat,  // sub-zMat aligned with unit chars
    const std::vector<double>& gammaE,
    const std::vector<double>& rates,
    double rateLoss
) {
  int nEdge = parent.size();
  int nCat = (int)rates.size();
  int maxNode = 2 * nTip - 1;
  const int* parPtr = INTEGER(parent);
  const int* chPtr  = INTEGER(child);
  const double* elPtr = REAL(edgeLen);
  const double* wPtr  = REAL(wEdge);
  const int* zPtr     = INTEGER(zMat);
  const double* phiPtr = REAL(phi);

  std::vector<uint8_t> initFlg(maxNode + 1, 0);
  std::vector<double> psFactor(3 * kEco), pdFactor(3 * kEco);
  std::vector<double> Pfactor;
  if (unit.isMkN) Pfactor.assign(3 * kEco * 4, 0.0);

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];
    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;
    for (int t = 1; t <= nTip; ++t) initFlg[t] = 1;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parPtr[e], ch = chPtr[e];
      double tBase = elPtr[e] * rate * unit.rateScale;

      if (unit.isMkN) {
        eco_compute_mkn_factors(kEco, refEcology, mode, phiPtr, gammaE,
                                rateLoss, tBase, Pfactor.data());
      } else {
        eco_compute_jc_factors(unit.kStates, kEco, refEcology, mode, phiPtr,
                               gammaE, tBase, psFactor.data(), pdFactor.data());
      }

      double* dst = unit.CL(cat, par);
      const double* src = unit.CL(cat, ch);
      bool firstChild = !initFlg[par];

      if (unit.isMkN) {
        eco_apply_mixture_mkn(unit, src, dst, e, nEdge, kEco,
                              refEcology, mode, Pfactor.data(),
                              wPtr, zPtr, firstChild);
      } else {
        eco_apply_mixture_jc(unit, src, dst, e, nEdge, kEco,
                             refEcology, mode,
                             psFactor.data(), pdFactor.data(),
                             wPtr, zPtr, firstChild);
      }
      initFlg[par] = 1;
    }
  }
}

// ---------------------------------------------------------------------------
// eco_recompute_dirty_unit: rebuild CLs only at `dirtyNodes` (in postorder)
// for one EcoCacheUnit using TreeNav child pointers.
//
// Each dirty node `node`'s CL is rebuilt by:
//   for cat in 0..nCat-1:
//     for each child slot in (ch0, ch1, ch2):
//       compute contribution from child CL (which is current — either it
//       was clean or already updated earlier in postorder)
//       first child: write; subsequent: multiply
//
// The per-edge wEdge / edgeLen come from the caller's *current* arrays
// (proposed topology for an NNI move).
// ---------------------------------------------------------------------------
[[maybe_unused]] static void eco_recompute_dirty_unit(
    EcoCacheUnit& unit,
    const EcoCLCache& cache,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    int kEco, int refEcology, int mode,
    const NumericVector& phi,
    const NumericMatrix& wEdge,
    const IntegerMatrix& zMat,
    const std::vector<double>& gammaE,
    const std::vector<double>& rates,
    double rateLoss,
    const std::vector<int>& dirtyNodes
) {
  int nEdge = parent.size();
  int nCat = (int)rates.size();
  const int* parPtr = INTEGER(parent);
  const int* chPtr  = INTEGER(child);
  const double* elPtr = REAL(edgeLen);
  const double* wPtr  = REAL(wEdge);
  const int* zPtr    = INTEGER(zMat);
  const double* phiPtr = REAL(phi);

  std::vector<double> psFactor(3 * kEco), pdFactor(3 * kEco);
  std::vector<double> Pfactor;
  if (unit.isMkN) Pfactor.assign(3 * kEco * 4, 0.0);

  // Map child node -> edge index (so we can grab edge length & wEdge row).
  // Build once per call.  Cheap (O(nEdge)).
  std::vector<int> childEdgeOf(cache.maxNode + 1, -1);
  for (int e = 0; e < nEdge; ++e) childEdgeOf[chPtr[e]] = e;
  (void)parPtr;  // unused; we use topo's parent pointers

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    for (int node : dirtyNodes) {
      if (node <= cache.nTip) continue;
      double* dst = unit.CL(cat, node);
      bool first = true;

      int children[3] = { cache.topo.ch0[node],
                          cache.topo.ch1[node],
                          cache.topo.ch2[node] };
      for (int ci = 0; ci < 3; ++ci) {
        int ch = children[ci];
        if (ch < 0) continue;
        int e = childEdgeOf[ch];
        if (e < 0) continue;  // shouldn't happen
        double tBase = elPtr[e] * rate * unit.rateScale;

        if (unit.isMkN) {
          eco_compute_mkn_factors(kEco, refEcology, mode, phiPtr, gammaE,
                                  rateLoss, tBase, Pfactor.data());
        } else {
          eco_compute_jc_factors(unit.kStates, kEco, refEcology, mode, phiPtr,
                                 gammaE, tBase, psFactor.data(), pdFactor.data());
        }

        const double* src = unit.CL(cat, ch);

        if (unit.isMkN) {
          eco_apply_mixture_mkn(unit, src, dst, e, nEdge, kEco,
                                refEcology, mode, Pfactor.data(),
                                wPtr, zPtr, first);
        } else {
          eco_apply_mixture_jc(unit, src, dst, e, nEdge, kEco,
                               refEcology, mode,
                               psFactor.data(), pdFactor.data(),
                               wPtr, zPtr, first);
        }
        first = false;
      }
    }
  }
}


// ---------------------------------------------------------------------------
// eco_unit_root_loglik: pruning log-lik contribution for one EcoCacheUnit
// (no ascertainment correction, no relabel — those are applied externally).
// Bit-identical to the legacy summation order.
// ---------------------------------------------------------------------------
[[maybe_unused]] static double eco_unit_root_loglik(
    const EcoCacheUnit& unit, int root, int nCat
) {
  int k = unit.kStates;
  int nChar = unit.nChar;
  const double* rf = unit.rootFreqs.data();

  if (nCat == 1) {
    const double* clRoot = unit.CL(0, root);
    double logLik = 0.0;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * k;
      double sl = 0.0;
      for (int s = 0; s < k; ++s) sl += rf[s] * clRoot[offset + s];
      if (sl <= 0.0) return R_NegInf;
      logLik += std::log(sl);
    }
    return logLik;
  }
  std::vector<double> siteLikSum(nChar, 0.0);
  for (int cat = 0; cat < nCat; ++cat) {
    const double* clRoot = unit.CL(cat, root);
    for (int c = 0; c < nChar; ++c) {
      int offset = c * k;
      double sl = 0.0;
      for (int s = 0; s < k; ++s) sl += rf[s] * clRoot[offset + s];
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


// ---------------------------------------------------------------------------
// find_dirty_nni_eco: same dirty-node walk as find_dirty_nni in
// node_cl_cache.h, but using EcoCLCache.topo (which is updated symmetrically
// after the proposal).  Returns nodes in postorder (v deepest first, root
// last).  When wEdge changes affect non-path nodes, the caller unions the
// result with `cache.invalidate_ancestors_of_dirty_edges`-derived nodes.
// ---------------------------------------------------------------------------
[[maybe_unused]] static std::vector<int> find_dirty_nni_eco(
    const TreeNav& topo, int v, int u
) {
  std::vector<int> path;
  path.reserve(16);
  path.push_back(v);
  int cur = u;
  while (cur >= 0) {
    path.push_back(cur);
    cur = topo.parentNode[cur];
  }
  return path;
}


// ---------------------------------------------------------------------------
// save_dirty_eco_cls / restore_dirty_eco_cls: M-158-style backup + restore
// of CL data for all (unit, cat, dirtyNode) cells.
// ---------------------------------------------------------------------------
struct EcoDirtyScratch {
  std::vector<double> savedCL;
  std::vector<int> dirtyNodes;
  // T-013: NNI rollback also restores TreeNav.  Track the NNI nodes so the
  // reject path knows which (v, u, cNode, wNode) tuple to reverse.
  int nniV = -1, nniU = -1, nniC = -1, nniW = -1;
  bool topoUpdated = false;
};

[[maybe_unused]] static void save_dirty_eco_cls(
    EcoCLCache& cache,
    EcoDirtyScratch& scratch,
    const std::vector<int>& dirtyNodes
) {
  size_t total = 0;
  for (const auto& u : cache.units)
    total += (size_t)u.nCat * dirtyNodes.size() * u.stride;
  scratch.savedCL.resize(total);
  scratch.dirtyNodes = dirtyNodes;
  double* dst = scratch.savedCL.data();
  for (const auto& u : cache.units) {
    for (int cat = 0; cat < u.nCat; ++cat) {
      for (int node : dirtyNodes) {
        const double* src = u.CL(cat, node);
        std::memcpy(dst, src, u.stride * sizeof(double));
        dst += u.stride;
      }
    }
  }
}

[[maybe_unused]] static void restore_dirty_eco_cls(
    EcoCLCache& cache,
    const EcoDirtyScratch& scratch
) {
  const double* src = scratch.savedCL.data();
  for (auto& u : cache.units) {
    for (int cat = 0; cat < u.nCat; ++cat) {
      for (int node : scratch.dirtyNodes) {
        double* dst = u.CL(cat, node);
        std::memcpy(dst, src, u.stride * sizeof(double));
        src += u.stride;
      }
    }
  }
}


#endif  // MKPRIME_ECOLOGY_CL_CACHE_H
