// Ecology partition reconstruction for the ecology-aware NT model.
//
// The ecology-aware extension treats edge ecology as latent and marginalises
// over it when evaluating the character likelihood. The ecology covariate
// itself evolves on the tree under a standard JC-K Mk process: a forward
// (Felsenstein, postorder) pass plus a backward (preorder, upward-message)
// pass yields per-node posterior marginals P(node = e | tip ecology data).
//
// Per-edge weights for the rate-modifier mixture are then computed as
//
//   w_{p->c}(e) = 0.5 * (P(parent = e | data) + P(child = e | data)).
//
// Using the mean of parent and child marginals (rather than just the parent
// marginal, or a single-time-point midpoint) is a cheap symmetric blend that
// reuses both marginals coming out of the same backward sweep. It still
// collapses path uncertainty along the edge — a known limitation of any
// marginal-reconstruction approximation; data augmentation MCMC over edge
// ecology paths is the principled upgrade and is deferred to a later phase.
//
// Single ecology character, no ACRV: the partition is reconstructed at one
// rate. Multiple ecology axes are out of scope here.

#include "mcmc_state.h"
#include "gibbs_z_workspace.h"
#include "ecology_cl_cache.h"
#include "fast_exp.h"
#include <cmath>
#include <vector>
#include <map>
#include <cstdlib>     // getenv
#include <cstdio>      // fprintf for env-gated diag
#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;

// Forward declarations from mcmc_likelihood.cpp.
double mk_prime_relabel_log(int kPrime, int kObs);
NumericVector cpp_acrv_rates(double rateLogSd, int nCat,
                              const std::vector<double>& acrvZ);
static NumericVector mkn_stationary_local(double rateLoss) {
  NumericVector f(2);
  f[0] = rateLoss / (1.0 + rateLoss);
  f[1] = 1.0 / (1.0 + rateLoss);
  return f;
}


// Fill `marg` with per-node posterior marginals under JC-K.
// Layout: marg has size (maxNode + 1) * kStates, 1-indexed by node.
// Tips with observed states get one-hot rows; missing tips (-1) get the
// posterior predictive from the rest of the tree.
static void compute_ecology_node_marginals(
    const int* parent, const int* child, int nEdge,
    const double* edgeLen,
    const int* tipStates, int nTip,
    int kStates,
    double rateMultiplier,
    std::vector<double>& marg
) {
  int maxNode = 2 * nTip - 1;
  int root    = nTip + 1;
  double inv_k = 1.0 / kStates;
  double km1   = static_cast<double>(kStates) - 1.0;

  // Forward CLs: fwd[node * kStates + s] = P(D_below(node) | node = s).
  std::vector<double>  fwd((maxNode + 1) * kStates, 0.0);
  std::vector<uint8_t> initFwd(maxNode + 1, 0u);

  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = fwd.data() + tip * kStates;
    int state = tipStates[tip - 1];
    if (state < 0) {
      for (int s = 0; s < kStates; ++s) cl[s] = 1.0;
    } else {
      cl[state] = 1.0;
    }
    initFwd[tip] = 1u;
  }

  // Forward pass: postorder traversal = reverse storage order.
  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e];
    int ch  = child[e];
    double t        = edgeLen[e] * rateMultiplier;
    double exp_term = MKP_EXP(-kStates * t / km1);
    double p_same   = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff   = inv_k - inv_k * exp_term;
    double diff_coeff = p_same - p_diff;

    const double* clCh  = fwd.data() + ch  * kStates;
    double*       clPar = fwd.data() + par * kStates;

    double sum_cl = 0.0;
    for (int j = 0; j < kStates; ++j) sum_cl += clCh[j];

    if (!initFwd[par]) {
      for (int i = 0; i < kStates; ++i)
        clPar[i] = p_diff * sum_cl + diff_coeff * clCh[i];
      initFwd[par] = 1u;
    } else {
      for (int i = 0; i < kStates; ++i)
        clPar[i] *= p_diff * sum_cl + diff_coeff * clCh[i];
    }
  }

  // Backward CLs: bwd[node * kStates + s] = P(D_above(node), node = s).
  // Stored absolute (with prior π folded into root) so that
  // marg[n, s] ∝ fwd[n, s] * bwd[n, s] holds at every node.
  std::vector<double> bwd((maxNode + 1) * kStates, 0.0);
  double* bwdRoot = bwd.data() + root * kStates;
  for (int s = 0; s < kStates; ++s) bwdRoot[s] = inv_k;

  // Sibling lookup: each parent has exactly two children (binary rooted tree).
  std::vector<int> firstChildEdge(maxNode + 1, -1);
  std::vector<int> secondChildEdge(maxNode + 1, -1);
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    if (firstChildEdge[par] < 0) firstChildEdge[par]  = e;
    else                          secondChildEdge[par] = e;
  }

  std::vector<double> m(kStates);

  // Backward pass: preorder traversal = forward storage order.
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    int ch  = child[e];
    int sibEdge = (firstChildEdge[par] == e) ? secondChildEdge[par]
                                             : firstChildEdge[par];
    int sib = child[sibEdge];

    // sibling message: msg[t] = Σ_u P_sib_edge(t → u) * fwd[sib, u]
    double t_sib    = edgeLen[sibEdge] * rateMultiplier;
    double exp_s    = MKP_EXP(-kStates * t_sib / km1);
    double ps_sib   = inv_k + (1.0 - inv_k) * exp_s;
    double pd_sib   = inv_k - inv_k * exp_s;
    double diff_sib = ps_sib - pd_sib;

    const double* clSib  = fwd.data() + sib * kStates;
    const double* bwdPar = bwd.data() + par * kStates;

    double sum_sib = 0.0;
    for (int u = 0; u < kStates; ++u) sum_sib += clSib[u];

    // m[t] = bwd[par, t] * sibling_msg[t]
    double sum_m = 0.0;
    for (int t = 0; t < kStates; ++t) {
      double msg_t = pd_sib * sum_sib + diff_sib * clSib[t];
      m[t] = bwdPar[t] * msg_t;
      sum_m += m[t];
    }

    // bwd[ch, s] = Σ_t m[t] * P_ch_edge(t → s)
    //            = p_diff_ch * Σ_t m[t] + (p_same_ch - p_diff_ch) * m[s]
    double t_ch    = edgeLen[e] * rateMultiplier;
    double exp_c   = MKP_EXP(-kStates * t_ch / km1);
    double ps_ch   = inv_k + (1.0 - inv_k) * exp_c;
    double pd_ch   = inv_k - inv_k * exp_c;
    double diff_ch = ps_ch - pd_ch;

    double* bwdCh = bwd.data() + ch * kStates;
    for (int s = 0; s < kStates; ++s)
      bwdCh[s] = pd_ch * sum_m + diff_ch * m[s];
  }

  // Marginals: normalise fwd * bwd at each node.
  marg.assign((maxNode + 1) * kStates, 0.0);
  for (int n = 1; n <= maxNode; ++n) {
    const double* fn = fwd.data() + n * kStates;
    const double* bn = bwd.data() + n * kStates;
    double*       mn = marg.data() + n * kStates;
    double z = 0.0;
    for (int s = 0; s < kStates; ++s) {
      mn[s] = fn[s] * bn[s];
      z += mn[s];
    }
    if (z > 0.0) {
      double inv = 1.0 / z;
      for (int s = 0; s < kStates; ++s) mn[s] *= inv;
    }
  }
}


// Per-node posterior marginals P(node = e | tip ecology data) under JC-K.
// Returns an nNode-by-kEcology matrix where row n-1 corresponds to node n
// (ape 1-indexed: tips 1..nTip, internal nodes nTip+1..2*nTip-1).
//
// [[Rcpp::export(.EcologyNodeMarginals)]]
NumericMatrix EcologyNodeMarginals(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    IntegerVector tipStates,
    int kStates
) {
  int nEdge = parent.size();
  int nTip  = tipStates.size();
  int maxNode = 2 * nTip - 1;

  if (kStates < 2)
    stop("kStates must be >= 2");
  if (child.size() != nEdge || edgeLen.size() != nEdge)
    stop("parent / child / edgeLen must have matching length");

  std::vector<double> marg;
  compute_ecology_node_marginals(
    INTEGER(parent), INTEGER(child), nEdge,
    REAL(edgeLen),
    INTEGER(tipStates), nTip,
    kStates, 1.0,
    marg
  );

  NumericMatrix out(maxNode, kStates);
  for (int n = 1; n <= maxNode; ++n) {
    for (int s = 0; s < kStates; ++s)
      out(n - 1, s) = marg[n * kStates + s];
  }
  return out;
}


// Per-edge ecology weights: mean of parent- and child-node marginals.
// nodeMarginals is the matrix returned by EcologyNodeMarginals (rows indexed
// 0-based as node n -> row n - 1, ape 1-indexed nodes).
//
// [[Rcpp::export(.EcologyEdgeWeights)]]
NumericMatrix EcologyEdgeWeights(
    NumericMatrix nodeMarginals,
    IntegerVector parent, IntegerVector child
) {
  int nEdge   = parent.size();
  int kStates = nodeMarginals.ncol();

  if (child.size() != nEdge)
    stop("parent and child must have matching length");

  NumericMatrix out(nEdge, kStates);
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    int ch  = child[e];
    for (int s = 0; s < kStates; ++s) {
      out(e, s) = 0.5 * (nodeMarginals(par - 1, s) + nodeMarginals(ch - 1, s));
    }
  }
  return out;
}


// v2: compute gamma_e (prior-expected per-cell rate) for non-reference ecology e.
//   gamma_e = pi0 + (1-pi0) * (theta_e * phi + (1-theta_e) / phi).
// For the reference ecology, gamma_e = 1 by construction (rate factor = 1
// always). `phi_e` is phi[0] under global mode and phi[ecoState] otherwise.
static inline double gamma_e_compute(double pi0, double theta_e, double phi_e) {
  return pi0 + (1.0 - pi0) * (theta_e * phi_e + (1.0 - theta_e) / phi_e);
}

// v2: resolve the rate factor for a transformational character.
// Reference ecology: factor = 1 regardless of z. Non-reference:
//   z == 0 (none):        factor = 1 / gamma_e
//   z == 1 (encouraged):  factor = phi / gamma_e
//   z == 2 (discouraged): factor = (1/phi) / gamma_e
// phi_e is the per-ecology phi (or the global phi under mode 0).
// gamma_e is precomputed.
static inline double trans_rate_factor(int z, int ecoState, int refEcology,
                                       const double* phi, int mode,
                                       double gamma_e) {
  if (ecoState == refEcology) return 1.0;
  double p = (mode == 0) ? phi[0] : phi[ecoState];
  double mu = (z == 0) ? 1.0 : (z == 1) ? p : 1.0 / p;
  return mu / gamma_e;
}


// Transformational ecology-aware pruning: JC-K with per-(edge, character)
// mixture transition matrix. The mixture preserves JC symmetry — the off-
// diagonal entry is state-independent — so the per-character propagation
// remains O(k) using mixed (p_same, p_diff) values.
//
// nCat = rate_multipliers.size() (ACRV categories).
// wEdge: nEdge x kEco, edge-ecology weights.
// zMat:  nChar x kEco, integer entries in {0, 1, 2}.
// phi:   length 1 (mode == 0) or kEco (mode == 1).
//
// Buffers buf / initFlg / stride follow the convention used by
// pruning_jc_acrv_flat: buf[node * stride + c * kStates + s] (1-indexed).
//
// Returns total log-likelihood across characters, or R_NegInf if any per-
// character average likelihood is non-positive.
// T-017: raw-pointer variant. The original Rcpp-argument signature is kept as
// a thin wrapper below for callers that already have Rcpp objects in hand.
// This variant takes only POD inputs and is safe to invoke from inside an
// OpenMP parallel region — no R-side state is touched.
//
// Pointer layouts (column-major):
//   parPtr[e], chPtr[e]    : 1-indexed node IDs, length nEdge
//   elPtr[e]               : edge lengths, length nEdge
//   tsPtr[(tip-1) + c*nTip]: tip state, c=0..nChar-1, tip=1..nTip
//   rfPtr[s]               : root frequencies, length kStates
//   rmPtr[cat]             : ACRV rate multipliers, length nCat
//   wPtr[e + s*nEdge]      : per-edge ecology weights, kEco columns
//   zPtr[c + j*nChar]      : per-character z-vector, kEco-1 columns
//   phiPtr[s_or_0]         : phi (length 1 if mode=0 else kEco)
static double pruning_jc_acrv_flat_ecology_raw(
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,
    const int* tsPtr, int nTip, int nChar,
    int kStates, const double* rfPtr,
    const double* rmPtr, int nCat,
    const double* wPtr, int kEco,
    const int* zPtr,
    const double* phiPtr, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {

  int maxNode = 2 * nTip - 1;
  int root    = nTip + 1;
  int clCols  = nChar * kStates;
  double inv_k = 1.0 / kStates;
  double km1   = static_cast<double>(kStates) - 1.0;

  std::vector<double> siteLik(nChar, 0.0);

  // Initialise tips (constant across rate categories).
  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state = tsPtr[(tip - 1) + c * nTip];
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1u;
  }

  // Per-(edge, ecology) base (p_same, p_diff) for the three rate factors
  // {none, encouraged, discouraged}. Stored as small stack buffers re-used
  // each (cat, edge) iteration.
  // ps[3 * kEco + s], indexed by (factor_idx, eco_state):
  //   factor_idx 0 = none (rate factor 1.0, ecology state irrelevant)
  //   factor_idx 1 = encouraged
  //   factor_idx 2 = discouraged
  std::vector<double> psFactor(3 * kEco), pdFactor(3 * kEco);

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rmPtr[cat];

    // Reset internal-node init flags (tips stay initialised).
    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parPtr[e];
      int ch  = chPtr[e];
      double tBase = elPtr[e] * rate;

      // v2: per-(z, s) rate factor incorporates phi and gamma_e normalisation.
      //   ref ecology: factor = 1 for all z (rate unchanged).
      //   non-ref:     z=0 -> 1/gamma_e; z=1 -> phi/gamma_e; z=2 -> (1/phi)/gamma_e.
      {
        for (int s = 0; s < kEco; ++s) {
          double gE = gammaE[s];
          for (int z = 0; z < 3; ++z) {
            double factor;
            if (s == refEcology) {
              factor = 1.0;
            } else {
              double p = (mode == 0) ? phiPtr[0] : phiPtr[s];
              double mu = (z == 0) ? 1.0 : (z == 1) ? p : 1.0 / p;
              factor = mu / gE;
            }
            double tEff = tBase * factor;
            double exV  = MKP_EXP(-kStates * tEff / km1);
            psFactor[z * kEco + s] = inv_k + (1.0 - inv_k) * exV;
            pdFactor[z * kEco + s] = inv_k - inv_k * exV;
          }
        }
      }

      double* clPar = buf + par * stride;
      double* clCh  = buf + ch  * stride;

      // T-009 (2026-05-21): tip-edge fast path.  When ch is a tip, the child
      // CL is one-hot (known state s: clCh[s]=1, rest 0) or all-ones (missing).
      // sum_cl = 1 (known) or kStates (missing) by construction; the per-state
      // update collapses:
      //   init + known s:   clPar[i] = pdMix (i≠s), psMix (i==s)
      //   init + missing:   clPar[i] = 1.0  (JC identity: Σ w_s·(K·pdF+(psF−pdF))=1)
      //   multiply + known: clPar[i] *= pdMix (i≠s), psMix (i==s)
      //   multiply + missing: no-op (× 1.0)
      // psMix/pdMix are computed in the same accumulation order as the internal
      // path to preserve FP results for those quantities.
      const bool tipChild = (ch <= nTip);

      if (tipChild) {
        if (!initFlg[par]) {
          for (int c = 0; c < nChar; ++c) {
            double psMix = 0.0, pdMix = 0.0;
            for (int s = 0; s < kEco; ++s) {
              int z;
              if (s == refEcology) {
                z = 0;
              } else {
                int zCol = (s < refEcology) ? s : (s - 1);
                z = zPtr[c + zCol * nChar];
              }
              double w = wPtr[e + s * nEdge];
              psMix += w * psFactor[z * kEco + s];
              pdMix += w * pdFactor[z * kEco + s];
            }
            int state  = tsPtr[(ch - 1) + c * nTip];
            int offset = c * kStates;
            if (state < 0) {
              // Missing: full marginalisation → clPar[i] = 1.
              for (int i = 0; i < kStates; ++i) clPar[offset + i] = 1.0;
            } else {
              for (int i = 0; i < kStates; ++i) clPar[offset + i] = pdMix;
              clPar[offset + state] = psMix;
            }
          }
          initFlg[par] = 1u;
        } else {
          for (int c = 0; c < nChar; ++c) {
            double psMix = 0.0, pdMix = 0.0;
            for (int s = 0; s < kEco; ++s) {
              int z;
              if (s == refEcology) {
                z = 0;
              } else {
                int zCol = (s < refEcology) ? s : (s - 1);
                z = zPtr[c + zCol * nChar];
              }
              double w = wPtr[e + s * nEdge];
              psMix += w * psFactor[z * kEco + s];
              pdMix += w * pdFactor[z * kEco + s];
            }
            int state  = tsPtr[(ch - 1) + c * nTip];
            if (state < 0) {
              // Missing: multiply by 1 — no-op.
            } else {
              int offset = c * kStates;
              for (int i = 0; i < kStates; ++i)
                clPar[offset + i] *= (i == state) ? psMix : pdMix;
            }
          }
        }
        continue;
      }

      // Internal-child path (unchanged).
      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          double psMix = 0.0, pdMix = 0.0;
          for (int s = 0; s < kEco; ++s) {
            // v2: zMat has kEco-1 cols (no column for ref ecology).
            // For ref ecology, z is treated as 0 (factor is 1 anyway).
            int z;
            if (s == refEcology) {
              z = 0;
            } else {
              // Column index among non-ref ecologies, ascending order.
              int zCol = (s < refEcology) ? s : (s - 1);
              z = zPtr[c + zCol * nChar];
            }
            double w = wPtr[e + s * nEdge];
            psMix += w * psFactor[z * kEco + s];
            pdMix += w * pdFactor[z * kEco + s];
          }
          double diff_coeff = psMix - pdMix;
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] = pdMix * sum_cl + diff_coeff * clCh[offset + i];
        }
        initFlg[par] = 1u;
      } else {
        for (int c = 0; c < nChar; ++c) {
          double psMix = 0.0, pdMix = 0.0;
          for (int s = 0; s < kEco; ++s) {
            // v2: zMat has kEco-1 cols (no column for ref ecology).
            // For ref ecology, z is treated as 0 (factor is 1 anyway).
            int z;
            if (s == refEcology) {
              z = 0;
            } else {
              // Column index among non-ref ecologies, ascending order.
              int zCol = (s < refEcology) ? s : (s - 1);
              z = zPtr[c + zCol * nChar];
            }
            double w = wPtr[e + s * nEdge];
            psMix += w * psFactor[z * kEco + s];
            pdMix += w * pdFactor[z * kEco + s];
          }
          double diff_coeff = psMix - pdMix;
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] *= pdMix * sum_cl + diff_coeff * clCh[offset + i];
        }
      }
    }

    double* clRoot = buf + root * stride;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double sl = 0.0;
      for (int s = 0; s < kStates; ++s)
        sl += rfPtr[s] * clRoot[offset + s];
      siteLik[c] += sl;
    }
  }

  double logLik   = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = siteLik[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// Thin Rcpp-arg wrapper around the raw-pointer pruner.  Kept for callers
// that already hold Rcpp objects and run serially.  Const-ref args avoid
// the Rcpp Vector copy-constructor (which touches the precious-object list
// and is NOT thread-safe).
static inline double pruning_jc_acrv_flat_ecology(
    const Rcpp::IntegerVector& parent, const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edge_length,
    const Rcpp::IntegerMatrix& tip_states,
    int kStates, const Rcpp::NumericVector& root_freqs,
    const Rcpp::NumericVector& rate_multipliers,
    const Rcpp::NumericMatrix& wEdge,
    const Rcpp::IntegerMatrix& zMat,
    const Rcpp::NumericVector& phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {
  return pruning_jc_acrv_flat_ecology_raw(
    INTEGER(parent), INTEGER(child), parent.size(),
    REAL(edge_length),
    INTEGER(tip_states), tip_states.nrow(), tip_states.ncol(),
    kStates, REAL(root_freqs),
    REAL(rate_multipliers), rate_multipliers.size(),
    REAL(wEdge), wEdge.ncol(),
    INTEGER(zMat),
    REAL(phi), mode, refEcology, gammaE,
    buf, initFlg, stride);
}


// Resolve the (rate01, rate10) pair for a neomorphic character with
// influence category z under magnitude mode `mode`. Unlike the
// transformational case, the modification is per-direction:
//   z == 0 (none):        (rate01,        rate10)
//   z == 1 (encouraged):  (rate01 * phi,  rate10 / phi)
//   z == 2 (discouraged): (rate01 / phi,  rate10 * phi)
// phi has length 1 in global mode and kEco in per_ecology mode.
// v2: reference ecology has rate factor = 1 (base rates unchanged).
// Non-reference: scale both directions by 1/gamma_e (z = 0), phi/gamma_e on
// rate01 and (1/phi)/gamma_e on rate10 (z = 1), or vice versa (z = 2).
static inline void mkn_rates_for_state(
    int z, int ecoState, int refEcology,
    double rate01_base, double rate10_base,
    const double* phi, int mode,
    double gamma_e,
    double& r01, double& r10) {
  if (ecoState == refEcology) {
    r01 = rate01_base;
    r10 = rate10_base;
    return;
  }
  double p = (mode == 0) ? phi[0] : phi[ecoState];
  double mu01, mu10;
  if (z == 0) {
    mu01 = 1.0; mu10 = 1.0;
  } else if (z == 1) {
    mu01 = p;   mu10 = 1.0 / p;
  } else {
    mu01 = 1.0 / p; mu10 = p;
  }
  r01 = rate01_base * mu01 / gamma_e;
  r10 = rate10_base * mu10 / gamma_e;
}


// Neomorphic ecology-aware pruning: 2-state asymmetric MkN with per-(edge,
// character) mixture transition matrix. Each ecology state contributes a
// full 2x2 P^{(s)}(t) and the mixture is a weighted sum of these. The
// mixture is not in general of MkN form (different stationary frequencies
// across ecology states), so the propagation is a full 2x2 matvec.
//
// kStates is hard-coded to 2.
// root_freqs is the prior at the root — typically the base MkN stationary
// distribution under the unmodified rates (rate10/lambda, rate01/lambda).
//
// Returns total log-likelihood across characters, or R_NegInf if any per-
// character average likelihood is non-positive.
// T-017: raw-pointer variant.  See header comment for pruning_jc_acrv_flat_ecology_raw.
static double pruning_mkn_acrv_flat_ecology_raw(
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,
    const int* tsPtr, int nTip, int nChar,
    double rate_loss, const double* rfPtr,
    const double* rmPtr, int nCat,
    const double* wPtr, int kEco,
    const int* zPtr,
    const double* phiPtr, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {

  const int kStates = 2;

  int maxNode = 2 * nTip - 1;
  int root    = nTip + 1;
  int clCols  = nChar * kStates;

  // Base rates (unmodified MkN with stationary pi_0 = rate_loss/(1+rate_loss)).
  double sum_rl = 1.0 + rate_loss;
  double r01_base = 2.0 / sum_rl;
  double r10_base = 2.0 * rate_loss / sum_rl;

  std::vector<double> siteLik(nChar, 0.0);

  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state = tsPtr[(tip - 1) + c * nTip];
      int offset = c * kStates;
      if (state < 0) { cl[offset] = 1.0; cl[offset + 1] = 1.0; }
      else            cl[offset + state] = 1.0;
    }
    initFlg[tip] = 1u;
  }

  // Per-(edge, factor variant, ecology state) precomputed P matrices.
  // P[v * kEco * 4 + s * 4 + ij], variant v in {0=none, 1=enc, 2=dis}.
  std::vector<double> Pfactor(3 * kEco * 4);

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rmPtr[cat];

    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parPtr[e];
      int ch  = chPtr[e];
      double tBase = elPtr[e] * rate;

      // v2: precompute P matrices for each (z, s). Reference ecology uses base
      // rates (no scaling); non-reference uses (r01 * mu01, r10 * mu10) / gamma_e
      // per mkn_rates_for_state.
      for (int s = 0; s < kEco; ++s) {
        double gE = gammaE[s];
        for (int z = 0; z < 3; ++z) {
          double r01, r10;
          mkn_rates_for_state(z, s, refEcology, r01_base, r10_base,
                              phiPtr, mode, gE, r01, r10);
          double lam = r01 + r10;
          double pi0_ = r10 / lam;
          double pi1_ = r01 / lam;
          double ex  = MKP_EXP(-lam * tBase);
          double P00 = pi0_ + pi1_ * ex;
          double P01 = pi1_ - pi1_ * ex;
          double P10 = pi0_ - pi0_ * ex;
          double P11 = pi1_ + pi0_ * ex;
          double* P = Pfactor.data() + (z * kEco + s) * 4;
          P[0] = P00; P[1] = P01; P[2] = P10; P[3] = P11;
        }
      }

      double* clPar = buf + par * stride;
      double* clCh  = buf + ch  * stride;

      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          double P00m = 0.0, P01m = 0.0, P10m = 0.0, P11m = 0.0;
          for (int s = 0; s < kEco; ++s) {
            int z;
            if (s == refEcology) {
              z = 0;
            } else {
              int zCol = (s < refEcology) ? s : (s - 1);
              z = zPtr[c + zCol * nChar];
            }
            double w = wPtr[e + s * nEdge];
            const double* P = Pfactor.data() + (z * kEco + s) * 4;
            P00m += w * P[0]; P01m += w * P[1];
            P10m += w * P[2]; P11m += w * P[3];
          }
          int offset = c * kStates;
          double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
          clPar[offset]     = P00m * cl0 + P01m * cl1;
          clPar[offset + 1] = P10m * cl0 + P11m * cl1;
        }
        initFlg[par] = 1u;
      } else {
        for (int c = 0; c < nChar; ++c) {
          double P00m = 0.0, P01m = 0.0, P10m = 0.0, P11m = 0.0;
          for (int s = 0; s < kEco; ++s) {
            int z;
            if (s == refEcology) {
              z = 0;
            } else {
              int zCol = (s < refEcology) ? s : (s - 1);
              z = zPtr[c + zCol * nChar];
            }
            double w = wPtr[e + s * nEdge];
            const double* P = Pfactor.data() + (z * kEco + s) * 4;
            P00m += w * P[0]; P01m += w * P[1];
            P10m += w * P[2]; P11m += w * P[3];
          }
          int offset = c * kStates;
          double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
          clPar[offset]     *= P00m * cl0 + P01m * cl1;
          clPar[offset + 1] *= P10m * cl0 + P11m * cl1;
        }
      }
    }

    double* clRoot = buf + root * stride;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double sl = rfPtr[0] * clRoot[offset] + rfPtr[1] * clRoot[offset + 1];
      siteLik[c] += sl;
    }
  }

  double logLik   = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = siteLik[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// Thin Rcpp-arg wrapper around the raw-pointer MkN pruner.  See note on
// pruning_jc_acrv_flat_ecology above re const-ref / thread-safety.
static inline double pruning_mkn_acrv_flat_ecology(
    const Rcpp::IntegerVector& parent, const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edge_length,
    const Rcpp::IntegerMatrix& tip_states,
    double rate_loss, const Rcpp::NumericVector& root_freqs,
    const Rcpp::NumericVector& rate_multipliers,
    const Rcpp::NumericMatrix& wEdge,
    const Rcpp::IntegerMatrix& zMat,
    const Rcpp::NumericVector& phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {
  return pruning_mkn_acrv_flat_ecology_raw(
    INTEGER(parent), INTEGER(child), parent.size(),
    REAL(edge_length),
    INTEGER(tip_states), tip_states.nrow(), tip_states.ncol(),
    rate_loss, REAL(root_freqs),
    REAL(rate_multipliers), rate_multipliers.size(),
    REAL(wEdge), wEdge.ncol(),
    INTEGER(zMat),
    REAL(phi), mode, refEcology, gammaE,
    buf, initFlg, stride);
}


// R-callable wrapper around pruning_mkn_acrv_flat_ecology for testing.
//
// [[Rcpp::export(.PruningMknEcology)]]
double PruningMknEcology(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, IntegerMatrix tipStates,
    double rateLoss, NumericVector rootFreqs,
    NumericVector rateMultipliers,
    NumericMatrix wEdge,
    IntegerMatrix zMat,
    NumericVector phi, int mode,
    int refEcology = -1,
    Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue,
    double pi0 = 0.0
) {
  int nTip    = tipStates.nrow();
  int nChar   = tipStates.ncol();
  int maxNode = 2 * nTip - 1;
  int stride  = nChar * 2;
  int kEco    = wEdge.ncol();

  if (rootFreqs.size() != 2)
    stop("rootFreqs length must equal 2");
  if (wEdge.nrow() != parent.size())
    stop("wEdge nrow must equal nEdge");
  if (zMat.nrow() != nChar)
    stop("zMat nrow must equal nChar");
  // v2: zMat has kEco-1 columns (one per non-ref ecology).
  if (zMat.ncol() != kEco - 1)
    stop("zMat ncol must equal kEco - 1 (v2)");
  if (mode != 0 && mode != 1)
    stop("mode must be 0 (global) or 1 (per_ecology)");
  if (mode == 0 && phi.size() != 1 && phi.size() != kEco)
    stop("global mode requires phi length 1 or kEco");
  if (mode == 1 && phi.size() != kEco)
    stop("per_ecology mode requires phi length kEco");
  if (rateLoss <= 0)
    stop("rateLoss must be positive");
  if (refEcology < 0 || refEcology >= kEco)
    stop("refEcology must be in [0, kEco)");

  // Build gammaE from theta + pi0. gammaE[refEcology] = 1.
  std::vector<double> gammaE(kEco, 1.0);
  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(kEco - 1, 0.5);
  if (thetaVec.size() != kEco - 1)
    stop("theta length must equal kEco - 1");
  for (int s = 0; s < kEco; ++s) {
    if (s == refEcology) { gammaE[s] = 1.0; continue; }
    int j = (s < refEcology) ? s : (s - 1);
    double phi_s = (mode == 0) ? phi[0] : phi[s];
    gammaE[s] = gamma_e_compute(pi0, thetaVec[j], phi_s);
  }

  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  return pruning_mkn_acrv_flat_ecology(
    parent, child, edgeLen, tipStates,
    rateLoss, rootFreqs, rateMultipliers,
    wEdge, zMat, phi, mode,
    refEcology, gammaE,
    buf.data(), initFlg.data(), stride
  );
}


// R-callable wrapper around pruning_jc_acrv_flat_ecology for testing.
// Allocates buffers locally; not for hot-path use.
//
// [[Rcpp::export(.PruningJcEcology)]]
double PruningJcEcology(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, IntegerMatrix tipStates,
    int kStates, NumericVector rootFreqs,
    NumericVector rateMultipliers,
    NumericMatrix wEdge,
    IntegerMatrix zMat,
    NumericVector phi, int mode,
    int refEcology = -1,
    Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue,
    double pi0 = 0.0
) {
  int nTip    = tipStates.nrow();
  int nChar   = tipStates.ncol();
  int maxNode = 2 * nTip - 1;
  int stride  = nChar * kStates;
  int kEco    = wEdge.ncol();

  if (kStates < 2)
    stop("kStates must be >= 2");
  if (rootFreqs.size() != kStates)
    stop("rootFreqs length must equal kStates");
  if (wEdge.nrow() != parent.size())
    stop("wEdge nrow must equal nEdge");
  if (zMat.nrow() != nChar)
    stop("zMat nrow must equal nChar");
  if (zMat.ncol() != kEco - 1)
    stop("zMat ncol must equal kEco - 1 (v2)");
  if (mode != 0 && mode != 1)
    stop("mode must be 0 (global) or 1 (per_ecology)");
  if (mode == 0 && phi.size() != 1 && phi.size() != kEco)
    stop("global mode requires phi length 1 or kEco");
  if (mode == 1 && phi.size() != kEco)
    stop("per_ecology mode requires phi length kEco");
  if (refEcology < 0 || refEcology >= kEco)
    stop("refEcology must be in [0, kEco)");

  std::vector<double> gammaE(kEco, 1.0);
  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(kEco - 1, 0.5);
  if (thetaVec.size() != kEco - 1)
    stop("theta length must equal kEco - 1");
  for (int s = 0; s < kEco; ++s) {
    if (s == refEcology) { gammaE[s] = 1.0; continue; }
    int j = (s < refEcology) ? s : (s - 1);
    double phi_s = (mode == 0) ? phi[0] : phi[s];
    gammaE[s] = gamma_e_compute(pi0, thetaVec[j], phi_s);
  }

  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  return pruning_jc_acrv_flat_ecology(
    parent, child, edgeLen, tipStates,
    kStates, rootFreqs, rateMultipliers,
    wEdge, zMat, phi, mode,
    refEcology, gammaE,
    buf.data(), initFlg.data(), stride
  );
}


// ---------------------------------------------------------------------------
// Ascertainment under the ecology mixture (variable coding)
// ---------------------------------------------------------------------------
//
// Per-character constant-site probabilities. Each character has its own
// z[c, :] vector and therefore its own mixture transition matrix, so the
// constant-site probability is computed per character via a single
// pseudo-character pruning pass.

// T-017: EcoWorkItem is one unit of pruning work on a single
// (partition, kPrime-subgroup) pair.  For neomorphic / known partitions
// there is exactly one item per partition (subgroupKp = -1).  For
// transformational partitions there is one item per kPrime subgroup
// (subgroupKp = the kPrime value, which is also the state-space size).
struct EcoWorkItem {
  int partIdx;       // index into data.parts[]
  int subgroupKp;    // -1 = neo/known; else kPrime value
};

// T-017 Phase 1b: per-thread scratch workspace.  Eliminates the per-call
// malloc/free that Phase 1 left in `compute_eco_work_item_ll`.
// T-016: extended with tipStatesZero/tipStatesOne so const_site_prob_*_raw
// can also use caller-supplied scratch instead of allocating per call.
struct EcoThreadScratch {
  std::vector<double>  buf;          // (maxNode + 1) * maxStride doubles
  std::vector<uint8_t> initFlg;      // (maxNode + 1) bytes
  std::vector<int>     zPartLocal;   // maxNCharPart * zCols
  std::vector<int>     subStates;    // nTip * maxNSub
  std::vector<int>     subZ;         // maxNSub * zCols
  std::vector<int>     cols;         // up to maxNCharPart
  std::vector<int>     zVec;         // zCols
  std::vector<double>  neoEl;        // nEdge doubles
  std::vector<double>  rootFreqs;    // up to maxKStates doubles
  std::vector<int>     tipStatesZero; // nTip zeros  — for const_site_prob single-char passes
  std::vector<int>     tipStatesOne;  // nTip ones   — for const_site_prob_mkn second pass
  int maxNode  = 0;
  int maxStride = 0;
  int maxNTip  = 0;

  void allocate(int maxNode_, int maxStride_, int maxNCharPart,
                int maxNSub, int zCols, int nTip, int nEdge,
                int maxKStates) {
    maxNode   = maxNode_;
    maxStride = maxStride_;
    maxNTip   = nTip;
    buf.assign(static_cast<size_t>(maxNode_ + 1) * maxStride_, 0.0);
    initFlg.assign(maxNode_ + 1, 0u);
    zPartLocal.assign(static_cast<size_t>(maxNCharPart) * std::max(1, zCols), 0);
    subStates.assign(static_cast<size_t>(nTip) * std::max(1, maxNSub), 0);
    subZ.assign(static_cast<size_t>(maxNSub) * std::max(1, zCols), 0);
    cols.assign(maxNCharPart, 0);
    zVec.assign(std::max(1, zCols), 0);
    neoEl.assign(nEdge, 0.0);
    rootFreqs.assign(std::max(2, maxKStates), 0.0);
    tipStatesZero.assign(nTip, 0);
    tipStatesOne.assign(nTip, 1);
  }
};

// T-017: raw-pointer variants of const_site_prob_*_eco_single.  All Rcpp
// temporaries replaced with std::vector so these are safe to call from
// inside an OpenMP parallel region.
//
// Raw-pointer layout matches the inner pruner:
//   parPtr[e], chPtr[e]       length nEdge
//   elPtr[e]                  length nEdge
//   wPtr[e + s*nEdge]         column-major (nEdge x kEco)
//   zVecPtr[j]                length (kEco - 1)
//   phiPtr[s_or_0]            length 1 (mode=0) or kEco (mode=1)
//   rmPtr[cat]                length nCat
static double const_site_prob_jc_eco_single_raw(
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,
    int nTip, int kStates,
    const double* rmPtr, int nCat,
    const double* wPtr, int kEco,
    const int* zVecPtr,           // length (kEco - 1); caller-owned, passed to pruner
    const double* phiPtr, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    EcoThreadScratch& scratch) {  // T-016: caller-owned scratch; no per-call alloc

  double* rootFreqs = scratch.rootFreqs.data();
  for (int k = 0; k < kStates; ++k) rootFreqs[k] = 1.0 / kStates;

  int maxNode = 2 * nTip - 1;
  int stride  = kStates;
  size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
  std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
  std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);

  // JC symmetry: P(constant in any state) = kStates * P(all-tips-0).
  double ll = pruning_jc_acrv_flat_ecology_raw(
    parPtr, chPtr, nEdge,
    elPtr,
    scratch.tipStatesZero.data(), nTip, /*nChar=*/1,
    kStates, rootFreqs,
    rmPtr, nCat,
    wPtr, kEco,
    zVecPtr,
    phiPtr, mode, refEcology, gammaE,
    scratch.buf.data(), scratch.initFlg.data(), stride);
  return kStates * std::exp(ll);
}


static double const_site_prob_mkn_eco_single_raw(
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,          // caller must pass already-scaled neoEl
    int nTip,
    double rateLoss,
    const double* rmPtr, int nCat,
    const double* wPtr, int kEco,
    const int* zVecPtr,           // length (kEco - 1); caller-owned, passed to pruner
    const double* phiPtr, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    EcoThreadScratch& scratch) {  // T-016: caller-owned scratch; no per-call alloc

  double* rootFreqs = scratch.rootFreqs.data();
  rootFreqs[0] = rateLoss / (1.0 + rateLoss);
  rootFreqs[1] = 1.0 / (1.0 + rateLoss);

  int maxNode = 2 * nTip - 1;
  int stride  = 2;
  size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;

  // Pseudo-char "all 0"
  std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
  std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);
  double ll0 = pruning_mkn_acrv_flat_ecology_raw(
    parPtr, chPtr, nEdge,
    elPtr,
    scratch.tipStatesZero.data(), nTip, /*nChar=*/1,
    rateLoss, rootFreqs,
    rmPtr, nCat,
    wPtr, kEco,
    zVecPtr,
    phiPtr, mode, refEcology, gammaE,
    scratch.buf.data(), scratch.initFlg.data(), stride);

  // Pseudo-char "all 1"
  std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
  std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);
  double ll1 = pruning_mkn_acrv_flat_ecology_raw(
    parPtr, chPtr, nEdge,
    elPtr,
    scratch.tipStatesOne.data(), nTip, /*nChar=*/1,
    rateLoss, rootFreqs,
    rmPtr, nCat,
    wPtr, kEco,
    zVecPtr,
    phiPtr, mode, refEcology, gammaE,
    scratch.buf.data(), scratch.initFlg.data(), stride);

  return std::exp(ll0) + std::exp(ll1);
}


// Thin Rcpp-arg wrappers around the raw variants.  Kept so serial callers
// (per_char_log_lik_ecology in the Gibbs z sweep) need no signature change.
// T-016: each wrapper carries a function-static EcoThreadScratch; the zVec
// is copied once from the Rcpp IntegerVector into scratch.zVec, then the raw
// variant receives scratch.zVec.data() — no allocation inside _raw.
static inline double const_site_prob_jc_eco_single(
    const Rcpp::IntegerVector& parent, const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen, int nTip, int kStates,
    const Rcpp::NumericVector& rates,
    const Rcpp::NumericMatrix& wEdge, const Rcpp::IntegerVector& zVec,
    const Rcpp::NumericVector& phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE) {
  static EcoThreadScratch jcScratch;
  int nEdge  = parent.size();
  int kEco   = wEdge.ncol();
  int zCols  = std::max(1, kEco - 1);
  int maxNode = 2 * nTip - 1;
  if (jcScratch.maxNode < maxNode || jcScratch.maxStride < kStates ||
      (int)jcScratch.rootFreqs.size() < kStates ||
      (int)jcScratch.zVec.size() < zCols ||
      jcScratch.maxNTip < nTip) {
    jcScratch.allocate(maxNode, kStates, 1, 1, zCols, nTip, nEdge, kStates);
  }
  int* zPtr = jcScratch.zVec.data();
  for (int j = 0; j < kEco - 1; ++j) zPtr[j] = zVec[j];
  return const_site_prob_jc_eco_single_raw(
    INTEGER(parent), INTEGER(child), nEdge,
    REAL(edgeLen),
    nTip, kStates,
    REAL(rates), rates.size(),
    REAL(wEdge), kEco,
    zPtr,
    REAL(phi), mode, refEcology, gammaE,
    jcScratch);
}

static inline double const_site_prob_mkn_eco_single(
    const Rcpp::IntegerVector& parent, const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen, int nTip,
    double rateLoss, const Rcpp::NumericVector& rates,
    const Rcpp::NumericMatrix& wEdge, const Rcpp::IntegerVector& zVec,
    const Rcpp::NumericVector& phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE) {
  static EcoThreadScratch mknScratch;
  int nEdge  = parent.size();
  int kEco   = wEdge.ncol();
  int zCols  = std::max(1, kEco - 1);
  int maxNode = 2 * nTip - 1;
  if (mknScratch.maxNode < maxNode || mknScratch.maxStride < 2 ||
      (int)mknScratch.zVec.size() < zCols ||
      mknScratch.maxNTip < nTip) {
    mknScratch.allocate(maxNode, 2, 1, 1, zCols, nTip, nEdge, 2);
  }
  int* zPtr = mknScratch.zVec.data();
  for (int j = 0; j < kEco - 1; ++j) zPtr[j] = zVec[j];
  return const_site_prob_mkn_eco_single_raw(
    INTEGER(parent), INTEGER(child), nEdge,
    REAL(edgeLen),
    nTip,
    rateLoss,
    REAL(rates), rates.size(),
    REAL(wEdge), kEco,
    zPtr,
    REAL(phi), mode, refEcology, gammaE,
    mknScratch);
}


// ---------------------------------------------------------------------------
// cpp_partition_log_likelihood_ecology: per-partition variant. (T-010)
// ---------------------------------------------------------------------------
//
// Mirrors the inner-loop body of cpp_log_likelihood_ecology for a single
// partition `partIdx`. Accepts pre-computed wEdge, gammaE, and ACRV rates
// so the orchestrator (or move handler) can hoist them out of the
// per-partition loop. Enables partition-level caching for non-tree moves
// in the ecology-aware path.
//
// Behaviour is bit-identical to one inner-loop iteration of the orchestrator
// (same per-partition summation order: transformational characters are
// sub-grouped by kPrime via std::map so ascending-k order is preserved).

// T-017: thread-safe partition-likelihood core.  All inputs are raw POD;
// no Rcpp constructors are invoked, so this can be called from inside an
// OpenMP parallel region.  All locally-needed scratch (zPart-equivalent,
// subset states, root freqs, neoEl, dirichlet/buffer/initFlg) is held in
// std::vector — glibc and Windows allocator are both thread-safe.
//
// Pointer layouts:
//   parPtr, chPtr      length nEdge
//   elPtr              length nEdge
//   wPtr               nEdge x kEco column-major (e + s*nEdge)
//   tipStatesGlobal    nTip x part.tipStates.ncol() column-major, length
//                      nTip * nCharPart, the global per-partition tip-state
//                      matrix; the function selects sub-columns as needed.
//   kPrimePtr          length nCharsTotal (global indices), only the entries
//                      part.globalCharIdx[c] are read.
//   globalCharIdxPtr   length nCharPart, maps local c -> global char index
//   phiPtr             length 1 (mode 0) or kEco (mode 1)
//   ratesPtr           length nCat
//   zMatrixGlobalPtr   nCharsTotal x (kEco-1) column-major; selects rows by
//                      globalCharIdx[c].  Pointer layout: zPtr[gi + j*nCharsTotal].
//   thinking-cap layout note (zPart-equivalent constructed locally):
//     local zPart is built column-major as nCharPart x (kEco-1):
//       zPart[c + j*nCharPart] = zMatrixGlobalPtr[globalCharIdxPtr[c] + j*nCharsTotal]
//     This matches the pruner's expected layout for INTEGER(zMat) where
//     zPtr[c + j*nChar] is read in the per-char loops.
static double cpp_partition_log_likelihood_ecology_raw(
    const McmcData& data, int partIdx,
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,
    const int* kPrimePtr,
    double rateLoss, double rateNeo,
    const double* phiPtr, int phiLen,
    const int* zMatrixGlobalPtr, int nCharsTotal,
    const double* wPtr,
    const std::vector<double>& gammaE,
    const double* ratesPtr, int nCat) {

  int nTip = data.nTip;
  int maxNode = 2 * nTip - 1;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int zCols = kEco - 1;

  const PartInfo& part = data.parts[partIdx];
  int nCharPart = part.tipStates.ncol();
  const int* partTipsPtr = INTEGER(part.tipStates);   // nTip x nCharPart, col-major
  const int* gciPtr      = INTEGER(part.globalCharIdx);

  // T-016: function-static scratch — eliminates per-call malloc/free for all
  // local vectors.  Serial-only (all callers are in the serial chain loop;
  // no OpenMP region wraps this function — confirmed by grep of mcmc.cpp).
  // Grow-only resize: re-allocate if any dimension exceeds current capacity.
  int maxKStates = 2, maxNSub = 1;
  if (part.type == 2) {
    maxKStates = part.k;
  } else if (part.type == 1) {
    for (int c = 0; c < nCharPart; ++c) {
      int kp = kPrimePtr[gciPtr[c]];
      if (kp > maxKStates) maxKStates = kp;
    }
    maxNSub = nCharPart;  // upper bound
  }
  int maxStride = nCharPart * std::max(2, maxKStates);

  static EcoThreadScratch orcScratch;
  if (orcScratch.maxNode < maxNode || orcScratch.maxStride < maxStride ||
      (int)orcScratch.zPartLocal.size() < nCharPart * std::max(1, zCols) ||
      (int)orcScratch.subStates.size() < nTip * maxNSub ||
      (int)orcScratch.neoEl.size() < nEdge ||
      (int)orcScratch.rootFreqs.size() < maxKStates ||
      orcScratch.maxNTip < nTip) {
    orcScratch.allocate(maxNode, maxStride, std::max(nCharPart, maxNSub),
                        maxNSub, std::max(1, zCols), nTip, nEdge, maxKStates);
  }

  // Build local zPart: nCharPart x zCols column-major, layout zPartLocal[c + j*nCharPart].
  int* zPartLocal = orcScratch.zPartLocal.data();
  for (int c = 0; c < nCharPart; ++c) {
    int gi = gciPtr[c];
    for (int j = 0; j < zCols; ++j)
      zPartLocal[c + j * nCharPart] = zMatrixGlobalPtr[gi + j * nCharsTotal];
  }

  // Suppress unused warnings (phiLen captured for documentation only).
  (void)phiLen;

  double ll = 0.0;

  if (part.type == 0) {
    // Neomorphic.  Build neoEl = edgeLen * rateNeo.
    double* neoEl = orcScratch.neoEl.data();
    for (int i = 0; i < nEdge; ++i) neoEl[i] = elPtr[i] * rateNeo;
    double* rootFreqs = orcScratch.rootFreqs.data();
    rootFreqs[0] = rateLoss / (1.0 + rateLoss);
    rootFreqs[1] = 1.0 / (1.0 + rateLoss);
    int stride = nCharPart * 2;
    size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
    std::fill_n(orcScratch.buf.begin(), bufBytes, 0.0);
    std::fill_n(orcScratch.initFlg.begin(), maxNode + 1, 0u);
    ll = pruning_mkn_acrv_flat_ecology_raw(
      parPtr, chPtr, nEdge,
      neoEl,
      partTipsPtr, nTip, nCharPart,
      rateLoss, rootFreqs,
      ratesPtr, nCat,
      wPtr, kEco,
      zPartLocal,
      phiPtr, mode, refE, gammaE,
      orcScratch.buf.data(), orcScratch.initFlg.data(), stride);
    if (data.codingType == 1) {  // variable
      // zVec for character c is column c of zPartLocal.
      int* zVec = orcScratch.zVec.data();
      for (int c = 0; c < nCharPart; ++c) {
        for (int j = 0; j < zCols; ++j)
          zVec[j] = zPartLocal[c + j * nCharPart];
        double pConst = const_site_prob_mkn_eco_single_raw(
          parPtr, chPtr, nEdge,
          neoEl,
          nTip,
          rateLoss,
          ratesPtr, nCat,
          wPtr, kEco,
          zVec,
          phiPtr, mode, refE, gammaE,
          orcScratch);
        ll -= std::log(1.0 - pConst);
      }
    }
  } else if (part.type == 2) {
    // Known state space.
    int kStates = part.k;
    double* rootFreqs = orcScratch.rootFreqs.data();
    for (int k = 0; k < kStates; ++k) rootFreqs[k] = 1.0 / kStates;
    int stride = nCharPart * kStates;
    size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
    std::fill_n(orcScratch.buf.begin(), bufBytes, 0.0);
    std::fill_n(orcScratch.initFlg.begin(), maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology_raw(
      parPtr, chPtr, nEdge,
      elPtr,
      partTipsPtr, nTip, nCharPart,
      kStates, rootFreqs,
      ratesPtr, nCat,
      wPtr, kEco,
      zPartLocal,
      phiPtr, mode, refE, gammaE,
      orcScratch.buf.data(), orcScratch.initFlg.data(), stride);
    if (data.codingType == 1) {
      int* zVec = orcScratch.zVec.data();
      for (int c = 0; c < nCharPart; ++c) {
        for (int j = 0; j < zCols; ++j)
          zVec[j] = zPartLocal[c + j * nCharPart];
        double pConst = const_site_prob_jc_eco_single_raw(
          parPtr, chPtr, nEdge,
          elPtr,
          nTip, kStates,
          ratesPtr, nCat,
          wPtr, kEco,
          zVec,
          phiPtr, mode, refE, gammaE,
          orcScratch);
        ll -= std::log(1.0 - pConst);
      }
    }
  } else {
    // Transformational: subgroup by kPrime (per-character).
    // NOTE: when one char's kPrime changes, the subgroup composition
    // changes — easiest correct option is to recompute the whole
    // transformational partition. Optimisation deferred (T-010 caveat).
    std::map<int, std::vector<int>> byKp;
    for (int c = 0; c < nCharPart; ++c) {
      int gi = gciPtr[c];
      byKp[kPrimePtr[gi]].push_back(c);
    }
    for (auto& kv : byKp) {
      int kp = kv.first;
      const std::vector<int>& cols = kv.second;
      int nSub = static_cast<int>(cols.size());

      // subStates: nTip x nSub column-major.  subZ: nSub x zCols column-major.
      int* subStates = orcScratch.subStates.data();
      int* subZ      = orcScratch.subZ.data();
      for (int c = 0; c < nSub; ++c) {
        for (int t = 0; t < nTip; ++t)
          subStates[t + c * nTip] = partTipsPtr[t + cols[c] * nTip];
        for (int j = 0; j < zCols; ++j)
          subZ[c + j * nSub] = zPartLocal[cols[c] + j * nCharPart];
      }
      double* rootFreqs = orcScratch.rootFreqs.data();
      for (int k = 0; k < kp; ++k) rootFreqs[k] = 1.0 / kp;
      int stride = nSub * kp;
      size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
      std::fill_n(orcScratch.buf.begin(), bufBytes, 0.0);
      std::fill_n(orcScratch.initFlg.begin(), maxNode + 1, 0u);

      double subLl = pruning_jc_acrv_flat_ecology_raw(
        parPtr, chPtr, nEdge,
        elPtr,
        subStates, nTip, nSub,
        kp, rootFreqs,
        ratesPtr, nCat,
        wPtr, kEco,
        subZ,
        phiPtr, mode, refE, gammaE,
        orcScratch.buf.data(), orcScratch.initFlg.data(), stride);

      if (data.codingType == 1) {
        int* zVec = orcScratch.zVec.data();
        for (int c = 0; c < nSub; ++c) {
          for (int j = 0; j < zCols; ++j)
            zVec[j] = subZ[c + j * nSub];
          double pConst = const_site_prob_jc_eco_single_raw(
            parPtr, chPtr, nEdge,
            elPtr,
            nTip, kp,
            ratesPtr, nCat,
            wPtr, kEco,
            zVec,
            phiPtr, mode, refE, gammaE,
            orcScratch);
          subLl -= std::log(1.0 - pConst);
        }
      }

      if (data.relabel) {
        for (int c = 0; c < nSub; ++c) {
          int gi = gciPtr[cols[c]];
          subLl += mk_prime_relabel_log(kp, data.kObs[gi]);
        }
      }
      ll += subLl;
    }
  }

  if (data.codingType == 2) {
    // Cannot call Rcpp::stop from a parallel region (longjmp out of a worker
    // thread is UB).  cpp_log_likelihood_ecology guards this case before the
    // parallel region; if reached here outside that path, return R_NegInf as
    // a safe sentinel — the existing tests verify codingType==1 paths.
    return R_NegInf;
  }

  return ll;
}


// T-017: existing Rcpp-arg entry point is preserved as a thin wrapper.  All
// move-handler callers in mcmc.cpp pass through here, untouched.  Const-ref
// args avoid the Vector copy-constructor (precious-object list mutation).
double cpp_partition_log_likelihood_ecology(
    const McmcData& data, int partIdx,
    const Rcpp::IntegerVector& parent, const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss, double rateNeo,
    const Rcpp::NumericVector& phi,
    const Rcpp::IntegerMatrix& zMatrix,
    const Rcpp::NumericMatrix& wEdge,
    const std::vector<double>& gammaE,
    const Rcpp::NumericVector& rates) {
  if (data.codingType == 2) {
    // Preserve the legacy diagnostic for the Rcpp entry (serial only).
    stop("informative coding is not yet supported under ecologyAware");
  }
  int nCharsTotal = zMatrix.nrow();  // global #characters
  return cpp_partition_log_likelihood_ecology_raw(
    data, partIdx,
    INTEGER(parent), INTEGER(child), parent.size(),
    REAL(edgeLen),
    INTEGER(kPrime),
    rateLoss, rateNeo,
    REAL(phi), phi.size(),
    INTEGER(zMatrix), nCharsTotal,
    REAL(wEdge),
    gammaE,
    REAL(rates), rates.size());
}


// Compute the kEco-vector gammaE[s] = pi0 + (1-pi0) * (theta_e * phi_e +
// (1-theta_e) / phi_e) with gammaE[refE] = 1. Cheap (O(kEco)) so it is
// recomputed per call rather than cached. Exposed so move handlers can
// build gammaE without going through the full orchestrator.
void compute_gamma_e_ecology(
    const McmcData& data,
    const Rcpp::NumericVector& phi,
    double pi0,
    const Rcpp::NumericVector& theta,
    std::vector<double>& gammaE) {
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  gammaE.assign(kEco, 1.0);
  for (int s = 0; s < kEco; ++s) {
    if (s == refE) { gammaE[s] = 1.0; continue; }
    int j = (s < refE) ? s : (s - 1);
    double phi_s = (mode == 0) ? phi[0] : phi[s];
    double th = (j >= 0 && j < theta.size()) ? theta[j] : 0.5;
    gammaE[s] = gamma_e_compute(pi0, th, phi_s);
  }
}


// ---------------------------------------------------------------------------
// cpp_log_likelihood_ecology: full-data orchestrator under the ecology mixture
// ---------------------------------------------------------------------------
//
// Mirrors the R-side .MkpEcologyLogLikelihood (R/likelihood.R): recompute
// ecology marginals + edge weights from the current tree, then iterate
// partitions invoking the ecology-aware pruning variants.
//
// T-010: refactored to call cpp_partition_log_likelihood_ecology per
// partition, preserving the per-partition summation order. The total is
// bit-identical to the pre-refactor monolithic body.

// (EcoWorkItem and EcoThreadScratch moved earlier — see above the
// const_site_prob_*_raw functions — so those functions can take
// EcoThreadScratch& parameters.)


// T-017: per-work-item core.  Bit-identical to the corresponding branch in
// cpp_partition_log_likelihood_ecology_raw, but operating on a single
// neo/known partition or a single transformational subgroup.  Pure
// std::vector scratch; no Rcpp constructors, no R-API calls.  Safe to invoke
// from inside an OpenMP parallel region.
static double compute_eco_work_item_ll(
    const EcoWorkItem& w,
    const McmcData& data,
    const int* parPtr, const int* chPtr, int nEdge,
    const double* elPtr,
    const int* kPrimePtr,
    double rateLoss, double rateNeo,
    const double* phiPtr, int phiLen,
    const int* zMatrixGlobalPtr, int nCharsTotal,
    const double* wPtr,
    const std::vector<double>& gammaE,
    const double* ratesPtr, int nCat,
    EcoThreadScratch& scratch) {

  (void)phiLen;

  int nTip = data.nTip;
  int maxNode = 2 * nTip - 1;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int zCols = kEco - 1;

  const PartInfo& part = data.parts[w.partIdx];
  int nCharPart = part.tipStates.ncol();
  const int* partTipsPtr = INTEGER(part.tipStates);
  const int* gciPtr      = INTEGER(part.globalCharIdx);

  // T-017 Phase 1b: reuse pre-allocated scratch buffers, zero-fill the
  // prefixes we'll touch.  Eliminates the per-call malloc/free that was
  // contending on glibc's global lock with 4 threads.

  if (w.subgroupKp < 0) {
    // Whole partition (neo or known).  Build local zPart from the global
    // zMatrix using globalCharIdx.
    int* zPartLocal = scratch.zPartLocal.data();
    for (int c = 0; c < nCharPart; ++c) {
      int gi = gciPtr[c];
      for (int j = 0; j < zCols; ++j)
        zPartLocal[c + j * nCharPart] = zMatrixGlobalPtr[gi + j * nCharsTotal];
    }
    double ll = 0.0;
    if (part.type == 0) {
      // Neomorphic.
      double* neoEl = scratch.neoEl.data();
      for (int i = 0; i < nEdge; ++i) neoEl[i] = elPtr[i] * rateNeo;
      double* rootFreqs = scratch.rootFreqs.data();
      rootFreqs[0] = rateLoss / (1.0 + rateLoss);
      rootFreqs[1] = 1.0 / (1.0 + rateLoss);
      int stride = nCharPart * 2;
      size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
      std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
      std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);
      ll = pruning_mkn_acrv_flat_ecology_raw(
        parPtr, chPtr, nEdge,
        neoEl,
        partTipsPtr, nTip, nCharPart,
        rateLoss, rootFreqs,
        ratesPtr, nCat,
        wPtr, kEco,
        zPartLocal,
        phiPtr, mode, refE, gammaE,
        scratch.buf.data(), scratch.initFlg.data(), stride);
      if (data.codingType == 1) {
        int* zVec = scratch.zVec.data();
        for (int c = 0; c < nCharPart; ++c) {
          for (int j = 0; j < zCols; ++j)
            zVec[j] = zPartLocal[c + j * nCharPart];
          double pConst = const_site_prob_mkn_eco_single_raw(
            parPtr, chPtr, nEdge,
            neoEl,
            nTip,
            rateLoss,
            ratesPtr, nCat,
            wPtr, kEco,
            zVec,
            phiPtr, mode, refE, gammaE,
            scratch);
          ll -= std::log(1.0 - pConst);
        }
      }
    } else { // part.type == 2
      int kStates = part.k;
      double* rootFreqs = scratch.rootFreqs.data();
      for (int k = 0; k < kStates; ++k) rootFreqs[k] = 1.0 / kStates;
      int stride = nCharPart * kStates;
      size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
      std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
      std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);
      ll = pruning_jc_acrv_flat_ecology_raw(
        parPtr, chPtr, nEdge,
        elPtr,
        partTipsPtr, nTip, nCharPart,
        kStates, rootFreqs,
        ratesPtr, nCat,
        wPtr, kEco,
        zPartLocal,
        phiPtr, mode, refE, gammaE,
        scratch.buf.data(), scratch.initFlg.data(), stride);
      if (data.codingType == 1) {
        int* zVec = scratch.zVec.data();
        for (int c = 0; c < nCharPart; ++c) {
          for (int j = 0; j < zCols; ++j)
            zVec[j] = zPartLocal[c + j * nCharPart];
          double pConst = const_site_prob_jc_eco_single_raw(
            parPtr, chPtr, nEdge,
            elPtr,
            nTip, kStates,
            ratesPtr, nCat,
            wPtr, kEco,
            zVec,
            phiPtr, mode, refE, gammaE,
            scratch);
          ll -= std::log(1.0 - pConst);
        }
      }
    }
    return ll;
  }

  // Transformational subgroup: one kPrime value (= w.subgroupKp).
  int kp = w.subgroupKp;
  // Collect the local column indices for this subgroup.  The order MUST
  // match the std::map iteration order in the serial reference (ascending
  // kPrime, then ascending insertion order of cols within a kp); the byKp
  // map ordering in the raw routine sorts by kp.  Within a kp the cols
  // vector preserves the order chars were visited (ascending c), so we
  // iterate c = 0..nCharPart-1 here too.
  int* cols = scratch.cols.data();
  int nSub = 0;
  for (int c = 0; c < nCharPart; ++c) {
    int gi = gciPtr[c];
    if (kPrimePtr[gi] == kp) cols[nSub++] = c;
  }

  int* subStates = scratch.subStates.data();
  int* subZ      = scratch.subZ.data();
  for (int c = 0; c < nSub; ++c) {
    for (int t = 0; t < nTip; ++t)
      subStates[t + c * nTip] = partTipsPtr[t + cols[c] * nTip];
    int gi = gciPtr[cols[c]];
    for (int j = 0; j < zCols; ++j)
      subZ[c + j * nSub] = zMatrixGlobalPtr[gi + j * nCharsTotal];
  }
  double* rootFreqs = scratch.rootFreqs.data();
  for (int k = 0; k < kp; ++k) rootFreqs[k] = 1.0 / kp;
  int stride = nSub * kp;
  size_t bufBytes = static_cast<size_t>(maxNode + 1) * stride;
  std::fill_n(scratch.buf.begin(), bufBytes, 0.0);
  std::fill_n(scratch.initFlg.begin(), maxNode + 1, 0u);

  double subLl = pruning_jc_acrv_flat_ecology_raw(
    parPtr, chPtr, nEdge,
    elPtr,
    subStates, nTip, nSub,
    kp, rootFreqs,
    ratesPtr, nCat,
    wPtr, kEco,
    subZ,
    phiPtr, mode, refE, gammaE,
    scratch.buf.data(), scratch.initFlg.data(), stride);

  if (data.codingType == 1) {
    int* zVec = scratch.zVec.data();
    for (int c = 0; c < nSub; ++c) {
      for (int j = 0; j < zCols; ++j)
        zVec[j] = subZ[c + j * nSub];
      double pConst = const_site_prob_jc_eco_single_raw(
        parPtr, chPtr, nEdge,
        elPtr,
        nTip, kp,
        ratesPtr, nCat,
        wPtr, kEco,
        zVec,
        phiPtr, mode, refE, gammaE,
        scratch);
      subLl -= std::log(1.0 - pConst);
    }
  }

  if (data.relabel) {
    for (int c = 0; c < nSub; ++c) {
      int gi = gciPtr[cols[c]];
      subLl += mk_prime_relabel_log(kp, data.kObs[gi]);
    }
  }
  return subLl;
}


double cpp_log_likelihood_ecology(
    const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0,
    const NumericVector& theta) {

  int nTip = data.nTip;
  int nEdge = parent.size();
  int kEco = data.ecology.kEcology;

  // Informative-coding guard: cannot longjmp from a parallel worker thread,
  // so reject here in the serial preamble (no parallel partition currently
  // supports informative coding).
  if (data.codingType == 2) {
    stop("informative coding is not yet supported under ecologyAware");
  }

  // v2: compute gammaE per ecology state (1.0 at refEcology).
  std::vector<double> gammaE;
  compute_gamma_e_ecology(data, phi, pi0, theta, gammaE);

  bool useAcrv = (rateLogSd > 0.0);
  NumericVector rates = useAcrv
    ? cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ)
    : NumericVector(1, 1.0);

  // Compute ecology marginals + per-edge weights once for this call.
  std::vector<double> margFlat;
  compute_ecology_node_marginals(
    INTEGER(parent), INTEGER(child), nEdge,
    REAL(edgeLen),
    INTEGER(data.ecology.tipStates), nTip,
    kEco, 1.0,
    margFlat
  );
  NumericMatrix wEdge(nEdge, kEco);
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e], ch = child[e];
    for (int s = 0; s < kEco; ++s) {
      wEdge(e, s) = 0.5 * (margFlat[par * kEco + s] + margFlat[ch * kEco + s]);
    }
  }

  // T-017: build the flat work-item list serially.  One item per neo/known
  // partition; one item per kPrime subgroup inside each transformational
  // partition.  Order in the list determines the post-parallel serial merge
  // order; this is fixed regardless of thread schedule so the final sum is
  // bit-identical to the serial reference.
  //
  // The serial reference summed partitions in index order, and within a
  // transformational partition summed kPrime subgroups in ascending-kPrime
  // order (std::map iteration order).  We reproduce that exact ordering
  // here.
  std::vector<EcoWorkItem> work;
  work.reserve(data.parts.size() * 2);
  for (int piIdx = 0; piIdx < (int)data.parts.size(); ++piIdx) {
    const PartInfo& part = data.parts[piIdx];
    if (part.type == 0 || part.type == 2) {
      work.push_back({piIdx, -1});
    } else {
      // Enumerate kPrime values present in this partition's chars.
      std::map<int, int> kpSet;  // kp -> count (count is unused; map sorts).
      int nCharPart = part.tipStates.ncol();
      const int* gciPtr = INTEGER(part.globalCharIdx);
      const int* kpPtr  = INTEGER(kPrime);
      for (int c = 0; c < nCharPart; ++c) {
        int gi = gciPtr[c];
        ++kpSet[kpPtr[gi]];
      }
      for (auto& kv : kpSet) {
        work.push_back({piIdx, kv.first});
      }
    }
  }

  // Pre-allocated per-item result buffer.  Each parallel iteration writes
  // exactly one slot (no contention).  Serial sum after the parallel region
  // → bit-identical to the serial reference.
  int nWork = static_cast<int>(work.size());
  std::vector<double> partialLL(nWork, 0.0);

  // T-017 Phase 1b: per-thread scratch workspace, persisted across calls
  // via function-static.  Eliminates the per-call malloc that was the
  // remaining bottleneck after Phase 1a (~300 KB `buf` × 4 threads
  // contending on glibc's global lock).
  //
  // Function-static is thread-safe for this function's current use
  // pattern (chain dispatch is serial — only one chain calls
  // cpp_log_likelihood_ecology at a time, even under PT).  Phase 2
  // (threaded move dispatch) would need to revisit this — e.g., move
  // scratch onto McmcState so it becomes per-chain.
  int maxNCharPart = 0, maxNSub = 0, maxKStates = 2, maxStride = 0;
  for (const EcoWorkItem& w : work) {
    const PartInfo& part = data.parts[w.partIdx];
    int nCharPart = part.tipStates.ncol();
    if (nCharPart > maxNCharPart) maxNCharPart = nCharPart;
    int stride;
    if (w.subgroupKp < 0) {
      int kStates = (part.type == 0) ? 2 : part.k;
      if (kStates > maxKStates) maxKStates = kStates;
      stride = nCharPart * kStates;
    } else {
      int kp = w.subgroupKp;
      if (kp > maxKStates) maxKStates = kp;
      int nSub = 0;
      const int* gciPtr = INTEGER(part.globalCharIdx);
      const int* kpPtr  = INTEGER(kPrime);
      for (int c = 0; c < nCharPart; ++c)
        if (kpPtr[gciPtr[c]] == kp) ++nSub;
      if (nSub > maxNSub) maxNSub = nSub;
      stride = nSub * kp;
    }
    if (stride > maxStride) maxStride = stride;
  }
  if (maxNSub < 1) maxNSub = 1;
  int maxNodeLocal = 2 * nTip - 1;
  int zColsLocal = std::max(1, kEco - 1);
  int nThreads = 1;
#ifdef _OPENMP
  nThreads = omp_get_max_threads();
#endif

  static std::vector<EcoThreadScratch> scratch;
  if ((int)scratch.size() < nThreads) scratch.resize(nThreads);
  // Grow-only resize: if current capacity insufficient for ANY dimension,
  // re-allocate.  Otherwise reuse.
  for (auto& s : scratch) {
    if (s.maxNode < maxNodeLocal || s.maxStride < maxStride ||
        (int)s.zPartLocal.size() < std::max(maxNCharPart, maxNSub) * zColsLocal ||
        (int)s.subStates.size() < nTip * maxNSub ||
        (int)s.neoEl.size() < nEdge ||
        (int)s.rootFreqs.size() < maxKStates ||
        s.maxNTip < nTip) {
      s.allocate(maxNodeLocal, maxStride, std::max(maxNCharPart, maxNSub),
                 maxNSub, zColsLocal, nTip, nEdge, maxKStates);
    }
  }

  const int* parPtr      = INTEGER(parent);
  const int* chPtr       = INTEGER(child);
  const double* elPtr    = REAL(edgeLen);
  const int* kPrimePtr   = INTEGER(kPrime);
  const double* phiPtr   = REAL(phi);
  const int* zMatrixPtr  = INTEGER(zMatrix);
  int nCharsTotal        = zMatrix.nrow();
  const double* wPtr     = REAL(wEdge);
  const double* ratesPtr = REAL(rates);
  int nCat               = rates.size();
  int phiLen             = phi.size();

  // T-017: parallelise the flat list with dynamic scheduling.  Per-item
  // result writes into a unique slot of partialLL, so no contention.  Each
  // worker only calls thread-safe std::vector operations + the inner pruner
  // raw entry points.  No Rcpp constructors, no R-API calls, no RNG.
  //
  // DO NOT replace this with reduction(+:totalLoglik): FP add isn't
  // associative and the reduction tree order depends on thread schedule,
  // which would break bit-identity vs the serial reference.

  // Env-gated diagnostic (T-017): prints thread + work-item info ONCE per
  // process when MKPRIME_T017_DIAG=1.  Helps verify the parallel region
  // is engaging.  Inert when env unset.
#ifdef _OPENMP
  {
    static bool t017_diag_done = false;
    if (!t017_diag_done) {
      const char* env = std::getenv("MKPRIME_T017_DIAG");
      if (env != nullptr && env[0] == '1') {
        t017_diag_done = true;
        std::fprintf(stderr,
          "[T017] cpp_log_likelihood_ecology: nWork=%d, "
          "omp_get_max_threads=%d, _OPENMP=%d\n",
          nWork, omp_get_max_threads(), _OPENMP);
      }
    }
  }
#endif

#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
#endif
  for (int wi = 0; wi < nWork; ++wi) {
    int tid = 0;
#ifdef _OPENMP
    tid = omp_get_thread_num();
#endif
    partialLL[wi] = compute_eco_work_item_ll(
      work[wi], data,
      parPtr, chPtr, nEdge,
      elPtr,
      kPrimePtr,
      rateLoss, rateNeo,
      phiPtr, phiLen,
      zMatrixPtr, nCharsTotal,
      wPtr,
      gammaE,
      ratesPtr, nCat,
      scratch[tid]);
  }

  // Serial merge.  Order = work-item list order = serial reference order →
  // bit-identical to the pre-T-017 serial sum.
  double totalLoglik = 0.0;
  for (int wi = 0; wi < nWork; ++wi) totalLoglik += partialLL[wi];
  return totalLoglik;
}


// R-callable wrapper for cpp_log_likelihood_ecology, used for tests that
// compare C++ to the R-level .MkpEcologyLogLikelihood orchestrator.
//
// [[Rcpp::export(.CppLogLikelihoodEcology)]]
double CppLogLikelihoodEcology(
    SEXP dataPtr,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    IntegerVector kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix,
    double pi0 = 0.0,
    Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  int kEco = data.ecology.kEcology;
  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(std::max(0, kEco - 1), 0.5);
  return cpp_log_likelihood_ecology(
    data, parent, child, edgeLen, kPrime,
    rateLoss, rateLogSd, rateNeo, phi, zMatrix,
    pi0, thetaVec);
}


// ---------------------------------------------------------------------------
// per_char_log_lik_ecology: per-character ecology mixture log-likelihood.
//
// Used by the z Gibbs sweep (Phase 3f) to evaluate the conditional for a
// single (character, z_{c,s}) cell without re-running the full data
// orchestrator.  All "expensive" intermediates — ACRV rates, per-edge
// ecology weights — are taken as inputs so the caller can hoist them out
// of the per-cell loop.
//
// Returns the SAME contribution that cpp_log_likelihood_ecology adds for
// this character: the pruning log-likelihood, plus the variable-coding
// correction when codingType == 1, plus the Mk' relabel correction when
// the character is transformational and data.relabel is on.  The sum
// across all characters equals cpp_log_likelihood_ecology(...).
//
// `globalCharIdx` is 0-based.
// ---------------------------------------------------------------------------

double per_char_log_lik_ecology(
    const McmcData& data,
    int globalCharIdx,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    int kp,
    double rateLoss, double rateNeo,
    const NumericVector& rates,
    NumericVector phi,
    IntegerVector zRow,        // v2: length (kEco - 1)
    const NumericMatrix& wEdge,// nEdge x kEco
    int refEcology,
    const std::vector<double>& gammaE,
    GibbsZWorkspace& ws        // caller-owned workspace; no per-call alloc
) {
  int nTip   = data.nTip;
  int nEdge  = parent.size();
  int maxNode = 2 * nTip - 1;
  int kEco   = data.ecology.kEcology;
  int mode   = data.magnitudeMode;
  int pi     = data.charToPartition[globalCharIdx];
  if (pi < 0) return 0.0;
  const PartInfo& part = data.parts[pi];

  // Find the local column index for this global char.
  int localCol = -1;
  for (int c = 0; c < part.globalCharIdx.size(); ++c) {
    if (part.globalCharIdx[c] == globalCharIdx) { localCol = c; break; }
  }
  if (localCol < 0) return 0.0;

  // Single-column tipStates and single-row zPart.
  IntegerMatrix tipStates(nTip, 1);
  for (int t = 0; t < nTip; ++t) tipStates(t, 0) = part.tipStates(t, localCol);
  int zCols = kEco - 1;
  IntegerMatrix zPart(1, zCols);
  for (int j = 0; j < zCols; ++j) zPart(0, j) = zRow[j];

  double ll;

  if (part.type == 0) {
    // Neomorphic (MkN)
    NumericVector neoEl(nEdge);
    for (int i = 0; i < nEdge; ++i) neoEl[i] = edgeLen[i] * rateNeo;
    NumericVector rootFreqs = mkn_stationary_local(rateLoss);
    int stride = 1 * 2;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_mkn_acrv_flat_ecology(
      parent, child, neoEl, tipStates,
      rateLoss, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      refEcology, gammaE,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_mkn_eco_single(
        parent, child, neoEl, nTip,
        rateLoss, rates, wEdge, zRow, phi, mode,
        refEcology, gammaE);
      ll -= std::log(1.0 - pConst);
    }
  } else if (part.type == 2) {
    int kStates = part.k;
    NumericVector rootFreqs(kStates, 1.0 / kStates);
    int stride = 1 * kStates;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology(
      parent, child, edgeLen, tipStates,
      kStates, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      refEcology, gammaE,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_jc_eco_single(
        parent, child, edgeLen, nTip, kStates,
        rates, wEdge, zRow, phi, mode,
        refEcology, gammaE);
      ll -= std::log(1.0 - pConst);
    }
  } else {
    NumericVector rootFreqs(kp, 1.0 / kp);
    int stride = 1 * kp;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology(
      parent, child, edgeLen, tipStates,
      kp, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      refEcology, gammaE,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_jc_eco_single(
        parent, child, edgeLen, nTip, kp,
        rates, wEdge, zRow, phi, mode,
        refEcology, gammaE);
      ll -= std::log(1.0 - pConst);
    }
    if (data.relabel) {
      ll += mk_prime_relabel_log(kp, data.kObs[globalCharIdx]);
    }
  }

  return ll;
}


// Recompute the nEdge x kEco mean-of-marginals weight matrix from the
// current tree + ecology tip states.  Exposed so the Gibbs sweep can
// hoist it out of the per-cell loop.
void recompute_w_edge(
    const McmcData& data,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    NumericMatrix& wEdgeOut
) {
  int nTip = data.nTip;
  int nEdge = parent.size();
  int kEco = data.ecology.kEcology;
  std::vector<double> margFlat;
  compute_ecology_node_marginals(
    INTEGER(parent), INTEGER(child), nEdge,
    REAL(edgeLen),
    INTEGER(data.ecology.tipStates), nTip,
    kEco, 1.0,
    margFlat
  );
  if (wEdgeOut.nrow() != nEdge || wEdgeOut.ncol() != kEco) {
    wEdgeOut = NumericMatrix(nEdge, kEco);
  }
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e], ch = child[e];
    for (int s = 0; s < kEco; ++s) {
      wEdgeOut(e, s) = 0.5 * (margFlat[par * kEco + s] +
                              margFlat[ch * kEco + s]);
    }
  }
}


// R-callable: per-character ecology log-likelihoods.  Invariant tested
// by the Gibbs scaffolding: sum equals cpp_log_likelihood_ecology(...).
//
// [[Rcpp::export(.CppLogLikelihoodEcologyPerChar)]]
NumericVector CppLogLikelihoodEcologyPerChar(
    SEXP dataPtr,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    IntegerVector kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix,
    double pi0 = 0.0,
    Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  if (data.codingType == 2) {
    stop("informative coding is not yet supported under ecologyAware");
  }

  int nChar = data.nChar;
  int kEco  = data.ecology.kEcology;
  int refE  = data.ecology.refEcology;
  int mode  = data.magnitudeMode;
  int zCols = kEco - 1;
  NumericVector out(nChar);

  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(std::max(0, zCols), 0.5);

  // gammaE
  std::vector<double> gammaE(kEco, 1.0);
  for (int s = 0; s < kEco; ++s) {
    if (s == refE) { gammaE[s] = 1.0; continue; }
    int j = (s < refE) ? s : (s - 1);
    double phi_s = (mode == 0) ? phi[0] : phi[s];
    double th = (j >= 0 && j < thetaVec.size()) ? thetaVec[j] : 0.5;
    gammaE[s] = gamma_e_compute(pi0, th, phi_s);
  }

  NumericMatrix wEdge(parent.size(), kEco);
  recompute_w_edge(data, parent, child, edgeLen, wEdge);

  NumericVector rates = (rateLogSd > 0.0)
    ? cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ)
    : NumericVector(1, 1.0);

  GibbsZWorkspace ws;
  IntegerVector zRow(zCols);
  for (int c = 0; c < nChar; ++c) {
    for (int j = 0; j < zCols; ++j) zRow[j] = zMatrix(c, j);
    int kp = (data.charToPartition[c] >= 0 &&
              data.parts[data.charToPartition[c]].type == 1)
              ? kPrime[c] : 2;
    out[c] = per_char_log_lik_ecology(
      data, c, parent, child, edgeLen,
      kp, rateLoss, rateNeo, rates, phi, zRow, wEdge,
      refE, gammaE, ws);
  }
  return out;
}


// ===========================================================================
// T-012: cached ecology pruner (populate + cached eval)
// ===========================================================================
//
// populate_eco_cache_full builds the per-unit storage and runs a full
// downpass writing CLs at every internal node, snapshots the (wEdge,
// scalar param) baseline, and marks topoValid / structureValid / nodeValid
// = true for every internal node.
//
// eco_cache_total_loglik computes the full ecology log-likelihood from the
// cached root CLs, applying ascertainment correction and (for transformational
// units) the Mk' relabel correction.  Bit-identical to
// cpp_log_likelihood_ecology when the cache was populated against the same
// inputs.

static IntegerMatrix make_sub_zmat_eco(const IntegerMatrix& zMatrix,
                                       const std::vector<int>& globalCharIdx,
                                       int zCols) {
  int nSub = (int)globalCharIdx.size();
  IntegerMatrix sub(nSub, zCols);
  for (int c = 0; c < nSub; ++c) {
    int gi = globalCharIdx[c];
    for (int j = 0; j < zCols; ++j) sub(c, j) = zMatrix(gi, j);
  }
  return sub;
}


// Populate the entire cache (units + tip CLs + full downpass + snapshot).
void populate_eco_cache_full(
    EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0, const NumericVector& theta,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE
) {
  int nTip = data.nTip;
  int maxNode = 2 * nTip - 1;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int nCat = (rateLogSd > 0.0) ? data.nCat : 1;

  cache.nTip = nTip;
  cache.maxNode = maxNode;
  cache.useAcrv = (rateLogSd > 0.0);
  cache.cachedRateLogSd = rateLogSd;

  NumericVector ratesV = (rateLogSd > 0.0)
    ? cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ)
    : NumericVector(1, 1.0);
  cache.rates.assign(ratesV.begin(), ratesV.end());

  cache.topo.build(parent, child, edgeLen, nTip);
  cache.nodeValid.assign(maxNode + 1, 0u);
  for (int t = 1; t <= nTip; ++t) cache.nodeValid[t] = 1u;

  build_eco_cache_units(cache, data, kPrime);
  int zCols = kEco - 1;

  for (auto& unit : cache.units) {
    unit.allocate(maxNode, nCat);
    if (unit.isMkN) {
      unit.rootFreqs = { rateLoss / (1.0 + rateLoss), 1.0 / (1.0 + rateLoss) };
      unit.rateScale = rateNeo;
    } else {
      unit.rateScale = 1.0;
      unit.rootFreqs.assign(unit.kStates, 1.0 / unit.kStates);
    }
    unit.init_tips(nTip);

    IntegerMatrix unitZ = make_sub_zmat_eco(zMatrix, unit.globalCharIdx, zCols);
    eco_full_downpass_unit(
      unit, parent, child, edgeLen, nTip,
      kEco, refE, mode, phi, wEdge, unitZ,
      gammaE, cache.rates, rateLoss
    );
  }

  for (int n = nTip + 1; n <= maxNode; ++n) cache.nodeValid[n] = 1u;
  cache.topoValid = true;
  cache.structureValid = true;
  snapshot_eco_inputs(cache, wEdge, edgeLen, phi, pi0, theta, gammaE);
}


// Compute total log-likelihood from cached CLs (incl. ascertainment +
// relabel correction).
double eco_cache_total_loglik(
    const EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const IntegerVector& kPrime, double rateLoss, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE
) {
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int nTip = data.nTip;
  int root = nTip + 1;
  int zCols = kEco - 1;
  int nCat = (int)cache.rates.size();
  if (nCat < 1) nCat = 1;
  NumericVector ratesV(cache.rates.size());
  for (size_t i = 0; i < cache.rates.size(); ++i) ratesV[i] = cache.rates[i];

  double total = 0.0;
  for (const auto& unit : cache.units) {
    double ll = eco_unit_root_loglik(unit, root, nCat);
    if (!R_FINITE(ll)) return R_NegInf;

    if (data.codingType == 1) {
      IntegerVector zVec(zCols);
      const int nSub = unit.nChar;
      for (int c = 0; c < nSub; ++c) {
        int gi = unit.globalCharIdx[c];
        for (int j = 0; j < zCols; ++j) zVec[j] = zMatrix(gi, j);
        double pConst;
        if (unit.isMkN) {
          NumericVector neoEl(parent.size());
          for (int e = 0; e < parent.size(); ++e) neoEl[e] = edgeLen[e] * rateNeo;
          pConst = const_site_prob_mkn_eco_single(
            parent, child, neoEl, nTip, rateLoss, ratesV,
            wEdge, zVec, phi, mode, refE, gammaE);
        } else {
          pConst = const_site_prob_jc_eco_single(
            parent, child, edgeLen, nTip, unit.kStates, ratesV,
            wEdge, zVec, phi, mode, refE, gammaE);
        }
        if (pConst > 0.0 && pConst < 1.0)
          ll -= std::log(1.0 - pConst);
      }
    } else if (data.codingType == 2) {
      stop("informative coding is not yet supported under ecologyAware");
    }

    if (unit.doRelabel) {
      for (int c = 0; c < unit.nChar; ++c) {
        int gi = unit.globalCharIdx[c];
        ll += mk_prime_relabel_log(unit.kStates, data.kObs[gi]);
      }
    }

    total += ll;
  }

  (void)kPrime;
  return total;
}


// Partial eval for an NNI move.  Assumes cache.topo + (parent, child) have
// already been updated to the proposed topology.  wEdge has been computed
// against the proposed topology too.  Returns the partial-eval logLik;
// `fallback` is set to true if the dirty set is too large (caller should
// then use a full eval instead).
//
// `dirtyFracThreshold` (default 0.70) is the fraction-of-internals at and
// beyond which partial-eval gives up.  Caller can raise this experimentally.
double eco_cache_partial_eval_nni(
    EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0, const NumericVector& theta,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE,
    int v, int u,
    EcoDirtyScratch& scratch,
    int& dirtyCount,
    bool& fallback,
    int* nWedgeDirtyEdges = nullptr,
    double dirtyFracThreshold = 0.70
) {
  fallback = false;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int nTip = data.nTip;
  int maxNode = cache.maxNode;
  int zCols = kEco - 1;

  std::vector<int> topoPath = find_dirty_nni_eco(cache.topo, v, u);
  std::vector<int> wDirty = detect_dirty_child_nodes(cache, child, wEdge, edgeLen);
  if (nWedgeDirtyEdges) *nWedgeDirtyEdges = (int)wDirty.size();

  std::vector<uint8_t> mark(maxNode + 1, 0);
  for (int n : topoPath) if (n >= 0) mark[n] = 1;
  for (int ch : wDirty) {
    int p = cache.topo.parentNode[ch];
    while (p >= 0) {
      if (mark[p]) break;
      mark[p] = 1;
      p = cache.topo.parentNode[p];
    }
  }

  std::vector<std::pair<int,int>> depthNode;
  for (int n = nTip + 1; n <= maxNode; ++n) {
    if (!mark[n]) continue;
    int d = 0;
    for (int x = n; x >= 0; x = cache.topo.parentNode[x]) ++d;
    depthNode.emplace_back(d, n);
  }
  std::sort(depthNode.begin(), depthNode.end(),
            [](const auto& a, const auto& b) { return a.first > b.first; });
  std::vector<int> dirty;
  dirty.reserve(depthNode.size());
  for (auto& dn : depthNode) dirty.push_back(dn.second);
  dirtyCount = (int)dirty.size();

  int nInternal = maxNode - nTip;
  if (nInternal > 0 && (double)dirtyCount > dirtyFracThreshold * (double)nInternal) {
    fallback = true;
    return R_NegInf;
  }

  save_dirty_eco_cls(cache, scratch, dirty);

  for (auto& unit : cache.units) {
    IntegerMatrix unitZ = make_sub_zmat_eco(zMatrix, unit.globalCharIdx, zCols);
    eco_recompute_dirty_unit(
      unit, cache,
      parent, child, edgeLen,
      kEco, refE, mode,
      phi, wEdge, unitZ, gammaE, cache.rates, rateLoss,
      dirty
    );
  }

  double ll = eco_cache_total_loglik(
    cache, data, parent, child, edgeLen, kPrime,
    rateLoss, rateNeo, phi, zMatrix, wEdge, gammaE
  );

  cache.diagPartialEvalCount++;
  cache.diagCachedSubtreeCount += dirtyCount;

  (void)pi0; (void)theta;
  return ll;
}


// R-callable wrapper: build a cache, populate it, return total logLik.
// Used in tests to verify bit-identity vs the legacy pruner.
//
// [[Rcpp::export(.CppLogLikelihoodEcologyCached)]]
double CppLogLikelihoodEcologyCached(
    SEXP dataPtr,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    IntegerVector kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix,
    double pi0 = 0.0,
    Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  int kEco = data.ecology.kEcology;
  int zCols = kEco - 1;
  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(std::max(0, zCols), 0.5);

  std::vector<double> gammaE;
  compute_gamma_e_ecology(data, phi, pi0, thetaVec, gammaE);

  NumericMatrix wEdge(parent.size(), kEco);
  recompute_w_edge(data, parent, child, edgeLen, wEdge);

  EcoCLCache cache;
  populate_eco_cache_full(cache, data, parent, child, edgeLen, kPrime,
                          rateLoss, rateLogSd, rateNeo,
                          phi, zMatrix, pi0, thetaVec, wEdge, gammaE);

  return eco_cache_total_loglik(cache, data, parent, child, edgeLen,
                                 kPrime, rateLoss, rateNeo,
                                 phi, zMatrix, wEdge, gammaE);
}


// T-013: R-callable partial-eval-NNI wrapper.
//
// Populate a fresh cache on (parentOld, childOld), then drive the partial-eval
// path for an NNI swap that produces (parentNew, childNew).  Returns a list
// with:
//   partial_ll  — log-likelihood from eco_cache_partial_eval_nni
//   dirty_count — # nodes recomputed
//   fallback    — heuristic flag if dirty set was too large
//   restored_ll — after restore_dirty_eco_cls, total loglik against OLD
//                 topology (rollback regression: should equal a fresh full
//                 eval on OLD).
//
// [[Rcpp::export(.CppPartialEvalEcologyNNI)]]
Rcpp::List CppPartialEvalEcologyNNI(
    SEXP dataPtr,
    IntegerVector parentOld, IntegerVector childOld,
    IntegerVector parentNew, IntegerVector childNew,
    NumericVector edgeLen,
    IntegerVector kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix,
    double pi0,
    Rcpp::Nullable<Rcpp::NumericVector> theta,
    int vNode, int uNode,
    double dirtyFracThreshold = 0.70
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  int kEco  = data.ecology.kEcology;
  int zCols = kEco - 1;
  int nTip  = data.nTip;
  Rcpp::NumericVector thetaVec = theta.isNotNull()
    ? Rcpp::NumericVector(theta) : Rcpp::NumericVector(std::max(0, zCols), 0.5);

  std::vector<double> gammaE;
  compute_gamma_e_ecology(data, phi, pi0, thetaVec, gammaE);

  // 1. Populate cache on OLD topology (this snapshots wEdge_old etc.)
  NumericMatrix wEdgeOld(parentOld.size(), kEco);
  recompute_w_edge(data, parentOld, childOld, edgeLen, wEdgeOld);

  EcoCLCache cache;
  populate_eco_cache_full(cache, data, parentOld, childOld, edgeLen, kPrime,
                          rateLoss, rateLogSd, rateNeo,
                          phi, zMatrix, pi0, thetaVec, wEdgeOld, gammaE);

  // 2. Build wEdge for NEW topology.
  NumericMatrix wEdgeNew(parentNew.size(), kEco);
  recompute_w_edge(data, parentNew, childNew, edgeLen, wEdgeNew);

  // 3. Update cache.topo to reflect NEW topology (mirrors what production does
  //    in the partial-eval branch before calling eco_cache_partial_eval_nni).
  cache.topo.build(parentNew, childNew, edgeLen, nTip);

  // 4. Partial eval.
  EcoDirtyScratch scratch;
  int dirtyCount = 0;
  int nWedgeDirtyEdges = 0;
  bool fallback = false;
  double partial_ll = eco_cache_partial_eval_nni(
    cache, data, parentNew, childNew, edgeLen, kPrime,
    rateLoss, rateNeo, phi, zMatrix, pi0, thetaVec,
    wEdgeNew, gammaE, vNode, uNode,
    scratch, dirtyCount, fallback,
    &nWedgeDirtyEdges, dirtyFracThreshold);
  int nInternal = (data.nTip - 1);  // for nTip tips, internal-node count = nTip-1

  // 5. Rollback: restore dirty CLs, then evaluate against OLD topology.
  double restored_ll = R_NaReal;
  if (!fallback) {
    restore_dirty_eco_cls(cache, scratch);
    // Restore cache.topo to OLD topology so eco_cache_total_loglik traverses
    // the correct tree.
    cache.topo.build(parentOld, childOld, edgeLen, nTip);
    restored_ll = eco_cache_total_loglik(
      cache, data, parentOld, childOld, edgeLen, kPrime,
      rateLoss, rateNeo, phi, zMatrix, wEdgeOld, gammaE);
  }

  return Rcpp::List::create(
    Rcpp::_["partial_ll"] = partial_ll,
    Rcpp::_["dirty_count"] = dirtyCount,
    Rcpp::_["n_wedge_dirty_edges"] = nWedgeDirtyEdges,
    Rcpp::_["n_internal"] = nInternal,
    Rcpp::_["fallback"] = fallback,
    Rcpp::_["restored_ll"] = restored_ll
  );
}
