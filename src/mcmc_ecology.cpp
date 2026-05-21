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
#include "fast_exp.h"
#include <cmath>
#include <vector>
#include <map>

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
static double pruning_jc_acrv_flat_ecology(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, NumericVector root_freqs,
    NumericVector rate_multipliers,
    NumericMatrix wEdge,
    IntegerMatrix zMat,
    NumericVector phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();
  int kEco  = wEdge.ncol();

  const int* parPtr = INTEGER(parent);
  const int* chPtr  = INTEGER(child);
  const double* elPtr = REAL(edge_length);
  const int* tsPtr   = INTEGER(tip_states);
  const double* rfPtr = REAL(root_freqs);
  const double* rmPtr = REAL(rate_multipliers);
  const double* wPtr  = REAL(wEdge);          // column-major: wPtr[e + s * nEdge]
  const int* zPtr     = INTEGER(zMat);         // column-major: zPtr[c + s * nChar]
  const double* phiPtr = REAL(phi);

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
static double pruning_mkn_acrv_flat_ecology(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    double rate_loss, NumericVector root_freqs,
    NumericVector rate_multipliers,
    NumericMatrix wEdge,
    IntegerMatrix zMat,
    NumericVector phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();
  int kEco  = wEdge.ncol();
  const int kStates = 2;

  const int* parPtr = INTEGER(parent);
  const int* chPtr  = INTEGER(child);
  const double* elPtr = REAL(edge_length);
  const int* tsPtr   = INTEGER(tip_states);
  const double* rfPtr = REAL(root_freqs);
  const double* rmPtr = REAL(rate_multipliers);
  const double* wPtr  = REAL(wEdge);
  const int* zPtr     = INTEGER(zMat);
  const double* phiPtr = REAL(phi);

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

static double const_site_prob_jc_eco_single(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, int nTip, int kStates,
    NumericVector rates,
    NumericMatrix wEdge, const IntegerVector& zVec,
    NumericVector phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE) {

  IntegerMatrix tipStates(nTip, 1);  // all zeros (default-constructed)
  // zVec is length (kEco - 1) — one entry per non-ref ecology in col order.
  IntegerMatrix zMat(1, zVec.size());
  for (int s = 0; s < zVec.size(); ++s) zMat(0, s) = zVec[s];
  NumericVector rootFreqs(kStates, 1.0 / kStates);

  int maxNode = 2 * nTip - 1;
  int stride  = kStates;
  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  double ll = pruning_jc_acrv_flat_ecology(
    parent, child, edgeLen, tipStates,
    kStates, rootFreqs, rates,
    wEdge, zMat, phi, mode,
    refEcology, gammaE,
    buf.data(), initFlg.data(), stride);
  // JC symmetry: P(constant in any state) = kStates * P(all-tips-0).
  return kStates * std::exp(ll);
}


static double const_site_prob_mkn_eco_single(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, int nTip,
    double rateLoss, NumericVector rates,
    NumericMatrix wEdge, const IntegerVector& zVec,
    NumericVector phi, int mode,
    int refEcology,
    const std::vector<double>& gammaE) {

  IntegerMatrix zMat(1, zVec.size());
  for (int s = 0; s < zVec.size(); ++s) zMat(0, s) = zVec[s];
  NumericVector rootFreqs(2);
  rootFreqs[0] = rateLoss / (1.0 + rateLoss);
  rootFreqs[1] = 1.0 / (1.0 + rateLoss);

  int maxNode = 2 * nTip - 1;
  int stride  = 2;
  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  // Pseudo-char "all 0"
  IntegerMatrix tipStates0(nTip, 1);  // zeros
  double ll0 = pruning_mkn_acrv_flat_ecology(
    parent, child, edgeLen, tipStates0,
    rateLoss, rootFreqs, rates,
    wEdge, zMat, phi, mode,
    refEcology, gammaE,
    buf.data(), initFlg.data(), stride);

  // Pseudo-char "all 1"
  IntegerMatrix tipStates1(nTip, 1);
  for (int t = 0; t < nTip; ++t) tipStates1(t, 0) = 1;
  // Re-zero buffers for the second pass
  std::fill(buf.begin(), buf.end(), 0.0);
  std::fill(initFlg.begin(), initFlg.end(), 0u);
  double ll1 = pruning_mkn_acrv_flat_ecology(
    parent, child, edgeLen, tipStates1,
    rateLoss, rootFreqs, rates,
    wEdge, zMat, phi, mode,
    refEcology, gammaE,
    buf.data(), initFlg.data(), stride);

  return std::exp(ll0) + std::exp(ll1);
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

double cpp_partition_log_likelihood_ecology(
    const McmcData& data, int partIdx,
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss, double rateNeo,
    Rcpp::NumericVector phi,
    Rcpp::IntegerMatrix zMatrix,
    const Rcpp::NumericMatrix& wEdge,
    const std::vector<double>& gammaE,
    const Rcpp::NumericVector& rates) {

  int nTip = data.nTip;
  int nEdge = parent.size();
  int maxNode = 2 * nTip - 1;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
  int refE = data.ecology.refEcology;
  int zCols = kEco - 1;

  const PartInfo& part = data.parts[partIdx];
  int nCharPart = part.tipStates.ncol();

  // Subset zMatrix to this partition's characters (global → local indices).
  IntegerMatrix zPart(nCharPart, zCols);
  for (int c = 0; c < nCharPart; ++c) {
    int gi = part.globalCharIdx[c];
    for (int j = 0; j < zCols; ++j) zPart(c, j) = zMatrix(gi, j);
  }

  double ll = 0.0;

  if (part.type == 0) {
    // Neomorphic
    NumericVector neoEl(nEdge);
    for (int i = 0; i < nEdge; ++i) neoEl[i] = edgeLen[i] * rateNeo;
    NumericVector rootFreqs = mkn_stationary_local(rateLoss);
    int stride = nCharPart * 2;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_mkn_acrv_flat_ecology(
      parent, child, neoEl, part.tipStates,
      rateLoss, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      refE, gammaE,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {  // variable
      for (int c = 0; c < nCharPart; ++c) {
        IntegerVector zVec(zCols);
        for (int j = 0; j < zCols; ++j) zVec[j] = zPart(c, j);
        double pConst = const_site_prob_mkn_eco_single(
          parent, child, neoEl, nTip,
          rateLoss, rates, wEdge, zVec, phi, mode,
          refE, gammaE);
        ll -= std::log(1.0 - pConst);
      }
    }
  } else if (part.type == 2) {
    // Known state space
    int kStates = part.k;
    NumericVector rootFreqs(kStates, 1.0 / kStates);
    int stride = nCharPart * kStates;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology(
      parent, child, edgeLen, part.tipStates,
      kStates, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      refE, gammaE,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      for (int c = 0; c < nCharPart; ++c) {
        IntegerVector zVec(zCols);
        for (int j = 0; j < zCols; ++j) zVec[j] = zPart(c, j);
        double pConst = const_site_prob_jc_eco_single(
          parent, child, edgeLen, nTip, kStates,
          rates, wEdge, zVec, phi, mode,
          refE, gammaE);
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
      int gi = part.globalCharIdx[c];
      byKp[kPrime[gi]].push_back(c);
    }
    for (auto& kv : byKp) {
      int kp = kv.first;
      const std::vector<int>& cols = kv.second;
      int nSub = static_cast<int>(cols.size());

      IntegerMatrix subStates(nTip, nSub);
      IntegerMatrix subZ(nSub, zCols);
      for (int c = 0; c < nSub; ++c) {
        for (int t = 0; t < nTip; ++t)
          subStates(t, c) = part.tipStates(t, cols[c]);
        for (int j = 0; j < zCols; ++j)
          subZ(c, j) = zPart(cols[c], j);
      }
      NumericVector rootFreqs(kp, 1.0 / kp);
      int stride = nSub * kp;
      std::vector<double>  buf((maxNode + 1) * stride, 0.0);
      std::vector<uint8_t> initFlg(maxNode + 1, 0u);

      double subLl = pruning_jc_acrv_flat_ecology(
        parent, child, edgeLen, subStates,
        kp, rootFreqs, rates,
        wEdge, subZ, phi, mode,
        refE, gammaE,
        buf.data(), initFlg.data(), stride);

      if (data.codingType == 1) {
        for (int c = 0; c < nSub; ++c) {
          IntegerVector zVec(zCols);
          for (int j = 0; j < zCols; ++j) zVec[j] = subZ(c, j);
          double pConst = const_site_prob_jc_eco_single(
            parent, child, edgeLen, nTip, kp,
            rates, wEdge, zVec, phi, mode,
            refE, gammaE);
          subLl -= std::log(1.0 - pConst);
        }
      }

      if (data.relabel) {
        for (int c = 0; c < nSub; ++c) {
          int gi = part.globalCharIdx[cols[c]];
          subLl += mk_prime_relabel_log(kp, data.kObs[gi]);
        }
      }
      ll += subLl;
    }
  }

  if (data.codingType == 2) {
    stop("informative coding is not yet supported under ecologyAware");
  }

  return ll;
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

double cpp_log_likelihood_ecology(
    const McmcData& data,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix,
    double pi0,
    NumericVector theta) {

  int nTip = data.nTip;
  int nEdge = parent.size();
  int kEco = data.ecology.kEcology;

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

  double totalLoglik = 0.0;
  for (int pi = 0; pi < (int)data.parts.size(); ++pi) {
    totalLoglik += cpp_partition_log_likelihood_ecology(
      data, pi, parent, child, edgeLen, kPrime,
      rateLoss, rateNeo, phi, zMatrix,
      wEdge, gammaE, rates);
  }
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
