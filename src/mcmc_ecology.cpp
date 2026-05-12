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


// Resolve the rate factor for a transformational character with influence
// category z under magnitude mode `mode` (0 = global, 1 = per_ecology).
// phi has length 1 in global mode and kEco in per_ecology mode.
//   z == 0 (none):        factor = 1
//   z == 1 (encouraged):  factor = phi[mode == 0 ? 0 : ecoState]
//   z == 2 (discouraged): factor = 1 / phi[mode == 0 ? 0 : ecoState]
static inline double trans_rate_factor(int z, int ecoState,
                                       const double* phi, int mode) {
  if (z == 0) return 1.0;
  double p = (mode == 0) ? phi[0] : phi[ecoState];
  return (z == 1) ? p : 1.0 / p;
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

      // Precompute (ps, pd) for each (z-factor, ecology state) at this edge.
      // factor 0 (none): identical across ecology states — compute once and
      // duplicate. factor 1 (encouraged) / 2 (discouraged): depend on phi
      // which varies with mode.
      {
        double exp0 = MKP_EXP(-kStates * tBase / km1);
        double ps0  = inv_k + (1.0 - inv_k) * exp0;
        double pd0  = inv_k - inv_k * exp0;
        for (int s = 0; s < kEco; ++s) {
          psFactor[0 * kEco + s] = ps0;
          pdFactor[0 * kEco + s] = pd0;
        }
        if (mode == 0) {
          // Global phi: same encouraged/discouraged values for every ecology.
          double p = phiPtr[0];
          double tE = tBase * p;
          double tD = tBase / p;
          double expE = MKP_EXP(-kStates * tE / km1);
          double expD = MKP_EXP(-kStates * tD / km1);
          double psE = inv_k + (1.0 - inv_k) * expE;
          double pdE = inv_k - inv_k * expE;
          double psD = inv_k + (1.0 - inv_k) * expD;
          double pdD = inv_k - inv_k * expD;
          for (int s = 0; s < kEco; ++s) {
            psFactor[1 * kEco + s] = psE;
            pdFactor[1 * kEco + s] = pdE;
            psFactor[2 * kEco + s] = psD;
            pdFactor[2 * kEco + s] = pdD;
          }
        } else {
          // Per-ecology phi.
          for (int s = 0; s < kEco; ++s) {
            double p = phiPtr[s];
            double tE = tBase * p;
            double tD = tBase / p;
            double expE = MKP_EXP(-kStates * tE / km1);
            double expD = MKP_EXP(-kStates * tD / km1);
            psFactor[1 * kEco + s] = inv_k + (1.0 - inv_k) * expE;
            pdFactor[1 * kEco + s] = inv_k - inv_k * expE;
            psFactor[2 * kEco + s] = inv_k + (1.0 - inv_k) * expD;
            pdFactor[2 * kEco + s] = inv_k - inv_k * expD;
          }
        }
      }

      double* clPar = buf + par * stride;
      double* clCh  = buf + ch  * stride;

      // Per-character propagation with the mixture (p_same, p_diff).
      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          double psMix = 0.0, pdMix = 0.0;
          for (int s = 0; s < kEco; ++s) {
            int z   = zPtr[c + s * nChar];
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
            int z   = zPtr[c + s * nChar];
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
static inline void mkn_rates_for_state(
    int z, int ecoState,
    double rate01_base, double rate10_base,
    const double* phi, int mode,
    double& r01, double& r10) {
  if (z == 0) {
    r01 = rate01_base;
    r10 = rate10_base;
    return;
  }
  double p = (mode == 0) ? phi[0] : phi[ecoState];
  if (z == 1) {
    r01 = rate01_base * p;
    r10 = rate10_base / p;
  } else {
    r01 = rate01_base / p;
    r10 = rate10_base * p;
  }
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

      // Precompute P for each variant (none / enc / dis) on this (cat, edge).
      // 'none' is ecology-independent; enc/dis depend on phi which may be
      // per-ecology.
      {
        double lam0 = r01_base + r10_base;
        double pi0_0 = r10_base / lam0;
        double pi1_0 = r01_base / lam0;
        double exp0 = MKP_EXP(-lam0 * tBase);
        double P00_0 = pi0_0 + pi1_0 * exp0;
        double P01_0 = pi1_0 - pi1_0 * exp0;
        double P10_0 = pi0_0 - pi0_0 * exp0;
        double P11_0 = pi1_0 + pi0_0 * exp0;
        for (int s = 0; s < kEco; ++s) {
          double* P = Pfactor.data() + (0 * kEco + s) * 4;
          P[0] = P00_0; P[1] = P01_0; P[2] = P10_0; P[3] = P11_0;
        }
      }
      for (int v = 1; v <= 2; ++v) {  // 1 = encouraged, 2 = discouraged
        int sEnd = (mode == 0) ? 1 : kEco;
        for (int s = 0; s < sEnd; ++s) {
          double r01, r10;
          mkn_rates_for_state(v, s, r01_base, r10_base, phiPtr, mode, r01, r10);
          double lam  = r01 + r10;
          double pi0  = r10 / lam;
          double pi1  = r01 / lam;
          double ex   = MKP_EXP(-lam * tBase);
          double P00 = pi0 + pi1 * ex;
          double P01 = pi1 - pi1 * ex;
          double P10 = pi0 - pi0 * ex;
          double P11 = pi1 + pi0 * ex;
          if (mode == 0) {
            for (int ss = 0; ss < kEco; ++ss) {
              double* P = Pfactor.data() + (v * kEco + ss) * 4;
              P[0] = P00; P[1] = P01; P[2] = P10; P[3] = P11;
            }
          } else {
            double* P = Pfactor.data() + (v * kEco + s) * 4;
            P[0] = P00; P[1] = P01; P[2] = P10; P[3] = P11;
          }
        }
      }

      double* clPar = buf + par * stride;
      double* clCh  = buf + ch  * stride;

      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          double P00m = 0.0, P01m = 0.0, P10m = 0.0, P11m = 0.0;
          for (int s = 0; s < kEco; ++s) {
            int z   = zPtr[c + s * nChar];
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
            int z   = zPtr[c + s * nChar];
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
    NumericVector phi, int mode
) {
  int nTip    = tipStates.nrow();
  int nChar   = tipStates.ncol();
  int maxNode = 2 * nTip - 1;
  int stride  = nChar * 2;

  if (rootFreqs.size() != 2)
    stop("rootFreqs length must equal 2");
  if (wEdge.nrow() != parent.size())
    stop("wEdge nrow must equal nEdge");
  if (zMat.nrow() != nChar)
    stop("zMat nrow must equal nChar");
  if (zMat.ncol() != wEdge.ncol())
    stop("zMat ncol must equal kEco");
  if (mode != 0 && mode != 1)
    stop("mode must be 0 (global) or 1 (per_ecology)");
  if (mode == 0 && phi.size() != 1)
    stop("global mode requires phi length 1");
  if (mode == 1 && phi.size() != wEdge.ncol())
    stop("per_ecology mode requires phi length kEco");
  if (rateLoss <= 0)
    stop("rateLoss must be positive");

  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  return pruning_mkn_acrv_flat_ecology(
    parent, child, edgeLen, tipStates,
    rateLoss, rootFreqs, rateMultipliers,
    wEdge, zMat, phi, mode,
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
    NumericVector phi, int mode
) {
  int nTip    = tipStates.nrow();
  int nChar   = tipStates.ncol();
  int maxNode = 2 * nTip - 1;
  int stride  = nChar * kStates;

  if (kStates < 2)
    stop("kStates must be >= 2");
  if (rootFreqs.size() != kStates)
    stop("rootFreqs length must equal kStates");
  if (wEdge.nrow() != parent.size())
    stop("wEdge nrow must equal nEdge");
  if (zMat.nrow() != nChar)
    stop("zMat nrow must equal nChar");
  if (zMat.ncol() != wEdge.ncol())
    stop("zMat ncol must equal kEco (wEdge ncol)");
  if (mode != 0 && mode != 1)
    stop("mode must be 0 (global) or 1 (per_ecology)");
  if (mode == 0 && phi.size() != 1)
    stop("global mode requires phi length 1");
  if (mode == 1 && phi.size() != wEdge.ncol())
    stop("per_ecology mode requires phi length kEco");

  std::vector<double>  buf((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> initFlg(maxNode + 1, 0u);

  return pruning_jc_acrv_flat_ecology(
    parent, child, edgeLen, tipStates,
    kStates, rootFreqs, rateMultipliers,
    wEdge, zMat, phi, mode,
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
    NumericVector phi, int mode) {

  IntegerMatrix tipStates(nTip, 1);  // all zeros (default-constructed)
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
    buf.data(), initFlg.data(), stride);
  // JC symmetry: P(constant in any state) = kStates * P(all-tips-0).
  return kStates * std::exp(ll);
}


static double const_site_prob_mkn_eco_single(
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen, int nTip,
    double rateLoss, NumericVector rates,
    NumericMatrix wEdge, const IntegerVector& zVec,
    NumericVector phi, int mode) {

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
    buf.data(), initFlg.data(), stride);

  return std::exp(ll0) + std::exp(ll1);
}


// ---------------------------------------------------------------------------
// cpp_log_likelihood_ecology: full-data orchestrator under the ecology mixture
// ---------------------------------------------------------------------------
//
// Mirrors the R-side .MkpEcologyLogLikelihood (R/likelihood.R): recompute
// ecology marginals + edge weights from the current tree, then iterate
// partitions invoking the ecology-aware pruning variants.
//
// No ascertainment correction yet — corresponds to coding = "none". A
// mixture-aware constant-/singleton-site pseudo-character path is the next
// addition (Phase 3h).

double cpp_log_likelihood_ecology(
    const McmcData& data,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    NumericVector phi,
    IntegerMatrix zMatrix) {

  int nTip = data.nTip;
  int nEdge = parent.size();
  int maxNode = 2 * nTip - 1;
  int kEco = data.ecology.kEcology;
  int mode = data.magnitudeMode;
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
    const PartInfo& part = data.parts[pi];
    int nCharPart = part.tipStates.ncol();

    // Subset zMatrix to this partition's characters (global → local indices).
    IntegerMatrix zPart(nCharPart, kEco);
    for (int c = 0; c < nCharPart; ++c) {
      int gi = part.globalCharIdx[c];
      for (int s = 0; s < kEco; ++s) zPart(c, s) = zMatrix(gi, s);
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
        buf.data(), initFlg.data(), stride);
      if (data.codingType == 1) {  // variable
        for (int c = 0; c < nCharPart; ++c) {
          IntegerVector zVec(kEco);
          for (int s = 0; s < kEco; ++s) zVec[s] = zPart(c, s);
          double pConst = const_site_prob_mkn_eco_single(
            parent, child, neoEl, nTip,
            rateLoss, rates, wEdge, zVec, phi, mode);
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
        buf.data(), initFlg.data(), stride);
      if (data.codingType == 1) {
        for (int c = 0; c < nCharPart; ++c) {
          IntegerVector zVec(kEco);
          for (int s = 0; s < kEco; ++s) zVec[s] = zPart(c, s);
          double pConst = const_site_prob_jc_eco_single(
            parent, child, edgeLen, nTip, kStates,
            rates, wEdge, zVec, phi, mode);
          ll -= std::log(1.0 - pConst);
        }
      }
    } else {
      // Transformational: subgroup by kPrime (per-character)
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
        IntegerMatrix subZ(nSub, kEco);
        for (int c = 0; c < nSub; ++c) {
          for (int t = 0; t < nTip; ++t)
            subStates(t, c) = part.tipStates(t, cols[c]);
          for (int s = 0; s < kEco; ++s)
            subZ(c, s) = zPart(cols[c], s);
        }
        NumericVector rootFreqs(kp, 1.0 / kp);
        int stride = nSub * kp;
        std::vector<double>  buf((maxNode + 1) * stride, 0.0);
        std::vector<uint8_t> initFlg(maxNode + 1, 0u);

        double subLl = pruning_jc_acrv_flat_ecology(
          parent, child, edgeLen, subStates,
          kp, rootFreqs, rates,
          wEdge, subZ, phi, mode,
          buf.data(), initFlg.data(), stride);

        if (data.codingType == 1) {
          for (int c = 0; c < nSub; ++c) {
            IntegerVector zVec(kEco);
            for (int s = 0; s < kEco; ++s) zVec[s] = subZ(c, s);
            double pConst = const_site_prob_jc_eco_single(
              parent, child, edgeLen, nTip, kp,
              rates, wEdge, zVec, phi, mode);
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

    totalLoglik += ll;
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
    IntegerMatrix zMatrix
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  return cpp_log_likelihood_ecology(
    data, parent, child, edgeLen, kPrime,
    rateLoss, rateLogSd, rateNeo, phi, zMatrix);
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
    IntegerVector zRow,        // length kEco
    const NumericMatrix& wEdge // nEdge x kEco
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
  IntegerMatrix zPart(1, kEco);
  for (int s = 0; s < kEco; ++s) zPart(0, s) = zRow[s];

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
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_mkn_eco_single(
        parent, child, neoEl, nTip,
        rateLoss, rates, wEdge, zRow, phi, mode);
      ll -= std::log(1.0 - pConst);
    }
  } else if (part.type == 2) {
    // Known state space
    int kStates = part.k;
    NumericVector rootFreqs(kStates, 1.0 / kStates);
    int stride = 1 * kStates;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology(
      parent, child, edgeLen, tipStates,
      kStates, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_jc_eco_single(
        parent, child, edgeLen, nTip, kStates,
        rates, wEdge, zRow, phi, mode);
      ll -= std::log(1.0 - pConst);
    }
  } else {
    // Transformational: kp drives the state count.
    NumericVector rootFreqs(kp, 1.0 / kp);
    int stride = 1 * kp;
    std::vector<double>  buf((maxNode + 1) * stride, 0.0);
    std::vector<uint8_t> initFlg(maxNode + 1, 0u);
    ll = pruning_jc_acrv_flat_ecology(
      parent, child, edgeLen, tipStates,
      kp, rootFreqs, rates,
      wEdge, zPart, phi, mode,
      buf.data(), initFlg.data(), stride);
    if (data.codingType == 1) {
      double pConst = const_site_prob_jc_eco_single(
        parent, child, edgeLen, nTip, kp,
        rates, wEdge, zRow, phi, mode);
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
    IntegerMatrix zMatrix
) {
  const McmcData& data = *Rcpp::XPtr<McmcData>(dataPtr).get();
  if (!data.ecologyAware) stop("McmcData was not built with ecologyAware = TRUE");
  if (data.codingType == 2) {
    stop("informative coding is not yet supported under ecologyAware");
  }

  int nChar = data.nChar;
  int kEco  = data.ecology.kEcology;
  NumericVector out(nChar);

  NumericMatrix wEdge(parent.size(), kEco);
  recompute_w_edge(data, parent, child, edgeLen, wEdge);

  NumericVector rates = (rateLogSd > 0.0)
    ? cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ)
    : NumericVector(1, 1.0);

  IntegerVector zRow(kEco);
  for (int c = 0; c < nChar; ++c) {
    for (int s = 0; s < kEco; ++s) zRow[s] = zMatrix(c, s);
    int kp = (data.charToPartition[c] >= 0 &&
              data.parts[data.charToPartition[c]].type == 1)
              ? kPrime[c] : 2;
    out[c] = per_char_log_lik_ecology(
      data, c, parent, child, edgeLen,
      kp, rateLoss, rateNeo, rates, phi, zRow, wEdge);
  }
  return out;
}
