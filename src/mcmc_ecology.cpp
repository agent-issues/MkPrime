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

using namespace Rcpp;


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
