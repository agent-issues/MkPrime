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
