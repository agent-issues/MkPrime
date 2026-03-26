#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

// Felsenstein pruning for a single partition under the JC(k) model.
//
// Computes the log-likelihood of a set of characters on a tree, where all
// characters share the same number of states (kStates) and use the JC model.
//
// The tree is given in ape's postorder format: parent[i] and child[i] define
// edge i, with edges ordered so children are processed before parents.
//
// Tip states are 0-indexed integers. A value of -1 indicates missing data
// (all states equally likely at that tip for that character).
//
// Parameters:
//   parent: integer vector of parent node indices (1-indexed, R convention)
//   child: integer vector of child node indices (1-indexed)
//   edge_length: branch lengths for each edge
//   tip_states: nTip x nChar matrix (0-indexed states, -1 for missing)
//   kStates: number of states for JC model
//   root_freqs: equilibrium frequencies at root (length kStates)
//
// Returns: total log-likelihood summed over all characters in the partition.

// [[Rcpp::export]]
double pruning_jc(Rcpp::IntegerVector parent,
                  Rcpp::IntegerVector child,
                  Rcpp::NumericVector edge_length,
                  Rcpp::IntegerMatrix tip_states,
                  int kStates,
                  Rcpp::NumericVector root_freqs) {
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nNode = nTip + nEdge / 2 + 1; // approx; use max node index instead
  // Find actual nNode from max parent index
  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  nNode = maxNode; // nodes are 1-indexed

  // Conditional likelihoods: CL[node][char * kStates + state]
  // Node indices are 1-based (matching R), so allocate nNode+1
  int clSize = nChar * kStates;
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  // Track whether a node's CL has been initialized (for the product)
  std::vector<bool> initialized(nNode + 1, false);

  // Initialize tip CLs
  for (int tip = 1; tip <= nTip; ++tip) {
    for (int c = 0; c < nChar; ++c) {
      int state = tip_states(tip - 1, c); // R matrix is 0-indexed in C++
      int offset = c * kStates;
      if (state < 0) {
        // Missing data: all states equally likely
        for (int s = 0; s < kStates; ++s) {
          CL[tip][offset + s] = 1.0;
        }
      } else {
        // Observed state
        CL[tip][offset + state] = 1.0;
      }
    }
    initialized[tip] = true;
  }

  // Post-order traversal: process edges in the given order
  // For each edge (parent -> child), compute the contribution of this child
  // to the parent's conditional likelihood
  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e]; // 1-indexed
    int ch = child[e];   // 1-indexed
    double t = edge_length[e];

    // Compute JC transition probabilities for this branch
    double inv_k = 1.0 / kStates;
    double exp_term = std::exp(-kStates * t / (kStates - 1.0));
    double p_same = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff = inv_k - inv_k * exp_term;

    if (!initialized[par]) {
      // First child: initialize parent CL
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        for (int i = 0; i < kStates; ++i) {
          double sum = 0.0;
          for (int j = 0; j < kStates; ++j) {
            double pij = (i == j) ? p_same : p_diff;
            sum += pij * CL[ch][offset + j];
          }
          CL[par][offset + i] = sum;
        }
      }
      initialized[par] = true;
    } else {
      // Subsequent children: multiply into parent CL
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        for (int i = 0; i < kStates; ++i) {
          double sum = 0.0;
          for (int j = 0; j < kStates; ++j) {
            double pij = (i == j) ? p_same : p_diff;
            sum += pij * CL[ch][offset + j];
          }
          CL[par][offset + i] *= sum;
        }
      }
    }
  }

  // Find root (the node that is a parent but never a child, or nTip+1)
  int root = nTip + 1; // ape convention

  // Compute log-likelihood at root
  double logLik = 0.0;
  for (int c = 0; c < nChar; ++c) {
    int offset = c * kStates;
    double site_lik = 0.0;
    for (int s = 0; s < kStates; ++s) {
      site_lik += root_freqs[s] * CL[root][offset + s];
    }
    if (site_lik <= 0.0) {
      return R_NegInf;
    }
    logLik += std::log(site_lik);
  }

  return logLik;
}


// Felsenstein pruning for MkN (asymmetric binary) model.
//
// Same structure as pruning_jc but uses asymmetric 2-state P(t).

// [[Rcpp::export]]
double pruning_mkn(Rcpp::IntegerVector parent,
                   Rcpp::IntegerVector child,
                   Rcpp::NumericVector edge_length,
                   Rcpp::IntegerMatrix tip_states,
                   double rate_loss,
                   Rcpp::NumericVector root_freqs) {
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  const int kStates = 2;

  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;

  int clSize = nChar * kStates;
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  // Initialize tip CLs
  for (int tip = 1; tip <= nTip; ++tip) {
    for (int c = 0; c < nChar; ++c) {
      int state = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        CL[tip][offset + 0] = 1.0;
        CL[tip][offset + 1] = 1.0;
      } else {
        CL[tip][offset + state] = 1.0;
      }
    }
    initialized[tip] = true;
  }

  // Precompute MkN rate parameters
  double sum_rl = 1.0 + rate_loss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rate_loss / sum_rl;
  double lambda = rate01 + rate10; // = 2.0

  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    int ch = child[e];
    double t = edge_length[e];

    double exp_term = std::exp(-lambda * t);
    double inv_lam_01 = rate01 / lambda;
    double inv_lam_10 = rate10 / lambda;

    // P matrix:
    // P[0][0] = inv_lam_10 + inv_lam_01 * exp
    // P[0][1] = inv_lam_01 - inv_lam_01 * exp
    // P[1][0] = inv_lam_10 - inv_lam_10 * exp
    // P[1][1] = inv_lam_01 + inv_lam_10 * exp
    double P00 = inv_lam_10 + inv_lam_01 * exp_term;
    double P01 = inv_lam_01 - inv_lam_01 * exp_term;
    double P10 = inv_lam_10 - inv_lam_10 * exp_term;
    double P11 = inv_lam_01 + inv_lam_10 * exp_term;

    if (!initialized[par]) {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double cl0 = CL[ch][offset + 0];
        double cl1 = CL[ch][offset + 1];
        CL[par][offset + 0] = P00 * cl0 + P01 * cl1;
        CL[par][offset + 1] = P10 * cl0 + P11 * cl1;
      }
      initialized[par] = true;
    } else {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double cl0 = CL[ch][offset + 0];
        double cl1 = CL[ch][offset + 1];
        CL[par][offset + 0] *= P00 * cl0 + P01 * cl1;
        CL[par][offset + 1] *= P10 * cl0 + P11 * cl1;
      }
    }
  }

  int root = nTip + 1;
  double logLik = 0.0;
  for (int c = 0; c < nChar; ++c) {
    int offset = c * kStates;
    double site_lik = root_freqs[0] * CL[root][offset + 0] +
                      root_freqs[1] * CL[root][offset + 1];
    if (site_lik <= 0.0) {
      return R_NegInf;
    }
    logLik += std::log(site_lik);
  }

  return logLik;
}
