#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>

// Felsenstein pruning with ACRV (Among-Character Rate Variation).
//
// For each rate category, performs a full tree traversal with scaled branch
// lengths. Per-site likelihoods are averaged across categories, then logged.
//
// This is the JC(k) version. The MkN version follows the same pattern.
//
// Parameters:
//   parent, child, edge_length: tree in canonical preorder (1-indexed)
//   tip_states: nTip x nChar (0-indexed, -1 for missing)
//   kStates: number of states
//   root_freqs: equilibrium frequencies (length kStates)
//   rate_multipliers: numeric vector of rate multipliers (length nCat)
//
// Returns: total log-likelihood summed over characters, accounting for ACRV.

// [[Rcpp::export]]
double pruning_jc_acrv(Rcpp::IntegerVector parent,
                       Rcpp::IntegerVector child,
                       Rcpp::NumericVector edge_length,
                       Rcpp::IntegerMatrix tip_states,
                       int kStates,
                       Rcpp::NumericVector root_freqs,
                       Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat = rate_multipliers.size();

  int maxNode = 0;
  for (int e = nEdge - 1; e >= 0; --e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  int clSize = nChar * kStates;

  // Accumulator for per-site likelihoods averaged across categories
  std::vector<double> site_lik_sum(nChar, 0.0);

  // Reusable CL storage for one rate category
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double inv_k = 1.0 / kStates;
  double km1 = kStates - 1.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    // Reset CL and initialized flags
    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    // Initialize tips
    for (int tip = 1; tip <= nTip; ++tip) {
      for (int c = 0; c < nChar; ++c) {
        int state = tip_states(tip - 1, c);
        int offset = c * kStates;
        if (state < 0) {
          for (int s = 0; s < kStates; ++s) {
            CL[tip][offset + s] = 1.0;
          }
        } else {
          CL[tip][offset + state] = 1.0;
        }
      }
      initialized[tip] = true;
    }

    // Post-order traversal with scaled branch lengths
    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch = child[e];
      double t = edge_length[e] * rate;

      double exp_term = std::exp(-kStates * t / km1);
      double p_same = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff = inv_k - inv_k * exp_term;

      // OPP-1: JC symmetry → O(k) product: new_cl[i] = p_diff*sum + (p_same-p_diff)*cl[i]
      double diff_coeff = p_same - p_diff;
      if (!initialized[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += CL[ch][offset + j];
          for (int i = 0; i < kStates; ++i)
            CL[par][offset + i] = p_diff * sum_cl + diff_coeff * CL[ch][offset + i];
        }
        initialized[par] = true;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += CL[ch][offset + j];
          for (int i = 0; i < kStates; ++i)
            CL[par][offset + i] *= p_diff * sum_cl + diff_coeff * CL[ch][offset + i];
        }
      }
    }

    // Add per-site likelihoods for this category
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double site_lik = 0.0;
      for (int s = 0; s < kStates; ++s) {
        site_lik += root_freqs[s] * CL[root][offset + s];
      }
      site_lik_sum[c] += site_lik;
    }
  }

  // Compute total log-likelihood: sum_c log(site_lik_sum[c] / nCat)
  double logLik = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg_lik = site_lik_sum[c] * inv_nCat;
    if (avg_lik <= 0.0) {
      return R_NegInf;
    }
    logLik += std::log(avg_lik);
  }

  return logLik;
}


// MkN version of ACRV-integrated pruning.

// [[Rcpp::export]]
double pruning_mkn_acrv(Rcpp::IntegerVector parent,
                        Rcpp::IntegerVector child,
                        Rcpp::NumericVector edge_length,
                        Rcpp::IntegerMatrix tip_states,
                        double rate_loss,
                        Rcpp::NumericVector root_freqs,
                        Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat = rate_multipliers.size();
  const int kStates = 2;

  int maxNode = 0;
  for (int e = nEdge - 1; e >= 0; --e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  int clSize = nChar * kStates;
  std::vector<double> site_lik_sum(nChar, 0.0);
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double sum_rl = 1.0 + rate_loss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rate_loss / sum_rl;
  double lambda = rate01 + rate10;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

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

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch = child[e];
      double t = edge_length[e] * rate;
      double exp_term = std::exp(-lambda * t);
      double inv_lam_01 = rate01 / lambda;
      double inv_lam_10 = rate10 / lambda;

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

    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double site_lik = root_freqs[0] * CL[root][offset + 0] +
                        root_freqs[1] * CL[root][offset + 1];
      site_lik_sum[c] += site_lik;
    }
  }

  double logLik = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg_lik = site_lik_sum[c] * inv_nCat;
    if (avg_lik <= 0.0) return R_NegInf;
    logLik += std::log(avg_lik);
  }

  return logLik;
}
