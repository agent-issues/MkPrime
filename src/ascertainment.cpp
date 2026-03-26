#include <Rcpp.h>
#include <cmath>
#include <vector>

// Compute the probability of a constant (invariant) site under JC(k).
//
// This is needed for the ascertainment bias correction (Lewis 2001):
//   logL_corrected = logL_uncorrected - nChar * log(1 - P_constant)
//
// P_constant = sum_{s=0}^{k-1} π_s * L(all_tips = s | tree, model)
//
// With ACRV:
//   P_constant = (1/nCat) * sum_cat sum_s π_s * L(all_tips=s | tree, rate_cat)
//
// Implementation: creates k "pseudo-characters" where all tips have state s,
// then runs a pruning pass for each rate category. Returns P_constant.
//
// Parameters:
//   parent, child, edge_length: tree in postorder (1-indexed)
//   nTip: number of tips
//   kStates: number of states
//   root_freqs: equilibrium frequencies
//   rate_multipliers: ACRV rate multipliers (length nCat)
//
// Returns: P(constant site), a scalar in [0, 1].

// [[Rcpp::export]]
double constant_site_prob_jc(Rcpp::IntegerVector parent,
                             Rcpp::IntegerVector child,
                             Rcpp::NumericVector edge_length,
                             int nTip,
                             int kStates,
                             Rcpp::NumericVector root_freqs,
                             Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nCat = rate_multipliers.size();

  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  // We have kStates "pseudo-characters" (one per constant state pattern)
  int nChar = kStates;
  int clSize = nChar * kStates;

  double inv_k = 1.0 / kStates;
  double km1 = kStates - 1.0;

  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double total_const_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    // Initialize tips: pseudo-char s has all tips in state s
    for (int tip = 1; tip <= nTip; ++tip) {
      for (int s = 0; s < kStates; ++s) {
        int offset = s * kStates;
        // This tip is in state s for pseudo-char s
        CL[tip][offset + s] = 1.0;
      }
      initialized[tip] = true;
    }

    for (int e = 0; e < nEdge; ++e) {
      int par = parent[e];
      int ch = child[e];
      double t = edge_length[e] * rate;

      double exp_term = std::exp(-kStates * t / km1);
      double p_same = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff = inv_k - inv_k * exp_term;

      if (!initialized[par]) {
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

    // Sum site likelihoods for all constant patterns
    for (int s = 0; s < kStates; ++s) {
      int offset = s * kStates;
      double site_lik = 0.0;
      for (int i = 0; i < kStates; ++i) {
        site_lik += root_freqs[i] * CL[root][offset + i];
      }
      total_const_prob += site_lik;
    }
  }

  // Average across rate categories
  total_const_prob /= nCat;

  return total_const_prob;
}


// Compute the probability of a singleton (autapomorphic) site under JC(k).
//
// A singleton pattern has exactly one tip in a unique state (parsimony score 1).
// Under JC(k) symmetry with uniform root frequencies, all k*(k-1) singleton
// patterns at a given tip have the same probability, so we only need n
// pseudo-characters (one per tip, with that tip in state 1 and all others in
// state 0).
//
// P(singleton) = k*(k-1) * sum_j P(all_tips=0, tip_j=1)
//
// With ACRV, averaged across rate categories.

// [[Rcpp::export]]
double singleton_site_prob_jc(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              int kStates,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nCat = rate_multipliers.size();

  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  // n pseudo-characters: pseudo-char j has tip j = state 1, all others = state 0
  int nChar = nTip;
  int clSize = nChar * kStates;

  double inv_k = 1.0 / kStates;
  double km1 = kStates - 1.0;

  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    // Initialize tip CLs
    for (int tip = 1; tip <= nTip; ++tip) {
      for (int j = 0; j < nTip; ++j) {
        int offset = j * kStates;
        if (j == tip - 1) {
          // This tip is the singleton: state 1
          CL[tip][offset + 1] = 1.0;
        } else {
          // Background: state 0
          CL[tip][offset + 0] = 1.0;
        }
      }
      initialized[tip] = true;
    }

    // Postorder traversal
    for (int e = 0; e < nEdge; ++e) {
      int par = parent[e];
      int ch = child[e];
      double t = edge_length[e] * rate;

      double exp_term = std::exp(-kStates * t / km1);
      double p_same = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff = inv_k - inv_k * exp_term;

      if (!initialized[par]) {
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

    // Sum site likelihoods for all singleton pseudo-characters
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double site_lik = 0.0;
      for (int s = 0; s < kStates; ++s) {
        site_lik += root_freqs[s] * CL[root][offset + s];
      }
      total_singleton_prob += site_lik;
    }
  }

  // Average across rate categories, multiply by k*(k-1) for JC symmetry
  total_singleton_prob /= nCat;
  total_singleton_prob *= kStates * (kStates - 1);

  return total_singleton_prob;
}


// MkN version: constant site probability for 2-state asymmetric model.

// [[Rcpp::export]]
double constant_site_prob_mkn(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              double rate_loss,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nCat = rate_multipliers.size();
  const int kStates = 2;

  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  int nChar = kStates; // 2 pseudo-characters
  int clSize = nChar * kStates;

  double sum_rl = 1.0 + rate_loss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rate_loss / sum_rl;
  double lambda = rate01 + rate10;

  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double total_const_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    for (int tip = 1; tip <= nTip; ++tip) {
      for (int s = 0; s < kStates; ++s) {
        CL[tip][s * kStates + s] = 1.0;
      }
      initialized[tip] = true;
    }

    for (int e = 0; e < nEdge; ++e) {
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

    for (int s = 0; s < kStates; ++s) {
      int offset = s * kStates;
      double site_lik = root_freqs[0] * CL[root][offset + 0] +
                        root_freqs[1] * CL[root][offset + 1];
      total_const_prob += site_lik;
    }
  }

  total_const_prob /= nCat;
  return total_const_prob;
}


// MkN singleton site probability for 2-state asymmetric model.
//
// Unlike JC, the two states are NOT symmetric under MkN, so we need
// 2n pseudo-characters: for each tip j, one with bg=0/single=1 and
// one with bg=1/single=0.
//
// P(singleton) = sum_j [P(all=0, j=1) + P(all=1, j=0)]

// [[Rcpp::export]]
double singleton_site_prob_mkn(Rcpp::IntegerVector parent,
                               Rcpp::IntegerVector child,
                               Rcpp::NumericVector edge_length,
                               int nTip,
                               double rate_loss,
                               Rcpp::NumericVector root_freqs,
                               Rcpp::NumericVector rate_multipliers) {
  int nEdge = parent.size();
  int nCat = rate_multipliers.size();
  const int kStates = 2;

  int maxNode = 0;
  for (int e = 0; e < nEdge; ++e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  // 2n pseudo-characters: first n have bg=0/single=1, next n have bg=1/single=0
  int nChar = 2 * nTip;
  int clSize = nChar * kStates;

  double sum_rl = 1.0 + rate_loss;
  double rate01 = 2.0 / sum_rl;
  double rate10 = 2.0 * rate_loss / sum_rl;
  double lambda = rate01 + rate10;

  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    // Initialize tip CLs
    for (int tip = 1; tip <= nTip; ++tip) {
      // First n pseudo-chars: bg=0, singleton at tip j has state 1
      for (int j = 0; j < nTip; ++j) {
        int offset = j * kStates;
        if (j == tip - 1) {
          CL[tip][offset + 1] = 1.0; // singleton: state 1
        } else {
          CL[tip][offset + 0] = 1.0; // background: state 0
        }
      }
      // Next n pseudo-chars: bg=1, singleton at tip j has state 0
      for (int j = 0; j < nTip; ++j) {
        int offset = (nTip + j) * kStates;
        if (j == tip - 1) {
          CL[tip][offset + 0] = 1.0; // singleton: state 0
        } else {
          CL[tip][offset + 1] = 1.0; // background: state 1
        }
      }
      initialized[tip] = true;
    }

    for (int e = 0; e < nEdge; ++e) {
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

    // Sum site likelihoods for all singleton pseudo-characters
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double site_lik = root_freqs[0] * CL[root][offset + 0] +
                        root_freqs[1] * CL[root][offset + 1];
      total_singleton_prob += site_lik;
    }
  }

  total_singleton_prob /= nCat;
  return total_singleton_prob;
}
