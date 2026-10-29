// Frozen copy of MkPrime's likelihood code as of origin/main de315cf
// (src/likelihood.cpp, src/acrv.cpp, src/mkn_rates.h), before the
// kernels moved to inst/include/MkPrime. test-header-legacy.R compiles it
// with Rcpp::sourceCpp() and requires identical() results. Do not edit.
#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>


// Neomorphic (asymmetric binary) rates for rate_loss r, with stationary
// frequencies pi0 = r / (1 + r), pi1 = 1 / (1 + r). Scaled so the
// pi-weighted mean rate pi0 * rate01 + pi1 * rate10 is 1, as RevBayes'
// fnFreeK(rescaled = TRUE) does; branch lengths are then in expected changes.
inline void mkn_rates(double rateLoss, double& rate01, double& rate10) {
  const double half = 0.5 * (1.0 + rateLoss);
  rate01 = half / rateLoss;  // gain, 0 -> 1
  rate10 = half;             // loss, 1 -> 0
}



// Felsenstein pruning for a single partition under the JC(k) model.
//
// Computes the log-likelihood of a set of characters on a tree, where all
// characters share the same number of states (kStates) and use the JC model.
//
// The tree is given in TreeTools canonical preorder: parent[i] / child[i]
// define edge i with root-to-tip ordering (children sorted by smallest
// descendant).  Traversal iterates edges in reverse so children are
// processed before parents (equivalent to postorder for Felsenstein).
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
double legacy_pruning_jc(Rcpp::IntegerVector parent,
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
  for (int e = nEdge - 1; e >= 0; --e) {
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
  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e]; // 1-indexed
    int ch = child[e];   // 1-indexed
    double t = edge_length[e];

    // Compute JC transition probabilities for this branch
    double inv_k = 1.0 / kStates;
    // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
    double arg = -kStates * t / (kStates - 1.0);
    double neg_expm1 = -std::expm1(arg);
    double exp_term  = 1.0 - neg_expm1;
    double p_same = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff = inv_k * neg_expm1;

    // OPP-1: JC symmetry → O(k) product: new_cl[i] = p_diff*sum + (p_same-p_diff)*cl[i]
    double diff_coeff = p_same - p_diff;
    if (!initialized[par]) {
      // First child: initialize parent CL
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double sum_cl = 0.0;
        for (int j = 0; j < kStates; ++j) sum_cl += CL[ch][offset + j];
        for (int i = 0; i < kStates; ++i)
          CL[par][offset + i] = p_diff * sum_cl + diff_coeff * CL[ch][offset + i];
      }
      initialized[par] = true;
    } else {
      // Subsequent children: multiply into parent CL
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double sum_cl = 0.0;
        for (int j = 0; j < kStates; ++j) sum_cl += CL[ch][offset + j];
        for (int i = 0; i < kStates; ++i)
          CL[par][offset + i] *= p_diff * sum_cl + diff_coeff * CL[ch][offset + i];
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


// Collapsed-state Felsenstein pruning for JC(kFull) when only kObs < kFull
// distinct states are observed across all tips.
//
// Exploits strong lumpability of JC(kFull) under the partition
// {{0}, {1}, ..., {kObs-1}, U}, where U is the set of n_U = kFull - kObs
// states that never appear at any tip. The pruning runs on
// kEff = kObs + 1 conditional-likelihood columns:
//   - columns 0..kObs-1: observed singleton classes
//   - column kObs: the lumped class (carries multiplicity n_U)
//
// Transition probabilities use kFull in the JC analytic formula (so branch
// lengths retain their full-model semantics: substitutions within U still
// consume branch length). The Felsenstein O(k) update uses a weighted
// row-sum that gives the lumped column weight n_U:
//   sum_eff = sum_{j<kObs} CL[j] + n_U * CL[kObs]
//   new_CL[i] = p_diff * sum_eff + (p_same - p_diff) * CL[i]   for all i
// Root frequencies are (1/kFull) on observed columns, (n_U/kFull) on the
// lumped column, so the root site-lik reduces to (1/kFull) * sum_eff_root.
//
// CL convention: average-over-class form
//   cl_lump[I] := (1 / |I|) * sum_{s in I} cl_full[s]
// This gives the uniform Felsenstein update above (no per-row class-size
// factor). Tip initialisation (data is already 0..kObs-1 contiguous per
// .PhyDatToIntMatrix; -1 = missing):
//   observed state s: CL[s] = 1, others 0
//   missing:          CL[0..kObs-1] = 1, CL[kObs] = 1
// (Under the average convention a missing tip is consistent with every
//  class, so cl_lump[I] = 1 for all I; the lumped column's multiplicity
//  enters only via the weighted Σ_eff above, not via the tip value.)
//
// Returns the identical log-likelihood as legacy_pruning_jc(... kStates=kFull ...)
// would on the equivalent uncollapsed tip data.

// [[Rcpp::export]]
double legacy_pruning_jc_collapsed(Rcpp::IntegerVector parent,
                            Rcpp::IntegerVector child,
                            Rcpp::NumericVector edge_length,
                            Rcpp::IntegerMatrix tip_states,
                            int kFull,
                            int kObs) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("pruning_jc_collapsed requires 1 <= kObs < kFull.");
  }
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  int kEff = kObs + 1;
  double n_U = static_cast<double>(kFull - kObs);

  int maxNode = 0;
  for (int e = nEdge - 1; e >= 0; --e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;

  int clSize = nChar * kEff;
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  for (int tip = 1; tip <= nTip; ++tip) {
    for (int c = 0; c < nChar; ++c) {
      int state = tip_states(tip - 1, c);
      int offset = c * kEff;
      if (state < 0) {
        for (int s = 0; s < kEff; ++s) CL[tip][offset + s] = 1.0;
      } else {
        // state is guaranteed to be in 0..kObs-1 by the data-remap invariant
        CL[tip][offset + state] = 1.0;
      }
    }
    initialized[tip] = true;
  }

  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e];
    int ch = child[e];
    double t = edge_length[e];

    double inv_k = 1.0 / kFull;
    // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
    double arg = -kFull * t / (kFull - 1.0);
    double neg_expm1 = -std::expm1(arg);
    double exp_term  = 1.0 - neg_expm1;
    double p_same = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff = inv_k * neg_expm1;
    double diff_coeff = p_same - p_diff;

    if (!initialized[par]) {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kEff;
        double sum_eff = 0.0;
        for (int j = 0; j < kObs; ++j) sum_eff += CL[ch][offset + j];
        sum_eff += n_U * CL[ch][offset + kObs];
        for (int i = 0; i < kEff; ++i)
          CL[par][offset + i] = p_diff * sum_eff + diff_coeff * CL[ch][offset + i];
      }
      initialized[par] = true;
    } else {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kEff;
        double sum_eff = 0.0;
        for (int j = 0; j < kObs; ++j) sum_eff += CL[ch][offset + j];
        sum_eff += n_U * CL[ch][offset + kObs];
        for (int i = 0; i < kEff; ++i)
          CL[par][offset + i] *= p_diff * sum_eff + diff_coeff * CL[ch][offset + i];
      }
    }
  }

  int root = nTip + 1;
  double inv_kFull = 1.0 / kFull;
  double logLik = 0.0;
  for (int c = 0; c < nChar; ++c) {
    int offset = c * kEff;
    double sum_eff = 0.0;
    for (int s = 0; s < kObs; ++s) sum_eff += CL[root][offset + s];
    sum_eff += n_U * CL[root][offset + kObs];
    double site_lik = inv_kFull * sum_eff;
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
double legacy_pruning_mkn(Rcpp::IntegerVector parent,
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
  for (int e = nEdge - 1; e >= 0; --e) {
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
  double rate01, rate10;
  mkn_rates(rate_loss, rate01, rate10);
  double lambda = rate01 + rate10;

  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e];
    int ch = child[e];
    double t = edge_length[e];

    // FAST-EXP-001: expm1 form avoids cancellation in P01/P10 at small lambda*t.
    double arg = -lambda * t;
    double neg_expm1 = -std::expm1(arg);
    double exp_term  = 1.0 - neg_expm1;
    double inv_lam_01 = rate01 / lambda;
    double inv_lam_10 = rate10 / lambda;

    // P matrix:
    // P[0][0] = inv_lam_10 + inv_lam_01 * exp
    // P[0][1] = inv_lam_01 * (1 - exp)
    // P[1][0] = inv_lam_10 * (1 - exp)
    // P[1][1] = inv_lam_01 + inv_lam_10 * exp
    double P00 = inv_lam_10 + inv_lam_01 * exp_term;
    double P01 = inv_lam_01 * neg_expm1;
    double P10 = inv_lam_10 * neg_expm1;
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
double legacy_pruning_jc_acrv(Rcpp::IntegerVector parent,
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

      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      double arg = -kStates * t / km1;
      double neg_expm1 = -std::expm1(arg);
      double exp_term  = 1.0 - neg_expm1;
      double p_same = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff = inv_k * neg_expm1;

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


// Collapsed-state ACRV variant of pruning_jc_acrv. See pruning_jc_collapsed
// in src/likelihood.cpp for the math; structure is identical to
// pruning_jc_acrv but with kEff = kObs + 1 columns, weighted Σ_eff using n_U
// on the lumped column, JC analytic formula in terms of kFull, and root
// frequencies (1/kFull) for observed columns and (n_U/kFull) for the lumped.

// [[Rcpp::export]]
double legacy_pruning_jc_acrv_collapsed(Rcpp::IntegerVector parent,
                                 Rcpp::IntegerVector child,
                                 Rcpp::NumericVector edge_length,
                                 Rcpp::IntegerMatrix tip_states,
                                 int kFull,
                                 int kObs,
                                 Rcpp::NumericVector rate_multipliers) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("pruning_jc_acrv_collapsed requires 1 <= kObs < kFull.");
  }
  int nEdge = parent.size();
  int nTip = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat = rate_multipliers.size();
  int kEff = kObs + 1;
  double n_U = static_cast<double>(kFull - kObs);

  int maxNode = 0;
  for (int e = nEdge - 1; e >= 0; --e) {
    if (parent[e] > maxNode) maxNode = parent[e];
    if (child[e] > maxNode) maxNode = child[e];
  }
  int nNode = maxNode;
  int root = nTip + 1;

  int clSize = nChar * kEff;
  std::vector<double> site_lik_sum(nChar, 0.0);
  std::vector<std::vector<double>> CL(nNode + 1, std::vector<double>(clSize, 0.0));
  std::vector<bool> initialized(nNode + 1, false);

  double inv_k = 1.0 / kFull;
  double km1 = kFull - 1.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = 0; n <= nNode; ++n) {
      std::fill(CL[n].begin(), CL[n].end(), 0.0);
      initialized[n] = false;
    }

    for (int tip = 1; tip <= nTip; ++tip) {
      for (int c = 0; c < nChar; ++c) {
        int state = tip_states(tip - 1, c);
        int offset = c * kEff;
        if (state < 0) {
          for (int s = 0; s < kEff; ++s) CL[tip][offset + s] = 1.0;
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

      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      double arg = -kFull * t / km1;
      double neg_expm1 = -std::expm1(arg);
      double exp_term  = 1.0 - neg_expm1;
      double p_same = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff = inv_k * neg_expm1;
      double diff_coeff = p_same - p_diff;

      if (!initialized[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kEff;
          double sum_eff = 0.0;
          for (int j = 0; j < kObs; ++j) sum_eff += CL[ch][offset + j];
          sum_eff += n_U * CL[ch][offset + kObs];
          for (int i = 0; i < kEff; ++i)
            CL[par][offset + i] = p_diff * sum_eff + diff_coeff * CL[ch][offset + i];
        }
        initialized[par] = true;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kEff;
          double sum_eff = 0.0;
          for (int j = 0; j < kObs; ++j) sum_eff += CL[ch][offset + j];
          sum_eff += n_U * CL[ch][offset + kObs];
          for (int i = 0; i < kEff; ++i)
            CL[par][offset + i] *= p_diff * sum_eff + diff_coeff * CL[ch][offset + i];
        }
      }
    }

    for (int c = 0; c < nChar; ++c) {
      int offset = c * kEff;
      double sum_eff = 0.0;
      for (int s = 0; s < kObs; ++s) sum_eff += CL[root][offset + s];
      sum_eff += n_U * CL[root][offset + kObs];
      site_lik_sum[c] += inv_k * sum_eff;
    }
  }

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
double legacy_pruning_mkn_acrv(Rcpp::IntegerVector parent,
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

  double rate01, rate10;
  mkn_rates(rate_loss, rate01, rate10);
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
      // FAST-EXP-001: expm1 form avoids cancellation in P01/P10 at small lambda*t.
      double arg = -lambda * t;
      double neg_expm1 = -std::expm1(arg);
      double exp_term  = 1.0 - neg_expm1;
      double inv_lam_01 = rate01 / lambda;
      double inv_lam_10 = rate10 / lambda;

      double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      double P01 = inv_lam_01 * neg_expm1;
      double P10 = inv_lam_10 * neg_expm1;
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

// Kernel expressions exactly as written at the call sites that now use
// MkPrime::JCRates / MkNRates. Form A: likelihood.cpp, acrv.cpp,
// ascertainment.cpp, mcmc_likelihood.cpp, rate_matrix.cpp. Form B
// (diff_coeff = exp_term): gibbs_partial_cl.h jc_trans_params().
// [[Rcpp::export]]
Rcpp::NumericMatrix legacy_jc_probs(int kStates, Rcpp::NumericVector tt) {
  Rcpp::NumericMatrix out(tt.size(), 4);
  double inv_k = 1.0 / kStates;
  double km1 = kStates - 1.0;
  for (R_xlen_t i = 0; i < tt.size(); ++i) {
    double t = tt[i];
    double arg = -kStates * t / km1;
    double neg_expm1 = -std::expm1(arg);
    double exp_term  = 1.0 - neg_expm1;
    double p_same = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff = inv_k * neg_expm1;
    double diff_coeff = p_same - p_diff;
    out(i, 0) = p_same; out(i, 1) = p_diff; out(i, 2) = diff_coeff; out(i, 3) = exp_term;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix legacy_mkn_probs(double rate_loss, Rcpp::NumericVector tt) {
  Rcpp::NumericMatrix out(tt.size(), 4);
  double rate01, rate10;
  mkn_rates(rate_loss, rate01, rate10);
  double lambda = rate01 + rate10;
  double inv_lam_01 = rate01 / lambda;
  double inv_lam_10 = rate10 / lambda;
  for (R_xlen_t i = 0; i < tt.size(); ++i) {
    double t = tt[i];
    double arg = -lambda * t;
    double neg_expm1 = -std::expm1(arg);
    double exp_term  = 1.0 - neg_expm1;
    out(i, 0) = inv_lam_10 + inv_lam_01 * exp_term;
    out(i, 1) = inv_lam_01 * neg_expm1;
    out(i, 2) = inv_lam_10 * neg_expm1;
    out(i, 3) = inv_lam_01 + inv_lam_10 * exp_term;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector legacy_acrv_rates(double rateLogSd, int nCat, std::vector<double> acrvZ) {
  if (rateLogSd <= 0.0) return Rcpp::NumericVector(nCat, 1.0);
  double mu = -rateLogSd * rateLogSd / 2.0;
  Rcpp::NumericVector rates(nCat);
  double total = 0.0;
  for (int i = 0; i < nCat; ++i) {
    rates[i] = std::exp(mu + rateLogSd * acrvZ[i]);
    total += rates[i];
  }
  for (int i = 0; i < nCat; ++i) rates[i] *= nCat / total;
  return rates;
}
