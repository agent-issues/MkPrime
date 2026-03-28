#include <Rcpp.h>
#include <cmath>
#include <vector>

// Ascertainment bias correction functions.
//
// Compute P(constant site) and P(singleton site) under JC(k) and MkN models,
// used for the Lewis (2001) variable-coding correction:
//   logL_corrected = logL_raw - nChar * log(1 - P_excluded)
//
// OPP-1: JC matrix-vector product reduced from O(k²) to O(k) per site per
//        edge, exploiting the two-value JC structure:
//        new_cl[i] = p_diff * sum_cl + (p_same - p_diff) * cl[i]
// OPP-2: maxNode = 2*nTip - 1 (covers both rooted and unrooted binary trees;
//        hot path in mcmc_likelihood.cpp uses 2*nTip-2 since trees are always unrooted there)
// OPP-4-lite: single flat std::vector<double> CL replaces vector<vector<double>>
//        (one heap allocation instead of nNode+1 allocations per call)


// ---------------------------------------------------------------------------
// constant_site_prob_jc
//
// P(constant site) for JC(k): sum over k constant patterns.
// Uses k pseudo-characters, one per constant state.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double constant_site_prob_jc(Rcpp::IntegerVector parent,
                             Rcpp::IntegerVector child,
                             Rcpp::NumericVector edge_length,
                             int nTip,
                             int kStates,
                             Rcpp::NumericVector root_freqs,
                             Rcpp::NumericVector rate_multipliers) {
  const int nEdge = parent.size();
  const int nCat  = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;  // OPP-2
  const int root    = nTip + 1;
  const int nChar   = kStates;        // one pseudo-char per constant state
  const int stride  = nChar * kStates;

  const double inv_k = 1.0 / kStates;
  const double km1   = kStates - 1.0;

  // OPP-4-lite: one flat allocation reused across all nCat categories
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  double total_const_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    std::fill(cl_flat.begin(), cl_flat.end(), 0.0);
    std::fill(cl_init.begin(), cl_init.end(), 0u);

    // Tips: pseudo-char s has all tips in state s
    for (int tip = 1; tip <= nTip; ++tip) {
      double* cl = cl_flat.data() + tip * stride;
      for (int s = 0; s < kStates; ++s)
        cl[s * kStates + s] = 1.0;
      cl_init[tip] = 1;
    }

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      const double exp_term = std::exp(-kStates * t / km1);
      const double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff   = inv_k - inv_k * exp_term;
      const double diff_coeff = p_same - p_diff;  // OPP-1
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      if (!cl_init[par]) {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
        cl_init[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] *= p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    for (int s = 0; s < kStates; ++s) {
      const int offset = s * kStates;
      double site_lik = 0.0;
      for (int i = 0; i < kStates; ++i)
        site_lik += root_freqs[i] * clRoot[offset + i];
      total_const_prob += site_lik;
    }
  }

  return total_const_prob / nCat;
}


// ---------------------------------------------------------------------------
// singleton_site_prob_jc
//
// P(singleton site) for JC(k): under JC symmetry, k*(k-1) singleton patterns
// at a given tip all have the same probability, so we need only nTip pseudo-
// characters (one per tip, that tip in state 1, all others in state 0).
// P(singleton) = k*(k-1) * sum_j P(all=0, tip_j=1)
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double singleton_site_prob_jc(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              int kStates,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;  // OPP-2
  const int root    = nTip + 1;
  const int nChar   = nTip;
  const int stride  = nChar * kStates;  // nTip * kStates

  const double inv_k = 1.0 / kStates;
  const double km1   = kStates - 1.0;

  // OPP-4-lite: one flat allocation reused across all nCat categories
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    std::fill(cl_flat.begin(), cl_flat.end(), 0.0);
    std::fill(cl_init.begin(), cl_init.end(), 0u);

    // Tips: pseudo-char j has tip j in state 1, all others in state 0
    for (int tip = 1; tip <= nTip; ++tip) {
      double* cl = cl_flat.data() + tip * stride;
      for (int j = 0; j < nTip; ++j) {
        const int offset = j * kStates;
        if (j == tip - 1)
          cl[offset + 1] = 1.0;  // singleton: state 1
        else
          cl[offset + 0] = 1.0;  // background: state 0
      }
      cl_init[tip] = 1;
    }

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      const double exp_term = std::exp(-kStates * t / km1);
      const double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff   = inv_k - inv_k * exp_term;
      const double diff_coeff = p_same - p_diff;  // OPP-1
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      if (!cl_init[par]) {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
        cl_init[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] *= p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    for (int c = 0; c < nChar; ++c) {
      const int offset = c * kStates;
      double site_lik = 0.0;
      for (int s = 0; s < kStates; ++s)
        site_lik += root_freqs[s] * clRoot[offset + s];
      total_singleton_prob += site_lik;
    }
  }

  total_singleton_prob /= nCat;
  return total_singleton_prob * kStates * (kStates - 1);
}


// ---------------------------------------------------------------------------
// constant_site_prob_mkn
//
// P(constant site) for the MkN 2-state asymmetric model.
// Two pseudo-characters: one for all-0, one for all-1.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double constant_site_prob_mkn(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              double rate_loss,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  const int nEdge    = parent.size();
  const int nCat     = rate_multipliers.size();
  const int maxNode  = 2 * nTip - 1;  // OPP-2
  const int root     = nTip + 1;
  const int kStates  = 2;
  const int nChar    = kStates;
  const int stride   = nChar * kStates;  // 4

  const double sum_rl     = 1.0 + rate_loss;
  const double rate01     = 2.0 / sum_rl;
  const double rate10     = 2.0 * rate_loss / sum_rl;
  const double lambda     = rate01 + rate10;
  const double inv_lam_01 = rate01 / lambda;
  const double inv_lam_10 = rate10 / lambda;

  // OPP-4-lite: single flat allocation
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  double total_const_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    std::fill(cl_flat.begin(), cl_flat.end(), 0.0);
    std::fill(cl_init.begin(), cl_init.end(), 0u);

    // Tips: pseudo-char s has all tips in state s (s=0 or s=1)
    for (int tip = 1; tip <= nTip; ++tip) {
      double* cl = cl_flat.data() + tip * stride;
      for (int s = 0; s < kStates; ++s)
        cl[s * kStates + s] = 1.0;
      cl_init[tip] = 1;
    }

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      const double exp_term = std::exp(-lambda * t);
      const double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      const double P01 = inv_lam_01 - inv_lam_01 * exp_term;
      const double P10 = inv_lam_10 - inv_lam_10 * exp_term;
      const double P11 = inv_lam_01 + inv_lam_10 * exp_term;
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      if (!cl_init[par]) {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          const double cl0 = clCh[offset];
          const double cl1 = clCh[offset + 1];
          clPar[offset]     = P00 * cl0 + P01 * cl1;
          clPar[offset + 1] = P10 * cl0 + P11 * cl1;
        }
        cl_init[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          const double cl0 = clCh[offset];
          const double cl1 = clCh[offset + 1];
          clPar[offset]     *= P00 * cl0 + P01 * cl1;
          clPar[offset + 1] *= P10 * cl0 + P11 * cl1;
        }
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    for (int s = 0; s < kStates; ++s) {
      const int offset = s * kStates;
      total_const_prob += root_freqs[0] * clRoot[offset]
                        + root_freqs[1] * clRoot[offset + 1];
    }
  }

  return total_const_prob / nCat;
}


// ---------------------------------------------------------------------------
// singleton_site_prob_mkn
//
// P(singleton site) for MkN 2-state asymmetric model.
// Unlike JC, states are not symmetric, so we need 2*nTip pseudo-characters:
// first nTip have bg=0/singleton=1; next nTip have bg=1/singleton=0.
// P(singleton) = sum_j [P(all=0, j=1) + P(all=1, j=0)]
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double singleton_site_prob_mkn(Rcpp::IntegerVector parent,
                               Rcpp::IntegerVector child,
                               Rcpp::NumericVector edge_length,
                               int nTip,
                               double rate_loss,
                               Rcpp::NumericVector root_freqs,
                               Rcpp::NumericVector rate_multipliers) {
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;  // OPP-2
  const int root    = nTip + 1;
  const int kStates = 2;
  const int nChar   = 2 * nTip;   // 2n pseudo-characters
  const int stride  = nChar * kStates;

  const double sum_rl     = 1.0 + rate_loss;
  const double rate01     = 2.0 / sum_rl;
  const double rate10     = 2.0 * rate_loss / sum_rl;
  const double lambda     = rate01 + rate10;
  const double inv_lam_01 = rate01 / lambda;
  const double inv_lam_10 = rate10 / lambda;

  // OPP-4-lite: single flat allocation
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    std::fill(cl_flat.begin(), cl_flat.end(), 0.0);
    std::fill(cl_init.begin(), cl_init.end(), 0u);

    // Tips: first nTip pseudo-chars bg=0/single=1; next nTip bg=1/single=0
    for (int tip = 1; tip <= nTip; ++tip) {
      double* cl = cl_flat.data() + tip * stride;
      for (int j = 0; j < nTip; ++j) {
        cl[j * kStates + (j == tip - 1 ? 1 : 0)] = 1.0;           // bg=0, single=1
        cl[(nTip + j) * kStates + (j == tip - 1 ? 0 : 1)] = 1.0;  // bg=1, single=0
      }
      cl_init[tip] = 1;
    }

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      const double exp_term = std::exp(-lambda * t);
      const double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      const double P01 = inv_lam_01 - inv_lam_01 * exp_term;
      const double P10 = inv_lam_10 - inv_lam_10 * exp_term;
      const double P11 = inv_lam_01 + inv_lam_10 * exp_term;
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      if (!cl_init[par]) {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          const double cl0 = clCh[offset];
          const double cl1 = clCh[offset + 1];
          clPar[offset]     = P00 * cl0 + P01 * cl1;
          clPar[offset + 1] = P10 * cl0 + P11 * cl1;
        }
        cl_init[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kStates;
          const double cl0 = clCh[offset];
          const double cl1 = clCh[offset + 1];
          clPar[offset]     *= P00 * cl0 + P01 * cl1;
          clPar[offset + 1] *= P10 * cl0 + P11 * cl1;
        }
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    for (int c = 0; c < nChar; ++c) {
      const int offset = c * kStates;
      total_singleton_prob += root_freqs[0] * clRoot[offset]
                            + root_freqs[1] * clRoot[offset + 1];
    }
  }

  return total_singleton_prob / nCat;
}
