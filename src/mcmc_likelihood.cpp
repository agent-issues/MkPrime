// C++ likelihood orchestration for the MCMC hot path.
//
// prepare_mcmc_data() converts the R MkPrimeData and MkPrimeModel objects
// into a persistent C++ McmcData struct (stored via XPtr).
//
// cpp_log_likelihood() replicates the R .MkpLogLikelihood() orchestration
// entirely in C++, calling the existing C++ pruning functions directly.
// This eliminates the per-iteration R overhead of partition looping and
// Rcpp list construction.
//
// M-065: cpp_log_likelihood now accepts parent/child vectors directly,
// avoiding edge matrix construction/decomposition round-trips.

#include "mcmc_state.h"
#include "fast_exp.h"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <set>

using namespace Rcpp;


// ---------------------------------------------------------------------------
// Forward declarations: pruning functions defined in other TUs
// ---------------------------------------------------------------------------

double pruning_jc(IntegerVector parent, IntegerVector child,
                  NumericVector edge_length, IntegerMatrix tip_states,
                  int kStates, NumericVector root_freqs);

double pruning_jc_acrv(IntegerVector parent, IntegerVector child,
                       NumericVector edge_length, IntegerMatrix tip_states,
                       int kStates, NumericVector root_freqs,
                       NumericVector rate_multipliers);

double pruning_mkn(IntegerVector parent, IntegerVector child,
                   NumericVector edge_length, IntegerMatrix tip_states,
                   double rate_loss, NumericVector root_freqs);

double pruning_mkn_acrv(IntegerVector parent, IntegerVector child,
                        NumericVector edge_length, IntegerMatrix tip_states,
                        double rate_loss, NumericVector root_freqs,
                        NumericVector rate_multipliers);

double constant_site_prob_jc(IntegerVector parent, IntegerVector child,
                             NumericVector edge_length, int nTip,
                             int kStates, NumericVector root_freqs,
                             NumericVector rate_multipliers);

double singleton_site_prob_jc(IntegerVector parent, IntegerVector child,
                              NumericVector edge_length, int nTip,
                              int kStates, NumericVector root_freqs,
                              NumericVector rate_multipliers);

double constant_site_prob_mkn(IntegerVector parent, IntegerVector child,
                              NumericVector edge_length, int nTip,
                              double rate_loss, NumericVector root_freqs,
                              NumericVector rate_multipliers);

double singleton_site_prob_mkn(IntegerVector parent, IntegerVector child,
                               NumericVector edge_length, int nTip,
                               double rate_loss, NumericVector root_freqs,
                               NumericVector rate_multipliers);

double mk_prime_relabel_log(int kPrime, int kObs);


// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// OPP-3: acrvZ = precomputed qnorm((i+0.5)/nCat) stored in McmcData — no
// transcendental calls per iteration; only exp() + scaling remain.
NumericVector cpp_acrv_rates(double rateLogSd, int nCat,
                             const std::vector<double>& acrvZ) {
  if (rateLogSd <= 0.0) return NumericVector(nCat, 1.0);
  double mu = -rateLogSd * rateLogSd / 2.0;
  NumericVector rates(nCat);
  double total = 0.0;
  for (int i = 0; i < nCat; ++i) {
    rates[i] = std::exp(mu + rateLogSd * acrvZ[i]);
    total += rates[i];
  }
  for (int i = 0; i < nCat; ++i) rates[i] *= nCat / total;
  return rates;
}

static NumericVector mkn_stationary(double rateLoss) {
  NumericVector f(2);
  f[0] = rateLoss / (1.0 + rateLoss);   // π₀ = rate_loss / (1 + rate_loss)
  f[1] = 1.0 / (1.0 + rateLoss);        // π₁ = 1 / (1 + rate_loss)
  return f;
}


// ---------------------------------------------------------------------------
// Flat-buffer pruning helpers (M-063)
//
// Drop-in alternatives to pruning_jc / pruning_jc_acrv / pruning_mkn /
// pruning_mkn_acrv that use a caller-supplied workspace buffer instead of
// allocating std::vector<std::vector<double>> CL on the heap each call.
//
// buf layout: buf[node * stride + c * kStates + s]  (node 1-indexed)
// initFlg:    uint8_t[nNodeMax+1], reset to 0 for each traversal.
//
// TRAVERSAL DIRECTION: edges must be in any valid preorder (root → tips).
// Iterating in reverse (nEdge-1 → 0) gives a valid bottom-up Felsenstein
// pass: every child CL is fully computed before its parent edge is visited.
// Canonical preorder (children sorted by smallest descendant) is the
// default, but in-place NNI (OPP-6b) may produce non-canonical preorder
// — still valid because NNI only modifies parent assignments when it is
// safe to do so (wRow > edgeRow guard in mcmc.cpp case 5).
// ---------------------------------------------------------------------------

#ifndef NDEBUG
// Debug-only check: verify edges are in valid preorder.
static void assert_valid_preorder(const IntegerVector& parent,
                                  const IntegerVector& child,
                                  int nTip) {
  const int nEdge = parent.size();
  const int root = nTip + 1;
  const int maxNode = 2 * nTip;
  std::vector<bool> seen(maxNode + 1, false);
  seen[root] = true;
  for (int i = 0; i < nEdge; ++i) {
    if (!seen[parent[i]])
      Rcpp::stop("pruning: parent %d not introduced at edge %d "
                 "(preorder invariant violated)", (int)parent[i], i);
    seen[child[i]] = true;
  }
}
#endif

static double pruning_jc_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, NumericVector root_freqs,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();

#ifndef NDEBUG
  assert_valid_preorder(parent, child, nTip);
#endif

  int maxNode = 2 * nTip - 1;  // OPP-2: +1 covers rooted trees (root = 2*nTip-1)
  int clCols = nChar * kStates;

  for (int n = 0; n <= maxNode; ++n) {
    std::fill(buf + n * stride, buf + n * stride + clCols, 0.0);
    initFlg[n] = 0;
  }
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  double inv_k = 1.0 / kStates;
  double km1   = kStates - 1.0;

  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e];
    int ch  = child[e];
    double t        = edge_length[e];
    double exp_term = MKP_EXP(-kStates * t / km1);
    double p_same   = inv_k + (1.0 - inv_k) * exp_term;
    double p_diff   = inv_k - inv_k * exp_term;
    double* clPar   = buf + par * stride;
    double* clCh    = buf + ch  * stride;

    // OPP-1: JC symmetry → O(k) product: new_cl[i] = p_diff*sum + (p_same-p_diff)*cl[i]
    double diff_coeff = p_same - p_diff;
    if (!initFlg[par]) {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double sum_cl = 0.0;
        for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
        for (int i = 0; i < kStates; ++i)
          clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
      }
      initFlg[par] = 1;
    } else {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double sum_cl = 0.0;
        for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
        for (int i = 0; i < kStates; ++i)
          clPar[offset + i] *= p_diff * sum_cl + diff_coeff * clCh[offset + i];
      }
    }
  }

  int root = nTip + 1;
  double* clRoot = buf + root * stride;
  double logLik  = 0.0;
  for (int c = 0; c < nChar; ++c) {
    int offset = c * kStates;
    double sl = 0.0;
    for (int s = 0; s < kStates; ++s)
      sl += root_freqs[s] * clRoot[offset + s];
    if (sl <= 0.0) return R_NegInf;
    logLik += std::log(sl);
  }
  return logLik;
}


static double pruning_jc_acrv_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, NumericVector root_freqs,
    NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();

#ifndef NDEBUG
  assert_valid_preorder(parent, child, nTip);
#endif

  int maxNode = 2 * nTip - 1;  // OPP-2
  int root   = nTip + 1;
  int clCols = nChar * kStates;

  std::vector<double> site_lik_sum(nChar, 0.0);
  double inv_k = 1.0 / kStates;
  double km1   = kStates - 1.0;

  // OPP-CL: tips are constant across rate categories — init once, not per cat.
  // Zero only tip regions, set observed states; internal nodes overwritten by
  // the first-child '=' branch so their stale values are harmless.
  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    // Reset only internal-node init flags (tips stay initialised).
    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch  = child[e];
      double t        = edge_length[e] * rate;
      double exp_term = MKP_EXP(-kStates * t / km1);
      double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff   = inv_k - inv_k * exp_term;
      double* clPar   = buf + par * stride;
      double* clCh    = buf + ch  * stride;

      // OPP-1: JC symmetry → O(k) product
      double diff_coeff = p_same - p_diff;
      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
        initFlg[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] *= p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
      }
    }

    double* clRoot = buf + root * stride;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double sl = 0.0;
      for (int s = 0; s < kStates; ++s)
        sl += root_freqs[s] * clRoot[offset + s];
      site_lik_sum[c] += sl;
    }
  }

  double logLik   = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = site_lik_sum[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// Per-site variant of pruning_jc_acrv_flat: fills siteLL[0..nChar-1] with
// per-character log(avg_lik) instead of returning the sum.
// Used by Gibbs kPrime sweep (M-155) to precompute likelihoods for all
// characters at each candidate k in a single batched traversal.
void pruning_jc_acrv_persite(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates,
    NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride,
    double* siteLL) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();

  int maxNode = 2 * nTip - 1;
  int root   = nTip + 1;
  int clCols = nChar * kStates;

  std::vector<double> site_lik_sum(nChar, 0.0);
  double inv_k = 1.0 / kStates;
  double km1   = kStates - 1.0;

  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];
    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch  = child[e];
      double t        = edge_length[e] * rate;
      double exp_term = std::exp(-kStates * t / km1);
      double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff   = inv_k - inv_k * exp_term;
      double* clPar   = buf + par * stride;
      double* clCh    = buf + ch  * stride;
      double diff_coeff = p_same - p_diff;

      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] = p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
        initFlg[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sum_cl = 0.0;
          for (int j = 0; j < kStates; ++j) sum_cl += clCh[offset + j];
          for (int i = 0; i < kStates; ++i)
            clPar[offset + i] *= p_diff * sum_cl + diff_coeff * clCh[offset + i];
        }
      }
    }

    double* clRoot = buf + root * stride;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double sl = 0.0;
      for (int s = 0; s < kStates; ++s)
        sl += inv_k * clRoot[offset + s];
      site_lik_sum[c] += sl;
    }
  }

  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = site_lik_sum[c] * inv_nCat;
    siteLL[c] = (avg > 0.0) ? std::log(avg) : R_NegInf;
  }
}


static double pruning_mkn_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    double rate_loss, NumericVector root_freqs,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  const int kStates = 2;

#ifndef NDEBUG
  assert_valid_preorder(parent, child, nTip);
#endif

  int maxNode = 2 * nTip - 1;  // OPP-2: +1 covers rooted trees (root = 2*nTip-1)
  int clCols = nChar * kStates;

  for (int n = 0; n <= maxNode; ++n) {
    std::fill(buf + n * stride, buf + n * stride + clCols, 0.0);
    initFlg[n] = 0;
  }
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        cl[offset + 0] = 1.0;
        cl[offset + 1] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  double sum_rl     = 1.0 + rate_loss;
  double rate01     = 2.0 / sum_rl;
  double rate10     = 2.0 * rate_loss / sum_rl;
  double lambda     = rate01 + rate10;
  double inv_lam_01 = rate01 / lambda;
  double inv_lam_10 = rate10 / lambda;

  for (int e = nEdge - 1; e >= 0; --e) {
    int par = parent[e];
    int ch  = child[e];
    double t        = edge_length[e];
    double exp_term = MKP_EXP(-lambda * t);
    double P00 = inv_lam_10 + inv_lam_01 * exp_term;
    double P01 = inv_lam_01 - inv_lam_01 * exp_term;
    double P10 = inv_lam_10 - inv_lam_10 * exp_term;
    double P11 = inv_lam_01 + inv_lam_10 * exp_term;
    double* clPar = buf + par * stride;
    double* clCh  = buf + ch  * stride;

    if (!initFlg[par]) {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
        clPar[offset]     = P00 * cl0 + P01 * cl1;
        clPar[offset + 1] = P10 * cl0 + P11 * cl1;
      }
      initFlg[par] = 1;
    } else {
      for (int c = 0; c < nChar; ++c) {
        int offset = c * kStates;
        double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
        clPar[offset]     *= P00 * cl0 + P01 * cl1;
        clPar[offset + 1] *= P10 * cl0 + P11 * cl1;
      }
    }
  }

  int root = nTip + 1;
  double* clRoot = buf + root * stride;
  double logLik  = 0.0;
  for (int c = 0; c < nChar; ++c) {
    int offset = c * kStates;
    double sl = root_freqs[0] * clRoot[offset] + root_freqs[1] * clRoot[offset + 1];
    if (sl <= 0.0) return R_NegInf;
    logLik += std::log(sl);
  }
  return logLik;
}


static double pruning_mkn_acrv_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    double rate_loss, NumericVector root_freqs,
    NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();
  const int kStates = 2;

#ifndef NDEBUG
  assert_valid_preorder(parent, child, nTip);
#endif

  int maxNode = 2 * nTip - 1;  // OPP-2
  int root   = nTip + 1;
  int clCols = nChar * kStates;

  std::vector<double> site_lik_sum(nChar, 0.0);
  double sum_rl     = 1.0 + rate_loss;
  double rate01     = 2.0 / sum_rl;
  double rate10     = 2.0 * rate_loss / sum_rl;
  double lambda     = rate01 + rate10;
  double inv_lam_01 = rate01 / lambda;
  double inv_lam_10 = rate10 / lambda;

  // OPP-CL: init tips once (constant across rate categories)
  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        cl[offset] = 1.0; cl[offset + 1] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

    for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e];
      int ch  = child[e];
      double t        = edge_length[e] * rate;
      double exp_term = MKP_EXP(-lambda * t);
      double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      double P01 = inv_lam_01 - inv_lam_01 * exp_term;
      double P10 = inv_lam_10 - inv_lam_10 * exp_term;
      double P11 = inv_lam_01 + inv_lam_10 * exp_term;
      double* clPar = buf + par * stride;
      double* clCh  = buf + ch  * stride;

      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
          clPar[offset]     = P00 * cl0 + P01 * cl1;
          clPar[offset + 1] = P10 * cl0 + P11 * cl1;
        }
        initFlg[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double cl0 = clCh[offset]; double cl1 = clCh[offset + 1];
          clPar[offset]     *= P00 * cl0 + P01 * cl1;
          clPar[offset + 1] *= P10 * cl0 + P11 * cl1;
        }
      }
    }

    double* clRoot = buf + root * stride;
    for (int c = 0; c < nChar; ++c) {
      int offset = c * kStates;
      double sl = root_freqs[0] * clRoot[offset] + root_freqs[1] * clRoot[offset + 1];
      site_lik_sum[c] += sl;
    }
  }

  double logLik   = 0.0;
  double inv_nCat = 1.0 / nCat;
  for (int c = 0; c < nChar; ++c) {
    double avg = site_lik_sum[c] * inv_nCat;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// ---------------------------------------------------------------------------
// M-052: Dirichlet-marginal Q-matrix heterogeneity (Het)
//
// compute_het_bins(): discretize Beta(alpha, (k-1)*alpha) into B bins.
//   Each bin is the conditional mean within an equal-probability interval.
//
// pruning_f81_het_acrv_flat(): Felsenstein pruning under a mixture of
//   F81 Q-matrices × ACRV rate categories.  For each (rate, bin, rotation)
//   combination, constructs P(t) analytically and runs a tree traversal.
//   Averages over all R × B × k components.
//
// For k=2, rotations are redundant (Beta(α,α) symmetry) so only one
// rotation is used, halving the cost.
//
// The F81 transition probability is O(1) per entry (no eigendecomposition):
//   P_ij(t) = π_j + (δ_ij − π_j) × exp(−μt)
//   μ = 1 / (1 − Σ π_j²)
// ---------------------------------------------------------------------------

// Discretize Beta(alpha, (k-1)*alpha) into nBins equal-probability bins.
// bins must point to nBins doubles.
static void compute_het_bins(double alpha, int k, int nBins, double* bins) {
  double a = alpha;
  double b = (k - 1.0) * alpha;
  for (int i = 0; i < nBins; ++i) {
    double lo = R::qbeta((double)i / nBins, a, b, 1, 0);
    double hi = R::qbeta((double)(i + 1) / nBins, a, b, 1, 0);
    // Guard: degenerate bin (lo ≈ hi) at extreme alpha.
    // Use midpoint as the bin representative; falls back to 1/k for large α.
    if (hi - lo < 1e-15) {
      bins[i] = 0.5 * (lo + hi);
      continue;
    }
    double p_lo = R::pbeta(lo, a + 1.0, b, 1, 0);
    double p_hi = R::pbeta(hi, a + 1.0, b, 1, 0);
    double denom = R::pbeta(hi, a, b, 1, 0) - R::pbeta(lo, a, b, 1, 0);
    if (denom < 1e-300) {
      bins[i] = 0.5 * (lo + hi);
    } else {
      bins[i] = (a / (a + b)) * (p_hi - p_lo) / denom;
    }
  }
}


// F81 Het + ACRV pruning (flat-buffer).
//
// Parameters:
//   kStates       – number of character states for this group
//   baseRL        – base rate_loss for neomorphic (1.0 for symmetric chars)
//   betaBins      – array of nBetaCat discretized Beta(α,(k-1)α) values
//   nBetaCat      – number of discretization bins (B)
//   rate_multipliers – ACRV rate categories (length nCat)
//
// For each (ACRV cat, Beta bin, rotation) triple, runs one full pruning
// pass using F81 transition probabilities.
//
// For k=2: one rotation suffices (Beta(α,α) symmetry); effective
//   components = nCat × nBetaCat.
// For k≥3: k rotations per bin; effective components = nCat × nBetaCat × k.
//
// Returns raw log-likelihood (no ascertainment correction).
static double pruning_f81_het_acrv_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, double baseRL,
    const double* betaBins, int nBetaCat,
    NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();

#ifndef NDEBUG
  assert_valid_preorder(parent, child, nTip);
#endif

  int maxNode = 2 * nTip - 1;
  int root    = nTip + 1;
  int clCols  = nChar * kStates;

  // For k=2 with symmetric Beta(α,α), rotations are redundant.
  int nRot = (kStates == 2) ? 1 : kStates;

  // Total components: nCat × nBetaCat × nRot
  int totalComp = nCat * nBetaCat * nRot;
  std::vector<double> site_lik_sum(nChar, 0.0);

  // Precompute base gain/loss for neomorphic composition (k=2 only).
  // For symmetric (baseRL = 1.0): gain_base = loss_base = 0.5.
  double gain_base, loss_base;
  if (kStates == 2) {
    double sum_rl = 1.0 + baseRL;
    gain_base = 1.0 / sum_rl;
    loss_base = baseRL / sum_rl;
  } else {
    gain_base = loss_base = 0.0;  // unused for k≥3
  }

  // OPP-CL: init tips once (constant across all cat × bin × rot combos)
  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  for (int cat = 0; cat < nCat; ++cat) {
    double acrvRate = rate_multipliers[cat];

    for (int bi = 0; bi < nBetaCat; ++bi) {
      double beta_val = betaBins[bi];

      for (int rot = 0; rot < nRot; ++rot) {

        // --- Build frequency vector π^(bi,rot) ---
        double pi[16];
        double sumPiSq = 0.0;

        if (kStates == 2) {
          double gain_b = gain_base * 2.0 * beta_val;
          double loss_b = loss_base * 2.0 * (1.0 - beta_val);
          double total_rate = gain_b + loss_b;
          pi[1] = gain_b / total_rate;
          pi[0] = 1.0 - pi[1];
          sumPiSq = pi[0] * pi[0] + pi[1] * pi[1];
        } else {
          double r = (1.0 - beta_val) / (kStates - 1.0);
          for (int s = 0; s < kStates; ++s) pi[s] = r;
          pi[rot] = beta_val;
          sumPiSq = beta_val * beta_val +
                    (kStates - 1.0) * r * r;
        }

        double mu = 1.0 / (1.0 - sumPiSq);

        // Reset internal-node init flags only
        for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

        // --- Tree traversal using F81 P(t) ---
        // P_ij(t) = π_j × (1 − d) + δ_ij × d
        // where d = exp(−μ × acrvRate × t)
        for (int e = nEdge - 1; e >= 0; --e) {
          int par = parent[e];
          int ch  = child[e];
          double t = edge_length[e] * acrvRate;
          double d = MKP_EXP(-mu * t);
          double one_minus_d = 1.0 - d;

          double* clPar = buf + par * stride;
          double* clCh  = buf + ch  * stride;

          // OPP-1 (F81): hoist Σ_j π_j·cl_j out of the i-loop → O(k) per char
          if (!initFlg[par]) {
            for (int c = 0; c < nChar; ++c) {
              int offset = c * kStates;
              double sum_pi_cl = 0.0;
              for (int j = 0; j < kStates; ++j)
                sum_pi_cl += pi[j] * clCh[offset + j];
              double base = one_minus_d * sum_pi_cl;
              for (int i = 0; i < kStates; ++i)
                clPar[offset + i] = base + d * clCh[offset + i];
            }
            initFlg[par] = 1;
          } else {
            for (int c = 0; c < nChar; ++c) {
              int offset = c * kStates;
              double sum_pi_cl = 0.0;
              for (int j = 0; j < kStates; ++j)
                sum_pi_cl += pi[j] * clCh[offset + j];
              double base = one_minus_d * sum_pi_cl;
              for (int i = 0; i < kStates; ++i)
                clPar[offset + i] *= base + d * clCh[offset + i];
            }
          }
        }

        // --- Accumulate site likelihoods (root weighted by π) ---
        double* clRoot = buf + root * stride;
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sl = 0.0;
          for (int s = 0; s < kStates; ++s)
            sl += pi[s] * clRoot[offset + s];
          site_lik_sum[c] += sl;
        }

      }  // rot
    }  // bi
  }  // cat

  // Average and log
  double logLik   = 0.0;
  double inv_comp = 1.0 / totalComp;
  for (int c = 0; c < nChar; ++c) {
    double avg = site_lik_sum[c] * inv_comp;
    if (avg <= 0.0) return R_NegInf;
    logLik += std::log(avg);
  }
  return logLik;
}


// Per-site variant of pruning_f81_het_acrv_flat: fills siteLL[0..nChar-1]
// with per-character log(avg_lik) instead of returning the sum.
// Used by Gibbs kPrime sweep (M-155).
void pruning_f81_het_acrv_persite(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, double baseRL,
    const double* betaBins, int nBetaCat,
    NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride,
    double* siteLL) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  int nCat  = rate_multipliers.size();

  int maxNode = 2 * nTip - 1;
  int root    = nTip + 1;
  int clCols  = nChar * kStates;

  int nRot = (kStates == 2) ? 1 : kStates;
  int totalComp = nCat * nBetaCat * nRot;
  std::vector<double> site_lik_sum(nChar, 0.0);

  double gain_base, loss_base;
  if (kStates == 2) {
    double sum_rl = 1.0 + baseRL;
    gain_base = 1.0 / sum_rl;
    loss_base = baseRL / sum_rl;
  } else {
    gain_base = loss_base = 0.0;
  }

  std::memset(initFlg, 0, (maxNode + 1) * sizeof(uint8_t));
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = buf + tip * stride;
    std::fill(cl, cl + clCols, 0.0);
    for (int c = 0; c < nChar; ++c) {
      int state  = tip_states(tip - 1, c);
      int offset = c * kStates;
      if (state < 0) {
        for (int s = 0; s < kStates; ++s) cl[offset + s] = 1.0;
      } else {
        cl[offset + state] = 1.0;
      }
    }
    initFlg[tip] = 1;
  }

  for (int cat = 0; cat < nCat; ++cat) {
    double acrvRate = rate_multipliers[cat];

    for (int bi = 0; bi < nBetaCat; ++bi) {
      double beta_val = betaBins[bi];

      for (int rot = 0; rot < nRot; ++rot) {

        double pi[16];
        double sumPiSq = 0.0;

        if (kStates == 2) {
          double gain_b = gain_base * 2.0 * beta_val;
          double loss_b = loss_base * 2.0 * (1.0 - beta_val);
          double total_rate = gain_b + loss_b;
          pi[1] = gain_b / total_rate;
          pi[0] = 1.0 - pi[1];
          sumPiSq = pi[0] * pi[0] + pi[1] * pi[1];
        } else {
          double r = (1.0 - beta_val) / (kStates - 1.0);
          for (int s = 0; s < kStates; ++s) pi[s] = r;
          pi[rot] = beta_val;
          sumPiSq = beta_val * beta_val +
                    (kStates - 1.0) * r * r;
        }

        double mu = 1.0 / (1.0 - sumPiSq);

        for (int n = nTip + 1; n <= maxNode; ++n) initFlg[n] = 0;

        for (int e = nEdge - 1; e >= 0; --e) {
          int par = parent[e];
          int ch  = child[e];
          double t = edge_length[e] * acrvRate;
          double d = std::exp(-mu * t);
          double one_minus_d = 1.0 - d;

          double* clPar = buf + par * stride;
          double* clCh  = buf + ch  * stride;

          if (!initFlg[par]) {
            for (int c = 0; c < nChar; ++c) {
              int offset = c * kStates;
              double sum_pi_cl = 0.0;
              for (int j = 0; j < kStates; ++j)
                sum_pi_cl += pi[j] * clCh[offset + j];
              double base = one_minus_d * sum_pi_cl;
              for (int i = 0; i < kStates; ++i)
                clPar[offset + i] = base + d * clCh[offset + i];
            }
            initFlg[par] = 1;
          } else {
            for (int c = 0; c < nChar; ++c) {
              int offset = c * kStates;
              double sum_pi_cl = 0.0;
              for (int j = 0; j < kStates; ++j)
                sum_pi_cl += pi[j] * clCh[offset + j];
              double base = one_minus_d * sum_pi_cl;
              for (int i = 0; i < kStates; ++i)
                clPar[offset + i] *= base + d * clCh[offset + i];
            }
          }
        }

        double* clRoot = buf + root * stride;
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          double sl = 0.0;
          for (int s = 0; s < kStates; ++s)
            sl += pi[s] * clRoot[offset + s];
          site_lik_sum[c] += sl;
        }

      }  // rot
    }  // bi
  }  // cat

  double inv_comp = 1.0 / totalComp;
  for (int c = 0; c < nChar; ++c) {
    double avg = site_lik_sum[c] * inv_comp;
    siteLL[c] = (avg > 0.0) ? std::log(avg) : R_NegInf;
  }
}


// Het-aware ascertainment correction helpers.
// Compute P(constant site) and P(singleton site) averaged over the
// Het × ACRV mixture, for use in variable/informative coding correction.
static double het_constant_site_prob(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, int nTip,
    int kStates, double baseRL,
    const double* betaBins, int nBetaCat,
    NumericVector rate_multipliers) {

  int nCat = rate_multipliers.size();
  int nRot = (kStates == 2) ? 1 : kStates;
  int totalComp = nCat * nBetaCat * nRot;

  double gain_base = 0.0, loss_base = 0.0;
  if (kStates == 2) {
    double sum_rl = 1.0 + baseRL;
    gain_base = 1.0 / sum_rl;
    loss_base = baseRL / sum_rl;
  }

  double totalP = 0.0;
  for (int cat = 0; cat < nCat; ++cat) {
    double acrvRate = rate_multipliers[cat];
    for (int bi = 0; bi < nBetaCat; ++bi) {
      double beta_val = betaBins[bi];
      for (int rot = 0; rot < nRot; ++rot) {
        double pi[16];
        double sumPiSq = 0.0;
        if (kStates == 2) {
          double g = gain_base * 2.0 * beta_val;
          double l = loss_base * 2.0 * (1.0 - beta_val);
          double tot = g + l;
          pi[1] = g / tot; pi[0] = 1.0 - pi[1];
          sumPiSq = pi[0]*pi[0] + pi[1]*pi[1];
        } else {
          double r = (1.0 - beta_val) / (kStates - 1.0);
          for (int s = 0; s < kStates; ++s) pi[s] = r;
          pi[rot] = beta_val;
          sumPiSq = beta_val*beta_val + (kStates-1.0)*r*r;
        }
        double mu = 1.0 / (1.0 - sumPiSq);

        // P(constant site in state s) = π_s × Π_edges P_ss(t)
        // P_ss(t) = π_s + (1 − π_s) exp(−μt)
        // Sum over all states s.
        double compP = 0.0;
        for (int s = 0; s < kStates; ++s) {
          double prod = pi[s];  // root frequency
          for (int e = 0; e < parent.size(); ++e) {
            double t = edge_length[e] * acrvRate;
            double Pss = pi[s] + (1.0 - pi[s]) * MKP_EXP(-mu * t);
            prod *= Pss;
          }
          compP += prod;
        }
        totalP += compP;
      }
    }
  }
  return totalP / totalComp;
}


static double het_singleton_site_prob(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, int nTip,
    int kStates, double baseRL,
    const double* betaBins, int nBetaCat,
    NumericVector rate_multipliers) {

  // For singletons, we use the existing singleton_site_prob_jc/mkn
  // functions by weight-averaging across mixture components.
  // This is correct because the singleton prob is linear in the
  // per-component transition matrices.
  //
  // However, the existing functions hardcode JC/MkN P(t), not F81.
  // For now, use the JC/MkN functions as an approximation.
  // TODO(M-052): implement exact F81 singleton prob if needed.
  //
  // For most datasets with coding="variable", only constant_site_prob
  // matters. Singleton correction ("informative") is Phase 7.
  // Return 0 for now — safe because informative coding is not yet supported.
  return 0.0;
}


// ---------------------------------------------------------------------------
// Single-character JC likelihood helpers (for Gibbs kPrime sweep)
// ---------------------------------------------------------------------------

// Inline JC pruning for one character with ACRV.
// Returns raw (not log) site likelihood averaged over rate categories.
// tipCol: pointer to nTip ints (0-indexed states, -1 = missing).
static double jc1_acrv(
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const int* tipCol, int nTip, int kStates,
    const NumericVector& rates, int nCat) {

  int nEdge = parent.size();
  int maxNode = 2 * nTip - 1;
  int root = nTip + 1;

  std::vector<double> cl((maxNode + 1) * kStates);
  std::vector<uint8_t> flg(maxNode + 1);

  // Init tips once (constant across rate categories)
  std::memset(cl.data(), 0, cl.size() * sizeof(double));
  std::memset(flg.data(), 0, flg.size() * sizeof(uint8_t));
  for (int t = 1; t <= nTip; ++t) {
    double* p = cl.data() + t * kStates;
    int st = tipCol[t - 1];
    if (st < 0) {
      for (int s = 0; s < kStates; ++s) p[s] = 1.0;
    } else {
      p[st] = 1.0;
    }
    flg[t] = 1;
  }

  double inv_k = 1.0 / kStates;
  double km1 = kStates - 1.0;
  double siteLikSum = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rates[cat];

    // Reset internal node init flags (tips stay init'd)
    for (int n = nTip + 1; n <= maxNode; ++n) flg[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      int par = parent[e], ch = child[e];
      double t = edgeLen[e] * rate;
      double ex = std::exp(-kStates * t / km1);
      double ps = inv_k + (1.0 - inv_k) * ex;
      double pd = inv_k - inv_k * ex;
      double dc = ps - pd;

      double* cp = cl.data() + par * kStates;
      double* cc = cl.data() + ch * kStates;

      double sum = 0.0;
      for (int j = 0; j < kStates; ++j) sum += cc[j];

      if (!flg[par]) {
        for (int i = 0; i < kStates; ++i)
          cp[i] = pd * sum + dc * cc[i];
        flg[par] = 1;
      } else {
        for (int i = 0; i < kStates; ++i)
          cp[i] *= pd * sum + dc * cc[i];
      }
    }

    double* clR = cl.data() + root * kStates;
    double sl = 0.0;
    for (int s = 0; s < kStates; ++s) sl += inv_k * clR[s];
    siteLikSum += sl;
  }

  return siteLikSum / nCat;
}


// Constant-site probability for a given kStates, handling JC and Het paths.
double const_site_prob_for_k(
    const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    int kStates, double betaScale,
    const NumericVector& acrvRates) {

  if (data.codingType == 0) return 0.0;

  if (data.qHeterogeneity) {
    double hetBins[16];
    compute_het_bins(betaScale, kStates, data.nBetaCat, hetBins);
    NumericVector rates = acrvRates;
    if (rates.size() == 0) rates = NumericVector(1, 1.0);
    double p = het_constant_site_prob(
      parent, child, edgeLen, data.nTip,
      kStates, 1.0, hetBins, data.nBetaCat, rates);
    // Phase 7: add het_singleton_site_prob when coding == 2
    return p;
  } else {
    NumericVector rootFreqs(kStates, 1.0 / kStates);
    double p = constant_site_prob_jc(parent, child, edgeLen, data.nTip,
                                      kStates, rootFreqs, acrvRates);
    // Phase 7: add singleton_site_prob_jc when coding == 2
    return p;
  }
}


// Full log-likelihood for one transformational character under JC(kStates).
// Handles JC, ACRV, Het/F81, ascertainment correction, and relabeling.
double single_char_loglik_jc(
    const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const int* tipCol,
    int kStates, int kObs,
    double betaScale,
    const NumericVector& acrvRates,
    double constSiteProb) {

  double ll;

  if (data.qHeterogeneity) {
    // F81 Het path: create 1-column matrix and use existing pruning function
    IntegerMatrix sub(data.nTip, 1);
    for (int t = 0; t < data.nTip; ++t) sub(t, 0) = tipCol[t];

    double hetBins[16];
    compute_het_bins(betaScale, kStates, data.nBetaCat, hetBins);

    int maxNode = 2 * data.nTip - 1;
    int tmpStride = kStates;  // 1 character x kStates states
    std::vector<double> tmpBuf((maxNode + 1) * tmpStride, 0.0);
    std::vector<uint8_t> tmpInit(maxNode + 1, 0);

    NumericVector rates = acrvRates;
    if (rates.size() == 0) rates = NumericVector(1, 1.0);

    ll = pruning_f81_het_acrv_flat(
      parent, child, edgeLen, sub,
      kStates, 1.0, hetBins, data.nBetaCat, rates,
      tmpBuf.data(), tmpInit.data(), tmpStride);
  } else {
    // JC path: inline single-character pruning
    double rawLik = jc1_acrv(parent, child, edgeLen, tipCol,
                              data.nTip, kStates, acrvRates,
                              acrvRates.size());
    if (rawLik <= 0.0) return R_NegInf;
    ll = std::log(rawLik);
  }

  // Ascertainment correction
  if (data.codingType != 0) {
    if (constSiteProb >= 1.0) return R_NegInf;
    ll -= std::log(1.0 - constSiteProb);
  }

  // Relabeling correction
  if (data.relabel) {
    ll += mk_prime_relabel_log(kStates, kObs);
  }

  return ll;
}


// ---------------------------------------------------------------------------
// C++ log-likelihood orchestration (mirrors .MkpLogLikelihood in R)
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// cpp_partition_log_likelihood: per-partition log-lik (M-064)
//
// Computes log-likelihood for partition partIdx only, allowing partial
// recomputation when only a subset of partitions is affected by a move.
// ---------------------------------------------------------------------------

double cpp_partition_log_likelihood(
    const McmcData& data, int partIdx,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    double betaScale,
    ClWorkspace* ws) {

  int nTip = data.nTip;
  NumericVector rates = cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ);  // OPP-3
  bool useAcrv = (rateLogSd > 0.0);
  int coding = data.codingType;
  bool useHet = data.qHeterogeneity;

  const PartInfo& part = data.parts[partIdx];
  double ll = 0.0;

  // Determine max node index for workspace fitness check
  int maxNode = 2 * data.nTip - 1;  // OPP-2

  // M-052: precompute Het bins for this partition's k if needed.
  // Stack-allocated for small B (typically 4).
  int nBC = data.nBetaCat;
  double hetBins[16];  // max nBetaCat we'd ever use
  // We'll compute bins lazily below when needed.

  if (part.type == 0) {
    // Neomorphic (kStates = 2): use flat-buffer variant when workspace fits.
    int neededStride = part.tipStates.ncol() * 2;
    bool useWs = ws && ws->fits(maxNode, neededStride);

    NumericVector neoEl(edgeLen.size());
    for (int i = 0; i < edgeLen.size(); ++i) neoEl[i] = edgeLen[i] * rateNeo;

    if (useHet) {
      // M-052: Het path — F81 mixture with rate_loss composition.
      compute_het_bins(betaScale, 2, nBC, hetBins);
      if (!useAcrv) rates = NumericVector(1, 1.0);
      if (useWs) {
        ll = pruning_f81_het_acrv_flat(
          parent, child, neoEl, part.tipStates,
          2, rateLoss, hetBins, nBC, rates,
          ws->buf.data(), ws->init.data(), ws->strideMax);
      } else {
        // Fallback: allocate temporary workspace
        int nNode = maxNode;
        int tmpStride = neededStride;
        std::vector<double> tmpBuf((nNode + 1) * tmpStride, 0.0);
        std::vector<uint8_t> tmpInit(nNode + 1, 0);
        ll = pruning_f81_het_acrv_flat(
          parent, child, neoEl, part.tipStates,
          2, rateLoss, hetBins, nBC, rates,
          tmpBuf.data(), tmpInit.data(), tmpStride);
      }
      if (coding != 0) {
        double p = het_constant_site_prob(
          parent, child, neoEl, nTip, 2, rateLoss, hetBins, nBC, rates);
        if (coding == 2) p += het_singleton_site_prob(
          parent, child, neoEl, nTip, 2, rateLoss, hetBins, nBC, rates);
        ll -= part.tipStates.ncol() * std::log(1.0 - p);
      }
    } else {
      // Homogeneous path (original).
      NumericVector rootFreqs = mkn_stationary(rateLoss);
      if (useWs) {
        ll = useAcrv
          ? pruning_mkn_acrv_flat(parent, child, neoEl, part.tipStates,
                                   rateLoss, rootFreqs, rates,
                                   ws->buf.data(), ws->init.data(), ws->strideMax)
          : pruning_mkn_flat(parent, child, neoEl, part.tipStates,
                              rateLoss, rootFreqs,
                              ws->buf.data(), ws->init.data(), ws->strideMax);
      } else {
        ll = useAcrv ? pruning_mkn_acrv(parent, child, neoEl, part.tipStates,
                                         rateLoss, rootFreqs, rates)
                     : pruning_mkn(parent, child, neoEl, part.tipStates,
                                    rateLoss, rootFreqs);
      }
      if (coding != 0) {
        double p = constant_site_prob_mkn(parent, child, neoEl, nTip,
                                          rateLoss, rootFreqs, rates);
        if (coding == 2) p += singleton_site_prob_mkn(parent, child, neoEl, nTip,
                                                       rateLoss, rootFreqs, rates);
        ll -= part.tipStates.ncol() * std::log(1.0 - p);
      }
    }
  } else if (part.type == 2) {
    // Known state space (kStates = part.k): flat-buffer when workspace fits.
    int kStates = part.k;
    int neededStride = part.tipStates.ncol() * kStates;
    bool useWs = ws && ws->fits(maxNode, neededStride);

    if (useHet) {
      // M-052: Het path for known-k partition (any k).
      compute_het_bins(betaScale, kStates, nBC, hetBins);
      if (!useAcrv) rates = NumericVector(1, 1.0);
      if (useWs) {
        ll = pruning_f81_het_acrv_flat(
          parent, child, edgeLen, part.tipStates,
          kStates, 1.0, hetBins, nBC, rates,
          ws->buf.data(), ws->init.data(), ws->strideMax);
      } else {
        int nNode = maxNode;
        std::vector<double> tmpBuf((nNode + 1) * neededStride, 0.0);
        std::vector<uint8_t> tmpInit(nNode + 1, 0);
        ll = pruning_f81_het_acrv_flat(
          parent, child, edgeLen, part.tipStates,
          kStates, 1.0, hetBins, nBC, rates,
          tmpBuf.data(), tmpInit.data(), neededStride);
      }
      if (coding != 0) {
        double p = het_constant_site_prob(
          parent, child, edgeLen, nTip, kStates, 1.0, hetBins, nBC, rates);
        if (coding == 2) p += het_singleton_site_prob(
          parent, child, edgeLen, nTip, kStates, 1.0, hetBins, nBC, rates);
        ll -= part.tipStates.ncol() * std::log(1.0 - p);
      }
    } else {
      // Homogeneous path (original).
      NumericVector rootFreqs(kStates, 1.0 / kStates);
      if (useWs) {
        ll = useAcrv
          ? pruning_jc_acrv_flat(parent, child, edgeLen, part.tipStates,
                                  kStates, rootFreqs, rates,
                                  ws->buf.data(), ws->init.data(), ws->strideMax)
          : pruning_jc_flat(parent, child, edgeLen, part.tipStates,
                             kStates, rootFreqs,
                             ws->buf.data(), ws->init.data(), ws->strideMax);
      } else {
        ll = useAcrv ? pruning_jc_acrv(parent, child, edgeLen, part.tipStates,
                                        kStates, rootFreqs, rates)
                     : pruning_jc(parent, child, edgeLen, part.tipStates,
                                   kStates, rootFreqs);
      }
      if (coding != 0) {
        double p = constant_site_prob_jc(parent, child, edgeLen, nTip,
                                         kStates, rootFreqs, rates);
        if (coding == 2) p += singleton_site_prob_jc(parent, child, edgeLen, nTip,
                                                      kStates, rootFreqs, rates);
        ll -= part.tipStates.ncol() * std::log(1.0 - p);
      }
    }
  } else {
    // Transformational: sub-group by kPrime value when heterogeneous.
    int nCharPart = part.tipStates.ncol();

    // Fast path: check if all kPrime in this partition are identical.
    // Common when k' = kObs for most characters (prior penalizes large k').
    // Avoids sort_unique, column-scan, and sub-matrix allocation+copy.
    int kp0 = kPrime[part.globalCharIdx[0]];
    bool allSame = true;
    for (int ci = 1; ci < nCharPart; ++ci) {
      if (kPrime[part.globalCharIdx[ci]] != kp0) { allSame = false; break; }
    }

    if (allSame) {
      // All characters share the same k' — use full partition tipStates directly.
      bool wsOk = ws && ws->fits(maxNode, nCharPart * kp0);
      if (useHet) {
        double hetBinsSub[16];
        compute_het_bins(betaScale, kp0, nBC, hetBinsSub);
        if (!useAcrv) rates = NumericVector(1, 1.0);
        if (wsOk) {
          ll += pruning_f81_het_acrv_flat(
            parent, child, edgeLen, part.tipStates,
            kp0, 1.0, hetBinsSub, nBC, rates,
            ws->buf.data(), ws->init.data(), ws->strideMax);
        } else {
          int tmpStride = nCharPart * kp0;
          std::vector<double> tmpBuf((maxNode + 1) * tmpStride, 0.0);
          std::vector<uint8_t> tmpInit(maxNode + 1, 0);
          ll += pruning_f81_het_acrv_flat(
            parent, child, edgeLen, part.tipStates,
            kp0, 1.0, hetBinsSub, nBC, rates,
            tmpBuf.data(), tmpInit.data(), tmpStride);
        }
        if (coding != 0) {
          double p = het_constant_site_prob(
            parent, child, edgeLen, nTip, kp0, 1.0, hetBinsSub, nBC, rates);
          if (coding == 2) p += het_singleton_site_prob(
            parent, child, edgeLen, nTip, kp0, 1.0, hetBinsSub, nBC, rates);
          ll -= nCharPart * std::log(1.0 - p);
        }
      } else {
        NumericVector rootFreqs(kp0, 1.0 / kp0);
        if (wsOk) {
          ll += useAcrv
            ? pruning_jc_acrv_flat(parent, child, edgeLen, part.tipStates,
                                    kp0, rootFreqs, rates,
                                    ws->buf.data(), ws->init.data(), ws->strideMax)
            : pruning_jc_flat(parent, child, edgeLen, part.tipStates,
                               kp0, rootFreqs,
                               ws->buf.data(), ws->init.data(), ws->strideMax);
        } else {
          ll += useAcrv ? pruning_jc_acrv(parent, child, edgeLen, part.tipStates,
                                           kp0, rootFreqs, rates)
                        : pruning_jc(parent, child, edgeLen, part.tipStates,
                                      kp0, rootFreqs);
        }
        if (coding != 0) {
          double p = constant_site_prob_jc(parent, child, edgeLen, nTip,
                                           kp0, rootFreqs, rates);
          if (coding == 2) p += singleton_site_prob_jc(parent, child, edgeLen, nTip,
                                                        kp0, rootFreqs, rates);
          ll -= nCharPart * std::log(1.0 - p);
        }
      }
      if (data.relabel) {
        for (int ci = 0; ci < nCharPart; ++ci)
          ll += mk_prime_relabel_log(kp0, part.kObsLocal[ci]);
      }
    } else {
      // Heterogeneous kPrime: sub-group by value.
      IntegerVector kPrimePart(nCharPart);
      for (int ci = 0; ci < nCharPart; ++ci)
        kPrimePart[ci] = kPrime[part.globalCharIdx[ci]];

      int maxKp = *std::max_element(kPrimePart.begin(), kPrimePart.end());
      bool wsOk = ws && ws->fits(maxNode, nCharPart * maxKp);

      IntegerVector uniqKp = sort_unique(kPrimePart);
      for (int ui = 0; ui < uniqKp.size(); ++ui) {
        int kp = uniqKp[ui];
        std::vector<int> cols;
        for (int ci = 0; ci < nCharPart; ++ci)
          if (kPrimePart[ci] == kp) cols.push_back(ci);
        int nSub = (int)cols.size();
        IntegerMatrix sub(nTip, nSub);
        for (int c = 0; c < nSub; ++c)
          for (int t = 0; t < nTip; ++t) sub(t, c) = part.tipStates(t, cols[c]);

        double subLl;
        if (useHet) {
          double hetBinsSub[16];
          compute_het_bins(betaScale, kp, nBC, hetBinsSub);
          if (!useAcrv) rates = NumericVector(1, 1.0);
          if (wsOk) {
            subLl = pruning_f81_het_acrv_flat(
              parent, child, edgeLen, sub,
              kp, 1.0, hetBinsSub, nBC, rates,
              ws->buf.data(), ws->init.data(), ws->strideMax);
          } else {
            int tmpStride = nSub * kp;
            std::vector<double> tmpBuf((maxNode + 1) * tmpStride, 0.0);
            std::vector<uint8_t> tmpInit(maxNode + 1, 0);
            subLl = pruning_f81_het_acrv_flat(
              parent, child, edgeLen, sub,
              kp, 1.0, hetBinsSub, nBC, rates,
              tmpBuf.data(), tmpInit.data(), tmpStride);
          }
          if (coding != 0) {
            double p = het_constant_site_prob(
              parent, child, edgeLen, nTip, kp, 1.0, hetBinsSub, nBC, rates);
            if (coding == 2) p += het_singleton_site_prob(
              parent, child, edgeLen, nTip, kp, 1.0, hetBinsSub, nBC, rates);
            subLl -= nSub * std::log(1.0 - p);
          }
        } else {
          NumericVector rootFreqs(kp, 1.0 / kp);
          if (wsOk) {
            subLl = useAcrv
              ? pruning_jc_acrv_flat(parent, child, edgeLen, sub,
                                      kp, rootFreqs, rates,
                                      ws->buf.data(), ws->init.data(), ws->strideMax)
              : pruning_jc_flat(parent, child, edgeLen, sub,
                                 kp, rootFreqs,
                                 ws->buf.data(), ws->init.data(), ws->strideMax);
          } else {
            subLl = useAcrv ? pruning_jc_acrv(parent, child, edgeLen, sub,
                                               kp, rootFreqs, rates)
                            : pruning_jc(parent, child, edgeLen, sub,
                                          kp, rootFreqs);
          }
          if (coding != 0) {
            double p = constant_site_prob_jc(parent, child, edgeLen, nTip,
                                             kp, rootFreqs, rates);
            if (coding == 2) p += singleton_site_prob_jc(parent, child, edgeLen, nTip,
                                                          kp, rootFreqs, rates);
            subLl -= nSub * std::log(1.0 - p);
          }
        }
        ll += subLl;
      }
      if (data.relabel) {
        for (int ci = 0; ci < nCharPart; ++ci)
          ll += mk_prime_relabel_log(kPrimePart[ci], part.kObsLocal[ci]);
      }
    }
  }
  return ll;
}


// M-065: accepts parent/child vectors directly — no edge matrix decomposition.
// M-063: optional ClWorkspace* threads flat-buffer workspace through all partitions.
// M-052: betaScale threads the Het parameter (ignored when data.qHeterogeneity is false).
double cpp_log_likelihood(
    const McmcData& data,
    IntegerVector parent,
    IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss,
    double rateLogSd,
    double rateNeo,
    double betaScale,
    ClWorkspace* ws) {

  double totalLoglik = 0.0;
  for (int pi = 0; pi < (int)data.parts.size(); ++pi) {
    totalLoglik += cpp_partition_log_likelihood(
      data, pi, parent, child, edgeLen, kPrime,
      rateLoss, rateLogSd, rateNeo, betaScale, ws);
  }
  return totalLoglik;
}


// ---------------------------------------------------------------------------
// prepare_mcmc_data: convert R mkd + model -> XPtr<McmcData>
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP prepare_mcmc_data(List partitions_r,
                       IntegerVector kObs_r,
                       CharacterVector charTypes_r,
                       bool hasNeo,
                       int nCat,
                       std::string codingStr,
                       bool relabelFlag,
                       double treeLengthShape, double treeLengthRate,
                       double rateLossMeanlog, double rateLossSdlog,
                       double rateLogSdShape,  double rateLogSdRate,
                       double rateNeoMeanlog,  double rateNeoSdlog,
                       double kprimeHyperA,    double kprimeHyperB,
                       bool   kPriorLogseries, double kprimeLogseriesC,
                       bool   kPriorBetaGeometric = false,
                       bool   qHeterogeneity = false,
                       int    nBetaCat = 4,
                       double betaScaleShape = 1.0,
                       double betaScaleRate = 1.0) {
  McmcData* d = new McmcData();
  d->hasNeo = hasNeo;
  d->nCat = nCat;
  // OPP-3: precompute qnorm midpoints once — eliminates nCat qnorm() calls per iteration
  d->acrvZ.resize(nCat);
  for (int i = 0; i < nCat; ++i)
    d->acrvZ[i] = R::qnorm((i + 0.5) / nCat, 0.0, 1.0, 1, 0);
  d->codingType = (codingStr == "none") ? 0 :
                  (codingStr == "variable") ? 1 : 2;
  d->relabel = relabelFlag;
  d->treeLengthShape = treeLengthShape;
  d->treeLengthRate  = treeLengthRate;
  d->rateLossMeanlog = rateLossMeanlog;
  d->rateLossSdlog   = rateLossSdlog;
  d->rateLogSdShape  = rateLogSdShape;
  d->rateLogSdRate   = rateLogSdRate;
  d->rateNeoMeanlog  = rateNeoMeanlog;
  d->rateNeoSdlog    = rateNeoSdlog;
  d->kprimeHyperA    = kprimeHyperA;
  d->kprimeHyperB    = kprimeHyperB;
  d->kPriorLogseries     = kPriorLogseries;
  d->kPriorBetaGeometric = kPriorBetaGeometric;
  d->kprimeLogseriesC    = kprimeLogseriesC;
  d->kObs = kObs_r;
  d->nChar = kObs_r.size();
  d->nTip = 0;

  for (int i = 0; i < charTypes_r.size(); ++i) {
    if (charTypes_r[i] == "transformational") {
      d->transIdxGlobal.push_back(i);
    }
  }

  int nPart = partitions_r.size();
  d->parts.resize(nPart);
  d->charToPartition.assign(d->nChar, -1);
  for (int pi = 0; pi < nPart; ++pi) {
    List p_r = partitions_r[pi];
    PartInfo& pinfo = d->parts[pi];

    std::string ptype = as<std::string>(p_r["type"]);
    pinfo.type = (ptype == "neomorphic") ? 0 :
                 (ptype == "transformational") ? 1 : 2;
    if (pinfo.type == 0) d->neoPartIndices.push_back(pi);

    SEXP k_sexp = p_r["k"];
    pinfo.k = (Rf_isNull(k_sexp) || IntegerVector::is_na(as<int>(k_sexp))) ? 0
                                                                            : as<int>(k_sexp);

    pinfo.tipStates = as<IntegerMatrix>(p_r["tip_states"]);
    if (d->nTip == 0) d->nTip = pinfo.tipStates.nrow();

    IntegerVector ci_r = as<IntegerVector>(p_r["char_indices"]);
    pinfo.globalCharIdx = IntegerVector(ci_r.size());
    for (int ci = 0; ci < ci_r.size(); ++ci) {
      pinfo.globalCharIdx[ci] = ci_r[ci] - 1;
      d->charToPartition[ci_r[ci] - 1] = pi;
    }

    int nCharPart = pinfo.tipStates.ncol();
    pinfo.kObsLocal = IntegerVector(nCharPart);
    for (int ci = 0; ci < nCharPart; ++ci) {
      pinfo.kObsLocal[ci] = kObs_r[pinfo.globalCharIdx[ci]];
    }
  }

  // Initialize branchBins with the default so weighted/block-Gibbs moves
  // work even if set_branch_bins() is never called (e.g. in unit tests).
  d->branchBins.init(d->nBranchBins);

  // M-052: Q-matrix heterogeneity parameters.
  d->qHeterogeneity = qHeterogeneity;
  d->nBetaCat       = nBetaCat;
  d->betaScaleShape = betaScaleShape;
  d->betaScaleRate  = betaScaleRate;

  // Collect distinct k values across all partitions for bin precomputation.
  if (qHeterogeneity) {
    std::set<int> kSet;
    for (int pi = 0; pi < nPart; ++pi) {
      const PartInfo& pinfo = d->parts[pi];
      if (pinfo.type == 0) {
        kSet.insert(2);
      } else if (pinfo.type == 2) {
        kSet.insert(pinfo.k);
      }
      // type 1 (trans): k' varies per character, handled at eval time
    }
    d->hetKValues.assign(kSet.begin(), kSet.end());
  }

  return Rcpp::XPtr<McmcData>(d, true);
}


// M-090: setter for nBranchBins (avoids changing prepare_mcmc_data signature).
// Also precomputes the BranchBins breakpoints so no lazy init is needed later.
// [[Rcpp::export]]
void set_branch_bins(SEXP dataPtr, int nBins) {
  McmcData* d = Rcpp::XPtr<McmcData>(dataPtr);
  d->nBranchBins = nBins;
  d->branchBins.init(nBins);
}
