#include <Rcpp.h>
#include "fast_exp.h"
#include "mkn_rates.h"
#include "ascertainment.h"
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
// OPP-2: maxNode = 2*nTip - 1 (covers both rooted and unrooted binary trees)
// OPP-4-lite: single flat std::vector<double> CL replaces vector<vector<double>>
//        (one heap allocation instead of nNode+1 allocations per call)


// ---------------------------------------------------------------------------
// constant_site_prob_jc
//
// P(constant site) for JC(k): sum over k constant patterns.
//
// M-170: JC symmetry — P(all tips = s) is identical for all s, so only
// one pseudo-character is needed (all tips in state 0).  Multiply result
// by kStates to recover the full sum.  This replaces the prior k pseudo-
// character traversals with a single traversal, giving a k× speedup.
// ---------------------------------------------------------------------------

void constant_site_probs_jc_impl(const Rcpp::IntegerVector& parent,
                                 const Rcpp::IntegerVector& child,
                                 const Rcpp::NumericVector& edge_length,
                                 int nTip,
                                 int kStates,
                                 const Rcpp::NumericVector& rate_multipliers,
                                 const std::vector<const uint8_t*>& masks,
                                 double* out) {
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;  // OPP-2
  const int root    = nTip + 1;
  // M-170: one pseudo-character (all observed tips in state 0) per mask
  const int nMask   = masks.empty() ? 1 : (int)masks.size();
  const int stride  = nMask * kStates;

  const double inv_k = 1.0 / kStates;
  const double km1   = kStates - 1.0;

  // OPP-4-lite: one flat allocation reused across all nCat categories
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  // Tip-init hoist: observed tips in state 0 (CL[0] = 1, rest 0 from init);
  // missing tips are marginalised.
  for (int tip = 1; tip <= nTip; ++tip) {
    for (int m = 0; m < nMask; ++m) {
      double* cl = cl_flat.data() + tip * stride + m * kStates;
      if (!masks.empty() && masks[m] && masks[m][tip - 1]) {
        std::fill_n(cl, kStates, 1.0);
      } else {
        cl[0] = 1.0;
      }
    }
    cl_init[tip] = 1;
  }

  std::fill_n(out, nMask, 0.0);

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    // Only reset internal nodes; tips are pre-initialized.
    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t          = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      const double arg        = -kStates * t / km1;
      const double neg_expm1  = -std::expm1(arg);
      const double exp_term   = 1.0 - neg_expm1;
      const double p_same     = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff     = inv_k * neg_expm1;
      const double diff_coeff = p_same - p_diff;  // OPP-1

      for (int m = 0; m < nMask; ++m) {
        double* clPar = cl_flat.data() + par * stride + m * kStates;
        double* clCh  = cl_flat.data() + ch  * stride + m * kStates;

        double sum_cl = 0.0;
        for (int j = 0; j < kStates; ++j) sum_cl += clCh[j];

        if (!cl_init[par]) {
          for (int i = 0; i < kStates; ++i)
            clPar[i] = p_diff * sum_cl + diff_coeff * clCh[i];
        } else {
          for (int i = 0; i < kStates; ++i)
            clPar[i] *= p_diff * sum_cl + diff_coeff * clCh[i];
        }
      }
      cl_init[par] = 1;
    }

    // Root: accumulate site likelihood for each pseudo-character.
    for (int m = 0; m < nMask; ++m) {
      const double* clRoot = cl_flat.data() + root * stride + m * kStates;
      double site_lik = 0.0;
      for (int i = 0; i < kStates; ++i)
        site_lik += inv_k * clRoot[i];  // root_freqs[i] = inv_k
      out[m] += site_lik;
    }
  }

  // Multiply by kStates: compensates for using 1 pseudo-char instead of k.
  for (int m = 0; m < nMask; ++m) out[m] = out[m] * kStates / nCat;
}

// [[Rcpp::export]]
double constant_site_prob_jc(Rcpp::IntegerVector parent,
                             Rcpp::IntegerVector child,
                             Rcpp::NumericVector edge_length,
                             int nTip,
                             int kStates,
                             Rcpp::NumericVector root_freqs,
                             Rcpp::NumericVector rate_multipliers) {
  double p;
  constant_site_probs_jc_impl(parent, child, edge_length, nTip, kStates,
                              rate_multipliers, {}, &p);
  // Return:
  return p;
}


// ---------------------------------------------------------------------------
// singleton_site_prob_jc
//
// P(singleton site) for JC(k): under JC symmetry, k*(k-1) singleton patterns
// at a given tip all have the same probability, so we need only nTip pseudo-
// characters (one per tip, that tip in state 1, all others in state 0).
// P(singleton) = k*(k-1) * sum_j P(all=0, tip_j=1)
// ---------------------------------------------------------------------------

double singleton_site_prob_jc_impl(const Rcpp::IntegerVector& parent,
                                   const Rcpp::IntegerVector& child,
                                   const Rcpp::NumericVector& edge_length,
                                   int nTip,
                                   int kStates,
                                   const Rcpp::NumericVector& root_freqs,
                                   const Rcpp::NumericVector& rate_multipliers,
                                   const uint8_t* missing) {
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

  // Tip-init hoist: tips are identical across rate categories — init once.
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = cl_flat.data() + tip * stride;
    if (missing && missing[tip - 1]) {
      std::fill_n(cl, stride, 1.0);
      cl_init[tip] = 1;
      continue;
    }
    for (int j = 0; j < nTip; ++j) {
      const int offset = j * kStates;
      if (j == tip - 1)
        cl[offset + 1] = 1.0;  // singleton: state 1
      else
        cl[offset + 0] = 1.0;  // background: state 0
    }
    cl_init[tip] = 1;
  }

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    // Only reset internal nodes; tips are pre-initialized.
    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      const double arg      = -kStates * t / km1;
      const double neg_expm1 = -std::expm1(arg);
      const double exp_term = 1.0 - neg_expm1;
      const double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff   = inv_k * neg_expm1;
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
      if (missing && missing[c]) continue;  // a missing tip is no singleton
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

// [[Rcpp::export]]
double singleton_site_prob_jc(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              int kStates,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  return singleton_site_prob_jc_impl(parent, child, edge_length, nTip,
                                     kStates, root_freqs, rate_multipliers,
                                     nullptr);
}


// ---------------------------------------------------------------------------
// constant_site_prob_jc_collapsed
//
// P(constant site) under JC(kFull) collapsed to kEff = kObs + 1 columns.
// By JC symmetry P(all tips = s) is identical for all s, so we evaluate the
// pseudo-character "all tips in observed state 0" on the collapsed tree and
// multiply by kFull (not kEff): the result accounts for both observed and
// lumped constant patterns.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double constant_site_prob_jc_collapsed(Rcpp::IntegerVector parent,
                                       Rcpp::IntegerVector child,
                                       Rcpp::NumericVector edge_length,
                                       int nTip,
                                       int kFull,
                                       int kObs,
                                       Rcpp::NumericVector rate_multipliers) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("constant_site_prob_jc_collapsed requires 1 <= kObs < kFull.");
  }
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;
  const int root    = nTip + 1;
  const int kEff    = kObs + 1;
  const double n_U  = static_cast<double>(kFull - kObs);
  const int stride  = kEff;  // 1 pseudo-character

  const double inv_k = 1.0 / kFull;
  const double km1   = kFull - 1.0;

  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  // All tips in observed state 0: CL[0] = 1, rest 0.
  for (int tip = 1; tip <= nTip; ++tip) {
    cl_flat[tip * stride] = 1.0;
    cl_init[tip] = 1;
  }

  double total_const_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t          = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      const double arg        = -kFull * t / km1;
      const double neg_expm1  = -std::expm1(arg);
      const double exp_term   = 1.0 - neg_expm1;
      const double p_same     = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff     = inv_k * neg_expm1;
      const double diff_coeff = p_same - p_diff;
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      double sum_eff = 0.0;
      for (int j = 0; j < kObs; ++j) sum_eff += clCh[j];
      sum_eff += n_U * clCh[kObs];

      if (!cl_init[par]) {
        for (int i = 0; i < kEff; ++i)
          clPar[i] = p_diff * sum_eff + diff_coeff * clCh[i];
        cl_init[par] = 1;
      } else {
        for (int i = 0; i < kEff; ++i)
          clPar[i] *= p_diff * sum_eff + diff_coeff * clCh[i];
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    double sum_eff_root = 0.0;
    for (int i = 0; i < kObs; ++i) sum_eff_root += clRoot[i];
    sum_eff_root += n_U * clRoot[kObs];
    total_const_prob += inv_k * sum_eff_root;
  }

  // total_const_prob is the per-category mean of P(all tips in state 0).
  // Multiply by kFull to recover sum over all kFull constant patterns.
  return total_const_prob * kFull / nCat;
}


// ---------------------------------------------------------------------------
// singleton_site_prob_jc_collapsed
//
// P(singleton site) under JC(kFull) collapsed to kEff = kObs + 1 columns.
// Constructs nTip pseudo-characters: tip j in column 1 (= observed state 1
// when kObs >= 2, else lumped column when kObs == 1), all others in state 0.
//
// Multiplier:
//   kObs >= 2: kFull * (kFull - 1)
//     Tip j sits in a singleton class; pruning gives P(bg=0, tip_j=1) and by
//     JC symmetry this equals every other (s, s') with s != s' singleton.
//   kObs == 1: kFull
//     Tip j sits in the lumped class (size n_U = kFull - 1); pruning yields
//     n_U * P(specific singleton), so the residual multiplier is
//     kFull*(kFull-1)/n_U = kFull.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
double singleton_site_prob_jc_collapsed(Rcpp::IntegerVector parent,
                                        Rcpp::IntegerVector child,
                                        Rcpp::NumericVector edge_length,
                                        int nTip,
                                        int kFull,
                                        int kObs,
                                        Rcpp::NumericVector rate_multipliers) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("singleton_site_prob_jc_collapsed requires 1 <= kObs < kFull.");
  }
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;
  const int root    = nTip + 1;
  const int nChar   = nTip;
  const int kEff    = kObs + 1;
  const double n_U  = static_cast<double>(kFull - kObs);
  const int stride  = nChar * kEff;

  const double inv_k = 1.0 / kFull;
  const double km1   = kFull - 1.0;

  // singleton tip occupies column 1: observed state 1 if kObs >= 2,
  // otherwise the lumped column (also index 1 because kEff = 2).
  const int singleton_col = 1;

  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = cl_flat.data() + tip * stride;
    for (int j = 0; j < nTip; ++j) {
      const int offset = j * kEff;
      if (j == tip - 1)
        cl[offset + singleton_col] = 1.0;
      else
        cl[offset + 0] = 1.0;
    }
    cl_init[tip] = 1;
  }

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in p_diff at small rt.
      const double arg      = -kFull * t / km1;
      const double neg_expm1 = -std::expm1(arg);
      const double exp_term = 1.0 - neg_expm1;
      const double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      const double p_diff   = inv_k * neg_expm1;
      const double diff_coeff = p_same - p_diff;
      double* clPar = cl_flat.data() + par * stride;
      double* clCh  = cl_flat.data() + ch  * stride;

      if (!cl_init[par]) {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kEff;
          double sum_eff = 0.0;
          for (int j = 0; j < kObs; ++j) sum_eff += clCh[offset + j];
          sum_eff += n_U * clCh[offset + kObs];
          for (int i = 0; i < kEff; ++i)
            clPar[offset + i] = p_diff * sum_eff + diff_coeff * clCh[offset + i];
        }
        cl_init[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          const int offset = c * kEff;
          double sum_eff = 0.0;
          for (int j = 0; j < kObs; ++j) sum_eff += clCh[offset + j];
          sum_eff += n_U * clCh[offset + kObs];
          for (int i = 0; i < kEff; ++i)
            clPar[offset + i] *= p_diff * sum_eff + diff_coeff * clCh[offset + i];
        }
      }
    }

    const double* clRoot = cl_flat.data() + root * stride;
    for (int c = 0; c < nChar; ++c) {
      const int offset = c * kEff;
      double sum_eff = 0.0;
      for (int s = 0; s < kObs; ++s) sum_eff += clRoot[offset + s];
      sum_eff += n_U * clRoot[offset + kObs];
      total_singleton_prob += inv_k * sum_eff;
    }
  }

  total_singleton_prob /= nCat;
  if (kObs >= 2) {
    return total_singleton_prob * kFull * (kFull - 1);
  } else {
    // kObs == 1: singleton tip used the lumped column, which contributed
    // multiplicity n_U = kFull - 1 already. Residual factor = kFull.
    return total_singleton_prob * kFull;
  }
}


// ---------------------------------------------------------------------------
// constant_site_prob_mkn
//
// P(constant site) for the MkN 2-state asymmetric model.
// Two pseudo-characters: one for all-0, one for all-1.
// ---------------------------------------------------------------------------

void constant_site_probs_mkn_impl(const Rcpp::IntegerVector& parent,
                                  const Rcpp::IntegerVector& child,
                                  const Rcpp::NumericVector& edge_length,
                                  int nTip,
                                  double rate_loss,
                                  const Rcpp::NumericVector& root_freqs,
                                  const Rcpp::NumericVector& rate_multipliers,
                                  const std::vector<const uint8_t*>& masks,
                                  double* out) {
  const int nEdge    = parent.size();
  const int nCat     = rate_multipliers.size();
  const int maxNode  = 2 * nTip - 1;  // OPP-2
  const int root     = nTip + 1;
  const int kStates  = 2;
  const int nMask    = masks.empty() ? 1 : (int)masks.size();
  const int nChar    = nMask * kStates;  // all-0 and all-1 per mask
  const int stride   = nChar * kStates;

  double rate01, rate10;
  mkn_rates(rate_loss, rate01, rate10);
  const double lambda     = rate01 + rate10;
  const double inv_lam_01 = rate01 / lambda;
  const double inv_lam_10 = rate10 / lambda;

  // OPP-4-lite: single flat allocation
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  // Tip-init hoist: tips are identical across rate categories — init once.
  for (int tip = 1; tip <= nTip; ++tip) {
    for (int m = 0; m < nMask; ++m) {
      double* cl = cl_flat.data() + tip * stride + m * kStates * kStates;
      if (!masks.empty() && masks[m] && masks[m][tip - 1]) {
        std::fill_n(cl, kStates * kStates, 1.0);
      } else {
        for (int s = 0; s < kStates; ++s)
          cl[s * kStates + s] = 1.0;
      }
    }
    cl_init[tip] = 1;
  }

  std::fill_n(out, nMask, 0.0);

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    // Only reset internal nodes; tips are pre-initialized.
    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in P01/P10 at small lambda*t.
      const double arg      = -lambda * t;
      const double neg_expm1 = -std::expm1(arg);
      const double exp_term = 1.0 - neg_expm1;
      const double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      const double P01 = inv_lam_01 * neg_expm1;
      const double P10 = inv_lam_10 * neg_expm1;
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
    for (int m = 0; m < nMask; ++m) {
      for (int s = 0; s < kStates; ++s) {
        const int offset = (m * kStates + s) * kStates;
        out[m] += root_freqs[0] * clRoot[offset]
                + root_freqs[1] * clRoot[offset + 1];
      }
    }
  }

  for (int m = 0; m < nMask; ++m) out[m] /= nCat;
}

// [[Rcpp::export]]
double constant_site_prob_mkn(Rcpp::IntegerVector parent,
                              Rcpp::IntegerVector child,
                              Rcpp::NumericVector edge_length,
                              int nTip,
                              double rate_loss,
                              Rcpp::NumericVector root_freqs,
                              Rcpp::NumericVector rate_multipliers) {
  double p;
  constant_site_probs_mkn_impl(parent, child, edge_length, nTip, rate_loss,
                               root_freqs, rate_multipliers, {}, &p);
  // Return:
  return p;
}


// ---------------------------------------------------------------------------
// singleton_site_prob_mkn
//
// P(singleton site) for MkN 2-state asymmetric model.
// Unlike JC, states are not symmetric, so we need 2*nTip pseudo-characters:
// first nTip have bg=0/singleton=1; next nTip have bg=1/singleton=0.
// P(singleton) = sum_j [P(all=0, j=1) + P(all=1, j=0)]
// ---------------------------------------------------------------------------

double singleton_site_prob_mkn_impl(const Rcpp::IntegerVector& parent,
                                    const Rcpp::IntegerVector& child,
                                    const Rcpp::NumericVector& edge_length,
                                    int nTip,
                                    double rate_loss,
                                    const Rcpp::NumericVector& root_freqs,
                                    const Rcpp::NumericVector& rate_multipliers,
                                    const uint8_t* missing) {
  const int nEdge   = parent.size();
  const int nCat    = rate_multipliers.size();
  const int maxNode = 2 * nTip - 1;  // OPP-2
  const int root    = nTip + 1;
  const int kStates = 2;
  const int nChar   = 2 * nTip;   // 2n pseudo-characters
  const int stride  = nChar * kStates;

  double rate01, rate10;
  mkn_rates(rate_loss, rate01, rate10);
  const double lambda     = rate01 + rate10;
  const double inv_lam_01 = rate01 / lambda;
  const double inv_lam_10 = rate10 / lambda;

  // OPP-4-lite: single flat allocation
  std::vector<double>  cl_flat((maxNode + 1) * stride, 0.0);
  std::vector<uint8_t> cl_init(maxNode + 1, 0u);

  // Tip-init hoist: tips are identical across rate categories — init once.
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = cl_flat.data() + tip * stride;
    if (missing && missing[tip - 1]) {
      std::fill_n(cl, stride, 1.0);
      cl_init[tip] = 1;
      continue;
    }
    for (int j = 0; j < nTip; ++j) {
      cl[j * kStates + (j == tip - 1 ? 1 : 0)] = 1.0;           // bg=0, single=1
      cl[(nTip + j) * kStates + (j == tip - 1 ? 0 : 1)] = 1.0;  // bg=1, single=0
    }
    cl_init[tip] = 1;
  }

  double total_singleton_prob = 0.0;

  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rate_multipliers[cat];

    // Only reset internal nodes; tips are pre-initialized.
    for (int n = nTip + 1; n <= maxNode; ++n) cl_init[n] = 0;

    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e];
      const int ch  = child[e];
      const double t        = edge_length[e] * rate;
      // FAST-EXP-001: expm1 form avoids cancellation in P01/P10 at small lambda*t.
      const double arg      = -lambda * t;
      const double neg_expm1 = -std::expm1(arg);
      const double exp_term = 1.0 - neg_expm1;
      const double P00 = inv_lam_10 + inv_lam_01 * exp_term;
      const double P01 = inv_lam_01 * neg_expm1;
      const double P10 = inv_lam_10 * neg_expm1;
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
      if (missing && missing[c % nTip]) continue;  // a missing tip is no singleton
      const int offset = c * kStates;
      total_singleton_prob += root_freqs[0] * clRoot[offset]
                            + root_freqs[1] * clRoot[offset + 1];
    }
  }

  return total_singleton_prob / nCat;
}

// [[Rcpp::export]]
double singleton_site_prob_mkn(Rcpp::IntegerVector parent,
                               Rcpp::IntegerVector child,
                               Rcpp::NumericVector edge_length,
                               int nTip,
                               double rate_loss,
                               Rcpp::NumericVector root_freqs,
                               Rcpp::NumericVector rate_multipliers) {
  return singleton_site_prob_mkn_impl(parent, child, edge_length, nTip,
                                      rate_loss, root_freqs, rate_multipliers,
                                      nullptr);
}


// Ascertainment probability for one character with the tips flagged in
// `missing` marginalised: the constant-site probability, plus the singleton
// probability if `informative`. MkN (rate_loss) if `neomorphic`, else JC(k).

// [[Rcpp::export]]
double asc_site_prob_missing(Rcpp::IntegerVector parent,
                             Rcpp::IntegerVector child,
                             Rcpp::NumericVector edge_length,
                             int nTip, int kStates, bool neomorphic,
                             double rate_loss,
                             Rcpp::NumericVector rate_multipliers,
                             Rcpp::LogicalVector missing, bool informative) {
  if (missing.size() != nTip) Rcpp::stop("`missing` must have length nTip");
  std::vector<uint8_t> miss(nTip);
  for (int t = 0; t < nTip; ++t) miss[t] = missing[t] == TRUE;
  double p;
  if (neomorphic) {
    const Rcpp::NumericVector rootFreqs = Rcpp::NumericVector::create(
      rate_loss / (1.0 + rate_loss), 1.0 / (1.0 + rate_loss));
    constant_site_probs_mkn_impl(parent, child, edge_length, nTip, rate_loss,
                                 rootFreqs, rate_multipliers, {miss.data()},
                                 &p);
    if (informative)
      p += singleton_site_prob_mkn_impl(parent, child, edge_length, nTip,
                                        rate_loss, rootFreqs,
                                        rate_multipliers, miss.data());
  } else {
    constant_site_probs_jc_impl(parent, child, edge_length, nTip, kStates,
                                rate_multipliers, {miss.data()}, &p);
    if (informative) {
      const Rcpp::NumericVector rootFreqs(kStates, 1.0 / kStates);
      p += singleton_site_prob_jc_impl(parent, child, edge_length, nTip,
                                       kStates, rootFreqs, rate_multipliers,
                                       miss.data());
    }
  }
  // Return:
  return p;
}
