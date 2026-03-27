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
#include <cmath>
#include <algorithm>

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
static NumericVector cpp_acrv_rates(double rateLogSd, int nCat,
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
  f[0] = 1.0 / (1.0 + rateLoss);
  f[1] = rateLoss / (1.0 + rateLoss);
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
// ---------------------------------------------------------------------------

static double pruning_jc_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    int kStates, NumericVector root_freqs,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();

  int maxNode = 2 * nTip - 2;  // OPP-2: unrooted binary tree; eliminates O(nEdge) scan
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

  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    int ch  = child[e];
    double t        = edge_length[e];
    double exp_term = std::exp(-kStates * t / km1);
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

  int maxNode = 2 * nTip - 2;  // OPP-2
  int root   = nTip + 1;
  int clCols = nChar * kStates;

  std::vector<double> site_lik_sum(nChar, 0.0);
  double inv_k = 1.0 / kStates;
  double km1   = kStates - 1.0;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

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

    for (int e = 0; e < nEdge; ++e) {
      int par = parent[e];
      int ch  = child[e];
      double t        = edge_length[e] * rate;
      double exp_term = std::exp(-kStates * t / km1);
      double p_same   = inv_k + (1.0 - inv_k) * exp_term;
      double p_diff   = inv_k - inv_k * exp_term;
      double* clPar   = buf + par * stride;
      double* clCh    = buf + ch  * stride;

      if (!initFlg[par]) {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          for (int i = 0; i < kStates; ++i) {
            double sum = 0.0;
            for (int j = 0; j < kStates; ++j)
              sum += ((i == j) ? p_same : p_diff) * clCh[offset + j];
            clPar[offset + i] = sum;
          }
        }
        initFlg[par] = 1;
      } else {
        for (int c = 0; c < nChar; ++c) {
          int offset = c * kStates;
          for (int i = 0; i < kStates; ++i) {
            double sum = 0.0;
            for (int j = 0; j < kStates; ++j)
              sum += ((i == j) ? p_same : p_diff) * clCh[offset + j];
            clPar[offset + i] *= sum;
          }
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


static double pruning_mkn_flat(
    IntegerVector parent, IntegerVector child,
    NumericVector edge_length, IntegerMatrix tip_states,
    double rate_loss, NumericVector root_freqs,
    double* buf, uint8_t* initFlg, int stride) {

  int nEdge = parent.size();
  int nTip  = tip_states.nrow();
  int nChar = tip_states.ncol();
  const int kStates = 2;

  int maxNode = 2 * nTip - 2;  // OPP-2: unrooted binary tree; eliminates O(nEdge) scan
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

  for (int e = 0; e < nEdge; ++e) {
    int par = parent[e];
    int ch  = child[e];
    double t        = edge_length[e];
    double exp_term = std::exp(-lambda * t);
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

  int maxNode = 2 * nTip - 2;  // OPP-2
  int root   = nTip + 1;
  int clCols = nChar * kStates;

  std::vector<double> site_lik_sum(nChar, 0.0);
  double sum_rl     = 1.0 + rate_loss;
  double rate01     = 2.0 / sum_rl;
  double rate10     = 2.0 * rate_loss / sum_rl;
  double lambda     = rate01 + rate10;
  double inv_lam_01 = rate01 / lambda;
  double inv_lam_10 = rate10 / lambda;

  for (int cat = 0; cat < nCat; ++cat) {
    double rate = rate_multipliers[cat];

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
          cl[offset] = 1.0; cl[offset + 1] = 1.0;
        } else {
          cl[offset + state] = 1.0;
        }
      }
      initFlg[tip] = 1;
    }

    for (int e = 0; e < nEdge; ++e) {
      int par = parent[e];
      int ch  = child[e];
      double t        = edge_length[e] * rate;
      double exp_term = std::exp(-lambda * t);
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
    ClWorkspace* ws) {

  int nTip = data.nTip;
  NumericVector rates = cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ);  // OPP-3
  bool useAcrv = (rateLogSd > 0.0);
  int coding = data.codingType;

  const PartInfo& part = data.parts[partIdx];
  double ll = 0.0;

  // Determine max node index for workspace fitness check
  int maxNode = 2 * data.nTip - 2;  // OPP-2

  if (part.type == 0) {
    // Neomorphic (kStates = 2): use flat-buffer variant when workspace fits.
    int neededStride = part.tipStates.ncol() * 2;
    bool useWs = ws && ws->fits(maxNode, neededStride);

    NumericVector neoEl(edgeLen.size());
    for (int i = 0; i < edgeLen.size(); ++i) neoEl[i] = edgeLen[i] * rateNeo;
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
  } else if (part.type == 2) {
    // Known state space (kStates = part.k): flat-buffer when workspace fits.
    int kStates = part.k;
    int neededStride = part.tipStates.ncol() * kStates;
    bool useWs = ws && ws->fits(maxNode, neededStride);

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
  } else {
    // Transformational: loop over sub-groups sharing the same kPrime value.
    // Workspace fitness: conservative upper bound nCharPart * max(kPrimePart).
    int nCharPart = part.tipStates.ncol();
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
      NumericVector rootFreqs(kp, 1.0 / kp);

      // Sub-call stride: nSub * kp <= nCharPart * maxKp, so wsOk implies it fits.
      double subLl;
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
      ll += subLl;
    }
    if (data.relabel) {
      for (int ci = 0; ci < nCharPart; ++ci)
        ll += mk_prime_relabel_log(kPrimePart[ci], part.kObsLocal[ci]);
    }
  }
  return ll;
}


// M-065: accepts parent/child vectors directly — no edge matrix decomposition.
// M-063: optional ClWorkspace* threads flat-buffer workspace through all partitions.
double cpp_log_likelihood(
    const McmcData& data,
    IntegerVector parent,
    IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss,
    double rateLogSd,
    double rateNeo,
    ClWorkspace* ws) {

  double totalLoglik = 0.0;
  for (int pi = 0; pi < (int)data.parts.size(); ++pi) {
    totalLoglik += cpp_partition_log_likelihood(
      data, pi, parent, child, edgeLen, kPrime,
      rateLoss, rateLogSd, rateNeo, ws);
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
                       double kprimeHyperA,    double kprimeHyperB) {
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

  return Rcpp::XPtr<McmcData>(d, true);
}
