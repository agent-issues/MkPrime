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

static NumericVector cpp_acrv_rates(double rateLogSd, int nCat) {
  if (rateLogSd <= 0.0) return NumericVector(nCat, 1.0);
  double mu = -rateLogSd * rateLogSd / 2.0;
  NumericVector rates(nCat);
  double total = 0.0;
  for (int i = 0; i < nCat; ++i) {
    double mid = ((double)i + 0.5) / nCat;
    double z = R::qnorm(mid, 0.0, 1.0, 1, 0);
    rates[i] = std::exp(mu + rateLogSd * z);
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
    double rateLoss, double rateLogSd, double rateNeo) {

  int nTip = data.nTip;
  NumericVector rates = cpp_acrv_rates(rateLogSd, data.nCat);
  bool useAcrv = (rateLogSd > 0.0);
  int coding = data.codingType;

  const PartInfo& part = data.parts[partIdx];
  double ll = 0.0;

  if (part.type == 0) {
    NumericVector neoEl(edgeLen.size());
    for (int i = 0; i < edgeLen.size(); ++i) neoEl[i] = edgeLen[i] * rateNeo;
    NumericVector rootFreqs = mkn_stationary(rateLoss);
    ll = useAcrv ? pruning_mkn_acrv(parent, child, neoEl, part.tipStates,
                                     rateLoss, rootFreqs, rates)
                 : pruning_mkn(parent, child, neoEl, part.tipStates,
                                rateLoss, rootFreqs);
    if (coding != 0) {
      double p = constant_site_prob_mkn(parent, child, neoEl, nTip,
                                        rateLoss, rootFreqs, rates);
      if (coding == 2) p += singleton_site_prob_mkn(parent, child, neoEl, nTip,
                                                     rateLoss, rootFreqs, rates);
      ll -= part.tipStates.ncol() * std::log(1.0 - p);
    }
  } else if (part.type == 2) {
    int kStates = part.k;
    NumericVector rootFreqs(kStates, 1.0 / kStates);
    ll = useAcrv ? pruning_jc_acrv(parent, child, edgeLen, part.tipStates,
                                    kStates, rootFreqs, rates)
                 : pruning_jc(parent, child, edgeLen, part.tipStates,
                               kStates, rootFreqs);
    if (coding != 0) {
      double p = constant_site_prob_jc(parent, child, edgeLen, nTip,
                                       kStates, rootFreqs, rates);
      if (coding == 2) p += singleton_site_prob_jc(parent, child, edgeLen, nTip,
                                                    kStates, rootFreqs, rates);
      ll -= part.tipStates.ncol() * std::log(1.0 - p);
    }
  } else {
    int nCharPart = part.tipStates.ncol();
    IntegerVector kPrimePart(nCharPart);
    for (int ci = 0; ci < nCharPart; ++ci)
      kPrimePart[ci] = kPrime[part.globalCharIdx[ci]];
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
      double subLl = useAcrv ? pruning_jc_acrv(parent, child, edgeLen, sub,
                                                 kp, rootFreqs, rates)
                              : pruning_jc(parent, child, edgeLen, sub,
                                            kp, rootFreqs);
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
double cpp_log_likelihood(
    const McmcData& data,
    IntegerVector parent,
    IntegerVector child,
    NumericVector edgeLen,
    const IntegerVector& kPrime,
    double rateLoss,
    double rateLogSd,
    double rateNeo) {

  double totalLoglik = 0.0;
  for (int pi = 0; pi < (int)data.parts.size(); ++pi) {
    totalLoglik += cpp_partition_log_likelihood(
      data, pi, parent, child, edgeLen, kPrime,
      rateLoss, rateLogSd, rateNeo);
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
