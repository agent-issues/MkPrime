# #18 (T-002): does singleton_site_prob_jc's per-call heap buffer cost enough
# to be worth threading a workspace through?
#
# Isolated micro-benchmark: the body of singleton_site_prob_jc_impl
# (src/ascertainment.cpp) compiled twice -- as-is (fresh zero-filled
# std::vector per call) and with caller-owned scratch (tips re-zeroed, internal
# nodes overwritten on first visit, as the real code already does). Timed
# in-C++ over many calls to keep R call overhead out. Then scaled by the number
# of calls per likelihood evaluation on Sun2018 under coding = "informative",
# and compared against one full likelihood evaluation.
#
# Usage: Rscript dev/profiling/drivers/18_singleton_alloc.R [lib]

args <- commandArgs(TRUE)
if (length(args)) .libPaths(c(args[1], .libPaths()))
suppressPackageStartupMessages({
  library(MkPrime)
  library(TreeTools)
})

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
#include <chrono>
using namespace Rcpp;

template <bool reuse>
double Singleton(const IntegerVector& parent, const IntegerVector& child,
                 const NumericVector& edge_length, int nTip, int kStates,
                 const NumericVector& root_freqs, const NumericVector& rates,
                 std::vector<double>& wsFlat, std::vector<uint8_t>& wsInit) {
  const int nEdge = parent.size(), nCat = rates.size();
  const int maxNode = 2 * nTip - 1, root = nTip + 1, nChar = nTip;
  const int stride = nChar * kStates;
  const double inv_k = 1.0 / kStates, km1 = kStates - 1.0;
  std::vector<double> localFlat;
  std::vector<uint8_t> localInit;
  double* flat;
  uint8_t* init;
  if (reuse) {
    const size_t need = (size_t)(maxNode + 1) * stride;
    if (wsFlat.size() < need) wsFlat.resize(need);
    if (wsInit.size() < (size_t)maxNode + 1) wsInit.resize(maxNode + 1);
    flat = wsFlat.data();
    init = wsInit.data();
    std::fill_n(flat + stride, (size_t)nTip * stride, 0.0);
  } else {
    localFlat.assign((size_t)(maxNode + 1) * stride, 0.0);
    localInit.assign(maxNode + 1, 0u);
    flat = localFlat.data();
    init = localInit.data();
  }
  for (int tip = 1; tip <= nTip; ++tip) {
    double* cl = flat + tip * stride;
    for (int j = 0; j < nTip; ++j) cl[j * kStates + (j == tip - 1)] = 1.0;
    init[tip] = 1;
  }
  double total = 0.0;
  for (int cat = 0; cat < nCat; ++cat) {
    const double rate = rates[cat];
    for (int n = nTip + 1; n <= maxNode; ++n) init[n] = 0;
    for (int e = nEdge - 1; e >= 0; --e) {
      const int par = parent[e], ch = child[e];
      const double t = edge_length[e] * rate;
      const double neg_expm1 = -std::expm1(-kStates * t / km1);
      const double p_same = inv_k + (1.0 - inv_k) * (1.0 - neg_expm1);
      const double p_diff = inv_k * neg_expm1;
      const double diff_coeff = p_same - p_diff;
      double* clPar = flat + par * stride;
      double* clCh = flat + ch * stride;
      for (int c = 0; c < nChar; ++c) {
        const int off = c * kStates;
        double s = 0.0;
        for (int j = 0; j < kStates; ++j) s += clCh[off + j];
        for (int i = 0; i < kStates; ++i) {
          const double v = p_diff * s + diff_coeff * clCh[off + i];
          if (init[par]) clPar[off + i] *= v; else clPar[off + i] = v;
        }
      }
      init[par] = 1;
    }
    const double* clRoot = flat + root * stride;
    for (int c = 0; c < nChar; ++c)
      for (int s = 0; s < kStates; ++s)
        total += root_freqs[s] * clRoot[c * kStates + s];
  }
  return total / nCat * kStates * (kStates - 1);
}

// [[Rcpp::export]]
NumericVector BenchSingleton(IntegerVector parent, IntegerVector child,
                             NumericVector edgeLen, int nTip, int k,
                             NumericVector rates, int nRep) {
  NumericVector rf(k, 1.0 / k);
  std::vector<double> wsF;
  std::vector<uint8_t> wsI;
  double a = 0, b = 0;
  auto t0 = std::chrono::steady_clock::now();
  for (int r = 0; r < nRep; ++r)
    a += Singleton<false>(parent, child, edgeLen, nTip, k, rf, rates, wsF, wsI);
  auto t1 = std::chrono::steady_clock::now();
  for (int r = 0; r < nRep; ++r)
    b += Singleton<true>(parent, child, edgeLen, nTip, k, rf, rates, wsF, wsI);
  auto t2 = std::chrono::steady_clock::now();
  if (a != b) stop("workspace variant changed the result");
  return NumericVector::create(
    std::chrono::duration<double, std::micro>(t1 - t0).count() / nRep,
    std::chrono::duration<double, std::micro>(t2 - t1).count() / nRep,
    Singleton<false>(parent, child, edgeLen, nTip, k, rf, rates, wsF, wsI));
}
')

nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
if (!nzchar(nexFile)) {
  nexFile <- file.path(Sys.getenv("TREESEARCH_SRC"),
                       "inst/datasets/Sun2018.nex")
}
sun <- ReadAsPhyDat(nexFile)
nTip <- length(sun)

tree <- Preorder(RenumberTips(NJTree(sun, edgeLengths = TRUE), names(sun)))
tree$edge.length[tree$edge.length <= 0] <- 1e-8
parent <- tree$edge[, 1]
child <- tree$edge[, 2]
nCat <- 4L
rates <- MkPrime:::DiscreteLognormalRates(0.5, nCat)

# Sanity: the kernel copy must reproduce the package function.
ref <- MkPrime:::singleton_site_prob_jc(parent, child, tree$edge.length,
                                        nTip, 2L, c(0.5, 0.5), rates)
nRep <- 2000L
timings <- replicate(21, BenchSingleton(parent, child, tree$edge.length,
                                       nTip, 2L, rates, nRep))
stopifnot(abs(timings[3, 1] / ref - 1) < 1e-12)
perCall <- apply(timings[1:2, ], 1, median)
spread <- apply(timings[1:2, ], 1, mad)
pairedDelta <- timings[1, ] - timings[2, ]
cat(sprintf("per call, fresh buffer : %.2f us (mad %.2f)\n",
            perCall[1], spread[1]))
cat(sprintf("per call, reused buffer: %.2f us (mad %.2f)\n",
            perCall[2], spread[2]))
cat(sprintf("paired delta: median %.2f us, IQR %.2f to %.2f us\n",
            median(pairedDelta), quantile(pairedDelta, 0.25),
            quantile(pairedDelta, 0.75)))

# Calls per likelihood evaluation (k' = kObs): one per binary
# transformational partition for the complete-data term, plus one per
# distinct missing-data mask among those characters (asc_probs_masked,
# memoised per (k, mask) within an evaluation).
mkd <- suppressWarnings(MkPrimeData(sun))
binParts <- Filter(function(p) p$type != "neomorphic" && p$kObs == 2L,
                   mkd$partitions)
masks <- unlist(lapply(binParts, function(p) {
  apply(is.na(p$tip_states) | p$tip_states < 0, 2, paste, collapse = "")
}))
nCalls <- length(binParts) +
  length(setdiff(unique(masks), strrep("FALSE", nTip)))
cat(sprintf("binary trans chars: %d; singleton calls per eval: %d\n",
            length(masks), nCalls))

# One full likelihood evaluation, same tree, for scale.
dataPtr <- MkPrime:::prepare_mcmc_data(
  partitions_r = mkd$partitions, kObs_r = mkd$kObs, charTypes_r = mkd$type,
  hasNeo = FALSE, nCat = nCat, codingStr = "informative",
  relabelFlag = TRUE, treeLengthShape = 1.5, treeLengthRate = 1,
  rateLossMeanlog = 0, rateLossSdlog = 1, rateLogSdShape = 1,
  rateLogSdRate = 1, rateNeoMeanlog = 0, rateNeoSdlog = 2,
  kprimeHyperA = 1, kprimeHyperB = 1, kPriorLogseries = TRUE,
  kprimeLogseriesC = 0.7
)
EvalOnce <- function() {
  MkPrime:::cpp_log_likelihood_xptr(dataPtr, parent, child,
                                    tree$edge.length,
                                    as.integer(mkd$kObs), rateLoss = 1,
                                    rateLogSd = 0.5, rateNeo = 1)
}
invisible(EvalOnce())
evalMs <- median(replicate(7, system.time(for (i in 1:20) EvalOnce())[[3]]))
evalUs <- evalMs / 20 * 1e6
saving <- nCalls * diff(rev(perCall))
cat(sprintf("full likelihood eval: %.0f us\n", evalUs))
cat(sprintf("saving per eval: %.2f us = %.3f %% of an eval\n",
            saving, 100 * saving / evalUs))
