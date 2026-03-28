#ifndef MKPRIME_MCMC_STATE_H
#define MKPRIME_MCMC_STATE_H

// Shared data structures for C++ MCMC engine.
// Included by mcmc_likelihood.cpp (which defines prepare_mcmc_data and
// cpp_log_likelihood) and by mcmc.cpp (which calls them).

#include <Rcpp.h>
#include <vector>
#include <string>

// Per-partition data: tip states and character mapping.
struct PartInfo {
  int type;          // 0 = neomorphic, 1 = transformational, 2 = known
  int k;             // state space for known partitions; else 0
  Rcpp::IntegerMatrix tipStates; // nTip x nChar (0-indexed, -1 = missing)
  Rcpp::IntegerVector kObsLocal; // kObs per character in partition (for relabelling)
  Rcpp::IntegerVector globalCharIdx; // 0-based map: local char → global kPrime index
};

// Bin breakpoints for weighted branch-length moves (M-087/088/089/054).
// Uses Beta(alpha, beta) quantile breakpoints to tile [0,1].
struct BranchBins {
  static constexpr double kBinAlpha = 0.25;
  static constexpr double kBinBeta  = 0.25;

  int nBins = 0;
  double concentration = 0.0;           // 2 * nBins (Beta concentration for draw)
  std::vector<double> breaks;           // length nBins+1: breaks[0]=0, breaks[nBins]=1
  std::vector<double> mids;             // length nBins: midpoint of each bin

  void init(int n) {
    nBins = n;
    concentration = 2.0 * n;
    breaks.resize(n + 1);
    mids.resize(n);
    breaks[0] = 0.0;
    breaks[n] = 1.0;
    for (int b = 1; b < n; ++b)
      breaks[b] = R::qbeta(static_cast<double>(b) / n, kBinAlpha, kBinBeta, 1, 0);
    for (int b = 0; b < n; ++b)
      mids[b] = 0.5 * (breaks[b] + breaks[b + 1]);
  }
};

// Pre-processed MCMC data: created once before the MCMC loop.
struct McmcData {
  int nTip;
  int nChar;
  bool hasNeo;
  std::vector<PartInfo> parts;
  Rcpp::IntegerVector kObs;           // global kObs, length nChar
  Rcpp::IntegerVector transIdxGlobal; // 0-based indices of trans chars in kPrime

  // Partition index maps for partial likelihood recomputation (M-064)
  std::vector<int> neoPartIndices;    // partition indices where type == 0
  std::vector<int> charToPartition;   // global char idx -> partition index

  // Model/prior parameters
  int nCat;
  int codingType;  // 0 = none, 1 = variable, 2 = informative
  bool relabel;
  std::vector<double> acrvZ;  // OPP-3: precomputed qnorm((i+0.5)/nCat), length nCat

  double treeLengthShape, treeLengthRate;
  double rateLossMeanlog,  rateLossSdlog;
  double rateLogSdShape,   rateLogSdRate;
  double rateNeoMeanlog,   rateNeoSdlog;
  double kprimeHyperA,     kprimeHyperB;

  // Logseries prior for k' (opt-in alternative to hierarchical geometric)
  bool   kPriorLogseries;   // true = log-series; false = hierarchical geometric
  double kprimeLogseriesC;  // c parameter (only used when kPriorLogseries = true)

  // Weighted-move configuration (M-090)
  int nBranchBins = 10;     // number of branch-fraction bins for weighted moves
  BranchBins branchBins;    // precomputed bin breakpoints (init by set_branch_bins)
};

// Pre-allocated flat CL workspace (M-063): eliminates per-call heap
// allocations inside the pruning hot path.
//
// Layout: buf[node * strideMax + c * kStates + s]  (node 1-indexed)
// Sized once at MCMC init with headroom for kPrime growth.
struct ClWorkspace {
  std::vector<double>  buf;
  std::vector<uint8_t> init;  // 0=unset, 1=set; uint8_t avoids vector<bool>
  int nNodeMax  = 0;
  int strideMax = 0;

  bool ready() const { return !buf.empty(); }

  // True iff workspace can serve a call needing nNode nodes and stride cols.
  bool fits(int nNode, int stride) const {
    return ready() && nNode <= nNodeMax && stride <= strideMax;
  }

  void allocate(int nNode, int stride) {
    nNodeMax  = nNode;
    strideMax = stride;
    buf.assign(static_cast<size_t>(nNode + 1) * stride, 0.0);
    init.assign(nNode + 1, 0u);
  }
};


// C++ log-likelihood orchestration — declared here, defined in
// mcmc_likelihood.cpp, called from mcmc.cpp.
// M-065: accepts parent/child vectors directly (no edge matrix round-trip).
// M-063: optional ClWorkspace* eliminates per-call heap allocation.
double cpp_log_likelihood(
    const McmcData& data,
    Rcpp::IntegerVector parent,
    Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss,
    double rateLogSd,
    double rateNeo,
    ClWorkspace* ws = nullptr);


// Per-partition log-likelihood for partial recomputation (M-064).
// M-063: optional ClWorkspace* eliminates per-call heap allocation.
double cpp_partition_log_likelihood(
    const McmcData& data, int partIdx,
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    ClWorkspace* ws = nullptr);
#endif  // MKPRIME_MCMC_STATE_H
