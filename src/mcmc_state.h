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
  int classIdx = 1;  // partition-API user class membership (1-based);
                     // 1 for everything when partition = NULL (legacy path).
                     // Read from the R-side partitions list by prepare_mcmc_data.
  Rcpp::IntegerMatrix tipStates; // nTip x nChar (0-indexed, -1 = missing)
  Rcpp::IntegerVector kObsLocal; // kObs per character in partition (for relabelling)
  Rcpp::IntegerVector globalCharIdx; // 0-based map: local char → global kPrime index

  // M-172: pattern-compressed tip states for Gibbs sweep deduplication.
  // uniqueTipStates holds one column per unique tip-pattern; patternIndex[c]
  // is the 0-based index into uniqueTipStates for the c-th local character.
  // All characters sharing a pattern have identical likelihoods under any
  // (k, tree, params) under the symmetric JC model.
  Rcpp::IntegerMatrix uniqueTipStates; // nTip x nUniquePatterns
  Rcpp::IntegerVector patternIndex;    // length nChar, 0-based
  int nUniquePatterns = 0;
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
  // Character counts per partition type (populated by prepare_mcmc_data).
  //   nNeo   = total chars in type==0 (neomorphic) partitions
  //   nTrans = total chars in type==1 (transformational) + type==2 (known) partitions
  // Used by compute_partition_scales() to apply RB-style partition-rate
  // normalisation so the nChar-weighted mean rate is 1 (audit Issue 1).
  int nNeo   = 0;
  int nTrans = 0;
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

  // kPrime prior selection (exactly one of these is true; when all false:
  // hierarchical geometric with shared p)
  bool   kPriorLogseries;          // log-series with fixed c
  bool   kPriorBetaGeometric;      // per-character Beta-Geometric (α, β)
  bool   kPriorEmpiricalGeometric; // convolution of empirical(N_obs) and Geom(p)
  double kprimeLogseriesC;         // c parameter (only used when kPriorLogseries)

  // Empirical prior parameters (only used when kPriorEmpiricalGeometric).
  // Pre-computed log P_emp(k) for k = 2, ..., empBodyLastK; entry idx = k - 2.
  // The tail extends k >= empTailStartK with log mass
  //   empLogTailStartP + (k - empTailStartK) * log(empTailDecay).
  std::vector<double> empLogBody;  // length empBodyLastK - 1
  int    empBodyLastK = 1;         // largest k with explicit body mass
  int    empTailStartK = 0;        // smallest k in geometric tail (0 = no tail)
  double empTailDecay = 0.0;       // q in (0, 1); 0 means no tail
  double empLogTailStartP = R_NegInf;  // log mass at empTailStartK

  // Weighted-move configuration (M-090)
  int nBranchBins = 10;     // number of branch-fraction bins for weighted moves
  BranchBins branchBins;    // precomputed bin breakpoints (init by set_branch_bins)

  // Partition-API state (Layer 1, plan v4).
  // nClasses == 1 corresponds to the legacy partition = NULL path: every
  // PartInfo carries classIdx = 1 and the user-class layer collapses into a
  // no-op. nClasses > 1 only when a user partition was supplied; in that
  // case classIdx on each PartInfo is read from the R-side partitions list.
  int nClasses = 1;

  // Partition-API (plan v4 §5.1): Dirichlet(α) concentration on class_w.
  // Default 1.0 = flat Dirichlet. Tunable via MkPrimeModel(classRateConcentration).
  double classRateConcentration = 1.0;

  // Q-matrix heterogeneity (M-052): Dirichlet-marginal discretization.
  // When enabled, characters evolve under a mixture of F81 Q-matrices with
  // equilibrium frequencies drawn from a discretized symmetric Dirichlet.
  // The marginal of one component of Dir(α,...,α) with k components is
  // Beta(α, (k-1)α).  Bins are quantile midpoints; rotations over k states
  // enforce labelling symmetry.  Single parameter: betaScale (= α).
  bool qHeterogeneity = false;
  int  nBetaCat = 4;        // number of discretization bins (B)
  double betaScaleShape = 1.0;  // Gamma prior shape for beta_scale
  double betaScaleRate  = 1.0;  // Gamma prior rate  for beta_scale
  // Distinct k values present in the dataset (populated at init).
  // Used to precompute per-k bins when beta_scale changes.
  std::vector<int> hetKValues;  // e.g., {2, 3, 5}
};

// ---------------------------------------------------------------------------
// Partition-rate scales (audit Issue 1: RB-style nChar-weighted-mean-1)
//
// MkPrime's rate_neo originally entered only as a one-sided neo multiplier:
//   neoEl = edgeLen * rate_neo;  transEl = edgeLen
// which gave nChar-weighted mean partition rate (n_neo * r + n_trans) / nChar
// ≠ 1 unless r = 1. This made `tree_length` lose its "expected substitutions
// per character" interpretation under data asymmetry.
//
// Symmetric RB-style formula (Mk' on RB side; identity-1 by construction):
//   neoScale   = r/(1+r) * nTotal / nNeo
//   transScale = 1/(1+r) * nTotal / nTrans
// Weighted mean: (nNeo*neoScale + nTrans*transScale)/nTotal == 1 exactly.
//
// Degenerate cases (no free parameter): if nNeo == 0 or nTrans == 0, the
// rate_neo parameter has no effect on the likelihood and both scales are 1.
struct PartitionScales {
  double neo;    // multiplier for type==0 (neomorphic) edge lengths
  double trans;  // multiplier for type==1 / type==2 edge lengths
};

static inline PartitionScales compute_partition_scales(
    double rateNeo, int nNeo, int nTrans) {
  PartitionScales s;
  if (nNeo == 0 || nTrans == 0) {
    s.neo = 1.0;
    s.trans = 1.0;
    return s;
  }
  const double denom = 1.0 + rateNeo;
  const double nTotal = static_cast<double>(nNeo + nTrans);
  s.neo   = rateNeo / denom * nTotal / static_cast<double>(nNeo);
  s.trans = 1.0     / denom * nTotal / static_cast<double>(nTrans);
  return s;
}


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

  // M-162B: pre-allocated scratch buffers to avoid per-call heap allocation
  std::vector<double>  ascBuf;      // fused ascertainment pseudo-CL
  std::vector<uint8_t> ascInit;     // fused ascertainment init flags
  int ascStrideMax = 0;             // max pseudo-char stride
  std::vector<double>  siteLikSum;  // per-character site likelihood accumulator
  int siteLikMax = 0;

  bool ready() const { return !buf.empty(); }

  // True iff workspace can serve a call needing nNode nodes and stride cols.
  bool fits(int nNode, int stride) const {
    return ready() && nNode <= nNodeMax && stride <= strideMax;
  }

  // True iff ascertainment scratch fits (nNode nodes, ascStride per node).
  bool ascFits(int nNode, int ascStride) const {
    return !ascBuf.empty() && nNode <= nNodeMax && ascStride <= ascStrideMax;
  }

  // True iff siteLikSum fits nChar entries.
  bool siteLikFits(int nChar) const {
    return !siteLikSum.empty() && nChar <= siteLikMax;
  }

  void allocate(int nNode, int stride) {
    nNodeMax  = nNode;
    strideMax = stride;
    buf.assign(static_cast<size_t>(nNode + 1) * stride, 0.0);
    init.assign(nNode + 1, 0u);
  }

  // M-162B: allocate scratch buffers for fused ascertainment + site_lik_sum.
  void allocateScratch(int nNode, int ascStride, int nSiteLik) {
    if (ascStride > 0) {
      ascStrideMax = ascStride;
      ascBuf.assign(static_cast<size_t>(nNode + 1) * ascStride, 0.0);
      ascInit.assign(nNode + 1, 0u);
    }
    if (nSiteLik > 0) {
      siteLikMax = nSiteLik;
      siteLikSum.assign(nSiteLik, 0.0);
    }
  }
};


// C++ log-likelihood orchestration — declared here, defined in
// mcmc_likelihood.cpp, called from mcmc.cpp.
// M-065: accepts parent/child vectors directly (no edge matrix round-trip).
// M-063: optional ClWorkspace* eliminates per-call heap allocation.
// M-052: betaScale controls the Dirichlet-marginal Q-matrix mixture.
// When data.qHeterogeneity is false, betaScale is ignored by the dispatch
// logic, so callers don't need to conditionally omit it.
double cpp_log_likelihood(
    const McmcData& data,
    Rcpp::IntegerVector parent,
    Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss,
    double rateLogSd,
    double rateNeo,
    double betaScale = 1.0,
    ClWorkspace* ws = nullptr);


// Partition-aware sibling (plan v4 §6, Layer 1). Adds per-class state for
// the "shape" (rateLogSd) and "ratemultiplier" (classRate) unlink tokens.
// Length-1 inputs collapse to the legacy scalar path; passing
// rateLogSd = NumericVector::create(rateLogSdScalar) and
// classRate  = NumericVector::create(1.0) reproduces cpp_log_likelihood
// to within fp tolerance (the §7b contract).
//
// Layer 1 only handles brlens-LINKED partitions (edgeLen is one vector
// shared across all classes); Layer 2 widens to nEdge × nClasses for the
// "brlens" unlink token (T3 treatment).
//
// etaNeo is the geometric-mean asymmetry parameter (plan §5.2) that
// replaces rate_neo in the partitioned path. It is only consulted when
// data.hasNeo; at etaNeo = 1 (default) the per-partition call matches
// the legacy rateNeo = 1 path bit-for-bit.
double cpp_log_likelihood_partitioned(
    const McmcData& data,
    Rcpp::IntegerVector parent,
    Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss,
    Rcpp::NumericVector rateLogSd,    // length 1 (linked) or data.nClasses
    Rcpp::NumericVector classRate,    // length 1 (linked) or data.nClasses
    double etaNeo,
    double betaScale = 1.0,
    ClWorkspace* ws = nullptr);


// Per-partition log-likelihood for partial recomputation (M-064).
// M-063: optional ClWorkspace* eliminates per-call heap allocation.
double cpp_partition_log_likelihood(
    const McmcData& data, int partIdx,
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    const Rcpp::IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    double betaScale = 1.0,
    ClWorkspace* ws = nullptr);

// ACRV rate computation (exposed for Gibbs kPrime sweep).
Rcpp::NumericVector cpp_acrv_rates(double rateLogSd, int nCat,
                                    const std::vector<double>& acrvZ);

// Per-site batched pruning for Gibbs kPrime sweep (M-155).
// Fill siteLL[0..nChar-1] with per-character log(avg_lik).
void pruning_jc_acrv_persite(
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edge_length, Rcpp::IntegerMatrix tip_states,
    int kStates,
    Rcpp::NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride,
    double* siteLL);

void pruning_f81_het_acrv_persite(
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edge_length, Rcpp::IntegerMatrix tip_states,
    int kStates, double baseRL,
    const double* betaBins, int nBetaCat,
    Rcpp::NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride,
    double* siteLL);

// Lumped-state variant of pruning_jc_acrv_persite for kObs < kFull.
// Uses kEff = kObs + 1 columns per character (observed + 1 lumped) with the
// JC(k) lumpability identity; result is mathematically identical to the
// uncollapsed kernel, with inner-loop arithmetic scaling as kEff/kFull.
// Caller must ensure stride >= nChar * (kObs + 1) and dispatch only when
// the saving is real (kObs + 1 < kFull).
void pruning_jc_acrv_persite_collapsed(
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edge_length, Rcpp::IntegerMatrix tip_states,
    int kFull, int kObs,
    Rcpp::NumericVector rate_multipliers,
    double* buf, uint8_t* initFlg, int stride,
    double* siteLL);


// Constant-site probability for a given kStates (cache helper for Gibbs sweep).
double const_site_prob_for_k(
    const McmcData& data,
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen,
    int kStates, double betaScale,
    const Rcpp::NumericVector& acrvRates);

#endif  // MKPRIME_MCMC_STATE_H
