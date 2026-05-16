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

  // M-172: pattern-compressed tip states for Gibbs sweep deduplication.
  // uniqueTipStates holds one column per unique tip-pattern; patternIndex[c]
  // is the 0-based index into uniqueTipStates for the c-th local character.
  // All characters sharing a pattern have identical likelihoods under any
  // (k, tree, params) under the symmetric JC model.
  Rcpp::IntegerMatrix uniqueTipStates; // nTip x nUniquePatterns
  Rcpp::IntegerVector patternIndex;    // length nChar, 0-based
  int nUniquePatterns = 0;
};

// Ecology partition for the ecology-aware NT model.
// Tip ecologies are integer-coded, contiguous in [0, kEcology).
// Phase 1: declarations only — likelihood/proposal wiring lands in Phase 2/3.
struct EcologyInfo {
  int kEcology = 0;                      // number of ecology states (0 = disabled)
  Rcpp::IntegerVector tipStates;         // length nTip, 0-indexed
  int refEcology = -1;                   // v2: 0-indexed; -1 if disabled
  Rcpp::NumericVector ecoEdgeMass;       // v2: length kEcology, edge mass per ecology
  bool active() const { return kEcology >= 2 && tipStates.size() > 0; }
};

// Per-iteration state of the ecology-aware extension.
// phi:     magnitude factor. Length 1 (global) or kEcology (per_ecology).
// pi0:     spike-and-slab sparsity (prob of z == none).
// theta:   v2: slab balance per non-reference ecology; length kEcology - 1.
//          P(z = enc) = (1-pi0) * theta_e; P(z = disc) = (1-pi0) * (1 - theta_e).
// zMatrix: v2: nChar x (kEcology - 1), one column per non-reference ecology in
//          ascending ecology-state order (skipping refEcology); values in
//          {0=none, 1=encouraged, 2=discouraged}.
// wEdge:   nEdge x kEcology cache of edge-ecology marginal weights, recomputed
//          once per likelihood evaluation.  Reference column carries weight as
//          usual; rate factor on reference-ecology edges is fixed at 1.
struct EcologyState {
  Rcpp::NumericVector phi;
  double pi0 = 0.0;
  Rcpp::NumericVector theta;   // v2: length kEcology - 1
  Rcpp::IntegerMatrix zMatrix;
  Rcpp::NumericMatrix wEdge;
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

  // Ecology-aware NT model (Phase 1: scaffolding only).
  // ecologyAware mirrors model$ecologyAware; ecology holds the per-tip
  // ecology codes and kEcology. magnitudeMode: 0 = global, 1 = per_ecology.
  // rho0Alpha/Beta parameterise the Beta prior on pi0; sigmaPhi is the
  // LogNormal sdlog for phi. gibbsZEvery is the sweep period over z.
  bool ecologyAware = false;
  EcologyInfo ecology;
  int    magnitudeMode = 0;
  // R5-1: ESS raised from 100 to 480 to match sim3 z-cell count (nChar=480).
  // R5-2: sigmaPhi raised from 1.0 to 1.5 so phi=4 is within 5x of phi=1.
  double rho0Alpha = 360.0;  // Beta(360, 120): mode 0.75, ESS 480
  double rho0Beta  = 120.0;
  double sigmaPhi  = 1.5;
  int    gibbsZEvery = 50;
  // v2: Beta hyperprior on theta_e (slab balance). R5-3: raised from (1,1)
  // uniform to (2,2) to resist theta drifting to the 0.5 saddle where
  // z=1 and z=2 become exchangeable and accelerate pi0 collapse.
  double thetaAlpha = 2.0;
  double thetaBeta  = 2.0;
  // v2: mapping from ecology state s -> column j in zMatrix/theta.
  // Length kEcology; -1 at the reference slot, ascending non-ref index
  // otherwise. Populated by prepare_mcmc_data (and refined by
  // edge-mass refEcology selection in mcmc_ecology).
  std::vector<int> ecologyToZCol;
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

// Single-character JC log-likelihood for Gibbs kPrime sweep.
// Handles JC + ACRV + Het + ascertainment correction + relabeling.
// constSiteProb should be pre-computed via const_site_prob_for_k().
double single_char_loglik_jc(
    const McmcData& data,
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen,
    const int* tipCol,
    int kStates, int kObs,
    double betaScale,
    const Rcpp::NumericVector& acrvRates,
    double constSiteProb);

// Constant-site probability for a given kStates (cache helper for Gibbs sweep).
double const_site_prob_for_k(
    const McmcData& data,
    const Rcpp::IntegerVector& parent,
    const Rcpp::IntegerVector& child,
    const Rcpp::NumericVector& edgeLen,
    int kStates, double betaScale,
    const Rcpp::NumericVector& acrvRates);

#endif  // MKPRIME_MCMC_STATE_H
