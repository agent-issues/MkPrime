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
  // Pre-computed log P_emp(k) for k = 2, ..., 1 + empLogBody.size();
  // entry idx = k - 2.  The tail extends k >= empTailStartK with log mass
  //   empLogTailStartP + (k - empTailStartK) * log(empTailDecay).
  std::vector<double> empLogBody;
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

  // Likelihood mode (v1 marginal-k landing, plan §2). When false, k'_i is
  // a sampled MCMC state variable (legacy sampled-k path). When true,
  // k'_i is analytically marginalised out of the likelihood at every
  // evaluation; the chain only carries (tree, mu, sigma, p). v1 supports
  // the geometric arm only — prepare_mcmc_data enforces this.
  bool marginalK = false;

  // Prior parameterisation for the geometric k' prior under marginal-k.
  //   false (default) = Model B (conditional): k'_i ~ kObs_i + Geometric(p),
  //                     marginal weight p (1-p)^(k - kObs_i).
  //   true            = Model A (unconditional): k'_i ~ 2 + Geometric(p),
  //                     marginal weight p (1-p)^(k - 2); the weight exponent
  //                     base shifts from kObs to 2. The sum is capped at
  //                     k <= kprimeTruncK and renormalised by log Z(p)
  //                     (MARGINAL-K-TRUNC-001). Model A and Model B differ by
  //                     a per-character factor (1-p)^(kObs_i - 2).
  bool unconditionalPrior = false;

  // MARGINAL-K-TRUNC-001: declared truncation cap K on k' for the geometric
  // prior under BOTH likelihoodModes (sampled_k truncates the prior; marginal_k
  // truncates the marginal sum). The prior is a truncated geometric on
  // k' in [2, K], renormalised by log Z(p) (Model A: Z = 1-(1-p)^(K-1), shared;
  // Model B: Z_i = 1-(1-p)^(K-kObs_i+1), per kObs). MUST equal the SBC forward's
  // K_MAX_PRIOR for calibration. Setup must enforce K >= max(kObs) (else a
  // character has empty support -> -Inf).
  //
  // The compile-time default MATCHES the MkPrimeModel default (200L) so that a
  // dataPtr built directly via prepare_mcmc_data (bypassing .InitMcmcData, which
  // wires the model's K via set_kprime_trunc_k) still agrees with the R-side
  // LogPrior. SBC and the truncation tests pin K=30 explicitly via the model.
  int kprimeTruncK = 200;
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


// The parameters partition `partIdx` is scored under by the partition-API
// likelihood (cpp_log_likelihood_partitioned): its class's ACRV shape and
// rate multiplier, with rate_neo fixed at 1. An evaluator that takes them from
// anywhere else scores a different model from the one state->logLik holds.
struct PartEvalParams {
  double rateLogSd;
  double classRate;  // multiplies every edge length
  double rateNeo;
};

static inline PartEvalParams partitioned_eval_params(
    const McmcData& data, int partIdx,
    const Rcpp::NumericVector& rateLogSd,
    const Rcpp::NumericVector& classRate) {
  const int ci = data.parts[partIdx].classIdx - 1;
  PartEvalParams pe;
  pe.rateLogSd = (rateLogSd.size() == 1) ? rateLogSd[0] : rateLogSd[ci];
  pe.classRate = (classRate.size() == 1) ? classRate[0] : classRate[ci];
  pe.rateNeo   = 1.0;
  return pe;
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


// ---------------------------------------------------------------------------
// Case-25 phase-1 helper: batched per-(char, k') log-weight precomputation.
//
// Shared across the sampled-k Gibbs path (case 25 in mcmc.cpp) and — once
// PR-B lands — the marginal-k evaluator in mcmc_likelihood.cpp. The helper
// owns the M-155 / M-164 / M-172 batching, dedup, JC-lumpability collapse,
// and prior-ceiling early termination that case 25 used to do inline.
//
// Output `charLogW` stores `β · LL(k = kObs_i + ko) + logPrior_k` flat in
// (ti, ko) with stride `kMaxKprimeCand`. `charNCand[ti]` is the number of
// evaluated ko slots for char ti before early termination. `charMaxLogW[ti]`
// is the running max of charLogW[ti, *] (consumed by phase-2 categorical
// sampling). Phase-2 scratch (`charMaxLL`, `terminated`, the per-partition
// active-pattern bookkeeping) is internal to the helper.
//
// PR-B note: a sibling helper will reuse this infrastructure but return raw
// per-(char, ko) log-likelihoods (no prior, no β) so the marginal evaluator
// can do its own logSumExp with arbitrary `P(u | hyperparams)` weights.
// ---------------------------------------------------------------------------

// Stage 1b (MARGINAL-K-TRUNC-001 coupling): the marginal numerator sums
// candidates ko with k' = kObs_i + ko, capped at min(nCand, K - kObs_i + 1)
// where K = McmcData.kprimeTruncK. For the sum to actually reach K the
// candidate cap must satisfy kMaxKprimeCand >= K - 1 (worst case kObs = 2);
// otherwise the numerator caps below K while Z_A renormalises [2, K],
// reintroducing the truncation bias. set_kprime_trunc_k() enforces
// K <= kMaxKprimeCand at the boundary. Raised 50 -> 256 to admit the
// real-data default K = 200 with margin. Cost: charLogW/charLLCache grow to
// nTrans * 256 doubles (~0.6 MB at 300 chars) and wBuf to 256 doubles (2 KB
// stack); the heavy per-(node, k') PerKpClCache is dormant in v1 (forward
// decl only), so this does NOT scale that allocation. Any K <= 30 result is
// bit-identical to the old cap (nEff = min(nCand, K - kObs + 1) is unchanged).
constexpr int    kMaxKprimeCand   = 256;    // absolute cap on ko candidates
constexpr double kKprimeLogCutoff = -25.0;  // M-164 prior-ceiling cutoff

struct McmcState;  // defined in mcmc.cpp; forward-declared so the helper
                   // declaration here only depends on pointer semantics.

struct KprimeCharWeights {
  // Flat, length nTrans * kMaxKprimeCand. logW[ti * kMaxKprimeCand + ko]
  // holds β · LL_ko + logPrior_ko (or R_NegInf for un-evaluated slots).
  std::vector<double> charLogW;
  // Number of evaluated ko slots per char (≤ kMaxKprimeCand).
  std::vector<int>    charNCand;
  // Running max of charLogW[ti, *]; phase-2 subtracts this before exp().
  std::vector<double> charMaxLogW;

  void resize(int nTrans) {
    charLogW.assign(static_cast<size_t>(nTrans) * kMaxKprimeCand, R_NegInf);
    charNCand.assign(nTrans, 0);
    charMaxLogW.assign(nTrans, R_NegInf);
  }
};

// Populate `out` with per-(char, ko) Gibbs weights for all transformational
// characters. Pure (no RNG); see `mcmc.cpp:gibbs_kprime_sweep_impl` for the
// canonical caller and the phase-2 categorical sampler that consumes `out`.
//
// NOTE: state->gibbsWs is grown if needed — that's workspace, not logical
// state, so callers can treat this as state-non-modifying. Topology is taken
// from the caller-supplied parent/child/edgeLen (MARGINAL-K-FREEZE-003: a
// proposed tree must be evaluable before it is committed to state); reads
// state->kPrime / p / betaScale; never writes them.
void compute_per_kprime_log_lik(
    McmcData* data, McmcState* state, double beta,
    Rcpp::IntegerVector parent, Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    Rcpp::NumericVector acrvRates,
    KprimeCharWeights& out);

// ---------------------------------------------------------------------------
// Marginal-k evaluator (v1: geometric arm only).
//
// Sums out per-character k'_i analytically:
//   L_marg(y_i | tree, mu, p) =
//     logSumExp_{u=0..u_max} [ LL(y_i | tree, mu, kObs_i + u)
//                              + log P(u | p) ]
// For neomorphic partitions, falls through to cpp_partition_log_likelihood
// (those are not marginalised — k is fixed).
//
// Reads data->parts, data->kObs, data->transIdxGlobal, data->codingType,
// data->relabel, the passed parent/child/edgeLen (NOT state's topology —
// MARGINAL-K-FREEZE-003), plus state->p/rateLoss/rateLogSd/rateNeo/
// betaScale. Does NOT modify state (apart from state->gibbsWs allocation
// as noted on compute_per_kprime_log_lik).
//
// REQUIREMENTS (enforced in MkPrimeModel.R):
//   - data->marginalK == true
//   - data->kPriorLogseries / kPriorBetaGeometric / kPriorEmpiricalGeometric
//     all false (geometric arm only in v1)
//   - data->qHeterogeneity == false
//
// `state` is non-const because the call mutates `state->gibbsWs` (workspace
// allocation) inside `compute_per_kprime_log_lik`, and may populate
// `state->charLLCache` if the marginal-k charLL cache is active.
//
// fillCharLLCache (MARGINAL-K-FREEZE-003 follow-up): when false, this is a
// SCRATCH eval — the per-(char,k') charLLCache is neither READ (forced cold,
// raw LLs always recomputed) nor WRITTEN (no fill, charLLCacheReady left
// untouched). Multi-config moves (weighted_*/block_gibbs_branch, which evaluate
// many topology/branch configs per call) MUST pass false: the cache is valid
// only across a pure-p change, so reusing it across configs reads stale LLs.
// Scratch evals also leave the cache exactly as do_move_impl's entry
// invalidation left it (cold), so no rejected/self-accepting move can poison a
// later mh_logit_p. (Tier-2 perKpCl cache is dormant — never read — so tier-1
// is the only persistent cache the cold path touches.)
double cpp_log_likelihood_marginal(
    McmcData& data, McmcState& state,
    Rcpp::IntegerVector parent,
    Rcpp::IntegerVector child,
    Rcpp::NumericVector edgeLen,
    double rateLoss,
    double rateLogSd,
    double rateNeo,
    ClWorkspace* ws = nullptr,
    bool fillCharLLCache = true);


// ---------------------------------------------------------------------------
// Marginal-k Option A cache (FU-3, plan §5).
//
// Two-tier cache for the marginal-k evaluator:
//
//   Tier 1 — per-(char, ko) `charLLCache` (Option A', shipped in PR-B):
//     Stores the raw per-(transformational char, k-offset) log-likelihood
//     after ascertainment + relabelling corrections. Invalidated by every
//     move except case 30 (mh_logit_p) — a p-move only changes the
//     `log P(u | p)` weights consumed by the per-character logSumExp.
//     Lives on `McmcState::charLLCache` (vector<double>) +
//     `McmcState::charLLNCand` (vector<int>) +
//     `McmcState::charLLCacheReady` (bool).
//
//   Tier 2 — per-(partition, node, ko) Felsenstein partial CL cache
//     (Option A, FU-3 structural landing): A sibling cache that stores
//     the Felsenstein conditional likelihoods at every internal node for
//     every candidate k value (kStates = kObs_partition + ko). Used by
//     the marginal evaluator to skip the full per-ko downpass when only
//     a subtree is dirty (NNI / beta_simplex / Dirichlet / SPR). Lives
//     on `McmcState::perKpCl` (see PerKpClCache below).
//
// Invariance under p-moves (proof: dev/red-team/proofs/marginal-k-
// geometric.md §5):
//   L(y_i | tree, μ, k) does NOT depend on p. Therefore neither Tier 1
//   nor Tier 2 are invalidated by case 30 (mh_logit_p). The per-character
//   `log P(u | p)` weights are recomputed on each p-move from scratch
//   (cheap: nTrans × kMaxKprimeCand scalar ops, no pruning).
//
// Invalidation rules (Tier 2 currently mirrors Tier 1):
//   | Move                        | Tier 1 (charLL)  | Tier 2 (per-(node,k)) |
//   | --------------------------- | ---------------- | --------------------- |
//   | case 30 (mh_logit_p)        | valid            | valid                 |
//   | tree topology (5, 6, ...)   | invalidate_all   | invalidate_all        |
//   | branch length (0, 4, 23..)  | invalidate_all   | invalidate_all        |
//   | rate_log_sd (2, 31, 34)     | invalidate_all   | invalidate_all        |
//   | rate_loss (1)               | invalidate_all   | invalidate_neo_only   |
//   | rate_neo (3)                | invalidate_all   | invalidate_all        |
//
// Memory ceiling: Lazy allocation per ko slice. Worst case at 150 tips,
// 300 chars, kMaxKprimeCand = 50: ~600 MB. PerKpClCache::ensure_capacity
// enforces a 4 GB cap (plan §12 fail-loud) BEFORE allocation, throwing
// with a clear message that names the offending slot count and a tuning
// hint (lower kMaxKprimeCand at the build level, or reduce nTrans).
//
// Status (this PR — FU-3 structural landing):
//   PerKpClCache fields are reserved on McmcState; invalidation hooks are
//   wired into do_move_impl in lockstep with charLLCacheReady. The cache
//   is NOT yet populated/consumed by the marginal evaluator (Tier 1 still
//   handles the only runtime win — p-move acceleration). Wiring Tier 2
//   into the per-move partial-CL dispatchers (lines 5083-5256 in mcmc.cpp)
//   is FU-3b. The structure here defines the contract that FU-3b must
//   honour, locked in by tests in
//   tests/testthat/test-marginal-k-cache-option-a.R.
// ---------------------------------------------------------------------------

// Forward declaration: defined in mcmc.cpp alongside McmcState.
struct PerKpClCache;

// Memory ceiling for the Tier 2 cache (plan §12). Exceeding this throws
// with a clear message at allocation time rather than OOMing.
constexpr size_t kMarginalKCacheMaxBytes =
    static_cast<size_t>(4) * 1024 * 1024 * 1024;  // 4 GB

#endif  // MKPRIME_MCMC_STATE_H
