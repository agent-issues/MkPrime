// C++ MCMC engine — state management, prior, and propose/accept cycle.
//
// McmcState (mutable per chain) holds tree + parameters.
// do_move_cpp() performs a full MH step, modifying state in-place.
// R creates states via init_mcmc_state(), reads via get_mcmc_state().
//
// M-065: Eliminated IntegerMatrix round-trips. do_move_impl now calls
// *_proposal_impl (parent/child vectors) and cpp_log_likelihood (vectors)
// directly, avoiding repeated matrix construction/decomposition.

#include "mcmc_state.h"
#include "gibbs_partial_cl.h"
#include "fitch.h"
#include "node_cl_cache.h"
#include <TreeTools/edge_to_splits.h>
#include <TreeTools/renumber_tree.h>
#include <cmath>
#include <cstring>
#include <chrono>
#include <cstdio>
#include <algorithm>

using namespace Rcpp;

// Forward declaration for relabelling correction (corrections.cpp)
double mk_prime_relabel_log(int kPrime, int kObs);

// Forward declarations for proposals in other TUs
// M-065: vector-based _impl versions (no edge matrix)
List nni_proposal_impl(IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
List spr_proposal_impl(IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
// M-084: subtree swap helpers (tree_moves.cpp)
List swap_subtrees_impl(IntegerVector parent, IntegerVector child,
                        int nTip, double treeLength,
                        NumericVector relBrLengths, int nodeA, int nodeB);
std::vector<int> get_valid_swap_partners_impl(const IntegerVector& parent,
                                              const IntegerVector& child,
                                              int nTip, int pruneNode);
// M-053: TBR proposal (tree_moves.cpp)
List tbr_proposal_impl(IntegerVector parent, IntegerVector child,
                        int nTip, double treeLength,
                        NumericVector relBrLengths);
List beta_simplex_proposal(NumericVector x, int index, double tuning);
bool beta_simplex_impl(NumericVector& x, int index, double tuning,
                       double& logHastings, int& outOther,
                       double& outOldIdx, double& outOldOther);  // OPP-5

// M-125: block Dirichlet simplex proposal (defined in proposals.cpp)
bool dirichlet_simplex_impl(NumericVector& x, int nCats, double alpha,
                            double& logHastings, NumericVector& snapshot,
                            std::vector<int>& modifiedEdges);
// M-127: localized Dirichlet simplex (neighborhood selection)
bool local_dirichlet_impl(NumericVector& x,
                          const IntegerVector& parent,
                          const IntegerVector& child,
                          int nCats, double alpha,
                          double& logHastings, NumericVector& snapshot,
                          std::vector<int>& modifiedEdges);


// ---------------------------------------------------------------------------
// Bactrian perturbation kernel  (M-118, Yang & Rodríguez 2013)
//
// Bimodal mixture 0.5*N(-m, 1-m²) + 0.5*N(+m, 1-m²) with m = 0.95.
// Replaces Uniform(-0.5, 0.5) in scale proposals.  The distribution is
// symmetric about zero, so the Hastings ratio for scale moves is unchanged
// (log(mult)).  Avoids near-zero perturbations, improving ESS/iter by
// ~15-50% for scalar parameters at zero computational overhead.
// ---------------------------------------------------------------------------
static constexpr double BACTRIAN_M  = 0.95;
static const     double BACTRIAN_SD = std::sqrt(1.0 - BACTRIAN_M * BACTRIAN_M);
// Scale so overall sd matches Uniform(-0.5, 0.5), i.e. 1/√12.
// The benefit is bimodal *shape* (avoids near-zero), not larger variance.
static const     double BACTRIAN_SCALE = 1.0 / std::sqrt(12.0);

static inline double bactrian_perturbation() {
  double z = R::rnorm(0.0, BACTRIAN_SD);
  double raw = (R::unif_rand() < 0.5) ? (BACTRIAN_M + z) : (-BACTRIAN_M + z);
  return raw * BACTRIAN_SCALE;
}

// Exported for unit testing (test-bactrian.R)
// [[Rcpp::export]]
NumericVector bactrian_draws(int n) {
  NumericVector out(n);
  for (int i = 0; i < n; ++i)
    out[i] = bactrian_perturbation();
  return out;
}


// ---------------------------------------------------------------------------
// M-120: 2D correlated Bactrian kernel for joint proposals.
// Both components share the same mode (±M) with correlated Gaussian noise.
// Correlation ρ is learned during warmup from posterior sample correlations.
// ---------------------------------------------------------------------------
static inline void bactrian_2d_perturbation(double rho,
                                            double& z1, double& z2) {
  // Correlated Gaussian noise
  double n1 = R::rnorm(0.0, 1.0);
  double n2 = R::rnorm(0.0, 1.0);
  double e1 = BACTRIAN_SD * n1;
  double e2 = BACTRIAN_SD * (rho * n1 + std::sqrt(1.0 - rho * rho) * n2);
  // Mode coupling: p_same = (1+rho)/2 gives Cor(z1,z2) = rho exactly.
  // Each marginal stays standard 1D Bactrian regardless of coupling.
  double s1 = (R::unif_rand() < 0.5) ? BACTRIAN_M : -BACTRIAN_M;
  double pSame = 0.5 * (1.0 + rho);
  double s2 = (R::unif_rand() < pSame) ? s1 : -s1;
  z1 = (s1 + e1) * BACTRIAN_SCALE;
  z2 = (s2 + e2) * BACTRIAN_SCALE;
}

// Exported for unit testing (test-joint-2d.R)
// [[Rcpp::export]]
NumericMatrix bactrian_2d_draws(int n, double rho) {
  NumericMatrix out(n, 2);
  for (int i = 0; i < n; ++i) {
    double z1, z2;
    bactrian_2d_perturbation(rho, z1, z2);
    out(i, 0) = z1;
    out(i, 1) = z2;
  }
  return out;
}


// Exported for unit testing (test-pspr.R)
// [[Rcpp::export]]
int fitch_score_r(IntegerVector parent, IntegerVector child,
                  IntegerMatrix tipStates, int nTip, int kStates) {
  std::vector<std::pair<IntegerMatrix, int>> parts = {{tipStates, kStates}};
  return fitch_score_all(INTEGER(parent), INTEGER(child),
                         parent.size(), nTip, parts);
}


// Exported for unit testing.  Row arguments are 1-based edge indices.
// [[Rcpp::export]]
IntegerVector fitch_score_candidates_r(
    IntegerVector parent, IntegerVector child,
    IntegerMatrix tipStates, int nTip, int kStates,
    int pruneRow, int parentRow, int sibRow,
    int u, int v, int sibNode, IntegerVector candidates) {
  const int nEdge = parent.size();
  // Every row argument indexes a raw int* inside fitch_score_candidates, so an
  // out-of-range one would write past the vector rather than fail.
  auto requireRow = [nEdge](int row, const char* name) {
    if (row < 1 || row > nEdge) Rcpp::stop("%s out of range", name);
  };
  if (child.size() != nEdge)
    Rcpp::stop("parent and child must have equal length");
  requireRow(pruneRow, "pruneRow");
  requireRow(parentRow, "parentRow");
  requireRow(sibRow, "sibRow");
  for (int i = 0; i < candidates.size(); ++i) requireRow(candidates[i], "candidates");
  if (nTip < 1 || nTip > tipStates.nrow())
    Rcpp::stop("nTip exceeds the rows of tipStates");
  if (kStates < 1 || kStates > 26)
    Rcpp::stop("kStates must lie in 1:26");
  std::vector<std::pair<IntegerMatrix, int>> parts = {{tipStates, kStates}};
  std::vector<int> cand(candidates.size());
  for (int i = 0; i < candidates.size(); ++i) cand[i] = candidates[i] - 1;
  IntegerVector wp = clone(parent), wc = clone(child);
  std::vector<int> scores;
  fitch_score_candidates(INTEGER(wp), INTEGER(wc), nEdge, nTip, parts,
                         pruneRow - 1, parentRow - 1, sibRow - 1,
                         u, v, sibNode, cand, scores);
  return wrap(scores);
}


// ---------------------------------------------------------------------------
// McmcState: mutable per-chain state
// ---------------------------------------------------------------------------

struct McmcState {
  IntegerVector parent;
  IntegerVector child;
  NumericVector relBrLengths;
  double treeLength;
  double rateLoss;
  double rateLogSd;
  double rateNeo;
  double p;
  double kprimeAlpha;   // Beta-Geometric hyperparameter α
  double kprimeBeta;    // Beta-Geometric hyperparameter β
  IntegerVector kPrime;
  double logLik;
  double logPrior;
  // Per-partition log-likelihood cache (M-064)
  std::vector<double> partLogLik;
  // Pre-allocated CL workspace (M-063): eliminates per-call heap allocations
  ClWorkspace clWs;
  // M-052: Q-matrix heterogeneity — Dirichlet-marginal beta_scale parameter.
  double betaScale = 1.0;
  // M-155: dedicated workspace for Gibbs kPrime batched pruning
  ClWorkspace gibbsWs;
  // M-121: persistent node-level CL cache for partial evaluation
  NodeCLCache nodeCL;
  // M-125: snapshot for block Dirichlet branch-length rollback
  NumericVector brSnapshot;
  // M-127: which edges the Dirichlet proposal modified (for partial CL eval)
  std::vector<int> dirEdges;
  // DIAG counters
  int diagDirPartialCount = 0;
  int diagDirFullbackCount = 0;
  int diagNniPartialCount = 0;
  int diagBsPartialCount = 0;
  int diagCachePopCount = 0;
  // Whether the next move selection uses cache-boosted weights. Set from the
  // type of the move just proposed, never from its acceptance or from
  // nodeCL.ready(): either would make the mixture weights depend on the state,
  // and the chain would no longer leave the posterior invariant.
  bool cacheBoostNext = false;

  // Partition-API (Layer 1, plan v4 §4.2).
  // When usePartitioned is false (legacy path), these are ignored and the
  // legacy scalar fields (rateLogSd, rateNeo) govern the MCMC loop exactly
  // as before. When usePartitioned is true, cpp_log_likelihood_partitioned
  // is called with the per-class vectors below.
  //
  // Invariant: classRateLogSd[0] == rateLogSd and etaNeo == 1.0 (frozen in
  // Layer 1) so that rateNeo == 1.0. They are kept in lockstep so the partial-
  // CL paths (which take scalars) can still use state->rateLogSd safely.
  bool usePartitioned = false;
  NumericVector classRateLogSd;   // length 1 (linked) or nClasses (unlinked shape)
  NumericVector classW;           // length nClasses (simplex; 1 when trivial)
  NumericVector classRate;        // length 1 (linked) or nClasses; derived from w
  IntegerVector nCharPerClass;    // length nClasses
  double etaNeo = 1.0;           // frozen at 1.0 in Layer 1

  // Pooled half-normal hyperprior on per-class σ_c (plan: hyperprior task).
  // Active iff useHyperpriorOnSigma is true. Under the non-centred
  // parameterisation σ_c = hyperTau * classZ[c]; classZ holds the sampled
  // z_c values, hyperTau holds the sampled population scale.
  // classRateLogSd is the derived σ_c that the likelihood consumes; it is
  // kept in sync after every move that touches classZ or hyperTau.
  //
  // When useHyperpriorOnSigma is false, classZ is empty and hyperTau is
  // unused; classRateLogSd is sampled directly via the legacy Gamma path.
  bool          useHyperpriorOnSigma = false;
  double        hyperTau = 1.0;
  NumericVector classZ;           // length classRateLogSd.size() when active

  // -----------------------------------------------------------------------
  // Marginal-k (v1) charLL cache (plan §5, sub-option A').
  //
  // When data->marginalK is true, the marginal evaluator stores per-(char,
  // ko) raw log-likelihoods (no prior, no β) flat in `charLLCache` with
  // stride kMaxKprimeCand. `charLLNCand[ti]` is the number of valid ko
  // slots for char ti (after early termination). `charLLCacheReady` is
  // true iff the cache is consistent with the current
  // (parent, child, edgeLen, rateLoss, rateLogSd, rateNeo) — that is, the
  // cache is invalidated on every move except case 30 (mh_logit_p), since
  // a p-move only changes the per-character weights `log P(u | p)`, not
  // the raw per-(char, k') LLs.
  //
  // The cache lets the marginal evaluator skip the helper call inside
  // case 30, recomputing only the per-character logSumExp against the new
  // weights — O(nTrans × u_max) ops, no pruning.
  std::vector<double> charLLCache;
  std::vector<int>    charLLNCand;
  bool                charLLCacheReady = false;

  // -----------------------------------------------------------------------
  // Marginal-k Tier 2 cache (Option A — per-(partition, node, ko) partial
  // Felsenstein CLs). FU-3 structural landing — see mcmc_state.h §Marginal-
  // k Option A cache for the full contract.
  //
  // Field layout matches Tier 1's invalidation rhythm:
  //   - perKpClSlots[ko] (when allocated) holds a NodeCLCache-style unit
  //     array for that k-offset. Lazy allocated on first marginal eval
  //     that touches that slot.
  //   - perKpClReady[ko] mirrors NodeCLCache::ready() for that slot.
  //   - perKpClReadyAny is the union — true iff at least one slot has been
  //     populated and is still valid.
  //
  // Invalidation rules: any non-p move (moveType != 30) sets
  // perKpClReadyAny = false and clears all per-slot ready flags. This
  // matches charLLCacheReady's invalidation contract (see do_move_impl).
  // Per-character partial-CL paths (NNI, beta_simplex, Dirichlet, SPR) in
  // the marginal-k branch would mutate these slots in-place — wiring is
  // FU-3b; for now this PR lands the structure + contract + tests.
  //
  // Memory: enforced via PerKpClCache::ensure_capacity (4 GB ceiling per
  // kMarginalKCacheMaxBytes in mcmc_state.h).
  std::vector<NodeCLCache> perKpClSlots;
  std::vector<bool>        perKpClReady;
  bool                     perKpClReadyAny = false;

  // Reset all Tier 2 slots to invalid. Cheap O(K_MAX_CAND) flag-flip; the
  // CL buffers themselves are not freed — lazy reuse keeps the next
  // populate path allocation-free when the same ko slots come back.
  void invalidate_per_kp_cl_all() {
    for (auto& slot : perKpClSlots) slot.invalidate_all();
    std::fill(perKpClReady.begin(), perKpClReady.end(), false);
    perKpClReadyAny = false;
  }
};


// ---------------------------------------------------------------------------
// Log prior (mirrors LogPrior in MkPrimeModel.R)
// ---------------------------------------------------------------------------

static double cpp_log_prior(
    const McmcData& data,
    double treeLength, const NumericVector& relBrLengths,
    double rateLoss, double rateLogSd, double rateNeo,
    double p, const IntegerVector& kPrime,
    double betaScale = 1.0,
    double kprimeAlpha = 1.0, double kprimeBeta = 1.0) {

  // Hard floor: prevents Mk_v singularity (corrected likelihood → +∞ at zero)
  if (treeLength < 1e-6) return R_NegInf;
  if (rateLogSd < 0.0)   return R_NegInf;
  for (int i = 0; i < relBrLengths.size(); ++i) {
    if (relBrLengths[i] <= 0.0) return R_NegInf;
  }
  if (data.hasNeo) {
    if (rateLoss <= 0.0) return R_NegInf;
    if (rateNeo <= 0.0)  return R_NegInf;
  }

  bool hasTrans = (data.transIdxGlobal.size() > 0);
  if (hasTrans) {
    // p boundary check applies to hierarchical geometric and empirical-geometric
    if (!data.kPriorLogseries && !data.kPriorBetaGeometric &&
        (p <= 0.0 || p >= 1.0)) return R_NegInf;
    // kprimeAlpha, kprimeBeta must be positive for beta_geometric
    if (data.kPriorBetaGeometric &&
        (kprimeAlpha <= 0.0 || kprimeBeta <= 0.0)) return R_NegInf;
    // The (plain) geometric arm is truncated at K (k' in [2, K],
    // MARGINAL-K-TRUNC-001): k' > K carries zero prior mass under sampled_k so
    // the joint marginalises to the marginal-k value (RB-consistency). Other
    // arms (logseries / beta_geometric / empirical_geometric) are NOT
    // K-truncated; under marginal_k k'_i is pinned to kObs_i and the evaluator
    // applies the cap itself, so the upper guard is scoped to sampled_k.
    const bool isPlainGeom = !data.kPriorLogseries &&
                             !data.kPriorBetaGeometric &&
                             !data.kPriorEmpiricalGeometric;
    for (int i = 0; i < data.transIdxGlobal.size(); ++i) {
      int gi = data.transIdxGlobal[i];
      if (kPrime[gi] < data.kObs[gi]) return R_NegInf;
      if (isPlainGeom && !data.marginalK && kPrime[gi] > data.kprimeTruncK)
        return R_NegInf;
    }
  }

  double lp = 0.0;

  // Tree length: Gamma(shape, rate)
  lp += R::dgamma(treeLength, data.treeLengthShape,
                   1.0 / data.treeLengthRate, 1);

  // Dirichlet(1,...,1) = log((n-1)!) = lgamma(n)
  lp += std::lgamma(static_cast<double>(relBrLengths.size()));

  if (data.hasNeo) {
    lp += R::dlnorm(rateLoss, data.rateLossMeanlog, data.rateLossSdlog, 1);
    lp += R::dlnorm(rateNeo,  data.rateNeoMeanlog,  data.rateNeoSdlog,  1);
  }

  if (rateLogSd > 0.0) {
    lp += R::dgamma(rateLogSd, data.rateLogSdShape,
                     1.0 / data.rateLogSdRate, 1);
  } else if (data.rateLogSdShape > 1.0) {
    return R_NegInf;
  }

  if (hasTrans) {
    int nTrans = data.transIdxGlobal.size();
    if (data.kPriorLogseries) {
      // Logseries: log P(k'_i; c) = k'_i*log(c) - log(k'_i) - log(-log(1-c))
      double c = data.kprimeLogseriesC;
      double logC   = std::log(c);
      double logNorm = std::log(-std::log1p(-c));  // log(-log(1-c))
      for (int i = 0; i < nTrans; ++i) {
        int kp = kPrime[data.transIdxGlobal[i]];
        lp += kp * logC - std::log(static_cast<double>(kp)) - logNorm;
      }
      // No p / Beta term
    } else if (data.kPriorBetaGeometric) {
      // Per-character Beta-Geometric: P(u | α, β) = B(α+1, β+u) / B(α, β)
      double a = kprimeAlpha;
      double b = kprimeBeta;
      double lbAB = R::lbeta(a, b);
      for (int i = 0; i < nTrans; ++i) {
        int gi = data.transIdxGlobal[i];
        int u = kPrime[gi] - data.kObs[gi];
        lp += R::lbeta(a + 1.0, b + static_cast<double>(u)) - lbAB;
      }
      // Hyperprior: Exp(1) on α and β
      lp += R::dexp(a, 1.0, 1);
      lp += R::dexp(b, 1.0, 1);
    } else if (data.kPriorEmpiricalGeometric) {
      // Convolution prior: k' = N_obs + N_unobs, with N_obs ~ empirical pmf
      // and N_unobs ~ Geometric(p).  For each character,
      //   log P(k'_i = m) = logSumExp_{j=2..m} [ log P_emp(j) + log p
      //                                          + (m - j) * log(1 - p) ]
      // Body of P_emp has explicit log values in data.empLogBody[];
      // beyond it the pmf decays geometrically with log-mass
      // `empLogTailStartP + (k - empTailStartK) * log(empTailDecay)`.
      double logP    = std::log(p);
      double log1mP  = std::log1p(-p);
      double logQ    = (data.empTailDecay > 0.0) ? std::log(data.empTailDecay)
                                                 : R_NegInf;
      int    bodyLen = static_cast<int>(data.empLogBody.size());
      // Scratch buffer for logSumExp (reused per call)
      std::vector<double> terms;
      terms.reserve(64);
      // log P(k' = m | p) = logSumExp_{j=2..m} [ logP_emp(j) + logP + (m-j) log1mP ]
      auto logPconv = [&](int m) -> double {
        terms.clear();
        double mx = R_NegInf;
        for (int j = 2; j <= m; ++j) {
          double logPemp;
          int bodyIdx = j - 2;
          if (bodyIdx < bodyLen) {
            logPemp = data.empLogBody[bodyIdx];
          } else if (data.empTailStartK > 0 && j >= data.empTailStartK &&
                     std::isfinite(data.empLogTailStartP) &&
                     std::isfinite(logQ)) {
            logPemp = data.empLogTailStartP +
                      (j - data.empTailStartK) * logQ;
          } else {
            continue;  // no mass at this j
          }
          if (!std::isfinite(logPemp)) continue;
          double term = logPemp + logP + (m - j) * log1mP;
          terms.push_back(term);
          if (term > mx) mx = term;
        }
        if (terms.empty() || !std::isfinite(mx)) return R_NegInf;
        double sumExp = 0.0;
        for (double t : terms) sumExp += std::exp(t - mx);
        return mx + std::log(sumExp);
      };
      for (int i = 0; i < nTrans; ++i) {
        int gi = data.transIdxGlobal[i];
        int m = kPrime[gi];
        if (m < 2) return R_NegInf;
        double numer = logPconv(m);             // untruncated convolution log-pmf
        if (!std::isfinite(numer)) return R_NegInf;
        lp += numer;
        // EG-001: under the conditional variant (Model B), LogPrior enforces
        // k'_i >= kObs_i, so renormalise over that truncated support.
        //   Z_i(p) = sum_{k>=kObs_i} P(k|p) = 1 - sum_{m'=2}^{kObs_i-1} P(m'|p)
        // (total mass over k>=2 is 1; proofs/kprime-priors.md s4.3). No-op when
        // kObs_i <= 2. Z_i depends on p, so it is NOT absorbed by p-varying MH.
        // Under the unconditional variant (Model A) the prior lives on the full
        // support k' >= 2 with Z_i == 1; the likelihood enforces the k' >= kObs
        // floor. Mirrors the geometric arm's Model A/B split above.
        if (!data.unconditionalPrior) {
          int kobs_i = data.kObs[gi];
          if (kobs_i > 2) {
            double belowMass = 0.0;
            for (int mm = 2; mm < kobs_i; ++mm) {
              double lpmm = logPconv(mm);
              if (std::isfinite(lpmm)) belowMass += std::exp(lpmm);
            }
            double Zi = 1.0 - belowMass;
            if (!(Zi > 0.0) || !std::isfinite(Zi)) return R_NegInf;
            lp -= std::log(Zi);
          }
        }
      }
      // p: Beta hyperprior (same as plain geometric)
      lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
    } else {
      // Hierarchical geometric, TRUNCATED at the declared cap K (k' in [2, K])
      // and renormalised by Z(p) — MARGINAL-K-TRUNC-001. This is the sampled-k
      // counterpart of the marginal-k normaliser (cpp_log_likelihood_marginal,
      // ~lines 4334-4338 and 4369-4373): summing exp(this per-character prior
      // term + that character's LL) over k'_i in [kObs_i, K] reproduces the
      // marginal-k per-character value, so likelihoodMode "sampled_k" and
      // "marginal_k" target the SAME posterior (Rao-Blackwell consistency).
      // The two Model A/B branches and the exact log1p/exp forms below mirror
      // the marginal evaluator so the deterministic logSumExp bit-check matches.
      //
      // Under marginal-k mode the per-character P(k'_i | p) mass (incl. the
      // truncation normaliser) is consumed by cpp_log_likelihood_marginal, so
      // only the hyperprior on p is added here in that mode.
      if (!data.marginalK) {
        const int    K      = data.kprimeTruncK;
        const double logP   = std::log(p);
        const double log1mP = std::log1p(-p);
        if (data.unconditionalPrior) {
          // Model A: P(k'_i = k | p) propto p (1-p)^(k-2), k in [2, K].
          // logZA = log(1 - (1-p)^(K-1)) is shared across all characters.
          const double logZA = std::log1p(-std::exp((K - 1) * log1mP));
          for (int i = 0; i < nTrans; ++i) {
            int gi = data.transIdxGlobal[i];
            lp += logP + (kPrime[gi] - 2) * log1mP - logZA;
          }
        } else {
          // Model B: P(k'_i = kObs_i + u | p) propto p (1-p)^u,
          // u in [0, K - kObs_i]; logZB_i = log(1 - (1-p)^(K - kObs_i + 1)).
          for (int i = 0; i < nTrans; ++i) {
            int gi = data.transIdxGlobal[i];
            int u  = kPrime[gi] - data.kObs[gi];
            lp += logP + u * log1mP
                - std::log1p(-std::exp((K - data.kObs[gi] + 1) * log1mP));
          }
        }
      }
      lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
    }
  }

  // M-052: beta_scale prior — Gamma(shape, rate)
  if (data.qHeterogeneity) {
    if (betaScale <= 0.0) return R_NegInf;
    lp += R::dgamma(betaScale, data.betaScaleShape,
                     1.0 / data.betaScaleRate, 1);
  }

  return lp;
}


// ---------------------------------------------------------------------------
// cpp_log_prior_partitioned — plan v4 §5.1 + §5.4
//
// Adds three per-class prior contributions on top of the legacy log-prior:
//
//   1. Dirichlet(α) on classW (unit simplex) — plan §5.1.
//      α == 1 gives a flat Dirichlet whose log-density is a constant;
//      at K == 1 (trivial spec) the Dirichlet degenerates and its
//      contribution is 0 — matching the legacy scalar path (§7b contract).
//
//   2. Gamma(rateLogSdShape, rateLogSdRate) i.i.d. on each
//      classRateLogSd[c] — plan §5.4. Skipped when classRateLogSd has
//      length 1 (linked shape) because the legacy rateLogSd prior already
//      covers that value — no double-counting.
//
//   3. LogNormal(0, rateNeoSdlog) on etaNeo — §5.2. Included for
//      forward compatibility; in Layer 1 etaNeo is frozen at 1.0 so
//      the contribution is the constant LogNormal-mode density.
//      Skipped entirely when !data.hasNeo.
//
// All other legacy contributions are identical to cpp_log_prior — the
// same scalar state variables are consumed.
// ---------------------------------------------------------------------------

static double cpp_log_prior_partitioned(
    const McmcData& data,
    double treeLength, const NumericVector& relBrLengths,
    double rateLoss, double rateLogSd, double rateNeo,
    double p, const IntegerVector& kPrime,
    const NumericVector& classRateLogSd,  // length 1 (linked) or nClasses
    const NumericVector& classW,          // length 1 or nClasses (simplex)
    double etaNeo,
    bool useHyperpriorOnSigma,
    double hyperTau,
    const NumericVector& classZ,
    double betaScale = 1.0,
    double kprimeAlpha = 1.0, double kprimeBeta = 1.0) {

  // 1. Legacy contributions — delegate to the existing static function.
  //    This covers tree length, Dirichlet(1) on rel branch lengths, rate_loss,
  //    rate_log_sd (the scalar shared across all classes), rate_neo (scalar),
  //    kPrime prior, and beta_scale.
  double lp = cpp_log_prior(
    data, treeLength, relBrLengths, rateLoss, rateLogSd, rateNeo,
    p, kPrime, betaScale, kprimeAlpha, kprimeBeta);

  if (!R_FINITE(lp)) return R_NegInf;

  // 2. Dirichlet(α) on classW — per-class rate contribution (§5.1).
  //    K == 1 (trivial spec): w = (1.0) is the only point on the degenerate
  //    simplex. The Dirichlet log-density is 0 (normalising constant = 0 for
  //    K=1 after the lgamma(1*α)-1*lgamma(α) = 0 identity). Skip the loop.
  int K = classW.size();
  if (K > 1) {
    double alpha = data.classRateConcentration;
    // log Dir(w; α) = lgamma(K*α) - K*lgamma(α) + (α-1)*sum(log(w))
    double logDirConst = std::lgamma(K * alpha) - K * std::lgamma(alpha);
    double sumLogW = 0.0;
    for (int c = 0; c < K; ++c) {
      if (classW[c] <= 0.0) return R_NegInf;
      sumLogW += std::log(classW[c]);
    }
    lp += logDirConst + (alpha - 1.0) * sumLogW;
  }

  // 3. Per-class ACRV-shape prior. Two regimes:
  //    K == 1 (length-1 classRateLogSd): linked or degenerate; the scalar
  //      rateLogSd prior in cpp_log_prior already covers it — nothing extra.
  //    K >= 2: branch on the hyperprior flag.
  //      (a) hyperprior_pooled (default): non-centred half-normal hierarchy
  //          σ_c = τ · z_c, z_c ~ HN(1) i.i.d., τ ~ HN(1). The Gamma on σ_0
  //          added by cpp_log_prior is the WRONG prior for this regime, so
  //          we SUBTRACT it and ADD HN(1) terms on each z_c (c = 0..K-1)
  //          and on τ.
  //      (b) gamma_independent (legacy): each σ_c independently ~ Gamma;
  //          σ_0 already counted, loop adds σ_1..σ_{K-1}.
  if (classRateLogSd.size() > 1) {
    if (useHyperpriorOnSigma) {
      // Subtract legacy Gamma on σ_0 added by cpp_log_prior. Boundary
      // semantics mirror cpp_log_prior (rateLogSd == 0 with shape > 1 has
      // already short-circuited via cpp_log_prior returning -Inf).
      if (rateLogSd > 0.0) {
        lp -= R::dgamma(rateLogSd, data.rateLogSdShape,
                         1.0 / data.rateLogSdRate, 1);
      }
      // HN(1) on every z_c: log p(z; 1) = log(2) + log φ(z; 0, 1) for z>=0.
      // Boundary z == 0 is finite (HN density at 0 equals √(2/π)).
      if (classZ.size() != classRateLogSd.size()) return R_NegInf;
      for (int c = 0; c < classZ.size(); ++c) {
        if (classZ[c] < 0.0) return R_NegInf;
        lp += std::log(2.0) + R::dnorm(classZ[c], 0.0, 1.0, 1);
      }
      // HN(1) on τ.
      if (hyperTau < 0.0) return R_NegInf;
      lp += std::log(2.0) + R::dnorm(hyperTau, 0.0, 1.0, 1);
    } else {
      for (int c = 1; c < classRateLogSd.size(); ++c) {
        double sd_c = classRateLogSd[c];
        if (sd_c < 0.0) return R_NegInf;
        if (sd_c > 0.0) {
          lp += R::dgamma(sd_c, data.rateLogSdShape,
                           1.0 / data.rateLogSdRate, 1);
        } else if (data.rateLogSdShape > 1.0) {
          return R_NegInf;
        }
        // sd_c == 0, shape == 1: dgamma at 0 = rateLogSdRate; the density is
        // finite and is added here. Consistent with the legacy boundary
        // treatment in cpp_log_prior.
      }
    }
  }

  // 4. etaNeo prior: LogNormal(0, rateNeoSdlog) — plan §5.2.
  //    Only when hasNeo; at etaNeo == 1.0 this is the mode of the LogNormal,
  //    but the prior is still finite and correct for the MH ratio.
  //    Note: the legacy rateNeo prior (inside cpp_log_prior) is already added
  //    for the scalar rateNeo. In the partitioned path rateNeo == etaNeo
  //    (Layer 1 sets both to 1.0 and keeps them in lockstep), so the legacy
  //    term already covers etaNeo. No extra term is needed here; this comment
  //    is a forward-compatibility marker for when etaNeo becomes free.

  return lp;
}


// ---------------------------------------------------------------------------
// init_mcmc_state: create XPtr<McmcState>
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP init_mcmc_state(IntegerVector parent, IntegerVector child,
                     NumericVector relBrLengths, double treeLength,
                     double rateLoss, double rateLogSd, double rateNeo,
                     double p, IntegerVector kPrime,
                     double logLik, double logPrior,
                     double betaScale = 1.0,
                     double kprimeAlpha = 1.0,
                     double kprimeBeta = 1.0,
                     NumericVector classRateLogSd = NumericVector(0),
                     NumericVector classW         = NumericVector(0),
                     NumericVector classRate      = NumericVector(0),
                     IntegerVector nCharPerClass  = IntegerVector(0),
                     double etaNeo = 1.0,
                     bool useHyperpriorOnSigma   = false,
                     double hyperTau             = 1.0,
                     NumericVector classZ        = NumericVector(0)) {
  McmcState* s = new McmcState();
  s->parent       = clone(parent);
  s->child        = clone(child);
  s->relBrLengths = clone(relBrLengths);
  s->treeLength   = treeLength;
  s->rateLoss     = rateLoss;
  s->rateLogSd    = rateLogSd;
  s->rateNeo      = rateNeo;
  s->p            = p;
  s->kprimeAlpha  = kprimeAlpha;
  s->kprimeBeta   = kprimeBeta;
  s->kPrime       = clone(kPrime);
  s->logLik       = logLik;
  s->logPrior     = logPrior;
  s->betaScale    = betaScale;
  s->brSnapshot   = NumericVector(relBrLengths.size());
  s->etaNeo       = etaNeo;
  // Partition-API (Layer 1): when per-class vectors are supplied, activate
  // the partitioned likelihood path. Legacy callers pass no per-class args;
  // their behaviour is unchanged (usePartitioned stays false).
  if (classRate.size() > 0) {
    s->usePartitioned  = true;
    s->classRateLogSd  = clone(classRateLogSd);
    s->classW          = clone(classW);
    s->classRate       = clone(classRate);
    s->nCharPerClass   = clone(nCharPerClass);
    s->etaNeo          = etaNeo;
    // Invariant: keep legacy scalar in lockstep with classRateLogSd[0].
    if (classRateLogSd.size() > 0)
      s->rateLogSd = classRateLogSd[0];
    // Layer 1 freezes etaNeo at 1.0 → rateNeo stays 1.0 (no change needed).

    // Hyperprior on per-class σ_c: only activates when there is more than
    // one σ (i.e. shape unlinked with K >= 2). Degenerate at K == 1, in
    // which case the legacy single-σ Gamma prior applies via cpp_log_prior.
    if (useHyperpriorOnSigma && classRateLogSd.size() >= 2) {
      s->useHyperpriorOnSigma = true;
      s->hyperTau             = hyperTau;
      s->classZ               = clone(classZ);
    }
  }
  return Rcpp::XPtr<McmcState>(s, true);
}


// Forward declaration (defined below, M-083) — fill_partition_cache needs it to
// set a marginal-aware initial logLik under marginal_k (MARGINAL-K-INIT-001).
static double compute_full_loglik(const McmcData& data, McmcState& state);

// [[Rcpp::export]]
void fill_partition_cache(SEXP dataPtr, SEXP statePtr) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  int nParts = (int)data->parts.size();
  int nEdge  = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];
  state->partLogLik.resize(nParts);
  double totalLogLik = 0.0;
  for (int pi = 0; pi < nParts; ++pi) {
    state->partLogLik[pi] = cpp_partition_log_likelihood(
      *data, pi, state->parent, state->child, edgeLen,
      state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
      state->betaScale,
      state->clWs.ready() ? &state->clWs : nullptr);
    totalLogLik += state->partLogLik[pi];
  }
  // Sync state->logLik with the C++ partition sum.  The R-computed initial
  // value may differ (e.g. ascertainment correction edge-cases returning -Inf).
  // A self-consistent logLik is required for the MH acceptance ratio.
  state->logLik = totalLogLik;

  // MARGINAL-K-INIT-001: under marginal_k the partition sum above is the
  // FIXED-kPrime likelihood (kPrime pinned to kObs), which sits ~+9.68 nats
  // ABOVE the true marginal-over-k likelihood. Left as the initial MH baseline
  // it freezes the WHOLE chain at the init tree for datasets where no early move
  // overcomes the inflation — the residual cause of the post-MARGINAL-K-SLICE-001
  // whole-chain freeze (~57% of SBC sims pinned at init tl=0.1*nEdge, p=0.5,
  // sigma=0.5). Recompute the init logLik via the marginal-aware dispatcher so
  // the baseline is correct; charLLCacheReady=false so that rebuild leaves a
  // coherent per-char cache.
  //
  // partLogLik is emptied because nothing refreshes it under marginal_k: it
  // would hold fixed-kPrime partition sums, and any `hasPLC` fast path reading
  // one commits a fixed-kPrime total as the marginal logLik. An empty cache is
  // what makes that class of fast path unselectable.
  if (data->marginalK) {
    state->charLLCacheReady = false;
    state->logLik = compute_full_loglik(*data, *state);
    state->partLogLik.clear();
  }
}


// ---------------------------------------------------------------------------
// allocate_cl_workspace: size and allocate the CL workspace in McmcState.
//
// Called once from R after fill_partition_cache(). Sizes the workspace to
// accommodate the largest per-partition pruning call needed given the current
// tree topology and kPrime values, with headroom for kPrime growth.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
void allocate_cl_workspace(SEXP dataPtr, SEXP statePtr) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();

  // nNode = max 1-indexed node in tree.
  // Flat pruning functions use maxNode = 2*nTip-1 (covers both rooted and
  // unrooted topologies).  Ensure the workspace is at least that large so
  // the fits() check succeeds and the workspace is actually used.
  int maxNode = 2 * data->nTip - 1;
  for (int i = 0; i < (int)state->parent.size(); ++i) {
    if (state->parent[i] > maxNode) maxNode = state->parent[i];
    if (state->child[i]  > maxNode) maxNode = state->child[i];
  }

  // kPrimeMax: current maximum kPrime across transformational characters,
  // with +4 headroom so reallocations are infrequent during MCMC.
  int kPrimeMax = 2;
  for (int gi = 0; gi < (int)data->transIdxGlobal.size(); ++gi) {
    int kp = state->kPrime[data->transIdxGlobal[gi]];
    if (kp > kPrimeMax) kPrimeMax = kp;
  }
  int kPrimeWithHeadroom = kPrimeMax + 4;

  // maxStride = max over partitions of (nCharPart * kMax_part).
  int maxStride = 0;
  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& pinfo = data->parts[pi];
    int nCharPart = pinfo.tipStates.ncol();
    int kMax = (pinfo.type == 0) ? 2 :
               (pinfo.type == 2) ? pinfo.k : kPrimeWithHeadroom;
    int stride = nCharPart * kMax;
    if (stride > maxStride) maxStride = stride;
  }
  if (maxStride < 2) maxStride = 2;

  state->clWs.allocate(maxNode, maxStride);

  // M-162B: pre-allocate scratch buffers for fused ascertainment + site_lik_sum.
  // ascStride: max pseudo-character stride across all partitions.
  //   JC/known: kMax (1 pseudo-char × kMax states, JC symmetry)
  //   MkN:      4   (2 pseudo-chars × 2 states)
  //   F81 het:  kMax² (k pseudo-chars × k states)
  // siteLikMax: max nChar across all partitions.
  int ascStride = 0;
  int maxNChar  = 0;
  bool useHet   = data->qHeterogeneity;
  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& pinfo = data->parts[pi];
    int nCharPart = pinfo.tipStates.ncol();
    if (nCharPart > maxNChar) maxNChar = nCharPart;

    int kMax = (pinfo.type == 0) ? 2 :
               (pinfo.type == 2) ? pinfo.k : kPrimeWithHeadroom;
    int partAsc;
    if (useHet) {
      partAsc = kMax * kMax;  // F81: k pseudo-chars × k states
    } else if (pinfo.type == 0) {
      partAsc = 4;  // MkN: 2 pseudo-chars × 2 states
    } else {
      partAsc = kMax;  // JC: 1 pseudo-char × k states (symmetry)
    }
    if (partAsc > ascStride) ascStride = partAsc;
  }

  state->clWs.allocateScratch(maxNode, ascStride, maxNChar);
}


// ---------------------------------------------------------------------------
// get_mcmc_state: extract R-accessible values from XPtr
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List get_mcmc_state(SEXP statePtr) {
  McmcState* s = Rcpp::XPtr<McmcState>(statePtr).get();
  IntegerMatrix edge(s->parent.size(), 2);
  for (int i = 0; i < s->parent.size(); ++i) {
    edge(i, 0) = s->parent[i];
    edge(i, 1) = s->child[i];
  }
  return List::create(
    _["edge"]          = edge,
    _["relBrLengths"]  = s->relBrLengths,
    _["treeLength"]    = s->treeLength,
    _["rateLoss"]      = s->rateLoss,
    _["rateLogSd"]     = s->rateLogSd,
    _["rateNeo"]       = s->rateNeo,
    _["p"]             = s->p,
    _["kprimeAlpha"]   = s->kprimeAlpha,
    _["kprimeBeta"]    = s->kprimeBeta,
    _["kPrime"]        = s->kPrime,
    _["logLik"]        = s->logLik,
    _["logPrior"]      = s->logPrior,
    _["logPost"]       = s->logLik + s->logPrior,
    _["betaScale"]     = s->betaScale,
    _["diagDirPartial"] = s->diagDirPartialCount,
    _["diagDirFullback"] = s->diagDirFullbackCount,
    _["diagNniPartial"] = s->diagNniPartialCount,
    _["diagBsPartial"]  = s->diagBsPartialCount,
    _["diagSelectivePop"] = s->nodeCL.diagSelectivePopCount,
    // Partition-API extras (zero-length / scalar defaults on legacy state).
    _["usePartitioned"]     = s->usePartitioned,
    _["classRateLogSd"]     = s->classRateLogSd,
    _["classW"]             = s->classW,
    _["classRate"]          = s->classRate,
    _["nCharPerClass"]      = s->nCharPerClass,
    _["etaNeo"]             = s->etaNeo,
    _["useHyperpriorOnSigma"] = s->useHyperpriorOnSigma,
    _["hyperTau"]           = s->hyperTau,
    _["classZ"]             = s->classZ
  );
}


// Fingerprint of the unrooted topology, cut to 53 bits so a double holds it
// exactly.  TreeTools needs parents ahead of children, which every move
// preserves even where it leaves the edges out of canonical preorder.
static double fnv_topo_hash(const IntegerVector& parent,
                            const IntegerVector& child, int nTip) {
  return static_cast<double>(TreeTools::topology_hash(
      parent.begin(), child.begin(), parent.size(), nTip) >> 11);
}

// [[Rcpp::export]]
double compute_topo_hash(IntegerVector parent, IntegerVector child, int nTip) {
  return fnv_topo_hash(parent, child, nTip);
}

// Exported for unit testing.  Case 6 takes the TreeNav SPR branch only while
// this cache is live.
// [[Rcpp::export]]
bool node_cl_ready(SEXP statePtr) {
  return Rcpp::XPtr<McmcState>(statePtr).get()->nodeCL.ready();
}


// [[Rcpp::export]]
double get_state_log_lik(SEXP statePtr) {
  return Rcpp::XPtr<McmcState>(statePtr).get()->logLik;
}


// Inspection helper for the marginal-k cache tiers (FU-3 tests).
// Exposes per-tier validity flags so tests can assert cache invariance
// under p-moves and invalidation under tree/branch/rate moves.
// [[Rcpp::export]]
List get_marginal_cache_state(SEXP statePtr) {
  McmcState* s = Rcpp::XPtr<McmcState>(statePtr).get();
  int nReady = 0;
  for (auto v : s->perKpClReady) if (v) ++nReady;
  return List::create(
    // Tier 1 (per-(char, k') cache shipped in PR-B):
    _["charLLReady"]   = s->charLLCacheReady,
    _["charLLSize"]    = static_cast<int>(s->charLLCache.size()),
    _["charLLNCand"]   = s->charLLNCand,
    // Tier 2 (per-(partition, node, k_offset) — FU-3 structural landing):
    _["perKpReadyAny"] = s->perKpClReadyAny,
    _["perKpSlots"]    = static_cast<int>(s->perKpClSlots.size()),
    _["perKpReady"]    = static_cast<int>(nReady)
  );
}


// Force-invalidate both marginal-k cache tiers (FU-3 tests).
// Useful for "baseline" cache-cold comparisons in the test harness.
// [[Rcpp::export]]
void invalidate_marginal_cache(SEXP statePtr) {
  McmcState* s = Rcpp::XPtr<McmcState>(statePtr).get();
  s->charLLCacheReady = false;
  s->invalidate_per_kp_cl_all();
}


// Recompute log prior for the current state (test helper for R↔C++ cross-check).
// [[Rcpp::export]]
double eval_log_prior_cpp(SEXP dataPtr, SEXP statePtr) {
  McmcData* d = Rcpp::XPtr<McmcData>(dataPtr);
  McmcState* s = Rcpp::XPtr<McmcState>(statePtr);
  return cpp_log_prior(
    *d, s->treeLength, s->relBrLengths,
    s->rateLoss, s->rateLogSd, s->rateNeo,
    s->p, s->kPrime, s->betaScale,
    s->kprimeAlpha, s->kprimeBeta);
}


// Partition-API: R-callable wrapper around cpp_log_prior_partitioned.
// Mirrors eval_log_prior_cpp (legacy scalar surface) for the partitioned path.
// classRateLogSd: length 1 (linked) or nClasses (unlinked).
// classW: length 1 (trivial) or nClasses (simplex).
// When classRateLogSd and classW each have length 1, the result must equal
// eval_log_prior_cpp to ~1e-10 (§7b analogue for the prior).
// [[Rcpp::export]]
double eval_log_prior_partitioned_cpp(
    SEXP dataPtr, SEXP statePtr,
    Rcpp::NumericVector classRateLogSd,
    Rcpp::NumericVector classW,
    double etaNeo,
    bool useHyperpriorOnSigma = false,
    double hyperTau = 1.0,
    Rcpp::NumericVector classZ = Rcpp::NumericVector::create()) {
  McmcData*  d = Rcpp::XPtr<McmcData>(dataPtr);
  McmcState* s = Rcpp::XPtr<McmcState>(statePtr);
  return cpp_log_prior_partitioned(
    *d, s->treeLength, s->relBrLengths,
    s->rateLoss, s->rateLogSd, s->rateNeo,
    s->p, s->kPrime,
    classRateLogSd, classW, etaNeo,
    useHyperpriorOnSigma, hyperTau, classZ,
    s->betaScale, s->kprimeAlpha, s->kprimeBeta);
}


// Setter for classRateConcentration (plan v4 §5.1). Avoids changing the
// prepare_mcmc_data signature; pattern mirrors set_branch_bins (M-090).
// [[Rcpp::export]]
void set_class_rate_concentration(SEXP dataPtr, double concentration) {
  McmcData* d = Rcpp::XPtr<McmcData>(dataPtr);
  d->classRateConcentration = concentration;
}


// ---------------------------------------------------------------------------
// compute_log_prior: dispatch to legacy or partitioned prior function.
//
// Mirrors the compute_full_loglik_at dispatch pattern from commit 2ba52bd.
// Two overloads:
//
//   compute_log_prior(data, state)
//     — reads all prior-relevant fields from state; use when MCMC is
//       evaluating the current state (no proposal).
//
//   compute_log_prior_at(data, state, relBrLengths)
//     — overrides the branch-length simplex (for Dirichlet move proposals
//       where only the topology / relBrLengths change but other fields stay
//       at state values).  All other fields read from state.
//
// ---------------------------------------------------------------------------

static double compute_log_prior_at(
    const McmcData& data, const McmcState& state,
    const NumericVector& relBrLengths) {
  if (state.usePartitioned) {
    return cpp_log_prior_partitioned(
      data, state.treeLength, relBrLengths,
      state.rateLoss, state.rateLogSd, state.rateNeo,
      state.p, state.kPrime,
      state.classRateLogSd, state.classW, state.etaNeo,
      state.useHyperpriorOnSigma, state.hyperTau, state.classZ,
      state.betaScale, state.kprimeAlpha, state.kprimeBeta);
  }
  return cpp_log_prior(
    data, state.treeLength, relBrLengths,
    state.rateLoss, state.rateLogSd, state.rateNeo,
    state.p, state.kPrime, state.betaScale,
    state.kprimeAlpha, state.kprimeBeta);
}

static double compute_log_prior(const McmcData& data, const McmcState& state) {
  return compute_log_prior_at(data, state, state.relBrLengths);
}


// ---------------------------------------------------------------------------
// preorder_into  (M-109)
//
// Lightweight preorder traversal into pre-allocated output buffers.
// Produces a valid (not necessarily canonical) preorder: parents before
// children.  Iterating the output backwards gives a valid postorder for
// Felsenstein pruning.  No Rcpp allocation — all output goes into caller-
// owned int*/double* arrays.
//
// For committing a chosen topology to state (where canonical ordering
// matters for subsequent in-place NNI), use TreeTools::preorder_weighted_impl
// instead.
// ---------------------------------------------------------------------------
static int preorder_into(const IntegerVector& parent,
                         const IntegerVector& child,
                         const NumericVector& edgeLen,
                         int nTip,
                         int* outParent, int* outChild, double* outLen) {
  const int nEdge = parent.size();
  const int root = nTip + 1;
  const int maxNode = 2 * nTip;

  // Build child-edge linked list: head[node] → first edge, nxt[e] → next
  std::vector<int> head(maxNode + 2, -1);
  std::vector<int> nxt(nEdge, -1);
  for (int e = nEdge - 1; e >= 0; --e) {
    nxt[e] = head[parent[e]];
    head[parent[e]] = e;
  }

  // DFS preorder from root
  int pos = 0;
  std::vector<int> stk;
  stk.reserve(nEdge);

  // Seed with root's children (reversed so first child is popped first)
  int cnt = 0;
  int rootEdges[4]; // root has ≤3 children (unrooted trifurcating)
  for (int e = head[root]; e >= 0; e = nxt[e])
    if (cnt < 4) rootEdges[cnt++] = e;
  for (int i = cnt - 1; i >= 0; --i)
    stk.push_back(rootEdges[i]);

  while (!stk.empty()) {
    int e = stk.back(); stk.pop_back();
    outParent[pos] = parent[e];
    outChild[pos]  = child[e];
    outLen[pos]    = edgeLen[e];
    ++pos;

    int ch = child[e];
    if (ch > nTip) {
      // Internal node: push children in reverse linked-list order
      int nc = 0;
      int ce[3]; // binary tree: ≤2 children per internal node
      for (int x = head[ch]; x >= 0; x = nxt[x])
        if (nc < 3) ce[nc++] = x;
      for (int i = nc - 1; i >= 0; --i)
        stk.push_back(ce[i]);
    }
  }
  return pos;
}


// ---------------------------------------------------------------------------
// compute_full_loglik / compute_full_loglik_at  (M-083)
//
// Pure C++ helpers for evaluating the total log-likelihood from within
// proposal code (GibbsSPR, GibbsSubtreeSwap, Weighted moves).  No R
// boundary crossing.  Both variants use the state's CL workspace if it
// has been allocated (allocate_cl_workspace was called).
//
// compute_full_loglik_at — evaluate at an arbitrary (parent, child, edgeLen),
//   but using state's kPrime / rateLoss / rateLogSd / rateNeo unchanged.
//   Does NOT update state->logLik or state->partLogLik.
//
// compute_full_loglik — convenience wrapper: evaluate at state's current tree.
// ---------------------------------------------------------------------------

static double compute_full_loglik_at(
    const McmcData& data, McmcState& state,
    const IntegerVector& parent,
    const IntegerVector& child,
    const NumericVector& edgeLen,
    bool fillCharLLCache = true) {
  // Marginal-k dispatch (v1: geometric arm only; known-k partitions and the
  // partition API are deferred and rejected by .RequireMarginalKSupported()).
  // cast away const on data: the marginal evaluator takes a non-const reference
  // because the underlying PR-A helper may grow state->gibbsWs and the cache
  // lives on state too.
  //
  // fillCharLLCache=false => SCRATCH eval (no charLLCache read/write): required
  // by multi-config moves so per-config marginal LLs are recomputed coherently
  // and the cache is left cold (FREEZE-003 follow-up). Ignored by the
  // non-marginal branches (they do not use charLLCache).
  if (data.marginalK) {
    return cpp_log_likelihood_marginal(
      const_cast<McmcData&>(data), state,
      parent, child, edgeLen,
      state.rateLoss, state.rateLogSd, state.rateNeo,
      state.clWs.ready() ? &state.clWs : nullptr,
      fillCharLLCache);
  }
  if (state.usePartitioned) {
    return cpp_log_likelihood_partitioned(
      data, parent, child, edgeLen,
      state.kPrime, state.rateLoss,
      state.classRateLogSd, state.classRate,
      state.etaNeo, state.betaScale,
      state.clWs.ready() ? &state.clWs : nullptr);
  }
  return cpp_log_likelihood(
    data, parent, child, edgeLen,
    state.kPrime, state.rateLoss, state.rateLogSd, state.rateNeo,
    state.betaScale,
    state.clWs.ready() ? &state.clWs : nullptr);
}

static double compute_full_loglik(const McmcData& data, McmcState& state) {
  int nEdge = state.relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state.treeLength * state.relBrLengths[i];
  return compute_full_loglik_at(data, state, state.parent, state.child, edgeLen);
}


// Rcpp-exported wrappers: used by R-level tests (test-full-loglik.R).
// [[Rcpp::export]]
double eval_full_loglik_cpp(SEXP dataPtr, SEXP statePtr) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  return compute_full_loglik(*data, *state);
}

// [[Rcpp::export]]
double eval_full_loglik_at_cpp(SEXP dataPtr, SEXP statePtr,
                                IntegerVector parent, IntegerVector child,
                                NumericVector edgeLen) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  return compute_full_loglik_at(*data, *state, parent, child, edgeLen);
}

// FREEZE-003 Phase-2 candidate-weight check (advisor-requested deterministic
// settler). The weighted topology moves SCORE candidates via preorder_into ->
// compute_full_loglik_at (the selection weight), but COMMIT via
// preorder_weighted_impl. `committed==cold` only exercises the commit path; a
// discrepancy between the two preorder canonicalisations would skew candidate
// SELECTION without showing up in any committed-LL test. This returns the
// marginal LL of the SAME (parent,child,edgeLen) computed through BOTH paths as
// SCRATCH evals; the test asserts bit-equality, so the LL used for selection is
// exactly the LL of the corresponding landed tree -> no proposal-selection skew.
// [[Rcpp::export]]
NumericVector eval_preorder_paths_cpp(SEXP dataPtr, SEXP statePtr,
                                      IntegerVector parent, IntegerVector child,
                                      NumericVector edgeLen) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  const int nTip  = data->nTip;
  const int nEdge = parent.size();
  // Candidate-evaluation path (M-109): in-place edits then preorder_into.
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);
  preorder_into(parent, child, edgeLen, nTip,
                INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
  double ll_into = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs,
                                          /*fillCharLLCache=*/false);
  // Commit path: TreeTools::preorder_weighted_impl.
  auto po = TreeTools::preorder_weighted_impl(parent, child, edgeLen);
  IntegerMatrix oe = po.first;
  NumericVector oaf = po.second;
  IntegerVector op = oe(_, 0);
  IntegerVector oc = oe(_, 1);
  double ll_weighted = compute_full_loglik_at(*data, *state, op, oc, oaf,
                                              /*fillCharLLCache=*/false);
  return NumericVector::create(_["preorder_into"]     = ll_into,
                               _["preorder_weighted"] = ll_weighted);
}


// ---------------------------------------------------------------------------
// gibbs_spr_impl  (M-085, M-105 partial CL; corrected kernel GSPR-001/004)
//
// Likelihood-weighted SPR with a full MH acceptance step.  The original
// "Gibbs" version committed the chosen regraft with a deterministic
// tau = 1/2 split and no accept step (GSPR-001: maps positive-measure sets
// into a pi-null set), and its candidate filter excluded the subtree's own
// position so the selection normaliser did not cancel (GSPR-004).  The
// corrected kernel (design: dev/plans/2026-08-12-gibbs-spr-fix-design.md
// §3a):
//
//   1. Prune a uniformly chosen eligible edge (u -> v).  Both endpoints of
//      the move prune to the SAME residual tree R, so a selection
//      distribution that is a function of R alone has the same normaliser
//      in both directions and it cancels exactly.
//   2. Enumerate ALL edges of R as regraft candidates: the surviving
//      original edges PLUS the merged (parentRow, sibRow) pair -- the
//      subtree's own position.  Weight each by exp(beta * logLik) with the
//      subtree reattached at a FIXED reference fraction tau = 1/2.  The
//      fixed reference is load-bearing: weights at a drawn tau would leave
//      direction-dependent normalisers that do not cancel.
//   3. Select an edge proportionally to weight, then draw the committed
//      split fraction tau ~ U(0,1) and evaluate the proposal there.
//   4. Accept by full MH:
//        beta * (logLik_y - logLik_x) + (logPrior_y - logPrior_x)
//          + log w_{e_x} - log w_{e_y} + log(lReg) - log(lMerge)
//      where the w-ratio is the selection-probability ratio (the shared
//      normaliser cancelled) and log(lReg) - log(lMerge) is the SPR
//      merge/split Jacobian, exactly as in spr_proposal_impl
//      (src/proposals.cpp:161).
//
// Drawing the merged edge is a legitimate branch-fraction move at unchanged
// topology (Jacobian 0), not a no-op: the old `if (rnd < wOrig) return
// false` early-out is gone.
//
// Three evaluation paths share gibbs_spr_plan / gibbs_spr_finish below and
// differ only in how the tau = 1/2 reference weights are computed: partial
// CL reuse (M-105), streaming Q-heterogeneity (M-114), or the
// full-evaluation fallback (M-083/M-109, also used under
// coding = "informative" per LIKE-001).
// ---------------------------------------------------------------------------

// Old full-evaluation path (used as fallback and for validation)
static bool gibbs_spr_impl_full(McmcData* data, McmcState* state, double beta);
// M-114: partial CL path for Q-heterogeneity
static bool gibbs_spr_impl_het(McmcData* data, McmcState* state, double beta);

// Shared plan: prune-edge choice and candidate enumeration.
struct GibbsSprPlan {
  int pruneRow = -1;             // edge u -> v
  int parentRow = -1;            // edge g -> u
  int sibRow = -1;               // edge u -> sibNode
  int u = -1, v = -1, sibNode = -1;
  std::vector<int> cands;        // regraft edge rows: E(R) minus the merged edge
  NumericVector absLen;          // absolute edge lengths of the current tree
  double lMerge = 0.0;           // absLen[parentRow] + absLen[sibRow]
  double lPrune = 0.0;           // absLen[pruneRow]
};

// Deterministic part: fill the plan for a GIVEN prune edge.  Returns false
// when the prune edge is degenerate (malformed rows or no regular
// candidates).
static bool gibbs_spr_plan_at(const McmcState* state, int nTip, int pruneRow,
                              GibbsSprPlan& plan) {
  const int nEdge = state->parent.size();
  plan.pruneRow = pruneRow;
  plan.u = state->parent[pruneRow];
  plan.v = state->child[pruneRow];

  plan.parentRow = plan.sibRow = plan.sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == plan.u) plan.parentRow = i;
    if (state->parent[i] == plan.u && state->child[i] != plan.v) {
      plan.sibRow = i; plan.sibNode = state->child[i];
    }
  }
  if (plan.parentRow < 0 || plan.sibRow < 0) return false;

  // BFS: mark descendants of v
  std::vector<bool> isDesc(2 * nTip + 2, false);
  isDesc[plan.v] = true;
  if (plan.v > nTip) {
    std::vector<int> queue = {plan.v};
    while (!queue.empty()) {
      int cur = queue.back(); queue.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (state->parent[i] == cur) {
          int c = state->child[i];
          isDesc[c] = true;
          if (c > nTip) queue.push_back(c);
        }
      }
    }
  }

  // Regraft candidates: every edge of the residual tree R.  The three edges
  // incident to u collapse in R to the single merged (parentRow, sibRow)
  // pair, appended as the extra candidate index nCand by the evaluation
  // paths; edges within the pruned subtree are not in R.
  plan.cands.clear();
  plan.cands.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[state->child[i]]) continue;
    if (state->parent[i] == plan.u || state->child[i] == plan.u) continue;
    plan.cands.push_back(i);
  }
  if (plan.cands.empty()) return false;

  plan.absLen = NumericVector(nEdge);
  for (int i = 0; i < nEdge; ++i)
    plan.absLen[i] = state->treeLength * state->relBrLengths[i];
  plan.lMerge = plan.absLen[plan.parentRow] + plan.absLen[plan.sibRow];
  plan.lPrune = plan.absLen[plan.pruneRow];
  return true;
}

// Random part: choose the prune edge uniformly among eligible edges.  The
// eligible count is nEdge - 3 independent of topology (exactly the root's
// three edges are excluded), so the choice probability cancels between the
// two directions of a move and needs no Hastings term.
static bool gibbs_spr_plan(const McmcData* data, const McmcState* state,
                           GibbsSprPlan& plan) {
  const int nEdge = state->parent.size();
  const int root  = data->nTip + 1;

  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  int pickIdx = (int)(R::unif_rand() * (double)eligible.size());
  if (pickIdx >= (int)eligible.size()) pickIdx = (int)eligible.size() - 1;
  return gibbs_spr_plan_at(state, data->nTip, eligible[pickIdx], plan);
}

// Shared selection + MH + commit for the three gibbs_spr evaluation paths.
//
// candLL[0 .. nCand-1] hold the log-likelihood of regrafting the pruned
// subtree at the tau = 1/2 midpoint of each plan.cands edge; candLL[nCand]
// is the same reference evaluation on the merged edge.  All entries must be
// produced by the same evaluation machinery (including any relabelling /
// ascertainment corrections) so the selection weights are a function of the
// residual tree alone.
//
// Deviation from §3a's cost table: the chosen-at-tau evaluation goes
// through compute_full_loglik_at (one full pruning) rather than the
// per-path candidate machinery.  That keeps the committed state->logLik
// canonical and lets all three paths share this finisher unchanged.
static bool gibbs_spr_finish(McmcData* data, McmcState* state, double beta,
                             const GibbsSprPlan& plan,
                             const std::vector<double>& candLL) {
  const int nCand = (int)plan.cands.size();
  const int nEdge = state->parent.size();

  double maxLL = R_NegInf;
  for (double ll : candLL) if (R_FINITE(ll) && ll > maxLL) maxLL = ll;
  if (!R_FINITE(maxLL)) return false;

  std::vector<double> ws(nCand + 1);
  double sumW = 0.0;
  for (int ci = 0; ci <= nCand; ++ci) {
    ws[ci] = R_FINITE(candLL[ci]) ? std::exp(beta * (candLL[ci] - maxLL)) : 0.0;
    sumW += ws[ci];
  }
  if (sumW <= 0.0) return false;

  // Select an edge of R: plan.cands first, the merged edge last.
  double rnd = R::unif_rand() * sumW;
  int chosen = nCand;
  for (int ci = 0; ci < nCand; ++ci) {
    if (rnd < ws[ci]) { chosen = ci; break; }
    rnd -= ws[ci];
  }
  const bool ontoMerged = (chosen == nCand);

  // Draw the committed split fraction (the tau = 1/2 reference above is only
  // for the selection weights).
  const double tau  = R::unif_rand();
  const double lReg = ontoMerged ? plan.lMerge
                                 : plan.absLen[plan.cands[chosen]];

  // Build the proposed tree.  Regrafting onto the merged edge reproduces the
  // current topology with the (parentRow, sibRow) pair re-split at tau.
  IntegerVector propPar = clone(state->parent);
  IntegerVector propCh  = clone(state->child);
  NumericVector propAbs = clone(plan.absLen);
  if (ontoMerged) {
    propAbs[plan.parentRow] = tau * plan.lMerge;
    propAbs[plan.sibRow]    = (1.0 - tau) * plan.lMerge;
  } else {
    const int rr = plan.cands[chosen];
    const int b  = state->child[rr];
    propCh[plan.parentRow] = plan.sibNode;
    propAbs[plan.parentRow] = plan.lMerge;
    propCh[rr]              = plan.u;
    propAbs[rr]             = tau * lReg;
    propPar[plan.sibRow]    = plan.u;
    propCh[plan.sibRow]     = b;
    propAbs[plan.sibRow]    = (1.0 - tau) * lReg;
  }

  auto po = TreeTools::preorder_weighted_impl(propPar, propCh, propAbs);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);

  // Evaluate the proposal at the drawn tau (NOT candLL[chosen], which is the
  // tau = 1/2 selection reference).
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbs,
                                            /*fillCharLLCache=*/false);
  if (!R_FINITE(newLogLik)) return false;

  NumericVector propRelBr(nEdge);
  for (int k = 0; k < nEdge; ++k)
    propRelBr[k] = ordAbs[k] / state->treeLength;
  double newLogPrior = compute_log_prior_at(*data, *state, propRelBr);
  if (!R_FINITE(newLogPrior)) return false;

  // MH: selection ratio (shared normaliser cancelled; reduces to the weight
  // log-ratio) + SPR merge/split Jacobian.  Both terms vanish when the
  // merged edge is drawn (lReg == lMerge), leaving a pure branch-fraction
  // MH step.
  double logAlpha = beta * (newLogLik - state->logLik)
                  + (newLogPrior - state->logPrior)
                  + beta * (candLL[nCand] - candLL[chosen])
                  + std::log(lReg) - std::log(plan.lMerge);
  if (!(R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha))
    return false;

  for (int k = 0; k < nEdge; ++k) {
    state->parent[k]       = op[k];
    state->child[k]        = oc[k];
    state->relBrLengths[k] = propRelBr[k];
  }
  state->logLik   = newLogLik;
  state->logPrior = newLogPrior;
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  return true;
}

// tau = 1/2 reference weights via partial CL reuse (M-105): one caching
// downpass per CLGroup, then an O(depth) update per candidate.  Fills
// candLL[0 .. nCand] (merged edge last), all through evaluate_candidate /
// evaluate_const_prob so every weight -- including the subtree's own
// position -- is computed on an identical footing (GSPR-004).
static void gibbs_spr_eval_partial(McmcData* data, McmcState* state,
                                   const GibbsSprPlan& plan,
                                   std::vector<double>& candLL) {
  const int nTip  = data->nTip;
  const int nCand = (int)plan.cands.size();

  // ===== M-105: Partial CL cache setup =====

  // Build tree navigation from current topology
  TreeNav topo;
  topo.build(state->parent, state->child, plan.absLen, nTip);

  // Compute ACRV rates
  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding = data->codingType;
  int maxNode = topo.maxNode;

  // Audit Issue 1: RB-style partition-rate normalisation (nChar-weighted
  // mean rate = 1). Both neo and trans/known unit rateScales are derived
  // from rate_neo + per-partition character counts; legacy "rateScale = 1
  // for trans" behaviour is recovered when nNeo == 0 or nTrans == 0.
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

  // Build CLGroups: one per (partition, kStates) evaluation unit
  std::vector<CLGroup> groups;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];

    if (part.type == 0) {
      // Neomorphic: k=2, MkN model
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));

    } else if (part.type == 2) {
      // Known state space: fixed k
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));

    } else {
      // Transformational: group by kPrime
      int nCharPart = part.tipStates.ncol();
      IntegerVector kPrimePart(nCharPart);
      for (int ci = 0; ci < nCharPart; ++ci)
        kPrimePart[ci] = state->kPrime[part.globalCharIdx[ci]];
      IntegerVector uniqKp = sort_unique(kPrimePart);

      for (int ui = 0; ui < uniqKp.size(); ++ui) {
        int kp = uniqKp[ui];
        std::vector<int> cols;
        for (int ci = 0; ci < nCharPart; ++ci)
          if (kPrimePart[ci] == kp) cols.push_back(ci);
        int nSub = (int)cols.size();

        IntegerMatrix sub(nTip, nSub);
        for (int c = 0; c < nSub; ++c)
          for (int t = 0; t < nTip; ++t)
            sub(t, c) = part.tipStates(t, cols[c]);

        CLGroup g;
        g.isMkN     = false;
        g.rateLoss  = 1.0;
        g.rateScale = pScales.trans;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
      }
    }
  }

  // Run caching downpass for each group
  for (auto& grp : groups)
    caching_downpass(grp, topo, state->parent, state->child, rates);

  // Compute residual CLs for each group (detach v from u)
  std::vector<ResidualCL> residuals(groups.size());
  for (size_t gi = 0; gi < groups.size(); ++gi)
    compute_residual_cl(residuals[gi], groups[gi], topo, rates,
                        plan.u, plan.sibNode, plan.lMerge);

  // Ascertainment correction: create pseudo-character groups (constant-site
  // patterns) and process through the same partial CL pipeline.
  std::vector<CLGroup> pseudoGroups;
  std::vector<ResidualCL> pseudoResiduals;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    pseudoResiduals.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, topo.maxNode, nCat);
      caching_downpass(pseudoGroups[gi], topo, state->parent, state->child, rates);
      compute_residual_cl(pseudoResiduals[gi], pseudoGroups[gi], topo, rates,
                          plan.u, plan.sibNode, plan.lMerge);
    }
  }

  // ===== Evaluate candidates (merged edge last) using partial CLs =====

  const int g = topo.parentNode[plan.u];

  candLL.assign(nCand + 1, 0.0);
  for (int ci = 0; ci <= nCand; ++ci) {
    const bool merged = (ci == nCand);
    const int rr = merged ? -1 : plan.cands[ci];
    const int a  = merged ? g : state->parent[rr];
    const int b  = merged ? plan.sibNode : state->child[rr];
    const double lHalf = 0.5 * (merged ? plan.lMerge : plan.absLen[rr]);

    double totalLL = 0.0;

    for (size_t gi = 0; gi < groups.size(); ++gi) {
      double grpLL = evaluate_candidate(
        groups[gi], topo, residuals[gi], rates,
        plan.v, plan.u, plan.sibNode, plan.lMerge, a, b, lHalf, plan.lPrune,
        nullptr, merged);

      // Ascertainment correction via pseudo-character partial CLs.
      // Note: only the constant-site term is computed here; coding == 2
      // (informative) additionally requires a singleton-site term that
      // this partial-CL path does not yet evaluate. The caller short-
      // circuits to gibbs_spr_impl_full when codingType == 2 (LIKE-001
      // interim), so we only reach this branch under coding == 1.
      if (coding != 0 && groups[gi].nChar > 0) {
        double constP = evaluate_const_prob(
          pseudoGroups[gi], topo, pseudoResiduals[gi], rates,
          plan.v, plan.u, plan.sibNode, plan.lMerge, a, b, lHalf, plan.lPrune,
          nullptr, merged);
        if (constP < 1.0)
          grpLL -= groups[gi].nChar * std::log(1.0 - constP);
      }

      totalLL += grpLL;
    }

    candLL[ci] = totalLL;
  }

  // Add relabelling correction to candLL — it's a topology-independent
  // constant that's included in state->logLik but not in evaluate_candidate.
  // Applied to every entry (including the merged edge) so all selection
  // weights share the same convention.
  if (data->relabel) {
    double relabelCorr = 0.0;
    for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
      const PartInfo& part = data->parts[pi];
      if (part.type == 1) {  // transformational only
        int nCharPart = part.tipStates.ncol();
        for (int ci = 0; ci < nCharPart; ++ci)
          relabelCorr += mk_prime_relabel_log(
            state->kPrime[part.globalCharIdx[ci]], part.kObsLocal[ci]);
      }
    }
    for (int ci = 0; ci <= nCand; ++ci)
      candLL[ci] += relabelCorr;
  }
}

static bool gibbs_spr_impl(McmcData* data, McmcState* state, double beta) {
  // M-114: Q-heterogeneity uses streaming partial CL
  if (data->qHeterogeneity)
    return gibbs_spr_impl_het(data, state, beta);

  // LIKE-001 interim (option 3): the pseudo-character partial-CL path
  // computes only the constant-site ascertainment term, not the singleton
  // term. Under coding="informative" (codingType == 2) that omission
  // biases the selection weights. Fall back to the full evaluator,
  // which routes through cpp_partition_log_likelihood with the singleton
  // correction applied. Restore partial-CL once evaluate_singleton_prob
  // (math-prover option 2) is implemented.
  if (data->codingType == 2)
    return gibbs_spr_impl_full(data, state, beta);

  GibbsSprPlan plan;
  if (!gibbs_spr_plan(data, state, plan)) return false;

  std::vector<double> candLL;
  gibbs_spr_eval_partial(data, state, plan, candLL);
  return gibbs_spr_finish(data, state, beta, plan, candLL);
}

// Deterministic enumeration hook for the GSPR-004 candidate-set symmetry
// test (tests/testthat/test-gibbs-spr-candidates.R).  Runs the same plan +
// tau = 1/2 reference evaluation as gibbs_spr_impl for a GIVEN prune edge
// (1-based row of the state's edge matrix) and returns the enumerated
// candidate edges — the merged (g, sibNode) pair last — with their
// reference log-likelihoods.  No RNG, no state mutation.
// [[Rcpp::export]]
List gibbs_spr_enumerate_cpp(SEXP dataPtr, SEXP statePtr, int pruneRow) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  if (data->qHeterogeneity || data->codingType == 2)
    Rcpp::stop("gibbs_spr_enumerate_cpp drives the partial-CL path only");
  const int nEdge = state->parent.size();
  if (pruneRow < 1 || pruneRow > nEdge)
    Rcpp::stop("pruneRow out of range");
  if (state->parent[pruneRow - 1] == data->nTip + 1)
    Rcpp::stop("prune edge must not hang off the root");

  GibbsSprPlan plan;
  if (!gibbs_spr_plan_at(state, data->nTip, pruneRow - 1, plan))
    Rcpp::stop("degenerate prune edge");  // # nocov: needs nTip < 4, but the root guard above fires first

  std::vector<double> candLL;
  gibbs_spr_eval_partial(data, state, plan, candLL);

  const int nCand = (int)plan.cands.size();
  IntegerMatrix edges(nCand + 1, 2);
  for (int ci = 0; ci < nCand; ++ci) {
    edges(ci, 0) = state->parent[plan.cands[ci]];
    edges(ci, 1) = state->child[plan.cands[ci]];
  }
  edges(nCand, 0) = state->parent[plan.parentRow];  // g
  edges(nCand, 1) = plan.sibNode;                   // merged edge (g, sibNode)

  return List::create(
    _["edges"]  = edges,
    _["logLik"] = NumericVector(candLL.begin(), candLL.end()),
    _["u"]      = plan.u,
    _["v"]      = plan.v,
    _["lMerge"] = plan.lMerge,
    _["lPrune"] = plan.lPrune);
}


// ---------------------------------------------------------------------------
// M-114: gibbs_spr with streaming partial CL for Q-heterogeneity.
//
// Same plan/selection/MH/commit logic as gibbs_spr_impl, but evaluates the
// tau = 1/2 reference weights under a mixture of F81 components by streaming
// over (betaBin, rotation) and accumulating per-site raw likelihoods.
// ---------------------------------------------------------------------------
static bool gibbs_spr_impl_het(McmcData* data, McmcState* state,
                                double beta) {
  const int nTip = data->nTip;

  GibbsSprPlan plan;
  if (!gibbs_spr_plan(data, state, plan)) return false;
  const int nCand = (int)plan.cands.size();
  // Evaluation slots: plan.cands first, the merged edge last.
  const int nEval = nCand + 1;

  TreeNav topo;
  topo.build(state->parent, state->child, plan.absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat,
                                          data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding  = data->codingType;
  int maxNode = topo.maxNode;
  int nBC     = data->nBetaCat;

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

  // Build CLGroups (same structure as non-het, but with useF81 flag)
  std::vector<CLGroup> groups;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
    } else {
      int nCharPart = part.tipStates.ncol();
      IntegerVector kPrimePart(nCharPart);
      for (int ci = 0; ci < nCharPart; ++ci)
        kPrimePart[ci] = state->kPrime[part.globalCharIdx[ci]];
      IntegerVector uniqKp = sort_unique(kPrimePart);
      for (int ui = 0; ui < uniqKp.size(); ++ui) {
        int kp = uniqKp[ui];
        std::vector<int> cols;
        for (int ci = 0; ci < nCharPart; ++ci)
          if (kPrimePart[ci] == kp) cols.push_back(ci);
        int nSub = (int)cols.size();
        IntegerMatrix sub(nTip, nSub);
        for (int c = 0; c < nSub; ++c)
          for (int t = 0; t < nTip; ++t)
            sub(t, c) = part.tipStates(t, cols[c]);
        CLGroup g;
        g.isMkN     = false;
        g.rateLoss  = 1.0;
        g.rateScale = pScales.trans;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
      }
    }
  }

  // Pseudo-character groups for ascertainment correction
  std::vector<CLGroup> pseudoGroups;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi)
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, maxNode, nCat);
  }

  const int g = topo.parentNode[plan.u];

  // ===== M-114: Streaming evaluation over (betaBin, rotation) =====
  //
  // Per-group, per-candidate accumulators for raw site likelihoods
  // siteLikAccums[gi][ci * nChar_gi .. (ci+1)*nChar_gi - 1]
  std::vector<std::vector<double>> siteLikAccums(groups.size());
  std::vector<std::vector<double>> constProbAccums(groups.size());
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    siteLikAccums[gi].assign((size_t)nEval * groups[gi].nChar, 0.0);
    if (coding != 0)
      constProbAccums[gi].assign(
        (size_t)nEval * pseudoGroups[gi].nChar, 0.0);
  }


  ResidualCL res;
  ResidualCL pseudoRes;

  // Process each group independently (different k → different bins/rotations)
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    CLGroup& grp = groups[gi];
    int k = grp.kStates;
    double baseRL = grp.isMkN ? state->rateLoss : 1.0;
    int nRot = (k == 2) ? 1 : k;

    double hetBins[16];
    gibbs_compute_het_bins(state->betaScale, k, nBC, hetBins);

    for (int bi = 0; bi < nBC; ++bi) {
      for (int rot = 0; rot < nRot; ++rot) {
        // Set F81 parameters for this component
        set_f81_component(grp, hetBins[bi], rot, baseRL);

        // Caching downpass + residual
        caching_downpass(grp, topo, state->parent, state->child, rates);
        compute_residual_cl(res, grp, topo, rates,
                            plan.u, plan.sibNode, plan.lMerge);

        // Same for pseudo-group (ascertainment)
        if (coding != 0) {
          set_f81_component(pseudoGroups[gi], hetBins[bi], rot, baseRL);
          caching_downpass(pseudoGroups[gi], topo, state->parent,
                           state->child, rates);
          compute_residual_cl(pseudoRes, pseudoGroups[gi], topo, rates,
                              plan.u, plan.sibNode, plan.lMerge);
        }

        // Evaluate all candidates (merged edge last) for this component
        for (int ci = 0; ci < nEval; ++ci) {
          const bool merged = (ci == nCand);
          const int rr = merged ? -1 : plan.cands[ci];
          const int a  = merged ? g : state->parent[rr];
          const int b  = merged ? plan.sibNode : state->child[rr];
          const double lHalf =
            0.5 * (merged ? plan.lMerge : plan.absLen[rr]);

          evaluate_candidate(
            grp, topo, res, rates,
            plan.v, plan.u, plan.sibNode, plan.lMerge, a, b, lHalf,
            plan.lPrune,
            siteLikAccums[gi].data() + (size_t)ci * grp.nChar, merged);

          if (coding != 0) {
            evaluate_const_prob(
              pseudoGroups[gi], topo, pseudoRes, rates,
              plan.v, plan.u, plan.sibNode, plan.lMerge, a, b, lHalf,
              plan.lPrune,
              constProbAccums[gi].data() +
                (size_t)ci * pseudoGroups[gi].nChar, merged);
          }
        }
      }
    }
  }

  // Convert accumulators to per-candidate log-likelihoods
  std::vector<double> candLL(nEval, 0.0);
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    const CLGroup& grp = groups[gi];
    int k     = grp.kStates;
    int nRot  = (k == 2) ? 1 : k;
    int totalComp = nCat * nBC * nRot;
    int nChar_gi  = grp.nChar;

    for (int ci = 0; ci < nEval; ++ci) {
      double grpLL = siteLikAccum_to_logLik(
        siteLikAccums[gi].data() + (size_t)ci * nChar_gi,
        nChar_gi, totalComp);

      if (coding != 0 && nChar_gi > 0) {
        int nPseudo = pseudoGroups[gi].nChar;  // = k
        double constP = 0.0;
        const double* cpa =
          constProbAccums[gi].data() + (size_t)ci * nPseudo;
        for (int c = 0; c < nPseudo; ++c) constP += cpa[c];
        constP /= totalComp;
        if (constP < 1.0)
          grpLL -= nChar_gi * std::log(1.0 - constP);
      }

      candLL[ci] += grpLL;
    }
  }

  // Relabelling correction (every entry, merged edge included)
  if (data->relabel) {
    double relabelCorr = 0.0;
    for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
      const PartInfo& part = data->parts[pi];
      if (part.type == 1) {
        int nCharPart = part.tipStates.ncol();
        for (int ci = 0; ci < nCharPart; ++ci)
          relabelCorr += mk_prime_relabel_log(
            state->kPrime[part.globalCharIdx[ci]], part.kObsLocal[ci]);
      }
    }
    for (int ci = 0; ci < nEval; ++ci)
      candLL[ci] += relabelCorr;
  }

  return gibbs_spr_finish(data, state, beta, plan, candLL);
}


// Old full-evaluation fallback (coding="informative" or validation), M-109
// in-place: tau = 1/2 reference weights via one full pruning per candidate.
static bool gibbs_spr_impl_full(McmcData* data, McmcState* state, double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  GibbsSprPlan plan;
  if (!gibbs_spr_plan(data, state, plan)) return false;
  const int nCand = (int)plan.cands.size();

  // Working copies: cloned ONCE, reused for all candidates (plan.absLen is
  // left pristine for gibbs_spr_finish)
  IntegerVector workPar = clone(state->parent);
  IntegerVector workCh  = clone(state->child);
  NumericVector workAbs = clone(plan.absLen);
  // Save original values for the 3 modified rows
  const int origPar_sibRow = workPar[plan.sibRow];
  const int origCh_parentRow = workCh[plan.parentRow];
  const int origCh_sibRow = workCh[plan.sibRow];
  const double origAbs_parentRow = workAbs[plan.parentRow];
  const double origAbs_sibRow = workAbs[plan.sibRow];

  // Pre-allocate output buffers for preorder_into
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  // Evaluate candidates (merged edge last).  Regrafting onto the merged
  // edge at tau = 1/2 is the current topology with the (parentRow, sibRow)
  // pair split evenly.
  std::vector<double> candLL(nCand + 1);
  for (int ci = 0; ci <= nCand; ++ci) {
    const bool merged = (ci == nCand);
    const int rr    = merged ? -1 : plan.cands[ci];
    const double lReg = merged ? plan.lMerge : workAbs[rr];
    const int bNode = merged ? -1 : workCh[rr];  // original child of regraft edge

    // Apply SPR in-place
    if (merged) {
      workAbs[plan.parentRow] = 0.5 * plan.lMerge;
      workAbs[plan.sibRow]    = 0.5 * plan.lMerge;
    } else {
      workCh[plan.parentRow] = plan.sibNode;
      workAbs[plan.parentRow] = plan.lMerge;
      workCh[rr]              = plan.u;
      workAbs[rr]             = 0.5 * lReg;
      workPar[plan.sibRow]    = plan.u;
      workCh[plan.sibRow]     = bNode;
      workAbs[plan.sibRow]    = 0.5 * lReg;
    }

    preorder_into(workPar, workCh, workAbs, nTip,
                  INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
    candLL[ci] = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs,
                                        /*fillCharLLCache=*/false);  // FREEZE-003

    // Restore
    if (merged) {
      workAbs[plan.parentRow] = origAbs_parentRow;
      workAbs[plan.sibRow]    = origAbs_sibRow;
    } else {
      workCh[plan.parentRow] = origCh_parentRow;
      workAbs[plan.parentRow] = origAbs_parentRow;
      workCh[rr]              = bNode;
      workAbs[rr]             = lReg;
      workPar[plan.sibRow]    = origPar_sibRow;
      workCh[plan.sibRow]     = origCh_sibRow;
      workAbs[plan.sibRow]    = origAbs_sibRow;
    }
  }

  return gibbs_spr_finish(data, state, beta, plan, candLL);
}


// ---------------------------------------------------------------------------
// gibbs_subtree_swap_impl  (M-086, M-109 in-place, M-111 partial CL;
// corrected kernel GSWAP-001)
//
// GibbsSubtreeSwap: anchor on a uniformly chosen non-root node, enumerate
// every valid swap partner, weight by exp(β × logLik), sample
// proportionally, then accept by MH at min(1, Z_x / Z_y).  Unlike gibbs_spr,
// whose two endpoints prune to one shared residual tree, the swap
// neighbourhood differs between the directions and its selection normaliser
// does not cancel: proof in
// dev/red-team/proofs/gibbs-subtree-swap-hastings.md.
//
// Branch lengths stay with their slots, not their subtrees: the swap
// transposes two entries of relBrLengths, so the Jacobian is 1 and the
// exchangeable Dirichlet branch prior is unchanged.
//
// Three evaluation paths fill the same neighbourhood and share one driver:
// partial CL reuse (M-111), streaming Q-heterogeneity (M-114), or full
// evaluation (also used under coding = "informative" per LIKE-001).
// ---------------------------------------------------------------------------

// Local helper: find edge row where child[i] == node
static int find_child_row_gibbs(const IntegerVector& child, int node) {
  for (int i = 0; i < child.size(); ++i)
    if (child[i] == node) return i;
  return -1;
}

static bool swap_neighbourhood_het(McmcData* data, McmcState* state,
                                   const IntegerVector& parent,
                                   const IntegerVector& child,
                                   const NumericVector& absLen,
                                   int nodeA,
                                   std::vector<int>& partners,
                                   std::vector<double>& candLL);

static bool swap_neighbourhood_full(McmcData* data, McmcState* state,
                                    const IntegerVector& parent,
                                    const IntegerVector& child,
                                    const NumericVector& absLenIn,
                                    int nodeA,
                                    std::vector<int>& partners,
                                    std::vector<double>& candLL);

// Selection normaliser of the anchored swap neighbourhood: log of the sum of
// exp(beta * logLik) over {stay} and every swap of nodeA with a valid
// partner.
static double swap_neighbourhood_logz(double selfLL, double beta,
                                      const std::vector<double>& candLL) {
  double maxLL = selfLL;
  for (double ll : candLL)
    if (R_FINITE(ll)) maxLL = std::max(maxLL, ll);
  if (!R_FINITE(maxLL)) return R_NegInf;
  double s = std::exp(beta * (selfLL - maxLL));
  for (double ll : candLL)
    if (R_FINITE(ll)) s += std::exp(beta * (ll - maxLL));
  return beta * maxLL + std::log(s);
}

// Log-likelihood of every swap of nodeA within the tree (parent, child,
// absLen), by partial CL reuse (M-111): one caching downpass per CLGroup,
// then an O(depth) update per candidate.  Both directions of the move go
// through this, so their normalisers rest on an identical footing.
//
// Edges need only be in a valid preorder, not the canonical one: the
// neighbourhood and its likelihoods are properties of the tree, not of the
// node labelling that canonicalisation permutes.
static bool swap_neighbourhood_partial(McmcData* data, McmcState* state,
                                       const IntegerVector& parent,
                                       const IntegerVector& child,
                                       const NumericVector& absLen,
                                       int nodeA,
                                       std::vector<int>& partners,
                                       std::vector<double>& candLL) {
  const int nTip = data->nTip;

  partners = get_valid_swap_partners_impl(parent, child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  TreeNav topo;
  topo.build(parent, child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding = data->codingType;
  int maxNode = topo.maxNode;

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

  // Build CLGroups: one per (partition, kStates) evaluation unit
  std::vector<CLGroup> groups;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];

    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));

    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));

    } else {
      int nCharPart = part.tipStates.ncol();
      IntegerVector kPrimePart(nCharPart);
      for (int ci = 0; ci < nCharPart; ++ci)
        kPrimePart[ci] = state->kPrime[part.globalCharIdx[ci]];
      IntegerVector uniqKp = sort_unique(kPrimePart);

      for (int ui = 0; ui < uniqKp.size(); ++ui) {
        int kp = uniqKp[ui];
        std::vector<int> cols;
        for (int ci = 0; ci < nCharPart; ++ci)
          if (kPrimePart[ci] == kp) cols.push_back(ci);
        int nSub = (int)cols.size();

        IntegerMatrix sub(nTip, nSub);
        for (int c = 0; c < nSub; ++c)
          for (int t = 0; t < nTip; ++t)
            sub(t, c) = part.tipStates(t, cols[c]);

        CLGroup g;
        g.isMkN     = false;
        g.rateLoss  = 1.0;
        g.rateScale = pScales.trans;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
      }
    }
  }

  // Run caching downpass for each group
  for (auto& grp : groups)
    caching_downpass(grp, topo, parent, child, rates);

  // Ascertainment correction: pseudo-character groups
  std::vector<CLGroup> pseudoGroups;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, maxNode, nCat);
      caching_downpass(pseudoGroups[gi], topo, parent, child, rates);
    }
  }

  // Precompute nodeA-fixed data for evaluate_swap_impl
  int pA    = topo.parentNode[nodeA];
  int slotA = topo.childSlot(pA, nodeA);
  double lenA = topo.edgeLen[topo.edgeToPar[nodeA]];

  std::vector<int> pathA;
  pathA.reserve(16);
  for (int n = pA; n >= 1; n = topo.parentNode[n])
    pathA.push_back(n);

  std::vector<int> pathAIdx(maxNode + 1, -1);
  for (int i = 0; i < (int)pathA.size(); ++i)
    pathAIdx[pathA[i]] = i;

  candLL.assign(nPart, 0.0);
  for (int pi = 0; pi < nPart; ++pi) {
    double totalLL = 0.0;

    for (size_t gi = 0; gi < groups.size(); ++gi) {
      double grpLL = evaluate_swap_candidate(
        groups[gi], topo, rates, nodeA, partners[pi],
        pA, slotA, lenA, pathA, pathAIdx);

      // Ascertainment correction
      if (coding != 0 && groups[gi].nChar > 0) {
        double constP = evaluate_swap_const_prob(
          pseudoGroups[gi], topo, rates, nodeA, partners[pi],
          pA, slotA, lenA, pathA, pathAIdx);
        if (constP < 1.0)
          grpLL -= groups[gi].nChar * std::log(1.0 - constP);
      }

      totalLL += grpLL;
    }

    candLL[pi] = totalLL;
  }

  // Relabelling correction (topology-independent constant)
  if (data->relabel) {
    double relabelCorr = 0.0;
    for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
      const PartInfo& part = data->parts[pi];
      if (part.type == 1) {
        int nCharPart = part.tipStates.ncol();
        for (int ci = 0; ci < nCharPart; ++ci)
          relabelCorr += mk_prime_relabel_log(
            state->kPrime[part.globalCharIdx[ci]], part.kObsLocal[ci]);
      }
    }
    for (int ci = 0; ci < nPart; ++ci)
      candLL[ci] += relabelCorr;
  }

  return true;
}

// Which evaluator fills the neighbourhood; the kernel is identical for all.
enum class SwapEval { Partial, Het, Full };

static bool swap_neighbourhood(McmcData* data, McmcState* state, SwapEval ev,
                               const IntegerVector& parent,
                               const IntegerVector& child,
                               const NumericVector& absLen, int nodeA,
                               std::vector<int>& partners,
                               std::vector<double>& candLL) {
  switch (ev) {
    case SwapEval::Het:
      return swap_neighbourhood_het(data, state, parent, child, absLen,
                                    nodeA, partners, candLL);
    case SwapEval::Full:
      return swap_neighbourhood_full(data, state, parent, child, absLen,
                                     nodeA, partners, candLL);
    default:
      return swap_neighbourhood_partial(data, state, parent, child, absLen,
                                        nodeA, partners, candLL);
  }
}

static bool gibbs_subtree_swap_impl(McmcData* data, McmcState* state,
                                    double beta) {
  // The pseudo-character partial-CL path omits the singleton-site
  // ascertainment term required under coding = "informative" (LIKE-001).
  const SwapEval ev = data->qHeterogeneity ? SwapEval::Het
                    : (data->codingType == 2 ? SwapEval::Full
                                             : SwapEval::Partial);

  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  // 1. Anchor on a uniformly chosen node (every edge child is non-root)
  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // 2. Evaluate the anchored neighbourhood of the current state
  std::vector<int> partners;
  std::vector<double> candLL;
  if (!swap_neighbourhood(data, state, ev, state->parent, state->child,
                          absLen, nodeA, partners, candLL))
    return false;
  const int nPart = (int)partners.size();

  // 3. Sampling weights: exp(beta * logLik), current state included
  const double llOrig = state->logLik;
  double maxLL = llOrig;
  for (int pi = 0; pi < nPart; ++pi)
    if (R_FINITE(candLL[pi])) maxLL = std::max(maxLL, candLL[pi]);

  double wOrig = std::exp(beta * (llOrig - maxLL));
  std::vector<double> ws(nPart);
  double sumW = wOrig;
  for (int pi = 0; pi < nPart; ++pi) {
    ws[pi] = R_FINITE(candLL[pi]) ? std::exp(beta * (candLL[pi] - maxLL)) : 0.0;
    sumW  += ws[pi];
  }
  if (!(sumW > 0.0)) return false;

  // 4. Sample; a self-draw is a no-op
  double rnd = R::unif_rand() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nPart - 1;
  for (int pi = 0; pi < nPart - 1; ++pi) {
    if (rnd < ws[pi]) { chosen = pi; break; }
    rnd -= ws[pi];
  }
  if (!R_FINITE(candLL[chosen])) return false;

  // 5. Build the proposal.  Node labels are preserved so nodeA still anchors
  //    the reverse neighbourhood; canonicalisation permutes internal labels,
  //    so it waits until the move is committed.  The evaluator attaches each
  //    subtree with the stem length of the slot it moves into, so the lengths
  //    are transposed with the parents.
  const int rowA = find_child_row_gibbs(state->child, nodeA);
  const int rowB = find_child_row_gibbs(state->child, partners[chosen]);
  if (rowA < 0 || rowB < 0) return false;

  IntegerVector propPar = clone(state->parent);
  NumericVector propAbs = clone(absLen);
  propPar[rowA] = state->parent[rowB];
  propPar[rowB] = state->parent[rowA];
  propAbs[rowA] = absLen[rowB];
  propAbs[rowB] = absLen[rowA];

  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);
  preorder_into(propPar, state->child, propAbs, nTip,
                INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));

  // 6. Accept with min(1, Z_x / Z_y)
  std::vector<int> revPartners;
  std::vector<double> revLL;
  if (!swap_neighbourhood(data, state, ev, ordPar, ordCh, ordAbs,
                          nodeA, revPartners, revLL))
    return false;

  const double logAlpha = swap_neighbourhood_logz(llOrig, beta, candLL)
                        - swap_neighbourhood_logz(candLL[chosen], beta, revLL);
  if (!(R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha))
    return false;

  // 7. Commit, in canonical preorder (the in-place NNI invariant needs it)
  state->parent[rowA] = propPar[rowA];
  state->parent[rowB] = propPar[rowB];
  const double tmpRel = state->relBrLengths[rowA];
  state->relBrLengths[rowA] = state->relBrLengths[rowB];
  state->relBrLengths[rowB] = tmpRel;

  auto po = TreeTools::preorder_weighted_impl(
      state->parent, state->child, propAbs);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbsFinal = po.second;
  for (int k = 0; k < nEdge; ++k) {
    state->parent[k]       = ordEdge(k, 0);
    state->child[k]        = ordEdge(k, 1);
    state->relBrLengths[k] = ordAbsFinal[k] / state->treeLength;
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  return true;
}


// M-114 streaming partial CL: same neighbourhood under Q-heterogeneity.
static bool swap_neighbourhood_het(McmcData* data, McmcState* state,
                                   const IntegerVector& parent,
                                   const IntegerVector& child,
                                   const NumericVector& absLen,
                                   int nodeA,
                                   std::vector<int>& partners,
                                   std::vector<double>& candLL) {
  const int nTip = data->nTip;

  partners = get_valid_swap_partners_impl(parent, child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  TreeNav topo;
  topo.build(parent, child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat,
                                          data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding  = data->codingType;
  int maxNode = topo.maxNode;
  int nBC     = data->nBetaCat;

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

  // Build CLGroups (same as gibbs_spr_impl_het)
  std::vector<CLGroup> groups;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
    } else {
      int nCharPart = part.tipStates.ncol();
      IntegerVector kPrimePart(nCharPart);
      for (int ci = 0; ci < nCharPart; ++ci)
        kPrimePart[ci] = state->kPrime[part.globalCharIdx[ci]];
      IntegerVector uniqKp = sort_unique(kPrimePart);
      for (int ui = 0; ui < uniqKp.size(); ++ui) {
        int kp = uniqKp[ui];
        std::vector<int> cols;
        for (int ci = 0; ci < nCharPart; ++ci)
          if (kPrimePart[ci] == kp) cols.push_back(ci);
        int nSub = (int)cols.size();
        IntegerMatrix sub(nTip, nSub);
        for (int c = 0; c < nSub; ++c)
          for (int t = 0; t < nTip; ++t)
            sub(t, c) = part.tipStates(t, cols[c]);
        CLGroup g;
        g.isMkN     = false;
        g.rateLoss  = 1.0;
        g.rateScale = pScales.trans;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
      }
    }
  }

  std::vector<CLGroup> pseudoGroups;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi)
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, maxNode, nCat);
  }

  // Precompute nodeA-fixed data
  int pA    = topo.parentNode[nodeA];
  int slotA = topo.childSlot(pA, nodeA);
  double lenA = topo.edgeLen[topo.edgeToPar[nodeA]];

  std::vector<int> pathA;
  pathA.reserve(16);
  for (int n = pA; n >= 1; n = topo.parentNode[n])
    pathA.push_back(n);

  std::vector<int> pathAIdx(maxNode + 1, -1);
  for (int i = 0; i < (int)pathA.size(); ++i)
    pathAIdx[pathA[i]] = i;

  // Per-group, per-candidate accumulators
  std::vector<std::vector<double>> siteLikAccums(groups.size());
  std::vector<std::vector<double>> constProbAccums(groups.size());
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    siteLikAccums[gi].assign((size_t)nPart * groups[gi].nChar, 0.0);
    if (coding != 0)
      constProbAccums[gi].assign(
        (size_t)nPart * pseudoGroups[gi].nChar, 0.0);
  }

  // Stream over (betaBin, rotation) per group
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    CLGroup& grp = groups[gi];
    int k = grp.kStates;
    double baseRL = grp.isMkN ? state->rateLoss : 1.0;
    int nRot = (k == 2) ? 1 : k;

    double hetBins[16];
    gibbs_compute_het_bins(state->betaScale, k, nBC, hetBins);

    for (int bi = 0; bi < nBC; ++bi) {
      for (int rot = 0; rot < nRot; ++rot) {
        set_f81_component(grp, hetBins[bi], rot, baseRL);
        caching_downpass(grp, topo, parent, child, rates);

        if (coding != 0) {
          set_f81_component(pseudoGroups[gi], hetBins[bi], rot, baseRL);
          caching_downpass(pseudoGroups[gi], topo, parent, child, rates);
        }

        for (int pi2 = 0; pi2 < nPart; ++pi2) {
          evaluate_swap_impl(
            grp, topo, rates, nodeA, partners[pi2],
            pA, slotA, lenA, pathA, pathAIdx, false,
            siteLikAccums[gi].data() + (size_t)pi2 * grp.nChar);

          if (coding != 0) {
            evaluate_swap_impl(
              pseudoGroups[gi], topo, rates, nodeA, partners[pi2],
              pA, slotA, lenA, pathA, pathAIdx, true,
              nullptr,
              constProbAccums[gi].data() +
                (size_t)pi2 * pseudoGroups[gi].nChar);
          }
        }
      }
    }
  }

  // Convert accumulators to per-candidate log-likelihoods
  candLL.assign(nPart, 0.0);
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    const CLGroup& grp = groups[gi];
    int k     = grp.kStates;
    int nRot  = (k == 2) ? 1 : k;
    int totalComp = nCat * nBC * nRot;
    int nChar_gi  = grp.nChar;

    for (int pi2 = 0; pi2 < nPart; ++pi2) {
      double grpLL = siteLikAccum_to_logLik(
        siteLikAccums[gi].data() + (size_t)pi2 * nChar_gi,
        nChar_gi, totalComp);

      if (coding != 0 && nChar_gi > 0) {
        int nPseudo = pseudoGroups[gi].nChar;
        double constP = 0.0;
        const double* cpa =
          constProbAccums[gi].data() + (size_t)pi2 * nPseudo;
        for (int c = 0; c < nPseudo; ++c) constP += cpa[c];
        constP /= totalComp;
        if (constP < 1.0)
          grpLL -= nChar_gi * std::log(1.0 - constP);
      }

      candLL[pi2] += grpLL;
    }
  }

  // Relabelling correction
  if (data->relabel) {
    double relabelCorr = 0.0;
    for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
      const PartInfo& part = data->parts[pi];
      if (part.type == 1) {
        int nCharPart = part.tipStates.ncol();
        for (int ci = 0; ci < nCharPart; ++ci)
          relabelCorr += mk_prime_relabel_log(
            state->kPrime[part.globalCharIdx[ci]], part.kObsLocal[ci]);
      }
    }
    for (int pi2 = 0; pi2 < nPart; ++pi2)
      candLL[pi2] += relabelCorr;
  }

  return true;
}


// Full-evaluation fallback for Q-heterogeneity (M-109 in-place pattern)
// Reference evaluator: same neighbourhood without partial CL reuse, for
// coding = "informative", where the pseudo-character path omits the
// singleton-site ascertainment term (LIKE-001).
static bool swap_neighbourhood_full(McmcData* data, McmcState* state,
                                    const IntegerVector& parent,
                                    const IntegerVector& child,
                                    const NumericVector& absLenIn,
                                    int nodeA,
                                    std::vector<int>& partners,
                                    std::vector<double>& candLL) {
  const int nEdge = parent.size();
  const int nTip  = data->nTip;

  partners = get_valid_swap_partners_impl(parent, child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  const int rowA = find_child_row_gibbs(child, nodeA);
  if (rowA < 0) return false;

  IntegerVector workPar = clone(parent);
  NumericVector absLen  = clone(absLenIn);

  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  const int origParA = workPar[rowA];
  const double origAbsA = absLen[rowA];

  candLL.assign(nPart, 0.0);
  for (int pi = 0; pi < nPart; ++pi) {
    int rowB = find_child_row_gibbs(child, partners[pi]);
    if (rowB < 0 || rowA == rowB) {
      candLL[pi] = R_NegInf;
      continue;
    }

    int origParB = workPar[rowB];
    double origAbsB = absLen[rowB];

    workPar[rowA] = origParB;
    workPar[rowB] = origParA;
    absLen[rowA]  = origAbsB;
    absLen[rowB]  = origAbsA;

    preorder_into(workPar, child, absLen, nTip,
                  INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
    candLL[pi] = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs);

    workPar[rowA] = origParA;
    workPar[rowB] = origParB;
    absLen[rowA]  = origAbsA;
    absLen[rowB]  = origAbsB;
  }

  return true;
}


// BranchBins struct now lives in mcmc_state.h and is precomputed by
// set_branch_bins() at MCMC init — no static globals or lazy init.

// Log proposal density of a branch fraction drawn by selecting a bin with
// probability proportional to `weights`, then drawing from the Beta centred
// on that bin's midpoint.  The bin is auxiliary and not carried in the state,
// so it is integrated out: the Beta components have full support on (0, 1),
// so any bin could have produced the fraction drawn.  `weights` need not be
// normalised; its sum is the same in both directions and cancels.
static double bin_mixture_logdensity(const std::vector<double>& weights,
                                     double f, const BranchBins& bins) {
  const double conc = bins.concentration;
  double logMax = R_NegInf;
  std::vector<double> terms(bins.nBins, R_NegInf);
  for (int b = 0; b < bins.nBins; ++b) {
    if (!(weights[b] > 0.0)) continue;
    const double mid = bins.mids[b];
    terms[b] = std::log(weights[b])
             + R::dbeta(f, mid * conc + 1.0, (1.0 - mid) * conc + 1.0, 1);
    if (terms[b] > logMax) logMax = terms[b];
  }
  if (!R_FINITE(logMax)) return R_NegInf;
  double sum = 0.0;
  for (int b = 0; b < bins.nBins; ++b)
    if (R_FINITE(terms[b])) sum += std::exp(terms[b] - logMax);
  return logMax + std::log(sum);
}

// Validation hook: the same mixture the weighted moves' Hastings ratio uses.
// [[Rcpp::export]]
double bin_mixture_log_density(NumericVector weights, double f, int nBins) {
  BranchBins bins;
  bins.init(nBins);
  return bin_mixture_logdensity(
    std::vector<double>(weights.begin(), weights.end()), f, bins);
}


// ---------------------------------------------------------------------------
// weighted_branch_scale_impl  (M-087)
//
// WeightedBranchLengthScale: pick two branches, discretise the branch-
// fraction space into B bins (Beta(0.25,0.25) quantile breakpoints),
// evaluate log-likelihood at each bin midpoint, weight by exp(beta * LL),
// Gibbs-sample a bin, then draw a new fraction from a Beta centred on the
// selected midpoint.  Returns logHastings for standard MH acceptance.
//
// The bin weights depend only on the edges the move leaves alone, so they
// and their normaliser are the same in both directions; the bin is integrated
// out.  Derivation: dev/red-team/proofs/weighted-branch-hastings.md.
// ---------------------------------------------------------------------------
static bool weighted_branch_scale_impl(
    McmcData* data, McmcState* state, double beta,
    double& logHastings,
    int& outIdx1, int& outIdx2, double& outOld1, double& outOld2) {

  const int nEdge = state->relBrLengths.size();
  if (nEdge < 2) return false;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Pick two branches (same scheme as beta_simplex_impl)
  int index = static_cast<int>(R::unif_rand() * nEdge);
  if (index >= nEdge) index = nEdge - 1;
  int other = static_cast<int>(R::unif_rand() * (nEdge - 1));
  if (other >= index) ++other;
  if (other >= nEdge) other = nEdge - 1;
  if (other == index) other = (index + 1) % nEdge;

  const double oldRelA = state->relBrLengths[index];
  const double oldRelB = state->relBrLengths[other];

  // Output for O(1) rollback
  outIdx1 = index; outIdx2 = other;
  outOld1 = oldRelA; outOld2 = oldRelB;
  const double relTotal = oldRelA + oldRelB;
  if (relTotal <= 0.0) return false;
  const double oldF = oldRelA / relTotal;

  // 2. Absolute edge lengths (modify only [index] and [other] per midpoint)
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];
  const double absTotal = absLen[index] + absLen[other];

  // 3. Evaluate log-likelihood at each bin midpoint
  std::vector<double> midLL(nBins);
  NumericVector trialAbs = clone(absLen);
  for (int b = 0; b < nBins; ++b) {
    const double mid = bins.mids[b];
    trialAbs[index] = mid * absTotal;
    trialAbs[other] = (1.0 - mid) * absTotal;
    midLL[b] = compute_full_loglik_at(*data, *state,
                                       state->parent, state->child, trialAbs,
                                       /*fillCharLLCache=*/false);  // FREEZE-003
  }

  // 4. Compute weights: exp(beta * LL), offset for numerical stability
  double maxLL = midLL[0];
  for (int b = 1; b < nBins; ++b)
    if (R_FINITE(midLL[b]) && midLL[b] > maxLL) maxLL = midLL[b];
  if (!R_FINITE(maxLL)) return false;

  std::vector<double> weights(nBins);
  double sumW = 0.0;
  for (int b = 0; b < nBins; ++b) {
    weights[b] = R_FINITE(midLL[b]) ?
                   std::exp(beta * (midLL[b] - maxLL)) : 0.0;
    sumW += weights[b];
  }
  if (sumW <= 0.0) return false;

  // 5. Sample a bin
  double rnd = R::unif_rand() * sumW;
  int chosenBin = nBins - 1;
  {
    double cum = 0.0;
    for (int b = 0; b < nBins; ++b) {
      cum += weights[b];
      if (rnd < cum) { chosenBin = b; break; }
    }
  }

  // 6. Draw fraction from Beta centred on chosen bin's midpoint.
  const double conc = bins.concentration;
  const double chosenMid = bins.mids[chosenBin];
  const double alphaNew = chosenMid * conc + 1.0;
  const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
  double newF = R::rbeta(alphaNew, betaNew);
  if (newF < 1e-8) newF = 1e-8;
  if (newF > 1.0 - 1e-8) newF = 1.0 - 1e-8;

  // 7. Hastings ratio.  The weights depend only on edges the move leaves
  //    alone, so they are identical in both directions.
  logHastings = bin_mixture_logdensity(weights, oldF, bins)
              - bin_mixture_logdensity(weights, newF, bins);

  // 8. Apply proposed fraction
  state->relBrLengths[index] = newF * relTotal;
  state->relBrLengths[other] = (1.0 - newF) * relTotal;
  return true;
}


// ---------------------------------------------------------------------------
// block_gibbs_branch_sweep_impl  (M-054 reframed)
//
// Random-permutation-scan MH-within-Gibbs sweep over ALL edge pairs.
// For each pair, uses the same bin-based approximate-conditional sampling
// as weighted_branch_scale_impl (M-087): evaluate LL at B bin midpoints,
// weight by exp(beta * LL), sample a bin, draw a fraction from a Beta
// centred on the bin midpoint, and accept/reject via MH.
//
// Each pair is accepted/rejected independently (composition of valid MH
// kernels).  The sweep returns true if at least one pair was accepted.
//
// The Dirichlet(1,...,1) prior on relBrLengths is constant w.r.t. branch
// values (only requires positivity), so prior recomputation within the
// sweep is unnecessary — only the likelihood changes.
//
// Cost: nEdge * (nBins + 1) full likelihood evaluations per sweep.
// ---------------------------------------------------------------------------
static bool block_gibbs_branch_sweep_impl(
    McmcData* data, McmcState* state, double beta) {

  const int nEdge = state->relBrLengths.size();
  if (nEdge < 2) return false;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;
  const double conc = bins.concentration;

  // Fisher-Yates shuffle for random permutation scan
  std::vector<int> perm(nEdge);
  for (int i = 0; i < nEdge; ++i) perm[i] = i;
  for (int i = nEdge - 1; i > 0; --i) {
    int j = static_cast<int>(R::unif_rand() * (i + 1));
    if (j > i) j = i;
    std::swap(perm[i], perm[j]);
  }

  // Working copy of absolute edge lengths (updated in-place across sweep)
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // Trial vector — shares memory with absLen except for two modified entries
  NumericVector trialAbs = clone(absLen);

  int nAccepted = 0;
  double currentLL = state->logLik;

  std::vector<double> midLL(nBins);
  std::vector<double> weights(nBins);

  for (int pi = 0; pi < nEdge; ++pi) {
    int index = perm[pi];

    // Pick a random partner edge
    int other = static_cast<int>(R::unif_rand() * (nEdge - 1));
    if (other >= index) ++other;
    if (other >= nEdge) other = nEdge - 1;
    if (other == index) other = (index + 1) % nEdge;

    const double oldRelA = state->relBrLengths[index];
    const double oldRelB = state->relBrLengths[other];
    const double relTotal = oldRelA + oldRelB;
    if (relTotal <= 0.0) continue;
    const double oldF = oldRelA / relTotal;
    const double absTotal = absLen[index] + absLen[other];
    if (absTotal <= 0.0) continue;

    // Evaluate LL at each bin midpoint
    // Reset trial vector to current state for these two edges
    for (int i = 0; i < nEdge; ++i) trialAbs[i] = absLen[i];

    for (int b = 0; b < nBins; ++b) {
      const double mid = bins.mids[b];
      trialAbs[index] = mid * absTotal;
      trialAbs[other] = (1.0 - mid) * absTotal;
      midLL[b] = compute_full_loglik_at(*data, *state,
                                         state->parent, state->child, trialAbs,
                                         /*fillCharLLCache=*/false);  // FREEZE-003
    }

    // Weight bins: exp(beta * (LL - maxLL))
    double maxLL = midLL[0];
    for (int b = 1; b < nBins; ++b)
      if (R_FINITE(midLL[b]) && midLL[b] > maxLL) maxLL = midLL[b];
    if (!R_FINITE(maxLL)) continue;

    double sumW = 0.0;
    for (int b = 0; b < nBins; ++b) {
      weights[b] = R_FINITE(midLL[b]) ?
                     std::exp(beta * (midLL[b] - maxLL)) : 0.0;
      sumW += weights[b];
    }
    if (sumW <= 0.0) continue;

    // Sample a bin
    double rnd = R::unif_rand() * sumW;
    int chosenBin = nBins - 1;
    {
      double cum = 0.0;
      for (int b = 0; b < nBins; ++b) {
        cum += weights[b];
        if (rnd < cum) { chosenBin = b; break; }
      }
    }

    // Draw fraction from Beta centred on chosen bin's midpoint
    const double chosenMid = bins.mids[chosenBin];
    const double alphaNew = chosenMid * conc + 1.0;
    const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
    double newF = R::rbeta(alphaNew, betaNew);
    if (newF < 1e-8) newF = 1e-8;
    if (newF > 1.0 - 1e-8) newF = 1.0 - 1e-8;

    // Hastings ratio: as in weighted_branch_scale_impl
    double logHastings = bin_mixture_logdensity(weights, oldF, bins)
                       - bin_mixture_logdensity(weights, newF, bins);

    if (!R_FINITE(logHastings)) continue;

    // Evaluate LL at proposed fraction
    trialAbs[index] = newF * absTotal;
    trialAbs[other] = (1.0 - newF) * absTotal;
    double proposedLL = compute_full_loglik_at(
      *data, *state, state->parent, state->child, trialAbs,
      /*fillCharLLCache=*/false);  // FREEZE-003 scratch eval
    if (!R_FINITE(proposedLL)) continue;

    // MH accept/reject (prior is constant for relBrLengths)
    double logAlpha = beta * (proposedLL - currentLL) + logHastings;

    if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
      state->relBrLengths[index] = newF * relTotal;
      state->relBrLengths[other] = (1.0 - newF) * relTotal;
      absLen[index] = trialAbs[index];
      absLen[other] = trialAbs[other];
      currentLL = proposedLL;
      ++nAccepted;
    } else {
      // Restore trial vector for next iteration
      trialAbs[index] = absLen[index];
      trialAbs[other] = absLen[other];
    }
  }

  // Update state with final likelihood
  if (nAccepted > 0) {
    state->logLik = currentLL;
    // Invalidate partition cache (sweep touched multiple partitions)
    state->partLogLik.clear();
    state->nodeCL.invalidate_all();  // M-143/M-161: branch lengths changed
  }
  return nAccepted > 0;
}


// ---------------------------------------------------------------------------
// weighted_spr_impl  (M-088)
//
// WeightedSPR: Gibbs SPR extended to integrate over branch fractions at each
// candidate reattachment position.  For each candidate, marginalise over B
// branch-fraction bins (same Beta(0.25,0.25) discretisation as M-087) to
// produce a marginal weight M_i.  Self (current topology) included in the
// candidate set.  Sample topology from {self, cand_1, ..., cand_N}
// proportional to marginal weights; if self drawn, return false (no-op).
// For chosen candidate, sample a bin from its conditional distribution, draw
// a fraction from Beta centred on the bin midpoint, construct the final
// proposed topology, and accept/reject via MH.
//
// Hastings ratio: topology selection cancels (both directions enumerate the
// same residual tree), leaving the branch-fraction mixture density and the
// SPR Jacobian; see step 14.
//
// Cost: O(N × B) likelihood evaluations + 1 for the final proposed state.
// ---------------------------------------------------------------------------
static bool weighted_spr_impl(McmcData* data, McmcState* state,
                               double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Eligible prune edges (parent != root)
  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  // 2. Pick random prune edge
  int pickIdx = (int)(R::unif_rand() * (double)eligible.size());
  if (pickIdx >= (int)eligible.size()) pickIdx = (int)eligible.size() - 1;
  const int pruneRow = eligible[pickIdx];
  const int u = state->parent[pruneRow];
  const int v = state->child[pruneRow];

  // 3. Find parentRow (edge into u) and sibRow (u's other child)
  int parentRow = -1, sibRow = -1, sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == u) parentRow = i;
    if (state->parent[i] == u && state->child[i] != v) {
      sibRow = i; sibNode = state->child[i];
    }
  }
  if (parentRow < 0 || sibRow < 0) return false;

  // 4. BFS: mark descendants of v
  std::vector<bool> isDesc(2 * nTip + 2, false);
  isDesc[v] = true;
  if (v > nTip) {
    std::vector<int> queue = {v};
    while (!queue.empty()) {
      int cur = queue.back(); queue.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (state->parent[i] == cur) {
          int c = state->child[i];
          isDesc[c] = true;
          if (c > nTip) queue.push_back(c);
        }
      }
    }
  }

  // 5. Collect valid regraft candidate edges (same filter as GibbsSPR)
  std::vector<int> cands;
  cands.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[state->child[i]]) continue;
    if (state->parent[i] == u || state->child[i] == u) continue;
    cands.push_back(i);
  }
  if (cands.empty()) return false;
  const int nCand = (int)cands.size();

  // 6. Absolute edge lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];
  const double lMerge = absLen[parentRow] + absLen[sibRow];
  const double fOld = (lMerge > 0.0) ? absLen[parentRow] / lMerge : 0.5;

  // 7. Self marginal: evaluate original topology at each bin midpoint
  //    varying the branch fraction at parentRow/sibRow (no topology change,
  //    tree already in preorder → no reorder needed)
  std::vector<double> selfLL(nBins);
  double selfMax = R_NegInf;
  {
    NumericVector trialAbs = clone(absLen);
    for (int b = 0; b < nBins; ++b) {
      trialAbs[parentRow] = bins.mids[b] * lMerge;
      trialAbs[sibRow]    = (1.0 - bins.mids[b]) * lMerge;
      selfLL[b] = compute_full_loglik_at(*data, *state,
                                          state->parent, state->child,
                                          trialAbs,
                                          /*fillCharLLCache=*/false);  // FREEZE-003
      if (R_FINITE(selfLL[b]) && selfLL[b] > selfMax) selfMax = selfLL[b];
    }
  }

  // 8. Candidate marginals: in-place topology + preorder_into (M-109)
  std::vector<std::vector<double>> candLL(nCand,
                                           std::vector<double>(nBins));
  std::vector<double> candMax(nCand, R_NegInf);

  // Working copies: cloned ONCE, reused for all candidates
  IntegerVector workPar = clone(state->parent);
  IntegerVector workCh  = clone(state->child);
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  const int origPar_sibRow = workPar[sibRow];
  const int origCh_parentRow = workCh[parentRow];
  const int origCh_sibRow = workCh[sibRow];
  const double origAbs_parentRow = absLen[parentRow];
  const double origAbs_sibRow = absLen[sibRow];

  for (int ci = 0; ci < nCand; ++ci) {
    const int rr      = cands[ci];
    const double lReg = absLen[rr];
    const int bNode   = workCh[rr];

    // Apply SPR topology in-place (once per candidate)
    workCh[parentRow] = sibNode;
    workCh[rr]        = u;
    workPar[sibRow]   = u;
    workCh[sibRow]    = bNode;
    absLen[parentRow] = lMerge;

    // Evaluate at each bin midpoint (only absLen changes per bin)
    for (int b = 0; b < nBins; ++b) {
      absLen[rr]     = bins.mids[b] * lReg;
      absLen[sibRow] = (1.0 - bins.mids[b]) * lReg;

      preorder_into(workPar, workCh, absLen, nTip,
                    INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
      candLL[ci][b] = compute_full_loglik_at(*data, *state,
                                              ordPar, ordCh, ordAbs,
                                              /*fillCharLLCache=*/false);  // FREEZE-003
      if (R_FINITE(candLL[ci][b]) && candLL[ci][b] > candMax[ci])
        candMax[ci] = candLL[ci][b];
    }

    // Restore
    workCh[parentRow] = origCh_parentRow;
    workCh[rr]        = bNode;
    workPar[sibRow]   = origPar_sibRow;
    workCh[sibRow]    = origCh_sibRow;
    absLen[parentRow] = origAbs_parentRow;
    absLen[rr]        = lReg;
    absLen[sibRow]    = origAbs_sibRow;
  }

  // 9. Compute marginal weights with a single global offset for stability
  double globalMax = selfMax;
  for (int ci = 0; ci < nCand; ++ci)
    if (candMax[ci] > globalMax) globalMax = candMax[ci];
  if (!R_FINITE(globalMax)) return false;

  // Self marginal
  double mSelf = 0.0;
  std::vector<double> selfW(nBins);
  for (int b = 0; b < nBins; ++b) {
    selfW[b] = R_FINITE(selfLL[b]) ?
                 std::exp(beta * (selfLL[b] - globalMax)) : 0.0;
    mSelf += selfW[b];
  }

  // Candidate marginals
  std::vector<double> mCand(nCand);
  std::vector<std::vector<double>> candW(nCand,
                                          std::vector<double>(nBins));
  double sumM = mSelf;
  for (int ci = 0; ci < nCand; ++ci) {
    mCand[ci] = 0.0;
    for (int b = 0; b < nBins; ++b) {
      candW[ci][b] = R_FINITE(candLL[ci][b]) ?
                       std::exp(beta * (candLL[ci][b] - globalMax)) : 0.0;
      mCand[ci] += candW[ci][b];
    }
    sumM += mCand[ci];
  }
  if (sumM <= 0.0) return false;

  // 10. Sample topology: self or candidate
  double rnd = R::unif_rand() * sumM;
  if (rnd < mSelf) return false;  // self-draw → no-op
  rnd -= mSelf;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < mCand[ci]) { chosen = ci; break; }
    rnd -= mCand[ci];
  }

  // 11. Sample bin within chosen candidate
  int chosenBin = nBins - 1;
  {
    double rndBin = R::unif_rand() * mCand[chosen];
    double cum = 0.0;
    for (int b = 0; b < nBins; ++b) {
      cum += candW[chosen][b];
      if (rndBin < cum) { chosenBin = b; break; }
    }
  }

  // 12. Draw fraction from Beta centred on chosen bin's midpoint
  const double conc = bins.concentration;
  const double chosenMid = bins.mids[chosenBin];
  const double alphaNew = chosenMid * conc + 1.0;
  const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
  double fNew = R::rbeta(alphaNew, betaNew);
  if (fNew < 1e-8) fNew = 1e-8;
  if (fNew > 1.0 - 1e-8) fNew = 1.0 - 1e-8;

  // 13. Construct final proposed topology with fNew using working copies
  //     (state untouched until acceptance confirmed)
  const int rr      = cands[chosen];
  const double lReg = absLen[rr];
  const int bNode   = workCh[rr];

  workCh[parentRow]  = sibNode;  absLen[parentRow] = lMerge;
  workCh[rr]         = u;        absLen[rr]        = fNew * lReg;
  workPar[sibRow]    = u;        workCh[sibRow]    = bNode;
  absLen[sibRow]     = (1.0 - fNew) * lReg;

  auto po = TreeTools::preorder_weighted_impl(workPar, workCh, absLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbsFinal = po.second;

  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbsFinal,
                                             /*fillCharLLCache=*/false);  // FREEZE-003
  if (!R_FINITE(newLogLik)) return false;

  // 14. Bin-mixture proposal density + the SPR merge/split Jacobian
  //     lReg / lMerge, the same term as spr_proposal_impl
  //     (proposals.cpp:161) and tbr_proposal_impl (tree_moves.cpp:492).
  //     The selection normaliser cancels: both directions prune to the same
  //     residual tree and enumerate its edges against the same bins.
  //     Reaching x from y is a CANDIDATE draw onto the merged edge, not y's
  //     self-draw, which is why selfW carries the reverse weight.  Derivation
  //     term by term: dev/red-team/proofs/weighted-spr-hastings.md.
  double logHR = bin_mixture_logdensity(selfW, fOld, bins)
               - bin_mixture_logdensity(candW[chosen], fNew, bins)
               + std::log(lReg) - std::log(lMerge);
  if (!R_FINITE(logHR)) return false;

  // 15. Prior at proposed state
  NumericVector propRelBr(nEdge);
  for (int k = 0; k < nEdge; ++k)
    propRelBr[k] = ordAbsFinal[k] / state->treeLength;

  double newLogPrior = compute_log_prior_at(*data, *state, propRelBr);
  if (!R_FINITE(newLogPrior)) return false;

  // 16. MH acceptance
  double logAlpha = beta * (newLogLik - state->logLik)
                  + (newLogPrior - state->logPrior) + logHR;
  if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = propRelBr[k];
    }
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik.clear();
    state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
    return true;
  }
  return false;
}


// ---------------------------------------------------------------------------
// weighted_subtree_swap_impl  (M-089)
//
// WeightedSubtreeSwap: GibbsSubtreeSwap extended to integrate over branch
// fractions at each candidate swap partner.  For each candidate B_i, the
// total branch length (brA + brB_i) is held fixed and redistributed across
// B bins.  Self (current topology) included as a point weight.  Sample
// partner from {self, cand_0, ..., cand_N}, then sample bin and fraction
// for the chosen partner, and accept/reject via MH.
//
// Neither the selection normaliser nor the bin weights cancel: the anchored
// neighbourhood differs between the two endpoints (as for gibbs_subtree_swap),
// and the reverse bins are scored on the reverse topology.  The reverse
// neighbourhood is therefore enumerated in full, doubling the cost.
// Derivation: dev/red-team/proofs/weighted-subtree-swap-hastings.md.
//
// Cost: O(N × B) likelihood evaluations in each direction.
// ---------------------------------------------------------------------------

// (Uses find_child_row_gibbs defined above)

// The anchored neighbourhood of `nodeA`: log selection weight of each
// (partner, bin) configuration, and the log normaliser including staying put.
// Rows are located by child, so any edge order is accepted.
struct SwapNeighbourhood {
  std::vector<int> partners;
  std::vector<int> rowBs;
  std::vector<double> totals;               // brA + brB_i
  std::vector<std::vector<double>> logW;    // [partner][bin]
  double logZ = R_NegInf;
};

static void weighted_swap_neighbourhood(
    McmcData* data, McmcState* state,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& absLenIn, int nodeA, double logSelf, double beta,
    SwapNeighbourhood& nb) {
  const int nEdge = parent.size();
  const int nTip  = data->nTip;
  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  nb.partners = get_valid_swap_partners_impl(parent, child, nTip, nodeA);
  const int nPart = (int)nb.partners.size();
  nb.rowBs.assign(nPart, -1);
  nb.totals.assign(nPart, 0.0);
  nb.logW.assign(nPart, std::vector<double>(nBins, R_NegInf));
  nb.logZ = R_NegInf;

  const int rowA = find_child_row_gibbs(child, nodeA);
  if (rowA < 0) return;

  // Working copies; child is unchanged by a subtree swap (M-109)
  IntegerVector workPar = clone(parent);
  NumericVector absLen  = clone(absLenIn);
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  const int origParA = workPar[rowA];
  const double origAbsA = absLen[rowA];
  double logMax = logSelf;

  for (int pi = 0; pi < nPart; ++pi) {
    const int rowB = find_child_row_gibbs(child, nb.partners[pi]);
    if (rowB < 0) continue;
    nb.rowBs[pi]  = rowB;
    nb.totals[pi] = absLen[rowA] + absLen[rowB];

    const int origParB = workPar[rowB];
    const double origAbsB = absLen[rowB];
    workPar[rowA] = origParB;
    workPar[rowB] = origParA;

    for (int b = 0; b < nBins; ++b) {
      absLen[rowA] = bins.mids[b] * nb.totals[pi];
      absLen[rowB] = (1.0 - bins.mids[b]) * nb.totals[pi];
      preorder_into(workPar, child, absLen, nTip,
                    INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
      const double ll = compute_full_loglik_at(*data, *state,
                                               ordPar, ordCh, ordAbs,
                                               /*fillCharLLCache=*/false);  // FREEZE-003
      if (R_FINITE(ll)) {
        nb.logW[pi][b] = beta * ll;
        if (nb.logW[pi][b] > logMax) logMax = nb.logW[pi][b];
      }
    }

    workPar[rowA] = origParA;
    workPar[rowB] = origParB;
    absLen[rowA]  = origAbsA;
    absLen[rowB]  = origAbsB;
  }

  if (!R_FINITE(logMax)) return;
  double sum = R_FINITE(logSelf) ? std::exp(logSelf - logMax) : 0.0;
  for (int pi = 0; pi < nPart; ++pi)
    for (int b = 0; b < nBins; ++b)
      sum += std::exp(nb.logW[pi][b] - logMax);
  nb.logZ = logMax + std::log(sum);
}

// Log density of fraction `f` under the bin mixture with log weights `logW`.
static double log_bin_mixture(const std::vector<double>& logW, double f,
                              const BranchBins& bins) {
  double logMax = R_NegInf;
  for (double lw : logW) if (lw > logMax) logMax = lw;
  if (!R_FINITE(logMax)) return R_NegInf;
  std::vector<double> w(logW.size());
  for (size_t b = 0; b < logW.size(); ++b) w[b] = std::exp(logW[b] - logMax);
  return logMax + bin_mixture_logdensity(w, f, bins);
}

static bool weighted_subtree_swap_impl(McmcData* data, McmcState* state,
                                        double beta) {
  const int nEdge = state->parent.size();
  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Pick a random node (any edge child)
  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];
  const int rowA = find_child_row_gibbs(state->child, nodeA);
  if (rowA < 0) return false;

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // 2. Forward neighbourhood
  const double logSelf = beta * state->logLik;
  SwapNeighbourhood fwd;
  weighted_swap_neighbourhood(data, state, state->parent, state->child,
                              absLen, nodeA, logSelf, beta, fwd);
  const int nPart = (int)fwd.partners.size();
  if (nPart == 0 || !R_FINITE(fwd.logZ)) return false;

  // 3. Sample (partner, bin) jointly; self-draw → no-op
  double rnd = R::unif_rand() - std::exp(logSelf - fwd.logZ);
  if (rnd < 0.0) return false;
  int chosen = -1, chosenBin = -1;
  for (int pi = 0; pi < nPart && chosen < 0; ++pi) {
    for (int b = 0; b < nBins; ++b) {
      rnd -= std::exp(fwd.logW[pi][b] - fwd.logZ);
      if (rnd < 0.0) { chosen = pi; chosenBin = b; break; }
    }
  }
  if (chosen < 0) return false;  // rounding at the top of [0, 1)
  const int rowB  = fwd.rowBs[chosen];
  const int nodeB = fwd.partners[chosen];

  // 4. Draw fraction from Beta centred on chosen bin's midpoint
  const double conc = bins.concentration;
  const double chosenMid = bins.mids[chosenBin];
  double fNew = R::rbeta(chosenMid * conc + 1.0,
                         (1.0 - chosenMid) * conc + 1.0);
  if (fNew < 1e-8) fNew = 1e-8;
  if (fNew > 1.0 - 1e-8) fNew = 1.0 - 1e-8;

  // 5. Construct the proposed tree
  const double tot  = fwd.totals[chosen];
  const double fOld = absLen[rowA] / tot;
  IntegerVector workPar = clone(state->parent);
  workPar[rowA] = state->parent[rowB];
  workPar[rowB] = state->parent[rowA];
  NumericVector propAbs = clone(absLen);
  propAbs[rowA] = fNew * tot;
  propAbs[rowB] = (1.0 - fNew) * tot;

  auto po = TreeTools::preorder_weighted_impl(workPar, state->child, propAbs);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbsFinal = po.second;
  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbsFinal,
                                             /*fillCharLLCache=*/false);  // FREEZE-003
  if (!R_FINITE(newLogLik)) return false;

  // 6. Reverse neighbourhood: the same anchor returns to the current state by
  //    swapping back with nodeB at fraction fOld.  Enumerated on the
  //    un-reordered proposal, because preordering renumbers internal nodes.
  SwapNeighbourhood rev;
  weighted_swap_neighbourhood(data, state, workPar, state->child, propAbs,
                              nodeA, beta * newLogLik, beta, rev);
  int back = -1;
  for (size_t pi = 0; pi < rev.partners.size(); ++pi)
    if (rev.partners[pi] == nodeB) { back = (int)pi; break; }
  if (back < 0 || !R_FINITE(rev.logZ)) return false;

  // 7. Hastings ratio
  const double logHR = log_bin_mixture(rev.logW[back], fOld, bins) - rev.logZ
                     - log_bin_mixture(fwd.logW[chosen], fNew, bins) + fwd.logZ;

  // 8. Prior at proposed state
  NumericVector propRelBr(nEdge);
  for (int k = 0; k < nEdge; ++k)
    propRelBr[k] = ordAbsFinal[k] / state->treeLength;

  double newLogPrior = compute_log_prior_at(*data, *state, propRelBr);
  if (!R_FINITE(newLogPrior)) return false;

  // 9. MH acceptance
  double logAlpha = beta * (newLogLik - state->logLik)
                  + (newLogPrior - state->logPrior) + logHR;
  if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = propRelBr[k];
    }
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik.clear();
    state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
    return true;
  }
  return false;
}


// ---------------------------------------------------------------------------
// Slice sampler for scalar parameters.
// paramIdx: 0=treeLength, 1=rateLoss, 2=rateLogSd, 3=rateNeo, 4=betaScale
// ---------------------------------------------------------------------------

static double get_scalar(const McmcState* state, int paramIdx) {
  switch (paramIdx) {
    case 0: return state->treeLength;
    case 1: return state->rateLoss;
    case 2: return state->rateLogSd;
    case 3: return state->rateNeo;
    case 4: return state->betaScale;
    default: return 0.0;
  }
}

static void set_scalar(McmcState* state, int paramIdx, double val) {
  switch (paramIdx) {
    case 0: state->treeLength = val; break;
    case 1: state->rateLoss = val; break;
    case 2:
      state->rateLogSd = val;
      // Lockstep invariant: classRateLogSd[0] == rateLogSd when partitioned.
      // The partitioned prior (cpp_log_prior_partitioned) skips c==0 on the
      // assumption this invariant holds; without the mirror, classRateLogSd[0]
      // would be unconstrained and the likelihood (which reads classRateLogSd)
      // would diverge from the prior path.
      if (state->usePartitioned && state->classRateLogSd.size() > 0)
        state->classRateLogSd[0] = val;
      break;
    case 3: state->rateNeo = val; break;
    case 4: state->betaScale = val; break;
  }
}

// Evaluate beta * logLik + logPrior for current state, using partial cache
// when the parameter only affects a subset of partitions.
static double eval_slice_target(McmcData* data, McmcState* state,
                                int paramIdx, double beta) {
  double logPrior = compute_log_prior(*data, *state);
  if (!R_FINITE(logPrior)) return R_NegInf;

  double logLik;
  int nEdge = state->parent.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  // !marginalK mirrors do_move_impl: partLogLik is fixed-kPrime, so the partial
  // update below would evaluate the wrong target under marginal_k.
  bool hasPLC = !state->partLogLik.empty() && !data->marginalK;
  ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

  // rate_loss (1): only neomorphic partitions change. (rate_neo / paramIdx 3
  // no longer qualifies — audit Issue 1: rate_neo now shifts transScale too,
  // so trans partitions are stale and we must fall through to full eval.)
  if (hasPLC && paramIdx == 1) {
    logLik = state->logLik;
    for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
      int pi = data->neoPartIndices[ni];
      double oldPart = state->partLogLik[pi];
      double newPart = cpp_partition_log_likelihood(
        *data, pi, state->parent, state->child, edgeLen,
        state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale, wsPtr);
      logLik += (newPart - oldPart);
    }
  } else if (data->marginalK) {
    // MARGINAL-K-SLICE-001: under marginal_k the slice target must be the
    // marginal-over-k likelihood, NOT the fixed-kPrime evaluator (state->kPrime
    // is pinned to kObs in marginal mode). Reset the charLL cache first so the
    // marginal evaluator does the full rebuild against the proposed scalar
    // rather than taking the p-only fast-path on stale rawLL. Routed through
    // compute_full_loglik_at — the same dispatcher the do_move_impl MH path
    // uses, which is empirically correct under marginal_k.
    state->charLLCacheReady = false;
    logLik = compute_full_loglik_at(
        *data, *state, state->parent, state->child, edgeLen);
  } else {
    logLik = state->usePartitioned
      ? cpp_log_likelihood_partitioned(
          *data, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss,
          state->classRateLogSd, state->classRate,
          state->etaNeo, state->betaScale, wsPtr)
      : cpp_log_likelihood(
          *data, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateLogSd,
          state->rateNeo, state->betaScale, wsPtr);
  }
  if (!R_FINITE(logLik)) return R_NegInf;
  return beta * logLik + logPrior;
}

// Univariate stepping-out slice sampler.
// Returns true on success (always, barring degenerate cases).
// Updates state in place, including logLik, logPrior, and partLogLik cache.
static bool slice_scalar_impl(McmcData* data, McmcState* state,
                               int paramIdx, double width,
                               double beta, int maxSteps = 10,
                               int* nExpansionsOut = nullptr) {
  double x0 = get_scalar(state, paramIdx);

  // Work on log scale: u = log(x).
  // Target includes Jacobian: log f(u) = beta*logLik + logPrior + u.
  double u0 = std::log(x0);
  double logY0 = beta * state->logLik + state->logPrior + u0;

  // Slice height
  double logZ = logY0 + std::log(R::unif_rand());

  // Stepping out on log scale (count expansions for width adaptation)
  int nExp = 0;
  double L = u0 - width * R::unif_rand();
  double R_bound = L + width;

  for (int j = 0; j < maxSteps; ++j) {
    set_scalar(state, paramIdx, std::exp(L));
    if (eval_slice_target(data, state, paramIdx, beta) + L <= logZ) break;
    L -= width;
    ++nExp;
  }
  for (int j = 0; j < maxSteps; ++j) {
    set_scalar(state, paramIdx, std::exp(R_bound));
    if (eval_slice_target(data, state, paramIdx, beta) + R_bound <= logZ)
      break;
    R_bound += width;
    ++nExp;
  }
  if (nExpansionsOut) *nExpansionsOut = nExp;

  // Shrink in on log scale
  for (int iter = 0; iter < 100; ++iter) {
    double u1 = L + R::unif_rand() * (R_bound - L);
    double x1 = std::exp(u1);
    set_scalar(state, paramIdx, x1);
    double logTarget1 = eval_slice_target(data, state, paramIdx, beta) + u1;
    if (logTarget1 >= logZ) {
      // Accept — recompute and cache logLik / logPrior / partLogLik
      state->logPrior = compute_log_prior(*data, *state);

      int nEdge = state->parent.size();
      NumericVector edgeLen(nEdge);
      for (int i = 0; i < nEdge; ++i)
        edgeLen[i] = state->treeLength * state->relBrLengths[i];
      ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

      bool hasPLC = !state->partLogLik.empty() && !data->marginalK;
      if (hasPLC && paramIdx == 1) {
        // rate_loss: only neo partitions change. Partial cache update.
        // (paramIdx == 3 / rate_neo intentionally excluded — audit Issue 1:
        // rate_neo now shifts transScale and so makes every partition stale.)
        for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
          int pi = data->neoPartIndices[ni];
          state->partLogLik[pi] = cpp_partition_log_likelihood(
            *data, pi, state->parent, state->child, edgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd,
            state->rateNeo, state->betaScale, wsPtr);
        }
        state->logLik = 0.0;
        for (size_t pi = 0; pi < state->partLogLik.size(); ++pi)
          state->logLik += state->partLogLik[pi];
      } else if (data->marginalK) {
        // MARGINAL-K-SLICE-001: recompute the accepted logLik via the
        // marginal-over-k evaluator (not the fixed-kPrime path). Reset the
        // charLL cache so the rebuild honours the accepted scalar; the rebuild
        // re-sets charLLCacheReady=true, leaving a coherent cache for the next
        // p-move fast-path (whose rawLL is independent of p).
        state->charLLCacheReady = false;
        state->logLik = compute_full_loglik_at(
            *data, *state, state->parent, state->child, edgeLen);
        // Invalidate partition cache (full recompute was done)
        state->partLogLik.clear();
      } else {
        state->logLik = state->usePartitioned
          ? cpp_log_likelihood_partitioned(
              *data, state->parent, state->child, edgeLen,
              state->kPrime, state->rateLoss,
              state->classRateLogSd, state->classRate,
              state->etaNeo, state->betaScale, wsPtr)
          : cpp_log_likelihood(
              *data, state->parent, state->child, edgeLen,
              state->kPrime, state->rateLoss, state->rateLogSd,
              state->rateNeo, state->betaScale, wsPtr);
        // Invalidate partition cache (full recompute was done)
        state->partLogLik.clear();
      }
      // M-145/M-161: Invalidate node CL cache — slice changed a model
      // parameter.  Granular: only invalidate affected units.
      switch (paramIdx) {
        case 1:  // rate_loss: only neomorphic units (enters mkn stationary
                 // freqs + Q-matrix; trans/known unchanged)
          state->nodeCL.invalidate_neo_cls();
          break;
        case 3:  // rate_neo: audit Issue 1 — RB-style partition-rate
                 // normalisation makes rateNeo affect BOTH neo and trans
                 // unit rateScales, so all CLs must be invalidated.
          state->nodeCL.invalidate_all_cls();
          break;
        case 2:  // rateLogSd: ACRV rates change, all units
          state->nodeCL.invalidate_all_cls();
          break;
        default:  // tree_length (0), beta_scale (4)
          state->nodeCL.invalidate_all();
          break;
      }
      return true;
    }
    // Shrink bracket on log scale
    if (u1 < u0) L = u1; else R_bound = u1;
  }

  // Fallback: restore original value
  set_scalar(state, paramIdx, x0);
  return false;
}

// ---------------------------------------------------------------------------
// Prior-only slice sampler for Beta-Geometric hyperparameters
// — reparameterised in (s, r) coordinates.
//
// s = log(α + β)            (size / concentration; unconstrained)
// r = log(α / β) = logit(α / (α + β))   (shape; unconstrained)
//
// α = e^s · σ(r),  β = e^s · (1 − σ(r)),  where σ(r) = 1 / (1 + e^{-r}).
//
// The α-axis / β-axis univariate slice (the pre-2026-05-28 implementation)
// was tightly coupled to the (log α, log β) posterior ridge of the BG model
// — neither axis-aligned slice nor axis-aligned Bactrian could traverse it.
// (s, r) decorrelates the BG posterior across the whole parameter space.
//
// paramCode 0 = sample s with r held fixed
// paramCode 1 = sample r with s held fixed
//
// Target density on (s, r) (change of variables from (α, β) ~ π):
//   log π_{s,r}(s, r) = log π_{α,β}(α(s,r), β(s,r)) + log|J|
// where |J| = |∂(α, β) / ∂(s, r)| = α · β.  So target adds
// `+ log α + log β` to compute_log_prior() at the candidate point.
// ---------------------------------------------------------------------------
static inline double bg_sigma_stable(double r) {
  // numerically stable σ(r) = 1 / (1 + e^{-r})
  if (r >= 0.0) return 1.0 / (1.0 + std::exp(-r));
  double e = std::exp(r);
  return e / (1.0 + e);
}

static bool slice_kprime_hyper_impl(McmcData* data, McmcState* state,
                                     int paramCode, double width,
                                     int maxSteps = 10,
                                     int* nExpansionsOut = nullptr) {
  double alpha0 = state->kprimeAlpha;
  double beta0  = state->kprimeBeta;
  if (alpha0 <= 0.0 || beta0 <= 0.0) return false;

  // Current (s, r) and the held-fixed coordinate
  double s0 = std::log(alpha0 + beta0);
  double r0 = std::log(alpha0 / beta0);
  double s_fixed = s0;
  double ratio_fixed = alpha0 / (alpha0 + beta0);  // = σ(r0)

  // Helper: derive (α, β) from a candidate value of the active coord
  auto srToAlphaBeta = [&](double v, double& aOut, double& bOut) -> bool {
    if (paramCode == 0) {
      // v = candidate s; r fixed → ratio_fixed = σ(r) held constant
      double ev = std::exp(v);
      if (!R_FINITE(ev) || ev <= 0.0) return false;
      aOut = ev * ratio_fixed;
      bOut = ev * (1.0 - ratio_fixed);
    } else {
      // v = candidate r; s fixed
      double sig = bg_sigma_stable(v);
      double exp_s = std::exp(s_fixed);
      if (!R_FINITE(exp_s) || exp_s <= 0.0) return false;
      aOut = exp_s * sig;
      bOut = exp_s * (1.0 - sig);
    }
    return R_FINITE(aOut) && R_FINITE(bOut) && aOut > 0.0 && bOut > 0.0;
  };

  // Target at current point: logPrior(α0, β0) + log|J| = logPrior + log α + log β
  double logY0 = state->logPrior + std::log(alpha0) + std::log(beta0);
  double logZ = logY0 + std::log(R::unif_rand());

  // Evaluate target at a candidate value of the active coord
  auto evalTarget = [&](double v) -> double {
    double aCand, bCand;
    if (!srToAlphaBeta(v, aCand, bCand)) return R_NegInf;
    double oldA = state->kprimeAlpha, oldB = state->kprimeBeta;
    state->kprimeAlpha = aCand;
    state->kprimeBeta  = bCand;
    double lp = compute_log_prior(*data, *state);
    state->kprimeAlpha = oldA;
    state->kprimeBeta  = oldB;
    if (!R_FINITE(lp)) return R_NegInf;
    return lp + std::log(aCand) + std::log(bCand);
  };

  double v0 = (paramCode == 0) ? s0 : r0;

  // Stepping out
  int nExp = 0;
  double L = v0 - width * R::unif_rand();
  double R_bound = L + width;
  for (int j = 0; j < maxSteps; ++j) {
    if (evalTarget(L) <= logZ) break;
    L -= width;
    ++nExp;
  }
  for (int j = 0; j < maxSteps; ++j) {
    if (evalTarget(R_bound) <= logZ) break;
    R_bound += width;
    ++nExp;
  }
  if (nExpansionsOut) *nExpansionsOut = nExp;

  // Shrink in
  for (int iter = 0; iter < 100; ++iter) {
    double v1 = L + R::unif_rand() * (R_bound - L);
    double logTarget1 = evalTarget(v1);
    if (logTarget1 >= logZ) {
      double a1, b1;
      if (!srToAlphaBeta(v1, a1, b1)) return false;
      state->kprimeAlpha = a1;
      state->kprimeBeta  = b1;
      state->logPrior = compute_log_prior(*data, *state);
      return true;
    }
    if (v1 < v0) L = v1; else R_bound = v1;
  }

  // Fallback: no change
  return false;
}

// ---------------------------------------------------------------------------
// Parsimony-guided SPR (pSPR) — M-119
//
// Like regular SPR but weights candidate regraft edges by their Fitch
// parsimony score: w_i = exp(-alpha * (score_i - score_min)).
// Better-scoring positions are proposed more often.
//
// The Hastings ratio correction for asymmetric proposal:
//   log(H) = log(l_regraft) - log(l_merge)       [branch-length Jacobian]
//          + alpha * (score_new - score_old)       [parsimony bias correction]
//          - log(w_chosen / sum_w_forward)         [forward proposal]
//          + log(w_reverse / sum_w_reverse)        [reverse proposal]
//
// Since the residual tree is the same for both directions, the candidate
// set and parsimony scores are identical, simplifying to:
//   log(H) = log(l_regraft) - log(l_merge) + log(w_orig) - log(w_chosen)
//          = log(l_regraft) - log(l_merge) - alpha * (score_orig - score_chosen)
// ---------------------------------------------------------------------------

static constexpr double PSPR_ALPHA = 0.1;  // parsimony bias strength

static List pspr_proposal_impl(
    const IntegerVector& stateParent, const IntegerVector& stateChild,
    int nTip, double treeLength,
    const NumericVector& relBrLengths,
    const McmcData* data)
{
  const int nEdge = stateParent.size();
  const int root  = nTip + 1;

  // Writable copies for in-place Fitch scoring
  IntegerVector workParent = clone(stateParent);
  IntegerVector workChild  = clone(stateChild);
  int* wp = INTEGER(workParent);
  int* wc = INTEGER(workChild);

  // 1. Eligible prune edges (parent != root)
  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (stateParent[i] != root) eligible.push_back(i);
  if (eligible.empty())
    return List::create(_["logHastings"] = R_NegInf);

  int pick = (int)(R::unif_rand() * (double)eligible.size());
  if (pick >= (int)eligible.size()) pick = (int)eligible.size() - 1;
  const int pruneRow = eligible[pick];
  const int u = stateParent[pruneRow];
  const int v = stateChild[pruneRow];

  // 2. Find parentRow (edge → u) and sibRow (u → sibling of v)
  int parentRow = -1, sibRow = -1, sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (stateChild[i] == u) parentRow = i;
    if (stateParent[i] == u && stateChild[i] != v) {
      sibRow = i; sibNode = stateChild[i];
    }
  }
  if (parentRow < 0 || sibRow < 0)
    return List::create(_["logHastings"] = R_NegInf);

  // 3. BFS: mark descendants of v
  std::vector<bool> isDesc(2 * nTip + 2, false);
  isDesc[v] = true;
  if (v > nTip) {
    std::vector<int> queue = {v};
    while (!queue.empty()) {
      int cur = queue.back(); queue.pop_back();
      for (int i = 0; i < nEdge; ++i) {
        if (stateParent[i] == cur) {
          int c = stateChild[i];
          isDesc[c] = true;
          if (c > nTip) queue.push_back(c);
        }
      }
    }
  }

  // 4. Collect candidate regraft edges
  std::vector<int> candidates;
  candidates.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[stateChild[i]]) continue;
    if (stateParent[i] == u || stateChild[i] == u) continue;
    candidates.push_back(i);
  }
  if (candidates.empty())
    return List::create(_["logHastings"] = R_NegInf);
  const int nCand = (int)candidates.size();

  // 5. Build partition list for Fitch scoring
  std::vector<std::pair<IntegerMatrix, int>> fitchParts;
  for (const auto& part : data->parts) {
    int kEff = (part.type == 0) ? 2 :
               (part.type == 2) ? part.k : 0;
    if (kEff == 0) {
      // Transformational: use max observed k across characters
      int maxK = 2;
      for (int c = 0; c < part.tipStates.ncol(); ++c) {
        for (int t = 0; t < nTip; ++t) {
          int s = part.tipStates(t, c);
          if (s >= maxK) maxK = s + 1;
        }
      }
      kEff = maxK;
    }
    fitchParts.push_back({part.tipStates, kEff});
  }

  // 6. Score all candidates using Fitch parsimony
  std::vector<int> scores;
  fitch_score_candidates(wp, wc, nEdge, nTip, fitchParts,
                         pruneRow, parentRow, sibRow,
                         u, v, sibNode, candidates, scores);

  // Also score the original tree to get the "original position" score.
  // The original position in the residual tree is the merged edge at
  // parentRow. We identify it by checking which candidate, if regrafted,
  // would reproduce a topology equivalent to the original.
  // Actually: the original tree score = fitch_score_all of unmodified tree.
  int scoreOrig = fitch_score_all(wp, wc, nEdge, nTip, fitchParts);

  // 7. Compute weights: w_i = exp(-alpha * (score_i - score_min))
  int minScore = scoreOrig;
  for (int ci = 0; ci < nCand; ++ci)
    if (scores[ci] < minScore) minScore = scores[ci];

  std::vector<double> logW(nCand);
  std::vector<double> w(nCand);
  double sumW = 0.0;
  for (int ci = 0; ci < nCand; ++ci) {
    logW[ci] = -PSPR_ALPHA * (double)(scores[ci] - minScore);
    w[ci] = std::exp(logW[ci]);
    sumW += w[ci];
  }

  // 8. Sample from weighted distribution
  double rnd = R::unif_rand() * sumW;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < w[ci]) { chosen = ci; break; }
    rnd -= w[ci];
  }

  // 9. Apply the chosen SPR
  const int regraftRow = candidates[chosen];
  const int b = stateChild[regraftRow];
  const double tau = R::unif_rand();

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = treeLength * relBrLengths[i];

  const double lRegraft = absLen[regraftRow];
  const double lMerge   = absLen[parentRow] + absLen[sibRow];

  IntegerVector newParent = clone(stateParent);
  IntegerVector newChild  = clone(stateChild);
  NumericVector newAbsLen = clone(absLen);

  // Suppress u
  newChild[parentRow]  = sibNode;
  newAbsLen[parentRow] = lMerge;
  // Insert u on regraft edge
  newChild[regraftRow]  = u;
  newAbsLen[regraftRow] = tau * lRegraft;
  // Reuse sibRow for u → b
  newParent[sibRow] = u;
  newChild[sibRow]  = b;
  newAbsLen[sibRow] = (1.0 - tau) * lRegraft;

  // Canonical preorder reordering
  auto po = TreeTools::preorder_weighted_impl(newParent, newChild, newAbsLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector ordParent(nEdge), ordChild(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    ordParent[i] = ordEdge(i, 0);
    ordChild[i]  = ordEdge(i, 1);
  }
  NumericVector orderedRelBr = ordAbs / treeLength;

  // 10. Hastings ratio: branch-length Jacobian + parsimony bias correction
  //
  // Forward candidate set F excludes the original-position edge; reverse
  // candidate set Rev excludes the chosen-position edge.  The normalization
  // constants differ: sumW_fwd = sumW, sumW_rev = sumW + wOrig - w[chosen].
  double logH_brlen = std::log(lRegraft) - std::log(lMerge);
  double wOrig   = std::exp(-PSPR_ALPHA * (double)(scoreOrig - minScore));
  double sumWRev = sumW + wOrig - w[chosen];
  double logH_pars = std::log(wOrig) - std::log(w[chosen])
                   + std::log(sumW) - std::log(sumWRev);
  double logHastings = logH_brlen + logH_pars;

  return List::create(_["parent"] = ordParent,
                      _["child"] = ordChild,
                      _["rel_br_lengths"] = orderedRelBr,
                      _["logHastings"] = logHastings);
}


// ---------------------------------------------------------------------------
// compute_per_kprime_log_lik — case-25 phase-1 helper.
//
// Spine of the marginal-k landing (dev/notes/2026-05-28-marginal-k-plan.md §3
// architectural insight, §14 PR-A). Extracted verbatim from the phase-1 body
// of gibbs_kprime_sweep_impl: batched per-(char, k') likelihood evaluation
// with M-155 batching, M-164 prior-ceiling cutoff, and M-172 pattern dedup.
//
// Populates `out.charLogW[ti * kMaxKprimeCand + ko]` with
//   β · LL(y_i | tree, μ, k = kObs_i + ko) + logPrior_k_under_arm
// for ko ∈ 0..(charNCand[ti] - 1). The prior and β are baked into the
// weights here so the existing sampled-k phase-2 categorical sampler can
// consume `out` unchanged. PR-B will add a sibling helper that returns raw
// LLs (no prior, no β) for the marginal-k logSumExp evaluator.
//
// No RNG: caller's `state->kPrime` is read but not written; `state->gibbsWs`
// is grown if needed (deterministic allocation).
// ---------------------------------------------------------------------------

void compute_per_kprime_log_lik(
    McmcData* data, McmcState* state, double beta,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    NumericVector acrvRates,
    KprimeCharWeights& out) {
  // MARGINAL-K-FREEZE-003 (Bug B): topology comes from the CALLER (parent,
  // child), NOT from state->parent/child. A topology-changing MH move
  // evaluates its proposal BEFORE committing, so state->parent still holds the
  // OLD tree at eval time; reading state here evaluated the wrong topology and
  // froze / mis-sampled marginal_k under free topology.
  int nTrans = (int)data->transIdxGlobal.size();
  out.resize(nTrans);
  if (nTrans == 0) return;

  // Partition-rate normalisation (issue #25). Every partition this helper
  // prunes is transformational, so cpp_partition_log_likelihood would scale its
  // edges by compute_partition_scales(...).trans; pruning at the caller's raw
  // lengths made marginal_k target a different posterior from sampled_k, and
  // made the case-25 Gibbs sweep sample the unscaled conditional and always
  // accept it. Applied here rather than at the two call sites so they cannot
  // diverge again. Both scales are 1 when nNeo == 0 or nTrans == 0.
  const double transScale =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans).trans;
  if (transScale != 1.0) {
    NumericVector scaledEdge(edgeLen.size());
    for (int i = 0; i < edgeLen.size(); ++i)
      scaledEdge[i] = edgeLen[i] * transScale;
    edgeLen = scaledEdge;   // rebinds the local only; caller's vector untouched
  }

  // Cache for constant-site probability by kStates (shared across characters)
  std::vector<double> cspCache;
  int cspCacheSize = 0;
  auto getCSP = [&](int kStates) -> double {
    if (data->codingType == 0) return 0.0;
    if (kStates >= cspCacheSize) {
      int newSize = kStates + 20;
      cspCache.resize(newSize, -1.0);
      cspCacheSize = newSize;
    }
    if (cspCache[kStates] < 0.0) {
      cspCache[kStates] = const_site_prob_for_k(
        *data, parent, child, edgeLen,
        kStates, state->betaScale, acrvRates);
    }
    return cspCache[kStates];
  };

  // Prior components (fixed during sweep)
  bool isBetaGeometric = data->kPriorBetaGeometric;
  bool isEmpGeom       = data->kPriorEmpiricalGeometric;
  bool isGeometric = !data->kPriorLogseries && !isBetaGeometric && !isEmpGeom;
  double logP = 0.0, log1mP = 0.0;
  double lsLogC = 0.0, lsLogNorm = 0.0;
  if (isGeometric || isEmpGeom) {
    logP   = std::log(state->p);
    log1mP = std::log1p(-state->p);
  } else if (!isBetaGeometric) {
    double c = data->kprimeLogseriesC;
    lsLogC    = std::log(c);
    lsLogNorm = std::log(-std::log1p(-c));
  }

  // ---------------------------------------------------------------
  // M-155: Batched precomputation of per-site likelihoods
  //
  // Instead of calling a per-character JC likelihood ~300 times (one per
  // character per candidate k'), batch all characters sharing the same
  // candidate k into a single partition-level pruning call.
  // ---------------------------------------------------------------

  // Beta-Geometric: precompute incremental log-prior for ko = 0..kMaxKprimeCand-1
  // logPrior(u=0) = log(α) - log(α+β)
  // logPrior(u=k) = logPrior(u=k-1) + log(β+k-1) - log(α+β+k)
  std::vector<double> bgLogPrior;
  if (isBetaGeometric) {
    double a = state->kprimeAlpha;
    double b = state->kprimeBeta;
    bgLogPrior.resize(kMaxKprimeCand);
    bgLogPrior[0] = std::log(a) - std::log(a + b);
    for (int ko = 1; ko < kMaxKprimeCand; ++ko) {
      bgLogPrior[ko] = bgLogPrior[ko - 1]
                      + std::log(b + ko - 1)
                      - std::log(a + b + ko);
    }
  }

  // Empirical-geometric: declared here, populated after transParts is built.
  std::vector<double> egLogPriorByK;

  // Build globalCharIdx → transIdx map (position in transIdxGlobal)
  std::vector<int> globalToTransIdx(data->nChar, -1);
  for (int ti = 0; ti < nTrans; ++ti)
    globalToTransIdx[data->transIdxGlobal[ti]] = ti;

  // Identify transformational partitions and determine k range
  // M-172: also track nUniq (unique tip-patterns per partition)
  struct TransPartInfo {
    int partIdx;
    int kObs;
    int nChar;
    int nUniq;  // number of unique tip-state patterns
  };
  std::vector<TransPartInfo> transParts;
  int maxNUniqPart = 0;
  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& p = data->parts[pi];
    if (p.type != 1) continue;
    int kObs_p = data->kObs[p.globalCharIdx[0]];
    int nCh    = p.tipStates.ncol();
    int nUniq  = p.nUniquePatterns;
    transParts.push_back({pi, kObs_p, nCh, nUniq});
    if (nUniq > maxNUniqPart) maxNUniqPart = nUniq;
  }

  // MARGINAL-K-TRUNC-001: the plain geometric prior is truncated at K, so
  // candidates with k' > K carry zero mass and phase 2 discards them (#66).
  std::vector<int> partKoMax(transParts.size(), kMaxKprimeCand);
  if (isGeometric) {
    for (int pi = 0; pi < (int)transParts.size(); ++pi) {
      int nEff = data->kprimeTruncK - transParts[pi].kObs + 1;
      if (nEff < 0) nEff = 0;
      if (nEff < kMaxKprimeCand) partKoMax[pi] = nEff;
    }
  }

  // Ensure Gibbs workspace is large enough for the worst-case stride.
  // M-172: stride is nUniq × k (not nChar × k) — savings proportional to redundancy.
  int maxNode = 2 * data->nTip - 1;
  int gibbsMaxStride = 0;
  for (auto& tp : transParts) {
    int s = tp.nUniq * (tp.kObs + kMaxKprimeCand);
    if (s > gibbsMaxStride) gibbsMaxStride = s;
  }
  if (!state->gibbsWs.fits(maxNode, gibbsMaxStride))
    state->gibbsWs.allocate(maxNode, gibbsMaxStride);

  // Empirical-geometric: precompute log P(k' = m) for all m in the candidate
  // range across all partitions.  The convolution depends only on k', not on
  // kObs, so tabulate once and look up by k'.
  if (isEmpGeom) {
    int maxKAll = 0;
    for (auto& tp_ : transParts) {
      int hi = tp_.kObs + kMaxKprimeCand - 1;
      if (hi > maxKAll) maxKAll = hi;
    }
    egLogPriorByK.assign(maxKAll + 2, R_NegInf);
    double logQ = (data->empTailDecay > 0.0)
                  ? std::log(data->empTailDecay) : R_NegInf;
    int bodyLen = (int)data->empLogBody.size();
    std::vector<double> terms;
    terms.reserve(64);
    for (int m = 2; m <= maxKAll + 1; ++m) {
      terms.clear();
      double mx = R_NegInf;
      for (int j = 2; j <= m; ++j) {
        double logPemp;
        int bodyIdx = j - 2;
        if (bodyIdx < bodyLen) {
          logPemp = data->empLogBody[bodyIdx];
        } else if (data->empTailStartK > 0 && j >= data->empTailStartK &&
                   std::isfinite(data->empLogTailStartP) &&
                   std::isfinite(logQ)) {
          logPemp = data->empLogTailStartP +
                    (j - data->empTailStartK) * logQ;
        } else {
          continue;
        }
        if (!std::isfinite(logPemp)) continue;
        double term = logPemp + logP + (m - j) * log1mP;
        terms.push_back(term);
        if (term > mx) mx = term;
      }
      if (!terms.empty() && std::isfinite(mx)) {
        double s = 0.0;
        for (double t : terms) s += std::exp(t - mx);
        egLogPriorByK[m] = mx + std::log(s);
      }
    }
  }

  bool useHet = data->qHeterogeneity;
  int nBC = data->nBetaCat;
  int coding = data->codingType;
  double hetBins[16];
  int nTip = data->nTip;

  // Per-character sampling state. Outputs charLogW / charNCand / charMaxLogW
  // are written into `out`; charMaxLL / terminated are phase-1 scratch.
  std::vector<double>& charLogW    = out.charLogW;
  std::vector<int>&    charNCand   = out.charNCand;
  std::vector<double>& charMaxLogW = out.charMaxLogW;
  // M-164: track best corrected log-likelihood per character for pre-filter
  std::vector<double> charMaxLL(nTrans, R_NegInf);
  std::vector<bool> terminated(nTrans, false);
  int nActive = nTrans;

  // M-172: per-site output buffer sized for unique patterns (≤ nChar)
  std::vector<double> siteLL(maxNUniqPart);

  // M-172: Per-partition active unique-pattern tracking.
  // patTrans[localPatIdx] lists all trans indices (ti) sharing that pattern.
  // Characters with the same tip-state column have identical likelihoods under
  // any (k, tree, params) for the symmetric JC model — evaluate once, replicate.
  struct PartActive {
    std::vector<int> activePatterns;         // local unique-pattern indices still active
    std::vector<std::vector<int>> patTrans;  // patTrans[localPat] = trans indices
  };
  std::vector<PartActive> partAct(transParts.size());
  for (int pi = 0; pi < (int)transParts.size(); ++pi) {
    const PartInfo& part = data->parts[transParts[pi].partIdx];
    int nCh   = transParts[pi].nChar;
    int nUniq = transParts[pi].nUniq;
    partAct[pi].patTrans.resize(nUniq);
    for (int c = 0; c < nCh; ++c) {
      int ti       = globalToTransIdx[part.globalCharIdx[c]];
      int localPat = part.patternIndex[c];
      partAct[pi].patTrans[localPat].push_back(ti);
    }
    partAct[pi].activePatterns.resize(nUniq);
    for (int j = 0; j < nUniq; ++j) partAct[pi].activePatterns[j] = j;
  }

  // Progressive batched precomputation with early termination (M-155 + M-172)
  //
  // M-172: operates on unique tip-state patterns within each partition.
  // siteLL[ai] is the likelihood for the ai-th active unique pattern; results
  // are scattered to all characters sharing that pattern after each call.
  for (int ko = 0; ko < kMaxKprimeCand && nActive > 0; ++ko) {

    for (int pi = 0; pi < (int)transParts.size(); ++pi) {
      auto& tp = transParts[pi];
      auto& pa = partAct[pi];
      int nAct = (int)pa.activePatterns.size();
      if (nAct == 0) continue;

      if (ko >= partKoMax[pi]) {
        for (int localPat : pa.activePatterns)
          for (int ti : pa.patTrans[localPat])
            if (!terminated[ti]) { terminated[ti] = true; nActive--; }
        pa.activePatterns.clear();
        continue;
      }

      int k = tp.kObs + ko;
      double csp = getCSP(k);
      double logAscCorr = (coding != 0 && csp < 1.0)
        ? -std::log(1.0 - csp) : 0.0;
      if (coding != 0 && csp >= 1.0) logAscCorr = R_NegInf;

      const PartInfo& part = data->parts[tp.partIdx];

      // M-164: prior-ceiling pre-filter for ko ≥ 2.
      // Use representative (first) trans index per pattern — all sharing a
      // pattern have identical weights so terminate together.
      if (ko >= 2) {
        std::vector<int> survivePatterns;
        for (int ai = 0; ai < nAct; ++ai) {
          int localPat = pa.activePatterns[ai];
          int ti0 = pa.patTrans[localPat][0];
          double logPrior_k;
          if (isBetaGeometric) {
            logPrior_k = bgLogPrior[ko];
          } else if (isGeometric) {
            logPrior_k = logP + ko * log1mP;
          } else if (isEmpGeom) {
            int k2 = tp.kObs + ko;
            logPrior_k = (k2 >= 0 && k2 < (int)egLogPriorByK.size())
                         ? egLogPriorByK[k2] : R_NegInf;
          } else {
            int k2 = tp.kObs + ko;
            logPrior_k = k2 * lsLogC
                       - std::log(static_cast<double>(k2)) - lsLogNorm;
          }
          double optimisticW = beta * charMaxLL[ti0] + logPrior_k;
          if (optimisticW < charMaxLogW[ti0] + kKprimeLogCutoff) {
            for (int ti : pa.patTrans[localPat]) {
              if (!terminated[ti]) { terminated[ti] = true; nActive--; }
            }
          } else {
            survivePatterns.push_back(localPat);
          }
        }
        pa.activePatterns = std::move(survivePatterns);
        nAct = (int)pa.activePatterns.size();
        if (nAct == 0) continue;
      }

      // M-172: call pruning on unique-pattern matrix (nAct ≤ nUniq ≤ nChar).
      // When all unique patterns are still active, use uniqueTipStates directly;
      // otherwise build a sub-matrix of the still-active unique patterns.
      //
      // JC-COLLAPSE: when ko >= 2 and not Het, run the lumped-state kernel at
      // kEff = kObs + 1; pruning cost scales as kEff/k (~k/3 for typical
      // kObs=2). Het is deferred — JC lumpability holds only under equal
      // stationary frequencies within the lumped class. The workspace stride
      // gibbsMaxStride is already worst-case (nUniq*(kObs+kMaxKprimeCand)) so the
      // collapsed kernel always fits without reallocation.
      const bool useCollapse = (!useHet) && (ko >= 2);
      int neededStride = nAct * (useCollapse ? (tp.kObs + 1) : k);
      if (nAct == tp.nUniq) {
        // All unique patterns active — use uniqueTipStates directly
        if (useHet) {
          gibbs_compute_het_bins(state->betaScale, k, nBC, hetBins);
          NumericVector rates = acrvRates;
          if (rates.size() == 0) rates = NumericVector(1, 1.0);
          pruning_f81_het_acrv_persite(
            parent, child, edgeLen, part.uniqueTipStates,
            k, 1.0, hetBins, nBC, rates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else if (useCollapse) {
          pruning_jc_acrv_persite_collapsed(
            parent, child, edgeLen, part.uniqueTipStates,
            k, tp.kObs, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else {
          pruning_jc_acrv_persite(
            parent, child, edgeLen, part.uniqueTipStates,
            k, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        }
      } else {
        // Build sub-matrix of active unique patterns
        IntegerMatrix sub(nTip, nAct);
        for (int ai = 0; ai < nAct; ++ai)
          for (int t = 0; t < nTip; ++t)
            sub(t, ai) = part.uniqueTipStates(t, pa.activePatterns[ai]);

        if (useHet) {
          gibbs_compute_het_bins(state->betaScale, k, nBC, hetBins);
          NumericVector rates = acrvRates;
          if (rates.size() == 0) rates = NumericVector(1, 1.0);
          pruning_f81_het_acrv_persite(
            parent, child, edgeLen, sub,
            k, 1.0, hetBins, nBC, rates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else if (useCollapse) {
          pruning_jc_acrv_persite_collapsed(
            parent, child, edgeLen, sub,
            k, tp.kObs, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else {
          pruning_jc_acrv_persite(
            parent, child, edgeLen, sub,
            k, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        }
      }

      // M-172: scatter siteLL[ai] to all characters sharing the active pattern.
      // logAscCorr and relabelling depend only on (k, kObs) — same for all
      // characters sharing a pattern — so corrected LL is identical across the
      // group and termination can be decided from one representative.
      double logPrior_k;
      if (isBetaGeometric) {
        logPrior_k = bgLogPrior[ko];
      } else if (isGeometric) {
        logPrior_k = logP + ko * log1mP;
      } else if (isEmpGeom) {
        logPrior_k = (k >= 0 && k < (int)egLogPriorByK.size())
                     ? egLogPriorByK[k] : R_NegInf;
      } else {
        logPrior_k = k * lsLogC
                   - std::log(static_cast<double>(k)) - lsLogNorm;
      }

      std::vector<int> stillActive;
      for (int ai = 0; ai < nAct; ++ai) {
        int localPat = pa.activePatterns[ai];
        double ll = siteLL[ai];
        if (R_FINITE(ll) && R_FINITE(logAscCorr))
          ll += logAscCorr;
        else if (!R_FINITE(logAscCorr))
          ll = R_NegInf;
        if (data->relabel && R_FINITE(ll))
          ll += mk_prime_relabel_log(k, tp.kObs);

        double w = beta * ll + logPrior_k;

        // Update per-character tracking for every character sharing this pattern
        for (int ti : pa.patTrans[localPat]) {
          if (R_FINITE(ll) && ll > charMaxLL[ti]) charMaxLL[ti] = ll;
          charLogW[ti * kMaxKprimeCand + ko] = w;
          charNCand[ti]++;
          if (w > charMaxLogW[ti]) charMaxLogW[ti] = w;
        }

        // Termination check: representative ti0 holds the same charMaxLogW as all
        // others in the group (identical weights throughout), so one check suffices.
        int ti0 = pa.patTrans[localPat][0];
        if (!R_FINITE(ll) || w < charMaxLogW[ti0] + kKprimeLogCutoff) {
          for (int ti : pa.patTrans[localPat]) {
            if (!terminated[ti]) { terminated[ti] = true; nActive--; }
          }
        } else {
          stillActive.push_back(localPat);
        }
      }

      pa.activePatterns = std::move(stillActive);
    }
  }
}


// ---------------------------------------------------------------------------
// cpp_log_likelihood_marginal (v1: geometric arm only)
//
// Marginal-k evaluator. Replaces per-character k'_i sampling with an
// analytic sum over u_i ∈ {0, .., kMaxKprimeCand-1} weighted by the
// geometric pmf P(u_i | p) = p (1-p)^u.
//
// For neomorphic partitions, falls through to cpp_partition_log_likelihood
// (k = kObs is the fixed correct choice — no marginalisation needed).
//
// Strategy (plan §3):
//   1. Call compute_per_kprime_log_lik(beta = 1.0). The returned
//      charLogW[ti, ko] = LL_ko + log P(ko | p) — i.e. exactly what the
//      marginal logSumExp wants under geometric. So the result for char ti
//      is simply logSumExp_ko(charLogW[ti, ko]).
//   2. Cache strategy: for the charLL cache to accelerate p-only moves,
//      stash charLL[ti, ko] = charLogW[ti, ko] - log P(ko | p_at_eval).
//      On a subsequent p-move (case 30) the marginal evaluator can skip
//      step 1 entirely and recompute the per-character logSumExp against
//      log P(ko | p_new).
//
// State invariant: this is an evaluator, not a sampler. It does not modify
// state->kPrime, state->logLik, state->logPrior, state->partLogLik, or any
// MCMC state variable. It MAY write state->gibbsWs (workspace allocation)
// and state->charLLCache (cache fill — see field comment in McmcState).
// ---------------------------------------------------------------------------

double cpp_log_likelihood_marginal(
    McmcData& data, McmcState& state,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    double rateLoss, double rateLogSd, double rateNeo,
    ClWorkspace* ws, bool fillCharLLCache) {

  // Defensive — v1 enforces these at the R layer in MkPrimeModel(), but
  // duplicate the check here in case a caller bypasses the constructor.
  if (data.kPriorLogseries || data.kPriorBetaGeometric ||
      data.kPriorEmpiricalGeometric) {
    Rcpp::stop("cpp_log_likelihood_marginal: v1 supports the geometric arm "
               "only (kPriorLogseries / kPriorBetaGeometric / "
               "kPriorEmpiricalGeometric must all be false).");
  }
  if (data.qHeterogeneity) {
    Rcpp::stop("cpp_log_likelihood_marginal: Q-matrix heterogeneity is not "
               "supported in v1 marginal-k mode (plan §13).");
  }

  const int nTrans = (int)data.transIdxGlobal.size();
  const int nParts = (int)data.parts.size();
  double totalLL = 0.0;

  // ---- Neomorphic partitions: fall through, unchanged ----------------------
  // k = kObs is exact (state->kPrime[ neo char ] == kObs in this code path
  // — case 25 never touches neo chars; the chain's state->kPrime starts
  // there and stays there).
  for (int pi = 0; pi < nParts; ++pi) {
    if (data.parts[pi].type == 0) {
      totalLL += cpp_partition_log_likelihood(
        data, pi, parent, child, edgeLen,
        state.kPrime, rateLoss, rateLogSd, rateNeo,
        state.betaScale, ws);
    }
  }

  if (nTrans == 0) return totalLL;

  // ---- Transformational partitions: marginal sum --------------------------
  // Geometric prior weights P(u | p): precompute log P(u | p) for u in
  // 0..kMaxKprimeCand-1. Same series the helper uses internally.
  const double p     = state.p;
  const double logP  = std::log(p);
  const double log1mP = std::log1p(-p);
  std::vector<double> logPriorByU(kMaxKprimeCand);
  for (int ko = 0; ko < kMaxKprimeCand; ++ko) {
    logPriorByU[ko] = logP + ko * log1mP;
  }

  // Model A (unconditional) vs Model B (conditional) differ only by a
  // per-character constant factor (1-p)^(kObs_i - 2). Under Model B the
  // marginal weight for state count k is p (1-p)^(k - kObs_i); under Model A
  // it is p (1-p)^(k - 2). Since k = kObs_i + ko, the Model A weight equals
  // the Model B weight times (1-p)^(kObs_i - 2), and that factor is constant
  // across all candidates ko of a given character, so it factors straight out
  // of the per-character logSumExp. We therefore keep logPriorByU (and the
  // cached raw LL) as Model B and add (kObs_i - 2)·log(1-p) to each
  // character's marginal contribution when unconditionalPrior is true.
  //
  // MARGINAL-K-TRUNC-001 (geometric truncation normaliser; proof
  // dev/red-team/proofs/marginal-k-truncation-normaliser.md). The geometric
  // prior is truncated at the declared cap K = data.kprimeTruncK (matching the
  // forward). The marginal is therefore BOTH (a) summed only over candidates
  // with k = kObs_i + c <= K, AND (b) renormalised by the truncated-tail mass
  //   Model A: logZA      = log(1 - (1-p)^(K-1))            [shared]
  //   Model B: logZB(kObs) = log(1 - (1-p)^(K - kObs_i + 1)) [per kObs]
  // subtracted once per character. Both are required: -logZ alone over-corrects
  // (biases p low at small p); the cap alone leaves the p-high bias. (Proof §5.3
  // as corrected 2026-05-30; R spec-check cap-z-spec-check.R.)
  const bool uncond = data.unconditionalPrior;
  const int    K     = data.kprimeTruncK;
  const double logZA = std::log1p(-std::exp((K - 1) * log1mP));

  // Cache fast-path: if the charLL cache is valid, skip the helper call
  // and recompute the per-char logSumExp against the current p-weights.
  bool useCache = fillCharLLCache && data.marginalK && state.charLLCacheReady &&
                  (int)state.charLLNCand.size() == nTrans &&
                  (int)state.charLLCache.size() ==
                    nTrans * kMaxKprimeCand;

  if (!useCache) {
    // Compute per-(char, ko) Gibbs weights via the PR-A helper at β = 1.
    // For the geometric arm these weights are LL + log P(u|p) — exactly
    // what the marginal logSumExp wants.
    NumericVector acrvRates = (rateLogSd > 0.0)
      ? cpp_acrv_rates(rateLogSd, data.nCat, data.acrvZ)
      : NumericVector(1, 1.0);

    // Save state's rateLogSd/etc temporarily? No — helper reads p from
    // state but does not directly use rateLoss/rateLogSd/rateNeo (those
    // enter via the per-partition pruning calls inside the helper, which
    // we feed `edgeLen` to externally). The helper signature passes
    // edgeLen directly, so we honour the caller's rateLoss/rateLogSd/
    // rateNeo by pre-scaling edgeLen if the caller did — but the existing
    // case-25 path uses state's own values. The MCMC chain always
    // evaluates at state's current scalars, so we mirror that by using
    // them here. Topology (parent/child/edgeLen) comes from the caller; the
    // scalar params (p/rateLoss/rateLogSd/rateNeo/betaScale) are read from
    // state inside the helper. (MARGINAL-K-FREEZE-003: pass the caller's
    // parent/child so a proposed topology is evaluated, not state's old tree.)
    KprimeCharWeights kw;
    compute_per_kprime_log_lik(&data, &state, /*beta=*/1.0,
                               parent, child, edgeLen, acrvRates, kw);

    // Allocate cache lazily on first marginal eval (skip entirely for scratch
    // evals, which never write the cache).
    if (fillCharLLCache &&
        (int)state.charLLCache.size() != nTrans * kMaxKprimeCand) {
      state.charLLCache.assign(
        static_cast<size_t>(nTrans) * kMaxKprimeCand, R_NegInf);
      state.charLLNCand.assign(nTrans, 0);
    }

    // Fill cache: charLL = charLogW - log P(u | p_at_eval).
    // Per-char logSumExp on charLogW directly == marginal LL for char ti.
    for (int ti = 0; ti < nTrans; ++ti) {
      const int kObs_ti = data.kObs[data.transIdxGlobal[ti]];
      const int nCand = kw.charNCand[ti];
      // MARGINAL-K-TRUNC-001 cap: candidate c has k = kObs_ti + c; keep only
      // k <= K. nEff <= 0 means kObs_ti > K (empty support) -> -Inf char.
      const int nEff = std::min(nCand, K - kObs_ti + 1);
      if (fillCharLLCache) state.charLLNCand[ti] = (nEff > 0) ? nEff : 0;
      double mx = R_NegInf;
      for (int c = 0; c < nCand; ++c) {
        double w = kw.charLogW[ti * kMaxKprimeCand + c];
        // Stash raw LL = w - logPriorByU[c] in cache (all candidates, so a
        // later K change or audit can re-derive; only c < nEff are summed).
        // Scratch evals (fillCharLLCache=false) skip the write but still need
        // mx, so the loop body otherwise runs unchanged.
        if (fillCharLLCache) {
          double rawLL = R_FINITE(w) ? (w - logPriorByU[c]) : R_NegInf;
          state.charLLCache[ti * kMaxKprimeCand + c] = rawLL;
        }
        if (c < nEff && R_FINITE(w) && w > mx) mx = w;   // cap: only k <= K
      }
      if (!R_FINITE(mx)) {
        // No finite candidate (incl. kObs_ti > K) — character unsupportable;
        // total LL is -Inf and downstream MH ratio will reject.
        totalLL = R_NegInf;
        // Continue filling cache for the remaining chars so a subsequent
        // call (e.g. on rejection rollback) finds a coherent cache.
        continue;
      }
      double s = 0.0;
      for (int c = 0; c < nEff; ++c) {
        double w = kw.charLogW[ti * kMaxKprimeCand + c];
        if (R_FINITE(w)) s += std::exp(w - mx);
      }
      double charLL = mx + std::log(s);
      if (uncond) {
        charLL += (kObs_ti - 2) * log1mP - logZA;                       // Model A
      } else {
        charLL -= std::log1p(-std::exp((K - kObs_ti + 1) * log1mP));    // Model B
      }
      if (R_FINITE(totalLL)) totalLL += charLL;
    }
    // Scratch evals leave the cache cold (charLLCacheReady stays whatever the
    // entry invalidation set — false for every non-p move), so no rejected or
    // self-accepting multi-config move can poison a subsequent mh_logit_p.
    if (fillCharLLCache) state.charLLCacheReady = true;
    return totalLL;
  }

  // Cache fast-path: reuse stored raw charLL, only re-do the logSumExp
  // against the new logPriorByU weights.
  for (int ti = 0; ti < nTrans; ++ti) {
    const int nCand = state.charLLNCand[ti];
    double mx = R_NegInf;
    // Compute weights into a small local buffer to avoid double-exp cost.
    double wBuf[kMaxKprimeCand];
    for (int c = 0; c < nCand; ++c) {
      double rawLL = state.charLLCache[ti * kMaxKprimeCand + c];
      double w = R_FINITE(rawLL)
        ? (rawLL + logPriorByU[c])
        : R_NegInf;
      wBuf[c] = w;
      if (R_FINITE(w) && w > mx) mx = w;
    }
    if (!R_FINITE(mx)) { totalLL = R_NegInf; continue; }
    double s = 0.0;
    for (int c = 0; c < nCand; ++c)
      if (R_FINITE(wBuf[c])) s += std::exp(wBuf[c] - mx);
    // MARGINAL-K-TRUNC-001: the k<=K cap is INHERITED here via charLLNCand[ti]
    // (= nEff, set capped during the fill above), so the loops already stop at
    // k=K; only the -logZ renormaliser is applied here.
    double charLL = mx + std::log(s);
    const int kObs_ti = data.kObs[data.transIdxGlobal[ti]];
    if (uncond) {
      charLL += (kObs_ti - 2) * log1mP - logZA;                       // Model A
    } else {
      charLL -= std::log1p(-std::exp((K - kObs_ti + 1) * log1mP));    // Model B
    }
    if (R_FINITE(totalLL)) totalLL += charLL;
  }
  return totalLL;
}


// ---------------------------------------------------------------------------
// kprime_sweep_candidates — how many k' candidates phase 1 evaluates per
// transformational character at the current state. Diagnostic only: the sweep
// itself is unaffected. Used to assert that the enumerated range tracks
// kprimeTruncK (issue #66) without a wall-clock measurement.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
IntegerVector kprime_sweep_candidates(SEXP dataPtr, SEXP statePtr,
                                      double beta = 1.0) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();

  int nEdge = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  NumericVector acrvRates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);

  KprimeCharWeights kw;
  compute_per_kprime_log_lik(data, state, beta,
                             state->parent, state->child, edgeLen,
                             acrvRates, kw);
  return Rcpp::wrap(kw.charNCand);
}


// ---------------------------------------------------------------------------
// Gibbs kPrime sweep (moveType 25)
//
// Samples each k'_i from its full conditional in a single random-order scan
// of all transformational characters. Always accepts (Gibbs update).
//
// Phase 1 (per-(char, k') weight precomputation) lives in
// compute_per_kprime_log_lik above; phase 2 (categorical sampling + cache
// rebuild) lives here. Same RNG draws as the pre-refactor monolithic
// implementation — the permutation Fisher-Yates and per-character categorical
// inverse-CDF samples consume R::unif_rand in the same order.
// ---------------------------------------------------------------------------

static bool gibbs_kprime_sweep_impl(McmcData* data, McmcState* state,
                                     double beta) {
  int nTrans = (int)data->transIdxGlobal.size();
  if (nTrans == 0) return false;

  // Pre-compute absolute edge lengths
  int nEdge = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  // Pre-compute ACRV rates
  NumericVector acrvRates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);

  // Phase 1: batched per-(char, k') log-weight precomputation.
  // Gibbs kPrime sweep evaluates at state's CURRENT tree, so pass state's
  // own parent/child (MARGINAL-K-FREEZE-003: helper no longer reads state).
  KprimeCharWeights kw;
  compute_per_kprime_log_lik(data, state, beta,
                             state->parent, state->child, edgeLen, acrvRates, kw);
  const std::vector<double>& charLogW    = kw.charLogW;
  const std::vector<int>&    charNCand   = kw.charNCand;
  const std::vector<double>& charMaxLogW = kw.charMaxLogW;

  // ---------------------------------------------------------------
  // Phase 2 — sampling phase: for each character, sample k' from the
  // precomputed log-weights in charLogW[].
  // ---------------------------------------------------------------

  // MARGINAL-K-TRUNC-001: for the (plain) geometric arm the prior is truncated
  // at K, so candidates with k' = kObs_i + c > K carry zero mass and must be
  // dropped from the Gibbs draw. This move always-accepts, so an out-of-support
  // draw would NOT be MH-rejected downstream; the cap must happen here. Mirrors
  // the marginal evaluator's nEff cap (compute over [0, nEff), nEff = min(nCand,
  // K - kObs_i + 1)). Other arms are not K-truncated and keep the full range.
  const bool isPlainGeom = !data->kPriorLogseries &&
                           !data->kPriorBetaGeometric &&
                           !data->kPriorEmpiricalGeometric;
  const int  truncK = data->kprimeTruncK;

  // Random permutation of transformational character indices
  std::vector<int> perm(nTrans);
  for (int i = 0; i < nTrans; ++i) perm[i] = i;
  for (int i = nTrans - 1; i > 0; --i) {
    int j = static_cast<int>(R::unif_rand() * (i + 1));
    if (j > i) j = i;
    std::swap(perm[i], perm[j]);
  }

  for (int si = 0; si < nTrans; ++si) {
    int ti = perm[si];
    int gi = data->transIdxGlobal[ti];
    int kObs_i = data->kObs[gi];
    int nCand = charNCand[ti];
    if (nCand == 0) continue;

    double maxW = charMaxLogW[ti];
    const double* logW = &charLogW[ti * kMaxKprimeCand];

    // MARGINAL-K-TRUNC-001 cap (geometric only). When truncation actually bites
    // (nEff < nCand), recompute maxW over the retained range: the cached
    // charMaxLogW is the max over the FULL pre-cap range and may sit above K,
    // which would skew the logSumExp toward an empty tail. Untruncated and
    // not-geometric paths keep the original maxW for bit-reproducibility.
    if (isPlainGeom) {
      int nEff = truncK - kObs_i + 1;
      if (nEff <= 0) continue;            // empty support: kObs_i > K
      if (nEff < nCand) {
        nCand = nEff;
        maxW = R_NegInf;
        for (int c = 0; c < nCand; ++c)
          if (R_FINITE(logW[c]) && logW[c] > maxW) maxW = logW[c];
        if (!R_FINITE(maxW)) continue;
      }
    }

    // Sample from categorical (log-sum-exp)
    double sumExp = 0.0;
    for (int c = 0; c < nCand; ++c)
      sumExp += std::exp(logW[c] - maxW);

    double u = R::unif_rand() * sumExp;
    double cum = 0.0;
    int chosen = nCand - 1;
    for (int c = 0; c < nCand; ++c) {
      cum += std::exp(logW[c] - maxW);
      if (cum >= u) { chosen = c; break; }
    }

    state->kPrime[gi] = kObs_i + chosen;
  }

  // Rebuild logLik, logPrior, and partition cache after sweep
  ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
  int nParts = (int)data->parts.size();
  std::vector<double> newPLC(nParts);
  double newLL = 0.0;
  for (int pi = 0; pi < nParts; ++pi) {
    newPLC[pi] = cpp_partition_log_likelihood(
      *data, pi, state->parent, state->child, edgeLen,
      state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
      state->betaScale, wsPtr);
    newLL += newPLC[pi];
  }
  state->logLik = newLL;
  state->partLogLik = std::move(newPLC);

  state->logPrior = compute_log_prior(*data, *state);

  state->nodeCL.invalidate_structure();  // M-161: kPrime changed, unit structure may differ

  return true;  // Gibbs: always accept
}


// ---------------------------------------------------------------------------
// Block kPrime shift (moveType 26)
//
// Proposes shifting ALL transformational characters by the same integer delta.
// Standard MH acceptance with symmetric proposal.
// ---------------------------------------------------------------------------

static bool block_kprime_shift_impl(McmcData* data, McmcState* state,
                                     int intWalkWindow, double beta) {
  int nTrans = (int)data->transIdxGlobal.size();
  if (nTrans == 0) return false;

  // Propose delta ~ Uniform({-W, ..., W})
  int range = 2 * intWalkWindow + 1;
  int delta = static_cast<int>(R::unif_rand() * range) - intWalkWindow;
  if (delta == 0) return false;

  // Feasibility: all k'_i + delta >= kObs_i
  for (int i = 0; i < nTrans; ++i) {
    int gi = data->transIdxGlobal[i];
    if (state->kPrime[gi] + delta < data->kObs[gi]) return false;
  }

  // Save old values and apply shift
  std::vector<int> oldKPrime(nTrans);
  for (int i = 0; i < nTrans; ++i) {
    int gi = data->transIdxGlobal[i];
    oldKPrime[i] = state->kPrime[gi];
    state->kPrime[gi] += delta;
  }

  // Recompute log-prior
  double newLogPrior = compute_log_prior(*data, *state);

  if (!R_FINITE(newLogPrior)) {
    for (int i = 0; i < nTrans; ++i)
      state->kPrime[data->transIdxGlobal[i]] = oldKPrime[i];
    return false;
  }

  // Recompute log-likelihood (only transformational partitions affected)
  int nEdge = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
  bool hasPLC = !state->partLogLik.empty();
  int nParts = (int)data->parts.size();
  std::vector<double> newPC;
  double newLogLik;

  if (hasPLC) {
    newPC = state->partLogLik;
    newLogLik = state->logLik;
    for (int pi = 0; pi < nParts; ++pi) {
      if (data->parts[pi].type == 1) {  // transformational only
        double v = cpp_partition_log_likelihood(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
          state->betaScale, wsPtr);
        newLogLik += (v - newPC[pi]);
        newPC[pi] = v;
      }
    }
  } else {
    newPC.resize(nParts);
    newLogLik = 0.0;
    for (int pi = 0; pi < nParts; ++pi) {
      newPC[pi] = cpp_partition_log_likelihood(
        *data, pi, state->parent, state->child, edgeLen,
        state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
        state->betaScale, wsPtr);
      newLogLik += newPC[pi];
    }
  }

  // MH acceptance (symmetric proposal, logHastings = 0)
  double logAlpha = beta * (newLogLik - state->logLik) +
                    (newLogPrior - state->logPrior);

  if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
    state->logLik = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik = std::move(newPC);
    state->nodeCL.invalidate_structure();  // M-161: kPrime changed
    return true;
  }

  // Rollback
  for (int i = 0; i < nTrans; ++i)
    state->kPrime[data->transIdxGlobal[i]] = oldKPrime[i];
  return false;
}


// ---------------------------------------------------------------------------
// do_move_impl: internal propose/evaluate/accept (raw pointers, no SEXP).
// do_move_cpp:  Rcpp-exported SEXP wrapper — calls do_move_impl.
//
// moveType: 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
//           4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p (legacy),
//           9=gibbs_p, 10=gibbs_spr, 11=gibbs_subtree_swap,
//           12=weighted_br_scale, 13=weighted_spr, 14=weighted_subtree_swap,
//           15=block_gibbs_branch, 16=beta_scale, 17=tbr,
//           19=slice_scalar, 20=pspr,
//           21=joint_tl_rls, 22=joint_tl_rl,
//           23=dirichlet_branch, 24=local_dirichlet,
//           25=gibbs_kprime_sweep, 26=block_kprime_shift,
//           27=scale_kprime_alpha, 28=scale_kprime_beta,
//           29=slice_kprime_hyper,
//           30=mh_logit_p (logit-scale MH on p for empirical_geometric prior)
//           31=scale_class_rate_log_sd (per-class ACRV shape; charIdx=1-based classIdx;
//              acts on z_c when state->useHyperpriorOnSigma, else on σ_c directly)
//           32=dirichlet_simplex_class_w (Dirichlet simplex on class_w)
//           33=joint_tl_rn (tree_length × rate_neo 2D Bactrian; partition-rate
//              ridge from Issue 1 fix — see dev/notes/2026-05-27-rate-neo-ridge-and-joint-moves.md)
//           34=scale_hyper_tau (Bactrian scale on the population scale τ of
//              the half-normal hyperprior on σ_c = τ·z_c; only fired when
//              state->useHyperpriorOnSigma)
//           35=gibbs_p_marginal (data-augmentation Metropolis-within-Gibbs
//              p-update for the marginal_k geometric arm: imputes latent u_i
//              from the cached per-(char,k') weights, proposes p* from the
//              untruncated conjugate Beta, accepts with the truncation-
//              normaliser ratio. Like case 30, a p-only move that PRESERVES
//              the charLL cache. Derivation + math-prover verification:
//              dev/red-team/proofs/marginal-k-gibbs-p.md)
//
// M-065: NNI/SPR now call _impl versions directly with parent/child vectors.
// Likelihood calls use vectors directly (no IntegerMatrix construction).
// ---------------------------------------------------------------------------

static bool do_move_impl(McmcData* data, McmcState* state,
                         int moveType, int charIdx,
                         double scaleTuning, double betaSimplexTuning,
                         int intWalkWindow, double beta,
                         double jointRho = 0.0) {

  // Snapshot scalar state for rollback
  double oldTL   = state->treeLength;
  double oldRL   = state->rateLoss;
  double oldRLSD = state->rateLogSd;
  double oldRN   = state->rateNeo;
  double oldP    = state->p;
  double oldBS   = state->betaScale;  // M-052
  double oldKpA  = state->kprimeAlpha;
  double oldKpB  = state->kprimeBeta;

  // Marginal-k: the charLL cache holds raw per-(char, k') LLs derived
  // from (parent, child, edgeLen, rateLoss, rateLogSd, rateNeo, betaScale)
  // but NOT from p. Only case 30 (mh_logit_p) leaves it valid
  // post-proposal. Defensively invalidate at the top for every other
  // move; the marginal evaluator will rebuild on the next call. Cost on
  // infeasible early-return moves: at most one wasted rebuild on the
  // next call, which is acceptable.
  //
  // FU-3 (Option A Tier 2): the per-(node, ko) Felsenstein CL cache is
  // governed by the same invariance proof as charLLCache — see
  // dev/red-team/proofs/marginal-k-geometric.md §5. A p-move preserves
  // both tiers; any other move invalidates both.
  if (data->marginalK && moveType != 30 && moveType != 35) {
    state->charLLCacheReady = false;
    state->invalidate_per_kp_cl_all();
  }

  double logHastings  = 0.0;
  bool topologyChanged = false;
  int oldKPrimeVal = 0;     // single-element rollback for case 7 (kPrime)
  int kPrimeCharIdx = -1;   // which character was changed

  // Rollback storage for per-class moves (cases 31, 32, 33)
  int    classIdx31   = -1;   // 0-based class index for case 31 rollback
  double oldClassRLS  = 0.0;  // old classRateLogSd[classIdx31]
  double oldClassZ    = 0.0;  // old classZ[classIdx31] (hyperprior path)
  NumericVector classWSnapshot;  // full classW snapshot for case 32 rollback
  // Case 33 (scale_hyper_tau) rollback: all σ_c change together with τ.
  bool             case34Active   = false;
  double           oldTau34       = 0.0;
  NumericVector    oldClassRLS34;  // snapshot of classRateLogSd (all entries)
  double           oldRateLogSd34 = 0.0;

  // Lockstep invariant: classRateLogSd[0] == rateLogSd when partitioned.
  // Cases 2 / 21 mirror writes to classRateLogSd[0]; this snapshot lets the
  // rollback paths below restore it unconditionally without coupling to the
  // case-31 rollback (which only fires when classIdx31 >= 0).
  double oldClassRLS0 =
    (state->usePartitioned && state->classRateLogSd.size() > 0)
    ? state->classRateLogSd[0] : 0.0;

  // O(1) relBr rollback for cases 4 and 12 (save 2 modified elements)
  int bsIdx1 = -1, bsIdx2 = -1;
  double bsOldVal1 = 0.0, bsOldVal2 = 0.0;
  // M-158: nodeEdgeLen rollback for beta_simplex (2 values)
  double bsOldNodeEdgeLen1 = 0.0, bsOldNodeEdgeLen2 = 0.0;

  // OPP-6b: in-place NNI rollback (2 parent values)
  bool nniInPlace = false;
  int nniCRow = -1, nniWRow = -1;
  int nniSavedP_cRow = 0, nniSavedP_wRow = 0;
  // M-121: NNI node identities for partial CL
  int nniVNode = 0, nniUNode = 0, nniCNode = 0, nniWNode = 0;

  // OPP-6: proposed topology held separately for SPR/TBR; state->parent/child
  // not overwritten until acceptance → no pre-proposal clone, no rollback copy.
  IntegerVector proposedParent, proposedChild;
  NumericVector proposedRelBr;
  // M-158: SPR partial CL metadata
  SprMeta sprMeta;
  sprMeta.valid = false;
  bool sprPartialCL = false;  // true if SPR used partial CL path

  // M-121/M-158: pre-proposal cache population for partial CL.
  // Must happen BEFORE the proposal modifies state in-place.
  // Populate for NNI (5), beta_simplex (4), Dirichlet (23, 24), SPR (6).
  //
  // Skipped under marginal-k mode: the partial-CL paths assume sampled-k
  // (they call cpp_partition_log_likelihood which sums over a single
  // k = state->kPrime). The marginal evaluator dispatches through
  // compute_full_loglik_at which always recomputes the full marginal sum
  // from scratch (or from the charLL cache for p-only moves).
  if ((moveType == 5 || moveType == 4 || moveType == 23 || moveType == 24
       || moveType == 6) &&
      !data->qHeterogeneity && !data->marginalK && !state->nodeCL.ready()) {
    state->diagCachePopCount++;
    int nEdge = state->relBrLengths.size();
    NumericVector absLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      absLen[i] = state->treeLength * state->relBrLengths[i];
    populate_cache(state->nodeCL, *data,
                        state->parent, state->child, absLen,
                        state->kPrime, state->rateLoss,
                        state->rateLogSd, state->rateNeo);
  }

  switch (moveType) {
    case 0: { // scale tree_length (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->treeLength = oldTL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 1: { // scale rate_loss (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->rateLoss = oldRL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 2: { // scale rate_log_sd (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->rateLogSd = oldRLSD * mult;
      if (state->usePartitioned && state->classRateLogSd.size() > 0)
        state->classRateLogSd[0] = state->rateLogSd;  // lockstep
      logHastings = std::log(mult);
      break;
    }
    case 3: { // scale rate_neo (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->rateNeo = oldRN * mult;
      logHastings = std::log(mult);
      break;
    }
    case 4: { // beta_simplex — O(1) rollback: save 2 modified elements
      int n_br = state->relBrLengths.size();
      int idx = static_cast<int>(R::unif_rand() * n_br);
      if (idx >= n_br) idx = n_br - 1;
      bsIdx1 = idx;
      if (!beta_simplex_impl(state->relBrLengths, idx,
                             betaSimplexTuning, logHastings,
                             bsIdx2, bsOldVal1, bsOldVal2))
        return false;
      break;
    }
    case 5: { // NNI — OPP-6b: in-place when safe, full reorder otherwise.
      //
      // NNI swaps 2 parent assignments: one child of v moves to u, one
      // child of u (sibling w) moves to v.  The in-place modification
      // preserves valid preorder ONLY when wRow > edgeRow (v is already
      // introduced before wRow in the edge list).  When wRow <= edgeRow,
      // the parent assignment at wRow references v which hasn't appeared
      // as a child yet — breaking the preorder invariant that the
      // reverse-iteration Felsenstein pruning depends on.
      const int nEdge = state->parent.size();
      const int nTip  = data->nTip;

      // Find internal edges (both endpoints internal)
      std::vector<int> intRows;
      intRows.reserve(nEdge / 2);
      for (int i = 0; i < nEdge; ++i)
        if (state->parent[i] > nTip && state->child[i] > nTip)
          intRows.push_back(i);
      if (intRows.empty()) return false;

      int pick = (int)(R::unif_rand() * (double)intRows.size());
      if (pick >= (int)intRows.size()) pick = intRows.size() - 1;
      const int edgeRow = intRows[pick];
      const int u = state->parent[edgeRow];
      const int v = state->child[edgeRow];

      // Find v's children and u's other children (not v)
      std::vector<int> vCh, uSib;
      for (int i = 0; i < nEdge; ++i) {
        if (state->parent[i] == v)
          vCh.push_back(i);
        else if (state->parent[i] == u && state->child[i] != v)
          uSib.push_back(i);
      }
      if (vCh.empty() || uSib.empty()) return false;

      int pV = (int)(R::unif_rand() * (double)vCh.size());
      if (pV >= (int)vCh.size()) pV = vCh.size() - 1;
      int pU = (int)(R::unif_rand() * (double)uSib.size());
      if (pU >= (int)uSib.size()) pU = uSib.size() - 1;
      const int cRow = vCh[pV];
      const int wRow = uSib[pU];

      if (wRow > edgeRow) {
        // Safe: v is introduced at edgeRow, which is before wRow.
        // In-place swap preserves valid preorder.
        nniCRow = cRow; nniWRow = wRow;
        nniSavedP_cRow = state->parent[cRow];
        nniSavedP_wRow = state->parent[wRow];
        // M-121: save node identities for partial CL dirty path
        nniVNode = v; nniUNode = u;
        nniCNode = state->child[cRow];
        nniWNode = state->child[wRow];
        state->parent[cRow] = u;
        state->parent[wRow] = v;
        nniInPlace = true;
      } else {
        // Unsafe: wRow <= edgeRow -- v not yet introduced at wRow.
        // Apply the same NNI swap but canonicalise via reorder.
        IntegerVector newPar = clone(state->parent);
        newPar[cRow] = u;
        newPar[wRow] = v;
        NumericVector absLen(nEdge);
        for (int i = 0; i < nEdge; ++i)
          absLen[i] = state->treeLength * state->relBrLengths[i];
        auto po = TreeTools::preorder_weighted_impl(
          newPar, state->child, absLen);
        IntegerMatrix oe = po.first;
        NumericVector oa = po.second;
        proposedParent = IntegerVector(nEdge);
        proposedChild  = IntegerVector(nEdge);
        for (int i = 0; i < nEdge; ++i) {
          proposedParent[i] = oe(i, 0);
          proposedChild[i]  = oe(i, 1);
        }
        proposedRelBr = oa / state->treeLength;
        topologyChanged = true;
      }
      logHastings = 0.0;
      break;
    }
    case 6: { // SPR — M-158: partial CL when cache valid, else OPP-6 full eval
      if (state->nodeCL.ready() && !data->qHeterogeneity) {
        // M-158: TreeNav-based SPR with partial CL evaluation
        sprMeta = propose_spr_treenav(state->nodeCL.topo);
        if (!sprMeta.valid || !R_FINITE(sprMeta.logHastings)) return false;
        logHastings = sprMeta.logHastings;
        sprPartialCL = true;
        // TreeNav is updated + partial eval happens in the evaluation section
      } else {
        // Fallback: full eval path (original OPP-6 pattern)
        List prop = spr_proposal_impl(state->parent, state->child,
                                      data->nTip, state->treeLength,
                                      state->relBrLengths);
        logHastings = as<double>(prop["logHastings"]);
        if (!R_FINITE(logHastings)) return false;
        proposedParent  = as<IntegerVector>(prop["parent"]);
        proposedChild   = as<IntegerVector>(prop["child"]);
        proposedRelBr   = as<NumericVector>(prop["rel_br_lengths"]);
        topologyChanged = true;
      }
      break;
    }
    case 17: { // TBR — M-053: OPP-6 pattern (defer topology commit)
      List prop = tbr_proposal_impl(state->parent, state->child,
                                    data->nTip, state->treeLength,
                                    state->relBrLengths);
      logHastings = as<double>(prop["logHastings"]);
      if (!R_FINITE(logHastings)) return false;
      proposedParent  = as<IntegerVector>(prop["parent"]);
      proposedChild   = as<IntegerVector>(prop["child"]);
      proposedRelBr   = as<NumericVector>(prop["rel_br_lengths"]);
      topologyChanged = true;
      break;
    }
    case 7: { // int_walk kPrime — O(1) rollback (save single element, not full clone)
      kPrimeCharIdx = charIdx;
      oldKPrimeVal  = state->kPrime[charIdx];
      int oldK      = oldKPrimeVal;
      int lowerK = data->kObs[charIdx];
      int range  = 2 * intWalkWindow + 1;
      int delta  = static_cast<int>(R::unif_rand() * range) - intWalkWindow;
      int newK   = oldK + delta;
      if (newK < lowerK) return false;
      state->kPrime[charIdx] = newK;
      logHastings = 0.0;
      break;
    }
    case 8: { // scale p (legacy MH — kept for backward compat, not used by default)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->p = oldP * mult;
      logHastings = std::log(mult);
      break;
    }
    case 9: { // gibbs_p — conjugate Beta draw, acceptance = 1
      // Full conditional: p | k' ~ Beta(a + nTrans, b + sum(k'_i - kObs_i))
      // p does not enter the likelihood, only the prior on k' and p itself.
      // Empirical-geometric breaks conjugacy: p enters the prior via the
      // convolution, so the full conditional is no longer Beta.  Reject
      // any misrouted call rather than silently producing wrong samples.
      if (data->kPriorEmpiricalGeometric) return false;
      // MARGINAL-K-TRUNC-001: the plain geometric prior is now truncated at K,
      // whose normaliser Z(p) = 1 - (1-p)^(K - kObs_i + 1) (Model B) or
      // 1 - (1-p)^(K-1) (Model A) is p-dependent, so the p full-conditional is
      // likewise non-Beta. p is sampled via case 30 (mh_logit_p) instead; reject
      // any misrouted gibbs_p call. (logseries / beta_geometric carry no scalar
      // p and never schedule gibbs_p, so the plain geometric is the only arm
      // that can reach here past the empirical guard above.)
      {
        const bool isPlainGeom = !data->kPriorLogseries &&
                                 !data->kPriorBetaGeometric &&
                                 !data->kPriorEmpiricalGeometric;
        if (isPlainGeom) return false;
      }
      int nTrans = (int)data->transIdxGlobal.size();
      if (nTrans == 0) return false;  // no transformational chars: skip
      double sumU = 0.0;
      for (int i = 0; i < nTrans; ++i) {
        int gi = data->transIdxGlobal[i];
        sumU += static_cast<double>(state->kPrime[gi] - data->kObs[gi]);
      }
      double shape1 = data->kprimeHyperA + nTrans;
      double shape2 = data->kprimeHyperB + sumU;
      state->p = R::rbeta(shape1, shape2);
      // Recompute prior (p changed; likelihood unchanged)
      state->logPrior = compute_log_prior(*data, *state);
      state->logLik = state->logLik;  // unchanged
      return true;  // Gibbs: always accept
    }
    case 10: { // gibbs_spr — M-085
      return gibbs_spr_impl(data, state, beta);
    }
    case 11: { // gibbs_subtree_swap — M-086
      return gibbs_subtree_swap_impl(data, state, beta);
    }
    case 12: { // weighted_branch_scale — M-087, O(1) rollback
      if (!weighted_branch_scale_impl(data, state, beta, logHastings,
                                       bsIdx1, bsIdx2, bsOldVal1, bsOldVal2))
        return false;
      break;
    }
    case 13: { // weighted_spr — M-088
      return weighted_spr_impl(data, state, beta);
    }
    case 14: { // weighted_subtree_swap — M-089
      return weighted_subtree_swap_impl(data, state, beta);
    }
    case 15: { // block_gibbs_branch — M-054 reframed
      return block_gibbs_branch_sweep_impl(data, state, beta);
    }
    case 16: { // M-052: scale beta_scale (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->betaScale = oldBS * mult;
      logHastings = std::log(mult);
      break;
    }
    case 20: { // pSPR — M-119: parsimony-guided SPR
      List prop = pspr_proposal_impl(state->parent, state->child,
                                     data->nTip, state->treeLength,
                                     state->relBrLengths, data);
      logHastings = as<double>(prop["logHastings"]);
      if (!R_FINITE(logHastings)) return false;
      proposedParent  = as<IntegerVector>(prop["parent"]);
      proposedChild   = as<IntegerVector>(prop["child"]);
      proposedRelBr   = as<NumericVector>(prop["rel_br_lengths"]);
      topologyChanged = true;
      break;
    }
    case 21: { // M-120: joint_tl_rls (tree_length × rate_log_sd)
      double z1, z2;
      bactrian_2d_perturbation(jointRho, z1, z2);
      double mult1 = std::exp(scaleTuning * z1);
      double mult2 = std::exp(scaleTuning * z2);
      state->treeLength = oldTL * mult1;
      state->rateLogSd  = oldRLSD * mult2;
      if (state->usePartitioned && state->classRateLogSd.size() > 0)
        state->classRateLogSd[0] = state->rateLogSd;  // lockstep
      logHastings = std::log(mult1) + std::log(mult2);
      break;
    }
    case 22: { // M-120: joint_tl_rl (tree_length × rate_loss)
      double z1, z2;
      bactrian_2d_perturbation(jointRho, z1, z2);
      double mult1 = std::exp(scaleTuning * z1);
      double mult2 = std::exp(scaleTuning * z2);
      state->treeLength = oldTL * mult1;
      state->rateLoss   = oldRL * mult2;
      logHastings = std::log(mult1) + std::log(mult2);
      break;
    }
    case 23: { // M-125: dirichlet_branch (block Dirichlet on relBrLengths)
      // nCats passed via intWalkWindow for this move type (repurposed)
      int nCats = intWalkWindow;
      if (nCats < 2) nCats = 2;
      if (!dirichlet_simplex_impl(state->relBrLengths, nCats,
                                  scaleTuning, logHastings,
                                  state->brSnapshot,
                                  state->dirEdges)) {
        return false;
      }
      break;
    }
    case 24: { // M-127: local_dirichlet (neighborhood Dirichlet on relBrLengths)
      int nCats = intWalkWindow;
      if (nCats < 2) nCats = 2;
      if (!local_dirichlet_impl(state->relBrLengths,
                                state->parent, state->child,
                                nCats, scaleTuning, logHastings,
                                state->brSnapshot,
                                state->dirEdges)) {
        return false;
      }
      break;
    }
    // Slice samplers. run_mcmc_batch_cpp intercepts these before reaching
    // do_move_impl, so that path carries its own adapted per-move width;
    // here the caller supplies the width through scaleTuning.
    case 19:
      return slice_scalar_impl(data, state, charIdx, scaleTuning, beta);
    case 29:
      return slice_kprime_hyper_impl(data, state, charIdx, scaleTuning);

    case 25: { // gibbs_kprime_sweep — Gibbs update of all k'_i
      return gibbs_kprime_sweep_impl(data, state, beta);
    }
    case 26: { // block_kprime_shift — shift all k'_i by same delta
      return block_kprime_shift_impl(data, state, intWalkWindow, beta);
    }
    // Move types 27 (scale_kprime_alpha) and 28 (scale_kprime_beta) — the
    // axis-aligned Bactrian moves on raw α and β — were removed when the
    // BG hyperparameter sampler was reparameterised to (s, r) coordinates.
    // The (s, r) univariate slice (move type 29) decorrelates the BG ridge
    // that those axis moves could not traverse. R-side .BuildMoves no
    // longer registers either move, so these case labels would be
    // unreachable if added back.
    case 30: { // mh_logit_p — logit-scale MH on p (for empirical_geometric)
      // Multiplicative MH on p ∈ (0,1) overshoots when p is close to 1, which
      // is the typical posterior region under the empirical_geometric prior.
      // Propose on the unbounded logit scale instead, so no rejections from
      // boundary violations. Jacobian is |dp/dlogit(p)| = p (1 − p).
      if (oldP <= 0.0 || oldP >= 1.0) return false;
      double logitP = std::log(oldP / (1.0 - oldP));
      double logitPnew = logitP + scaleTuning * bactrian_perturbation();
      double newP;
      if (logitPnew >= 0.0) {
        newP = 1.0 / (1.0 + std::exp(-logitPnew));
      } else {
        double e = std::exp(logitPnew);
        newP = e / (1.0 + e);
      }
      if (newP <= 0.0 || newP >= 1.0) return false;
      state->p = newP;
      logHastings = std::log(newP) + std::log1p(-newP)
                  - std::log(oldP) - std::log1p(-oldP);
      break;
    }
    case 35: { // gibbs_p_marginal — data-augmentation Metropolis-within-Gibbs
      // p-update for the marginal_k geometric arm. Self-contained (does its own
      // accept/reject and returns), like the Gibbs cases — it does NOT fall
      // through to the generic MH machinery, because the accept ratio is the
      // truncation-normaliser ratio Σ_i[logZ_i(p) − logZ_i(p*)], not the
      // marginal-LL ratio. Derivation + math-prover verification + numerical
      // invariance check: dev/red-team/proofs/marginal-k-gibbs-p.md and
      // dev/red-team/numerical/gibbs-p-identity-check.R.
      if (!data->marginalK) return false;
      // Geometric arm only (the other priors carry no scalar conjugate p).
      {
        const bool isPlainGeom = !data->kPriorLogseries &&
                                 !data->kPriorBetaGeometric &&
                                 !data->kPriorEmpiricalGeometric;
        if (!isPlainGeom) return false;
      }
      const int nTrans = (int)data->transIdxGlobal.size();
      if (nTrans == 0) return false;
      // Precondition: a WARM cache (raw per-(char,k') LLs against the committed
      // tree/rates). If a preceding non-p move left it cold, skip this iteration
      // and let mh_logit_p (case 30) move p; never impute from a stale cache.
      if (!state->charLLCacheReady ||
          (int)state->charLLNCand.size() != nTrans ||
          (int)state->charLLCache.size() != nTrans * kMaxKprimeCand)
        return false;
      if (!(oldP > 0.0 && oldP < 1.0)) return false;

      const double logP   = std::log(oldP);
      const double log1mP = std::log1p(-oldP);
      const bool   uncond = data->unconditionalPrior;
      const int    K      = data->kprimeTruncK;

      // --- Step 1: impute u_i ~ Categorical(charLogW[i,·]) over the cached
      //     support {0..nEff_i-1}; accumulate S = Σ u_i and (Model A) the shift
      //     c_A = Σ (kObs_i - 2). charLogW[i,c] = charLLCache[i,c] + logP +
      //     c·log1mP (the evaluator's per-char logSumExp argument); g_i(p) and
      //     Z_i(p) are u-independent and cancel out of the categorical.
      double sumU = 0.0;
      double cA   = 0.0;
      for (int i = 0; i < nTrans; ++i) {
        const int gi    = data->transIdxGlobal[i];
        const int kObsi = data->kObs[gi];
        const int nEff  = state->charLLNCand[i];
        if (nEff <= 0) return false;  // empty support (kObs_i > K): degenerate
        const double* rawRow =
          &state->charLLCache[(size_t)i * kMaxKprimeCand];
        double mx = R_NegInf;
        for (int c = 0; c < nEff; ++c) {
          if (R_FINITE(rawRow[c])) {
            const double w = rawRow[c] + logP + c * log1mP;
            if (w > mx) mx = w;
          }
        }
        if (!R_FINITE(mx)) return false;  // no finite candidate for this char
        double sum = 0.0;
        for (int c = 0; c < nEff; ++c) {
          if (R_FINITE(rawRow[c]))
            sum += std::exp(rawRow[c] + logP + c * log1mP - mx);
        }
        const double targ = R::unif_rand() * sum;
        int ui = -1;  // chosen candidate; tracks the last FINITE one as fallback
        double acc = 0.0;
        for (int c = 0; c < nEff; ++c) {
          if (!R_FINITE(rawRow[c])) continue;   // zero-prob candidate: never select
          acc += std::exp(rawRow[c] + logP + c * log1mP - mx);
          ui = c;                                // last finite candidate (FP fallback)
          if (acc >= targ) break;                // inverse-CDF selection
        }
        // ui >= 0 guaranteed: mx finite => at least one finite candidate exists,
        // and only finite candidates are ever assigned to ui.
        sumU += (double)ui;
        if (uncond) cA += (double)(kObsi - 2);
      }

      // --- Step 2: propose p* from the untruncated conjugate Beta.
      const double shape1 = data->kprimeHyperA + (double)nTrans;
      const double shape2 = data->kprimeHyperB + sumU + cA;
      const double pStar  = R::rbeta(shape1, shape2);
      if (!(pStar > 0.0 && pStar < 1.0)) return false;  // boundary draw → reject

      // --- Accept ratio: log α = Σ_i[ logZ_i(oldP) − logZ_i(pStar) ].
      const double log1mPstar = std::log1p(-pStar);
      double sumLogZ_old, sumLogZ_new;
      if (uncond) {
        // Model A: Z = 1 − (1−p)^(K−1), shared across all trans chars.
        sumLogZ_old =
          (double)nTrans * std::log1p(-std::exp((K - 1) * log1mP));
        sumLogZ_new =
          (double)nTrans * std::log1p(-std::exp((K - 1) * log1mPstar));
      } else {
        // Model B: Z_i = 1 − (1−p)^(K − kObs_i + 1), per character.
        sumLogZ_old = 0.0;
        sumLogZ_new = 0.0;
        for (int i = 0; i < nTrans; ++i) {
          const int kObsi = data->kObs[data->transIdxGlobal[i]];
          const int e     = K - kObsi + 1;
          sumLogZ_old += std::log1p(-std::exp(e * log1mP));
          sumLogZ_new += std::log1p(-std::exp(e * log1mPstar));
        }
      }
      // Deep small-p tail: Z → 0, logZ → −∞. The non-finite guard is symmetric
      // in (oldP, pStar) so detailed balance is preserved; mh_logit_p covers
      // that region. (See proof §6.4.)
      if (!R_FINITE(sumLogZ_old) || !R_FINITE(sumLogZ_new)) return false;
      const double logAlpha = sumLogZ_old - sumLogZ_new;

      if (std::log(R::unif_rand()) < logAlpha) {
        state->p = pStar;
        // The charLL cache's candidate SUPPORT (nEff per char) is p-dependent:
        // it is the early-termination set chosen at FILL-TIME p. The raw LLs are
        // p-independent, but a large downward p-jump fattens the geometric tail
        // P(u|p)=p(1-p)^u, moving non-negligible mass onto candidates the warm
        // cache never stored — so the warm fast-path would UNDERCOUNT. Force a
        // COLD re-fill at p* (re-derives the support via pruning + M-164
        // early-termination at p*): this gives a correct committed logLik AND
        // leaves the cache appropriate for the NEXT move's imputation. Same
        // force-cold idiom as the FREEZE-003 weighted-move re-enable. (The warm
        // shortcut is only valid for mh_logit_p's small steps; case 30 keeps it.)
        state->charLLCacheReady = false;
        state->logLik   = compute_full_loglik(*data, *state);
        state->logPrior = compute_log_prior(*data, *state);
        return true;   // accept
      }
      return false;    // reject: p, logLik, logPrior, cache all untouched
    }
    case 31: { // scale_class_rate_log_sd - Bactrian multiplicative move.
      // charIdx carries the 1-based class index supplied by run_mcmc_batch_cpp.
      // Two arms:
      //   - useHyperpriorOnSigma: scale z_c by mult, recompute σ_c = τ · z_c.
      //   - legacy: scale σ_c directly.
      int cIdx0 = charIdx - 1;  // convert to 0-based
      if (cIdx0 < 0 || cIdx0 >= (int)state->classRateLogSd.size()) return false;
      classIdx31  = cIdx0;
      oldClassRLS = state->classRateLogSd[cIdx0];
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      if (state->useHyperpriorOnSigma) {
        if (cIdx0 >= (int)state->classZ.size()) return false;
        oldClassZ = state->classZ[cIdx0];
        state->classZ[cIdx0]         = oldClassZ * mult;
        state->classRateLogSd[cIdx0] = state->hyperTau * state->classZ[cIdx0];
        if (state->classZ[cIdx0] <= 0.0 ||
            state->classRateLogSd[cIdx0] <= 0.0) {
          state->classZ[cIdx0]         = oldClassZ;
          state->classRateLogSd[cIdx0] = oldClassRLS;
          return false;
        }
        // Lockstep: keep scalar rateLogSd in sync with σ_0 for partial-CL paths.
        if (cIdx0 == 0) state->rateLogSd = state->classRateLogSd[0];
      } else {
        state->classRateLogSd[cIdx0] = oldClassRLS * mult;
        if (state->classRateLogSd[cIdx0] <= 0.0) {
          state->classRateLogSd[cIdx0] = oldClassRLS;
          return false;
        }
        if (cIdx0 == 0) state->rateLogSd = state->classRateLogSd[0];
      }
      // Lockstep: keep state->rateLogSd in sync with classRateLogSd[0] so
      // the legacy prior term and the trace's `rate_log_sd` column stay
      // coherent with the value the partitioned likelihood actually uses.
      if (cIdx0 == 0) state->rateLogSd = state->classRateLogSd[0];
      logHastings = std::log(mult);
      break;
    }
    case 32: { // dirichlet_simplex_class_w - Dirichlet proposal on classW simplex
      int nC = (int)state->classW.size();
      if (nC < 2) return false;
      classWSnapshot = clone(state->classW);
      std::vector<int> dummyEdges;
      NumericVector tmpSnap = clone(state->classW);
      if (!dirichlet_simplex_impl(state->classW, nC,
                                  betaSimplexTuning, logHastings,
                                  tmpSnap, dummyEdges)) {
        state->classW = classWSnapshot;
        return false;
      }
      // Recompute classRate[c] = classW[c] * nChar / nChar_c  (section 5.1)
      int nChar = 0;
      for (int k = 0; k < (int)state->nCharPerClass.size(); ++k)
        nChar += state->nCharPerClass[k];
      for (int k = 0; k < nC; ++k) {
        int nk = (k < (int)state->nCharPerClass.size()) ? state->nCharPerClass[k] : 1;
        state->classRate[k] = (nk > 0) ? state->classW[k] * nChar / nk : 1.0;
      }
      // NOTE: Dirichlet prior on classW is pending the parallel agent
      // cpp_log_prior extension. Acceptance is LL-ratio only until reconciled.
      break;
    }
    case 33: { // joint_tl_rn (tree_length × rate_neo): partition-rate ridge
      // from Issue 1 fix. r' = r·exp(σ z2), T' = T·exp(σ z1) with correlated
      // 2D Bactrian; ρ adapted from warmup posterior cor(log T, log r).
      double z1, z2;
      bactrian_2d_perturbation(jointRho, z1, z2);
      double mult1 = std::exp(scaleTuning * z1);
      double mult2 = std::exp(scaleTuning * z2);
      state->treeLength = oldTL * mult1;
      state->rateNeo    = oldRN * mult2;
      logHastings = std::log(mult1) + std::log(mult2);
      break;
    }
    case 34: { // scale_hyper_tau — Bactrian scale on the population scale τ.
      // Only fired when state->useHyperpriorOnSigma is true. Scaling τ
      // rescales every σ_c = τ · z_c, so all per-class likelihoods need a
      // full reeval (CL cache invalidation handled in the accept block).
      if (!state->useHyperpriorOnSigma) return false;
      if (state->classZ.size() != state->classRateLogSd.size()) return false;
      case34Active   = true;
      oldTau34       = state->hyperTau;
      oldClassRLS34  = clone(state->classRateLogSd);
      oldRateLogSd34 = state->rateLogSd;
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      double newTau = oldTau34 * mult;
      if (newTau <= 0.0) return false;
      state->hyperTau = newTau;
      for (int k = 0; k < state->classRateLogSd.size(); ++k) {
        state->classRateLogSd[k] = newTau * state->classZ[k];
        if (state->classRateLogSd[k] <= 0.0) {
          // Roll back partial update before reporting rejection.
          state->hyperTau     = oldTau34;
          state->classRateLogSd = oldClassRLS34;
          state->rateLogSd    = oldRateLogSd34;
          return false;
        }
      }
      state->rateLogSd = state->classRateLogSd[0];
      logHastings = std::log(mult);
      break;
    }
    default:
      return false;
  }

  if (!R_FINITE(logHastings)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS; state->kprimeAlpha = oldKpA; state->kprimeBeta = oldKpB;
    if (bsIdx1 >= 0) {
      state->relBrLengths[bsIdx1] = bsOldVal1;
      state->relBrLengths[bsIdx2] = bsOldVal2;
    }
    if (moveType == 7 && kPrimeCharIdx >= 0) state->kPrime[kPrimeCharIdx] = oldKPrimeVal;
    if (moveType == 23 || moveType == 24) {
      const int nE = state->relBrLengths.size();
      for (int i = 0; i < nE; ++i) state->relBrLengths[i] = state->brSnapshot[i];
    }
    if (nniInPlace) { state->parent[nniCRow] = nniSavedP_cRow; state->parent[nniWRow] = nniSavedP_wRow; }
    if (classIdx31 >= 0) {
      state->classRateLogSd[classIdx31] = oldClassRLS;
      if (state->useHyperpriorOnSigma &&
          classIdx31 < (int)state->classZ.size()) {
        state->classZ[classIdx31] = oldClassZ;
      }
    }
    if (case34Active) {
      // case 34 (scale_hyper_tau) full restore — supersedes the
      // oldClassRLS0 lockstep mirror since this case rebuilds the entire
      // classRateLogSd vector.
      state->hyperTau     = oldTau34;
      state->classRateLogSd = oldClassRLS34;
      state->rateLogSd    = oldRateLogSd34;
    } else if (state->usePartitioned && state->classRateLogSd.size() > 0) {
      // Lockstep rollback: cases 2/21 mirror rateLogSd → classRateLogSd[0].
      state->classRateLogSd[0] = oldClassRLS0;
    }
    if (moveType == 32 && !classWSnapshot.isNULL()) {
      state->classW = classWSnapshot;
      int nChar32 = 0;
      for (int k = 0; k < (int)state->nCharPerClass.size(); ++k) nChar32 += state->nCharPerClass[k];
      for (int k = 0; k < (int)state->classW.size(); ++k) {
        int nk = (k < (int)state->nCharPerClass.size()) ? state->nCharPerClass[k] : 1;
        state->classRate[k] = (nk > 0) ? state->classW[k] * nChar32 / nk : 1.0;
      }
    }
    return false;
  }

  // NNI doesn't change any parameter → prior is unchanged; skip evaluation.
  double newLogPrior;
  if (nniInPlace) {
    newLogPrior = state->logPrior;
  } else {
    newLogPrior = compute_log_prior(*data, *state);
  }

  if (!R_FINITE(newLogPrior)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS; state->kprimeAlpha = oldKpA; state->kprimeBeta = oldKpB;
    if (bsIdx1 >= 0) {
      state->relBrLengths[bsIdx1] = bsOldVal1;
      state->relBrLengths[bsIdx2] = bsOldVal2;
    }
    if (moveType == 7 && kPrimeCharIdx >= 0) state->kPrime[kPrimeCharIdx] = oldKPrimeVal;
    if (moveType == 23 || moveType == 24) {
      const int nE = state->relBrLengths.size();
      for (int i = 0; i < nE; ++i) state->relBrLengths[i] = state->brSnapshot[i];
    }
    if (nniInPlace) { state->parent[nniCRow] = nniSavedP_cRow; state->parent[nniWRow] = nniSavedP_wRow; }
    if (classIdx31 >= 0) {
      state->classRateLogSd[classIdx31] = oldClassRLS;
      if (state->useHyperpriorOnSigma &&
          classIdx31 < (int)state->classZ.size()) {
        state->classZ[classIdx31] = oldClassZ;
      }
    }
    if (case34Active) {
      // case 34 (scale_hyper_tau) full restore — supersedes oldClassRLS0.
      state->hyperTau     = oldTau34;
      state->classRateLogSd = oldClassRLS34;
      state->rateLogSd    = oldRateLogSd34;
    } else if (state->usePartitioned && state->classRateLogSd.size() > 0) {
      // Lockstep rollback: cases 2/21 mirror rateLogSd → classRateLogSd[0].
      state->classRateLogSd[0] = oldClassRLS0;
    }
    if (moveType == 32 && !classWSnapshot.isNULL()) {
      state->classW = classWSnapshot;
      int nChar32 = 0;
      for (int k = 0; k < (int)state->nCharPerClass.size(); ++k) nChar32 += state->nCharPerClass[k];
      for (int k = 0; k < (int)state->classW.size(); ++k) {
        int nk = (k < (int)state->nCharPerClass.size()) ? state->nCharPerClass[k] : 1;
        state->classRate[k] = (nk > 0) ? state->classW[k] * nChar32 / nk : 1.0;
      }
    }
    return false;
  }

  // OPP-6: select evaluation topology — proposed values for NNI/SPR,
  // current state for all other moves.
  const IntegerVector& evalParent = topologyChanged ? proposedParent : state->parent;
  const IntegerVector& evalChild  = topologyChanged ? proposedChild  : state->child;
  const NumericVector& evalRelBr  = topologyChanged ? proposedRelBr  : state->relBrLengths;

  // ---- Likelihood evaluation (M-064: partial, M-065: vectors, M-121: node CL) ----
  // Moves that only touch `p` (case 8 legacy multiplicative, case 30 logit MH)
  // leave the likelihood untouched UNDER SAMPLED-K: there p enters only the
  // prior, since kPrime is an explicit state and L = L(y | tree, mu, kPrime).
  //
  // Under MARGINAL-K (FU-5 fix, 2026-05-28), this is WRONG: p enters the
  // per-character marginal LL via the P(u | p) = p (1 - p)^u weights
  // consumed inside the per-character logSumExp. A p-move therefore changes
  // the marginal LL even though it changes no other state. The marginal
  // evaluator (compute_full_loglik_at -> cpp_log_likelihood_marginal) is
  // p-aware; we just need to take the recompute path. Under the buggy
  // `likChanges = (moveType != 8 && moveType != 30)`, case-30 reused
  // state->logLik in the MH ratio and the chain random-walked on p
  // ignoring data (flat-prior + zero-LL → guaranteed acceptance mod
  // Jacobian). See dev/red-team/heavy-tests/marginal-k-* and the FU-1
  // side-finding for the diagnostic.
  bool likChanges = data->marginalK || (moveType != 8 && moveType != 30);
  // Under marginal-k mode the per-partition cache (partLogLik) is built
  // for sampled-k semantics (cpp_partition_log_likelihood with fixed
  // state->kPrime). The marginal evaluator does NOT update it. Force
  // hasPLC=false so we route through compute_full_loglik_at (the
  // marginal-aware dispatcher).
  bool hasPLC = !state->partLogLik.empty() && !data->marginalK;
  double newLogLik;
  std::vector<double> newPC;
  bool usedPartialCL = false;

  if (!likChanges) {
    newLogLik = state->logLik;
  } else if (nniInPlace && state->nodeCL.ready()) {
    // M-121: NNI with valid node CL cache → partial evaluation
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    // Update TreeNav topology for the NNI swap (before partial eval)
    update_topo_nni(state->nodeCL.topo, nniVNode, nniUNode,
                    nniCNode, nniWNode);

    auto dirty = find_dirty_nni(state->nodeCL.topo, nniVNode, nniUNode);
    newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                    evalParent, evalChild, propEdgeLen,
                                    state->rateLoss, state->rateNeo,
                                    state->rateLogSd, state->betaScale, dirty);
    usedPartialCL = true;
    state->diagNniPartialCount++;

  } else if (moveType == 4 && state->nodeCL.ready()) {
    // M-121: beta_simplex with valid node CL cache → partial evaluation
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    // M-158: sync nodeEdgeLen for the two modified edges (save for rollback)
    bsOldNodeEdgeLen1 = state->nodeCL.topo.nodeEdgeLen[evalChild[bsIdx1]];
    bsOldNodeEdgeLen2 = state->nodeCL.topo.nodeEdgeLen[evalChild[bsIdx2]];
    state->nodeCL.topo.nodeEdgeLen[evalChild[bsIdx1]] = propEdgeLen[bsIdx1];
    state->nodeCL.topo.nodeEdgeLen[evalChild[bsIdx2]] = propEdgeLen[bsIdx2];

    auto dirty = find_dirty_beta_simplex(state->nodeCL.topo,
                                          evalParent, bsIdx1, bsIdx2);
    newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                    evalParent, evalChild, propEdgeLen,
                                    state->rateLoss, state->rateNeo,
                                    state->rateLogSd, state->betaScale, dirty);
    usedPartialCL = true;
    state->diagBsPartialCount++;

  } else if ((moveType == 23 || moveType == 24) && state->nodeCL.ready()) {
    // M-127: Dirichlet (random or local) with valid node CL cache → partial eval
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    // M-158: sync nodeEdgeLen for all modified Dirichlet edges
    for (int idx : state->dirEdges)
      state->nodeCL.topo.nodeEdgeLen[evalChild[idx]] = propEdgeLen[idx];

    auto dirty = find_dirty_dirichlet(state->nodeCL.topo,
                                       evalParent, state->dirEdges);
    // Heuristic: if dirty set covers most of the tree, fall back to full eval
    int nInternal = nEdge + 1 - data->nTip;
    if ((int)dirty.size() > (int)(0.8 * (nInternal + data->nTip))) {
      state->diagDirFullbackCount++;
      { ClWorkspace* dWs = state->clWs.ready() ? &state->clWs : nullptr;
        newLogLik = state->usePartitioned
          ? cpp_log_likelihood_partitioned(
              *data, evalParent, evalChild, propEdgeLen,
              state->kPrime, state->rateLoss,
              state->classRateLogSd, state->classRate,
              state->etaNeo, state->betaScale, dWs)
          : cpp_log_likelihood(
              *data, evalParent, evalChild, propEdgeLen,
              state->kPrime, state->rateLoss, state->rateLogSd,
              state->rateNeo, state->betaScale, dWs);
      }
    } else {
      newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                      evalParent, evalChild, propEdgeLen,
                                      state->rateLoss, state->rateNeo,
                                      state->rateLogSd, state->betaScale, dirty);
      usedPartialCL = true;

      state->diagDirPartialCount++;
    }

  } else if (sprPartialCL) {
    // M-158: SPR with partial CL evaluation via TreeNav
    // 1. Apply SPR to TreeNav (updates topology + nodeEdgeLen)
    update_topo_spr(state->nodeCL.topo, sprMeta);

    // 2. Build proposed parent/child/relBr from updated TreeNav
    treenav_to_preorder(state->nodeCL.topo, state->treeLength,
                        proposedParent, proposedChild, proposedRelBr);
    topologyChanged = true;
    const IntegerVector& sprParent = proposedParent;
    const IntegerVector& sprChild  = proposedChild;

    int nEdge = proposedRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * proposedRelBr[i];

    // 3. Find dirty set and check if partial eval is worthwhile
    auto dirty = find_dirty_spr(state->nodeCL.topo,
                                 sprMeta.u, sprMeta.p, sprMeta.a);
    int nInternal = nEdge + 1 - data->nTip;
    if ((int)dirty.size() > (int)(0.8 * (nInternal + data->nTip))) {
      // Dirty set too large — fall back to full eval
      // (TreeNav already updated; will be reversed on rejection)
      { ClWorkspace* sWs = state->clWs.ready() ? &state->clWs : nullptr;
        newLogLik = state->usePartitioned
          ? cpp_log_likelihood_partitioned(
              *data, sprParent, sprChild, propEdgeLen,
              state->kPrime, state->rateLoss,
              state->classRateLogSd, state->classRate,
              state->etaNeo, state->betaScale, sWs)
          : cpp_log_likelihood(
              *data, sprParent, sprChild, propEdgeLen,
              state->kPrime, state->rateLoss, state->rateLogSd,
              state->rateNeo, state->betaScale, sWs);
      }
    } else {
      // 4. Partial CL evaluation
      newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                      sprParent, sprChild, propEdgeLen,
                                      state->rateLoss, state->rateNeo,
                                      state->rateLogSd, state->betaScale, dirty);
      usedPartialCL = true;
    }

  } else if (!hasPLC) {
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];
    // Dispatch through compute_full_loglik_at so marginal-k mode routes to
    // cpp_log_likelihood_marginal automatically.
    newLogLik = compute_full_loglik_at(
      *data, *state, evalParent, evalChild, propEdgeLen);
  } else {
    int nParts = (int)data->parts.size();
    int nEdge  = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];
    newPC = state->partLogLik;
    switch (moveType) {
      case 1: {
        // rate_loss: enters mkn stationary frequencies + Q-matrix; only
        // neomorphic partitions change. Partial recompute is safe.
        ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
        newLogLik = state->logLik;
        for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
          int pi = data->neoPartIndices[ni];
          double v = cpp_partition_log_likelihood(*data, pi,
            evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
            state->betaScale, wsPtr);
          newLogLik += (v - newPC[pi]);
          newPC[pi] = v;
        }
        break;
      }
      // case 3 (rate_neo scale) intentionally falls through to the default
      // full-recompute branch: under the RB-style partition-rate normalisation
      // (audit Issue 1), rate_neo shifts BOTH neoScale AND transScale, so
      // every partition's log-likelihood is stale — not just neo.
      case 7: {
        ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
        int ap = data->charToPartition[charIdx];
        newLogLik = state->logLik;
        if (ap >= 0) {
          double v = cpp_partition_log_likelihood(*data, ap,
            evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
            state->betaScale, wsPtr);
          newLogLik += (v - newPC[ap]);
          newPC[ap] = v;
        }
        break;
      }
      default: {
        ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;
        newLogLik = 0.0;
        for (int pi = 0; pi < nParts; ++pi) {
          newPC[pi] = cpp_partition_log_likelihood(*data, pi,
            evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
            state->betaScale, wsPtr);
          newLogLik += newPC[pi];
        }
        break;
      }
    }
  }

  double logAlpha = beta * (newLogLik - state->logLik) +
                    (newLogPrior - state->logPrior) + logHastings;
  if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    if (!newPC.empty()) state->partLogLik = std::move(newPC);
    // OPP-6: commit proposed topology to state on acceptance
    if (topologyChanged) {
      state->parent       = std::move(proposedParent);
      state->child        = std::move(proposedChild);
      state->relBrLengths = std::move(proposedRelBr);
    }
    // OPP-6b: in-place NNI — state->parent already modified in-place.
    // Partition cache handled by std::move(newPC) above (recomputed via
    // default branch of the partial-lik switch).

    // MARGINAL-K-FREEZE-003 coherence guard (opt-in: compile with
    // -DMKPRIME_CHECK_MARGINAL_COHERENCE). After an accepted move under
    // marginal_k the committed state->logLik MUST equal a fresh COLD marginal
    // recompute of the committed tree; a mismatch is the freeze / wrong-
    // posterior failure mode (a move left a stale baseline). Off by default
    // (zero production cost); the testthat free-topology test is the always-on
    // CI guard. Covers the MH-eval moves that reach this accept block; the
    // early-return Gibbs moves are disabled under marginal_k (see .BuildMoves).
#ifdef MKPRIME_CHECK_MARGINAL_COHERENCE
    if (data->marginalK) {
      bool savedReady = state->charLLCacheReady;
      state->charLLCacheReady = false;                       // force cold rebuild
      double freshLL = compute_full_loglik(*data, *state);
      state->charLLCacheReady = savedReady;
      if (std::abs(state->logLik - freshLL) > 1e-6)
        Rf_warning(
          "MARGINAL-K coherence: moveType %d logLik %.8g != fresh %.8g (gap %.4g)",
          moveType, state->logLik, freshLL, state->logLik - freshLL);
    }
#endif

    // M-121/M-161: cache management on acceptance.
    // Partial CL moves keep the cache valid (already updated).
    // Other moves: granular invalidation by move type.
    if (!usedPartialCL && likChanges) {
      switch (moveType) {
        case 1:
          // rate_loss: only neomorphic units affected (mkn stationary +
          // Q-matrix; trans/known unchanged)
          state->nodeCL.invalidate_neo_cls();
          break;
        case 3:
          // rate_neo: audit Issue 1 — partition-rate normalisation makes
          // rateNeo affect BOTH neo and trans unit rateScales, so all CLs
          // are stale.
          state->nodeCL.invalidate_all_cls();
          break;
        case 2:
        case 31:  // scale_class_rate_log_sd: σ_c change → ACRV rates change
        case 34:  // scale_hyper_tau: rescales all σ_c → ACRV rates change
          // rateLogSd: ACRV rates change, all units stale
          state->nodeCL.invalidate_all_cls();
          break;
        case 7:
          // kPrime int_walk: unit structure may change
          state->nodeCL.invalidate_structure();
          break;
        case 16:
          // beta_scale: Q-het parameter, all units stale
          state->nodeCL.invalidate_all_cls();
          break;
        default:
          // tree_length (0), topology, etc.: full invalidation
          state->nodeCL.invalidate_all();
          break;
      }
    }
    // When partial CL was used, clear partition-level cache
    // (it's not maintained by partial eval; will be rebuilt if needed)
    if (usedPartialCL && !state->partLogLik.empty())
      state->partLogLik.clear();

    return true;
  }

  // Reject: rollback
  state->treeLength   = oldTL;
  state->rateLoss     = oldRL;
  state->rateLogSd    = oldRLSD;
  state->rateNeo      = oldRN;
  state->p            = oldP;
  state->betaScale    = oldBS;
  state->kprimeAlpha  = oldKpA;
  state->kprimeBeta   = oldKpB;
  if (bsIdx1 >= 0) {
    state->relBrLengths[bsIdx1] = bsOldVal1;
    state->relBrLengths[bsIdx2] = bsOldVal2;
  }
  if (moveType == 7 && kPrimeCharIdx >= 0) state->kPrime[kPrimeCharIdx] = oldKPrimeVal;
  // M-125/M-127: Dirichlet simplex rollback — restore full vector from snapshot
  if (moveType == 23 || moveType == 24) {
    const int nE = state->relBrLengths.size();
    for (int i = 0; i < nE; ++i) state->relBrLengths[i] = state->brSnapshot[i];
  }
  // OPP-6b: in-place NNI rollback — restore 2 parent values
  if (nniInPlace) { state->parent[nniCRow] = nniSavedP_cRow; state->parent[nniWRow] = nniSavedP_wRow; }
  // Per-class rollback (cases 31, 32, 34) + lockstep mirror for cases 2/21.
  if (classIdx31 >= 0) {
    state->classRateLogSd[classIdx31] = oldClassRLS;
    if (state->useHyperpriorOnSigma &&
        classIdx31 < (int)state->classZ.size()) {
      state->classZ[classIdx31] = oldClassZ;
    }
  }
  if (case34Active) {
    // case 34 (scale_hyper_tau) full restore — supersedes oldClassRLS0.
    state->hyperTau     = oldTau34;
    state->classRateLogSd = oldClassRLS34;
    state->rateLogSd    = oldRateLogSd34;
  } else if (state->usePartitioned && state->classRateLogSd.size() > 0) {
    // Lockstep rollback: cases 2/21 mirror rateLogSd → classRateLogSd[0].
    state->classRateLogSd[0] = oldClassRLS0;
  }
  if (moveType == 32 && !classWSnapshot.isNULL()) {
    state->classW = classWSnapshot;
    // Restore classRate from snapshot
    int nChar = 0;
    for (int k = 0; k < (int)state->nCharPerClass.size(); ++k) nChar += state->nCharPerClass[k];
    for (int k = 0; k < (int)state->classW.size(); ++k) {
      int nk = (k < (int)state->nCharPerClass.size()) ? state->nCharPerClass[k] : 1;
      state->classRate[k] = (nk > 0) ? state->classW[k] * nChar / nk : 1.0;
    }
  }

  // M-121: rollback node CL cache on rejection of partial-eval moves
  if (usedPartialCL) {
    restore_dirty_cls(state->nodeCL, state->nodeCL.dirtyNodes);
    // Also rollback TreeNav for NNI (topology was updated before partial eval)
    if (nniInPlace) {
      update_topo_nni(state->nodeCL.topo, nniVNode, nniUNode,
                      nniWNode, nniCNode);  // reverse swap
    }
    // M-158: rollback nodeEdgeLen for beta_simplex
    if (moveType == 4 && bsIdx1 >= 0) {
      state->nodeCL.topo.nodeEdgeLen[state->child[bsIdx1]] = bsOldNodeEdgeLen1;
      state->nodeCL.topo.nodeEdgeLen[state->child[bsIdx2]] = bsOldNodeEdgeLen2;
    }
    // M-158: rollback nodeEdgeLen for Dirichlet (recompute from restored brSnapshot)
    if (moveType == 23 || moveType == 24) {
      for (int idx : state->dirEdges)
        state->nodeCL.topo.nodeEdgeLen[state->child[idx]] =
          state->treeLength * state->relBrLengths[idx];
    }
  }
  // M-158: rollback SPR TreeNav on rejection (whether or not partial CL was used)
  if (sprPartialCL) {
    reverse_topo_spr(state->nodeCL.topo, sprMeta);
    if (usedPartialCL) {
      restore_dirty_cls(state->nodeCL, state->nodeCL.dirtyNodes);
    }
  }

  // Marginal-k: on rejection of a non-30 move, the cache contents are
  // stale (they were filled against the proposed parent/child/edgeLen).
  // Invalidate so the next marginal eval rebuilds against the rolled-
  // back state. For case 30 (p-only) the cache was not written during
  // the proposal (cache path skipped re-pruning) and is independent of
  // p, so leave it valid.
  //
  // FU-3 Tier 2 mirrors charLLCacheReady on rejection: any subtree CLs
  // written during partial-CL evaluation (FU-3b future work) are
  // invalidated together with the per-(char, k') cache, since both were
  // computed against the proposed (parent, child, edgeLen).
  if (data->marginalK && moveType != 30 && moveType != 35) {
    state->charLLCacheReady = false;
    state->invalidate_per_kp_cl_all();
  }

  return false;
}


// [[Rcpp::export]]
bool do_move_cpp(SEXP dataPtr, SEXP statePtr,
                 int moveType, int charIdx,
                 double scaleTuning, double betaSimplexTuning,
                 int intWalkWindow, double beta) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  return do_move_impl(data, state, moveType, charIdx,
                      scaleTuning, betaSimplexTuning, intWalkWindow, beta);
}


// Moves that can score through the node CL cache: beta_simplex, NNI, SPR,
// dirichlet_branch, local_dirichlet.
static inline bool is_partial_cl_move(int moveType) {
  return moveType == 4 || moveType == 5 || moveType == 6 ||
         moveType == 23 || moveType == 24;
}

// ---------------------------------------------------------------------------
// run_mcmc_batch_cpp: C++ inner loop — runs nBatch iterations for one run.
//
// Handles weighted move selection, do_move_impl calls, chain swaps, and
// sample collection. R calls this per-batch (default 200 iters) and handles
// adaptation, convergence checks, progress, and checkpointing at boundaries.
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
List run_mcmc_batch_cpp(
    SEXP dataPtr,
    List stateXPtrs,
    NumericVector betas,
    IntegerVector moveTypeCodes,
    IntegerVector transIdxCpp,
    IntegerVector sliceParamCodes,
    NumericVector moveWeights,
    NumericMatrix chainScaleTunings,
    NumericVector chainBsmpTunings,
    IntegerVector chainIntWalkWins,
    IntegerVector moveIntParams,
    NumericMatrix sliceWidths,
    NumericMatrix jointRhos,
    int nBatch,
    int startIter,
    int warmup,
    int thin,
    bool hasNeo,
    int nEdge,
    double cacheBonus = 1.0
) {
  McmcData* data = Rcpp::XPtr<McmcData>(dataPtr).get();
  int nChains    = stateXPtrs.size();
  int nMoves     = moveTypeCodes.size();
  int nTrans     = transIdxCpp.size();

  // Extract raw state pointers
  std::vector<McmcState*> states(nChains);
  for (int ch = 0; ch < nChains; ++ch)
    states[ch] = Rcpp::XPtr<McmcState>(stateXPtrs[ch]).get();

  // Cumulative move weights for O(nMoves) weighted sampling
  // Base weights (used when node CL cache is invalid)
  std::vector<double> cumWeights(nMoves);
  double totalWeight = 0.0;
  for (int m = 0; m < nMoves; ++m) {
    totalWeight += moveWeights[m];
    cumWeights[m] = totalWeight;
  }

  // M-159: Cache-boosted weights, used on the iteration after a
  // partial-CL-eligible move (see McmcState::cacheBoostNext). The node CL
  // cache is never populated under qHeterogeneity or marginalK.
  bool haveCacheBoost = cacheBonus > 1.0 && !data->qHeterogeneity &&
                        !data->marginalK;
  std::vector<double> cumWeightsCached(nMoves);
  double totalWeightCached = 0.0;
  if (haveCacheBoost) {
    for (int m = 0; m < nMoves; ++m) {
      int mt = moveTypeCodes[m];
      double w = moveWeights[m];
      if (is_partial_cl_move(mt)) w *= cacheBonus;
      totalWeightCached += w;
      cumWeightsCached[m] = totalWeightCached;
    }
  }
  int cacheHits = 0, cacheMisses = 0;

  // Accept/propose counters (nChains x nMoves)
  IntegerMatrix acceptCounts(nChains, nMoves);
  IntegerMatrix proposeCounts(nChains, nMoves);
  // Per-move wall-time tracking for adaptive scheduler (M-092)
  NumericMatrix moveTimeNs(nChains, nMoves);
  // Slice sampler stepping-out expansion counts (for width adaptation)
  IntegerMatrix sliceExpansions(nChains, nMoves);
  int nSwapPairs = std::max(0, nChains - 1);
  IntegerVector swapAccept(nSwapPairs, 0);
  IntegerVector swapPropose(nSwapPairs, 0);

  // Sample storage
  // Hyperparameter columns depend on kPrime prior:
  //   geometric: "p" (1 col)
  //   beta_geometric: "kprime_alpha", "kprime_beta" (2 cols)
  //   logseries: none (0 cols)
  bool includeBG = data->kPriorBetaGeometric;
  bool includeP  = !data->kPriorLogseries && !includeBG;
  int nKpHyperCols = includeBG ? 2 : (includeP ? 1 : 0);
  bool includeBS = data->qHeterogeneity;  // M-052: beta_scale column
  // Per-class partition columns (Layer 1).
  // Emit classRateLogSd columns iff any move in the schedule has type 31
  // (scale_class_rate_log_sd), which indicates shape is unlinked.
  // Emit classW columns iff any move has type 32 (dirichlet_simplex_class_w).
  // Emit hyper_tau + class_z columns iff state->useHyperpriorOnSigma
  // (move 33 is the τ move; state flag drives column emission).
  bool hasMove31 = false, hasMove32 = false;
  for (int m = 0; m < nMoves; ++m) {
    if (moveTypeCodes[m] == 31) hasMove31 = true;
    if (moveTypeCodes[m] == 32) hasMove32 = true;
  }
  int nClassRLS = (hasMove31 && nChains > 0 && states[0]->usePartitioned)
                  ? (int)states[0]->classRateLogSd.size() : 0;
  int nClassW   = (hasMove32 && nChains > 0 && states[0]->usePartitioned)
                  ? (int)states[0]->classW.size() : 0;
  bool hyperOn = (nChains > 0 && states[0]->usePartitioned &&
                  states[0]->useHyperpriorOnSigma);
  int nHyperTauCol = hyperOn ? 1 : 0;
  int nClassZ      = hyperOn ? (int)states[0]->classZ.size() : 0;
  // Base columns: log_post, log_lik, tree_length, rate_log_sd (4).
  // rate_loss included only when hasNeo (like rate_neo, p, beta_scale).
  // +2 diagnostic columns: swap_cold (cold-chain swaps since last sample),
  // topo_hash (topology fingerprint for change detection).
  int nScalarCols = 4 + (hasNeo ? 2 : 0) + nKpHyperCols +
                    (includeBS ? 1 : 0) + 2 + nTrans + nEdge +
                    nClassRLS + nClassW + nHyperTauCol + nClassZ;
  int maxSaved    = nBatch / thin + 2;
  std::vector<std::vector<double>> scalarRows;
  scalarRows.reserve(maxSaved);
  List edgeSamples;

  // Diagnostic: count accepted swaps involving the cold chain (index 0)
  int coldSwapsSinceSample = 0;

  // PT-RT-001: round-trip-time tracking. Each "particle" is an MCMC state
  // identity; particleAtSlot[s] gives which particle currently sits at
  // ladder slot s (slot 0 = cold = beta 1). On every accepted swap we
  // swap the two slot entries. We count a "round trip" each time a
  // particle returns to the cold slot after having reached the hot slot
  // since its last cold visit. Particle 0 starts at the cold slot, so
  // its initial state is PS_AT_COLD; particles initially at heated slots
  // start with state PS_NONE so they don't trivially count a partial trip.
  std::vector<int> particleAtSlot(nChains);
  for (int s = 0; s < nChains; ++s) particleAtSlot[s] = s;
  enum ParticleState : char { PS_NONE = 0, PS_AT_COLD = 1, PS_AT_HOT = 2 };
  std::vector<char> particleState(nChains, PS_NONE);
  if (nChains > 0) particleState[0] = PS_AT_COLD;
  int roundTripCount = 0;

  // Main iteration loop
  for (int i = 0; i < nBatch; ++i) {
    // Check for user interrupt every 10 iterations (expensive moves can take
    // seconds each, so we want to stay responsive to Ctrl-C / ESC).
    if (i % 10 == 0) R_CheckUserInterrupt();

    int iter = startIter + i;

    // Advance each chain
    for (int ch = 0; ch < nChains; ++ch) {
      bool useBoost = haveCacheBoost && states[ch]->cacheBoostNext;
      double tw = useBoost ? totalWeightCached : totalWeight;
      const auto& cw = useBoost ? cumWeightsCached : cumWeights;
      if (useBoost) ++cacheHits; else ++cacheMisses;

      double u = R::unif_rand() * tw;
      int moveIdx = 0;
      while (moveIdx < nMoves - 1 && u > cw[moveIdx]) ++moveIdx;

      proposeCounts(ch, moveIdx)++;

      int moveType = moveTypeCodes[moveIdx];

      // charIdx: for int_walk → random trans character; for slice → paramIdx;
      //          for scale_class_rate_log_sd (31) → 1-based classIdx from moveIntParams.
      int charIdx = 0;
      if (moveType == 7 && nTrans > 0) {
        int r = static_cast<int>(R::unif_rand() * nTrans);
        if (r >= nTrans) r = nTrans - 1;
        charIdx = transIdxCpp[r];
      } else if (moveType == 19) {
        charIdx = sliceParamCodes[moveIdx];
      } else if (moveType == 31) {
        charIdx = moveIntParams[moveIdx];  // 1-based classIdx
      }

      auto t0 = std::chrono::steady_clock::now();
      bool accepted;
      if (moveType == 19) {
        // Slice sampling — self-contained, no MH accept/reject
        int nExp = 0;
        accepted = slice_scalar_impl(
          data, states[ch], charIdx,
          sliceWidths(ch, moveIdx), betas[ch], 10, &nExp);
        sliceExpansions(ch, moveIdx) += nExp;
      } else if (moveType == 29) {
        // Prior-only slice sampler for BG hyperparameters (M-163)
        int nExp = 0;
        accepted = slice_kprime_hyper_impl(
          data, states[ch], sliceParamCodes[moveIdx],
          sliceWidths(ch, moveIdx), 10, &nExp);
        sliceExpansions(ch, moveIdx) += nExp;
      } else {
        // Per-move int param overrides chain-level intWalkWindow
        int iww = moveIntParams[moveIdx] > 0
                    ? moveIntParams[moveIdx]
                    : chainIntWalkWins[ch];
        accepted = do_move_impl(
          data, states[ch],
          moveType, charIdx,
          chainScaleTunings(ch, moveIdx),
          chainBsmpTunings[ch],
          iww,
          betas[ch],
          jointRhos(ch, moveIdx)
        );
      }
      auto t1 = std::chrono::steady_clock::now();
      moveTimeNs(ch, moveIdx) +=
        (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
          t1 - t0).count();
      if (accepted) acceptCounts(ch, moveIdx)++;
      states[ch]->cacheBoostNext = is_partial_cl_move(moveType);
    }

    // Chain swap: propose one random adjacent pair per iteration
    if (nChains > 1) {
      int iPair = static_cast<int>(R::unif_rand() * nSwapPairs);
      if (iPair >= nSwapPairs) iPair = nSwapPairs - 1;
      int jPair  = iPair + 1;
      swapPropose[iPair]++;
      double logAlpha = (betas[iPair] - betas[jPair]) *
                        (states[jPair]->logLik - states[iPair]->logLik);
      if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
        std::swap(*states[iPair], *states[jPair]);
        swapAccept[iPair]++;
        if (iPair == 0) coldSwapsSinceSample++;

        // PT-RT-001: update particle-at-slot map and possibly count round trip.
        std::swap(particleAtSlot[iPair], particleAtSlot[jPair]);
        int slotHot = nChains - 1;
        for (int s_chk : { iPair, jPair }) {
          int p = particleAtSlot[s_chk];
          if (s_chk == 0) {
            if (particleState[p] == PS_AT_HOT) roundTripCount++;
            particleState[p] = PS_AT_COLD;
          } else if (s_chk == slotHot) {
            if (particleState[p] == PS_AT_COLD) particleState[p] = PS_AT_HOT;
          }
        }
      }
    }

    // Save cold chain (index 0) sample post-warmup on thinning interval
    if (iter > warmup && (iter - warmup) % thin == 0) {
      McmcState* s0 = states[0];
      std::vector<double> row(nScalarCols);
      int col = 0;
      row[col++] = s0->logLik + s0->logPrior;   // log_post
      row[col++] = s0->logLik;
      row[col++] = s0->treeLength;
      if (hasNeo) row[col++] = s0->rateLoss;
      row[col++] = s0->rateLogSd;
      if (includeBG) {
        row[col++] = s0->kprimeAlpha;
        row[col++] = s0->kprimeBeta;
      } else if (includeP) {
        row[col++] = s0->p;
      }
      if (hasNeo) row[col++] = s0->rateNeo;
      if (includeBS) row[col++] = s0->betaScale;  // M-052
      // Diagnostic: cold-chain swaps since last sample
      row[col++] = static_cast<double>(coldSwapsSinceSample);
      coldSwapsSinceSample = 0;
      // Diagnostic: unrooted topology fingerprint
      row[col++] = fnv_topo_hash(s0->parent, s0->child, data->nTip);
      for (int j = 0; j < nTrans; ++j)
        row[col++] = static_cast<double>(s0->kPrime[transIdxCpp[j]]);
      for (int k = 0; k < nEdge; ++k)
        row[col++] = s0->relBrLengths[k];
      // Per-class partition columns (Layer 1)
      for (int k = 0; k < nClassRLS; ++k)
        row[col++] = s0->classRateLogSd[k];
      for (int k = 0; k < nClassW; ++k)
        row[col++] = s0->classW[k];
      // Hyperprior columns (hyper_tau and per-class z_c).
      if (nHyperTauCol > 0) row[col++] = s0->hyperTau;
      for (int k = 0; k < nClassZ; ++k)
        row[col++] = s0->classZ[k];
      scalarRows.push_back(row);

      // Edge matrix for tree reconstruction in R
      IntegerMatrix edgeMat(nEdge, 2);
      for (int k = 0; k < nEdge; ++k) {
        edgeMat(k, 0) = s0->parent[k];
        edgeMat(k, 1) = s0->child[k];
      }
      edgeSamples.push_back(edgeMat);
    }
  }

  // Pack scalar rows into R matrix
  int nSaved = static_cast<int>(scalarRows.size());
  NumericMatrix scalarMat(nSaved, nScalarCols);
  for (int i = 0; i < nSaved; ++i)
    for (int j = 0; j < nScalarCols; ++j)
      scalarMat(i, j) = scalarRows[i][j];

  // DIAG: aggregate counters from cold chain (index 0)
  IntegerVector diagCounters = IntegerVector::create(
    _["dir_partial"] = states[0]->diagDirPartialCount,
    _["dir_fullback"] = states[0]->diagDirFullbackCount,
    _["nni_partial"] = states[0]->diagNniPartialCount,
    _["bs_partial"] = states[0]->diagBsPartialCount
  );

  return List::create(
    _["accept_counts"]    = acceptCounts,
    _["propose_counts"]   = proposeCounts,
    _["move_time_ns"]     = moveTimeNs,
    _["slice_expansions"] = sliceExpansions,
    _["swap_accept"]      = swapAccept,
    _["swap_propose"]     = swapPropose,
    _["round_trip_count"] = roundTripCount,
    _["scalar_samples"]   = scalarMat,
    _["edge_samples"]     = edgeSamples,
    _["n_saved"]          = nSaved,
    _["diag_counters"]    = diagCounters,
    _["cache_hits"]       = cacheHits,
    _["cache_misses"]     = cacheMisses
  );
}

// [[Rcpp::export]]
List debug_mcmc_data(SEXP dataPtr) {
  McmcData* data = Rcpp::XPtr<McmcData>(dataPtr).get();
  return List::create(
    _["kPriorLogseries"]      = data->kPriorLogseries,
    _["kPriorBetaGeometric"]  = data->kPriorBetaGeometric,
    _["kprimeLogseriesC"]     = data->kprimeLogseriesC,
    _["kprimeHyperA"]      = data->kprimeHyperA,
    _["kprimeHyperB"]      = data->kprimeHyperB,
    _["treeLengthShape"]   = data->treeLengthShape,
    _["treeLengthRate"]    = data->treeLengthRate,
    _["transIdxGlobal"]    = data->transIdxGlobal,
    _["kObs"]              = data->kObs,
    _["nCat"]              = data->nCat
  );
}


// M-111: Validation — compare partial CL vs full evaluation for all swap
// candidates of a given nodeA.  Returns a data.frame with columns:
//   nodeA, nodeB, ll_partial, ll_full
// [[Rcpp::export]]
DataFrame validate_swap_partial_cl(SEXP dataPtr, SEXP statePtr, int nodeA) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();

  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  std::vector<int> partners = get_valid_swap_partners_impl(
      state->parent, state->child, nTip, nodeA);
  const int nPart = (int)partners.size();

  // --- Full evaluation (same as gibbs_subtree_swap_impl_full) ---
  const int rowA = find_child_row_gibbs(state->child, nodeA);
  IntegerVector workPar = clone(state->parent);
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);
  const int origParA = workPar[rowA];
  const double origAbsA = absLen[rowA];

  NumericVector llFull(nPart);
  for (int pi = 0; pi < nPart; ++pi) {
    int rowB = find_child_row_gibbs(state->child, partners[pi]);
    if (rowB < 0 || rowA == rowB) { llFull[pi] = R_NegInf; continue; }

    int origParB = workPar[rowB];
    double origAbsB = absLen[rowB];

    workPar[rowA] = origParB;  workPar[rowB] = origParA;
    absLen[rowA]  = origAbsB;  absLen[rowB]  = origAbsA;

    preorder_into(workPar, state->child, absLen, nTip,
                  INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
    llFull[pi] = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs);

    workPar[rowA] = origParA;  workPar[rowB] = origParB;
    absLen[rowA]  = origAbsA;  absLen[rowB]  = origAbsB;
  }

  // --- Partial CL evaluation ---
  // Through the same entry point the move uses, so this validates the kernel
  // rather than a copy of it.
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  std::vector<int> partialPartners;
  std::vector<double> partialLL;
  swap_neighbourhood_partial(data, state, state->parent, state->child,
                             absLen, nodeA, partialPartners, partialLL);

  NumericVector llPartial(nPart);
  for (int pi = 0; pi < nPart; ++pi)
    llPartial[pi] = (pi < (int)partialLL.size() &&
                     partialPartners[pi] == partners[pi])
                      ? partialLL[pi] : R_NaReal;

  IntegerVector nodeAVec(nPart, nodeA);
  IntegerVector nodeBVec(nPart);
  for (int pi = 0; pi < nPart; ++pi) nodeBVec[pi] = partners[pi];

  return DataFrame::create(
    _["nodeA"]      = nodeAVec,
    _["nodeB"]      = nodeBVec,
    _["ll_partial"] = llPartial,
    _["ll_full"]    = llFull
  );
}


