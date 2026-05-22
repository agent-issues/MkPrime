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
#include "gibbs_z_workspace.h"
#include "chain_rng.h"
#include "gibbs_partial_cl.h"
#include "fitch.h"
#include "node_cl_cache.h"
#include "ecology_cl_cache.h"
#include <TreeTools/renumber_tree.h>
#include <cmath>
#include <cstring>
#include <chrono>
#include <cstdio>

using namespace Rcpp;

// Forward declaration for relabelling correction (corrections.cpp)
double mk_prime_relabel_log(int kPrime, int kObs);

// Forward declaration for ecology orchestrator (mcmc_ecology.cpp).
// T-017: const-ref args (was pass-by-value) to avoid the Rcpp Vector copy-
// constructor (precious-object list mutation, not thread-safe).
double cpp_log_likelihood_ecology(
    const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0,
    const NumericVector& theta);

// Forward declarations for Phase 3f Gibbs sweep helpers (mcmc_ecology.cpp).
double per_char_log_lik_ecology(
    const McmcData& data,
    int globalCharIdx,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    int kp,
    double rateLoss, double rateNeo,
    const NumericVector& rates,
    NumericVector phi,
    IntegerVector zRow,
    const NumericMatrix& wEdge,
    int refEcology,
    const std::vector<double>& gammaE,
    GibbsZWorkspace& ws);

void recompute_w_edge(
    const McmcData& data,
    IntegerVector parent, IntegerVector child,
    NumericVector edgeLen,
    NumericMatrix& wEdgeOut);

// T-010: per-partition ecology likelihood (mcmc_ecology.cpp).
// Accepts pre-computed wEdge, gammaE, rates so the caller can hoist them
// out of a multi-partition or move loop and reuse them across calls.
//
// T-017: signature changed to const-ref for all Rcpp args.  Pass-by-value
// invokes the Rcpp Vector copy-constructor which touches the precious-object
// list (Rcpp_PreciousPreserve/Release) — that's not thread-safe.  Const-ref
// avoids the copy entirely.  No call-site change required; all existing
// callers pass named objects.
double cpp_partition_log_likelihood_ecology(
    const McmcData& data, int partIdx,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const IntegerVector& kPrime,
    double rateLoss, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE,
    const NumericVector& rates);

// T-010: build the kEco-vector gammaE from current state (cheap O(kEco)).
void compute_gamma_e_ecology(
    const McmcData& data,
    const NumericVector& phi,
    double pi0,
    const NumericVector& theta,
    std::vector<double>& gammaE);

// T-013: ecology partial-CL cache management (defined in mcmc_ecology.cpp).
void populate_eco_cache_full(
    EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateLogSd, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0, const NumericVector& theta,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE
);

double eco_cache_total_loglik(
    const EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen,
    const IntegerVector& kPrime, double rateLoss, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE
);

double eco_cache_partial_eval_nni(
    EcoCLCache& cache, const McmcData& data,
    const IntegerVector& parent, const IntegerVector& child,
    const NumericVector& edgeLen, const IntegerVector& kPrime,
    double rateLoss, double rateNeo,
    const NumericVector& phi,
    const IntegerMatrix& zMatrix,
    double pi0, const NumericVector& theta,
    const NumericMatrix& wEdge,
    const std::vector<double>& gammaE,
    int v, int u,
    EcoDirtyScratch& scratch,
    int& dirtyCount,
    bool& fallback,
    int* nWedgeDirtyEdges = nullptr,
    double dirtyFracThreshold = 0.70
);

// cpp_acrv_rates lives in mcmc_likelihood.cpp.
NumericVector cpp_acrv_rates(double rateLogSd, int nCat,
                              const std::vector<double>& acrvZ);

// Forward declarations for proposals in other TUs
// M-065: vector-based _impl versions (no edge matrix)
// T-017-IIa: all proposal impls take ChainRng& rng as first argument
List nni_proposal_impl(ChainRng& rng,
                       IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
List spr_proposal_impl(ChainRng& rng,
                       IntegerVector parent, IntegerVector child,
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
List tbr_proposal_impl(ChainRng& rng,
                       IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
List beta_simplex_proposal(NumericVector x, int index, double tuning);
bool beta_simplex_impl(ChainRng& rng,
                       NumericVector& x, int index, double tuning,
                       double& logHastings, int& outOther,
                       double& outOldIdx, double& outOldOther);  // OPP-5

// M-125: block Dirichlet simplex proposal (defined in proposals.cpp)
bool dirichlet_simplex_impl(ChainRng& rng,
                            NumericVector& x, int nCats, double alpha,
                            double& logHastings, NumericVector& snapshot,
                            std::vector<int>& modifiedEdges);
// M-127: localized Dirichlet simplex (neighborhood selection)
bool local_dirichlet_impl(ChainRng& rng,
                          NumericVector& x,
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

// T-017-IIa: ChainRng overload for in-MCMC callers (15 sites).
// Zero-arg version above is kept for the R-exported bactrian_draws wrapper.
static inline double bactrian_perturbation(ChainRng& rng) {
  double z = rng.rnorm(0.0, BACTRIAN_SD);
  double raw = (rng.unif() < 0.5) ? (BACTRIAN_M + z) : (-BACTRIAN_M + z);
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

// T-017-IIa: ChainRng overload for in-MCMC callers (4 sites).
// Zero-arg version above is kept for the R-exported bactrian_2d_draws wrapper.
static inline void bactrian_2d_perturbation(ChainRng& rng, double rho,
                                            double& z1, double& z2) {
  double n1 = rng.rnorm(0.0, 1.0);
  double n2 = rng.rnorm(0.0, 1.0);
  double e1 = BACTRIAN_SD * n1;
  double e2 = BACTRIAN_SD * (rho * n1 + std::sqrt(1.0 - rho * rho) * n2);
  double s1 = (rng.unif() < 0.5) ? BACTRIAN_M : -BACTRIAN_M;
  double pSame = 0.5 * (1.0 + rho);
  double s2 = (rng.unif() < pSame) ? s1 : -s1;
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

  // Ecology-aware NT model state. Only meaningful when
  // McmcData::ecologyAware is true; default values are inert.
  NumericVector phi;           // length 1 (global) or kEcology (per_ecology)
  double        pi0 = 0.0;
  // v2: theta has length (kEcology - 1); one entry per non-reference ecology
  // in ascending ecology-state order (i.e. ecology states < refEcology come
  // first, then ecology states > refEcology). zMatrix matches columnwise.
  NumericVector theta;
  IntegerMatrix zMatrix;       // nChar x (kEcology - 1), entries in {0, 1, 2}
  NumericMatrix wEdge;         // nEdge x kEcology, cached per-likelihood
  bool          wEdgeDirty = true;  // recompute wEdge on next likelihood call
  // M-155: dedicated workspace for Gibbs kPrime batched pruning
  ClWorkspace gibbsWs;
  // M-121: persistent node-level CL cache for partial evaluation
  NodeCLCache nodeCL;
  // T-011: persistent ecology CL cache (parallel to nodeCL for the eco path).
  // Populated lazily; consulted only when ecologyAware mode is active.
  EcoCLCache ecoCL;
  // M-125: snapshot for block Dirichlet branch-length rollback
  NumericVector brSnapshot;
  // M-127: which edges the Dirichlet proposal modified (for partial CL eval)
  std::vector<int> dirEdges;
  // DIAG counters
  int diagDirPartialCount = 0;
  int diagDirMismatchCount = 0;
  int diagDirFullbackCount = 0;
  int diagNniPartialCount = 0;
  int diagNniMismatchCount = 0;
  int diagBsPartialCount = 0;
  int diagBsMismatchCount = 0;
  int diagDriftCount = 0;
  int diagCachePopCount = 0;
  double diagMaxDiff = 0.0;

};


// ---------------------------------------------------------------------------
// T-010: helpers for the ecology-aware partition cache
// ---------------------------------------------------------------------------
//
// `state->wEdge` is the per-edge ecology-state weight matrix (nEdge × kEco),
// expensive to recompute (full forward+backward sweep on the tree).
// `state->wEdgeDirty` is set true whenever any move modifies the tree topology
// or any branch length; cleared after a successful recompute. gammaE and ACRV
// rates are cheap and computed per call.
//
// `state->partLogLik` is reused as the ecology partition cache (the blind path
// uses it likewise). When ecologyAware is true the cache is initialised from
// the ecology orchestrator (see fill_partition_cache); when set the
// per-partition sum equals state->logLik.

// Recompute and cache state->wEdge if dirty. No-op otherwise.
static inline void eco_refresh_wedge(McmcData* data, McmcState* state,
                                     const NumericVector& edgeLen) {
  if (!state->wEdgeDirty &&
      state->wEdge.nrow() == edgeLen.size() &&
      state->wEdge.ncol() == data->ecology.kEcology) {
    return;
  }
  if (state->wEdge.nrow() != edgeLen.size() ||
      state->wEdge.ncol() != data->ecology.kEcology) {
    state->wEdge = NumericMatrix(edgeLen.size(), data->ecology.kEcology);
  }
  recompute_w_edge(*data, state->parent, state->child, edgeLen, state->wEdge);
  state->wEdgeDirty = false;
}

// Compute one partition's ecology log-lik, using the cached wEdge and a
// freshly-computed gammaE / rates. Convenience wrapper around the per-
// partition function for move-handler call sites.
static double eco_partition_loglik(McmcData* data, McmcState* state,
                                   int partIdx, const NumericVector& edgeLen) {
  eco_refresh_wedge(data, state, edgeLen);
  std::vector<double> gammaE;
  compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta, gammaE);
  NumericVector rates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);
  return cpp_partition_log_likelihood_ecology(
    *data, partIdx, state->parent, state->child, edgeLen,
    state->kPrime, state->rateLoss, state->rateNeo,
    state->phi, state->zMatrix,
    state->wEdge, gammaE, rates);
}

// Recompute ALL ecology partitions, refresh state->partLogLik and
// state->logLik. Cheaper than cpp_log_likelihood_ecology by ~tree-build cost
// when wEdge is already cached (e.g. phi/pi0/theta/rateLogSd moves).
static double eco_recompute_all_partitions(McmcData* data, McmcState* state,
                                           const NumericVector& edgeLen) {
  eco_refresh_wedge(data, state, edgeLen);
  std::vector<double> gammaE;
  compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta, gammaE);
  NumericVector rates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);
  int nParts = (int)data->parts.size();
  if ((int)state->partLogLik.size() != nParts) state->partLogLik.assign(nParts, 0.0);
  double total = 0.0;
  for (int pi = 0; pi < nParts; ++pi) {
    state->partLogLik[pi] = cpp_partition_log_likelihood_ecology(
      *data, pi, state->parent, state->child, edgeLen,
      state->kPrime, state->rateLoss, state->rateNeo,
      state->phi, state->zMatrix,
      state->wEdge, gammaE, rates);
    total += state->partLogLik[pi];
  }
  return total;
}


// ---------------------------------------------------------------------------
// Log prior (mirrors LogPrior in MkPrimeModel.R)
// ---------------------------------------------------------------------------

static double cpp_log_prior(
    const McmcData& data,
    double treeLength, const NumericVector& relBrLengths,
    double rateLoss, double rateLogSd, double rateNeo,
    double p, const IntegerVector& kPrime,
    double betaScale = 1.0,
    double kprimeAlpha = 1.0, double kprimeBeta = 1.0,
    const NumericVector* phiPtr   = nullptr,
    double pi0                    = 0.0,
    const IntegerMatrix* zMatPtr  = nullptr,
    const NumericVector* thetaPtr = nullptr) {

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
    for (int i = 0; i < data.transIdxGlobal.size(); ++i) {
      int gi = data.transIdxGlobal[i];
      if (kPrime[gi] < data.kObs[gi]) return R_NegInf;
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
      // beyond data.empBodyLastK the pmf decays geometrically with
      // log-mass `empLogTailStartP + (k - empTailStartK) * log(empTailDecay)`.
      double logP    = std::log(p);
      double log1mP  = std::log1p(-p);
      double logQ    = (data.empTailDecay > 0.0) ? std::log(data.empTailDecay)
                                                 : R_NegInf;
      int    bodyLen = static_cast<int>(data.empLogBody.size());
      // Scratch buffer for logSumExp (reused per character)
      std::vector<double> terms;
      terms.reserve(64);
      for (int i = 0; i < nTrans; ++i) {
        int gi = data.transIdxGlobal[i];
        int m = kPrime[gi];
        if (m < 2) return R_NegInf;
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
        lp += mx + std::log(sumExp);
      }
      // p: Beta hyperprior (same as plain geometric)
      lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
    } else {
      // Hierarchical geometric: P(k'_i = kObs_i + u) = p*(1-p)^u
      double sumU = 0.0;
      for (int i = 0; i < nTrans; ++i) {
        int gi = data.transIdxGlobal[i];
        sumU += (kPrime[gi] - data.kObs[gi]);
      }
      lp += nTrans * std::log(p) + sumU * std::log1p(-p);
      lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
    }
  }

  // M-052: beta_scale prior — Gamma(shape, rate)
  if (data.qHeterogeneity) {
    if (betaScale <= 0.0) return R_NegInf;
    lp += R::dgamma(betaScale, data.betaScaleShape,
                     1.0 / data.betaScaleRate, 1);
  }

  // Ecology-aware NT model (v2): phi, pi0, theta, z priors.
  // Mirrors LogPrior in R/MkPrimeModel.R. Asymmetric slab:
  //   P(z = none) = pi0
  //   P(z = enc)  = (1 - pi0) * theta_e
  //   P(z = disc) = (1 - pi0) * (1 - theta_e)
  // theta_e ~ Beta(thetaAlpha, thetaBeta).
  if (data.ecologyAware) {
    if (phiPtr == nullptr || zMatPtr == nullptr || thetaPtr == nullptr)
      return R_NegInf;
    const NumericVector& phi = *phiPtr;
    const IntegerMatrix& zMat = *zMatPtr;
    const NumericVector& theta = *thetaPtr;
    for (int i = 0; i < phi.size(); ++i) {
      if (phi[i] <= 0.0) return R_NegInf;
    }
    if (pi0 <= 0.0 || pi0 >= 1.0) return R_NegInf;
    int nCharZ = zMat.nrow();
    int kEcoZ  = zMat.ncol();
    if (theta.size() != kEcoZ) return R_NegInf;
    for (int e = 0; e < kEcoZ; ++e) {
      if (theta[e] <= 0.0 || theta[e] >= 1.0) return R_NegInf;
    }
    for (int i = 0; i < phi.size(); ++i) {
      lp += R::dlnorm(phi[i], 0.0, data.sigmaPhi, 1);
    }
    lp += R::dbeta(pi0, data.rho0Alpha, data.rho0Beta, 1);

    double logPi0     = std::log(pi0);
    double log1mPi0   = std::log1p(-pi0);
    for (int e = 0; e < kEcoZ; ++e) {
      long nNone = 0, nEnc = 0, nDisc = 0;
      for (int c = 0; c < nCharZ; ++c) {
        int zv = zMat(c, e);
        if      (zv == 0) ++nNone;
        else if (zv == 1) ++nEnc;
        else if (zv == 2) ++nDisc;
        else return R_NegInf;
      }
      double logTheta   = std::log(theta[e]);
      double log1mTheta = std::log1p(-theta[e]);
      if (nNone > 0) lp += static_cast<double>(nNone) * logPi0;
      if (nEnc  > 0) lp += static_cast<double>(nEnc)  * (log1mPi0 + logTheta);
      if (nDisc > 0) lp += static_cast<double>(nDisc) * (log1mPi0 + log1mTheta);
      lp += R::dbeta(theta[e], data.thetaAlpha, data.thetaBeta, 1);
    }
  }

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
                     Rcpp::Nullable<Rcpp::NumericVector> phi = R_NilValue,
                     double pi0 = 0.0,
                     Rcpp::Nullable<Rcpp::IntegerMatrix> zMatrix = R_NilValue,
                     Rcpp::Nullable<Rcpp::NumericVector> theta = R_NilValue) {
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
  if (phi.isNotNull()) {
    Rcpp::NumericVector phiVec(phi);
    if (phiVec.size() > 0) {
      s->phi      = clone(phiVec);
      s->pi0      = pi0;
      if (zMatrix.isNotNull()) {
        s->zMatrix = clone(Rcpp::IntegerMatrix(zMatrix));
      }
      if (theta.isNotNull()) {
        s->theta = clone(Rcpp::NumericVector(theta));
      }
      s->wEdgeDirty = true;
    }
  }
  return Rcpp::XPtr<McmcState>(s, true);
}


// [[Rcpp::export]]
void fill_partition_cache(SEXP dataPtr, SEXP statePtr) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();
  int nParts = (int)data->parts.size();
  int nEdge  = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  if (data->ecologyAware) {
    // T-010: populate the ecology partition cache. wEdge is recomputed and
    // cached on state->wEdge so subsequent non-tree moves can skip the
    // forward+backward sweep over the ecology tree.
    state->wEdgeDirty = true;
    state->logLik = eco_recompute_all_partitions(data, state, edgeLen);
    return;
  }

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
    _["diagDirMismatch"] = s->diagDirMismatchCount,
    _["diagDirFullback"] = s->diagDirFullbackCount,
    _["diagNniPartial"] = s->diagNniPartialCount,
    _["diagBsPartial"]  = s->diagBsPartialCount,
    _["diagDriftCount"] = s->diagDriftCount,
    _["diagMaxDiff"]    = s->diagMaxDiff,
    _["diagSelectivePop"] = s->nodeCL.diagSelectivePopCount,
    // Clone mutable ecology vectors so that R-side snapshots taken before a
    // move don't alias the underlying storage and silently reflect post-move
    // mutations.
    _["phi"]            = clone(s->phi),
    _["pi0"]            = s->pi0,
    _["theta"]          = clone(s->theta),
    _["zMatrix"]        = clone(s->zMatrix)
  );
}


// FNV-1a topology hash of the parent vector.
// Edges must be in canonical preorder (guaranteed by all tree moves).
static double fnv_topo_hash(const IntegerVector& parent) {
  uint64_t h = 0xcbf29ce484222325ULL;
  for (int k = 0; k < parent.size(); ++k) {
    h ^= static_cast<uint64_t>(parent[k]);
    h *= 0x100000001b3ULL;
  }
  return static_cast<double>(h >> 11);  // 53 significant bits
}

// [[Rcpp::export]]
double compute_topo_hash(IntegerVector parent) {
  return fnv_topo_hash(parent);
}

// [[Rcpp::export]]
double get_state_log_lik(SEXP statePtr) {
  return Rcpp::XPtr<McmcState>(statePtr).get()->logLik;
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
    s->kprimeAlpha, s->kprimeBeta,
    &s->phi, s->pi0, &s->zMatrix, &s->theta);
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
    const NumericVector& edgeLen) {
  if (data.ecologyAware) {
    return cpp_log_likelihood_ecology(
      data, parent, child, edgeLen, state.kPrime,
      state.rateLoss, state.rateLogSd, state.rateNeo,
      state.phi, state.zMatrix, state.pi0, state.theta);
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


// ---------------------------------------------------------------------------
// gibbs_spr_impl  (M-085, M-105 partial CL)
//
// GibbsSPR with partial CL reuse: enumerate all valid SPR reattachment
// positions for a randomly chosen subtree.  Instead of running a full-tree
// pruning per candidate, cache per-node CLs from one downpass and update
// only the O(depth) affected path per candidate.
//
// Falls back to the old full-evaluation path when Q-heterogeneity is
// enabled (M-052), since that changes the pruning model.
// ---------------------------------------------------------------------------

// Old full-evaluation path (used as fallback and for validation)
static bool gibbs_spr_impl_full(ChainRng& rng, McmcData* data, McmcState* state, double beta);
// M-114: partial CL path for Q-heterogeneity
static bool gibbs_spr_impl_het(ChainRng& rng, McmcData* data, McmcState* state, double beta);

static bool gibbs_spr_impl(ChainRng& rng, McmcData* data, McmcState* state, double beta) {
  // M-114: Q-heterogeneity uses streaming partial CL
  if (data->qHeterogeneity)
    return gibbs_spr_impl_het(rng, data, state, beta);

  // Ecology mode: the partial-CL machinery used to compute candLL below
  // does not include the per-character ecology rate modifiers (phi, z,
  // pi0, gamma_e). Letting this move fire would silently overwrite
  // state->logLik with a non-ecology value and bias subsequent MH
  // ratios. Skip until an eco-aware streaming candidate evaluator
  // exists. The chain falls back on spr/nni/tbr/pspr for topology
  // exploration in eco mode.
  if (data->ecologyAware) return false;

  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  // 1. Eligible prune edges (parent != root)
  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  // 2. Pick random prune edge
  int pickIdx = (int)(rng.unif() * (double)eligible.size());
  if (pickIdx >= (int)eligible.size()) pickIdx = (int)eligible.size() - 1;
  const int pruneRow = eligible[pickIdx];
  const int u = state->parent[pruneRow];
  const int v = state->child[pruneRow];

  // 3. Find parentRow, sibRow, sibNode
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

  // 5. Collect valid regraft candidate edges
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
  const double lPrune = absLen[pruneRow];

  // ===== M-105: Partial CL cache setup =====

  // Build tree navigation from current topology
  TreeNav topo;
  topo.build(state->parent, state->child, absLen, nTip);

  // Compute ACRV rates
  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding = data->codingType;
  int maxNode = topo.maxNode;

  // Build CLGroups: one per (partition, kStates) evaluation unit
  std::vector<CLGroup> groups;
  // Track which groups belong to which partition (for ascertainment)
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];

    if (part.type == 0) {
      // Neomorphic: k=2, MkN model
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = state->rateNeo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

    } else if (part.type == 2) {
      // Known state space: fixed k
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = 1.0;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

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
        g.rateScale = 1.0;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
        groupMeta.push_back({pi, nCharPart});
      }
    }
  }

  // Run caching downpass for each group
  for (auto& grp : groups)
    caching_downpass(grp, topo, state->parent, state->child, rates);

  // Compute residual CLs for each group (detach v from u)
  std::vector<ResidualCL> residuals(groups.size());
  for (size_t gi = 0; gi < groups.size(); ++gi)
    compute_residual_cl(residuals[gi], groups[gi], topo, rates, u, sibNode, lMerge);

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
                          u, sibNode, lMerge);
    }
  }

  // ===== Evaluate candidates using partial CLs =====

  std::vector<double> candLL(nCand);
  for (int ci = 0; ci < nCand; ++ci) {
    const int rr = cands[ci];
    const int a  = state->parent[rr];
    const int b  = state->child[rr];
    const double lReg = absLen[rr];
    const double lHalf = 0.5 * lReg;

    double totalLL = 0.0;

    for (size_t gi = 0; gi < groups.size(); ++gi) {
      double grpLL = evaluate_candidate(
        groups[gi], topo, residuals[gi], rates,
        v, u, sibNode, lMerge, a, b, lHalf, lPrune);

      // Ascertainment correction via pseudo-character partial CLs
      if (coding != 0 && groups[gi].nChar > 0) {
        double constP = evaluate_const_prob(
          pseudoGroups[gi], topo, pseudoResiduals[gi], rates,
          v, u, sibNode, lMerge, a, b, lHalf, lPrune);
        // TODO: coding == 2 (informative) needs singleton_site_prob too
        if (constP < 1.0)
          grpLL -= groups[gi].nChar * std::log(1.0 - constP);
      }

      totalLL += grpLL;
    }

    candLL[ci] = totalLL;
  }

  // Add relabelling correction to candLL — it's a topology-independent
  // constant that's included in state->logLik but not in evaluate_candidate.
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
    for (int ci = 0; ci < nCand; ++ci)
      candLL[ci] += relabelCorr;
  }

  // 8. Sampling weights: exp(β × logLik), current state included
  const double llOrig = state->logLik;
  double maxLL = llOrig;
  for (int ci = 0; ci < nCand; ++ci) maxLL = std::max(maxLL, candLL[ci]);

  double wOrig = std::exp(beta * (llOrig - maxLL));
  std::vector<double> ws(nCand);
  double sumW = wOrig;
  for (int ci = 0; ci < nCand; ++ci) {
    ws[ci] = std::exp(beta * (candLL[ci] - maxLL));
    sumW  += ws[ci];
  }

  // 9. Sample
  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < ws[ci]) { chosen = ci; break; }
    rnd -= ws[ci];
  }

  // 10. Apply chosen SPR: modify state in-place, canonical reorder
  {
    const int rr      = cands[chosen];
    const double lReg = absLen[rr];
    const int b = state->child[rr];

    state->child[parentRow]  = sibNode;  absLen[parentRow] = lMerge;
    state->child[rr]         = u;        absLen[rr]        = 0.5 * lReg;
    state->parent[sibRow]    = u;        state->child[sibRow] = b;
    absLen[sibRow]           = 0.5 * lReg;

    auto po = TreeTools::preorder_weighted_impl(
        state->parent, state->child, absLen);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbs  = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbs[k] / state->treeLength;
    }
  }

  // 11. Commit
  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// ---------------------------------------------------------------------------
// M-114: Gibbs SPR with streaming partial CL for Q-heterogeneity.
//
// Same prune/candidate/sampling logic as gibbs_spr_impl, but evaluates
// each CLGroup's likelihood under a mixture of F81 components by streaming
// over (betaBin, rotation) and accumulating per-site raw likelihoods.
// ---------------------------------------------------------------------------
static bool gibbs_spr_impl_het(ChainRng& rng, McmcData* data, McmcState* state,
                                double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  // 1. Eligible prune edges (parent != root)
  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  int pickIdx = (int)(rng.unif() * (double)eligible.size());
  if (pickIdx >= (int)eligible.size()) pickIdx = (int)eligible.size() - 1;
  const int pruneRow = eligible[pickIdx];
  const int u = state->parent[pruneRow];
  const int v = state->child[pruneRow];

  int parentRow = -1, sibRow = -1, sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == u) parentRow = i;
    if (state->parent[i] == u && state->child[i] != v) {
      sibRow = i; sibNode = state->child[i];
    }
  }
  if (parentRow < 0 || sibRow < 0) return false;

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

  std::vector<int> cands;
  cands.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[state->child[i]]) continue;
    if (state->parent[i] == u || state->child[i] == u) continue;
    cands.push_back(i);
  }
  if (cands.empty()) return false;
  const int nCand = (int)cands.size();

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];
  const double lMerge = absLen[parentRow] + absLen[sibRow];
  const double lPrune = absLen[pruneRow];

  TreeNav topo;
  topo.build(state->parent, state->child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat,
                                          data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding  = data->codingType;
  int maxNode = topo.maxNode;
  int nBC     = data->nBetaCat;

  // Build CLGroups (same structure as non-het, but with useF81 flag)
  std::vector<CLGroup> groups;
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = state->rateNeo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = 1.0;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
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
        g.rateScale = 1.0;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
        groupMeta.push_back({pi, nCharPart});
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

  // ===== M-114: Streaming evaluation over (betaBin, rotation) =====
  //
  // Per-group, per-candidate accumulators for raw site likelihoods
  // siteLikAccums[gi][ci * nChar_gi .. (ci+1)*nChar_gi - 1]
  std::vector<std::vector<double>> siteLikAccums(groups.size());
  std::vector<std::vector<double>> constProbAccums(groups.size());
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    siteLikAccums[gi].assign((size_t)nCand * groups[gi].nChar, 0.0);
    if (coding != 0)
      constProbAccums[gi].assign(
        (size_t)nCand * pseudoGroups[gi].nChar, 0.0);
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
        compute_residual_cl(res, grp, topo, rates, u, sibNode, lMerge);

        // Same for pseudo-group (ascertainment)
        if (coding != 0) {
          set_f81_component(pseudoGroups[gi], hetBins[bi], rot, baseRL);
          caching_downpass(pseudoGroups[gi], topo, state->parent,
                           state->child, rates);
          compute_residual_cl(pseudoRes, pseudoGroups[gi], topo, rates,
                              u, sibNode, lMerge);
        }

        // Evaluate all candidates for this component
        for (int ci = 0; ci < nCand; ++ci) {
          const int rr = cands[ci];
          const int a  = state->parent[rr];
          const int b  = state->child[rr];
          const double lHalf = 0.5 * absLen[rr];

          evaluate_candidate(
            grp, topo, res, rates,
            v, u, sibNode, lMerge, a, b, lHalf, lPrune,
            siteLikAccums[gi].data() + (size_t)ci * grp.nChar);

          if (coding != 0) {
            evaluate_const_prob(
              pseudoGroups[gi], topo, pseudoRes, rates,
              v, u, sibNode, lMerge, a, b, lHalf, lPrune,
              constProbAccums[gi].data() +
                (size_t)ci * pseudoGroups[gi].nChar);
          }
        }
      }
    }
  }

  // Convert accumulators to per-candidate log-likelihoods
  std::vector<double> candLL(nCand, 0.0);
  for (size_t gi = 0; gi < groups.size(); ++gi) {
    const CLGroup& grp = groups[gi];
    int k     = grp.kStates;
    int nRot  = (k == 2) ? 1 : k;
    int totalComp = nCat * nBC * nRot;
    int nChar_gi  = grp.nChar;

    for (int ci = 0; ci < nCand; ++ci) {
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
    for (int ci = 0; ci < nCand; ++ci)
      candLL[ci] += relabelCorr;
  }

  // Sampling and commit (identical to gibbs_spr_impl)
  const double llOrig = state->logLik;
  double maxLL = llOrig;
  for (int ci = 0; ci < nCand; ++ci)
    maxLL = std::max(maxLL, candLL[ci]);

  double wOrig = std::exp(beta * (llOrig - maxLL));
  std::vector<double> ws(nCand);
  double sumW = wOrig;
  for (int ci = 0; ci < nCand; ++ci) {
    ws[ci] = std::exp(beta * (candLL[ci] - maxLL));
    sumW  += ws[ci];
  }

  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < ws[ci]) { chosen = ci; break; }
    rnd -= ws[ci];
  }

  // Apply chosen SPR
  {
    const int rr      = cands[chosen];
    const double lReg = absLen[rr];
    const int b = state->child[rr];

    state->child[parentRow]  = sibNode;  absLen[parentRow] = lMerge;
    state->child[rr]         = u;        absLen[rr]        = 0.5 * lReg;
    state->parent[sibRow]    = u;        state->child[sibRow] = b;
    absLen[sibRow]           = 0.5 * lReg;

    auto po = TreeTools::preorder_weighted_impl(
        state->parent, state->child, absLen);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbs  = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbs[k] / state->treeLength;
    }
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// Old full-evaluation fallback (Q-het or validation), M-109 in-place
static bool gibbs_spr_impl_full(ChainRng& rng, McmcData* data, McmcState* state, double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  int pickIdx = (int)(rng.unif() * (double)eligible.size());
  if (pickIdx >= (int)eligible.size()) pickIdx = (int)eligible.size() - 1;
  const int pruneRow = eligible[pickIdx];
  const int u = state->parent[pruneRow];
  const int v = state->child[pruneRow];

  int parentRow = -1, sibRow = -1, sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == u) parentRow = i;
    if (state->parent[i] == u && state->child[i] != v) {
      sibRow = i; sibNode = state->child[i];
    }
  }
  if (parentRow < 0 || sibRow < 0) return false;

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

  std::vector<int> cands;
  cands.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    if (isDesc[state->child[i]]) continue;
    if (state->parent[i] == u || state->child[i] == u) continue;
    cands.push_back(i);
  }
  if (cands.empty()) return false;
  const int nCand = (int)cands.size();

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];
  const double lMerge = absLen[parentRow] + absLen[sibRow];

  // Working copies: cloned ONCE, reused for all candidates
  IntegerVector workPar = clone(state->parent);
  IntegerVector workCh  = clone(state->child);
  // Save original values for the 3 modified rows
  const int origPar_sibRow = workPar[sibRow];
  const int origCh_parentRow = workCh[parentRow];
  const int origCh_sibRow = workCh[sibRow];
  const double origAbs_parentRow = absLen[parentRow];
  const double origAbs_sibRow = absLen[sibRow];

  // Pre-allocate output buffers for preorder_into
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  std::vector<double> candLL(nCand);
  for (int ci = 0; ci < nCand; ++ci) {
    const int rr      = cands[ci];
    const double lReg = absLen[rr];
    const int bNode   = workCh[rr];  // original child of regraft edge

    // Apply SPR in-place
    workCh[parentRow] = sibNode;   absLen[parentRow] = lMerge;
    workCh[rr]        = u;         absLen[rr]        = 0.5 * lReg;
    workPar[sibRow]   = u;         workCh[sibRow]    = bNode;
    absLen[sibRow]    = 0.5 * lReg;

    preorder_into(workPar, workCh, absLen, nTip,
                  INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
    candLL[ci] = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs);

    // Restore
    workCh[parentRow] = origCh_parentRow;  absLen[parentRow] = origAbs_parentRow;
    workCh[rr]        = bNode;             absLen[rr]        = lReg;
    workPar[sibRow]   = origPar_sibRow;    workCh[sibRow]    = origCh_sibRow;
    absLen[sibRow]    = origAbs_sibRow;
  }

  const double llOrig = state->logLik;
  double maxLL = llOrig;
  for (int ci = 0; ci < nCand; ++ci) maxLL = std::max(maxLL, candLL[ci]);
  double wOrig = std::exp(beta * (llOrig - maxLL));
  std::vector<double> ws(nCand);
  double sumW = wOrig;
  for (int ci = 0; ci < nCand; ++ci) {
    ws[ci] = std::exp(beta * (candLL[ci] - maxLL));
    sumW  += ws[ci];
  }

  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < ws[ci]) { chosen = ci; break; }
    rnd -= ws[ci];
  }

  // Apply chosen SPR: modify state in-place, canonical reorder
  {
    const int rr      = cands[chosen];
    const double lReg = absLen[rr];  // absLen already restored to original
    const int bNode   = state->child[rr];

    state->child[parentRow]  = sibNode;  absLen[parentRow] = lMerge;
    state->child[rr]         = u;        absLen[rr]        = 0.5 * lReg;
    state->parent[sibRow]    = u;        state->child[sibRow] = bNode;
    absLen[sibRow]           = 0.5 * lReg;

    auto po = TreeTools::preorder_weighted_impl(
        state->parent, state->child, absLen);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbsFinal = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbsFinal[k] / state->treeLength;
    }
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// ---------------------------------------------------------------------------
// gibbs_subtree_swap_impl  (M-086, M-109 in-place, M-111 partial CL)
//
// GibbsSubtreeSwap: enumerate all valid subtree-swap partners for a randomly
// chosen node, weight by exp(β × logLik), sample proportionally, apply.
// Same Gibbs semantics and design choices as gibbs_spr_impl.
// Branch lengths swap with their subtrees (Jacobian = 1); see M-084.
//
// M-111: Partial CL reuse — cache per-node CLs from one downpass, then
// evaluate each candidate by updating only the O(depth) affected path
// (union of paths from pA and pB to root).  Falls back to full evaluation
// when Q-heterogeneity is enabled.
// ---------------------------------------------------------------------------

// Local helper: find edge row where child[i] == node
static int find_child_row_gibbs(const IntegerVector& child, int node) {
  for (int i = 0; i < child.size(); ++i)
    if (child[i] == node) return i;
  return -1;
}

// Full-evaluation fallback (Q-het or validation)
static bool gibbs_subtree_swap_impl_full(ChainRng& rng, McmcData* data, McmcState* state,
                                         double beta);
// M-114: partial CL path for Q-heterogeneity
static bool gibbs_subtree_swap_impl_het(ChainRng& rng, McmcData* data, McmcState* state,
                                         double beta);

static bool gibbs_subtree_swap_impl(ChainRng& rng, McmcData* data, McmcState* state,
                                    double beta) {
  // M-114: Q-heterogeneity uses streaming partial CL
  if (data->qHeterogeneity)
    return gibbs_subtree_swap_impl_het(rng, data, state, beta);

  // Ecology mode: same reason as gibbs_spr_impl — the partial-CL
  // candidate evaluation does not include phi/z/pi0/gamma_e, so the
  // move would silently corrupt state->logLik. Skip in eco mode.
  if (data->ecologyAware) return false;

  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  // 1. Pick a random node (any edge child is a valid non-root candidate)
  int pickIdx = (int)(rng.unif() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];

  // 2. Get valid swap partners
  std::vector<int> partners = get_valid_swap_partners_impl(
      state->parent, state->child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  // 3. Build absolute edge lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // ===== M-111: Partial CL cache setup =====

  TreeNav topo;
  topo.build(state->parent, state->child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding = data->codingType;
  int maxNode = topo.maxNode;

  // Build CLGroups: one per (partition, kStates) evaluation unit
  std::vector<CLGroup> groups;
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];

    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = state->rateNeo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = 1.0;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

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
        g.rateScale = 1.0;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
        groupMeta.push_back({pi, nCharPart});
      }
    }
  }

  // Run caching downpass for each group
  for (auto& grp : groups)
    caching_downpass(grp, topo, state->parent, state->child, rates);

  // Ascertainment correction: pseudo-character groups
  std::vector<CLGroup> pseudoGroups;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, maxNode, nCat);
      caching_downpass(pseudoGroups[gi], topo, state->parent, state->child,
                       rates);
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

  // ===== Evaluate candidates using partial CLs =====

  std::vector<double> candLL(nPart);
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

  // 4. Sampling weights: exp(β × logLik), current state included
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

  // 5. Sample: self-draw → no-op
  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nPart - 1;
  for (int pi = 0; pi < nPart - 1; ++pi) {
    if (rnd < ws[pi]) { chosen = pi; break; }
    rnd -= ws[pi];
  }
  if (!R_FINITE(candLL[chosen])) return false;

  // 6. Apply chosen swap: modify in-place then canonical reorder for state
  {
    const int rowA = find_child_row_gibbs(state->child, nodeA);
    const int rowB = find_child_row_gibbs(state->child, partners[chosen]);
    int swapParA = state->parent[rowB];
    int swapParB = state->parent[rowA];
    state->parent[rowA] = swapParA;
    state->parent[rowB] = swapParB;
    absLen[rowA] = state->treeLength * state->relBrLengths[rowB];
    absLen[rowB] = state->treeLength * state->relBrLengths[rowA];
    double tmpRel = state->relBrLengths[rowA];
    state->relBrLengths[rowA] = state->relBrLengths[rowB];
    state->relBrLengths[rowB] = tmpRel;

    // Canonical reorder (needed for NNI in-place invariant)
    auto po = TreeTools::preorder_weighted_impl(
        state->parent, state->child, absLen);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbsFinal = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbsFinal[k] / state->treeLength;
    }
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// ---------------------------------------------------------------------------
// M-114: Gibbs subtree swap with streaming partial CL for Q-heterogeneity.
// ---------------------------------------------------------------------------
static bool gibbs_subtree_swap_impl_het(ChainRng& rng, McmcData* data, McmcState* state,
                                         double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  int pickIdx = (int)(rng.unif() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];

  std::vector<int> partners = get_valid_swap_partners_impl(
      state->parent, state->child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  TreeNav topo;
  topo.build(state->parent, state->child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat,
                                          data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding  = data->codingType;
  int maxNode = topo.maxNode;
  int nBC     = data->nBetaCat;

  // Build CLGroups (same as gibbs_spr_impl_het)
  std::vector<CLGroup> groups;
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN     = true;
      g.rateLoss  = state->rateLoss;
      g.rateScale = state->rateNeo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = 1.0;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
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
        g.rateScale = 1.0;
        g.tipData   = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
        groupMeta.push_back({pi, nCharPart});
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
        caching_downpass(grp, topo, state->parent, state->child, rates);

        if (coding != 0) {
          set_f81_component(pseudoGroups[gi], hetBins[bi], rot, baseRL);
          caching_downpass(pseudoGroups[gi], topo, state->parent,
                           state->child, rates);
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
  std::vector<double> candLL(nPart, 0.0);
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

  // Sampling (identical to gibbs_subtree_swap_impl)
  const double llOrig = state->logLik;
  double maxLL = llOrig;
  for (int pi2 = 0; pi2 < nPart; ++pi2)
    maxLL = std::max(maxLL, candLL[pi2]);

  double wOrig = std::exp(beta * (llOrig - maxLL));
  std::vector<double> ws(nPart);
  double sumW = wOrig;
  for (int pi2 = 0; pi2 < nPart; ++pi2) {
    ws[pi2] = std::exp(beta * (candLL[pi2] - maxLL));
    sumW += ws[pi2];
  }

  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nPart - 1;
  for (int pi2 = 0; pi2 < nPart - 1; ++pi2) {
    if (rnd < ws[pi2]) { chosen = pi2; break; }
    rnd -= ws[pi2];
  }

  // Apply swap (same as gibbs_subtree_swap_impl)
  int nodeB = partners[chosen];
  int rowA = -1, rowB = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == nodeA) rowA = i;
    if (state->child[i] == nodeB) rowB = i;
  }
  state->child[rowA] = nodeB;
  state->child[rowB] = nodeA;
  absLen[rowA] = topo.edgeLen[topo.edgeToPar[nodeB]];
  absLen[rowB] = topo.edgeLen[topo.edgeToPar[nodeA]];

  auto po = TreeTools::preorder_weighted_impl(
      state->parent, state->child, absLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  for (int k2 = 0; k2 < nEdge; ++k2) {
    state->parent[k2]       = ordEdge(k2, 0);
    state->child[k2]        = ordEdge(k2, 1);
    state->relBrLengths[k2] = ordAbs[k2] / state->treeLength;
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// Full-evaluation fallback for Q-heterogeneity (M-109 in-place pattern)
static bool gibbs_subtree_swap_impl_full(ChainRng& rng, McmcData* data, McmcState* state,
                                         double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  int pickIdx = (int)(rng.unif() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];

  std::vector<int> partners = get_valid_swap_partners_impl(
      state->parent, state->child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  const int rowA = find_child_row_gibbs(state->child, nodeA);
  if (rowA < 0) return false;

  IntegerVector workPar = clone(state->parent);

  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  const int origParA = workPar[rowA];
  const double origAbsA = absLen[rowA];

  std::vector<double> candLL(nPart);
  for (int pi = 0; pi < nPart; ++pi) {
    int rowB = find_child_row_gibbs(state->child, partners[pi]);
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

    preorder_into(workPar, state->child, absLen, nTip,
                  INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
    candLL[pi] = compute_full_loglik_at(*data, *state, ordPar, ordCh, ordAbs);

    workPar[rowA] = origParA;
    workPar[rowB] = origParB;
    absLen[rowA]  = origAbsA;
    absLen[rowB]  = origAbsB;
  }

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

  double rnd = rng.unif() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nPart - 1;
  for (int pi = 0; pi < nPart - 1; ++pi) {
    if (rnd < ws[pi]) { chosen = pi; break; }
    rnd -= ws[pi];
  }
  if (!R_FINITE(candLL[chosen])) return false;

  {
    int rowB = find_child_row_gibbs(state->child, partners[chosen]);
    int swapParA = state->parent[rowB];
    int swapParB = state->parent[rowA];
    state->parent[rowA] = swapParA;
    state->parent[rowB] = swapParB;
    absLen[rowA] = state->treeLength * state->relBrLengths[rowB];
    absLen[rowB] = state->treeLength * state->relBrLengths[rowA];
    double tmpRel = state->relBrLengths[rowA];
    state->relBrLengths[rowA] = state->relBrLengths[rowB];
    state->relBrLengths[rowB] = tmpRel;

    auto po = TreeTools::preorder_weighted_impl(
        state->parent, state->child, absLen);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbsFinal = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbsFinal[k] / state->treeLength;
    }
  }

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
  state->ecoCL.invalidate_all();   // T-011: topology/branches changed
  state->wEdgeDirty = true;        // T-010: topology/branches changed
  return true;
}


// BranchBins struct now lives in mcmc_state.h and is precomputed by
// set_branch_bins() at MCMC init — no static globals or lazy init.


// ---------------------------------------------------------------------------
// weighted_branch_scale_impl  (M-087)
//
// WeightedBranchLengthScale: pick two branches, discretise the branch-
// fraction space into B bins (Beta(0.25,0.25) quantile breakpoints),
// evaluate log-likelihood at each bin midpoint, weight by exp(beta * LL),
// Gibbs-sample a bin, then draw a new fraction from a Beta centred on the
// selected midpoint.  Returns logHastings for standard MH acceptance.
//
// The bin-selection weights are symmetric (forward == reverse) because
// midpoint evaluations are independent of the current fraction.
// logHastings = log(w_oldBin) + logBeta(f_old|a_old,b_old)
//             - log(w_chosenBin) - logBeta(f_new|a_new,b_new).
// ---------------------------------------------------------------------------
static bool weighted_branch_scale_impl(
    ChainRng& rng, McmcData* data, McmcState* state, double beta,
    double& logHastings,
    int& outIdx1, int& outIdx2, double& outOld1, double& outOld2) {

  const int nEdge = state->relBrLengths.size();
  if (nEdge < 2) return false;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Pick two branches (same scheme as beta_simplex_impl)
  int index = static_cast<int>(rng.unif() * nEdge);
  if (index >= nEdge) index = nEdge - 1;
  int other = static_cast<int>(rng.unif() * (nEdge - 1));
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
                                       state->parent, state->child, trialAbs);
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
  double rnd = rng.unif() * sumW;
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
  double newF = rng.rbeta(alphaNew, betaNew);
  if (newF < 1e-8) newF = 1e-8;
  if (newF > 1.0 - 1e-8) newF = 1.0 - 1e-8;

  // 7. Hastings ratio
  //    Bin weights cancel (symmetric); within-bin densities remain.
  int oldBin = nBins - 1;
  for (int b = 0; b < nBins; ++b) {
    if (oldF <= bins.breaks[b + 1]) { oldBin = b; break; }
  }
  const double oldMid = bins.mids[oldBin];
  const double alphaOld = oldMid * conc + 1.0;
  const double betaOld  = (1.0 - oldMid) * conc + 1.0;

  logHastings = std::log(weights[oldBin])
              + R::dbeta(oldF, alphaOld, betaOld, 1)
              - std::log(weights[chosenBin])
              - R::dbeta(newF, alphaNew, betaNew, 1);

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
    ChainRng& rng, McmcData* data, McmcState* state, double beta) {

  const int nEdge = state->relBrLengths.size();
  if (nEdge < 2) return false;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;
  const double conc = bins.concentration;

  // Fisher-Yates shuffle for random permutation scan
  std::vector<int> perm(nEdge);
  for (int i = 0; i < nEdge; ++i) perm[i] = i;
  for (int i = nEdge - 1; i > 0; --i) {
    int j = static_cast<int>(rng.unif() * (i + 1));
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
    int other = static_cast<int>(rng.unif() * (nEdge - 1));
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
                                         state->parent, state->child, trialAbs);
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
    double rnd = rng.unif() * sumW;
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
    double newF = rng.rbeta(alphaNew, betaNew);
    if (newF < 1e-8) newF = 1e-8;
    if (newF > 1.0 - 1e-8) newF = 1.0 - 1e-8;

    // Hastings ratio: bin weights cancel; within-bin Beta densities remain
    int oldBin = nBins - 1;
    for (int b = 0; b < nBins; ++b) {
      if (oldF <= bins.breaks[b + 1]) { oldBin = b; break; }
    }
    const double oldMid = bins.mids[oldBin];
    const double alphaOld = oldMid * conc + 1.0;
    const double betaOld  = (1.0 - oldMid) * conc + 1.0;

    double logHastings = std::log(weights[oldBin])
                       + R::dbeta(oldF, alphaOld, betaOld, 1)
                       - std::log(weights[chosenBin])
                       - R::dbeta(newF, alphaNew, betaNew, 1);

    if (!R_FINITE(logHastings)) continue;

    // Evaluate LL at proposed fraction
    trialAbs[index] = newF * absTotal;
    trialAbs[other] = (1.0 - newF) * absTotal;
    double proposedLL = compute_full_loglik_at(
      *data, *state, state->parent, state->child, trialAbs);
    if (!R_FINITE(proposedLL)) continue;

    // MH accept/reject (prior is constant for relBrLengths)
    double logAlpha = beta * (proposedLL - currentLL) + logHastings;

    if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
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
    state->ecoCL.invalidate_all();   // T-011: branch lengths changed
    state->wEdgeDirty = true;        // T-010: branch lengths changed
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
// Hastings ratio: topology selection cancels (Z = Z' by symmetry of the
// candidate set).  Remaining branch-fraction component:
//   logHR = log(w_{self,b_old}) + logBeta(f_old | b_old)
//         - log(w_{chosen,b_new}) - logBeta(f_new | b_new)
//
// Cost: O(N × B) likelihood evaluations + 1 for the final proposed state.
// ---------------------------------------------------------------------------
static bool weighted_spr_impl(ChainRng& rng, McmcData* data, McmcState* state,
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
  int pickIdx = (int)(rng.unif() * (double)eligible.size());
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
                                          trialAbs);
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
                                              ordPar, ordCh, ordAbs);
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
  double rnd = rng.unif() * sumM;
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
    double rndBin = rng.unif() * mCand[chosen];
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
  double fNew = rng.rbeta(alphaNew, betaNew);
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
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbsFinal);
  if (!R_FINITE(newLogLik)) return false;

  // 14. Hastings ratio (branch-fraction component only; topology cancels)
  int oldBin = nBins - 1;
  for (int b = 0; b < nBins; ++b) {
    if (fOld <= bins.breaks[b + 1]) { oldBin = b; break; }
  }
  const double oldMid   = bins.mids[oldBin];
  const double alphaOld = oldMid * conc + 1.0;
  const double betaOld  = (1.0 - oldMid) * conc + 1.0;

  double logHR = std::log(std::max(selfW[oldBin], 1e-300))
               + R::dbeta(fOld, alphaOld, betaOld, 1)
               - std::log(std::max(candW[chosen][chosenBin], 1e-300))
               - R::dbeta(fNew, alphaNew, betaNew, 1);

  // 15. Prior at proposed state
  NumericVector propRelBr(nEdge);
  for (int k = 0; k < nEdge; ++k)
    propRelBr[k] = ordAbsFinal[k] / state->treeLength;

  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, propRelBr,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
  if (!R_FINITE(newLogPrior)) return false;

  // 16. MH acceptance
  double logAlpha = beta * (newLogLik - state->logLik)
                  + (newLogPrior - state->logPrior) + logHR;
  if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = propRelBr[k];
    }
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik.clear();
    state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
    state->ecoCL.invalidate_all();   // T-011: topology/branches changed
    state->wEdgeDirty = true;        // T-010: topology/branches changed
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
// for the chosen partner.  MH acceptance corrects the approximation.
//
// Cost: O(N × B) likelihood evaluations.
// ---------------------------------------------------------------------------

// (Uses find_child_row_gibbs defined above)

static bool weighted_subtree_swap_impl(ChainRng& rng, McmcData* data, McmcState* state,
                                        double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Pick a random node (any edge child)
  int pickIdx = (int)(rng.unif() * (double)nEdge);
  if (pickIdx >= nEdge) pickIdx = nEdge - 1;
  const int nodeA = state->child[pickIdx];

  // 2. Get valid swap partners
  std::vector<int> partners = get_valid_swap_partners_impl(
      state->parent, state->child, nTip, nodeA);
  if (partners.empty()) return false;
  const int nPart = (int)partners.size();

  // 3. Find rowA
  const int rowA = find_child_row_gibbs(state->child, nodeA);
  if (rowA < 0) return false;

  // 4. Absolute edge lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // 5. Candidate marginals: in-place topology + preorder_into (M-109)
  std::vector<std::vector<double>> candLL(nPart,
                                           std::vector<double>(nBins));
  std::vector<double> candMax(nPart, R_NegInf);
  std::vector<int> rowBs(nPart);       // edge row for each partner
  std::vector<double> totals(nPart);   // brA + brB_i

  // Working copy of parent (cloned ONCE); child unchanged in subtree swap
  IntegerVector workPar = clone(state->parent);
  IntegerVector ordPar(nEdge), ordCh(nEdge);
  NumericVector ordAbs(nEdge);

  const int origParA = workPar[rowA];
  const double origAbsA = absLen[rowA];

  for (int pi = 0; pi < nPart; ++pi) {
    int rowB = find_child_row_gibbs(state->child, partners[pi]);
    if (rowB < 0) { candMax[pi] = R_NegInf; rowBs[pi] = -1; continue; }
    rowBs[pi]  = rowB;
    totals[pi] = absLen[rowA] + absLen[rowB];

    int origParB = workPar[rowB];
    double origAbsB = absLen[rowB];

    // Apply swap in-place
    workPar[rowA] = origParB;
    workPar[rowB] = origParA;

    for (int b = 0; b < nBins; ++b) {
      absLen[rowA] = bins.mids[b] * totals[pi];
      absLen[rowB] = (1.0 - bins.mids[b]) * totals[pi];

      preorder_into(workPar, state->child, absLen, nTip,
                    INTEGER(ordPar), INTEGER(ordCh), REAL(ordAbs));
      candLL[pi][b] = compute_full_loglik_at(*data, *state,
                                              ordPar, ordCh, ordAbs);
      if (R_FINITE(candLL[pi][b]) && candLL[pi][b] > candMax[pi])
        candMax[pi] = candLL[pi][b];
    }

    // Restore
    workPar[rowA] = origParA;
    workPar[rowB] = origParB;
    absLen[rowA]  = origAbsA;
    absLen[rowB]  = origAbsB;
  }

  // 6. Compute marginal weights with global offset
  double globalMax = state->logLik;
  for (int pi = 0; pi < nPart; ++pi)
    if (candMax[pi] > globalMax) globalMax = candMax[pi];
  if (!R_FINITE(globalMax)) return false;

  double wOrig = std::exp(beta * (state->logLik - globalMax));
  std::vector<double> mCand(nPart);
  std::vector<std::vector<double>> candW(nPart,
                                          std::vector<double>(nBins));
  double sumM = wOrig;
  for (int pi = 0; pi < nPart; ++pi) {
    mCand[pi] = 0.0;
    for (int b = 0; b < nBins; ++b) {
      candW[pi][b] = R_FINITE(candLL[pi][b]) ?
                       std::exp(beta * (candLL[pi][b] - globalMax)) : 0.0;
      mCand[pi] += candW[pi][b];
    }
    sumM += mCand[pi];
  }
  if (sumM <= 0.0) return false;

  // 7. Sample: self-draw → no-op
  double rnd = rng.unif() * sumM;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nPart - 1;
  for (int pi = 0; pi < nPart - 1; ++pi) {
    if (rnd < mCand[pi]) { chosen = pi; break; }
    rnd -= mCand[pi];
  }
  if (rowBs[chosen] < 0) return false;

  // 8. Sample bin within chosen candidate
  int chosenBin = nBins - 1;
  {
    double rndBin = rng.unif() * mCand[chosen];
    double cum = 0.0;
    for (int b = 0; b < nBins; ++b) {
      cum += candW[chosen][b];
      if (rndBin < cum) { chosenBin = b; break; }
    }
  }

  // 9. Draw fraction from Beta centred on chosen bin's midpoint
  const double conc = bins.concentration;
  const double chosenMid = bins.mids[chosenBin];
  const double alphaNew = chosenMid * conc + 1.0;
  const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
  double fNew = rng.rbeta(alphaNew, betaNew);
  if (fNew < 1e-8) fNew = 1e-8;
  if (fNew > 1.0 - 1e-8) fNew = 1.0 - 1e-8;

  // 10. Construct final proposed topology at fNew (using working copy)
  const int rowB   = rowBs[chosen];
  const double tot = totals[chosen];
  workPar[rowA] = state->parent[rowB];
  workPar[rowB] = state->parent[rowA];
  absLen[rowA] = fNew * tot;
  absLen[rowB] = (1.0 - fNew) * tot;

  auto po = TreeTools::preorder_weighted_impl(workPar, state->child, absLen);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbsFinal = po.second;
  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbsFinal);
  if (!R_FINITE(newLogLik)) return false;

  // 11. Hastings ratio
  //     f_default = absLen[rowB] / tot (the default swap fraction)
  const double fDefault = (tot > 0.0) ? absLen[rowB] / tot : 0.5;
  int defaultBin = nBins - 1;
  for (int b = 0; b < nBins; ++b) {
    if (fDefault <= bins.breaks[b + 1]) { defaultBin = b; break; }
  }
  const double defaultMid  = bins.mids[defaultBin];
  const double alphaOld    = defaultMid * conc + 1.0;
  const double betaOld     = (1.0 - defaultMid) * conc + 1.0;

  double logHR = std::log(std::max(wOrig, 1e-300))
               + R::dbeta(fDefault, alphaOld, betaOld, 1)
               - std::log(std::max(candW[chosen][chosenBin], 1e-300))
               - R::dbeta(fNew, alphaNew, betaNew, 1);

  // 12. Prior at proposed state
  NumericVector propRelBr(nEdge);
  for (int k = 0; k < nEdge; ++k)
    propRelBr[k] = ordAbsFinal[k] / state->treeLength;

  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, propRelBr,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
  if (!R_FINITE(newLogPrior)) return false;

  // 13. MH acceptance
  double logAlpha = beta * (newLogLik - state->logLik)
                  + (newLogPrior - state->logPrior) + logHR;
  if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = propRelBr[k];
    }
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik.clear();
    state->nodeCL.invalidate_all();  // M-143/M-161: topology/branches changed
    state->ecoCL.invalidate_all();   // T-011: topology/branches changed
    state->wEdgeDirty = true;        // T-010: topology/branches changed
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
    case 2: state->rateLogSd = val; break;
    case 3: state->rateNeo = val; break;
    case 4: state->betaScale = val; break;
  }
}

// Evaluate beta * logLik + logPrior for current state, using partial cache
// when the parameter only affects a subset of partitions.
static double eval_slice_target(McmcData* data, McmcState* state,
                                int paramIdx, double beta) {
  double logPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
  if (!R_FINITE(logPrior)) return R_NegInf;

  double logLik;
  int nEdge = state->parent.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  bool hasPLC = !state->partLogLik.empty();
  ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

  if (data->ecologyAware) {
    // T-010: partition cache restored for ecology. Tree-length (paramIdx 0)
    // changes every edge length and therefore invalidates wEdge; rate_loss
    // and rate_neo touch only neomorphic partitions; rate_log_sd changes
    // ACRV rates for every partition. wEdge is independent of all four
    // slice params, so we just need to recompute the right partitions.
    // This is a probe — DO NOT mutate state->partLogLik (mirror blind path).
    if (paramIdx == 0) state->wEdgeDirty = true;
    bool hasPLCEco = !state->partLogLik.empty();
    eco_refresh_wedge(data, state, edgeLen);
    std::vector<double> gammaE;
    compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                            gammaE);
    NumericVector rates = (state->rateLogSd > 0.0)
      ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
      : NumericVector(1, 1.0);
    if (hasPLCEco && (paramIdx == 1 || paramIdx == 3)) {
      // rate_loss / rate_neo: neomorphic partitions only.
      logLik = state->logLik;
      for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
        int pi = data->neoPartIndices[ni];
        double newPart = cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
        logLik += (newPart - state->partLogLik[pi]);
      }
    } else {
      // Tree length, rate_log_sd, or stale cache: recompute every partition,
      // reusing cached wEdge when clean.
      int nParts = (int)data->parts.size();
      logLik = 0.0;
      for (int pi = 0; pi < nParts; ++pi) {
        logLik += cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
      }
    }
  } else if (hasPLC && (paramIdx == 1 || paramIdx == 3)) {
    // rate_loss (1), rate_neo (3): only neomorphic partitions change
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
  } else {
    logLik = cpp_log_likelihood(
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
static bool slice_scalar_impl(ChainRng& rng, McmcData* data, McmcState* state,
                               int paramIdx, double width,
                               double beta, int maxSteps = 10,
                               int* nExpansionsOut = nullptr) {
  double x0 = get_scalar(state, paramIdx);

  // Work on log scale: u = log(x).
  // Target includes Jacobian: log f(u) = beta*logLik + logPrior + u.
  double u0 = std::log(x0);
  double logY0 = beta * state->logLik + state->logPrior + u0;

  // Slice height
  double logZ = logY0 + std::log(rng.unif());

  // Stepping out on log scale (count expansions for width adaptation)
  int nExp = 0;
  double L = u0 - width * rng.unif();
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
    double u1 = L + rng.unif() * (R_bound - L);
    double x1 = std::exp(u1);
    set_scalar(state, paramIdx, x1);
    double logTarget1 = eval_slice_target(data, state, paramIdx, beta) + u1;
    if (logTarget1 >= logZ) {
      // Accept — recompute and cache logLik / logPrior / partLogLik
      state->logPrior = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);

      int nEdge = state->parent.size();
      NumericVector edgeLen(nEdge);
      for (int i = 0; i < nEdge; ++i)
        edgeLen[i] = state->treeLength * state->relBrLengths[i];
      ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

      bool hasPLC = !state->partLogLik.empty();
      if (data->ecologyAware) {
        // T-010: refresh partition cache + state->logLik. wEdge cache is
        // reused unless paramIdx == 0 (tree length) which invalidates it.
        if (paramIdx == 0) state->wEdgeDirty = true;
        bool hasPLCEco = !state->partLogLik.empty();
        if (hasPLCEco && (paramIdx == 1 || paramIdx == 3)) {
          eco_refresh_wedge(data, state, edgeLen);
          std::vector<double> gammaE;
          compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                                  gammaE);
          NumericVector rates = (state->rateLogSd > 0.0)
            ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
            : NumericVector(1, 1.0);
          for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
            int pi = data->neoPartIndices[ni];
            state->partLogLik[pi] = cpp_partition_log_likelihood_ecology(
              *data, pi, state->parent, state->child, edgeLen,
              state->kPrime, state->rateLoss, state->rateNeo,
              state->phi, state->zMatrix,
              state->wEdge, gammaE, rates);
          }
          state->logLik = 0.0;
          for (size_t pi = 0; pi < state->partLogLik.size(); ++pi)
            state->logLik += state->partLogLik[pi];
        } else {
          state->logLik = eco_recompute_all_partitions(data, state, edgeLen);
        }
      } else if (hasPLC && (paramIdx == 1 || paramIdx == 3)) {
        // Update only neo partitions in cache
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
      } else {
        state->logLik = cpp_log_likelihood(
          *data, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateLogSd,
          state->rateNeo, state->betaScale, wsPtr);
        // Invalidate partition cache (full recompute was done)
        state->partLogLik.clear();
      }
      // M-145/M-161: Invalidate node CL cache — slice changed a model
      // parameter.  Granular: only invalidate affected units.
      switch (paramIdx) {
        case 1: case 3:  // rate_loss, rate_neo: only neomorphic units
          state->nodeCL.invalidate_neo_cls();
          state->ecoCL.invalidate_neo_cls();   // T-011
          break;
        case 2:  // rateLogSd: ACRV rates change, all units
          state->nodeCL.invalidate_all_cls();
          state->ecoCL.invalidate_all_cls();   // T-011
          break;
        default:  // tree_length (0), beta_scale (4)
          state->nodeCL.invalidate_all();
          state->ecoCL.invalidate_all();       // T-011
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
// Prior-only slice sampler for Beta-Geometric hyperparameters (M-163)
//
// Samples kprimeAlpha (paramCode=0) or kprimeBeta (paramCode=1) using a
// univariate slice sampler on the log scale.  Target = logPrior + log(x)
// (the log(x) Jacobian arises from sampling u = log(x) and transforming).
// No likelihood evaluation needed — α/β only affect the prior.
// ---------------------------------------------------------------------------
static bool slice_kprime_hyper_impl(ChainRng& rng, McmcData* data, McmcState* state,
                                     int paramCode, double width,
                                     int maxSteps = 10,
                                     int* nExpansionsOut = nullptr) {
  double x0 = (paramCode == 0) ? state->kprimeAlpha : state->kprimeBeta;
  if (x0 <= 0.0) return false;

  // Work on log scale: u = log(x)
  double u0 = std::log(x0);

  // Target: logPrior(current) + u (Jacobian)
  double logY0 = state->logPrior + u0;
  double logZ = logY0 + std::log(rng.unif());

  // Helper lambda to evaluate target at a candidate u
  auto evalTarget = [&](double u) -> double {
    double xCand = std::exp(u);
    if (xCand <= 0.0 || !R_FINITE(xCand)) return R_NegInf;
    double oldVal = (paramCode == 0) ? state->kprimeAlpha : state->kprimeBeta;
    if (paramCode == 0) state->kprimeAlpha = xCand;
    else                state->kprimeBeta  = xCand;
    double lp = cpp_log_prior(
      *data, state->treeLength, state->relBrLengths,
      state->rateLoss, state->rateLogSd, state->rateNeo,
      state->p, state->kPrime, state->betaScale,
      state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
    // Restore
    if (paramCode == 0) state->kprimeAlpha = oldVal;
    else                state->kprimeBeta  = oldVal;
    return lp + u;  // logPrior + Jacobian
  };

  // Stepping out
  int nExp = 0;
  double L = u0 - width * rng.unif();
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
    double u1 = L + rng.unif() * (R_bound - L);
    double logTarget1 = evalTarget(u1);
    if (logTarget1 >= logZ) {
      // Accept
      double x1 = std::exp(u1);
      if (paramCode == 0) state->kprimeAlpha = x1;
      else                state->kprimeBeta  = x1;
      // Recompute and cache logPrior
      state->logPrior = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
        state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
      return true;
    }
    if (u1 < u0) L = u1; else R_bound = u1;
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
    ChainRng& rng,
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

  int pick = (int)(rng.unif() * (double)eligible.size());
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
  double rnd = rng.unif() * sumW;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < w[ci]) { chosen = ci; break; }
    rnd -= w[ci];
  }

  // 9. Apply the chosen SPR
  const int regraftRow = candidates[chosen];
  const int b = stateChild[regraftRow];
  const double tau = rng.unif();

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
// Gibbs kPrime sweep (moveType 25)
//
// Samples each k'_i from its full conditional in a single random-order scan
// of all transformational characters. Always accepts (Gibbs update).
//
// Red-team S-1 follow-up (2026-05-21): under ecologyAware==TRUE, the blind
// per-site pruners ignore phi/pi0/theta/zMatrix/wEdge and therefore draw
// from a wrong full conditional. The ecology branch below evaluates the
// candidate weights via per_char_log_lik_ecology, the same helper used by
// gibbs_z_sweep_impl, which produces fully-corrected per-character log-
// likelihoods (pruner + asc + relabel) under the ecology mixture model.
// ---------------------------------------------------------------------------

// Forward declaration of the ecology branch.
static bool gibbs_kprime_sweep_impl_ecology(ChainRng& rng, McmcData* data, McmcState* state,
                                            double beta);

static bool gibbs_kprime_sweep_impl(ChainRng& rng, McmcData* data, McmcState* state,
                                     double beta) {
  int nTrans = (int)data->transIdxGlobal.size();
  if (nTrans == 0) return false;

  // Ecology-aware path uses a dedicated implementation that evaluates the
  // candidate weights under the full ecology mixture conditional.
  if (data->ecologyAware) {
    return gibbs_kprime_sweep_impl_ecology(rng, data, state, beta);
  }

  // Pre-compute absolute edge lengths
  int nEdge = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  // Pre-compute ACRV rates
  bool useAcrv = (state->rateLogSd > 0.0);
  NumericVector acrvRates = useAcrv
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);

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
        *data, state->parent, state->child, edgeLen,
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

  static const int K_MAX_CAND = 50;   // absolute cap (fallback path)
  // M-164: tightened from -57.5 to -25.0.  exp(-25) ≈ 1.4e-11 relative
  // probability; even 50 such candidates contribute ~7e-10 total mass,
  // far below double-precision RNG resolution (~2.2e-16).
  static const double LOG_CUTOFF = -25.0;

  // Beta-Geometric: precompute incremental log-prior for ko = 0..K_MAX_CAND-1
  // logPrior(u=0) = log(α) - log(α+β)
  // logPrior(u=k) = logPrior(u=k-1) + log(β+k-1) - log(α+β+k)
  std::vector<double> bgLogPrior;
  if (isBetaGeometric) {
    double a = state->kprimeAlpha;
    double b = state->kprimeBeta;
    bgLogPrior.resize(K_MAX_CAND);
    bgLogPrior[0] = std::log(a) - std::log(a + b);
    for (int ko = 1; ko < K_MAX_CAND; ++ko) {
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

  // Ensure Gibbs workspace is large enough for the worst-case stride.
  // M-172: stride is nUniq × k (not nChar × k) — savings proportional to redundancy.
  int maxNode = 2 * data->nTip - 1;
  int gibbsMaxStride = 0;
  for (auto& tp : transParts) {
    int s = tp.nUniq * (tp.kObs + K_MAX_CAND);
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
      int hi = tp_.kObs + K_MAX_CAND - 1;
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

  // Per-character sampling state
  // logW stores weights for each candidate; nCand tracks how many
  std::vector<double> charMaxLogW(nTrans, R_NegInf);
  // M-164: track best corrected log-likelihood per character for pre-filter
  std::vector<double> charMaxLL(nTrans, R_NegInf);
  // Flat: logW[ti * K_MAX_CAND + ko]
  std::vector<double> charLogW(nTrans * K_MAX_CAND, R_NegInf);
  std::vector<int> charNCand(nTrans, 0);
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
  for (int ko = 0; ko < K_MAX_CAND && nActive > 0; ++ko) {

    for (int pi = 0; pi < (int)transParts.size(); ++pi) {
      auto& tp = transParts[pi];
      auto& pa = partAct[pi];
      int nAct = (int)pa.activePatterns.size();
      if (nAct == 0) continue;

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
          if (optimisticW < charMaxLogW[ti0] + LOG_CUTOFF) {
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
      // gibbsMaxStride is already worst-case (nUniq*(kObs+K_MAX_CAND)) so the
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
            state->parent, state->child, edgeLen, part.uniqueTipStates,
            k, 1.0, hetBins, nBC, rates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else if (useCollapse) {
          pruning_jc_acrv_persite_collapsed(
            state->parent, state->child, edgeLen, part.uniqueTipStates,
            k, tp.kObs, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else {
          pruning_jc_acrv_persite(
            state->parent, state->child, edgeLen, part.uniqueTipStates,
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
            state->parent, state->child, edgeLen, sub,
            k, 1.0, hetBins, nBC, rates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else if (useCollapse) {
          pruning_jc_acrv_persite_collapsed(
            state->parent, state->child, edgeLen, sub,
            k, tp.kObs, acrvRates,
            state->gibbsWs.buf.data(), state->gibbsWs.init.data(),
            neededStride, siteLL.data());
        } else {
          pruning_jc_acrv_persite(
            state->parent, state->child, edgeLen, sub,
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
          charLogW[ti * K_MAX_CAND + ko] = w;
          charNCand[ti]++;
          if (w > charMaxLogW[ti]) charMaxLogW[ti] = w;
        }

        // Termination check: representative ti0 holds the same charMaxLogW as all
        // others in the group (identical weights throughout), so one check suffices.
        int ti0 = pa.patTrans[localPat][0];
        if (!R_FINITE(ll) || w < charMaxLogW[ti0] + LOG_CUTOFF) {
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

  // ---------------------------------------------------------------
  // Sampling phase: for each character, sample k' from the
  // precomputed log-weights in charLogW[].
  // ---------------------------------------------------------------

  // Random permutation of transformational character indices
  std::vector<int> perm(nTrans);
  for (int i = 0; i < nTrans; ++i) perm[i] = i;
  for (int i = nTrans - 1; i > 0; --i) {
    int j = static_cast<int>(rng.unif() * (i + 1));
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
    double* logW = &charLogW[ti * K_MAX_CAND];

    // Sample from categorical (log-sum-exp)
    double sumExp = 0.0;
    for (int c = 0; c < nCand; ++c)
      sumExp += std::exp(logW[c] - maxW);

    double u = rng.unif() * sumExp;
    double cum = 0.0;
    int chosen = nCand - 1;
    for (int c = 0; c < nCand; ++c) {
      cum += std::exp(logW[c] - maxW);
      if (cum >= u) { chosen = c; break; }
    }

    state->kPrime[gi] = kObs_i + chosen;
  }

  // Rebuild logLik, logPrior, and partition cache after sweep.
  // T-010: ecology partition cache restored. kPrime changes affect only
  // transformational partitions, so we recompute just those (kPrime subgroup
  // composition for one trans char can change — easiest correct option is
  // to recompute the whole trans partition).
  if (data->ecologyAware) {
    bool hasPLCEco = !state->partLogLik.empty();
    if (hasPLCEco) {
      eco_refresh_wedge(data, state, edgeLen);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      NumericVector rates = (state->rateLogSd > 0.0)
        ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
        : NumericVector(1, 1.0);
      double newLL = 0.0;
      int nP = (int)data->parts.size();
      for (int pi = 0; pi < nP; ++pi) {
        if (data->parts[pi].type == 1) {  // transformational only
          state->partLogLik[pi] = cpp_partition_log_likelihood_ecology(
            *data, pi, state->parent, state->child, edgeLen,
            state->kPrime, state->rateLoss, state->rateNeo,
            state->phi, state->zMatrix,
            state->wEdge, gammaE, rates);
        }
        newLL += state->partLogLik[pi];
      }
      state->logLik = newLL;
    } else {
      state->logLik = eco_recompute_all_partitions(data, state, edgeLen);
    }
  } else {
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
  }

  state->logPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);

  state->nodeCL.invalidate_structure();  // M-161: kPrime changed, unit structure may differ
  state->ecoCL.invalidate_structure();   // T-011: kPrime changed

  return true;  // Gibbs: always accept
}


// ---------------------------------------------------------------------------
// Gibbs kPrime sweep — ecology-aware branch (red-team S-1 follow-up).
//
// Mirrors the blind sweep's structure (per-character categorical over
// ko = 0..K_MAX_CAND-1, log-sum-exp sampling, prior families, M-164
// early termination) but evaluates each candidate's likelihood via
// per_char_log_lik_ecology — the same helper used by gibbs_z_sweep_impl —
// which honours phi/pi0/theta/zMatrix/wEdge and the ecology mixture.
//
// Notes vs. the blind path:
//   * M-172 unique-pattern compression is DROPPED. Per-character ecology
//     weights depend on the character's zRow (per-character), so two
//     characters sharing a tip-state pattern can have different ecology
//     likelihoods. Correctness > speed; this is a known cost.
//   * No external ascertainment or relabel correction is applied here:
//     per_char_log_lik_ecology already adds both (see mcmc_ecology.cpp
//     :1300-1346). Adding them externally would double-count.
//   * useHet is not supported in ecology mode (the ecology pruner has no
//     het variant — same constraint as gibbs_z_sweep_impl).
//   * Post-sweep recompute keeps the partition cache path from the blind
//     impl (lines that begin "if (data->ecologyAware)").
// ---------------------------------------------------------------------------
static bool gibbs_kprime_sweep_impl_ecology(McmcData* data, McmcState* state,
                                            double beta) {
  int nTrans = (int)data->transIdxGlobal.size();
  if (nTrans == 0) return false;
  if (data->codingType == 2) return false;  // informative not supported under ecology

  int nEdge = state->relBrLengths.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  // Hoist wEdge, gammaE, ACRV rates out of the per-cell loop (same pattern
  // as gibbs_z_sweep_impl).
  eco_refresh_wedge(data, state, edgeLen);
  std::vector<double> gammaE;
  compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta, gammaE);
  NumericVector rates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);
  int refE  = data->ecology.refEcology;
  int kEco  = data->ecology.kEcology;
  int zCols = std::max(0, kEco - 1);
  int nTip  = data->nTip;
  int maxNode = 2 * nTip - 1;

  // Prior components (identical to blind path)
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

  static const int K_MAX_CAND = 50;
  static const double LOG_CUTOFF = -25.0;

  // Beta-Geometric incremental log-prior
  std::vector<double> bgLogPrior;
  if (isBetaGeometric) {
    double a = state->kprimeAlpha;
    double b = state->kprimeBeta;
    bgLogPrior.resize(K_MAX_CAND);
    bgLogPrior[0] = std::log(a) - std::log(a + b);
    for (int ko = 1; ko < K_MAX_CAND; ++ko) {
      bgLogPrior[ko] = bgLogPrior[ko - 1]
                      + std::log(b + ko - 1)
                      - std::log(a + b + ko);
    }
  }

  // Per-trans bookkeeping (kObs, globalCharIdx) and max candidate k for
  // workspace sizing.
  std::vector<int> tiKObs(nTrans, 0);
  std::vector<int> tiGi(nTrans, 0);
  int maxCandK = 0;
  for (int ti = 0; ti < nTrans; ++ti) {
    int gi = data->transIdxGlobal[ti];
    tiGi[ti]   = gi;
    tiKObs[ti] = data->kObs[gi];
    int hi = tiKObs[ti] + K_MAX_CAND - 1;
    if (hi > maxCandK) maxCandK = hi;
  }

  // Empirical-geometric: tabulate log P(k' = m) once.
  std::vector<double> egLogPriorByK;
  if (isEmpGeom) {
    egLogPriorByK.assign(maxCandK + 2, R_NegInf);
    double logQ = (data->empTailDecay > 0.0)
                  ? std::log(data->empTailDecay) : R_NegInf;
    int bodyLen = (int)data->empLogBody.size();
    std::vector<double> terms;
    terms.reserve(64);
    for (int m = 2; m <= maxCandK + 1; ++m) {
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

  // Lambda computing the prior contribution at candidate (kObs + ko).
  auto logPriorAt = [&](int kObs_, int ko) -> double {
    int k = kObs_ + ko;
    if (isBetaGeometric) {
      return bgLogPrior[ko];
    } else if (isGeometric) {
      return logP + ko * log1mP;
    } else if (isEmpGeom) {
      return (k >= 0 && k < (int)egLogPriorByK.size())
             ? egLogPriorByK[k] : R_NegInf;
    } else {
      return k * lsLogC - std::log(static_cast<double>(k)) - lsLogNorm;
    }
  };

  // Allocate one workspace large enough for any candidate k we'll evaluate.
  // per_char_log_lik_ecology still allocates an internal buf per call (see
  // mcmc_ecology.cpp ~1292), but the GibbsZWorkspace covers its zPart/tipStates.
  GibbsZWorkspace ws;
  ws.ensure(maxNode, maxCandK, nTip, zCols);

  IntegerVector zRow(zCols);

  // Per-character logW table (flat) and termination flags.
  std::vector<double> charLogW(nTrans * K_MAX_CAND, R_NegInf);
  std::vector<double> charMaxLogW(nTrans, R_NegInf);
  std::vector<double> charMaxLL(nTrans, R_NegInf);
  std::vector<int>    charNCand(nTrans, 0);
  std::vector<bool>   terminated(nTrans, false);

  // Outer loop: build per-character weights ko = 0..K_MAX_CAND-1, with M-164
  // early termination per-character (no pattern grouping).
  for (int ko = 0; ko < K_MAX_CAND; ++ko) {
    bool anyActive = false;
    for (int ti = 0; ti < nTrans; ++ti) {
      if (terminated[ti]) continue;
      anyActive = true;
      int gi   = tiGi[ti];
      int kObs = tiKObs[ti];
      int k    = kObs + ko;

      double logPrior_k = logPriorAt(kObs, ko);

      // M-164 pre-filter: if even the best plausible likelihood seen so far
      // would put this candidate below the current best by LOG_CUTOFF, drop
      // the character from further consideration. Identical logic to the
      // blind path but on a per-character basis (no pattern groups).
      if (ko >= 2) {
        double optimisticW = beta * charMaxLL[ti] + logPrior_k;
        if (optimisticW < charMaxLogW[ti] + LOG_CUTOFF) {
          terminated[ti] = true;
          continue;
        }
      }

      // Pull zRow for this character.
      for (int j = 0; j < zCols; ++j) zRow[j] = state->zMatrix(gi, j);

      // Per-character ecology log-lik at candidate k. This call already
      // applies the ascertainment correction (codingType == 1) and the
      // Mk' relabel correction (data.relabel) internally — do NOT add them
      // externally here or they will be double-counted.
      double ll = per_char_log_lik_ecology(
        *data, gi, state->parent, state->child, edgeLen,
        k, state->rateLoss, state->rateNeo,
        rates, state->phi, zRow, state->wEdge,
        refE, gammaE, ws);

      double w = (R_FINITE(ll)) ? (beta * ll + logPrior_k) : R_NegInf;
      charLogW[ti * K_MAX_CAND + ko] = w;
      charNCand[ti]++;
      if (R_FINITE(ll) && ll > charMaxLL[ti]) charMaxLL[ti] = ll;
      if (R_FINITE(w) && w > charMaxLogW[ti]) charMaxLogW[ti] = w;

      // Per-character M-164 termination on observed weight.
      if (!R_FINITE(ll) || w < charMaxLogW[ti] + LOG_CUTOFF) {
        terminated[ti] = true;
      }
    }
    if (!anyActive) break;
  }

  // Sampling phase (identical structure to the blind path).
  std::vector<int> perm(nTrans);
  for (int i = 0; i < nTrans; ++i) perm[i] = i;
  for (int i = nTrans - 1; i > 0; --i) {
    int j = static_cast<int>(R::unif_rand() * (i + 1));
    if (j > i) j = i;
    std::swap(perm[i], perm[j]);
  }

  for (int si = 0; si < nTrans; ++si) {
    int ti = perm[si];
    int gi = tiGi[ti];
    int kObs_i = tiKObs[ti];
    int nCand = charNCand[ti];
    if (nCand == 0) continue;

    double maxW = charMaxLogW[ti];
    if (!R_FINITE(maxW)) continue;  // no usable candidate
    double* logW = &charLogW[ti * K_MAX_CAND];

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

  // Post-sweep: refresh wEdge (kPrime doesn't change wEdge but harmless to
  // verify) and recompute partition cache + logLik via the ecology path.
  // Mirror the blind path's post-sweep ecology block.
  bool hasPLCEco = !state->partLogLik.empty();
  if (hasPLCEco) {
    eco_refresh_wedge(data, state, edgeLen);
    compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                            gammaE);
    double newLL = 0.0;
    int nP = (int)data->parts.size();
    for (int pi = 0; pi < nP; ++pi) {
      if (data->parts[pi].type == 1) {  // transformational only
        state->partLogLik[pi] = cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
      }
      newLL += state->partLogLik[pi];
    }
    state->logLik = newLL;
  } else {
    state->logLik = eco_recompute_all_partitions(data, state, edgeLen);
  }

  state->logPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);

  state->nodeCL.invalidate_structure();
  state->ecoCL.invalidate_structure();   // T-011: kPrime changed

  return true;
}


// ---------------------------------------------------------------------------
// Block kPrime shift (moveType 26)
//
// Proposes shifting ALL transformational characters by the same integer delta.
// Standard MH acceptance with symmetric proposal.
// ---------------------------------------------------------------------------

static bool block_kprime_shift_impl(ChainRng& rng, McmcData* data, McmcState* state,
                                     int intWalkWindow, double beta) {
  int nTrans = (int)data->transIdxGlobal.size();
  if (nTrans == 0) return false;

  // Propose delta ~ Uniform({-W, ..., W})
  int range = 2 * intWalkWindow + 1;
  int delta = static_cast<int>(rng.unif() * range) - intWalkWindow;
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
  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);

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

  if (data->ecologyAware) {
    // T-010: block kPrime shift — transformational partitions only;
    // wEdge unchanged. Reuse cached partLogLik for non-trans partitions.
    bool hasPLCEco = !state->partLogLik.empty();
    if (hasPLCEco) {
      newPC = state->partLogLik;
      newLogLik = state->logLik;
      eco_refresh_wedge(data, state, edgeLen);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      NumericVector rates = (state->rateLogSd > 0.0)
        ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
        : NumericVector(1, 1.0);
      for (int pi = 0; pi < nParts; ++pi) {
        if (data->parts[pi].type == 1) {
          double v = cpp_partition_log_likelihood_ecology(
            *data, pi, state->parent, state->child, edgeLen,
            state->kPrime, state->rateLoss, state->rateNeo,
            state->phi, state->zMatrix,
            state->wEdge, gammaE, rates);
          newLogLik += (v - newPC[pi]);
          newPC[pi] = v;
        }
      }
    } else {
      // No cached PLC — compute all partitions into a local vector
      // (don't mutate state in case the move is rejected).
      newPC.resize(nParts);
      newLogLik = 0.0;
      eco_refresh_wedge(data, state, edgeLen);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      NumericVector rates = (state->rateLogSd > 0.0)
        ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
        : NumericVector(1, 1.0);
      for (int pi = 0; pi < nParts; ++pi) {
        newPC[pi] = cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
        newLogLik += newPC[pi];
      }
    }
  } else if (hasPLC) {
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

  if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
    state->logLik = newLogLik;
    state->logPrior = newLogPrior;
    state->partLogLik = std::move(newPC);
    state->nodeCL.invalidate_structure();  // M-161: kPrime changed
    state->ecoCL.invalidate_structure();   // T-011: kPrime changed
    return true;
  }

  // Rollback
  for (int i = 0; i < nTrans; ++i)
    state->kPrime[data->transIdxGlobal[i]] = oldKPrime[i];
  return false;
}


// ---------------------------------------------------------------------------
// gibbs_z_sweep_impl (Phase 3f) — cell-by-cell Gibbs sweep over the
// per-(character, ecology) latent z_{c, s} ∈ {0=none, 1=encouraged,
// 2=discouraged} under the spike-and-slab prior P(z = 0) = pi0,
// P(z = 1) = P(z = 2) = (1 - pi0) / 2.
//
// For each cell:
//   1. Enumerate the three candidate values.
//   2. Evaluate the per-character log-likelihood under each
//      (everything else fixed), tempered by the chain's beta.
//   3. Sample from the resulting categorical (softmax of
//      lp_v + beta * ll_v).
//
// wEdge and ACRV rates are hoisted out of the per-cell loop.  After the
// sweep, recompute the full logLik / logPrior in one shot to avoid drift
// from any small numerical mismatch between the per-character helper and
// the orchestrator (advisor recommendation).
// ---------------------------------------------------------------------------

static bool gibbs_z_sweep_impl(McmcData* data, McmcState* state, double beta) {
  if (!data->ecologyAware) return false;
  if (state->zMatrix.nrow() == 0 || state->zMatrix.ncol() == 0) return false;
  if (data->codingType == 2) return false;  // informative not supported

  int nChar = data->nChar;
  int kEco  = data->ecology.kEcology;
  int refE  = data->ecology.refEcology;
  int mode  = data->magnitudeMode;
  int nEdge = state->relBrLengths.size();
  int zCols = kEco - 1;
  if (zCols <= 0) return false;

  // Absolute edge lengths (relBr × treeLength)
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  // Hoisted: per-edge ecology weights, ACRV rate categories.
  NumericMatrix wEdge(nEdge, kEco);
  recompute_w_edge(*data, state->parent, state->child, edgeLen, wEdge);

  NumericVector rates = (state->rateLogSd > 0.0)
    ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
    : NumericVector(1, 1.0);

  // v2: gammaE depends on theta and pi0; constant across z choices for fixed cell.
  std::vector<double> gammaE(kEco, 1.0);
  for (int s = 0; s < kEco; ++s) {
    if (s == refE) { gammaE[s] = 1.0; continue; }
    int j = (s < refE) ? s : (s - 1);
    double phi_s = (mode == 0) ? state->phi[0] : state->phi[s];
    double th = (j < state->theta.size()) ? state->theta[j] : 0.5;
    gammaE[s] = state->pi0 + (1.0 - state->pi0) *
                (th * phi_s + (1.0 - th) / phi_s);
  }

  // Allocate ONE workspace for the whole sweep — eliminates per-call heap allocs.
  // Size to the maximum stride (= kStates) needed across all partitions.
  int nTip    = data->nTip;
  int maxNode = 2 * nTip - 1;
  int maxStride = 0;
  for (const PartInfo& part : data->parts) {
    int s = (part.type == 0) ? 2 : part.k;
    if (s > maxStride) maxStride = s;
  }
  // Transformational partitions can have kPrime > part.k (which is 0 for type 1).
  // Use current kPrime max for type-1 parts.
  for (int c = 0; c < nChar; ++c) {
    int pi = data->charToPartition[c];
    if (pi < 0) continue;
    if (data->parts[pi].type == 1) {
      int kp = state->kPrime[c];
      if (kp > maxStride) maxStride = kp;
    }
  }
  GibbsZWorkspace ws;
  if (maxStride > 0) ws.ensure(maxNode, maxStride, nTip, zCols);

  IntegerVector zRow(zCols);
  double log_pi0    = std::log(state->pi0);
  double log_1m_pi0 = std::log1p(-state->pi0);

  for (int c = 0; c < nChar; ++c) {
    int pi = data->charToPartition[c];
    if (pi < 0) continue;
    int kp = (data->parts[pi].type == 1) ? state->kPrime[c] : 2;

    for (int j = 0; j < zCols; ++j) zRow[j] = state->zMatrix(c, j);

    for (int j = 0; j < zCols; ++j) {
      double th = state->theta[j];
      double log_none = log_pi0;
      double log_enc  = log_1m_pi0 + std::log(th);
      double log_disc = log_1m_pi0 + std::log1p(-th);

      double ll[3];
      for (int v = 0; v < 3; ++v) {
        zRow[j] = v;
        ll[v] = per_char_log_lik_ecology(
          *data, c, state->parent, state->child, edgeLen,
          kp, state->rateLoss, state->rateNeo,
          rates, state->phi, zRow, wEdge,
          refE, gammaE, ws);
      }
      double lp[3] = { log_none, log_enc, log_disc };

      double logits[3];
      double mx = R_NegInf;
      for (int v = 0; v < 3; ++v) {
        logits[v] = lp[v] + beta * ll[v];
        if (logits[v] > mx) mx = logits[v];
      }
      if (!std::isfinite(mx)) {
        zRow[j] = state->zMatrix(c, j);
        continue;
      }
      double sumExp = 0.0;
      double pCum[3];
      for (int v = 0; v < 3; ++v) {
        pCum[v] = std::exp(logits[v] - mx);
        sumExp += pCum[v];
      }
      double u = R::unif_rand() * sumExp;
      int chosen = 2;
      double acc = 0.0;
      for (int v = 0; v < 3; ++v) {
        acc += pCum[v];
        if (u < acc) { chosen = v; break; }
      }
      zRow[j] = chosen;
      state->zMatrix(c, j) = chosen;
    }
  }

  // Recompute logLik / logPrior fresh to avoid drift.
  // T-010: refresh every partition into the cache; wEdge unchanged by z.
  // T-011: every char's per-edge mixture changed → every cached node CL stale.
  state->ecoCL.invalidate_all_cls();
  state->logLik = eco_recompute_all_partitions(data, state, edgeLen);
  state->logPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
  return true;
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
//           18=neo_joint_scale, 19=slice_scalar, 20=pspr,
//           21=joint_tl_rls, 22=joint_tl_rl,
//           23=dirichlet_branch, 24=local_dirichlet,
//           25=gibbs_kprime_sweep, 26=block_kprime_shift,
//           27=scale_kprime_alpha, 28=scale_kprime_beta,
//           29=slice_kprime_hyper,
//           30=mh_logit_p (logit-scale MH on p for empirical_geometric prior),
//           34=scale_phi (ecology),
//           35=scale_pi0 (ecology, logit-Bactrian),
//           36=gibbs_z_sweep (ecology),
//           37=scale_theta (ecology, logit-Bactrian)
//
// M-065: NNI/SPR now call _impl versions directly with parent/child vectors.
// Likelihood calls use vectors directly (no IntegerMatrix construction).
// ---------------------------------------------------------------------------

static bool do_move_impl(ChainRng& rng, McmcData* data, McmcState* state,
                         int moveType, int charIdx,
                         double scaleTuning, double betaSimplexTuning,
                         int intWalkWindow, double beta,
                         double jointRho = 0.0) {

  // (Periodic eco resync now lives in run_mcmc_batch_cpp, so it also
  // covers slice_scalar_impl moves which bypass do_move_impl.)

  // DIAG: pre-proposal LL consistency check (every 100 iterations).
  // In ecology mode we MUST use cpp_log_likelihood_ecology — the non-eco
  // path drops phi/z/pi0/gamma_e and would always report drift even when
  // the accumulator is correct.
  {
    static int preCheckCount = 0;
    if (++preCheckCount % 100 == 0) {
      int nE = state->relBrLengths.size();
      NumericVector curEl(nE);
      for (int i = 0; i < nE; ++i)
        curEl[i] = state->treeLength * state->relBrLengths[i];
      double freshLL;
      if (data->ecologyAware) {
        freshLL = cpp_log_likelihood_ecology(
          *data, state->parent, state->child, curEl,
          state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
          state->phi, state->zMatrix, state->pi0, state->theta);
      } else {
        freshLL = cpp_log_likelihood(*data, state->parent, state->child,
          curEl, state->kPrime, state->rateLoss, state->rateLogSd,
          state->rateNeo, state->betaScale,
          state->clWs.ready() ? &state->clWs : nullptr);
      }
      double drift = std::abs(state->logLik - freshLL);
      if (drift > 1e-4) {
        state->diagDriftCount++;
      }
    }
  }

  // Snapshot scalar state for rollback
  double oldTL   = state->treeLength;
  double oldRL   = state->rateLoss;
  double oldRLSD = state->rateLogSd;
  double oldRN   = state->rateNeo;
  double oldP    = state->p;
  double oldBS   = state->betaScale;  // M-052
  double oldKpA  = state->kprimeAlpha;
  double oldKpB  = state->kprimeBeta;

  double logHastings  = 0.0;
  bool topologyChanged = false;
  int oldKPrimeVal = 0;     // single-element rollback for case 7 (kPrime)
  int kPrimeCharIdx = -1;   // which character was changed

  // O(1) relBr rollback for cases 4 and 12 (save 2 modified elements)
  int bsIdx1 = -1, bsIdx2 = -1;
  double bsOldVal1 = 0.0, bsOldVal2 = 0.0;
  // M-158: nodeEdgeLen rollback for beta_simplex (2 values)
  double bsOldNodeEdgeLen1 = 0.0, bsOldNodeEdgeLen2 = 0.0;

  // Ecology phi rollback (case 30): single element modified.
  int phiOldIdx = -1;
  double phiOldVal = 0.0;

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
  if ((moveType == 5 || moveType == 4 || moveType == 23 || moveType == 24
       || moveType == 6) &&
      !data->qHeterogeneity && !data->ecologyAware && !state->nodeCL.ready()) {
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

  // T-013: pre-NNI population of the ecology partial-CL cache.  Disabled by
  // default — opt in via env MKPRIME_ECO_PARTIAL_CL=1.  On the rodent matrix
  // (64 tips, kEco=4) the dirty fraction per NNI swap is ~91 % median (see
  // dev/profiling/findings.md T-013) so partial-eval is a net loss after
  // dirty-walk + save/restore overhead.  May benefit larger trees.
  static const bool ecoPartialClEnabled = []{
    const char* env = std::getenv("MKPRIME_ECO_PARTIAL_CL");
    return env != nullptr && env[0] == '1';
  }();
  if (ecoPartialClEnabled && moveType == 5 && data->ecologyAware
      && !state->ecoCL.ready()) {
    int nEdge = state->relBrLengths.size();
    NumericVector absLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      absLen[i] = state->treeLength * state->relBrLengths[i];
    NumericMatrix wEdgeNow(nEdge, data->ecology.kEcology);
    recompute_w_edge(*data, state->parent, state->child, absLen, wEdgeNow);
    std::vector<double> gammaENow;
    compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                            gammaENow);
    populate_eco_cache_full(
      state->ecoCL, *data, state->parent, state->child, absLen,
      state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
      state->phi, state->zMatrix, state->pi0, state->theta,
      wEdgeNow, gammaENow);
  }

  switch (moveType) {
    case 0: { // scale tree_length (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->treeLength = oldTL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 1: { // scale rate_loss (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->rateLoss = oldRL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 2: { // scale rate_log_sd (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->rateLogSd = oldRLSD * mult;
      logHastings = std::log(mult);
      break;
    }
    case 3: { // scale rate_neo (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->rateNeo = oldRN * mult;
      logHastings = std::log(mult);
      break;
    }
    case 4: { // beta_simplex — O(1) rollback: save 2 modified elements
      int n_br = state->relBrLengths.size();
      int idx = static_cast<int>(rng.unif() * n_br);
      if (idx >= n_br) idx = n_br - 1;
      bsIdx1 = idx;
      if (!beta_simplex_impl(rng, state->relBrLengths, idx,
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

      int pick = (int)(rng.unif() * (double)intRows.size());
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

      int pV = (int)(rng.unif() * (double)vCh.size());
      if (pV >= (int)vCh.size()) pV = vCh.size() - 1;
      int pU = (int)(rng.unif() * (double)uSib.size());
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
      if (state->nodeCL.ready() && !data->qHeterogeneity && !data->ecologyAware) {
        // M-158: TreeNav-based SPR with partial CL evaluation
        sprMeta = propose_spr_treenav(state->nodeCL.topo);
        if (!sprMeta.valid || !R_FINITE(sprMeta.logHastings)) return false;
        logHastings = sprMeta.logHastings;
        sprPartialCL = true;
        // TreeNav is updated + partial eval happens in the evaluation section
      } else {
        // Fallback: full eval path (original OPP-6 pattern)
        List prop = spr_proposal_impl(rng, state->parent, state->child,
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
      List prop = tbr_proposal_impl(rng, state->parent, state->child,
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
      int delta  = static_cast<int>(rng.unif() * range) - intWalkWindow;
      int newK   = oldK + delta;
      if (newK < lowerK) return false;
      state->kPrime[charIdx] = newK;
      logHastings = 0.0;
      break;
    }
    case 8: { // scale p (legacy MH — kept for backward compat, not used by default)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
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
      int nTrans = (int)data->transIdxGlobal.size();
      if (nTrans == 0) return false;  // no transformational chars: skip
      double sumU = 0.0;
      for (int i = 0; i < nTrans; ++i) {
        int gi = data->transIdxGlobal[i];
        sumU += static_cast<double>(state->kPrime[gi] - data->kObs[gi]);
      }
      double shape1 = data->kprimeHyperA + nTrans;
      double shape2 = data->kprimeHyperB + sumU;
      state->p = rng.rbeta(shape1, shape2);
      // Recompute prior (p changed; likelihood unchanged)
      state->logPrior = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
      state->logLik = state->logLik;  // unchanged
      return true;  // Gibbs: always accept
    }
    case 10: { // gibbs_spr — M-085
      return gibbs_spr_impl(rng, data, state, beta);
    }
    case 11: { // gibbs_subtree_swap — M-086
      return gibbs_subtree_swap_impl(rng, data, state, beta);
    }
    case 12: { // weighted_branch_scale — M-087, O(1) rollback
      if (!weighted_branch_scale_impl(rng, data, state, beta, logHastings,
                                       bsIdx1, bsIdx2, bsOldVal1, bsOldVal2))
        return false;
      break;
    }
    case 13: { // weighted_spr — M-088
      return weighted_spr_impl(rng, data, state, beta);
    }
    case 14: { // weighted_subtree_swap — M-089
      return weighted_subtree_swap_impl(rng, data, state, beta);
    }
    case 15: { // block_gibbs_branch — M-054 reframed
      return block_gibbs_branch_sweep_impl(rng, data, state, beta);
    }
    case 16: { // M-052: scale beta_scale (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->betaScale = oldBS * mult;
      logHastings = std::log(mult);
      break;
    }
    case 18: { // neo_joint_scale (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->rateLoss = oldRL * mult;
      state->rateNeo  = oldRN * mult;
      logHastings = 2.0 * std::log(mult);
      break;
    }
    case 20: { // pSPR — M-119: parsimony-guided SPR
      List prop = pspr_proposal_impl(rng, state->parent, state->child,
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
      bactrian_2d_perturbation(rng, jointRho, z1, z2);
      double mult1 = std::exp(scaleTuning * z1);
      double mult2 = std::exp(scaleTuning * z2);
      state->treeLength = oldTL * mult1;
      state->rateLogSd  = oldRLSD * mult2;
      logHastings = std::log(mult1) + std::log(mult2);
      break;
    }
    case 22: { // M-120: joint_tl_rl (tree_length × rate_loss)
      double z1, z2;
      bactrian_2d_perturbation(rng, jointRho, z1, z2);
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
      if (!dirichlet_simplex_impl(rng, state->relBrLengths, nCats,
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
      if (!local_dirichlet_impl(rng, state->relBrLengths,
                                state->parent, state->child,
                                nCats, scaleTuning, logHastings,
                                state->brSnapshot,
                                state->dirEdges)) {
        return false;
      }
      break;
    }
    case 25: { // gibbs_kprime_sweep — Gibbs update of all k'_i
      return gibbs_kprime_sweep_impl(rng, data, state, beta);
    }
    case 26: { // block_kprime_shift — shift all k'_i by same delta
      return block_kprime_shift_impl(rng, data, state, intWalkWindow, beta);
    }
    case 27: { // scale_kprime_alpha — Bactrian scale for Beta-Geometric α
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->kprimeAlpha = oldKpA * mult;
      if (state->kprimeAlpha <= 0.0) return false;
      // Prior-only: likelihood is unchanged, compute prior ratio directly
      double newLP = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
        state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
      if (!R_FINITE(newLP)) { state->kprimeAlpha = oldKpA; return false; }
      double logAlpha = (newLP - state->logPrior) + std::log(mult);
      if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
        state->logPrior = newLP;
        return true;
      }
      state->kprimeAlpha = oldKpA;
      return false;
    }
    case 28: { // scale_kprime_beta — Bactrian scale for Beta-Geometric β
      double mult = std::exp(scaleTuning * bactrian_perturbation(rng));
      state->kprimeBeta = oldKpB * mult;
      if (state->kprimeBeta <= 0.0) return false;
      double newLP = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
        state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
      if (!R_FINITE(newLP)) { state->kprimeBeta = oldKpB; return false; }
      double logAlpha = (newLP - state->logPrior) + std::log(mult);
      if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
        state->logPrior = newLP;
        return true;
      }
      state->kprimeBeta = oldKpB;
      return false;
    }
    case 30: { // mh_logit_p — logit-scale MH on p (for empirical_geometric)
      // Multiplicative MH on p ∈ (0,1) overshoots when p is close to 1, which
      // is the typical posterior region under the empirical_geometric prior.
      // Propose on the unbounded logit scale instead, so no rejections from
      // boundary violations. Jacobian is |dp/dlogit(p)| = p (1 − p).
      if (oldP <= 0.0 || oldP >= 1.0) return false;
      double logitP = std::log(oldP / (1.0 - oldP));
      double logitPnew = logitP + scaleTuning * bactrian_perturbation(rng);
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
    case 34: { // scale phi (Bactrian) — ecology-aware NT
      if (!data->ecologyAware || state->phi.size() == 0) return false;
      int nPhi = state->phi.size();
      phiOldIdx = (nPhi == 1) ? 0
                : static_cast<int>(R::unif_rand() * nPhi);
      if (phiOldIdx >= nPhi) phiOldIdx = nPhi - 1;
      phiOldVal = state->phi[phiOldIdx];
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->phi[phiOldIdx] = phiOldVal * mult;
      logHastings = std::log(mult);
      break;
    }
    case 36: { // gibbs_z_sweep (ecology) — Gibbs over z_{c, s}
      return gibbs_z_sweep_impl(data, state, beta);
    }
    case 35: { // logit-Bactrian on pi0 (ecology) — v2: full MH with likelihood
      if (!data->ecologyAware) return false;
      double pi0Old = state->pi0;
      if (pi0Old <= 0.0 || pi0Old >= 1.0) return false;
      double logitOld = std::log(pi0Old / (1.0 - pi0Old));
      double logitNew = logitOld + scaleTuning * bactrian_perturbation();
      double pi0New;
      if (logitNew >= 0.0) {
        double e = std::exp(-logitNew);
        pi0New = 1.0 / (1.0 + e);
      } else {
        double e = std::exp(logitNew);
        pi0New = e / (1.0 + e);
      }
      if (!(pi0New > 0.0 && pi0New < 1.0)) return false;
      state->pi0 = pi0New;
      double logHast = std::log(pi0New * (1.0 - pi0New)) -
                       std::log(pi0Old * (1.0 - pi0Old));
      // T-010: pi0 affects gamma_e only (not wEdge); all partitions need a
      // fresh per-partition log-lik but the expensive wEdge sweep is reused.
      int nEdge = state->relBrLengths.size();
      NumericVector edgeLen(nEdge);
      for (int i = 0; i < nEdge; ++i)
        edgeLen[i] = state->treeLength * state->relBrLengths[i];
      eco_refresh_wedge(data, state, edgeLen);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      NumericVector rates = (state->rateLogSd > 0.0)
        ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
        : NumericVector(1, 1.0);
      int nP = (int)data->parts.size();
      std::vector<double> newPC(nP);
      double newLL = 0.0;
      for (int pi = 0; pi < nP; ++pi) {
        newPC[pi] = cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
        newLL += newPC[pi];
      }
      double newLP = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
        state->kprimeAlpha, state->kprimeBeta,
        &state->phi, state->pi0, &state->zMatrix, &state->theta);
      if (!R_FINITE(newLP) || !R_FINITE(newLL)) {
        state->pi0 = pi0Old;
        return false;
      }
      double logAlpha = (newLL - state->logLik) + (newLP - state->logPrior) + logHast;
      if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
        state->logLik = newLL;
        state->logPrior = newLP;
        state->partLogLik = std::move(newPC);
        state->ecoCL.invalidate_all_cls();  // T-011: pi0 changed factor table
        return true;
      }
      state->pi0 = pi0Old;
      return false;
    }
    case 37: { // logit-Bactrian on theta_e (ecology v2) — full MH
      if (!data->ecologyAware) return false;
      if (state->theta.size() == 0) return false;
      int nT = state->theta.size();
      int idx = (nT == 1) ? 0 : static_cast<int>(R::unif_rand() * nT);
      if (idx >= nT) idx = nT - 1;
      double thOld = state->theta[idx];
      if (thOld <= 0.0 || thOld >= 1.0) return false;
      double logitOld = std::log(thOld / (1.0 - thOld));
      double logitNew = logitOld + scaleTuning * bactrian_perturbation();
      double thNew;
      if (logitNew >= 0.0) {
        double e = std::exp(-logitNew);
        thNew = 1.0 / (1.0 + e);
      } else {
        double e = std::exp(logitNew);
        thNew = e / (1.0 + e);
      }
      if (!(thNew > 0.0 && thNew < 1.0)) return false;
      state->theta[idx] = thNew;
      double logHast = std::log(thNew * (1.0 - thNew)) -
                       std::log(thOld * (1.0 - thOld));
      // T-010: theta affects gamma_e only; wEdge cached.
      int nEdge = state->relBrLengths.size();
      NumericVector edgeLen(nEdge);
      for (int i = 0; i < nEdge; ++i)
        edgeLen[i] = state->treeLength * state->relBrLengths[i];
      eco_refresh_wedge(data, state, edgeLen);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      NumericVector rates = (state->rateLogSd > 0.0)
        ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
        : NumericVector(1, 1.0);
      int nP = (int)data->parts.size();
      std::vector<double> newPC(nP);
      double newLL = 0.0;
      for (int pi = 0; pi < nP; ++pi) {
        newPC[pi] = cpp_partition_log_likelihood_ecology(
          *data, pi, state->parent, state->child, edgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          state->wEdge, gammaE, rates);
        newLL += newPC[pi];
      }
      double newLP = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale,
        state->kprimeAlpha, state->kprimeBeta,
        &state->phi, state->pi0, &state->zMatrix, &state->theta);
      if (!R_FINITE(newLP) || !R_FINITE(newLL)) {
        state->theta[idx] = thOld;
        return false;
      }
      double logAlpha = (newLL - state->logLik) + (newLP - state->logPrior) + logHast;
      if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
        state->logLik = newLL;
        state->logPrior = newLP;
        state->partLogLik = std::move(newPC);
        state->ecoCL.invalidate_all_cls();  // T-011: theta changed factor table
        return true;
      }
      state->theta[idx] = thOld;
      return false;
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
    if (phiOldIdx >= 0) state->phi[phiOldIdx] = phiOldVal;
    return false;
  }

  // NNI doesn't change any parameter → prior is unchanged; skip evaluation.
  double newLogPrior;
  if (nniInPlace) {
    newLogPrior = state->logPrior;
  } else {
    newLogPrior = cpp_log_prior(
      *data, state->treeLength, state->relBrLengths,
      state->rateLoss, state->rateLogSd, state->rateNeo,
      state->p, state->kPrime, state->betaScale,
    state->kprimeAlpha, state->kprimeBeta,
    &state->phi, state->pi0, &state->zMatrix, &state->theta);
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
    if (phiOldIdx >= 0) state->phi[phiOldIdx] = phiOldVal;
    return false;
  }

  // OPP-6: select evaluation topology — proposed values for NNI/SPR,
  // current state for all other moves.
  const IntegerVector& evalParent = topologyChanged ? proposedParent : state->parent;
  const IntegerVector& evalChild  = topologyChanged ? proposedChild  : state->child;
  const NumericVector& evalRelBr  = topologyChanged ? proposedRelBr  : state->relBrLengths;

  // ---- Likelihood evaluation (M-064: partial, M-065: vectors, M-121: node CL) ----
  // Moves that only touch `p` (case 8 legacy multiplicative, case 30 logit MH)
  // leave the likelihood untouched.
  bool likChanges = (moveType != 8 && moveType != 30);
  // T-010: ecology mode now uses partition cache for non-tree moves. Tree /
  // branch moves invalidate wEdge; non-tree moves keep wEdge valid.
  const bool eco = data->ecologyAware;
  bool hasPLC = !eco && !state->partLogLik.empty();
  double newLogLik;
  std::vector<double> newPC;
  bool usedPartialCL = false;

  // Move types that change tree topology or any edge length (-> wEdge dirty
  // for the next eco call). Note: case 7 (int_walk kPrime), case 34 (phi)
  // do NOT touch the tree.
  auto move_invalidates_wedge = [](int mt) {
    switch (mt) {
      case 0:   // scale_tl
      case 4:   // beta_simplex
      case 5:   // NNI
      case 6:   // SPR
      case 10:  // gibbs_spr (returns false in eco — keep here for completeness)
      case 11:  // gibbs_subtree_swap (returns false in eco)
      case 12:  // weighted_branch_scale
      case 13:  // weighted_spr
      case 14:  // weighted_subtree_swap
      case 15:  // block_gibbs_branch
      case 17:  // TBR
      case 20:  // pSPR
      case 21:  // joint_tl_rls
      case 22:  // joint_tl_rl
      case 23:  // dirichlet_branch
      case 24:  // local_dirichlet
        return true;
      default:
        return false;
    }
  };

  // T-013: track whether the eco partial-CL path was used (for accept/reject)
  bool usedEcoPartialCL = false;

  if (!likChanges) {
    newLogLik = state->logLik;
  } else if (eco) {
    // T-010: ecology partition cache.
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    // T-013: NNI partial-eval path.  Try only when NNI in-place AND cache is
    // ready.  Falls through to T-010 full-eval logic on fallback (dirty set
    // too large) or when partial-CL is inapplicable.
    bool ecoPartialAccepted = false;
    if (nniInPlace && state->ecoCL.ready()) {
      // Mirror the M-158 NNI pattern: update cache's TreeNav in-place
      // before partial eval, revert on fallback / reject.
      update_topo_nni(state->ecoCL.topo, nniVNode, nniUNode,
                      nniCNode, nniWNode);

      NumericMatrix evalWEdge(nEdge, data->ecology.kEcology);
      recompute_w_edge(*data, evalParent, evalChild, propEdgeLen, evalWEdge);
      std::vector<double> gammaE;
      compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                              gammaE);
      EcoDirtyScratch scratch;
      int dirtyCount = 0;
      bool fallback = false;
      double partial_ll = eco_cache_partial_eval_nni(
        state->ecoCL, *data, evalParent, evalChild, propEdgeLen,
        state->kPrime, state->rateLoss, state->rateNeo,
        state->phi, state->zMatrix, state->pi0, state->theta,
        evalWEdge, gammaE,
        nniVNode, nniUNode,
        scratch, dirtyCount, fallback);

      if (!fallback) {
        newLogLik = partial_ll;
        ecoPartialAccepted = true;
        usedEcoPartialCL = true;
        // Stash rollback state on the cache.
        state->ecoCL.rollbackSavedCL    = std::move(scratch.savedCL);
        state->ecoCL.rollbackDirtyNodes = std::move(scratch.dirtyNodes);
        state->ecoCL.rollbackNniV = nniVNode;
        state->ecoCL.rollbackNniU = nniUNode;
        state->ecoCL.rollbackNniC = nniCNode;
        state->ecoCL.rollbackNniW = nniWNode;
        state->ecoCL.rollbackTopoUpdated = true;
      } else {
        // Dirty set too large.  Revert cache TreeNav; fall through.
        update_topo_nni(state->ecoCL.topo, nniVNode, nniUNode,
                        nniWNode, nniCNode);
      }
    }

    if (ecoPartialAccepted) {
      // Partial-eval path took over; skip the legacy logic.
      // partLogLik is now stale relative to newLogLik (partial-eval doesn't
      // maintain it); will be cleared on accept.
    } else {

    // Determine which partitions need recomputation and whether to refresh
    // wEdge from the proposed tree.
    bool wEdgeChanges = move_invalidates_wedge(moveType);
    bool hasPLCEco = !state->partLogLik.empty();

    // Build wEdge for the *evaluation* topology (which is either the proposed
    // tree for topology moves, or state->parent/child otherwise).
    NumericMatrix evalWEdge;
    if (wEdgeChanges || state->wEdgeDirty ||
        state->wEdge.nrow() != nEdge ||
        state->wEdge.ncol() != data->ecology.kEcology) {
      // Recompute fresh against the proposed topology / edge lengths.
      evalWEdge = NumericMatrix(nEdge, data->ecology.kEcology);
      recompute_w_edge(*data, evalParent, evalChild, propEdgeLen, evalWEdge);
    } else {
      evalWEdge = state->wEdge;
    }

    std::vector<double> gammaE;
    compute_gamma_e_ecology(*data, state->phi, state->pi0, state->theta,
                            gammaE);
    NumericVector rates = (state->rateLogSd > 0.0)
      ? cpp_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ)
      : NumericVector(1, 1.0);

    if (hasPLCEco && moveType == 7) {
      // kPrime int_walk: only the partition containing charIdx changes
      // (wEdge unchanged; subgroup composition handled by recomputing
      // the whole trans partition).
      int ap = data->charToPartition[charIdx];
      newPC = state->partLogLik;
      newLogLik = state->logLik;
      if (ap >= 0) {
        double v = cpp_partition_log_likelihood_ecology(
          *data, ap, evalParent, evalChild, propEdgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          evalWEdge, gammaE, rates);
        newLogLik += (v - newPC[ap]);
        newPC[ap] = v;
      }
    } else if (hasPLCEco && (moveType == 1 || moveType == 3 || moveType == 18)) {
      // rate_loss / rate_neo / neo_joint: neomorphic partitions only.
      newPC = state->partLogLik;
      newLogLik = state->logLik;
      for (size_t ni = 0; ni < data->neoPartIndices.size(); ++ni) {
        int pi = data->neoPartIndices[ni];
        double v = cpp_partition_log_likelihood_ecology(
          *data, pi, evalParent, evalChild, propEdgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          evalWEdge, gammaE, rates);
        newLogLik += (v - newPC[pi]);
        newPC[pi] = v;
      }
    } else {
      // Tree moves, phi (34), rateLogSd (2), and any cold-start: every
      // partition. wEdge has been built for the eval topology above.
      int nP = (int)data->parts.size();
      newPC.assign(nP, 0.0);
      newLogLik = 0.0;
      for (int pi = 0; pi < nP; ++pi) {
        newPC[pi] = cpp_partition_log_likelihood_ecology(
          *data, pi, evalParent, evalChild, propEdgeLen,
          state->kPrime, state->rateLoss, state->rateNeo,
          state->phi, state->zMatrix,
          evalWEdge, gammaE, rates);
        newLogLik += newPC[pi];
      }
    }
    }  // close T-013 ecoPartialAccepted-else
  } else if (!eco && nniInPlace && state->nodeCL.ready()) {
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
    // DIAG: compare NNI partial-CL with full eval
    { double fullLL = cpp_log_likelihood(*data, evalParent, evalChild,
        propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale,
        state->clWs.ready() ? &state->clWs : nullptr);
      double diff = std::abs(newLogLik - fullLL);
      state->diagNniPartialCount++;
      if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
      if (diff > 1e-6) state->diagNniMismatchCount++;
    }

  } else if (!eco && moveType == 4 && state->nodeCL.ready()) {
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
    // DIAG: compare BS partial-CL with full eval
    { double fullLL = cpp_log_likelihood(*data, evalParent, evalChild,
        propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale,
        state->clWs.ready() ? &state->clWs : nullptr);
      double diff = std::abs(newLogLik - fullLL);
      state->diagBsPartialCount++;
      if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
      if (diff > 1e-6) state->diagBsMismatchCount++;
    }

  } else if (!eco && (moveType == 23 || moveType == 24) && state->nodeCL.ready()) {
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
      newLogLik = cpp_log_likelihood(*data, evalParent, evalChild,
        propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale,
        state->clWs.ready() ? &state->clWs : nullptr);
    } else {
      newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                      evalParent, evalChild, propEdgeLen,
                                      state->rateLoss, state->rateNeo,
                                      state->rateLogSd, state->betaScale, dirty);
      usedPartialCL = true;

      // DIAG: compare Dirichlet partial-CL with full eval
      { double fullLL = cpp_log_likelihood(*data, evalParent, evalChild,
          propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
          state->rateNeo, state->betaScale,
          state->clWs.ready() ? &state->clWs : nullptr);
        double diff = std::abs(newLogLik - fullLL);
        state->diagDirPartialCount++;
        if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
        if (diff > 1e-6) {
          state->diagDirMismatchCount++;
          if (state->diagDirMismatchCount <= 5) {
            FILE* f2 = std::fopen("C:/Users/pjjg18/GitHub/mkp/pcl_diag.txt", "a");
            if (f2) {
              std::fprintf(f2, "  DIR_MM #%d: partial=%.6f full=%.6f diff=%.6f dirty=%d nEdges=%d mt=%d\n",
                      state->diagDirMismatchCount, newLogLik, fullLL, diff,
                      (int)dirty.size(), nEdge, moveType);
              std::fclose(f2);
            }
          }
        }
      }
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
      newLogLik = cpp_log_likelihood(*data, sprParent, sprChild,
        propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale,
        state->clWs.ready() ? &state->clWs : nullptr);
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
    if (eco) {
      newLogLik = cpp_log_likelihood_ecology(*data, evalParent, evalChild,
        propEdgeLen, state->kPrime,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->phi, state->zMatrix, state->pi0, state->theta);
    } else {
      newLogLik = cpp_log_likelihood(*data, evalParent, evalChild,
        propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
        state->rateNeo, state->betaScale,
        state->clWs.ready() ? &state->clWs : nullptr);
    }
  } else {
    int nParts = (int)data->parts.size();
    int nEdge  = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];
    newPC = state->partLogLik;
    switch (moveType) {
      case 1:
      case 3:
      case 18: {
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
  if (R_FINITE(logAlpha) && std::log(rng.unif()) < logAlpha) {
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

    // T-010: invalidate cached wEdge when the accepted move changed the
    // tree topology or any edge length. Use the same predicate as eval.
    if (eco) {
      switch (moveType) {
        case 0:  case 4:  case 5:  case 6:
        case 10: case 11: case 12: case 13: case 14: case 15:
        case 17: case 20: case 21: case 22: case 23: case 24:
          state->wEdgeDirty = true;
          break;
        default:
          break;
      }
    }

    // M-121/M-161: cache management on acceptance.
    // Partial CL moves keep the cache valid (already updated).
    // Other moves: granular invalidation by move type.
    // T-011: mirror the same granular invalidation on state->ecoCL so the
    // ecology partial-CL cache stays in sync.  When eco mode is off the
    // ecoCL is empty and the invalidate_* calls are O(1) no-ops.
    if (!usedPartialCL && likChanges) {
      switch (moveType) {
        case 1: case 3: case 18:
          // rate_loss, rate_neo, neo_joint: only neomorphic units affected
          state->nodeCL.invalidate_neo_cls();
          state->ecoCL.invalidate_neo_cls();
          break;
        case 2:
          // rateLogSd: ACRV rates change, all units stale
          state->nodeCL.invalidate_all_cls();
          state->ecoCL.invalidate_all_cls();
          break;
        case 7:
          // kPrime int_walk: unit structure may change
          state->nodeCL.invalidate_structure();
          state->ecoCL.invalidate_structure();
          break;
        case 16:
          // beta_scale: Q-het parameter, all units stale
          state->nodeCL.invalidate_all_cls();
          state->ecoCL.invalidate_all_cls();
          break;
        // T-011: ecology scalar params — wEdge unchanged, factor table dirty.
        case 34: case 35: case 37:  // scale_phi, scale_pi0, scale_theta
          state->ecoCL.invalidate_all_cls();
          break;
        default:
          // tree_length (0), topology, etc.: full invalidation
          state->nodeCL.invalidate_all();
          state->ecoCL.invalidate_all();
          break;
      }
    }
    // When partial CL was used, clear partition-level cache
    // (it's not maintained by partial eval; will be rebuilt if needed)
    if (usedPartialCL && !state->partLogLik.empty())
      state->partLogLik.clear();

    // T-013: eco partial-CL accept — cache already reflects the NEW topology.
    // partLogLik isn't maintained on this path; clear so future non-tree-move
    // evals fall through to the T-010 full-recompute path which rebuilds it.
    if (usedEcoPartialCL) {
      state->partLogLik.clear();
      state->wEdgeDirty = true;
      state->ecoCL.rollbackTopoUpdated = false;
    }

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
  // Ecology phi rollback
  if (phiOldIdx >= 0) state->phi[phiOldIdx] = phiOldVal;

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

  // T-013: rollback eco partial-CL on rejection.  Restore the saved CLs at
  // dirty nodes and reverse the NNI swap on the cache's TreeNav.
  if (usedEcoPartialCL) {
    EcoCLCache& c = state->ecoCL;
    if (!c.rollbackSavedCL.empty() && !c.rollbackDirtyNodes.empty()) {
      const double* src = c.rollbackSavedCL.data();
      for (auto& u : c.units) {
        for (int cat = 0; cat < u.nCat; ++cat) {
          for (int node : c.rollbackDirtyNodes) {
            double* dst = u.CL(cat, node);
            std::memcpy(dst, src, u.stride * sizeof(double));
            src += u.stride;
          }
        }
      }
    }
    if (c.rollbackTopoUpdated && c.rollbackNniV >= 0) {
      update_topo_nni(c.topo, c.rollbackNniV, c.rollbackNniU,
                      c.rollbackNniW, c.rollbackNniC);  // reverse swap
      c.rollbackTopoUpdated = false;
    }
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
  // T-017-IIa: R-exported wrapper uses one R-RNG draw to seed a ChainRng
  uint64_t seed = static_cast<uint64_t>(
    R::unif_rand() * static_cast<double>(std::numeric_limits<uint32_t>::max()));
  ChainRng rng(seed);
  return do_move_impl(rng, data, state, moveType, charIdx,
                      scaleTuning, betaSimplexTuning, intWalkWindow, beta);
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

  // T-014: configurable drift-resync gate frequency. Default 200 (was 20
  // hard-coded; that fired per-chain at 5 % of iterations, dominating
  // ~25 % of the aware-PT amplification on rodent). Lower (e.g. 1 or 5)
  // when debugging a new partial-CL path; higher (or absent) trusts the
  // cached path. The 5000-iter drift driver (`11g_t011_drift.R`) confirmed
  // 0 events at the previous /20 cadence so /200 is a safe production default.
  int resyncEvery = 200;
  const char* resyncEnv = std::getenv("MKPRIME_ECO_RESYNC_EVERY");
  if (resyncEnv != nullptr) {
    int v = std::atoi(resyncEnv);
    if (v >= 1) resyncEvery = v;
  }

  // T-015 diagnostic: per-chain timing of the PT loop. Off unless
  // MKPRIME_T015_DIAG=1. Emits per-chain move-dispatch / drift-resync
  // accumulators and one swap accumulator at end of batch, to stderr.
  // Numbers cross-validate the throughput-arithmetic claim that aware
  // nChains=4 PT runs at near-linear scaling (no super-linear
  // amplification) vs blind's super-linear cache bonus.
  bool t015Diag = false;
  {
    const char* env = std::getenv("MKPRIME_T015_DIAG");
    t015Diag = (env != nullptr && env[0] == '1');
  }
  std::vector<double> t015MoveNs(nChains, 0.0);
  std::vector<double> t015ResyncNs(nChains, 0.0);
  std::vector<long long> t015MoveCalls(nChains, 0);
  std::vector<long long> t015ResyncCalls(nChains, 0);
  double t015SwapNs = 0.0;
  long long t015SwapCalls = 0;

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

  // M-159: Cache-boosted weights (used when nodeCL.ready() && !qHeterogeneity).
  // Partial-CL-eligible moves {4=beta_simplex, 5=NNI, 6=SPR, 23=dirichlet,
  // 24=local_dirichlet} get cacheBonus multiplier.
  bool haveCacheBoost = cacheBonus > 1.0 && !data->qHeterogeneity &&
                        !data->ecologyAware;
  std::vector<double> cumWeightsCached(nMoves);
  double totalWeightCached = 0.0;
  if (haveCacheBoost) {
    for (int m = 0; m < nMoves; ++m) {
      int mt = moveTypeCodes[m];
      double w = moveWeights[m];
      if (mt == 4 || mt == 5 || mt == 6 || mt == 23 || mt == 24)
        w *= cacheBonus;
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
  bool includeEco = data->ecologyAware;
  int nPhiCols   = includeEco ? (int)states[0]->phi.size()   : 0;
  int nPi0Cols   = includeEco ? 1 : 0;
  int nThetaCols = includeEco ? (int)states[0]->theta.size() : 0;
  // Base columns: log_post, log_lik, tree_length, rate_log_sd (4).
  // rate_loss included only when hasNeo (like rate_neo, p, beta_scale).
  // +2 diagnostic columns: swap_cold (cold-chain swaps since last sample),
  // topo_hash (topology fingerprint for change detection).
  // +ecology: phi (nPhiCols), pi0 (1), theta (nThetaCols) when ecologyAware.
  int nScalarCols = 4 + (hasNeo ? 2 : 0) + nKpHyperCols +
                    (includeBS ? 1 : 0) + nPhiCols + nPi0Cols + nThetaCols +
                    2 + nTrans + nEdge;
  int maxSaved    = nBatch / thin + 2;
  std::vector<std::vector<double>> scalarRows;
  scalarRows.reserve(maxSaved);
  List edgeSamples;
  List zSamples;  // nChar x kEco IntegerMatrix snapshots when ecologyAware

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

  // T-017-IIa: Per-chain RNG instances (one per chain, indexed by ch).
  std::vector<ChainRng> rngs;
  rngs.reserve(nChains);

  // T-017-IIa: Per-chain RNG seed stream.
  // Draw exactly ONE R-RNG value to establish the base seed, then use
  // Weyl-mix to spread nChains seeds far apart in the mt19937_64 state space.
  // This is the ONLY R-RNG draw consumed inside run_mcmc_batch_cpp; all
  // other random draws within the per-chain loop go through ChainRng.
  {
    uint64_t base_seed = static_cast<uint64_t>(
      R::unif_rand() * static_cast<double>(std::numeric_limits<uint32_t>::max()));
    // Weyl-mix constant: 2^64 / phi (golden ratio). Spreads chain seeds
    // to avoid low-bit correlation between adjacent-index seeds.
    constexpr uint64_t WEYL_MIX = 0x9E3779B97F4A7C15ULL;
    for (int ch = 0; ch < nChains; ++ch) {
      uint64_t chain_seed = base_seed ^ (static_cast<uint64_t>(ch) * WEYL_MIX);
      rngs.emplace_back(chain_seed);
    }
  }

  // Main iteration loop
  for (int i = 0; i < nBatch; ++i) {
    // Check for user interrupt every 10 iterations (expensive moves can take
    // seconds each, so we want to stay responsive to Ctrl-C / ESC).
    if (i % 10 == 0) R_CheckUserInterrupt();

    int iter = startIter + i;

    // Advance each chain
    for (int ch = 0; ch < nChains; ++ch) {
      // M-159: Cache-aware weighted move selection.
      // When the node CL cache is valid, boost partial-CL-eligible moves.
      bool useBoost = haveCacheBoost && states[ch]->nodeCL.ready();
      double tw = useBoost ? totalWeightCached : totalWeight;
      const auto& cw = useBoost ? cumWeightsCached : cumWeights;
      if (useBoost) ++cacheHits; else ++cacheMisses;

      double u = rngs[ch].unif() * tw;
      int moveIdx = 0;
      while (moveIdx < nMoves - 1 && u > cw[moveIdx]) ++moveIdx;

      proposeCounts(ch, moveIdx)++;

      int moveType = moveTypeCodes[moveIdx];

      // charIdx: for int_walk → random trans character; for slice → paramIdx
      int charIdx = 0;
      if (moveType == 7 && nTrans > 0) {
        int r = static_cast<int>(rngs[ch].unif() * nTrans);
        if (r >= nTrans) r = nTrans - 1;
        charIdx = transIdxCpp[r];
      } else if (moveType == 19) {
        charIdx = sliceParamCodes[moveIdx];
      }

      auto t0 = std::chrono::steady_clock::now();
      auto t015_t0 = t0;
      bool accepted;
      if (moveType == 19) {
        // Slice sampling — self-contained, no MH accept/reject
        int nExp = 0;
        accepted = slice_scalar_impl(
          rngs[ch], data, states[ch], charIdx,
          sliceWidths(ch, moveIdx), betas[ch], 10, &nExp);
        sliceExpansions(ch, moveIdx) += nExp;
      } else if (moveType == 29) {
        // Prior-only slice sampler for BG hyperparameters (M-163)
        int nExp = 0;
        accepted = slice_kprime_hyper_impl(
          rngs[ch], data, states[ch], sliceParamCodes[moveIdx],
          sliceWidths(ch, moveIdx), 10, &nExp);
        sliceExpansions(ch, moveIdx) += nExp;
      } else {
        // Per-move int param overrides chain-level intWalkWindow
        int iww = moveIntParams[moveIdx] > 0
                    ? moveIntParams[moveIdx]
                    : chainIntWalkWins[ch];
        accepted = do_move_impl(
          rngs[ch], data, states[ch],
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
      if (t015Diag) {
        t015MoveNs[ch] +=
          (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
            t1 - t015_t0).count();
        t015MoveCalls[ch]++;
      }
      if (accepted) acceptCounts(ch, moveIdx)++;

      // Periodic from-scratch resync for ecology mode. Every K=20 iter
      // recompute state->logLik / state->logPrior from
      // cpp_log_likelihood_ecology + cpp_log_prior and overwrite the
      // accumulators. Guarantees chain dynamics see a correct baseline
      // regardless of any remaining incremental-update leaks. Logs any
      // drift > 0.5 nats with iter / moveType so leaks can be located.
      // Placed AFTER each move so it catches drift from every code
      // path (including slice_scalar_impl, which bypasses
      // do_move_impl).
      if (data->ecologyAware && (iter % resyncEvery == 0)) {
        auto t015_r0 = std::chrono::steady_clock::now();
        int nE = states[ch]->relBrLengths.size();
        NumericVector curEl(nE);
        for (int e = 0; e < nE; ++e)
          curEl[e] = states[ch]->treeLength * states[ch]->relBrLengths[e];
        double freshLL = cpp_log_likelihood_ecology(
          *data, states[ch]->parent, states[ch]->child, curEl,
          states[ch]->kPrime, states[ch]->rateLoss,
          states[ch]->rateLogSd, states[ch]->rateNeo,
          states[ch]->phi, states[ch]->zMatrix,
          states[ch]->pi0, states[ch]->theta);
        double freshLP = cpp_log_prior(
          *data, states[ch]->treeLength, states[ch]->relBrLengths,
          states[ch]->rateLoss, states[ch]->rateLogSd,
          states[ch]->rateNeo, states[ch]->p, states[ch]->kPrime,
          states[ch]->betaScale,
          states[ch]->kprimeAlpha, states[ch]->kprimeBeta,
          &states[ch]->phi, states[ch]->pi0,
          &states[ch]->zMatrix, &states[ch]->theta);
        double driftLL = freshLL - states[ch]->logLik;
        double driftLP = freshLP - states[ch]->logPrior;
        if (R_FINITE(freshLL) && R_FINITE(freshLP) &&
            std::abs(driftLL) + std::abs(driftLP) > 0.5) {
          REprintf("[eco-resync iter=%d ch=%d lastMt=%d] "
                   "dLL=%+.3f dLP=%+.3f\n",
                   iter, ch, moveType, driftLL, driftLP);
          states[ch]->diagDriftCount++;
        }
        if (R_FINITE(freshLL)) states[ch]->logLik   = freshLL;
        if (R_FINITE(freshLP)) states[ch]->logPrior = freshLP;
        if (t015Diag) {
          auto t015_r1 = std::chrono::steady_clock::now();
          t015ResyncNs[ch] +=
            (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
              t015_r1 - t015_r0).count();
          t015ResyncCalls[ch]++;
        }
      }
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
        auto t015_s0 = std::chrono::steady_clock::now();
        std::swap(*states[iPair], *states[jPair]);
        if (t015Diag) {
          auto t015_s1 = std::chrono::steady_clock::now();
          t015SwapNs +=
            (double)std::chrono::duration_cast<std::chrono::nanoseconds>(
              t015_s1 - t015_s0).count();
          t015SwapCalls++;
        }
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
      if (includeEco) {
        for (int e = 0; e < nPhiCols; ++e) row[col++] = s0->phi[e];
        row[col++] = s0->pi0;
        for (int e = 0; e < nThetaCols; ++e) row[col++] = s0->theta[e];
      }
      // Diagnostic: cold-chain swaps since last sample
      row[col++] = static_cast<double>(coldSwapsSinceSample);
      coldSwapsSinceSample = 0;
      // Diagnostic: topology hash (FNV-1a of canonical-preorder parent vector)
      row[col++] = fnv_topo_hash(s0->parent);
      for (int j = 0; j < nTrans; ++j)
        row[col++] = static_cast<double>(s0->kPrime[transIdxCpp[j]]);
      for (int k = 0; k < nEdge; ++k)
        row[col++] = s0->relBrLengths[k];
      scalarRows.push_back(row);

      // Edge matrix for tree reconstruction in R
      IntegerMatrix edgeMat(nEdge, 2);
      for (int k = 0; k < nEdge; ++k) {
        edgeMat(k, 0) = s0->parent[k];
        edgeMat(k, 1) = s0->child[k];
      }
      edgeSamples.push_back(edgeMat);

      // Ecology z snapshot for the per-(c, s) posterior.
      if (includeEco) {
        zSamples.push_back(clone(s0->zMatrix));
      }
    }
  }

  // Pack scalar rows into R matrix
  int nSaved = static_cast<int>(scalarRows.size());
  NumericMatrix scalarMat(nSaved, nScalarCols);
  for (int i = 0; i < nSaved; ++i)
    for (int j = 0; j < nScalarCols; ++j)
      scalarMat(i, j) = scalarRows[i][j];

  // T-015: emit per-chain timing breakdown when MKPRIME_T015_DIAG=1.
  if (t015Diag) {
    double tot_move = 0.0, tot_resync = 0.0;
    long long tot_move_calls = 0, tot_resync_calls = 0;
    for (int ch = 0; ch < nChains; ++ch) {
      tot_move        += t015MoveNs[ch];
      tot_resync      += t015ResyncNs[ch];
      tot_move_calls  += t015MoveCalls[ch];
      tot_resync_calls += t015ResyncCalls[ch];
    }
    REprintf("[T015 batch nBatch=%d nChains=%d] move_total=%.3f s "
             "resync_total=%.3f s swap_total=%.3f s swap_calls=%lld\n",
             nBatch, nChains, tot_move * 1e-9, tot_resync * 1e-9,
             t015SwapNs * 1e-9, t015SwapCalls);
    for (int ch = 0; ch < nChains; ++ch) {
      double per_move = (t015MoveCalls[ch] > 0)
        ? t015MoveNs[ch] / (double)t015MoveCalls[ch] : 0.0;
      double per_resync = (t015ResyncCalls[ch] > 0)
        ? t015ResyncNs[ch] / (double)t015ResyncCalls[ch] : 0.0;
      REprintf("[T015 ch=%d] move_calls=%lld move_mean_us=%.3f "
               "resync_calls=%lld resync_mean_us=%.3f\n",
               ch, t015MoveCalls[ch], per_move * 1e-3,
               t015ResyncCalls[ch], per_resync * 1e-3);
    }
  }

  // DIAG: aggregate counters from cold chain (index 0)
  IntegerVector diagCounters = IntegerVector::create(
    _["dir_partial"] = states[0]->diagDirPartialCount,
    _["dir_mismatch"] = states[0]->diagDirMismatchCount,
    _["dir_fullback"] = states[0]->diagDirFullbackCount,
    _["nni_partial"] = states[0]->diagNniPartialCount,
    _["bs_partial"] = states[0]->diagBsPartialCount,
    _["drift"] = states[0]->diagDriftCount
  );
  // DIAG: write summary to file (Rprintf is swallowed in RStudio batch loops)
  {
    FILE* f = std::fopen("pcl_diag.txt", "a");
    if (f) {
      std::fprintf(f, "cachePop=%d nni=%d(mm=%d) bs=%d(mm=%d) dir=%d(mm=%d,fb=%d) drift=%d maxD=%.2e\n",
              states[0]->diagCachePopCount,
              states[0]->diagNniPartialCount, states[0]->diagNniMismatchCount,
              states[0]->diagBsPartialCount, states[0]->diagBsMismatchCount,
              states[0]->diagDirPartialCount, states[0]->diagDirMismatchCount,
              states[0]->diagDirFullbackCount, states[0]->diagDriftCount,
              states[0]->diagMaxDiff);
      std::fclose(f);
    }
  }

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
    _["z_samples"]        = zSamples,
    _["n_saved"]          = nSaved,
    _["diag_counters"]    = diagCounters,
    _["diag_max_diff"]    = states[0]->diagMaxDiff,
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
  // Rebuild absLen (may have been modified)
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  TreeNav topo;
  topo.build(state->parent, state->child, absLen, nTip);

  NumericVector rates = gibbs_acrv_rates(state->rateLogSd, data->nCat, data->acrvZ);
  bool useAcrv = (state->rateLogSd > 0.0);
  int nCat = useAcrv ? data->nCat : 1;
  if (!useAcrv) rates = NumericVector(1, 1.0);

  int coding  = data->codingType;
  int maxNode = topo.maxNode;

  std::vector<CLGroup> groups;
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN = true; g.rateLoss = state->rateLoss; g.rateScale = state->rateNeo;
      g.tipData = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN = false; g.rateLoss = 1.0; g.rateScale = 1.0;
      g.tipData = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), part.k);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
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
        g.isMkN = false; g.rateLoss = 1.0; g.rateScale = 1.0;
        g.tipData = sub;
        g.allocate(maxNode, nCat, nSub, kp);
        groups.push_back(std::move(g));
        groupMeta.push_back({pi, nCharPart});
      }
    }
  }

  for (auto& grp : groups)
    caching_downpass(grp, topo, state->parent, state->child, rates);

  std::vector<CLGroup> pseudoGroups;
  if (coding != 0) {
    pseudoGroups.resize(groups.size());
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      pseudoGroups[gi] = create_const_pseudo_group(
        groups[gi], nTip, maxNode, nCat);
      caching_downpass(pseudoGroups[gi], topo, state->parent, state->child,
                       rates);
    }
  }

  int pA    = topo.parentNode[nodeA];
  int slotA = topo.childSlot(pA, nodeA);
  double lenA_val = topo.edgeLen[topo.edgeToPar[nodeA]];

  std::vector<int> pathA;
  pathA.reserve(16);
  for (int n = pA; n >= 1; n = topo.parentNode[n])
    pathA.push_back(n);
  std::vector<int> pathAIdx(maxNode + 1, -1);
  for (int i = 0; i < (int)pathA.size(); ++i)
    pathAIdx[pathA[i]] = i;

  NumericVector llPartial(nPart);
  for (int pi = 0; pi < nPart; ++pi) {
    double totalLL = 0.0;
    for (size_t gi = 0; gi < groups.size(); ++gi) {
      double grpLL = evaluate_swap_candidate(
        groups[gi], topo, rates, nodeA, partners[pi],
        pA, slotA, lenA_val, pathA, pathAIdx);
      if (coding != 0 && groups[gi].nChar > 0) {
        double constP = evaluate_swap_const_prob(
          pseudoGroups[gi], topo, rates, nodeA, partners[pi],
          pA, slotA, lenA_val, pathA, pathAIdx);
        if (constP < 1.0)
          grpLL -= groups[gi].nChar * std::log(1.0 - constP);
      }
      totalLL += grpLL;
    }

    if (data->relabel) {
      for (int pii = 0; pii < (int)data->parts.size(); ++pii) {
        const PartInfo& part = data->parts[pii];
        if (part.type == 1) {
          int nCharPart = part.tipStates.ncol();
          for (int ci = 0; ci < nCharPart; ++ci)
            totalLL += mk_prime_relabel_log(
              state->kPrime[part.globalCharIdx[ci]], part.kObsLocal[ci]);
        }
      }
    }

    llPartial[pi] = totalLL;
  }

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





