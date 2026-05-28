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
#include <TreeTools/renumber_tree.h>
#include <cmath>
#include <cstring>
#include <chrono>
#include <cstdio>

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
  int diagDirMismatchCount = 0;
  int diagDirFullbackCount = 0;
  int diagNniPartialCount = 0;
  int diagNniMismatchCount = 0;
  int diagBsPartialCount = 0;
  int diagBsMismatchCount = 0;
  int diagDriftCount = 0;
  int diagCachePopCount = 0;
  double diagMaxDiff = 0.0;

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
    const NumericVector& edgeLen) {
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
static bool gibbs_spr_impl_full(McmcData* data, McmcState* state, double beta);
// M-114: partial CL path for Q-heterogeneity
static bool gibbs_spr_impl_het(McmcData* data, McmcState* state, double beta);

static bool gibbs_spr_impl(McmcData* data, McmcState* state, double beta) {
  // M-114: Q-heterogeneity uses streaming partial CL
  if (data->qHeterogeneity)
    return gibbs_spr_impl_het(data, state, beta);

  // LIKE-001 interim (option 3): the pseudo-character partial-CL path below
  // computes only the constant-site ascertainment term, not the singleton
  // term. Under coding="informative" (codingType == 2) that omission
  // biases the Gibbs sampling weights. Fall back to the full evaluator,
  // which routes through cpp_partition_log_likelihood with the singleton
  // correction applied. Restore partial-CL once evaluate_singleton_prob
  // (math-prover option 2) is implemented.
  if (data->codingType == 2)
    return gibbs_spr_impl_full(data, state, beta);

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
  int pickIdx = (int)(R::unif_rand() * (double)eligible.size());
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

  // Audit Issue 1: RB-style partition-rate normalisation (nChar-weighted
  // mean rate = 1). Both neo and trans/known unit rateScales are derived
  // from rate_neo + per-partition character counts; legacy "rateScale = 1
  // for trans" behaviour is recovered when nNeo == 0 or nTrans == 0.
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

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
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

    } else if (part.type == 2) {
      // Known state space: fixed k
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
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
        g.rateScale = pScales.trans;
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

      // Ascertainment correction via pseudo-character partial CLs.
      // Note: only the constant-site term is computed here; coding == 2
      // (informative) additionally requires a singleton-site term that
      // this partial-CL path does not yet evaluate. The caller short-
      // circuits to gibbs_spr_impl_full when codingType == 2 (LIKE-001
      // interim), so we only reach this branch under coding == 1.
      if (coding != 0 && groups[gi].nChar > 0) {
        double constP = evaluate_const_prob(
          pseudoGroups[gi], topo, pseudoResiduals[gi], rates,
          v, u, sibNode, lMerge, a, b, lHalf, lPrune);
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
  double rnd = R::unif_rand() * sumW;
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
  return true;
}


// ---------------------------------------------------------------------------
// M-114: Gibbs SPR with streaming partial CL for Q-heterogeneity.
//
// Same prune/candidate/sampling logic as gibbs_spr_impl, but evaluates
// each CLGroup's likelihood under a mixture of F81 components by streaming
// over (betaBin, rotation) and accumulating per-site raw likelihoods.
// ---------------------------------------------------------------------------
static bool gibbs_spr_impl_het(McmcData* data, McmcState* state,
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

  int pickIdx = (int)(R::unif_rand() * (double)eligible.size());
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

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

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
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
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
        g.rateScale = pScales.trans;
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

  double rnd = R::unif_rand() * sumW;
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
  return true;
}


// Old full-evaluation fallback (Q-het or validation), M-109 in-place
static bool gibbs_spr_impl_full(McmcData* data, McmcState* state, double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  std::vector<int> eligible;
  eligible.reserve(nEdge);
  for (int i = 0; i < nEdge; ++i)
    if (state->parent[i] != root) eligible.push_back(i);
  if (eligible.empty()) return false;

  int pickIdx = (int)(R::unif_rand() * (double)eligible.size());
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

  double rnd = R::unif_rand() * sumW;
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
static bool gibbs_subtree_swap_impl_full(McmcData* data, McmcState* state,
                                         double beta);
// M-114: partial CL path for Q-heterogeneity
static bool gibbs_subtree_swap_impl_het(McmcData* data, McmcState* state,
                                         double beta);

static bool gibbs_subtree_swap_impl(McmcData* data, McmcState* state,
                                    double beta) {
  // M-114: Q-heterogeneity uses streaming partial CL
  if (data->qHeterogeneity)
    return gibbs_subtree_swap_impl_het(data, state, beta);

  // LIKE-001 interim (option 3): see comment in gibbs_spr_impl. The pseudo-
  // character partial-CL path omits the singleton-site ascertainment term
  // required under coding="informative"; fall back to the full evaluator.
  if (data->codingType == 2)
    return gibbs_subtree_swap_impl_full(data, state, beta);

  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  // 1. Pick a random node (any edge child is a valid non-root candidate)
  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
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

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

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
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});

    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
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
        g.rateScale = pScales.trans;
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
  double rnd = R::unif_rand() * sumW;
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
  return true;
}


// ---------------------------------------------------------------------------
// M-114: Gibbs subtree swap with streaming partial CL for Q-heterogeneity.
// ---------------------------------------------------------------------------
static bool gibbs_subtree_swap_impl_het(McmcData* data, McmcState* state,
                                         double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
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

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

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
      g.rateScale = pScales.neo;
      g.tipData   = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN     = false;
      g.rateLoss  = 1.0;
      g.rateScale = pScales.trans;
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
        g.rateScale = pScales.trans;
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

  double rnd = R::unif_rand() * sumW;
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
  return true;
}


// Full-evaluation fallback for Q-heterogeneity (M-109 in-place pattern)
static bool gibbs_subtree_swap_impl_full(McmcData* data, McmcState* state,
                                         double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
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

  double rnd = R::unif_rand() * sumW;
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
// Hastings ratio: topology selection cancels (Z = Z' by symmetry of the
// candidate set).  Remaining branch-fraction component:
//   logHR = log(w_{self,b_old}) + logBeta(f_old | b_old)
//         - log(w_{chosen,b_new}) - logBeta(f_new | b_new)
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
// for the chosen partner.  MH acceptance corrects the approximation.
//
// Cost: O(N × B) likelihood evaluations.
// ---------------------------------------------------------------------------

// (Uses find_child_row_gibbs defined above)

static bool weighted_subtree_swap_impl(McmcData* data, McmcState* state,
                                        double beta) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  const BranchBins& bins = data->branchBins;
  const int nBins = bins.nBins;

  // 1. Pick a random node (any edge child)
  int pickIdx = (int)(R::unif_rand() * (double)nEdge);
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
  double rnd = R::unif_rand() * sumM;
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
    double rndBin = R::unif_rand() * mCand[chosen];
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
  double fNew = R::rbeta(alphaNew, betaNew);
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

  double newLogPrior = compute_log_prior_at(*data, *state, propRelBr);
  if (!R_FINITE(newLogPrior)) return false;

  // 13. MH acceptance
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

  bool hasPLC = !state->partLogLik.empty();
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

      bool hasPLC = !state->partLogLik.empty();
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
// Prior-only slice sampler for Beta-Geometric hyperparameters (M-163)
//
// Samples kprimeAlpha (paramCode=0) or kprimeBeta (paramCode=1) using a
// univariate slice sampler on the log scale.  Target = logPrior + log(x)
// (the log(x) Jacobian arises from sampling u = log(x) and transforming).
// No likelihood evaluation needed — α/β only affect the prior.
// ---------------------------------------------------------------------------
static bool slice_kprime_hyper_impl(McmcData* data, McmcState* state,
                                     int paramCode, double width,
                                     int maxSteps = 10,
                                     int* nExpansionsOut = nullptr) {
  double x0 = (paramCode == 0) ? state->kprimeAlpha : state->kprimeBeta;
  if (x0 <= 0.0) return false;

  // Work on log scale: u = log(x)
  double u0 = std::log(x0);

  // Target: logPrior(current) + u (Jacobian)
  double logY0 = state->logPrior + u0;
  double logZ = logY0 + std::log(R::unif_rand());

  // Helper lambda to evaluate target at a candidate u
  auto evalTarget = [&](double u) -> double {
    double xCand = std::exp(u);
    if (xCand <= 0.0 || !R_FINITE(xCand)) return R_NegInf;
    double oldVal = (paramCode == 0) ? state->kprimeAlpha : state->kprimeBeta;
    if (paramCode == 0) state->kprimeAlpha = xCand;
    else                state->kprimeBeta  = xCand;
    double lp = compute_log_prior(*data, *state);
    // Restore
    if (paramCode == 0) state->kprimeAlpha = oldVal;
    else                state->kprimeBeta  = oldVal;
    return lp + u;  // logPrior + Jacobian
  };

  // Stepping out
  int nExp = 0;
  double L = u0 - width * R::unif_rand();
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
    double u1 = L + R::unif_rand() * (R_bound - L);
    double logTarget1 = evalTarget(u1);
    if (logTarget1 >= logZ) {
      // Accept
      double x1 = std::exp(u1);
      if (paramCode == 0) state->kprimeAlpha = x1;
      else                state->kprimeBeta  = x1;
      // Recompute and cache logPrior
      state->logPrior = compute_log_prior(*data, *state);
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
// Gibbs kPrime sweep (moveType 25)
//
// Samples each k'_i from its full conditional in a single random-order scan
// of all transformational characters. Always accepts (Gibbs update).
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
    double* logW = &charLogW[ti * K_MAX_CAND];

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
  if ((moveType == 5 || moveType == 4 || moveType == 23 || moveType == 24
       || moveType == 6) &&
      !data->qHeterogeneity && !state->nodeCL.ready()) {
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
    case 25: { // gibbs_kprime_sweep — Gibbs update of all k'_i
      return gibbs_kprime_sweep_impl(data, state, beta);
    }
    case 26: { // block_kprime_shift — shift all k'_i by same delta
      return block_kprime_shift_impl(data, state, intWalkWindow, beta);
    }
    case 27: { // scale_kprime_alpha — Bactrian scale for Beta-Geometric α
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->kprimeAlpha = oldKpA * mult;
      if (state->kprimeAlpha <= 0.0) return false;
      // Prior-only: likelihood is unchanged, compute prior ratio directly
      double newLP = compute_log_prior(*data, *state);
      if (!R_FINITE(newLP)) { state->kprimeAlpha = oldKpA; return false; }
      double logAlpha = (newLP - state->logPrior) + std::log(mult);
      if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
        state->logPrior = newLP;
        return true;
      }
      state->kprimeAlpha = oldKpA;
      return false;
    }
    case 28: { // scale_kprime_beta — Bactrian scale for Beta-Geometric β
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->kprimeBeta = oldKpB * mult;
      if (state->kprimeBeta <= 0.0) return false;
      double newLP = compute_log_prior(*data, *state);
      if (!R_FINITE(newLP)) { state->kprimeBeta = oldKpB; return false; }
      double logAlpha = (newLP - state->logPrior) + std::log(mult);
      if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
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
  // leave the likelihood untouched.
  bool likChanges = (moveType != 8 && moveType != 30);
  bool hasPLC = !state->partLogLik.empty();
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
    // DIAG: compare NNI partial-CL with full eval
    { ClWorkspace* nniWs = state->clWs.ready() ? &state->clWs : nullptr;
      double fullLL = state->usePartitioned
        ? cpp_log_likelihood_partitioned(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss,
            state->classRateLogSd, state->classRate,
            state->etaNeo, state->betaScale, nniWs)
        : cpp_log_likelihood(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd,
            state->rateNeo, state->betaScale, nniWs);
      double diff = std::abs(newLogLik - fullLL);
      state->diagNniPartialCount++;
      if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
      if (diff > 1e-6) state->diagNniMismatchCount++;
    }

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
    // DIAG: compare BS partial-CL with full eval
    { ClWorkspace* bsWs = state->clWs.ready() ? &state->clWs : nullptr;
      double fullLL = state->usePartitioned
        ? cpp_log_likelihood_partitioned(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss,
            state->classRateLogSd, state->classRate,
            state->etaNeo, state->betaScale, bsWs)
        : cpp_log_likelihood(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd,
            state->rateNeo, state->betaScale, bsWs);
      double diff = std::abs(newLogLik - fullLL);
      state->diagBsPartialCount++;
      if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
      if (diff > 1e-6) state->diagBsMismatchCount++;
    }

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

      // DIAG: compare Dirichlet partial-CL with full eval
      { ClWorkspace* dWs2 = state->clWs.ready() ? &state->clWs : nullptr;
        double fullLL = state->usePartitioned
          ? cpp_log_likelihood_partitioned(
              *data, evalParent, evalChild, propEdgeLen,
              state->kPrime, state->rateLoss,
              state->classRateLogSd, state->classRate,
              state->etaNeo, state->betaScale, dWs2)
          : cpp_log_likelihood(
              *data, evalParent, evalChild, propEdgeLen,
              state->kPrime, state->rateLoss, state->rateLogSd,
              state->rateNeo, state->betaScale, dWs2);
        double diff = std::abs(newLogLik - fullLL);
        state->diagDirPartialCount++;
        if (diff > state->diagMaxDiff) state->diagMaxDiff = diff;
        if (diff > 1e-6) {
          state->diagDirMismatchCount++;
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
    { ClWorkspace* pWs = state->clWs.ready() ? &state->clWs : nullptr;
      newLogLik = state->usePartitioned
        ? cpp_log_likelihood_partitioned(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss,
            state->classRateLogSd, state->classRate,
            state->etaNeo, state->betaScale, pWs)
        : cpp_log_likelihood(
            *data, evalParent, evalChild, propEdgeLen,
            state->kPrime, state->rateLoss, state->rateLogSd,
            state->rateNeo, state->betaScale, pWs);
    }
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

  // M-159: Cache-boosted weights (used when nodeCL.ready() && !qHeterogeneity).
  // Partial-CL-eligible moves {4=beta_simplex, 5=NNI, 6=SPR, 23=dirichlet,
  // 24=local_dirichlet} get cacheBonus multiplier.
  bool haveCacheBoost = cacheBonus > 1.0 && !data->qHeterogeneity;
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
      // M-159: Cache-aware weighted move selection.
      // When the node CL cache is valid, boost partial-CL-eligible moves.
      bool useBoost = haveCacheBoost && states[ch]->nodeCL.ready();
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
      // Diagnostic: topology hash (FNV-1a of canonical-preorder parent vector)
      row[col++] = fnv_topo_hash(s0->parent);
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
    _["dir_mismatch"] = states[0]->diagDirMismatchCount,
    _["dir_fullback"] = states[0]->diagDirFullbackCount,
    _["nni_partial"] = states[0]->diagNniPartialCount,
    _["bs_partial"] = states[0]->diagBsPartialCount,
    _["drift"] = states[0]->diagDriftCount
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

  // Audit Issue 1: partition-rate normalisation (see comment at first call site).
  const PartitionScales pScales =
      compute_partition_scales(state->rateNeo, data->nNeo, data->nTrans);

  std::vector<CLGroup> groups;
  struct GroupMeta { int partIdx; int nCharInPart; };
  std::vector<GroupMeta> groupMeta;

  for (int pi = 0; pi < (int)data->parts.size(); ++pi) {
    const PartInfo& part = data->parts[pi];
    if (part.type == 0) {
      CLGroup g;
      g.isMkN = true; g.rateLoss = state->rateLoss; g.rateScale = pScales.neo;
      g.tipData = part.tipStates;
      g.allocate(maxNode, nCat, part.tipStates.ncol(), 2);
      groups.push_back(std::move(g));
      groupMeta.push_back({pi, part.tipStates.ncol()});
    } else if (part.type == 2) {
      CLGroup g;
      g.isMkN = false; g.rateLoss = 1.0; g.rateScale = pScales.trans;
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
        g.isMkN = false; g.rateLoss = 1.0; g.rateScale = pScales.trans;
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





