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
  IntegerVector kPrime;
  double logLik;
  double logPrior;
  // Per-partition log-likelihood cache (M-064)
  std::vector<double> partLogLik;
  // Pre-allocated CL workspace (M-063): eliminates per-call heap allocations
  ClWorkspace clWs;
  // M-052: Q-matrix heterogeneity — Dirichlet-marginal beta_scale parameter.
  double betaScale = 1.0;
  // M-121: persistent node-level CL cache for partial evaluation
  NodeCLCache nodeCL;
  // M-125: snapshot for block Dirichlet branch-length rollback
  NumericVector brSnapshot;
  // M-127: which edges the Dirichlet proposal modified (for partial CL eval)
  std::vector<int> dirEdges;
};


// ---------------------------------------------------------------------------
// Log prior (mirrors LogPrior in MkPrimeModel.R)
// ---------------------------------------------------------------------------

static double cpp_log_prior(
    const McmcData& data,
    double treeLength, const NumericVector& relBrLengths,
    double rateLoss, double rateLogSd, double rateNeo,
    double p, const IntegerVector& kPrime,
    double betaScale = 1.0) {

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
    // p boundary check only applies to hierarchical geometric
    if (!data.kPriorLogseries && (p <= 0.0 || p >= 1.0)) return R_NegInf;
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
// init_mcmc_state: create XPtr<McmcState>
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
SEXP init_mcmc_state(IntegerVector parent, IntegerVector child,
                     NumericVector relBrLengths, double treeLength,
                     double rateLoss, double rateLogSd, double rateNeo,
                     double p, IntegerVector kPrime,
                     double logLik, double logPrior,
                     double betaScale = 1.0) {
  McmcState* s = new McmcState();
  s->parent       = clone(parent);
  s->child        = clone(child);
  s->relBrLengths = clone(relBrLengths);
  s->treeLength   = treeLength;
  s->rateLoss     = rateLoss;
  s->rateLogSd    = rateLogSd;
  s->rateNeo      = rateNeo;
  s->p            = p;
  s->kPrime       = clone(kPrime);
  s->logLik       = logLik;
  s->logPrior     = logPrior;
  s->betaScale    = betaScale;
  s->brSnapshot   = NumericVector(relBrLengths.size());
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
    _["kPrime"]        = s->kPrime,
    _["logLik"]        = s->logLik,
    _["logPrior"]      = s->logPrior,
    _["logPost"]       = s->logLik + s->logPrior,
    _["betaScale"]     = s->betaScale
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

  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, propRelBr,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale);
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

  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, propRelBr,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale);
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
    state->p, state->kPrime, state->betaScale);
  if (!R_FINITE(logPrior)) return R_NegInf;

  double logLik;
  int nEdge = state->parent.size();
  NumericVector edgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    edgeLen[i] = state->treeLength * state->relBrLengths[i];

  bool hasPLC = !state->partLogLik.empty();
  ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

  // rate_loss (1), rate_neo (3): only neomorphic partitions change
  if (hasPLC && (paramIdx == 1 || paramIdx == 3)) {
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
      state->logPrior = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale);

      int nEdge = state->parent.size();
      NumericVector edgeLen(nEdge);
      for (int i = 0; i < nEdge; ++i)
        edgeLen[i] = state->treeLength * state->relBrLengths[i];
      ClWorkspace* wsPtr = state->clWs.ready() ? &state->clWs : nullptr;

      bool hasPLC = !state->partLogLik.empty();
      if (hasPLC && (paramIdx == 1 || paramIdx == 3)) {
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
// do_move_impl: internal propose/evaluate/accept (raw pointers, no SEXP).
// do_move_cpp:  Rcpp-exported SEXP wrapper — calls do_move_impl.
//
// moveType: 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
//           4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p (legacy),
//           9=gibbs_p, 10=gibbs_spr, 11=gibbs_subtree_swap,
//           12=weighted_br_scale, 13=weighted_spr, 14=weighted_subtree_swap,
//           15=block_gibbs_branch, 16=beta_scale, 17=tbr,
//           18=neo_joint_scale, 19=slice_scalar, 20=pspr,
//           21=joint_tl_rls, 22=joint_tl_rl
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

  double logHastings  = 0.0;
  bool topologyChanged = false;
  int oldKPrimeVal = 0;     // single-element rollback for case 7 (kPrime)
  int kPrimeCharIdx = -1;   // which character was changed

  // O(1) relBr rollback for cases 4 and 12 (save 2 modified elements)
  int bsIdx1 = -1, bsIdx2 = -1;
  double bsOldVal1 = 0.0, bsOldVal2 = 0.0;

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

  // M-121: pre-proposal cache population for partial CL.
  // Must happen BEFORE the proposal modifies state in-place.
  if ((moveType == 5 || moveType == 4) && !data->qHeterogeneity &&
      !state->nodeCL.valid) {
    int nEdge = state->relBrLengths.size();
    NumericVector absLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      absLen[i] = state->treeLength * state->relBrLengths[i];
    populate_cache_full(state->nodeCL, *data,
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
    case 6: { // SPR — OPP-6: no pre-proposal clone; defer state update to accept
      List prop = spr_proposal_impl(state->parent, state->child,
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
      state->logPrior = cpp_log_prior(
        *data, state->treeLength, state->relBrLengths,
        state->rateLoss, state->rateLogSd, state->rateNeo,
        state->p, state->kPrime, state->betaScale);
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
    case 18: { // neo_joint_scale (Bactrian, M-118)
      double mult = std::exp(scaleTuning * bactrian_perturbation());
      state->rateLoss = oldRL * mult;
      state->rateNeo  = oldRN * mult;
      logHastings = 2.0 * std::log(mult);
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
    default:
      return false;
  }

  if (!R_FINITE(logHastings)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS;
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
      state->p, state->kPrime, state->betaScale);
  }

  if (!R_FINITE(newLogPrior)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS;
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
    return false;
  }

  // OPP-6: select evaluation topology — proposed values for NNI/SPR,
  // current state for all other moves.
  const IntegerVector& evalParent = topologyChanged ? proposedParent : state->parent;
  const IntegerVector& evalChild  = topologyChanged ? proposedChild  : state->child;
  const NumericVector& evalRelBr  = topologyChanged ? proposedRelBr  : state->relBrLengths;

  // ---- Likelihood evaluation (M-064: partial, M-065: vectors, M-121: node CL) ----
  bool likChanges = (moveType != 8);
  bool hasPLC = !state->partLogLik.empty();
  double newLogLik;
  std::vector<double> newPC;
  bool usedPartialCL = false;

  if (!likChanges) {
    newLogLik = state->logLik;
  } else if (nniInPlace && state->nodeCL.valid) {
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

  } else if (moveType == 4 && state->nodeCL.valid) {
    // M-121: beta_simplex with valid node CL cache → partial evaluation
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    auto dirty = find_dirty_beta_simplex(state->nodeCL.topo,
                                          evalParent, bsIdx1, bsIdx2);
    newLogLik = partial_eval_dirty(state->nodeCL, *data,
                                    evalParent, evalChild, propEdgeLen,
                                    state->rateLoss, state->rateNeo,
                                    state->rateLogSd, state->betaScale, dirty);
    usedPartialCL = true;

  } else if ((moveType == 23 || moveType == 24) && state->nodeCL.valid) {
    // M-127: Dirichlet (random or local) with valid node CL cache → partial eval
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];

    auto dirty = find_dirty_dirichlet(state->nodeCL.topo,
                                       evalParent, state->dirEdges);
    // Heuristic: if dirty set covers most of the tree, fall back to full eval
    int nInternal = nEdge + 1 - data->nTip;
    if ((int)dirty.size() > (int)(0.8 * (nInternal + data->nTip))) {
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
    }

  } else if (!hasPLC) {
    int nEdge = evalRelBr.size();
    NumericVector propEdgeLen(nEdge);
    for (int i = 0; i < nEdge; ++i)
      propEdgeLen[i] = state->treeLength * evalRelBr[i];
    newLogLik = cpp_log_likelihood(*data, evalParent, evalChild,
      propEdgeLen, state->kPrime, state->rateLoss, state->rateLogSd,
      state->rateNeo, state->betaScale,
      state->clWs.ready() ? &state->clWs : nullptr);
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

    // M-121: cache management on acceptance.
    // Partial CL moves keep the cache valid (already updated).
    // All other moves that change CLs invalidate the cache.
    if (!usedPartialCL) {
      // Invalidate node CL cache for any accepted move that changes
      // topology, branch lengths, or model parameters.
      if (likChanges) state->nodeCL.valid = false;
    }
    // When partial CL was used, clear partition-level cache
    // (it's not maintained by partial eval; will be rebuilt if needed)
    if (usedPartialCL && !state->partLogLik.empty())
      state->partLogLik.clear();

    return true;
  }

  // Reject: rollback
  state->treeLength = oldTL;
  state->rateLoss   = oldRL;
  state->rateLogSd  = oldRLSD;
  state->rateNeo    = oldRN;
  state->p          = oldP;
  state->betaScale  = oldBS;
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

  // M-121: rollback node CL cache on rejection of partial-eval moves
  if (usedPartialCL) {
    restore_dirty_cls(state->nodeCL, state->nodeCL.dirtyNodes);
    // Also rollback TreeNav for NNI (topology was updated before partial eval)
    if (nniInPlace) {
      update_topo_nni(state->nodeCL.topo, nniVNode, nniUNode,
                      nniWNode, nniCNode);  // reverse swap
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
    int nEdge
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
  std::vector<double> cumWeights(nMoves);
  double totalWeight = 0.0;
  for (int m = 0; m < nMoves; ++m) {
    totalWeight += moveWeights[m];
    cumWeights[m] = totalWeight;
  }

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
  // p column is omitted when using the log-series prior (no hyperparameter)
  bool includeP  = !data->kPriorLogseries;
  bool includeBS = data->qHeterogeneity;  // M-052: beta_scale column
  // Base columns: log_post, log_lik, tree_length, rate_log_sd (4).
  // rate_loss included only when hasNeo (like rate_neo, p, beta_scale).
  // +2 diagnostic columns: swap_cold (cold-chain swaps since last sample),
  // topo_hash (topology fingerprint for change detection).
  int nScalarCols = 4 + (hasNeo ? 2 : 0) + (includeP ? 1 : 0) +
                    (includeBS ? 1 : 0) + 2 + nTrans + nEdge;
  int maxSaved    = nBatch / thin + 2;
  std::vector<std::vector<double>> scalarRows;
  scalarRows.reserve(maxSaved);
  List edgeSamples;

  // Diagnostic: count accepted swaps involving the cold chain (index 0)
  int coldSwapsSinceSample = 0;

  // Main iteration loop
  for (int i = 0; i < nBatch; ++i) {
    // Check for user interrupt every 10 iterations (expensive moves can take
    // seconds each, so we want to stay responsive to Ctrl-C / ESC).
    if (i % 10 == 0) R_CheckUserInterrupt();

    int iter = startIter + i;

    // Advance each chain
    for (int ch = 0; ch < nChains; ++ch) {
      // Weighted move selection
      double u = R::unif_rand() * totalWeight;
      int moveIdx = 0;
      while (moveIdx < nMoves - 1 && u > cumWeights[moveIdx]) ++moveIdx;

      proposeCounts(ch, moveIdx)++;

      int moveType = moveTypeCodes[moveIdx];

      // charIdx: for int_walk → random trans character; for slice → paramIdx
      int charIdx = 0;
      if (moveType == 7 && nTrans > 0) {
        int r = static_cast<int>(R::unif_rand() * nTrans);
        if (r >= nTrans) r = nTrans - 1;
        charIdx = transIdxCpp[r];
      } else if (moveType == 19) {
        charIdx = sliceParamCodes[moveIdx];
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
      if (includeP) row[col++] = s0->p;
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

  return List::create(
    _["accept_counts"]    = acceptCounts,
    _["propose_counts"]   = proposeCounts,
    _["move_time_ns"]     = moveTimeNs,
    _["slice_expansions"] = sliceExpansions,
    _["swap_accept"]      = swapAccept,
    _["swap_propose"]     = swapPropose,
    _["scalar_samples"]   = scalarMat,
    _["edge_samples"]     = edgeSamples,
    _["n_saved"]          = nSaved
  );
}

// [[Rcpp::export]]
List debug_mcmc_data(SEXP dataPtr) {
  McmcData* data = Rcpp::XPtr<McmcData>(dataPtr).get();
  return List::create(
    _["kPriorLogseries"]   = data->kPriorLogseries,
    _["kprimeLogseriesC"]  = data->kprimeLogseriesC,
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





