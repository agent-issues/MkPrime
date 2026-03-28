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
#include <TreeTools/renumber_tree.h>
#include <cmath>
#include <chrono>

using namespace Rcpp;

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
List beta_simplex_proposal(NumericVector x, int index, double tuning);
bool beta_simplex_impl(NumericVector& x, int index, double tuning,
                       double& logHastings);  // OPP-5: in-place, no List alloc


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

  if (treeLength <= 0.0) return R_NegInf;
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

  // nNode = max 1-indexed node in tree
  int maxNode = 0;
  for (int i = 0; i < state->parent.size(); ++i) {
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


// [[Rcpp::export]]
double get_state_log_lik(SEXP statePtr) {
  return Rcpp::XPtr<McmcState>(statePtr).get()->logLik;
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
// gibbs_spr_impl  (M-085)
//
// GibbsSPR: enumerate all valid SPR reattachment positions for a randomly
// chosen subtree, weight each by exp(β × logLik) (prior is topology-
// independent), sample one proportionally (including the current state
// to ensure proper discrete Gibbs semantics), and apply.
//
// Design choices:
//   tau = 0.5 (midpoint insertion) → no branch-length Jacobian.
//   Current state included in candidate set → pure Gibbs, acceptance = 1
//     on non-self draws; self-draw recorded as a no-op (returns false).
//   Dirichlet(1,...,1) prior on relBrLengths is uniform and equal for all
//     candidates → only likelihoods enter the weights.
// ---------------------------------------------------------------------------
static bool gibbs_spr_impl(McmcData* data, McmcState* state, double beta) {
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
  const int u = state->parent[pruneRow];   // pruned subtree's parent node
  const int v = state->child[pruneRow];    // pruned subtree root

  // 3. Find parentRow (edge into u) and sibRow (u's other child)
  int parentRow = -1, sibRow = -1, sibNode = -1;
  for (int i = 0; i < nEdge; ++i) {
    if (state->child[i] == u) parentRow = i;
    if (state->parent[i] == u && state->child[i] != v) {
      sibRow = i; sibNode = state->child[i];
    }
  }
  if (parentRow < 0 || sibRow < 0) return false;

  // 4. BFS: mark all descendants of v
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

  // 6. Absolute edge lengths (treeLength is preserved across topology change)
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];
  const double lMerge = absLen[parentRow] + absLen[sibRow];

  // 7. Compute logLik for each candidate at tau = 0.5
  std::vector<double> candLL(nCand);
  for (int ci = 0; ci < nCand; ++ci) {
    const int rr      = cands[ci];
    const double lReg = absLen[rr];

    IntegerVector np = clone(state->parent);
    IntegerVector nc = clone(state->child);
    NumericVector na = clone(absLen);
    const int b = nc[rr];

    nc[parentRow] = sibNode;  na[parentRow] = lMerge;
    nc[rr]        = u;        na[rr]        = 0.5 * lReg;
    np[sibRow]    = u;        nc[sibRow]    = b;
    na[sibRow]    = 0.5 * lReg;

    auto po = TreeTools::preorder_weighted_impl(np, nc, na);
    IntegerVector op = po.first(_, 0);
    IntegerVector oc = po.first(_, 1);
    candLL[ci] = compute_full_loglik_at(*data, *state, op, oc, po.second);
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

  // 9. Sample: self-draw → no-op (recorded as rejection)
  double rnd = R::unif_rand() * sumW;
  if (rnd < wOrig) return false;
  rnd -= wOrig;
  int chosen = nCand - 1;
  for (int ci = 0; ci < nCand - 1; ++ci) {
    if (rnd < ws[ci]) { chosen = ci; break; }
    rnd -= ws[ci];
  }

  // 10. Apply chosen SPR: reconstruct and commit topology to state
  {
    const int rr      = cands[chosen];
    const double lReg = absLen[rr];
    IntegerVector np  = clone(state->parent);
    IntegerVector nc  = clone(state->child);
    NumericVector na  = clone(absLen);
    const int b = nc[rr];

    nc[parentRow] = sibNode;  na[parentRow] = lMerge;
    nc[rr]        = u;        na[rr]        = 0.5 * lReg;
    np[sibRow]    = u;        nc[sibRow]    = b;
    na[sibRow]    = 0.5 * lReg;

    auto po = TreeTools::preorder_weighted_impl(np, nc, na);
    IntegerMatrix ordEdge = po.first;
    NumericVector ordAbs  = po.second;
    for (int k = 0; k < nEdge; ++k) {
      state->parent[k]       = ordEdge(k, 0);
      state->child[k]        = ordEdge(k, 1);
      state->relBrLengths[k] = ordAbs[k] / state->treeLength;
    }
  }

  // 11. Commit: logLik updated; prior/logPrior unchanged (topology-
  //     independent); clear per-partition cache (stale after topology change)
  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  return true;
}


// ---------------------------------------------------------------------------
// gibbs_subtree_swap_impl  (M-086)
//
// GibbsSubtreeSwap: enumerate all valid subtree-swap partners for a randomly
// chosen node, weight by exp(β × logLik), sample proportionally, apply.
// Same Gibbs semantics and design choices as gibbs_spr_impl.
// Branch lengths swap with their subtrees (Jacobian = 1); see M-084.
// ---------------------------------------------------------------------------
static bool gibbs_subtree_swap_impl(McmcData* data, McmcState* state,
                                    double beta) {
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

  // 3. Compute logLik for each candidate swap
  std::vector<double> candLL(nPart);
  for (int pi = 0; pi < nPart; ++pi) {
    List sw = swap_subtrees_impl(clone(state->parent), clone(state->child),
                                 nTip, state->treeLength,
                                 clone(state->relBrLengths),
                                 nodeA, partners[pi]);
    if (as<double>(sw["logHastings"]) == R_NegInf) {
      candLL[pi] = R_NegInf;
      continue;
    }
    IntegerVector sp = sw["parent"];
    IntegerVector sc = sw["child"];
    NumericVector sr = sw["rel_br_lengths"];
    NumericVector sa(nEdge);
    for (int k = 0; k < nEdge; ++k) sa[k] = state->treeLength * sr[k];
    candLL[pi] = compute_full_loglik_at(*data, *state, sp, sc, sa);
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

  // 6. Apply chosen swap
  List sw = swap_subtrees_impl(clone(state->parent), clone(state->child),
                               nTip, state->treeLength,
                               clone(state->relBrLengths),
                               nodeA, partners[chosen]);
  state->parent       = as<IntegerVector>(sw["parent"]);
  state->child        = as<IntegerVector>(sw["child"]);
  state->relBrLengths = as<NumericVector>(sw["rel_br_lengths"]);

  state->logLik = candLL[chosen];
  state->partLogLik.clear();
  return true;
}


// ---------------------------------------------------------------------------
// Bin structure for WeightedBranchLengthScale (M-087) — Beta(0.25, 0.25)
// quantile breakpoints, lazily initialised on first call.
// ---------------------------------------------------------------------------

struct BranchBins {
  int nBins = 0;
  std::vector<double> breaks;  // size nBins+1: [0]=0, [nBins]=1
  std::vector<double> mids;    // size nBins: midpoint of each bin
};

static BranchBins s_branchBins;

static const BranchBins& get_branch_bins(int nBins) {
  if (s_branchBins.nBins == nBins) return s_branchBins;
  s_branchBins.nBins = nBins;
  s_branchBins.breaks.resize(nBins + 1);
  s_branchBins.mids.resize(nBins);
  s_branchBins.breaks[0] = 0.0;
  s_branchBins.breaks[nBins] = 1.0;
  for (int b = 1; b < nBins; ++b)
    s_branchBins.breaks[b] = R::qbeta((double)b / nBins, 0.25, 0.25, 1, 0);
  for (int b = 0; b < nBins; ++b)
    s_branchBins.mids[b] =
      0.5 * (s_branchBins.breaks[b] + s_branchBins.breaks[b + 1]);
  return s_branchBins;
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
// The bin-selection weights are symmetric (forward == reverse) because
// midpoint evaluations are independent of the current fraction.
// logHastings = log(w_oldBin) + logBeta(f_old|a_old,b_old)
//             - log(w_chosenBin) - logBeta(f_new|a_new,b_new).
// ---------------------------------------------------------------------------
static bool weighted_branch_scale_impl(
    McmcData* data, McmcState* state, double beta, int nBins,
    double& logHastings) {

  const int nEdge = state->relBrLengths.size();
  if (nEdge < 2) return false;

  const BranchBins& bins = get_branch_bins(nBins);

  // 1. Pick two branches (same scheme as beta_simplex_impl)
  int index = static_cast<int>(R::unif_rand() * nEdge);
  if (index >= nEdge) index = nEdge - 1;
  int other = static_cast<int>(R::unif_rand() * (nEdge - 1));
  if (other >= index) ++other;
  if (other >= nEdge) other = nEdge - 1;
  if (other == index) other = (index + 1) % nEdge;

  const double oldRelA = state->relBrLengths[index];
  const double oldRelB = state->relBrLengths[other];
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
  //    Concentration = 2 * nBins gives moderate spread matching bin width.
  const double conc = 2.0 * nBins;
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
                               double beta, int nBins) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;
  const int root  = nTip + 1;

  const BranchBins& bins = get_branch_bins(nBins);

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

  // 8. Candidate marginals: for each candidate SPR, evaluate at each bin
  //    midpoint varying the branch fraction at the regraft point.
  //    candLL[ci][b] stores the log-likelihood.
  std::vector<std::vector<double>> candLL(nCand,
                                           std::vector<double>(nBins));
  std::vector<double> candMax(nCand, R_NegInf);

  for (int ci = 0; ci < nCand; ++ci) {
    const int rr      = cands[ci];
    const double lReg = absLen[rr];

    // Construct SPR topology once (same as GibbsSPR)
    IntegerVector np = clone(state->parent);
    IntegerVector nc = clone(state->child);
    const int b_node = nc[rr];
    nc[parentRow] = sibNode;
    nc[rr]        = u;
    np[sibRow]    = u;
    nc[sibRow]    = b_node;

    // Evaluate at each bin midpoint
    for (int b = 0; b < nBins; ++b) {
      NumericVector na = clone(absLen);
      na[parentRow] = lMerge;
      na[rr]        = bins.mids[b] * lReg;
      na[sibRow]    = (1.0 - bins.mids[b]) * lReg;

      auto po = TreeTools::preorder_weighted_impl(np, nc, na);
      IntegerVector op = po.first(_, 0);
      IntegerVector oc = po.first(_, 1);
      candLL[ci][b] = compute_full_loglik_at(*data, *state,
                                              op, oc, po.second);
      if (R_FINITE(candLL[ci][b]) && candLL[ci][b] > candMax[ci])
        candMax[ci] = candLL[ci][b];
    }
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
  const double conc = 2.0 * nBins;
  const double chosenMid = bins.mids[chosenBin];
  const double alphaNew = chosenMid * conc + 1.0;
  const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
  double fNew = R::rbeta(alphaNew, betaNew);
  if (fNew < 1e-8) fNew = 1e-8;
  if (fNew > 1.0 - 1e-8) fNew = 1.0 - 1e-8;

  // 13. Construct final proposed topology with fNew
  const int rr      = cands[chosen];
  const double lReg = absLen[rr];
  IntegerVector np  = clone(state->parent);
  IntegerVector nc  = clone(state->child);
  NumericVector na  = clone(absLen);
  const int b_node  = nc[rr];
  nc[parentRow] = sibNode;  na[parentRow] = lMerge;
  nc[rr]        = u;        na[rr]        = fNew * lReg;
  np[sibRow]    = u;        nc[sibRow]    = b_node;
  na[sibRow]    = (1.0 - fNew) * lReg;

  auto po = TreeTools::preorder_weighted_impl(np, nc, na);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;

  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbs);
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
    propRelBr[k] = ordAbs[k] / state->treeLength;

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

// Local helper: find edge row where child[i] == node (mirrors tree_moves.cpp)
static int find_child_row_local(const IntegerVector& child, int node) {
  for (int i = 0; i < child.size(); ++i)
    if (child[i] == node) return i;
  return -1;
}

static bool weighted_subtree_swap_impl(McmcData* data, McmcState* state,
                                        double beta, int nBins) {
  const int nEdge = state->parent.size();
  const int nTip  = data->nTip;

  const BranchBins& bins = get_branch_bins(nBins);

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
  const int rowA = find_child_row_local(state->child, nodeA);
  if (rowA < 0) return false;

  // 4. Absolute edge lengths
  NumericVector absLen(nEdge);
  for (int i = 0; i < nEdge; ++i)
    absLen[i] = state->treeLength * state->relBrLengths[i];

  // 5. Candidate marginals: for each partner B_i, construct swapped
  //    topology and evaluate at each bin midpoint
  std::vector<std::vector<double>> candLL(nPart,
                                           std::vector<double>(nBins));
  std::vector<double> candMax(nPart, R_NegInf);
  std::vector<int> rowBs(nPart);       // edge row for each partner
  std::vector<double> totals(nPart);   // brA + brB_i

  for (int pi = 0; pi < nPart; ++pi) {
    int rowB = find_child_row_local(state->child, partners[pi]);
    if (rowB < 0) { candMax[pi] = R_NegInf; rowBs[pi] = -1; continue; }
    rowBs[pi]  = rowB;
    totals[pi] = absLen[rowA] + absLen[rowB];

    // Construct swapped topology: swap parent assignments
    IntegerVector np = clone(state->parent);
    np[rowA] = state->parent[rowB];
    np[rowB] = state->parent[rowA];
    // child vector unchanged

    for (int b = 0; b < nBins; ++b) {
      NumericVector na = clone(absLen);
      na[rowA] = bins.mids[b] * totals[pi];
      na[rowB] = (1.0 - bins.mids[b]) * totals[pi];

      auto po = TreeTools::preorder_weighted_impl(np, state->child, na);
      IntegerVector op = po.first(_, 0);
      IntegerVector oc = po.first(_, 1);
      candLL[pi][b] = compute_full_loglik_at(*data, *state,
                                              op, oc, po.second);
      if (R_FINITE(candLL[pi][b]) && candLL[pi][b] > candMax[pi])
        candMax[pi] = candLL[pi][b];
    }
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
  const double conc = 2.0 * nBins;
  const double chosenMid = bins.mids[chosenBin];
  const double alphaNew = chosenMid * conc + 1.0;
  const double betaNew  = (1.0 - chosenMid) * conc + 1.0;
  double fNew = R::rbeta(alphaNew, betaNew);
  if (fNew < 1e-8) fNew = 1e-8;
  if (fNew > 1.0 - 1e-8) fNew = 1.0 - 1e-8;

  // 10. Construct final proposed topology at fNew
  const int rowB   = rowBs[chosen];
  const double tot = totals[chosen];
  IntegerVector np = clone(state->parent);
  np[rowA] = state->parent[rowB];
  np[rowB] = state->parent[rowA];
  NumericVector na = clone(absLen);
  na[rowA] = fNew * tot;
  na[rowB] = (1.0 - fNew) * tot;

  auto po = TreeTools::preorder_weighted_impl(np, state->child, na);
  IntegerMatrix ordEdge = po.first;
  NumericVector ordAbs  = po.second;
  IntegerVector op = ordEdge(_, 0);
  IntegerVector oc = ordEdge(_, 1);
  double newLogLik = compute_full_loglik_at(*data, *state, op, oc, ordAbs);
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
    propRelBr[k] = ordAbs[k] / state->treeLength;

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
// do_move_impl: internal propose/evaluate/accept (raw pointers, no SEXP).
// do_move_cpp:  Rcpp-exported SEXP wrapper — calls do_move_impl.
//
// moveType: 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
//           4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p (legacy),
//           9=gibbs_p, 10=gibbs_spr, 11=gibbs_subtree_swap,
//           12=weighted_br_scale, 13=weighted_spr, 14=weighted_subtree_swap
//
// M-065: NNI/SPR now call _impl versions directly with parent/child vectors.
// Likelihood calls use vectors directly (no IntegerMatrix construction).
// ---------------------------------------------------------------------------

static bool do_move_impl(McmcData* data, McmcState* state,
                         int moveType, int charIdx,
                         double scaleTuning, double betaSimplexTuning,
                         int intWalkWindow, double beta) {

  // Snapshot scalar state for rollback
  double oldTL   = state->treeLength;
  double oldRL   = state->rateLoss;
  double oldRLSD = state->rateLogSd;
  double oldRN   = state->rateNeo;
  double oldP    = state->p;
  double oldBS   = state->betaScale;  // M-052

  double logHastings  = 0.0;
  bool topologyChanged = false;
  NumericVector oldRelBr;   // rollback for case 4 only
  IntegerVector oldKPrime;

  // OPP-6: proposed topology held separately; state->parent/child not
  // overwritten until acceptance → no pre-proposal clone, no rollback copy.
  IntegerVector proposedParent, proposedChild;
  NumericVector proposedRelBr;

  switch (moveType) {
    case 0: { // scale tree_length
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->treeLength = oldTL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 1: { // scale rate_loss
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->rateLoss = oldRL * mult;
      logHastings = std::log(mult);
      break;
    }
    case 2: { // scale rate_log_sd
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->rateLogSd = oldRLSD * mult;
      logHastings = std::log(mult);
      break;
    }
    case 3: { // scale rate_neo
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->rateNeo = oldRN * mult;
      logHastings = std::log(mult);
      break;
    }
    case 4: { // beta_simplex — OPP-5: call impl directly, no Rcpp::List alloc
      oldRelBr = clone(state->relBrLengths);
      int n_br = state->relBrLengths.size();
      int idx = static_cast<int>(R::unif_rand() * n_br);
      if (idx >= n_br) idx = n_br - 1;
      if (!beta_simplex_impl(state->relBrLengths, idx,
                             betaSimplexTuning, logHastings))
        return false;
      break;
    }
    case 5: { // NNI — OPP-6: no pre-proposal clone; defer state update to accept
      List prop = nni_proposal_impl(state->parent, state->child,
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
    case 7: { // int_walk kPrime
      oldKPrime = clone(state->kPrime);
      int oldK   = state->kPrime[charIdx];
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
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
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
    case 12: { // weighted_branch_scale — M-087
      oldRelBr = clone(state->relBrLengths);
      if (!weighted_branch_scale_impl(data, state, beta,
                                       data->nBranchBins, logHastings))
        return false;
      break;
    }
    case 13: { // weighted_spr — M-088
      return weighted_spr_impl(data, state, beta, data->nBranchBins);
    }
    case 14: { // weighted_subtree_swap — M-089
      return weighted_subtree_swap_impl(data, state, beta, data->nBranchBins);
    }
    case 16: { // M-052: scale beta_scale
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->betaScale = oldBS * mult;
      logHastings = std::log(mult);
      break;
    }
    default:
      return false;
  }

  if (!R_FINITE(logHastings)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS;
    if (moveType == 4 || moveType == 12) state->relBrLengths = oldRelBr;
    if (moveType == 7) state->kPrime = oldKPrime;
    return false;
  }

  // Evaluate prior (relBrLengths prior = lgamma(n), value-independent)
  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime, state->betaScale);

  if (!R_FINITE(newLogPrior)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    state->betaScale = oldBS;
    if (moveType == 4 || moveType == 12) state->relBrLengths = oldRelBr;
    if (moveType == 7) state->kPrime = oldKPrime;
    return false;
  }

  // OPP-6: select evaluation topology — proposed values for NNI/SPR,
  // current state for all other moves.
  const IntegerVector& evalParent = topologyChanged ? proposedParent : state->parent;
  const IntegerVector& evalChild  = topologyChanged ? proposedChild  : state->child;
  const NumericVector& evalRelBr  = topologyChanged ? proposedRelBr  : state->relBrLengths;

  // ---- Likelihood evaluation (M-064: partial, M-065: vectors) ----
  bool likChanges = (moveType != 8);
  bool hasPLC = !state->partLogLik.empty();
  double newLogLik;
  std::vector<double> newPC;

  if (!likChanges) {
    newLogLik = state->logLik;
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
      case 3: {
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
    return true;
  }

  // Reject: rollback scalar state only
  // OPP-6: state->parent/child were never overwritten — no topology rollback needed
  state->treeLength = oldTL;
  state->rateLoss   = oldRL;
  state->rateLogSd  = oldRLSD;
  state->rateNeo    = oldRN;
  state->p          = oldP;
  state->betaScale  = oldBS;
  if (moveType == 4 || moveType == 12) state->relBrLengths = oldRelBr;
  if (moveType == 7) state->kPrime = oldKPrime;
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
    NumericVector moveWeights,
    NumericMatrix chainScaleTunings,
    NumericVector chainBsmpTunings,
    IntegerVector chainIntWalkWins,
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
  int nSwapPairs = std::max(0, nChains - 1);
  IntegerVector swapAccept(nSwapPairs, 0);
  IntegerVector swapPropose(nSwapPairs, 0);

  // Sample storage
  // p column is omitted when using the log-series prior (no hyperparameter)
  bool includeP  = !data->kPriorLogseries;
  bool includeBS = data->qHeterogeneity;  // M-052: beta_scale column
  int nScalarCols = 5 + (includeP ? 1 : 0) + (hasNeo ? 1 : 0) +
                    (includeBS ? 1 : 0) + nTrans + nEdge;
  int maxSaved    = nBatch / thin + 2;
  std::vector<std::vector<double>> scalarRows;
  scalarRows.reserve(maxSaved);
  List edgeSamples;

  // Main iteration loop
  for (int i = 0; i < nBatch; ++i) {
    int iter = startIter + i;

    // Advance each chain
    for (int ch = 0; ch < nChains; ++ch) {
      // Weighted move selection
      double u = R::unif_rand() * totalWeight;
      int moveIdx = 0;
      while (moveIdx < nMoves - 1 && u > cumWeights[moveIdx]) ++moveIdx;

      proposeCounts(ch, moveIdx)++;

      int moveType = moveTypeCodes[moveIdx];

      // charIdx for int_walk (kPrime) move
      int charIdx = 0;
      if (moveType == 7 && nTrans > 0) {
        int r = static_cast<int>(R::unif_rand() * nTrans);
        if (r >= nTrans) r = nTrans - 1;
        charIdx = transIdxCpp[r];
      }

      auto t0 = std::chrono::steady_clock::now();
      bool accepted = do_move_impl(
        data, states[ch],
        moveType, charIdx,
        chainScaleTunings(ch, moveIdx),
        chainBsmpTunings[ch],
        chainIntWalkWins[ch],
        betas[ch]
      );
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
      row[col++] = s0->rateLoss;
      row[col++] = s0->rateLogSd;
      if (includeP) row[col++] = s0->p;
      if (hasNeo) row[col++] = s0->rateNeo;
      if (includeBS) row[col++] = s0->betaScale;  // M-052
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
    _["accept_counts"]  = acceptCounts,
    _["propose_counts"] = proposeCounts,
    _["move_time_ns"]   = moveTimeNs,
    _["swap_accept"]    = swapAccept,
    _["swap_propose"]   = swapPropose,
    _["scalar_samples"] = scalarMat,
    _["edge_samples"]   = edgeSamples,
    _["n_saved"]        = nSaved
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





