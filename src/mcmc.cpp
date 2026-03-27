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

using namespace Rcpp;

// Forward declarations for proposals in other TUs
// M-065: vector-based _impl versions (no edge matrix)
List nni_proposal_impl(IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
List spr_proposal_impl(IntegerVector parent, IntegerVector child,
                       int nTip, double treeLength,
                       NumericVector relBrLengths);
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
};


// ---------------------------------------------------------------------------
// Log prior (mirrors LogPrior in MkPrimeModel.R)
// ---------------------------------------------------------------------------

static double cpp_log_prior(
    const McmcData& data,
    double treeLength, const NumericVector& relBrLengths,
    double rateLoss, double rateLogSd, double rateNeo,
    double p, const IntegerVector& kPrime) {

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
    if (p <= 0.0 || p >= 1.0) return R_NegInf;
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
    double sumU = 0.0;
    for (int i = 0; i < nTrans; ++i) {
      int gi = data.transIdxGlobal[i];
      sumU += (kPrime[gi] - data.kObs[gi]);
    }
    lp += nTrans * std::log(p) + sumU * std::log1p(-p);
    lp += R::dbeta(p, data.kprimeHyperA, data.kprimeHyperB, 1);
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
                     double logLik, double logPrior) {
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
  for (int pi = 0; pi < nParts; ++pi)
    state->partLogLik[pi] = cpp_partition_log_likelihood(
      *data, pi, state->parent, state->child, edgeLen,
      state->kPrime, state->rateLoss, state->rateLogSd, state->rateNeo,
      state->clWs.ready() ? &state->clWs : nullptr);
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
    _["logPost"]       = s->logLik + s->logPrior
  );
}


// [[Rcpp::export]]
double get_state_log_lik(SEXP statePtr) {
  return Rcpp::XPtr<McmcState>(statePtr).get()->logLik;
}


// ---------------------------------------------------------------------------
// do_move_impl: internal propose/evaluate/accept (raw pointers, no SEXP).
// do_move_cpp:  Rcpp-exported SEXP wrapper — calls do_move_impl.
//
// moveType: 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
//           4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p
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
    case 8: { // scale p
      double mult = std::exp(scaleTuning * (R::unif_rand() - 0.5));
      state->p = oldP * mult;
      logHastings = std::log(mult);
      break;
    }
    default:
      return false;
  }

  if (!R_FINITE(logHastings)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    if (moveType == 4) state->relBrLengths = oldRelBr;
    if (moveType == 7) state->kPrime = oldKPrime;
    return false;
  }

  // Evaluate prior (relBrLengths prior = lgamma(n), value-independent)
  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime);

  if (!R_FINITE(newLogPrior)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    if (moveType == 4) state->relBrLengths = oldRelBr;
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
      state->rateNeo, state->clWs.ready() ? &state->clWs : nullptr);
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
            wsPtr);
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
            wsPtr);
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
            wsPtr);
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
  if (moveType == 4) state->relBrLengths = oldRelBr;
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
  int nSwapPairs = std::max(0, nChains - 1);
  IntegerVector swapAccept(nSwapPairs, 0);
  IntegerVector swapPropose(nSwapPairs, 0);

  // Sample storage
  int nScalarCols = 6 + (hasNeo ? 1 : 0) + nTrans + nEdge;
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

      bool accepted = do_move_impl(
        data, states[ch],
        moveType, charIdx,
        chainScaleTunings(ch, moveIdx),
        chainBsmpTunings[ch],
        chainIntWalkWins[ch],
        betas[ch]
      );
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
      row[col++] = s0->p;
      if (hasNeo) row[col++] = s0->rateNeo;
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
    _["swap_accept"]    = swapAccept,
    _["swap_propose"]   = swapPropose,
    _["scalar_samples"] = scalarMat,
    _["edge_samples"]   = edgeSamples,
    _["n_saved"]        = nSaved
  );
}

