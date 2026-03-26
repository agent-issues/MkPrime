// C++ MCMC engine — state management, prior, and propose/accept cycle.
//
// McmcState (mutable per chain) holds tree + parameters.
// do_move_cpp() performs a full MH step, modifying state in-place.
// R creates states via init_mcmc_state(), reads via get_mcmc_state().

#include "mcmc_state.h"
#include <TreeTools/renumber_tree.h>
#include <cmath>

using namespace Rcpp;

// Forward declarations for proposals in other TUs
List nni_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths);
List spr_proposal(IntegerMatrix edge, int nTip, double treeLength,
                  NumericVector relBrLengths);
List beta_simplex_proposal(NumericVector x, int index, double tuning);


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

  // Tree length: Gamma(shape, rate) — R::dgamma uses scale = 1/rate
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

  double logHastings = 0.0;
  bool topologyChanged = false;
  IntegerVector oldParent, oldChild;
  NumericVector oldRelBr;
  IntegerVector oldKPrime;

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
    case 4: { // beta_simplex
      oldRelBr = clone(state->relBrLengths);
      int n_br = state->relBrLengths.size();
      int idx = static_cast<int>(R::unif_rand() * n_br);
      if (idx >= n_br) idx = n_br - 1;
      List prop = beta_simplex_proposal(state->relBrLengths, idx,
                                        betaSimplexTuning);
      state->relBrLengths = as<NumericVector>(prop["value"]);
      logHastings = as<double>(prop["logHastings"]);
      break;
    }
    case 5: { // NNI
      oldParent = clone(state->parent);
      oldChild  = clone(state->child);
      oldRelBr  = clone(state->relBrLengths);
      IntegerMatrix edge(state->parent.size(), 2);
      for (int i = 0; i < state->parent.size(); ++i) {
        edge(i, 0) = state->parent[i];
        edge(i, 1) = state->child[i];
      }
      List prop = nni_proposal(edge, data->nTip, state->treeLength,
                               state->relBrLengths);
      logHastings = as<double>(prop["logHastings"]);
      if (!R_FINITE(logHastings)) return false;
      IntegerMatrix newEdge = as<IntegerMatrix>(prop["edge"]);
      for (int i = 0; i < state->parent.size(); ++i) {
        state->parent[i] = newEdge(i, 0);
        state->child[i]  = newEdge(i, 1);
      }
      state->relBrLengths = as<NumericVector>(prop["rel_br_lengths"]);
      topologyChanged = true;
      break;
    }
    case 6: { // SPR
      oldParent = clone(state->parent);
      oldChild  = clone(state->child);
      oldRelBr  = clone(state->relBrLengths);
      IntegerMatrix edge(state->parent.size(), 2);
      for (int i = 0; i < state->parent.size(); ++i) {
        edge(i, 0) = state->parent[i];
        edge(i, 1) = state->child[i];
      }
      List prop = spr_proposal(edge, data->nTip, state->treeLength,
                               state->relBrLengths);
      logHastings = as<double>(prop["logHastings"]);
      if (!R_FINITE(logHastings)) return false;
      IntegerMatrix newEdge = as<IntegerMatrix>(prop["edge"]);
      for (int i = 0; i < state->parent.size(); ++i) {
        state->parent[i] = newEdge(i, 0);
        state->child[i]  = newEdge(i, 1);
      }
      state->relBrLengths = as<NumericVector>(prop["rel_br_lengths"]);
      topologyChanged = true;
      break;
    }
    case 7: { // int_walk kPrime
      oldKPrime = clone(state->kPrime);
      int oldK = state->kPrime[charIdx];
      int lowerK = data->kObs[charIdx];
      int range = 2 * intWalkWindow + 1;
      int delta = static_cast<int>(R::unif_rand() * range) - intWalkWindow;
      int newK = oldK + delta;
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
    if (topologyChanged) {
      state->parent = oldParent; state->child = oldChild;
      state->relBrLengths = oldRelBr;
    }
    if (moveType == 4) state->relBrLengths = oldRelBr;
    if (moveType == 7) state->kPrime = oldKPrime;
    return false;
  }

  // Evaluate prior
  double newLogPrior = cpp_log_prior(
    *data, state->treeLength, state->relBrLengths,
    state->rateLoss, state->rateLogSd, state->rateNeo,
    state->p, state->kPrime);

  if (!R_FINITE(newLogPrior)) {
    state->treeLength = oldTL; state->rateLoss = oldRL;
    state->rateLogSd = oldRLSD; state->rateNeo = oldRN; state->p = oldP;
    if (topologyChanged) {
      state->parent = oldParent; state->child = oldChild;
      state->relBrLengths = oldRelBr;
    }
    if (moveType == 4) state->relBrLengths = oldRelBr;
    if (moveType == 7) state->kPrime = oldKPrime;
    return false;
  }

  // Evaluate likelihood (edge matrix + absolute lengths)
  int nEdge = state->relBrLengths.size();
  IntegerMatrix propEdge(nEdge, 2);
  NumericVector propEdgeLen(nEdge);
  for (int i = 0; i < nEdge; ++i) {
    propEdge(i, 0) = state->parent[i];
    propEdge(i, 1) = state->child[i];
    propEdgeLen[i] = state->treeLength * state->relBrLengths[i];
  }
  double newLogLik = cpp_log_likelihood(
    *data, propEdge, propEdgeLen, state->kPrime,
    state->rateLoss, state->rateLogSd, state->rateNeo);

  // Accept/reject
  double logAlpha = beta * (newLogLik - state->logLik) +
                    (newLogPrior - state->logPrior) + logHastings;

  if (R_FINITE(logAlpha) && std::log(R::unif_rand()) < logAlpha) {
    state->logLik   = newLogLik;
    state->logPrior = newLogPrior;
    return true;
  }

  // Reject: rollback
  state->treeLength = oldTL;
  state->rateLoss   = oldRL;
  state->rateLogSd  = oldRLSD;
  state->rateNeo    = oldRN;
  state->p          = oldP;
  if (topologyChanged) {
    state->parent       = oldParent;
    state->child        = oldChild;
    state->relBrLengths = oldRelBr;
  }
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
//
// stateXPtrs          List of nChains XPtr<McmcState>; index 0 = cold chain
// betas               Per-chain temperature (cold = 1.0)
// moveTypeCodes       Move type codes 0–8, length nMoves
// transIdxCpp         0-based global kPrime indices for transformational chars
// moveWeights         Unnormalized proposal weights, length nMoves
// chainScaleTunings   NumericMatrix nChains × nMoves — scale tuning per move
// chainBsmpTunings    NumericVector nChains — beta-simplex tuning
// chainIntWalkWins    IntegerVector nChains — int-walk window
// nBatch              Number of iterations to run
// startIter           Global iteration counter at batch start (1-based)
// warmup              Total warmup iterations; samples saved for iter > warmup
// thin                Save cold-chain state every thin-th post-warmup iter
// hasNeo              Whether rate_neo is a model parameter
// nEdge               Number of edges (length of relBrLengths)
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

  // Accept/propose counters (nChains × nMoves)
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
