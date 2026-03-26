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
// do_move_cpp: full propose/evaluate/accept — updates state in-place
//
// moveType: 0=scale_tl, 1=scale_rl, 2=scale_rls, 3=scale_rn,
//           4=beta_simplex, 5=nni, 6=spr, 7=int_walk, 8=scale_p
// ---------------------------------------------------------------------------

// [[Rcpp::export]]
bool do_move_cpp(SEXP dataPtr, SEXP statePtr,
                 int moveType, int charIdx,
                 double scaleTuning, double betaSimplexTuning,
                 int intWalkWindow, double beta) {
  McmcData*  data  = Rcpp::XPtr<McmcData>(dataPtr).get();
  McmcState* state = Rcpp::XPtr<McmcState>(statePtr).get();

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
