// Engine extensions (#441): validation and attachment of the
// MkPrimeExtension structs described in inst/include/MkPrime/extension.h.
#include <Rcpp.h>
#include "mcmc_state.h"
#include "../inst/include/MkPrime/MkPrime.h"
#include <vector>
#include <cmath>
#include <string>

static const MkPrimeExtension* extension_from_sexp(SEXP x) {
  if (TYPEOF(x) != EXTPTRSXP ||
      R_ExternalPtrTag(x) != Rf_install("MkPrimeExtension")) {
    Rcpp::stop("Not an MkPrimeExtension pointer.");
  }
  const MkPrimeExtension* ext =
    static_cast<const MkPrimeExtension*>(R_ExternalPtrAddr(x));
  if (!ext) Rcpp::stop("MkPrimeExtension pointer is NULL.");
  if (ext->abi_version < MKPRIME_EXTENSION_ABI_MIN ||
      ext->abi_version > MKPRIME_EXTENSION_ABI_VERSION) {
    Rcpp::stop("Extension '%s' was built for extension ABI %d; MkPrime "
               "accepts %d to %d. Rebuild it against this MkPrime.",
               ext->name ? ext->name : "?", ext->abi_version,
               MKPRIME_EXTENSION_ABI_MIN, MKPRIME_EXTENSION_ABI_VERSION);
  }
  const bool v3 = ext->abi_version >= 3;
  const bool hasPrior = ext->log_prior || (v3 && ext->log_prior_ints);
  if (!ext->name || ext->n_params < 0 || !hasPrior ||
      (ext->n_params > 0 && (!ext->param_names || !ext->param_init)) ||
      ext->n_params + ext_n_ints(ext) < 1) {
    Rcpp::stop("Extension '%s' is incomplete: it needs a name, at least one "
               "parameter with names and initial values, and log_prior.",
               ext->name ? ext->name : "?");
  }
  const int nMoves = ext_n_moves(ext);
  const bool hasPropose = ext->propose || (v3 && ext->propose_ints);
  if (nMoves < 0 ||
      (nMoves > 0 && (!ext->move_names || !ext->move_weight || !hasPropose))) {
    Rcpp::stop("Extension '%s' declares %d moves; each needs a name, a "
               "weight and the propose callback.", ext->name, nMoves);
  }
  for (int m = 0; m < nMoves; ++m) {
    if (!ext->move_names[m] || !(ext->move_weight[m] > 0.0)) {
      Rcpp::stop("Extension '%s' move %d needs a name and a positive weight.",
                 ext->name, m + 1);
    }
  }
  const int nInts = ext_n_ints(ext);
  if (nInts < 0 ||
      (nInts > 0 && (!ext->int_names || !ext->int_init || !ext->int_levels))) {
    Rcpp::stop("Extension '%s' declares %d integer parameters; each needs a "
               "name, an initial value and a number of levels.", ext->name,
               nInts);
  }
  for (int i = 0; i < nInts; ++i) {
    if (!ext->int_names[i] || ext->int_levels[i] < 2 || ext->int_init[i] < 0 ||
        ext->int_init[i] >= ext->int_levels[i]) {
      Rcpp::stop("Extension '%s' integer parameter %d needs a name, at least "
                 "two levels and an initial value among them.", ext->name,
                 i + 1);
    }
  }
  if (ext_claims(ext) & ~MKPRIME_CLAIM_NEOMORPHIC) {
    Rcpp::stop("Extension '%s' makes an unknown claim (%d).", ext->name,
               ext_claims(ext));
  }
  if (ext_claims(ext) && !ext_has_log_lik(ext)) {
    Rcpp::stop("Extension '%s' claims characters but cannot score them: it "
               "needs a log-likelihood callback.", ext->name);
  }
  return ext;
}

// [[Rcpp::export]]
Rcpp::List extension_info_cpp(SEXP ext) {
  const MkPrimeExtension* e = extension_from_sexp(ext);
  Rcpp::CharacterVector names(e->n_params);
  Rcpp::NumericVector init(e->n_params);
  for (int i = 0; i < e->n_params; ++i) {
    if (!e->param_names[i]) Rcpp::stop("Extension parameter %d has no name.", i + 1);
    names[i] = e->param_names[i];
    init[i] = e->param_init[i];
  }
  const int nMoves = ext_n_moves(e);
  Rcpp::CharacterVector moveNames(nMoves);
  Rcpp::NumericVector moveWeights(nMoves);
  for (int m = 0; m < nMoves; ++m) {
    moveNames[m] = e->move_names[m];
    moveWeights[m] = e->move_weight[m];
  }
  const int nInts = ext_n_ints(e);
  Rcpp::CharacterVector intNames(nInts);
  Rcpp::IntegerVector intInit(nInts), intLevels(nInts);
  for (int i = 0; i < nInts; ++i) {
    intNames[i] = e->int_names[i];
    intInit[i] = e->int_init[i];
    intLevels[i] = e->int_levels[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("name") = std::string(e->name),
    Rcpp::Named("abi_version") = e->abi_version,
    Rcpp::Named("params") = names,
    Rcpp::Named("init") = init,
    Rcpp::Named("log_lik") = ext_has_log_lik(e),
    Rcpp::Named("moves") = moveNames,
    Rcpp::Named("move_weights") = moveWeights,
    Rcpp::Named("ints") = intNames,
    Rcpp::Named("int_init") = intInit,
    Rcpp::Named("int_levels") = intLevels,
    Rcpp::Named("claims_neomorphic") =
      (ext_claims(e) & MKPRIME_CLAIM_NEOMORPHIC) != 0
  );
}

// [[Rcpp::export]]
double extension_log_prior_cpp(SEXP ext, Rcpp::NumericVector params,
                               Rcpp::IntegerVector ints = Rcpp::IntegerVector::create()) {
  const MkPrimeExtension* e = extension_from_sexp(ext);
  if (params.size() != e->n_params) {
    Rcpp::stop("Extension '%s' has %d parameters; got %d values.", e->name,
               e->n_params, static_cast<int>(params.size()));
  }
  if (ints.size() != ext_n_ints(e)) {
    Rcpp::stop("Extension '%s' has %d integer parameters; got %d values.",
               e->name, ext_n_ints(e), static_cast<int>(ints.size()));
  }
  return ext_log_prior(e, params.begin(), ints.begin());
}

// Move the neomorphic partitions out of `d->parts` into the claiming
// extension's data. Their characters stay in nNeo, so rate_neo is
// normalised over them as before, and hasNeo is unchanged.
static void claim_neomorphic(McmcData* d, int owner) {
  if (d->nClasses > 1) {
    Rcpp::stop("Extensions cannot claim characters under a user partition.");
  }
  d->extClaimNeo = owner;
  std::vector<PartInfo> kept;
  std::vector<MaskTally> keptTally;
  for (size_t pi = 0; pi < d->parts.size(); ++pi) {
    const PartInfo& part = d->parts[pi];
    if (part.type != 0) {
      kept.push_back(part);
      if (pi < d->partMaskTally.size()) keptTally.push_back(d->partMaskTally[pi]);
      continue;
    }
    const int nCharPart = part.tipStates.ncol();
    for (int c = 0; c < nCharPart; ++c) {
      for (int t = 0; t < d->nTip; ++t) {
        d->extClaimStates.push_back(part.tipStates(t, c));
      }
      d->extClaimChar.push_back(part.globalCharIdx[c] + 1);
    }
  }
  d->parts.swap(kept);
  d->partMaskTally.swap(keptTally);
  d->neoPartIndices.clear();
  d->charToPartition.assign(d->nChar, -1);
  for (size_t pi = 0; pi < d->parts.size(); ++pi) {
    const PartInfo& part = d->parts[pi];
    for (R_xlen_t c = 0; c < part.globalCharIdx.size(); ++c) {
      d->charToPartition[part.globalCharIdx[c]] = static_cast<int>(pi);
    }
  }
}

// Attach extensions, in order, to an McmcData. Their parameters follow one
// another in McmcState::extParams (and extInts). Call once per McmcData:
// a claim removes partitions that a second call could not restore.
// [[Rcpp::export]]
void attach_extensions_cpp(SEXP dataPtr, Rcpp::List extensions) {
  McmcData* d = Rcpp::XPtr<McmcData>(dataPtr);
  if (d->extClaimNeo >= 0 || !d->extensions.empty()) {
    Rcpp::stop("Extensions are already attached to this data.");
  }
  std::vector<const MkPrimeExtension*> exts;
  for (R_xlen_t i = 0; i < extensions.size(); ++i) {
    exts.push_back(extension_from_sexp(extensions[i]));
  }
  int claimant = -1;
  for (size_t i = 0; i < exts.size(); ++i) {
    if (ext_claims(exts[i]) & MKPRIME_CLAIM_NEOMORPHIC) {
      if (claimant >= 0) {
        Rcpp::stop("Extensions '%s' and '%s' both claim the neomorphic "
                   "characters.", exts[claimant]->name, exts[i]->name);
      }
      claimant = static_cast<int>(i);
    }
  }
  d->nExtParams = 0;
  d->nExtInts = 0;
  d->extLik = false;
  for (size_t i = 0; i < exts.size(); ++i) {
    const MkPrimeExtension* e = exts[i];
    d->extensions.push_back(e);
    d->extKeep.push_back(Rcpp::RObject(extensions[i]));
    d->extOffset.push_back(d->nExtParams);
    d->nExtParams += e->n_params;
    d->extIntOffset.push_back(d->nExtInts);
    for (int j = 0; j < ext_n_ints(e); ++j) {
      d->extIntLevels.push_back(e->int_levels[j]);
    }
    d->nExtInts += ext_n_ints(e);
    d->extLik = d->extLik || ext_has_log_lik(e);
    for (int m = 0; m < ext_n_moves(e); ++m) {
      d->extMoves.push_back({static_cast<int>(i), m});
    }
  }
  if (claimant >= 0) claim_neomorphic(d, claimant);
}

// [[Rcpp::export]]
int extension_abi_version() {
  return MKPRIME_EXTENSION_ABI_VERSION;
}


// Test fixture: two parameters with a known prior and no likelihood role.
//   toy_sigma ~ LogNormal(0, 0.5),  toy_rate ~ Gamma(shape 2, rate 4).
static double toy_log_prior(const double* x, void*) {
  if (!(x[0] > 0.0) || !(x[1] > 0.0)) return R_NegInf;
  const double sdlog = 0.5, lx = std::log(x[0]);
  const double lnorm = -lx - std::log(sdlog) - 0.5 * std::log(2.0 * M_PI) -
    0.5 * lx * lx / (sdlog * sdlog);
  const double shape = 2.0, rate = 4.0;
  const double lgam = shape * std::log(rate) - std::lgamma(shape) +
    (shape - 1.0) * std::log(x[1]) - rate * x[1];
  return lnorm + lgam;
}

static const char* const kToyNames[] = {"toy_sigma", "toy_rate"};
static const double kToyInit[] = {1.0, 0.5};
static const MkPrimeExtension kToyExtension = {
  MKPRIME_EXTENSION_ABI_VERSION, "toy", 2, kToyNames, kToyInit,
  toy_log_prior, nullptr
};

// [[Rcpp::export]]
SEXP toy_extension_cpp() {
  return MkPrimeExtensionPtr(&kToyExtension);
}

// ABI v1 build of the same toy, to check that version 1 structs still load.
static const MkPrimeExtension kToyV1Extension = {
  1, "toy_v1", 2, kToyNames, kToyInit, toy_log_prior, nullptr
};

// [[Rcpp::export]]
SEXP toy_v1_extension_cpp() {
  return MkPrimeExtensionPtr(&kToyV1Extension);
}

// Test fixture with a likelihood (ABI v2): kToyEvents events of a Poisson
// process at rate toy_lambda along the whole tree, so
//   log L = N log(lambda T) - lambda T - log N!,  T = tree length,
// with toy_lambda ~ Gamma(2, 4); toy_lambda | T ~ Gamma(2 + N, 4 + T).
// Moves: "gibbs" draws from that conditional; "jump" is a log-scale
// random walk that declines a step beyond exp(+-0.45).
static const int kToyEvents = 3;

static double toy_tree_length(const MkPrimeTreeView* tree) {
  double t = 0.0;
  for (int i = 0; i < tree->n_edge; ++i) t += tree->edge_length[i];
  return t;
}

static double toy_lik_log_prior(const double* x, void*) {
  if (!(x[0] > 0.0)) return R_NegInf;
  return 2.0 * std::log(4.0) - std::lgamma(2.0) + std::log(x[0]) - 4.0 * x[0];
}

static double toy_lik_log_lik(const MkPrimeTreeView* tree, const double* x,
                              void*) {
  const double mu = x[0] * toy_tree_length(tree);
  if (!(mu > 0.0)) return R_NegInf;
  return kToyEvents * std::log(mu) - mu - std::lgamma(kToyEvents + 1.0);
}

// Gamma(shape, rate) log density, for an integer shape.
static double toy_log_dgamma(double x, int shape, double rate) {
  return shape * std::log(rate) - std::lgamma(static_cast<double>(shape)) +
    (shape - 1) * std::log(x) - rate * x;
}

static double toy_lik_propose(int move, const MkPrimeTreeView* tree,
                              double* x, void*) {
  const double old = x[0];
  if (move == 0) {
    const int shape = 2 + kToyEvents;
    const double rate = 4.0 + toy_tree_length(tree);
    double logProd = 0.0;
    for (int i = 0; i < shape; ++i) logProd += std::log(tree->unif_rand());
    x[0] = -logProd / rate;
    return toy_log_dgamma(old, shape, rate) - toy_log_dgamma(x[0], shape, rate);
  }
  const double step = tree->unif_rand() - 0.5;
  if (std::fabs(step) > 0.45) return R_NegInf;
  x[0] = old * std::exp(step);
  return step;
}

static const char* const kToyLikNames[] = {"toy_lambda"};
static const double kToyLikInit[] = {0.5};
static const char* const kToyLikMoves[] = {"gibbs", "jump"};
static const double kToyLikMoveWeight[] = {1.0, 0.5};
static const MkPrimeExtension kToyLikExtension = {
  MKPRIME_EXTENSION_ABI_VERSION, "toylik", 1, kToyLikNames, kToyLikInit,
  toy_lik_log_prior, nullptr,
  toy_lik_log_lik, 2, kToyLikMoves, kToyLikMoveWeight, toy_lik_propose
};

// [[Rcpp::export]]
SEXP toy_lik_extension_cpp() {
  return MkPrimeExtensionPtr(&kToyLikExtension);
}

// Test fixture (ABI v3): claims the neomorphic characters and scores them
// as MkPrime would without it (MkN at rate_loss, edge lengths scaled by
// neo_scale, the ACRV rates; no ascertainment correction, so use it with
// coding = "none"). The total likelihood is then unchanged. It adds an
// integer toy_switch (levels 0, 1; prior 0.3, 0.7), moved by "flip", and
// toy_mu ~ Gamma(2, 4), which enters nothing but its prior.
static double toyclaim_log_prior(const double* x, const int* z, void*) {
  if (!(x[0] > 0.0)) return R_NegInf;
  return 2.0 * std::log(4.0) + std::log(x[0]) - 4.0 * x[0] +
    std::log(z[0] == 1 ? 0.7 : 0.3);
}

static double toyclaim_log_lik(const MkPrimeTreeView* v, const double*,
                               const int*, void*) {
  if (v->n_claimed == 0) return 0.0;
  std::vector<double> len(v->edge_length, v->edge_length + v->n_edge);
  for (double& l : len) l *= v->neo_scale;
  const std::vector<double> rates(v->rate_cat, v->rate_cat + v->n_rate_cat);
  MkPrime::Pruner<MkPrime::MkNModel> pruner(
    MkPrime::MkNModel(v->core[MKPRIME_CORE_RATE_LOSS]), v->n_tip, v->n_claimed,
    v->claimed_states, rates);
  pruner.set_tree(v->parent, v->child, len.data(), v->n_edge);
  return pruner.log_likelihood();
}

static double toyclaim_propose(int, const MkPrimeTreeView*, double*, int* z,
                               void*) {
  z[0] = 1 - z[0];
  return 0.0;
}

static const char* const kToyClaimNames[] = {"toy_mu"};
static const double kToyClaimInit[] = {0.5};
static const char* const kToyClaimMoves[] = {"flip"};
static const double kToyClaimMoveWeight[] = {1.0};
static const char* const kToyClaimInts[] = {"toy_switch"};
static const int kToyClaimIntInit[] = {0};
static const int kToyClaimIntLevels[] = {2};
static const MkPrimeExtension kToyClaimExtension = {
  MKPRIME_EXTENSION_ABI_VERSION, "toyclaim", 1, kToyClaimNames, kToyClaimInit,
  nullptr, nullptr,
  nullptr, 1, kToyClaimMoves, kToyClaimMoveWeight, nullptr,
  MKPRIME_CLAIM_NEOMORPHIC, 1, kToyClaimInts, kToyClaimIntInit,
  kToyClaimIntLevels, toyclaim_log_prior, toyclaim_log_lik, toyclaim_propose
};

// [[Rcpp::export]]
SEXP toy_claim_extension_cpp() {
  return MkPrimeExtensionPtr(&kToyClaimExtension);
}
