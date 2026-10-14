// Engine extensions (#441): validation and attachment of the
// MkPrimeExtension structs described in inst/include/MkPrime/extension.h.
#include <Rcpp.h>
#include "mcmc_state.h"
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
  if (!ext->name || ext->n_params < 1 || !ext->param_names ||
      !ext->param_init || !ext->log_prior) {
    Rcpp::stop("Extension '%s' is incomplete: it needs a name, at least one "
               "parameter with names and initial values, and log_prior.",
               ext->name ? ext->name : "?");
  }
  const int nMoves = ext_n_moves(ext);
  if (nMoves < 0 ||
      (nMoves > 0 && (!ext->move_names || !ext->move_weight || !ext->propose))) {
    Rcpp::stop("Extension '%s' declares %d moves; each needs a name, a "
               "weight and the propose callback.", ext->name, nMoves);
  }
  for (int m = 0; m < nMoves; ++m) {
    if (!ext->move_names[m] || !(ext->move_weight[m] > 0.0)) {
      Rcpp::stop("Extension '%s' move %d needs a name and a positive weight.",
                 ext->name, m + 1);
    }
  }
  return ext;
}

// [[Rcpp::export]]
Rcpp::List extension_info_cpp(SEXP ext) {
  const MkPrimeExtension* e = extension_from_sexp(ext);
  Rcpp::CharacterVector names(e->n_params);
  Rcpp::NumericVector init(e->n_params);
  const int nMoves = ext_n_moves(e);
  Rcpp::CharacterVector moveNames(nMoves);
  Rcpp::NumericVector moveWeights(nMoves);
  for (int m = 0; m < nMoves; ++m) {
    moveNames[m] = e->move_names[m];
    moveWeights[m] = e->move_weight[m];
  }
  for (int i = 0; i < e->n_params; ++i) {
    if (!e->param_names[i]) Rcpp::stop("Extension parameter %d has no name.", i + 1);
    names[i] = e->param_names[i];
    init[i] = e->param_init[i];
  }
  return Rcpp::List::create(
    Rcpp::Named("name") = std::string(e->name),
    Rcpp::Named("abi_version") = e->abi_version,
    Rcpp::Named("params") = names,
    Rcpp::Named("init") = init,
    Rcpp::Named("log_lik") = ext_has_log_lik(e),
    Rcpp::Named("moves") = moveNames,
    Rcpp::Named("move_weights") = moveWeights
  );
}

// [[Rcpp::export]]
double extension_log_prior_cpp(SEXP ext, Rcpp::NumericVector params) {
  const MkPrimeExtension* e = extension_from_sexp(ext);
  if (params.size() != e->n_params) {
    Rcpp::stop("Extension '%s' has %d parameters; got %d values.", e->name,
               e->n_params, static_cast<int>(params.size()));
  }
  return e->log_prior(params.begin(), e->user);
}

// Attach extensions, in order, to an McmcData. Their parameters follow one
// another in McmcState::extParams.
// [[Rcpp::export]]
void attach_extensions_cpp(SEXP dataPtr, Rcpp::List extensions) {
  McmcData* d = Rcpp::XPtr<McmcData>(dataPtr);
  d->extensions.clear();
  d->extOffset.clear();
  d->nExtParams = 0;
  d->extLik = false;
  d->extMoves.clear();
  for (R_xlen_t i = 0; i < extensions.size(); ++i) {
    const MkPrimeExtension* e = extension_from_sexp(extensions[i]);
    d->extensions.push_back(e);
    d->extOffset.push_back(d->nExtParams);
    d->nExtParams += e->n_params;
    d->extLik = d->extLik || ext_has_log_lik(e);
    for (int m = 0; m < ext_n_moves(e); ++m) {
      d->extMoves.push_back({static_cast<int>(i), m});
    }
  }
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
