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
  if (ext->abi_version != MKPRIME_EXTENSION_ABI_VERSION) {
    Rcpp::stop("Extension '%s' was built for extension ABI %d; MkPrime "
               "provides %d. Rebuild it against this MkPrime.",
               ext->name ? ext->name : "?", ext->abi_version,
               MKPRIME_EXTENSION_ABI_VERSION);
  }
  if (!ext->name || ext->n_params < 1 || !ext->param_names ||
      !ext->param_init || !ext->log_prior) {
    Rcpp::stop("Extension '%s' is incomplete: it needs a name, at least one "
               "parameter with names and initial values, and log_prior.",
               ext->name ? ext->name : "?");
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
  return Rcpp::List::create(
    Rcpp::Named("name") = std::string(e->name),
    Rcpp::Named("abi_version") = e->abi_version,
    Rcpp::Named("params") = names,
    Rcpp::Named("init") = init
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
  for (R_xlen_t i = 0; i < extensions.size(); ++i) {
    const MkPrimeExtension* e = extension_from_sexp(extensions[i]);
    d->extensions.push_back(e);
    d->extOffset.push_back(d->nExtParams);
    d->nExtParams += e->n_params;
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
