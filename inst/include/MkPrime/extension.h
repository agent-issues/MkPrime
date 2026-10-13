#ifndef MKPRIME_EXTENSION_H
#define MKPRIME_EXTENSION_H

/* Engine extension interface (#441), ABI version 1.
 *
 * An extension package describes its parameters in a static
 * MkPrimeExtension and returns MkPrimeExtensionPtr(&ext) from a registered
 * .Call entry point taking no arguments. R code then passes
 * MkPrimeExtension(entry, package) to MkPrimeModel(extensions = ...).
 * MkPrime resolves the entry point afresh in every R process that builds a
 * chain, so the pointer itself is never serialised.
 *
 * Version 1 contributes parameters and their prior: each parameter is
 * logged, checkpointed, carried through tempering swaps and updated by a
 * Bactrian scale move. Parameters do not yet enter the likelihood.
 *
 * Plain C, so that C and C++ packages can use it. The struct and everything
 * it points to must outlive the R session (static storage).
 */

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>

#define MKPRIME_EXTENSION_ABI_VERSION 1

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MkPrimeExtension {
  /* MKPRIME_EXTENSION_ABI_VERSION as compiled by the extension. */
  int abi_version;
  /* Short identifier, e.g. the package name. */
  const char* name;
  /* Parameters, each a positive real. Names become log columns and must not
   * clash with MkPrime's own or another extension's. */
  int n_params;
  const char* const* param_names;
  const double* param_init;
  /* Log prior density of this extension's parameters (length n_params);
   * -INFINITY outside their support. `user` is passed through unchanged. */
  double (*log_prior)(const double* params, void* user);
  void* user;
} MkPrimeExtension;

/* External pointer handed back to R by the extension's entry point. */
static inline SEXP MkPrimeExtensionPtr(const MkPrimeExtension* ext) {
  return R_MakeExternalPtr((void*) ext, Rf_install("MkPrimeExtension"),
                           R_NilValue);
}

#ifdef __cplusplus
}
#endif

#endif /* MKPRIME_EXTENSION_H */
