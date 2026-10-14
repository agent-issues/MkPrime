#ifndef MKPRIME_EXTENSION_H
#define MKPRIME_EXTENSION_H

/* Engine extension interface (#441), ABI version 2.
 *
 * An extension package describes itself in a static MkPrimeExtension and
 * returns MkPrimeExtensionPtr(&ext) from a registered .Call entry point.
 * R code wraps a function returning that pointer,
 *   MkPrimeExtension(function() .Call(entryPoint)),
 * and passes it to MkPrimeModel(extensions = ...). MkPrime calls the
 * function afresh in every R process that builds a chain, so the pointer
 * itself is never serialised.
 *
 * Version 1: parameters and their prior. Each parameter is a positive real
 * that is logged, checkpointed, carried through tempering swaps and updated
 * by a Bactrian scale move.
 *
 * Version 2 adds, all optional (NULL / 0 for none):
 *  - log_lik: a log-likelihood term over the extension's own data, given
 *    the tree and the extension's parameters. It is added to MkPrime's
 *    likelihood and tempered with it. Characters an extension models are
 *    given to the extension, not to MkPrimeData().
 *  - moves: Metropolis-Hastings proposals on the extension's parameters.
 *
 * Callbacks are called once per move or likelihood evaluation, never per
 * edge or character, and from a single thread per chain. They must not
 * keep pointers into the MkPrimeTreeView after returning.
 *
 * A struct built against version 1 remains valid: MkPrime reads the
 * version 2 fields only when abi_version >= 2.
 *
 * Plain C, so that C and C++ packages can use it. The struct and everything
 * it points to must outlive the R session (static storage).
 */

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>

#define MKPRIME_EXTENSION_ABI_VERSION 2
/* Oldest ABI version this MkPrime still accepts. */
#define MKPRIME_EXTENSION_ABI_MIN 1

#ifdef __cplusplus
extern "C" {
#endif

/* The tree a callback sees. Edge lists use 1-based node numbers (tips
 * 1..n_tip, root n_tip + 1) in preorder; edge_length is absolute (expected
 * changes per character). */
typedef struct MkPrimeTreeView {
  int n_tip;
  int n_edge;
  const int* parent;
  const int* child;
  const double* edge_length;
  /* Uniform(0, 1) deviate from the chain's random number stream. Use it,
   * not another generator, so that runs are reproducible. */
  double (*unif_rand)(void);
} MkPrimeTreeView;

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

  /* ---- ABI version 2 ---- */

  /* Log-likelihood of the extension's data given the tree and its
   * parameters; -INFINITY if impossible. NULL: no likelihood term. */
  double (*log_lik)(const MkPrimeTreeView* tree, const double* params,
                    void* user);
  /* Custom moves, each drawn with relative weight move_weight[i] alongside
   * MkPrime's own (whose weights are of order 1 to nEdge). */
  int n_moves;
  const char* const* move_names;
  const double* move_weight;
  /* Propose move `move` (0-based) by changing `params` in place, and return
   * the log Hastings ratio, log q(old | new) - log q(new | old). Return
   * -INFINITY to decline. MkPrime scores the proposal (prior and likelihood,
   * tempered) and restores `params` if it is rejected. A draw from a full
   * conditional is a proposal whose Hastings ratio cancels the target
   * ratio; it is then always accepted in the cold chain. */
  double (*propose)(int move, const MkPrimeTreeView* tree, double* params,
                    void* user);
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
