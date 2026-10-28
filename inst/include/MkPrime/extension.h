#ifndef MKPRIME_EXTENSION_H
#define MKPRIME_EXTENSION_H

/* Engine extension interface (#441), ABI version 3.
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
 * Version 3 adds:
 *  - integer parameters (e.g. latent categories), each in 0..levels - 1,
 *    logged, checkpointed and carried through tempering like the real ones;
 *    the *_ints callbacks see them and, when given, replace the v1/v2 ones;
 *  - MkPrime's own scalars, the neomorphic edge-length scale and the ACRV
 *    rate multipliers in the tree view, so an extension can share them;
 *  - claims: an extension may score MkPrime's neomorphic characters itself
 *    (MKPRIME_CLAIM_NEOMORPHIC). MkPrime then leaves them out of its own
 *    likelihood, but still counts them when it normalises rate_neo, samples
 *    rate_loss and rate_neo, and logs them. Their tip states arrive in the
 *    tree view. The extension owns any ascertainment correction for them.
 *
 * A struct built against an earlier version remains valid: MkPrime reads
 * the fields of version v only when abi_version >= v.
 *
 * Plain C, so that C and C++ packages can use it. The struct and everything
 * it points to must outlive every chain that uses it: static storage, or
 * memory owned by the external pointer, which MkPrime keeps alive for as
 * long as it holds the pointer.
 */

#ifndef R_NO_REMAP
#define R_NO_REMAP
#endif
#include <Rinternals.h>

#define MKPRIME_EXTENSION_ABI_VERSION 3
/* Oldest ABI version this MkPrime still accepts. */
#define MKPRIME_EXTENSION_ABI_MIN 1

/* Indices into MkPrimeTreeView::core (ABI 3). */
#define MKPRIME_CORE_TREE_LENGTH 0
#define MKPRIME_CORE_RATE_LOSS 1
#define MKPRIME_CORE_RATE_LOG_SD 2
#define MKPRIME_CORE_RATE_NEO 3
#define MKPRIME_CORE_BETA_SCALE 4

/* MkPrimeExtension::claims (ABI 3). */
#define MKPRIME_CLAIM_NEOMORPHIC 1

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

  /* ---- ABI version 3; read only by extensions built for 3 or later ---- */

  /* MkPrime's scalars, indexed by MKPRIME_CORE_*; n_core entries. */
  int n_core;
  const double* core;
  /* Multiplier of edge_length for neomorphic characters: rate_neo's share
   * of the rate, normalised over all characters, claimed ones included. */
  double neo_scale;
  /* Among-character rate variation: n_rate_cat equally weighted rate
   * multipliers from rate_log_sd ({1} when rate_log_sd is 0). */
  int n_rate_cat;
  const double* rate_cat;
  /* Ascertainment coding: 0 none, 1 variable, 2 informative. */
  int coding;
  /* Characters claimed by this extension: tip states (n_tip x n_claimed,
   * column-major, -1 = missing) and their 1-based column in the data. */
  int n_claimed;
  const int* claimed_states;
  const int* claimed_char;
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

  /* ---- ABI version 3 ---- */

  /* MKPRIME_CLAIM_* flags; 0 for none. At most one extension may claim
   * a given kind of character. */
  int claims;
  /* Integer parameters, each in 0 .. int_levels[i] - 1. Names become log
   * columns. MkPrime gives them no moves of its own. */
  int n_ints;
  const char* const* int_names;
  const int* int_init;
  const int* int_levels;
  /* As log_prior, log_lik and propose, with the integer parameters. When
   * given they are used instead of those; propose_ints may change both
   * blocks, and MkPrime restores both on rejection. */
  double (*log_prior_ints)(const double* params, const int* ints, void* user);
  double (*log_lik_ints)(const MkPrimeTreeView* tree, const double* params,
                         const int* ints, void* user);
  double (*propose_ints)(int move, const MkPrimeTreeView* tree,
                         double* params, int* ints, void* user);
} MkPrimeExtension;

/* External pointer handed back to R by the extension's entry point. */
static inline SEXP MkPrimeExtensionPtr(const MkPrimeExtension* ext) {
  return R_MakeExternalPtr((void*) ext, Rf_install("MkPrimeExtension"),
                           R_NilValue);
}

/* As MkPrimeExtensionPtr(), keeping `prot` (e.g. the R object that owns a
 * heap-allocated struct and its data) alive as long as the pointer. */
static inline SEXP MkPrimeExtensionPtrProtect(const MkPrimeExtension* ext,
                                              SEXP prot) {
  return R_MakeExternalPtr((void*) ext, Rf_install("MkPrimeExtension"), prot);
}

#ifdef __cplusplus
}
#endif

#endif /* MKPRIME_EXTENSION_H */
