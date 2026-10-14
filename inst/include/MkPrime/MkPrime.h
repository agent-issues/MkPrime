#ifndef MKPRIME_MKPRIME_H
#define MKPRIME_MKPRIME_H

// MkPrime likelihood kernels for downstream packages.
//
// In DESCRIPTION, `LinkingTo: MkPrime`; then `#include <MkPrime/MkPrime.h>`.
// Everything is header-only C++11 in namespace MkPrime and needs neither R,
// Rcpp nor MkPrime's compiled code. MkPrime's own src/ uses these headers,
// so results agree with MkPrime bit for bit on a given build.
//
//   transition.h     JC(k) / Mk(k) / Mk', MkN and F81 transition kernels
//   models.h         JCModel, JCLumpedModel, MkNModel, F81Model
//   pruning.h        Pruner<Model>: cached pruning, checkpoint / reject
//   acrv.h           lognormal_rates()
//   corrections.h    relabel_log(), ascertainment_log_correction()
//   ascertainment.h  constant_patterns()
//
// MKPRIME_HEADER_API_VERSION increases on any change to these signatures or
// to their numerical results.

#define MKPRIME_HEADER_API_VERSION 1

#include "transition.h"
#include "models.h"
#include "pruning.h"
#include "acrv.h"
#include "corrections.h"
#include "ascertainment.h"

#endif  // MKPRIME_MKPRIME_H
