#ifndef MKPRIME_CORRECTIONS_H
#define MKPRIME_CORRECTIONS_H

// Corrections added to a Felsenstein log-likelihood.

#include <cmath>
#include <stdexcept>

namespace MkPrime {

// Mk' relabelling correction for a character with kObs observed states and
// k' >= kObs true states. Under JC(k') every injective map from observed
// labels to model states gives the same likelihood, and there are
// k'! / (k' - kObs)! such maps; the log of that count is added to the
// log-likelihood. At k' = kObs it is the constant log(kObs!).
inline double relabel_log(int kPrime, int kObs) {
  if (kObs < 1) throw std::invalid_argument("relabel_log: kObs must be >= 1");
  if (kPrime < kObs) throw std::invalid_argument("relabel_log: kPrime < kObs");
  return std::lgamma(kPrime + 1.0) - std::lgamma(kPrime - kObs + 1.0);
}

// Lewis (2001) ascertainment correction for nChar characters, each excluded
// from the data with probability pExcluded (e.g. P(constant) under
// variable coding): add -nChar * log(1 - pExcluded) to the log-likelihood.
inline double ascertainment_log_correction(int nChar, double pExcluded) {
  return -nChar * std::log(1.0 - pExcluded);
}

}  // namespace MkPrime

#endif  // MKPRIME_CORRECTIONS_H
