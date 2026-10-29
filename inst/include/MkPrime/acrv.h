#ifndef MKPRIME_ACRV_H
#define MKPRIME_ACRV_H

// Among-character rate variation: discrete lognormal rate categories.

#include <cmath>

namespace MkPrime {

// Rates for nCat equiprobable categories of a lognormal with log-scale
// standard deviation sd and mean 1: out[i] is proportional to
// exp(-sd^2 / 2 + sd * z[i]), rescaled so the mean of out is exactly 1.
// z[i] are standard-normal quantiles, z[i] = qnorm((i + 0.5) / nCat) in
// MkPrime (precompute once; R's qnorm() or R::qnorm()). sd <= 0 gives 1s.
inline void lognormal_rates(double sd, int nCat, const double* z, double* out) {
  if (sd <= 0.0) {
    for (int i = 0; i < nCat; ++i) out[i] = 1.0;
    return;
  }
  double mu = -sd * sd / 2.0;
  double total = 0.0;
  for (int i = 0; i < nCat; ++i) {
    out[i] = std::exp(mu + sd * z[i]);
    total += out[i];
  }
  for (int i = 0; i < nCat; ++i) out[i] *= nCat / total;
}

}  // namespace MkPrime

#endif  // MKPRIME_ACRV_H
