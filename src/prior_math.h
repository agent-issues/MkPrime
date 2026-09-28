#ifndef MKP_PRIOR_MATH_H
#define MKP_PRIOR_MATH_H

#include <cmath>

namespace mkp {

// log(1 - exp(x)) for x <= 0 (Maechler 2012, "Accurately computing
// log(1 - exp(-|a|))"). Neither one-branch form is safe across the range the
// truncation normalisers reach: log1p(-exp(x)) underflows to -Inf as x -> 0-
// (tiny p), and log(-expm1(x)) rounds to exactly 0 once exp(x) < eps (x far
// below 0), losing the whole value. The -log(2) switch keeps both accurate.
// Mirrored by .Log1mExp() in R/MkPrimeModel.R.
inline double log1m_exp(double x) {
  constexpr double kLn2 = 0.6931471805599453;
  return (x > -kLn2) ? std::log(-std::expm1(x)) : std::log1p(-std::exp(x));
}

// log of the logseries normaliser over k' >= 2:
//   sum_{k>=2} c^k / k = -log(1 - c) - c.
// Below c = 0.25 the difference cancels, so sum the series directly; its terms
// fall by at least 4x each, so 60 terms reach well below double precision.
// Mirrored by .LogseriesLogNorm() in R/MkPrimeModel.R.
inline double logseries_log_norm(double c) {
  double s;
  if (c < 0.25) {
    s = 0.0;
    double ck = c;
    for (int k = 2; k <= 60; ++k) {
      ck *= c;
      s += ck / k;
    }
  } else {
    s = -std::log1p(-c) - c;
  }
  return std::log(s);
}

}  // namespace mkp

#endif  // MKP_PRIOR_MATH_H
