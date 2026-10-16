#ifndef MKPRIME_SCALING_H
#define MKPRIME_SCALING_H

// Conditional-likelihood rescaling (header API v3).
//
// A character's CLs at a node are rescaled only once their sum falls below
// kRescaleBelow, so characters that never come near underflow keep their
// arithmetic unchanged. Rescaling multiplies by a power of two, which is
// exact: stored CLs are the true ones times 2^shift, and the integer shift
// is carried alongside them, per (rate category, node, character).
//
// kRescaleBelow is 2^-256. A binary node multiplies the CLs of two
// children, and the root of an unrooted tree three, none of which has been
// left below 2^-256, so a product stays above 2^-768 times the spread across
// states -- far from the subnormal range at 2^-1022. A polytomy of more than
// three children has less headroom. A character must reach a
// log-likelihood below about -177 before any of its CLs is rescaled.

#include <cmath>

namespace MkPrime {

const double kRescaleBelow = 8.6361685550944446e-78;  // 2^-256
const double kLn2 = 0.69314718055994530942;

// Rescales cl[0..k) so that its largest entry lies in [0.5, 1); returns the
// shift applied (stored = true * 2^shift), 0 if every entry is zero.
inline int rescale_cl(double* cl, int k) {
  double mx = 0.0;
  for (int i = 0; i < k; ++i) {
    if (cl[i] > mx) mx = cl[i];
  }
  if (!(mx > 0.0)) return 0;
  int e;
  std::frexp(mx, &e);
  for (int i = 0; i < k; ++i) cl[i] = std::ldexp(cl[i], -e);
  return -e;
}

// acc * 2^-accShift += x * 2^-xShift. The sum keeps the shift of its larger
// term, so neither term overflows and the smaller one can only underflow
// where it is negligible.
inline void add_scaled(double& acc, int& accShift, double x, int xShift) {
  if (xShift == accShift) {
    acc += x;
  } else if (x == 0.0) {
    return;
  } else if (acc == 0.0 || xShift < accShift) {
    acc = std::ldexp(acc, xShift - accShift) + x;
    accShift = xShift;
  } else {
    acc += std::ldexp(x, accShift - xShift);
  }
}

// log(x * 2^-shift); x <= 0 gives -Inf or NaN as std::log(x) would.
inline double log_scaled(double x, int shift) {
  return shift ? std::log(x) - shift * kLn2 : std::log(x);
}

}  // namespace MkPrime

#endif  // MKPRIME_SCALING_H
