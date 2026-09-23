#ifndef MKPRIME_F81_H
#define MKPRIME_F81_H

#include <cmath>

// F81 transition: (P × cl)_i = (1 - e^{-μt}) × dot(π, cl) + e^{-μt} × cl_i
// M-114: used by Q-heterogeneity partial CL
inline void f81_transition(const double* cl, double* result,
                            int nChar, int kStates,
                            const double* pi, double mu, double t) {
  // FAST-EXP-001: expm1 form is exact even at tiny mu*t.
  double arg = -mu * t;
  double one_minus_exp = -std::expm1(arg);
  double exp_t         = 1.0 - one_minus_exp;
  for (int c = 0; c < nChar; ++c) {
    int off = c * kStates;
    double piDotCl = 0.0;
    for (int s = 0; s < kStates; ++s)
      piDotCl += pi[s] * cl[off + s];
    double baseTerm = one_minus_exp * piDotCl;
    for (int s = 0; s < kStates; ++s)
      result[off + s] = baseTerm + exp_t * cl[off + s];
  }
}

#endif  // MKPRIME_F81_H
