#ifndef MKPRIME_TRANSITION_H
#define MKPRIME_TRANSITION_H

// Transition probabilities for the Mk family at branch length t, measured in
// expected changes per character. Header-only; needs neither R nor Rcpp.
//
// Mk(k) and Mk' (with k = k') are both JC(k): use JCRates. MkN is the
// asymmetric binary model: use MkNRates. F81 acts on CL vectors directly.

#include <cmath>

namespace MkPrime {

// Entries of the JC(k) transition matrix P(t).
// diffCoeff and expTerm are equal analytically (both pSame - pDiff), but are
// computed differently and can differ in the last bit; both are exposed so
// that callers keep the rounding they were validated with.
struct JCProbs {
  double pSame;      // P_ii(t)
  double pDiff;      // P_ij(t), i != j
  double diffCoeff;  // pSame - pDiff
  double expTerm;    // exp(-k t / (k - 1))
};

// JC(k) with Q normalised to one expected change per unit time:
// Q_ii = -1, Q_ij = 1 / (k - 1). Construct once per k; probs() per branch.
class JCRates {
 public:
  explicit JCRates(int k) : k_(k), invK_(1.0 / k), km1_(k - 1.0) {}

  int k() const { return k_; }

  JCProbs probs(double t) const {
    // expm1 avoids cancellation in pDiff at small k t (FAST-EXP-001).
    const double arg = -k_ * t / km1_;
    const double negExpm1 = -std::expm1(arg);
    JCProbs p;
    p.expTerm = 1.0 - negExpm1;
    p.pSame = invK_ + (1.0 - invK_) * p.expTerm;
    p.pDiff = invK_ * negExpm1;
    p.diffCoeff = p.pSame - p.pDiff;
    return p;
  }

 private:
  int k_;
  double invK_;
  double km1_;
};

inline JCProbs jc_probs(int k, double t) { return JCRates(k).probs(t); }

// Neomorphic (asymmetric binary) rates for rate_loss r, with stationary
// frequencies pi0 = r / (1 + r), pi1 = 1 / (1 + r). Scaled so the
// pi-weighted mean rate pi0 * rate01 + pi1 * rate10 is 1, as RevBayes'
// fnFreeK(rescaled = TRUE) does; branch lengths are then in expected changes.
inline void mkn_rates(double rateLoss, double& rate01, double& rate10) {
  const double half = 0.5 * (1.0 + rateLoss);
  rate01 = half / rateLoss;  // gain, 0 -> 1
  rate10 = half;             // loss, 1 -> 0
}

struct MkNProbs {
  double p00, p01, p10, p11;
};

// MkN(rate_loss). Construct once per rate_loss; probs() per branch.
class MkNRates {
 public:
  explicit MkNRates(double rateLoss) : rateLoss_(rateLoss) {
    mkn_rates(rateLoss, rate01_, rate10_);
    lambda_ = rate01_ + rate10_;
    w01_ = rate01_ / lambda_;
    w10_ = rate10_ / lambda_;
  }

  double rate_loss() const { return rateLoss_; }
  double rate01() const { return rate01_; }
  double rate10() const { return rate10_; }
  double lambda() const { return lambda_; }

  // Stationary frequencies, the root distribution of a reversible MkN tree.
  double pi0() const { return rateLoss_ / (1.0 + rateLoss_); }
  double pi1() const { return 1.0 / (1.0 + rateLoss_); }

  MkNProbs probs(double t) const {
    // expm1 avoids cancellation in p01 / p10 at small lambda t (FAST-EXP-001).
    const double arg = -lambda_ * t;
    const double negExpm1 = -std::expm1(arg);
    const double expTerm = 1.0 - negExpm1;
    MkNProbs p;
    p.p00 = w10_ + w01_ * expTerm;
    p.p01 = w01_ * negExpm1;
    p.p10 = w10_ * negExpm1;
    p.p11 = w01_ + w10_ * expTerm;
    return p;
  }

 private:
  double rateLoss_;
  double rate01_;
  double rate10_;
  double lambda_;
  double w01_;
  double w10_;
};

// Kernels mapping a child's conditional likelihoods (CLs) to their
// contribution at the parent: result[c * k + i] = sum_j P_ij(t) cl[c * k + j],
// for nChar characters stored contiguously, k values each.

// JC(k), O(k) per character.
inline void jc_transition(const double* cl, double* result,
                          int nChar, int kStates,
                          double pDiff, double diffCoeff) {
  for (int c = 0; c < nChar; ++c) {
    int off = c * kStates;
    double s = 0.0;
    for (int j = 0; j < kStates; ++j) s += cl[off + j];
    for (int i = 0; i < kStates; ++i)
      result[off + i] = pDiff * s + diffCoeff * cl[off + i];
  }
}

inline void mkn_transition(const double* cl, double* result, int nChar,
                           double p00, double p01, double p10, double p11) {
  for (int c = 0; c < nChar; ++c) {
    int off = c * 2;
    double c0 = cl[off], c1 = cl[off + 1];
    result[off]     = p00 * c0 + p01 * c1;
    result[off + 1] = p10 * c0 + p11 * c1;
  }
}

// F81 with stationary frequencies pi (length kStates) and rate mu:
// (P(t) cl)_i = (1 - e^{-mu t}) * dot(pi, cl) + e^{-mu t} * cl_i.
inline void f81_transition(const double* cl, double* result,
                           int nChar, int kStates,
                           const double* pi, double mu, double t) {
  // expm1 form is exact even at tiny mu * t (FAST-EXP-001).
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

}  // namespace MkPrime

#endif  // MKPRIME_TRANSITION_H
