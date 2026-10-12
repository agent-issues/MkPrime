#ifndef MKPRIME_MODELS_H
#define MKPRIME_MODELS_H

// Substitution models for Pruner (pruning.h). A model is any copyable type
// with these three members:
//
//   int n_states() const;
//     CL entries per character.
//   void transition(const double* src, double* dst, int nChar,
//                   double t) const;
//     dst = P(t) src, character by character (see transition.h).
//   double root_site_lik(const double* cl) const;
//     Likelihood of one character given its CL at the root.
//
// All four models here are time-reversible, so a rooted tree gives the same
// likelihood as the unrooted tree obtained by merging the root's two edges.

#include "transition.h"
#include <stdexcept>
#include <vector>

namespace MkPrime {

// JC(k): Mk(k), and Mk' with k = k'. Root frequencies default to uniform.
class JCModel {
 public:
  explicit JCModel(int k) : rates_(check_k(k)), rootFreqs_(k, 1.0 / k) {}
  JCModel(int k, const double* rootFreqs)
    : rates_(check_k(k)), rootFreqs_(rootFreqs, rootFreqs + k) {}

  int n_states() const { return rates_.k(); }

  void transition(const double* src, double* dst, int nChar, double t) const {
    const JCProbs p = rates_.probs(t);
    jc_transition(src, dst, nChar, rates_.k(), p.pDiff, p.diffCoeff);
  }

  double root_site_lik(const double* cl) const {
    double siteLik = 0.0;
    for (int s = 0; s < rates_.k(); ++s) siteLik += rootFreqs_[s] * cl[s];
    return siteLik;
  }

 private:
  static int check_k(int k) {
    if (k < 2) throw std::invalid_argument("JCModel: k must be >= 2");
    return k;
  }
  JCRates rates_;
  std::vector<double> rootFreqs_;
};

// JC(kFull) when only kObs < kFull states are observed at any tip. The
// kFull - kObs unobserved states are lumped into one CL column (index kObs),
// so pruning costs O(kObs + 1) rather than O(kFull) per character while
// giving the JC(kFull) likelihood. Tip states must lie in 0..kObs-1.
// CL convention: the lumped column holds the mean over the lumped states.
class JCLumpedModel {
 public:
  JCLumpedModel(int kFull, int kObs)
    : rates_(kFull), kObs_(kObs), nU_(static_cast<double>(kFull - kObs)),
      invKFull_(1.0 / kFull) {
    if (kObs < 1 || kObs >= kFull) {
      throw std::invalid_argument("JCLumpedModel: need 1 <= kObs < kFull");
    }
  }

  int n_states() const { return kObs_ + 1; }

  void transition(const double* src, double* dst, int nChar, double t) const {
    const JCProbs p = rates_.probs(t);
    const int kEff = kObs_ + 1;
    for (int c = 0; c < nChar; ++c) {
      const int off = c * kEff;
      double sumEff = 0.0;
      for (int j = 0; j < kObs_; ++j) sumEff += src[off + j];
      sumEff += nU_ * src[off + kObs_];
      for (int i = 0; i < kEff; ++i)
        dst[off + i] = p.pDiff * sumEff + p.diffCoeff * src[off + i];
    }
  }

  double root_site_lik(const double* cl) const {
    double sumEff = 0.0;
    for (int s = 0; s < kObs_; ++s) sumEff += cl[s];
    sumEff += nU_ * cl[kObs_];
    return invKFull_ * sumEff;
  }

 private:
  JCRates rates_;
  int kObs_;
  double nU_;
  double invKFull_;
};

// MkN (asymmetric binary). Root frequencies default to the stationary ones.
class MkNModel {
 public:
  explicit MkNModel(double rateLoss)
    : rates_(check_rate(rateLoss)), rf0_(rates_.pi0()), rf1_(rates_.pi1()) {}
  MkNModel(double rateLoss, const double* rootFreqs)
    : rates_(check_rate(rateLoss)), rf0_(rootFreqs[0]), rf1_(rootFreqs[1]) {}

  int n_states() const { return 2; }

  void transition(const double* src, double* dst, int nChar, double t) const {
    const MkNProbs p = rates_.probs(t);
    mkn_transition(src, dst, nChar, p.p00, p.p01, p.p10, p.p11);
  }

  double root_site_lik(const double* cl) const {
    return rf0_ * cl[0] + rf1_ * cl[1];
  }

 private:
  static double check_rate(double r) {
    if (!(r > 0.0)) throw std::invalid_argument("MkNModel: rate_loss must be > 0");
    return r;
  }
  MkNRates rates_;
  double rf0_;
  double rf1_;
};

// F81: stationary frequencies pi (also the root distribution), rate mu.
// mu = 1 / (1 - sum(pi^2)) gives one expected change per unit time.
class F81Model {
 public:
  F81Model(const double* pi, int k, double mu) : pi_(pi, pi + k), mu_(mu) {
    if (k < 2) throw std::invalid_argument("F81Model: k must be >= 2");
  }

  int n_states() const { return static_cast<int>(pi_.size()); }

  void transition(const double* src, double* dst, int nChar, double t) const {
    f81_transition(src, dst, nChar, n_states(), pi_.data(), mu_, t);
  }

  double root_site_lik(const double* cl) const {
    double siteLik = 0.0;
    for (int s = 0; s < n_states(); ++s) siteLik += pi_[s] * cl[s];
    return siteLik;
  }

 private:
  std::vector<double> pi_;
  double mu_;
};

}  // namespace MkPrime

#endif  // MKPRIME_MODELS_H
