#ifndef MKPRIME_MKN_RATES_H
#define MKPRIME_MKN_RATES_H

// Neomorphic (asymmetric binary) rates for rate_loss r, with stationary
// frequencies pi0 = r / (1 + r), pi1 = 1 / (1 + r). Scaled so the
// pi-weighted mean rate pi0 * rate01 + pi1 * rate10 is 1, as RevBayes'
// fnFreeK(rescaled = TRUE) does; branch lengths are then in expected changes.
inline void mkn_rates(double rateLoss, double& rate01, double& rate10) {
  const double half = 0.5 * (1.0 + rateLoss);
  rate01 = half / rateLoss;  // gain, 0 -> 1
  rate10 = half;             // loss, 1 -> 0
}

#endif  // MKPRIME_MKN_RATES_H
