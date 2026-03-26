#include <Rcpp.h>
#include <cmath>

// Mk' relabelling correction
//
// When the true number of states is k' but only kObs are observed in the
// data, we need to account for the number of ways to assign the k' labels
// to the kObs observed categories.
//
// The correction (to be ADDED to the log-likelihood) is:
//   log C(k', kObs) = log(kObs!) - log(k'!) + log((k' - kObs)!) + kObs * log(k')
//
// This is the log of:
//   C(k', kObs) = kObs! / k'! * (k' - kObs)! * k'^kObs
//               = k'^kObs / C(k', kObs)  [binomial coefficient]
//
// Derivation: There are C(k', kObs) ways to choose which kObs of the k'
// states are observed. Given that choice, there are kObs! ways to assign
// the observed state labels. So the total is:
//   kObs! * C(k', kObs) * (1/k'^kObs)^{-1}
// Wait — let's be precise. The correction comes from the Stirling number /
// relabelling argument in Mk' (see mkprime-model.md).
//
// Parameters:
//   kPrime: true number of states (k' >= kObs)
//   kObs: observed number of states (>= 2)
//
// Returns: log of the relabelling correction factor

// [[Rcpp::export]]
double mk_prime_relabel_log(int kPrime, int kObs) {
  if (kPrime < kObs) {
    Rcpp::stop("kPrime (%d) must be >= kObs (%d)", kPrime, kObs);
  }
  if (kObs < 1) {
    Rcpp::stop("kObs must be >= 1");
  }

  // Use lgamma(n+1) = log(n!)
  double log_correction =
    std::lgamma(kObs + 1.0) -      // log(kObs!)
    std::lgamma(kPrime + 1.0) +     // -log(k'!)
    std::lgamma(kPrime - kObs + 1.0) + // log((k' - kObs)!)
    kObs * std::log(static_cast<double>(kPrime));  // kObs * log(k')

  return log_correction;
}


// Vectorized version: compute relabelling correction for multiple characters
// with the same kObs but potentially different k' values.
//
// Parameters:
//   kPrime_vec: integer vector of k' values (one per character)
//   kObs: shared observed state count for this batch
//
// Returns: numeric vector of log corrections

// [[Rcpp::export]]
Rcpp::NumericVector mk_prime_relabel_log_batch(
    Rcpp::IntegerVector kPrime_vec, int kObs) {
  int n = kPrime_vec.size();
  Rcpp::NumericVector result(n);

  // Precompute the constant part: log(kObs!)
  double log_kObs_fact = std::lgamma(kObs + 1.0);

  for (int i = 0; i < n; ++i) {
    int kp = kPrime_vec[i];
    if (kp < kObs) {
      Rcpp::stop("kPrime[%d] = %d is less than kObs = %d", i + 1, kp, kObs);
    }
    result[i] = log_kObs_fact -
      std::lgamma(kp + 1.0) +
      std::lgamma(kp - kObs + 1.0) +
      kObs * std::log(static_cast<double>(kp));
  }

  return result;
}
