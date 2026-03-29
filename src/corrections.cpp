#include <Rcpp.h>
#include <cmath>

// Mk' relabelling correction
//
// When the true number of states is k' but only kObs are observed in the
// data, we need to account for the number of ways to assign the k' model-
// state labels to the kObs observed categories.
//
// Under JC(k'), all states are exchangeable, so every injective map from
// {observed labels 0..kObs-1} → {model states 0..k'-1} gives the same
// Felsenstein likelihood.  The number of such maps is the falling factorial:
//
//   P(k', kObs) = k'! / (k' - kObs)!
//
// The correction (ADDED to the Felsenstein log-likelihood) is therefore:
//
//   log P(k', kObs) = log(k'!) - log((k' - kObs)!)
//
// This INCREASES with k' (more states → more equivalent label assignments
// → higher total probability).
//
// When k' = kObs the correction reduces to log(kObs!), a constant that
// cancels in MH ratios but is included for correct absolute log-posteriors.
//
// Verified by direct brute-force enumeration on a star tree (see tests).
//
// Parameters:
//   kPrime: true number of states (k' >= kObs)
//   kObs:   observed number of states (>= 1)
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

  // log P(k', kObs) = lgamma(k'+1) - lgamma(k'-kObs+1)
  return std::lgamma(kPrime + 1.0) - std::lgamma(kPrime - kObs + 1.0);
}


// Vectorized version: compute relabelling correction for multiple characters
// with potentially different k' values (all sharing the same kObs).
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

  double log_denom_base = std::lgamma(1.0); // lgamma(0+1) for the kObs=kPrime case

  for (int i = 0; i < n; ++i) {
    int kp = kPrime_vec[i];
    if (kp < kObs) {
      Rcpp::stop("kPrime[%d] = %d is less than kObs = %d", i + 1, kp, kObs);
    }
    result[i] = std::lgamma(kp + 1.0) - std::lgamma(kp - kObs + 1.0);
  }

  return result;
}
