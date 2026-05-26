#include <Rcpp.h>
#include "fast_exp.h"
#include <cmath>
#include <vector>

// JC(k') transition probability matrix — analytical formula
//
// For a JC model with k states (equal rates), the transition probability is:
//   P_ij(t) = (1/k) + ((k-1)/k) * exp(-k*t/(k-1))   for i == j
//   P_ij(t) = (1/k) - (1/k)     * exp(-k*t/(k-1))   for i != j
//
// The rate matrix Q has Q_ii = -1, Q_ij = 1/(k-1), so the expected number
// of substitutions per unit time is 1 (rate is already normalized).
//
// Parameters:
//   k: number of states (>= 2)
//   t: branch length (expected substitutions)
//
// Returns: k x k matrix of transition probabilities (column-major)

// [[Rcpp::export]]
Rcpp::NumericMatrix jc_transition_probs(int k, double t) {
  if (k < 2) {
    Rcpp::stop("k must be >= 2");
  }
  if (t < 0.0) {
    Rcpp::stop("Branch length t must be >= 0");
  }

  Rcpp::NumericMatrix P(k, k);
  double inv_k = 1.0 / k;
  // FAST-EXP-001: use expm1 to avoid catastrophic cancellation in (1-exp) at small kt.
  double arg = -k * t / (k - 1.0);
  double neg_expm1 = -std::expm1(arg);     // = 1 - exp(arg)
  double exp_term  = 1.0 - neg_expm1;       // = exp(arg)
  double diag = inv_k + (1.0 - inv_k) * exp_term;
  double off_diag = inv_k * neg_expm1;

  for (int i = 0; i < k; ++i) {
    for (int j = 0; j < k; ++j) {
      P(i, j) = (i == j) ? diag : off_diag;
    }
  }
  return P;
}


// MkN (asymmetric binary) transition probability matrix
//
// Q-matrix parameterized by rate_loss:
//   rate01 = 2 / (1 + rate_loss)        (gain rate: 0 → 1)
//   rate10 = 2 * rate_loss / (1 + rate_loss)  (loss rate: 1 → 0)
//
// This parameterization normalizes: π0 * rate01 + π1 * rate10 = 1
// where π is the stationary distribution.
//
// Analytical P(t) for 2-state:
//   λ = rate01 + rate10
//   P_00(t) = rate10/λ + rate01/λ * exp(-λ*t)
//   P_01(t) = rate01/λ - rate01/λ * exp(-λ*t)
//   P_10(t) = rate10/λ - rate10/λ * exp(-λ*t)
//   P_11(t) = rate01/λ + rate10/λ * exp(-λ*t)
//
// Parameters:
//   rate_loss: asymmetry parameter (1.0 = symmetric Mk2)
//   t: branch length

// [[Rcpp::export]]
Rcpp::NumericMatrix mkn_transition_probs(double rate_loss, double t) {
  if (rate_loss <= 0.0) {
    Rcpp::stop("rate_loss must be > 0");
  }
  if (t < 0.0) {
    Rcpp::stop("Branch length t must be >= 0");
  }

  double sum_rl = 1.0 + rate_loss;
  double rate01 = 2.0 / sum_rl;      // gain
  double rate10 = 2.0 * rate_loss / sum_rl;  // loss
  double lambda = rate01 + rate10;    // = 2.0 always (by construction)
  // FAST-EXP-001: expm1 form avoids cancellation in P01/P10 at small lambda*t.
  double arg = -lambda * t;
  double neg_expm1 = -std::expm1(arg);
  double exp_term  = 1.0 - neg_expm1;

  double inv_lambda_01 = rate01 / lambda;
  double inv_lambda_10 = rate10 / lambda;

  Rcpp::NumericMatrix P(2, 2);
  P(0, 0) = inv_lambda_10 + inv_lambda_01 * exp_term;  // P_00
  P(0, 1) = inv_lambda_01 * neg_expm1;                 // P_01
  P(1, 0) = inv_lambda_10 * neg_expm1;                 // P_10
  P(1, 1) = inv_lambda_01 + inv_lambda_10 * exp_term;  // P_11

  return P;
}


// Stationary frequencies for the MkN model
//
// π_0 = rate10 / (rate01 + rate10) = rate_loss / (1 + rate_loss)
// π_1 = rate01 / (rate01 + rate10) = 1 / (1 + rate_loss)

// [[Rcpp::export]]
Rcpp::NumericVector mkn_stationary_freqs(double rate_loss) {
  if (rate_loss <= 0.0) {
    Rcpp::stop("rate_loss must be > 0");
  }
  double sum_rl = 1.0 + rate_loss;
  return Rcpp::NumericVector::create(rate_loss / sum_rl, 1.0 / sum_rl);
}
