// Fréchet correlation ESS — C++ inner loop
//
// Adapted from the treess package by Andrew F. Magee (GPL-3).
//   Source: https://github.com/afmagee/treess
//   Paper:  Magee et al. (2021), arXiv:2109.07629
//
// The algorithm generalises univariate autocorrelation ESS to
// non-Euclidean spaces via Fréchet variance.  This C++ version
// replaces the R-level loop that repeatedly subsets the distance
// matrix with in-place accumulation over the squared-distance matrix.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double frechet_correlation_ess_cpp(const NumericMatrix& dmat_sq,
                                   int min_nsamples) {
  const int n = dmat_sq.nrow();
  const int max_lag = n - min_nsamples - 1;
  if (max_lag < 1) return NA_REAL;

  // P stores paired sums of consecutive correlations (Geyer's initial
  // positive sequence estimator, as in Vehtari et al.)
  std::vector<double> P;
  P.reserve((max_lag + 1) / 2);

  double prev_cor = 1.0;  // cors[0] = 1.0 by definition

  for (int i = 1; i <= max_lag; ++i) {
    // Fréchet variance of samples [i .. n-1] (front-trimmed)
    // var1 = sum(dmat_sq[i:, i:]) / (2 * m * (m - 1))  where m = n - i
    double sum1 = 0.0;
    const int m = n - i;
    for (int r = i; r < n; ++r) {
      for (int c = i; c < n; ++c) {
        sum1 += dmat_sq(r, c);
      }
    }
    double var1 = sum1 / (2.0 * m * (m - 1));

    // Fréchet variance of samples [0 .. n-i-1] (back-trimmed)
    double sum2 = 0.0;
    for (int r = 0; r < m; ++r) {
      for (int c = 0; c < m; ++c) {
        sum2 += dmat_sq(r, c);
      }
    }
    double var2 = sum2 / (2.0 * m * (m - 1));

    // Mean squared distance at lag i (the i-th super-diagonal)
    double d12 = 0.0;
    for (int j = 0; j < n - i; ++j) {
      d12 += dmat_sq(j, j + i);
    }
    d12 /= (n - i);

    // Lower-bound covariance
    double covar = (var1 + var2 - d12) / 2.0;

    double cor_i;
    if (var1 == 0.0 || var2 == 0.0) {
      cor_i = 1.0;
    } else {
      cor_i = covar / std::sqrt(var1 * var2);
    }

    // Pair consecutive correlations: P_k = cors[2k-1] + cors[2k]
    if (i % 2 == 1) {
      double p_val = cor_i + prev_cor;
      if (p_val < 0.0) break;  // initial positive sequence truncation
      P.push_back(p_val);
    }
    prev_cor = cor_i;
  }

  if (P.empty()) {
    return static_cast<double>(n);  // no autocorrelation detected
  }

  // Monotone (smoothed) sequence: P'[k] = min(P[k], P[k-1])
  for (size_t k = 1; k < P.size(); ++k) {
    if (P[k] > P[k - 1]) P[k] = P[k - 1];
  }

  // The last P may be negative (it's the one that triggered the break
  // at the boundary).  If the final P is positive, include it.
  int K = static_cast<int>(P.size()) - 1;
  if (P.back() > 0.0) K = static_cast<int>(P.size());

  double tau_hat = -1.0;
  for (int k = 0; k < K; ++k) {
    tau_hat += 2.0 * P[k];
  }
  if (tau_hat < 1.0) tau_hat = 1.0;  // floor at 1

  return static_cast<double>(n) / tau_hat;
}
