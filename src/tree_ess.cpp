// Fréchet correlation ESS — C++ inner loop
//
// Adapted from the treess package by Andrew F. Magee (GPL-3).
//   Source: https://github.com/afmagee/treess
//   Paper:  Magee et al. (2021), arXiv:2109.07629
//
// The algorithm generalises univariate autocorrelation ESS to
// non-Euclidean spaces via Fréchet variance.  This C++ version
// replaces the R-level loop that repeatedly subsets the distance
// matrix.  Instead of recomputing O((n-i)^2) submatrix sums at each
// lag, it maintains running trailing and leading submatrix sums
// via O(n)-per-step incremental updates.
//
// Complexity: O(n^2) initial sum + O(n * L) for L lags, versus
// O(n^2 * L) in the naive implementation.

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double frechet_correlation_ess_cpp(const NumericMatrix& dmat_sq,
                                   int min_nsamples) {
  const int n = dmat_sq.nrow();
  const int max_lag = n - min_nsamples - 1;
  if (max_lag < 1) return NA_REAL;

  // R matrices are column-major: dmat_sq(r, c) = data[c * n + r].
  // Access with varying r (fixed c) is contiguous; varying c is stride-n.
  // For a symmetric matrix, row_sum(r, a..b) = col_sum(a..b, r),
  // so we prefer iterating over the column dimension.

  // Total matrix sum: O(n^2), column-contiguous access.
  double total_sum = 0.0;
  for (int c = 0; c < n; ++c) {
    for (int r = 0; r < n; ++r) {
      total_sum += dmat_sq(r, c);
    }
  }

  // Running submatrix sums, updated incrementally at each lag.
  //   trailing = T(i) = sum(dmat_sq[i..n-1, i..n-1])
  //   leading  = S(m) = sum(dmat_sq[0..m-1, 0..m-1])  where m = n - i
  double trailing = total_sum;  // T(0)
  double leading  = total_sum;  // S(n)

  std::vector<double> P;
  P.reserve((max_lag + 1) / 2);
  double prev_cor = 1.0;  // cors[0] = 1.0 by definition

  for (int i = 1; i <= max_lag; ++i) {
    const int m = n - i;

    // T(i) = T(i-1) - 2 * row_sum(i-1, i-1..n-1) + dmat_sq(i-1, i-1)
    // By symmetry, row_sum = col_sum(i-1..n-1, i-1): contiguous access.
    {
      double col_sum = 0.0;
      for (int r = i - 1; r < n; ++r) {
        col_sum += dmat_sq(r, i - 1);
      }
      trailing -= 2.0 * col_sum - dmat_sq(i - 1, i - 1);
    }

    // S(m) = S(m+1) - 2 * row_sum(m, 0..m) + dmat_sq(m, m)
    // By symmetry, row_sum = col_sum(0..m, m): contiguous access.
    {
      double col_sum = 0.0;
      for (int r = 0; r <= m; ++r) {
        col_sum += dmat_sq(r, m);
      }
      leading -= 2.0 * col_sum - dmat_sq(m, m);
    }

    // Super-diagonal sum at lag i (stride-(n+1) access; unavoidable)
    double d12_sum = 0.0;
    for (int j = 0; j < m; ++j) {
      d12_sum += dmat_sq(j, j + i);
    }

    const double denom = 2.0 * m * (m - 1);
    double var1 = trailing / denom;
    double var2 = leading / denom;
    double d12 = d12_sum / m;

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
