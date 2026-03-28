// Fréchet correlation ESS — C++ inner loop
//
// Adapted from the treess package by Andrew F. Magee (GPL-3).
//   Source: https://github.com/afmagee/treess
//   Paper:  Magee et al. (2021), arXiv:2109.07629
//
// The algorithm generalises univariate autocorrelation ESS to
// non-Euclidean spaces via Fréchet variance.
//
// Optimisation strategy:
//   1. Exploit matrix symmetry: read only the upper triangle, deriving
//      lower-triangle column sums via M[r,c] = M[c,r].  This halves
//      memory traffic (the n x n distance matrix is the largest object).
//   2. Precompute column partial sums and super-diagonal sums in one
//      fused O(n^2/2) pass with column-contiguous access.
//   3. Main lag loop is O(1) per iteration (lookup only).
//   Overall: O(n^2/2 + L).

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
double frechet_correlation_ess_cpp(const NumericMatrix& dmat_sq,
                                   int min_nsamples) {
  const int n = dmat_sq.nrow();
  const int max_lag = n - min_nsamples - 1;
  if (max_lag < 1) return NA_REAL;

  const double* const d = REAL(dmat_sq);

  // ── Upper-triangle precomputation: O(n²/2), column-contiguous ────
  //
  // For symmetric M, column c's lower triangle = column c's upper
  // triangle transposed.  We read only rows 0..c of each column and
  // accumulate the cross-contribution via lower_partial[r] += M[r,c]
  // (which, by symmetry, equals M[c,r] = the missing lower element).

  std::vector<double> upper_cs(n);             // sum(rows 0..c, col c)
  std::vector<double> lower_partial(n, 0.0);   // sum(rows c+1..n-1, col c)

  const int n_diags = std::min(max_lag + 1, n);
  std::vector<double> diag_sums(n_diags, 0.0); // super-diagonal sums

  for (int c = 0; c < n; ++c) {
    const double* col = d + static_cast<std::size_t>(c) * n;
    const int r_diag = std::max(0, c - n_diags + 1);

    double ucs = 0.0;

    // Region 1: rows where lag (c-r) >= n_diags — no diagonal contribution.
    // Vectorisable: scalar reduction + element-wise accumulation.
    for (int r = 0; r < r_diag; ++r) {
      const double val = col[r];
      ucs += val;
      lower_partial[r] += val;
    }

    // Region 2: rows where lag < n_diags — also accumulate diag_sums.
    for (int r = r_diag; r < c; ++r) {
      const double val = col[r];
      ucs += val;
      lower_partial[r] += val;
      diag_sums[c - r] += val;
    }

    ucs += col[c];  // diagonal element
    upper_cs[c] = ucs;
  }

  // Derive lower_cs and total_sum from the upper-triangle accumulators.
  //   lower_cs[c] = M[c,c] + lower_partial[c]
  //   total_sum    = Σ (upper_cs[c] + lower_partial[c])
  double total_sum = 0.0;
  std::vector<double> lower_cs(n);
  for (int c = 0; c < n; ++c) {
    const double diag = d[static_cast<std::size_t>(c) * n + c];
    lower_cs[c] = diag + lower_partial[c];
    total_sum += upper_cs[c] + lower_partial[c];
  }

  // ── Main lag loop: O(1) per iteration ────────────────────────────

  double trailing = total_sum;  // T(0) = full matrix sum
  double leading  = total_sum;  // S(n) = full matrix sum

  std::vector<double> P;
  P.reserve(static_cast<std::size_t>((max_lag + 1) / 2));
  double prev_cor = 1.0;

  for (int i = 1; i <= max_lag; ++i) {
    const int m = n - i;

    // T(i) = T(i-1) - 2·lower_cs[i-1] + M(i-1,i-1)
    trailing -= 2.0 * lower_cs[i - 1]
              - d[static_cast<std::size_t>(i - 1) * n + (i - 1)];

    // S(m) = S(m+1) - 2·upper_cs[m] + M(m,m)
    leading -= 2.0 * upper_cs[m]
             - d[static_cast<std::size_t>(m) * n + m];

    const double denom = 2.0 * m * (m - 1);
    const double var1 = trailing / denom;
    const double var2 = leading / denom;
    const double d12  = diag_sums[i] / m;

    const double covar = (var1 + var2 - d12) / 2.0;

    double cor_i;
    if (var1 == 0.0 || var2 == 0.0) {
      cor_i = 1.0;
    } else {
      cor_i = covar / std::sqrt(var1 * var2);
    }

    if (i % 2 == 1) {
      const double p_val = cor_i + prev_cor;
      if (p_val < 0.0) break;
      P.push_back(p_val);
    }
    prev_cor = cor_i;
  }

  if (P.empty()) {
    return static_cast<double>(n);
  }

  for (std::size_t k = 1; k < P.size(); ++k) {
    if (P[k] > P[k - 1]) P[k] = P[k - 1];
  }

  int K = static_cast<int>(P.size()) - 1;
  if (P.back() > 0.0) K = static_cast<int>(P.size());

  double tau_hat = -1.0;
  for (int k = 0; k < K; ++k) {
    tau_hat += 2.0 * P[k];
  }
  if (tau_hat < 1.0) tau_hat = 1.0;

  return static_cast<double>(n) / tau_hat;
}
