// Tree-topology ESS — C++ inner loops
//
// frechet_correlation_ess_cpp():
//   Adapted from the treess package by Andrew F. Magee (GPL-3).
//   Source: https://github.com/afmagee/treess
//   Paper:  Magee et al. (2021), arXiv:2109.07629
//
//   The algorithm generalises univariate autocorrelation ESS to
//   non-Euclidean spaces via Fréchet variance.
//
//   Optimisation strategy:
//     1. Exploit matrix symmetry: read only the upper triangle, deriving
//        lower-triangle column sums via M[r,c] = M[c,r].  This halves
//        memory traffic (the n x n distance matrix is the largest object).
//     2. Precompute column partial sums and super-diagonal sums in one
//        fused O(n^2/2) pass with column-contiguous access.
//     3. Main lag loop is O(1) per iteration (lookup only).
//     Overall: O(n^2/2 + L).
//
// median_pseudo_ess_cpp():
//   Computes the Lanfear et al. (2016) median pseudo-ESS using the
//   Geyer (1992) initial-monotone-sequence estimator instead of
//   coda::effectiveSize (which fits an AR model per row).  This gives
//   a ~100x speedup for large matrices.

#include <Rcpp.h>
#include <cmath>
#include <algorithm>
#include <vector>
#include <numeric>
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


// ── Univariate initial-monotone-sequence ESS (Geyer 1992) ─────────
//
// For a univariate time series x[0..n-1]:
//   1. Compute autocovariances gamma(k)
//   2. Form consecutive pairs P_m = gamma(2m) + gamma(2m+1)
//   3. Apply initial-positive-sequence: truncate when P_m < 0
//   4. Apply monotone-sequence: enforce P_m is non-increasing
//   5. tau = -gamma(0) + 2 * sum(P_m), clamped to >= 1
//   6. ESS = n / tau

// Compute a single autocovariance gamma(k) = (1/n) sum (x[i]-mu)(x[i+k]-mu)
static inline double autocovariance(const double* x, int n, double mean,
                                     int k) {
  double s = 0.0;
  for (int i = 0; i < n - k; ++i) {
    s += (x[i] - mean) * (x[i + k] - mean);
  }
  return s / n;
}

static double univariate_ess_geyer(const double* x, int n,
                                    int min_nsamples) {
  if (n < min_nsamples + 2) return NA_REAL;

  double mean = 0.0;
  for (int i = 0; i < n; ++i) mean += x[i];
  mean /= n;

  // A constant row means the chain never left one topology. Within one chain
  // that cannot be told apart from a posterior concentrated on one topology,
  // so it is unassessable, not maximally mixed.
  const double var0 = autocovariance(x, n, mean, 0);
  if (var0 == 0.0) return NA_REAL;

  // Lazy evaluation: compute autocovariances on demand and stop as soon as
  // the initial-positive-sequence criterion is met.  For chains with ESS ~E,
  // this evaluates ~2E lags instead of all n-min_nsamples lags.
  const int max_lag = n - min_nsamples - 1;

  std::vector<double> P;
  P.reserve(std::min(64, (max_lag + 1) / 2));

  for (int m = 0; 2 * m + 1 <= max_lag; ++m) {
    double rho_even = autocovariance(x, n, mean, 2 * m) / var0;
    double rho_odd  = autocovariance(x, n, mean, 2 * m + 1) / var0;
    double p = rho_even + rho_odd;
    if (p < 0.0) break;
    P.push_back(p);
  }

  if (P.empty()) return static_cast<double>(n);

  // Monotone-sequence correction
  for (std::size_t k = 1; k < P.size(); ++k) {
    if (P[k] > P[k - 1]) P[k] = P[k - 1];
  }

  int K = static_cast<int>(P.size()) - 1;
  if (P.back() > 0.0) K = static_cast<int>(P.size());

  double tau = -1.0;
  for (int k = 0; k < K; ++k) {
    tau += 2.0 * P[k];
  }
  if (tau < 1.0) tau = 1.0;

  return static_cast<double>(n) / tau;
}


// [[Rcpp::export]]
double median_pseudo_ess_cpp(const NumericMatrix& dmat,
                              int min_nsamples,
                              int max_rows) {
  const int n_rows = dmat.nrow();   // number of anchor trees (rows to eval)
  const int n_cols = dmat.ncol();   // chain length (series length per row)

  if (n_cols < min_nsamples + 2) return NA_REAL;

  // Determine which rows to use.  For rectangular matrices from the
  // cross-distance API (200 x n), n_rows is already small; use all.
  std::vector<int> rows;
  if (max_rows >= n_rows || max_rows <= 0) {
    rows.resize(n_rows);
    std::iota(rows.begin(), rows.end(), 0);
  } else {
    // Evenly-spaced subsample (deterministic, reproducible)
    rows.resize(max_rows);
    for (int i = 0; i < max_rows; ++i) {
      rows[i] = static_cast<int>(
        std::round(static_cast<double>(i) * (n_rows - 1) / (max_rows - 1)));
    }
  }

  // Compute per-row ESS.  Each row is a time series of length n_cols.
  const double* d = REAL(dmat);
  std::vector<double> ess_vals;
  ess_vals.reserve(rows.size());

  for (int row : rows) {
    // Column-major: element [row, j] is at d[row + j * n_rows]
    std::vector<double> row_data(n_cols);
    for (int j = 0; j < n_cols; ++j) {
      row_data[j] = d[row + static_cast<std::size_t>(j) * n_rows];
    }

    double e = univariate_ess_geyer(row_data.data(), n_cols, min_nsamples);
    if (!ISNA(e)) ess_vals.push_back(e);
  }

  if (ess_vals.empty()) return NA_REAL;

  // Median
  std::size_t mid = ess_vals.size() / 2;
  std::nth_element(ess_vals.begin(), ess_vals.begin() + mid, ess_vals.end());
  if (ess_vals.size() % 2 == 0) {
    double upper = ess_vals[mid];
    double lower = *std::max_element(ess_vals.begin(), ess_vals.begin() + mid);
    return (lower + upper) / 2.0;
  }
  return ess_vals[mid];
}
