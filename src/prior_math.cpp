#include <Rcpp.h>
#include "prior_math.h"
using namespace Rcpp;

// Test hooks: expose the C++ helpers so the R tests can pin them directly.
// [[Rcpp::export]]
NumericVector log1m_exp_cpp(NumericVector x) {
  NumericVector out(x.size());
  for (R_xlen_t i = 0; i < x.size(); ++i) out[i] = mkp::log1m_exp(x[i]);
  return out;
}

// [[Rcpp::export]]
double logseries_log_norm_cpp(double c) {
  return mkp::logseries_log_norm(c);
}
