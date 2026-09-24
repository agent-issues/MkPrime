#ifndef MKPRIME_ASCERTAINMENT_H
#define MKPRIME_ASCERTAINMENT_H

#include <Rcpp.h>
#include <cstdint>
#include <vector>

// Ascertainment probabilities of ascertainment.cpp. A missing-data mask
// (length nTip, or nullptr) flags tips whose state is unknown: they are
// marginalised, so the result is the probability of a constant (or
// singleton) pattern among the observed tips only. The constant-site
// functions take several masks, share one traversal between them and write
// one probability per mask to `out`; no masks means every tip is observed.

void constant_site_probs_jc_impl(const Rcpp::IntegerVector& parent,
                                 const Rcpp::IntegerVector& child,
                                 const Rcpp::NumericVector& edge_length,
                                 int nTip, int kStates,
                                 const Rcpp::NumericVector& rate_multipliers,
                                 const std::vector<const uint8_t*>& masks,
                                 double* out);

double singleton_site_prob_jc_impl(const Rcpp::IntegerVector& parent,
                                   const Rcpp::IntegerVector& child,
                                   const Rcpp::NumericVector& edge_length,
                                   int nTip, int kStates,
                                   const Rcpp::NumericVector& root_freqs,
                                   const Rcpp::NumericVector& rate_multipliers,
                                   const uint8_t* missing);

void constant_site_probs_mkn_impl(const Rcpp::IntegerVector& parent,
                                  const Rcpp::IntegerVector& child,
                                  const Rcpp::NumericVector& edge_length,
                                  int nTip, double rate_loss,
                                  const Rcpp::NumericVector& root_freqs,
                                  const Rcpp::NumericVector& rate_multipliers,
                                  const std::vector<const uint8_t*>& masks,
                                  double* out);

double singleton_site_prob_mkn_impl(const Rcpp::IntegerVector& parent,
                                    const Rcpp::IntegerVector& child,
                                    const Rcpp::NumericVector& edge_length,
                                    int nTip, double rate_loss,
                                    const Rcpp::NumericVector& root_freqs,
                                    const Rcpp::NumericVector& rate_multipliers,
                                    const uint8_t* missing);

#endif  // MKPRIME_ASCERTAINMENT_H
