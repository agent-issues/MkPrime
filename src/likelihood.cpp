#include <Rcpp.h>
#include "../inst/include/MkPrime/pruning.h"
#include "f81.h"
#include <vector>

// R-callable Felsenstein pruning for a single partition, as thin wrappers
// around MkPrime::Pruner (inst/include/MkPrime/pruning.h).
//
// parent, child: edges of the tree (1-based node numbers, as in phylo$edge);
//   the root is the one node that is never a child.
// edge_length: branch lengths, in expected changes per character.
// tip_states: nTip x nChar, 0-indexed states; -1 for missing.
// Returns the log-likelihood summed over characters, or -Inf.

template <class Model>
static double prune(const Model& model, const Rcpp::IntegerVector& parent,
                    const Rcpp::IntegerVector& child,
                    const Rcpp::NumericVector& edge_length,
                    const Rcpp::IntegerMatrix& tip_states,
                    const std::vector<double>& rates) {
  if (tip_states.ncol() == 0) return 0.0;
  MkPrime::Pruner<Model> pruner(model, tip_states.nrow(), tip_states.ncol(),
                                tip_states.begin(), rates);
  pruner.set_tree(parent.begin(), child.begin(), edge_length.begin(),
                  parent.size());
  // Return:
  return pruner.log_likelihood();
}

static std::vector<double> category_rates(const Rcpp::NumericVector& rates) {
  if (rates.size() == 0) Rcpp::stop("`rate_multipliers` must not be empty.");
  // Return:
  return std::vector<double>(rates.begin(), rates.end());
}

// JC(k), with root frequencies `root_freqs` (length kStates).
// [[Rcpp::export]]
double pruning_jc(Rcpp::IntegerVector parent,
                  Rcpp::IntegerVector child,
                  Rcpp::NumericVector edge_length,
                  Rcpp::IntegerMatrix tip_states,
                  int kStates,
                  Rcpp::NumericVector root_freqs) {
  // Return:
  return prune(MkPrime::JCModel(kStates, root_freqs.begin()), parent, child,
               edge_length, tip_states, std::vector<double>(1, 1.0));
}

// JC(kFull) when only states 0..kObs-1 are observed: the unobserved states
// are lumped into one CL column (MkPrime::JCLumpedModel). Gives the same
// log-likelihood as pruning_jc(..., kStates = kFull) with uniform root
// frequencies.
// [[Rcpp::export]]
double pruning_jc_collapsed(Rcpp::IntegerVector parent,
                            Rcpp::IntegerVector child,
                            Rcpp::NumericVector edge_length,
                            Rcpp::IntegerMatrix tip_states,
                            int kFull,
                            int kObs) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("pruning_jc_collapsed requires 1 <= kObs < kFull.");
  }
  // Return:
  return prune(MkPrime::JCLumpedModel(kFull, kObs), parent, child,
               edge_length, tip_states, std::vector<double>(1, 1.0));
}

// MkN (asymmetric binary), with root frequencies `root_freqs` (length 2).
// [[Rcpp::export]]
double pruning_mkn(Rcpp::IntegerVector parent,
                   Rcpp::IntegerVector child,
                   Rcpp::NumericVector edge_length,
                   Rcpp::IntegerMatrix tip_states,
                   double rate_loss,
                   Rcpp::NumericVector root_freqs) {
  // Return:
  return prune(MkPrime::MkNModel(rate_loss, root_freqs.begin()), parent,
               child, edge_length, tip_states, std::vector<double>(1, 1.0));
}

// ACRV variants: site likelihoods averaged over equiprobable rate
// categories, branch lengths scaled by each of `rate_multipliers` in turn.

// [[Rcpp::export]]
double pruning_jc_acrv(Rcpp::IntegerVector parent,
                       Rcpp::IntegerVector child,
                       Rcpp::NumericVector edge_length,
                       Rcpp::IntegerMatrix tip_states,
                       int kStates,
                       Rcpp::NumericVector root_freqs,
                       Rcpp::NumericVector rate_multipliers) {
  // Return:
  return prune(MkPrime::JCModel(kStates, root_freqs.begin()), parent, child,
               edge_length, tip_states, category_rates(rate_multipliers));
}

// [[Rcpp::export]]
double pruning_jc_acrv_collapsed(Rcpp::IntegerVector parent,
                                 Rcpp::IntegerVector child,
                                 Rcpp::NumericVector edge_length,
                                 Rcpp::IntegerMatrix tip_states,
                                 int kFull,
                                 int kObs,
                                 Rcpp::NumericVector rate_multipliers) {
  if (kObs < 1 || kObs >= kFull) {
    Rcpp::stop("pruning_jc_acrv_collapsed requires 1 <= kObs < kFull.");
  }
  // Return:
  return prune(MkPrime::JCLumpedModel(kFull, kObs), parent, child,
               edge_length, tip_states, category_rates(rate_multipliers));
}

// [[Rcpp::export]]
double pruning_mkn_acrv(Rcpp::IntegerVector parent,
                        Rcpp::IntegerVector child,
                        Rcpp::NumericVector edge_length,
                        Rcpp::IntegerMatrix tip_states,
                        double rate_loss,
                        Rcpp::NumericVector root_freqs,
                        Rcpp::NumericVector rate_multipliers) {
  // Return:
  return prune(MkPrime::MkNModel(rate_loss, root_freqs.begin()), parent,
               child, edge_length, tip_states,
               category_rates(rate_multipliers));
}


// Apply the F81 transition kernel to a matrix of conditional likelihoods.
//
// Exposes the kernel the Q-heterogeneity partial-CL path uses
// (gibbs_partial_cl.h) so a test can check it against an independently
// constructed P(t).  Column c of `cl` holds one character's kStates
// conditional likelihoods; the return has the same shape and holds
// (P(t) %*% cl[, c])_i.

// [[Rcpp::export]]
Rcpp::NumericMatrix f81_transition_cpp(Rcpp::NumericMatrix cl,
                                       Rcpp::NumericVector pi,
                                       double mu,
                                       double t) {
  int kStates = cl.nrow();
  int nChar = cl.ncol();
  if (pi.size() != kStates) {
    Rcpp::stop("`pi` must have one entry per row of `cl`.");
  }
  Rcpp::NumericMatrix result(kStates, nChar);
  f81_transition(cl.begin(), result.begin(), nChar, kStates,
                 pi.begin(), mu, t);
  return result;
}
