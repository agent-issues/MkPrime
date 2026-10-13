#ifndef MKPRIME_JOINT_ASC_H
#define MKPRIME_JOINT_ASC_H

// Joint ascertainment for transformational characters.
//
// A character's k' belongs to the character, and an invariant character goes
// unrecorded whatever its k'. So the recorded population is the pairs
// (K, y) ~ pi(k) L(y | k) that happen to be variable:
//
//   Pr(K = k, y | V) = pi(k) L(y | k) / D,   D = sum_j pi(j) Pr(V | K = j),
//
// where L(y | k) is the unconditioned pattern probability (relabelling
// factor included). The per-character likelihood therefore carries no
// ascertainment term at its own k'; each character instead contributes
// -log D, shared by every character with the same missing-data mask (and,
// under a prior that conditions on kObs, the same kObs). D depends on the
// tree, the rates and the k' prior, but not on k'.
//
// D needs Pr(V | k) at every k in the prior's support [lo, K]. Pr(V | k) is
// smooth in u = 1 / (k - 1), so k < kJointExactBelow is summed exactly and
// the tail by interpolation at Chebyshev nodes in u, with a check on the
// interpolant's coefficients that falls back to three times as many nodes
// (reusing the first set), then to the direct sum.

#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <memory>
#include <vector>

struct McmcData;
struct JointLayout;

// Hyperparameters of the k' prior: p (geometric, empirical_geometric) and
// alpha, beta (beta_geometric). The logseries prior carries none.
struct KprimeHyper {
  double p = 0.5;
  double alpha = 1.0;
  double beta = 1.0;
};

// k below this value is summed exactly.
constexpr int kJointExactBelow = 8;

// Pr(V | k) on one tree for the missing-data masks of the transformational
// partitions in one rate class, at the k values their normalisers need.
// Pr(V | k) depends on the mask, not the partition, so one traversal serves
// every partition in the class.
struct JointAscTable {
  int lo = 2;                    // smallest k any partition's prior supports
  int K = 2;                     // largest k in the prior's support
  int nMask = 0;
  std::vector<int> row;          // row of g by data.missingMasks id; -1 if none
  // Exact k: integers exactLo .. exactHi; g is nMask x nK, row-major by mask.
  int exactLo = 0, exactHi = -1;
  // Chebyshev nodes in u over [0, 1 / (kJointExactBelow - 1)]: nNode > 0
  // means the k >= kJointExactBelow tail is interpolated, nNode == 0 that
  // every k in [lo, K] is exact.
  int nNode = 0;
  std::vector<double> kVals;
  std::vector<double> g;
};

// Pr(all observed tips in one state | k) under JC(k), for each missing-data
// mask (1 = missing tip; nullptr = none missing) and each k > 1 in kVals,
// averaged over the ACRV rates. out is nMask x nK, row-major by mask.
void lumped_const_probs_jc(const Rcpp::IntegerVector& parent,
                           const Rcpp::IntegerVector& child,
                           const Rcpp::NumericVector& edgeLen, int nTip,
                           const Rcpp::NumericVector& rates,
                           const std::vector<const uint8_t*>& masks,
                           const std::vector<double>& kVals, double* out);

// The same under the Dirichlet-marginal F81 mixture of qHeterogeneity.
void lumped_const_probs_het(const Rcpp::IntegerVector& parent,
                            const Rcpp::IntegerVector& child,
                            const Rcpp::NumericVector& edgeLen, int nTip,
                            const Rcpp::NumericVector& rates,
                            double betaScale, int nBetaCat,
                            const std::vector<const uint8_t*>& masks,
                            const std::vector<double>& kVals, double* out);

// log pi(k) for k = *lo .. data.kprimeTruncK under the data's k' prior, for
// a character observed in kObs states, normalised over that support. *lo is
// 2 for a prior on the full support and kObs for one that conditions on it.
std::vector<double> kprime_log_prior_table(const McmcData& data, int kObs,
                                           const KprimeHyper& hyper, int* lo);

// Pr(V | k) table for the rate class of transformational partition partIdx;
// scaledEdge carries the partition scale and rates the ACRV rates (length
// >= 1). forceDirect evaluates every k exactly.
JointAscTable joint_asc_table(const McmcData& data, int partIdx,
                              const Rcpp::IntegerVector& parent,
                              const Rcpp::IntegerVector& child,
                              const Rcpp::NumericVector& scaledEdge,
                              const Rcpp::NumericVector& rates,
                              double betaScale, bool forceDirect = false);

// sum_m n_m log D_m over partition partIdx's masks, from its class's table
// and its log prior table over [lo, K]; +Inf when some D_m vanishes, so that
// subtracting it gives a -Inf likelihood.
double joint_asc_log(const McmcData& data, int partIdx,
                     const JointAscTable& table,
                     const std::vector<double>& logPrior, int lo);

// Pr(V | k) tables by rate class, for the last two trees each was scored on.
// A move that changes neither the tree nor the rates (a k' or hyperparameter
// update, or a return to the current tree after a rejected proposal) reuses
// a table rather than rebuilding it; one that keeps the topology reuses the
// way the masks share partials on it.
class JointAscCache {
 public:
  const JointAscTable& table(const McmcData& data, int partIdx,
                             const Rcpp::IntegerVector& parent,
                             const Rcpp::IntegerVector& child,
                             const Rcpp::NumericVector& scaledEdge,
                             const Rcpp::NumericVector& rates,
                             double betaScale);

 private:
  struct Slot {
    const McmcData* data = nullptr;
    std::vector<int> parent, child;
    std::vector<double> edge, rates;
    double betaScale = 0.0;
    std::shared_ptr<const JointLayout> layout;
    JointAscTable table;
  };
  std::vector<std::array<Slot, 2>> slots_;  // by classIdx
  std::vector<int> recent_;
};

// The amount partition partIdx's log-likelihood loses to ascertainment under
// the joint form: sum_m n_m log D_m at the given hyperparameters. cache may
// be null.
double trans_joint_asc_log(const McmcData& data, int partIdx,
                           const Rcpp::IntegerVector& parent,
                           const Rcpp::IntegerVector& child,
                           const Rcpp::NumericVector& scaledEdge,
                           const Rcpp::NumericVector& rates,
                           double betaScale, const KprimeHyper& hyper,
                           JointAscCache* cache = nullptr);

#endif  // MKPRIME_JOINT_ASC_H
