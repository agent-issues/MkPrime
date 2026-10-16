#ifndef MKPRIME_EDGE_SHIFTS_H
#define MKPRIME_EDGE_SHIFTS_H

// CL rescaling for the single-pass edge-list kernels (mcmc_likelihood.cpp).
//
// Those kernels visit edges in reverse preorder, so a node is complete when
// the edge above it is read: that is where its CLs are checked and, if
// small, rescaled in place (scaling.h). Shifts are allocated at the first
// rescale, so a traversal that never rescales carries no extra state and
// performs its original arithmetic.

#include "../inst/include/MkPrime/scaling.h"
#include <cstddef>
#include <vector>

class EdgeShifts {
 public:
  EdgeShifts(int maxNode, int nChar) : maxNode_(maxNode), nChar_(nChar) {}

  bool on() const { return !node_.empty(); }

  // Rescales character c of node `ch` (cl points at its k entries); returns
  // the shift, by which the caller scales any sum it holds over them.
  int rescale(double* cl, int k, int ch, int c) {
    if (node_.empty()) {
      node_.assign(static_cast<std::size_t>(maxNode_ + 1) * nChar_, 0);
      site_.assign(nChar_, 0);
    }
    const int s = MkPrime::rescale_cl(cl, k);
    node_[static_cast<std::size_t>(ch) * nChar_ + c] += s;
    return s;
  }

  // After the edge par -> ch: par's shifts are its children's sum so far.
  void propagate(int par, int ch, bool firstChild) {
    int* p = node_.data() + static_cast<std::size_t>(par) * nChar_;
    const int* q = node_.data() + static_cast<std::size_t>(ch) * nChar_;
    if (firstChild) {
      for (int c = 0; c < nChar_; ++c) p[c] = q[c];
    } else {
      for (int c = 0; c < nChar_; ++c) p[c] += q[c];
    }
  }

  int shift(int n, int c) const {
    return node_[static_cast<std::size_t>(n) * nChar_ + c];
  }

  // siteLik += root site likelihood `sl` of character c (shift `rootShift`).
  void add_site(double& siteLik, int c, double sl, int rootShift) {
    MkPrime::add_scaled(siteLik, site_[c], sl, rootShift);
  }

  int site_shift(int c) const { return site_[c]; }

 private:
  int maxNode_;
  int nChar_;
  std::vector<int> node_;  // [node * nChar + c]; empty until a rescale
  std::vector<int> site_;  // shifts of the per-character accumulators
};

#endif  // MKPRIME_EDGE_SHIFTS_H
