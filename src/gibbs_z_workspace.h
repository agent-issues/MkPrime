#ifndef MKPRIME_GIBBS_Z_WORKSPACE_H
#define MKPRIME_GIBBS_Z_WORKSPACE_H

// GibbsZWorkspace: pre-allocated scratch buffers for per_char_log_lik_ecology.
//
// Allocate ONE instance at the top of gibbs_z_sweep_impl (sized to the
// maximum stride across all partitions) and pass it by reference into every
// per_char_log_lik_ecology call.  Each call writes into the workspace in-place;
// no heap allocation occurs per call.
//
// Layout mirrors ClWorkspace (src/mcmc_state.h) used by the blind pruner path.

#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cstring>

struct GibbsZWorkspace {
  std::vector<double>  buf;
  std::vector<uint8_t> initFlg;
  Rcpp::IntegerMatrix  tipStates;  // nTip x 1 — reused across calls
  Rcpp::IntegerMatrix  zPart;      // 1 x maxZCols — reused across calls
  int maxNode   = 0;
  int maxStride = 0;
  int curNTip   = 0;
  int curZCols  = 0;

  // Ensure workspace can handle a call needing (maxNode_, stride, nTip, zCols).
  // Resizes only when the current allocation is insufficient.
  void ensure(int maxNode_, int stride, int nTip, int zCols) {
    if (maxNode_ > maxNode || stride > maxStride) {
      maxNode   = std::max(maxNode_, maxNode);
      maxStride = std::max(stride,   maxStride);
      buf.resize(static_cast<size_t>(maxNode + 1) * maxStride, 0.0);
      initFlg.resize(static_cast<size_t>(maxNode + 1), 0u);
    }
    if (nTip != curNTip) {
      curNTip   = nTip;
      tipStates = Rcpp::IntegerMatrix(nTip, 1);
    }
    if (zCols > curZCols) {
      curZCols = zCols;
      zPart    = Rcpp::IntegerMatrix(1, zCols);
    }
  }
};

#endif  // MKPRIME_GIBBS_Z_WORKSPACE_H
