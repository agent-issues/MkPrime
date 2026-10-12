#include "joint_asc.h"
#include "mcmc_state.h"
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <map>
#include <unordered_map>
#include <vector>

using namespace Rcpp;

namespace {

// Chebyshev node counts tried in turn before the direct sum.
constexpr int kChebLevels[] = {12, 32};
constexpr double kChebRelTol = 1e-12;
constexpr double kChebAbsTol = 1e-15;

inline double cheb_width() { return 1.0 / (kJointExactBelow - 1.0); }

// Chebyshev points of the first kind on [0, cheb_width()], in u = 1/(k - 1).
std::vector<double> cheb_nodes_u(int n) {
  std::vector<double> u(n);
  const double h = cheb_width();
  for (int j = 0; j < n; ++j) {
    u[j] = 0.5 * h * (1.0 + std::cos((2.0 * j + 1.0) * M_PI / (2.0 * n)));
  }
  return u;
}

// Masks that agree on which tips below a node are missing share one lumped
// partial there. Children fold into their parent one at a time, so a node may
// have any number of children: the first child's classes are copied, and
// each later child combines the parent's running class foldA[j] with its own
// class foldB[j] into running class j. A tip has two nominal classes
// (observed, missing) whether or not each occurs, so a fold can leave a node
// fewer classes than the copy before it: cap holds the most it ever has.
struct MaskClasses {
  int nTip = 0, maxNode = 0, root = 0;
  std::vector<int> nCls;                  // by node, after its last fold
  std::vector<int> cap;                   // by node
  std::vector<int> rootCls;               // by mask
  std::vector<uint8_t> first;             // by edge
  std::vector<int> nOut;                  // by edge: parent classes after it
  std::vector<std::vector<int>> foldA, foldB;  // by edge
};

MaskClasses build_mask_classes(const IntegerVector& parent,
                               const IntegerVector& child, int nTip,
                               const std::vector<const uint8_t*>& masks) {
  const int nEdge = parent.size();
  const int nMask = std::max<int>(1, masks.size());
  MaskClasses mc;
  mc.nTip = nTip;
  mc.root = nTip + 1;
  mc.maxNode = 2 * nTip - 1;
  for (int e = 0; e < nEdge; ++e) {
    mc.maxNode = std::max({mc.maxNode, (int)parent[e], (int)child[e]});
  }
  std::vector<std::vector<int>> cls(mc.maxNode + 1);
  mc.nCls.assign(mc.maxNode + 1, 0);
  mc.cap.assign(mc.maxNode + 1, 0);
  for (int t = 1; t <= nTip; ++t) {
    cls[t].resize(nMask);
    for (int m = 0; m < nMask; ++m) {
      cls[t][m] = (!masks.empty() && masks[m] && masks[m][t - 1]) ? 1 : 0;
    }
    mc.nCls[t] = 2;
    mc.cap[t] = 2;
  }
  mc.first.assign(nEdge, 0);
  mc.nOut.assign(nEdge, 0);
  mc.foldA.assign(nEdge, {});
  mc.foldB.assign(nEdge, {});
  std::vector<uint8_t> seen(mc.maxNode + 1, 0);
  std::vector<int> stamp, slot;
  int stampNow = 0;
  std::unordered_map<long long, int> ids;
  for (int e = nEdge - 1; e >= 0; --e) {
    const int par = parent[e], ch = child[e];
    if (!seen[par]) {
      seen[par] = 1;
      cls[par] = cls[ch];
      mc.nCls[par] = mc.nCls[ch];
      mc.cap[par] = mc.nCls[ch];
      mc.first[e] = 1;
      mc.nOut[e] = mc.nCls[ch];
      continue;
    }
    std::vector<int>& cp = cls[par];
    const std::vector<int>& cc = cls[ch];
    std::vector<int>& fa = mc.foldA[e];
    std::vector<int>& fb = mc.foldB[e];
    const long long nB = mc.nCls[ch];
    const long long nKey = (long long)mc.nCls[par] * nB;
    int nNew = 0;
    if (nKey <= 4096) {
      // Small key space: a stamped dense lookup avoids hashing.
      if ((long long)stamp.size() < nKey) {
        stamp.resize(nKey, 0);
        slot.resize(nKey, 0);
      }
      ++stampNow;
      for (int m = 0; m < nMask; ++m) {
        const long long key = cp[m] * nB + cc[m];
        if (stamp[key] != stampNow) {
          stamp[key] = stampNow;
          slot[key] = nNew++;
          fa.push_back(cp[m]);
          fb.push_back(cc[m]);
        }
        cp[m] = slot[key];
      }
    } else {
      ids.clear();
      for (int m = 0; m < nMask; ++m) {
        const long long key = cp[m] * nB + cc[m];
        auto found = ids.find(key);
        if (found == ids.end()) {
          found = ids.emplace(key, nNew++).first;
          fa.push_back(cp[m]);
          fb.push_back(cc[m]);
        }
        cp[m] = found->second;
      }
    }
    mc.nCls[par] = nNew;
    mc.cap[par] = std::max(mc.cap[par], nNew);
    mc.nOut[e] = nNew;
  }
  mc.rootCls = cls[mc.root];
  return mc;
}

// Offsets of each node's classes in a flat per-class buffer.
std::vector<int> class_offsets(const MaskClasses& mc, int* total) {
  std::vector<int> off(mc.maxNode + 2, 0);
  for (int v = 1; v <= mc.maxNode; ++v) off[v + 1] = off[v] + mc.cap[v];
  *total = off[mc.maxNode + 1];
  return off;
}

// JC(k), lumped: at each node, x is the likelihood of the observed tips below
// all being in state 0 given the node is in state 0, y the same given the
// node is in some other (fixed) state. A constant pattern then has
// probability x + (k - 1) y at the root, summed over its k states.
void lumped_jc_values(const MaskClasses& mc, const IntegerVector& parent,
                      const IntegerVector& child, const NumericVector& edgeLen,
                      const NumericVector& rates, const double* kVals, int nK,
                      int nMask, double* out, int outStride) {
  const int nEdge = parent.size();
  const int nCat = rates.size();
  int total;
  const std::vector<int> off = class_offsets(mc, &total);
  std::vector<double> x(total), y(total);
  int maxCls = 0;
  for (int v = 1; v <= mc.maxNode; ++v) maxCls = std::max(maxCls, mc.cap[v]);
  std::vector<double> mx(maxCls), my(maxCls), tx(maxCls), ty(maxCls);
  for (int t = 1; t <= mc.nTip; ++t) {
    x[off[t]] = 1.0;
    y[off[t]] = 0.0;
    x[off[t] + 1] = 1.0;
    y[off[t] + 1] = 1.0;
  }
  for (int ki = 0; ki < nK; ++ki) {
    const double k = kVals[ki];
    const double invK = 1.0 / k;
    const double lam = k / (k - 1.0);
    for (int m = 0; m < nMask; ++m) out[m * outStride + ki] = 0.0;
    for (int c = 0; c < nCat; ++c) {
      const double rate = rates[c];
      for (int e = nEdge - 1; e >= 0; --e) {
        const int ch = child[e], par = parent[e];
        const double omE = -std::expm1(-lam * edgeLen[e] * rate);
        const double a = (1.0 - invK) * omE;
        const double b = invK * omE;
        const double* xc = &x[off[ch]];
        const double* yc = &y[off[ch]];
        const int nc = mc.nCls[ch];
        for (int j = 0; j < nc; ++j) {
          const double d = yc[j] - xc[j];
          mx[j] = xc[j] + a * d;
          my[j] = yc[j] - b * d;
        }
        double* xp = &x[off[par]];
        double* yp = &y[off[par]];
        const int nOut = mc.nOut[e];
        if (mc.first[e]) {
          std::copy_n(mx.data(), nOut, xp);
          std::copy_n(my.data(), nOut, yp);
        } else {
          const int* fa = mc.foldA[e].data();
          const int* fb = mc.foldB[e].data();
          for (int j = 0; j < nOut; ++j) {
            tx[j] = xp[fa[j]] * mx[fb[j]];
            ty[j] = yp[fa[j]] * my[fb[j]];
          }
          std::copy_n(tx.data(), nOut, xp);
          std::copy_n(ty.data(), nOut, yp);
        }
      }
      const double* xr = &x[off[mc.root]];
      const double* yr = &y[off[mc.root]];
      for (int m = 0; m < nMask; ++m) {
        const int j = mc.rootCls[m];
        out[m * outStride + ki] += (xr[j] + (k - 1.0) * yr[j]) / nCat;
      }
    }
  }
}

// Bins of Beta(alpha, (k - 1) alpha), as compute_het_bins in
// mcmc_likelihood.cpp but at real k. Cached: the qbeta and pbeta calls cost
// more than the traversal they feed, and alpha changes only on its own move.
const std::vector<double>& het_bins(double alpha, double k, int nBins) {
  static double cacheAlpha = -1.0;
  static int cacheBins = -1;
  static std::map<double, std::vector<double>> cache;
  if (alpha != cacheAlpha || nBins != cacheBins) {
    cache.clear();
    cacheAlpha = alpha;
    cacheBins = nBins;
  }
  auto found = cache.find(k);
  if (found != cache.end()) return found->second;
  std::vector<double> bins(nBins);
  const double a = alpha;
  const double b = (k - 1.0) * alpha;
  for (int i = 0; i < nBins; ++i) {
    const double lo = R::qbeta((double)i / nBins, a, b, 1, 0);
    const double hi = R::qbeta((double)(i + 1) / nBins, a, b, 1, 0);
    if (hi - lo < 1e-15) {
      bins[i] = 0.5 * (lo + hi);
      continue;
    }
    const double pLo = R::pbeta(lo, a + 1.0, b, 1, 0);
    const double pHi = R::pbeta(hi, a + 1.0, b, 1, 0);
    const double denom = R::pbeta(hi, a, b, 1, 0) - R::pbeta(lo, a, b, 1, 0);
    bins[i] = denom < 1e-300 ? 0.5 * (lo + hi)
                             : (a / (a + b)) * (pHi - pLo) / denom;
  }
  return cache.emplace(k, std::move(bins)).first->second;
}

// F81 with one state ("the spike") at frequency beta and k - 1 at
// r = (1 - beta) / (k - 1), lumped. Constant patterns in the spike state (A)
// need its value at a spike node (a0) and at any other node (a1); those in a
// non-spike target state (B) need it at the spike (b0), the target (b1) and
// any of the k - 2 other states (b2). Averaged over ACRV rates and bins.
void lumped_het_values(const MaskClasses& mc, const IntegerVector& parent,
                       const IntegerVector& child, const NumericVector& edgeLen,
                       const NumericVector& rates, double betaScale,
                       int nBetaCat, const double* kVals, int nK, int nMask,
                       double* out, int outStride) {
  const int nEdge = parent.size();
  const int nCat = rates.size();
  int total;
  const std::vector<int> off = class_offsets(mc, &total);
  constexpr int nV = 5;
  std::vector<double> v(total * nV);
  int maxCls = 0;
  for (int n = 1; n <= mc.maxNode; ++n) maxCls = std::max(maxCls, mc.cap[n]);
  std::vector<double> msg(maxCls * nV), tmp(maxCls * nV);
  for (int t = 1; t <= mc.nTip; ++t) {
    double* obs = &v[off[t] * nV];
    obs[0] = 1.0; obs[1] = 0.0; obs[2] = 0.0; obs[3] = 1.0; obs[4] = 0.0;
    std::fill_n(obs + nV, nV, 1.0);
  }
  const double norm = 1.0 / (nCat * nBetaCat);
  for (int ki = 0; ki < nK; ++ki) {
    const double k = kVals[ki];
    const std::vector<double>& bins = het_bins(betaScale, k, nBetaCat);
    for (int m = 0; m < nMask; ++m) out[m * outStride + ki] = 0.0;
    for (int bi = 0; bi < nBetaCat; ++bi) {
      const double beta = bins[bi];
      const double r = (1.0 - beta) / (k - 1.0);
      const double mu = 1.0 / (1.0 - beta * beta - (k - 1.0) * r * r);
      for (int c = 0; c < nCat; ++c) {
        const double rate = rates[c];
        for (int e = nEdge - 1; e >= 0; --e) {
          const int ch = child[e], par = parent[e];
          const double omd = -std::expm1(-mu * edgeLen[e] * rate);
          const double d = 1.0 - omd;
          const double* vc = &v[off[ch] * nV];
          const int nc = mc.nCls[ch];
          for (int j = 0; j < nc; ++j) {
            const double* s = vc + j * nV;
            double* o = &msg[j * nV];
            const double sA = omd * (beta * s[0] + (1.0 - beta) * s[1]);
            const double sB = omd * (beta * s[2] + r * s[3]
                                     + (k - 2.0) * r * s[4]);
            o[0] = d * s[0] + sA;
            o[1] = d * s[1] + sA;
            o[2] = d * s[2] + sB;
            o[3] = d * s[3] + sB;
            o[4] = d * s[4] + sB;
          }
          double* vp = &v[off[par] * nV];
          const int nOut = mc.nOut[e];
          if (mc.first[e]) {
            std::copy_n(msg.data(), nOut * nV, vp);
          } else {
            const int* fa = mc.foldA[e].data();
            const int* fb = mc.foldB[e].data();
            for (int j = 0; j < nOut; ++j) {
              for (int q = 0; q < nV; ++q) {
                tmp[j * nV + q] = vp[fa[j] * nV + q] * msg[fb[j] * nV + q];
              }
            }
            std::copy_n(tmp.data(), nOut * nV, vp);
          }
        }
        const double* vr = &v[off[mc.root] * nV];
        for (int m = 0; m < nMask; ++m) {
          const double* s = vr + mc.rootCls[m] * nV;
          const double pA = beta * s[0] + (1.0 - beta) * s[1];
          const double pB = beta * s[2] + r * s[3] + (k - 2.0) * r * s[4];
          out[m * outStride + ki] += norm * (pA + (k - 1.0) * pB);
        }
      }
    }
  }
}

// Pr(V | k) = 1 - Pr(constant | k) at kVals, appended to g (nMask rows).
void append_variable_probs(const MaskClasses& mc, const McmcData& data,
                           const IntegerVector& parent,
                           const IntegerVector& child,
                           const NumericVector& scaledEdge,
                           const NumericVector& rates, double betaScale,
                           const std::vector<double>& kNew, int nMask,
                           JointAscTable& table) {
  if (kNew.empty()) return;
  const int nOld = table.kVals.size();
  const int nK = nOld + kNew.size();
  std::vector<double> g(nMask * nK);
  for (int m = 0; m < nMask; ++m) {
    std::copy_n(&table.g[m * nOld], nOld, &g[m * nK]);
  }
  if (data.qHeterogeneity) {
    lumped_het_values(mc, parent, child, scaledEdge, rates, betaScale,
                      data.nBetaCat, kNew.data(), kNew.size(), nMask,
                      &g[nOld], nK);
  } else {
    lumped_jc_values(mc, parent, child, scaledEdge, rates, kNew.data(),
                     kNew.size(), nMask, &g[nOld], nK);
  }
  for (int m = 0; m < nMask; ++m) {
    for (int i = nOld; i < nK; ++i) g[m * nK + i] = 1.0 - g[m * nK + i];
  }
  table.kVals.insert(table.kVals.end(), kNew.begin(), kNew.end());
  table.g = std::move(g);
}

// The smallest k' the prior supports for a character observed in kObs
// states: kObs under a prior that conditions on it, else 2.
int prior_lo(const McmcData& data, int kObs) {
  const bool conditional = data.kPriorBetaGeometric ||
    (!data.kPriorLogseries && !data.unconditionalPrior);
  return conditional ? kObs : 2;
}

// Whether every mask's interpolant through its last n values has converged:
// its two highest Chebyshev coefficients are negligible.
bool cheb_converged(const JointAscTable& table, int n) {
  const int nK = table.kVals.size();
  const int nMask = table.nMask;
  std::vector<double> basis(n * n);
  for (int q = 0; q < n; ++q) {
    for (int j = 0; j < n; ++j) {
      basis[q * n + j] = std::cos(q * (2.0 * j + 1.0) * M_PI / (2.0 * n));
    }
  }
  for (int m = 0; m < nMask; ++m) {
    const double* gv = &table.g[m * nK + nK - n];
    double scale = 0.0, tail = 0.0;
    for (int q = 0; q < n; ++q) {
      double cq = 0.0;
      for (int j = 0; j < n; ++j) cq += gv[j] * basis[q * n + j];
      cq *= (q == 0 ? 1.0 : 2.0) / n;
      scale = std::max(scale, std::abs(cq));
      if (q >= n - 2) tail += std::abs(cq);
    }
    if (tail > kChebRelTol * scale + kChebAbsTol) return false;
  }
  return true;
}

}  // namespace


void lumped_const_probs_jc(const IntegerVector& parent,
                           const IntegerVector& child,
                           const NumericVector& edgeLen, int nTip,
                           const NumericVector& rates,
                           const std::vector<const uint8_t*>& masks,
                           const std::vector<double>& kVals, double* out) {
  const MaskClasses mc = build_mask_classes(parent, child, nTip, masks);
  lumped_jc_values(mc, parent, child, edgeLen, rates, kVals.data(),
                   kVals.size(), std::max<int>(1, masks.size()), out,
                   kVals.size());
}

void lumped_const_probs_het(const IntegerVector& parent,
                            const IntegerVector& child,
                            const NumericVector& edgeLen, int nTip,
                            const NumericVector& rates,
                            double betaScale, int nBetaCat,
                            const std::vector<const uint8_t*>& masks,
                            const std::vector<double>& kVals, double* out) {
  const MaskClasses mc = build_mask_classes(parent, child, nTip, masks);
  lumped_het_values(mc, parent, child, edgeLen, rates, betaScale, nBetaCat,
                    kVals.data(), kVals.size(),
                    std::max<int>(1, masks.size()), out, kVals.size());
}


std::vector<double> kprime_log_prior_table(const McmcData& data, int kObs,
                                           const KprimeHyper& hyper, int* lo) {
  const int K = data.kprimeTruncK;
  *lo = prior_lo(data, kObs);
  if (*lo > K) return {};
  const int n = K - *lo + 1;
  std::vector<double> lw(n);
  if (data.kPriorLogseries) {
    const double logC = std::log(data.kprimeLogseriesC);
    for (int i = 0; i < n; ++i) {
      const int k = *lo + i;
      lw[i] = k * logC - std::log((double)k);
    }
  } else if (data.kPriorBetaGeometric) {
    const double a = hyper.alpha, b = hyper.beta;
    for (int i = 0; i < n; ++i) lw[i] = R::lbeta(a + 1.0, b + i);
  } else if (data.kPriorEmpiricalGeometric) {
    // P(k | p) = p P_emp(k) + (1 - p) P(k - 1 | p): the convolution of the
    // empirical pmf with Geometric(p), as cpp_log_prior.
    const double logP = std::log(hyper.p);
    const double log1mP = std::log1p(-hyper.p);
    const int bodyLen = data.empLogBody.size();
    const double logQ = data.empTailDecay > 0.0 ? std::log(data.empTailDecay)
                                                : R_NegInf;
    double prev = R_NegInf;
    for (int k = 2; k <= K; ++k) {
      double logPemp = R_NegInf;
      if (k - 2 < bodyLen) {
        logPemp = data.empLogBody[k - 2];
      } else if (data.empTailStartK > 0 && k >= data.empTailStartK &&
                 std::isfinite(data.empLogTailStartP) && std::isfinite(logQ)) {
        logPemp = data.empLogTailStartP + (k - data.empTailStartK) * logQ;
      }
      const double fresh = logP + logPemp;
      const double carry = log1mP + prev;
      const double mx = std::max(fresh, carry);
      prev = std::isfinite(mx)
        ? mx + std::log(std::exp(fresh - mx) + std::exp(carry - mx))
        : R_NegInf;
      if (k >= *lo) lw[k - *lo] = prev;
    }
  } else {
    const double log1mP = std::log1p(-hyper.p);
    for (int i = 0; i < n; ++i) lw[i] = i * log1mP;
  }
  double mx = R_NegInf;
  for (double w : lw) mx = std::max(mx, w);
  if (!std::isfinite(mx)) {
    std::fill(lw.begin(), lw.end(), R_NegInf);
    return lw;
  }
  double s = 0.0;
  for (double w : lw) s += std::exp(w - mx);
  const double lse = mx + std::log(s);
  for (double& w : lw) w -= lse;
  return lw;
}


JointAscTable joint_asc_table(const McmcData& data, int partIdx,
                              const IntegerVector& parent,
                              const IntegerVector& child,
                              const NumericVector& scaledEdge,
                              const NumericVector& rates,
                              double betaScale, bool forceDirect) {
  const int classIdx = data.parts[partIdx].classIdx;
  JointAscTable table;
  table.K = data.kprimeTruncK;
  table.lo = table.K + 1;
  std::vector<uint8_t> used(data.missingMasks.size(), 0);
  for (const PartInfo& part : data.parts) {
    const int nChar = part.tipStates.ncol();
    if (part.type != 1 || part.classIdx != classIdx || nChar == 0) continue;
    table.lo = std::min(table.lo,
                        prior_lo(data, data.kObs[part.globalCharIdx[0]]));
    for (int c = 0; c < nChar; ++c) {
      used[data.charMask[part.globalCharIdx[c]]] = 1;
    }
  }
  table.row.assign(used.size(), -1);
  std::vector<int> ids;
  for (int id = 0; id < (int)used.size(); ++id) {
    if (!used[id]) continue;
    table.row[id] = ids.size();
    ids.push_back(id);
  }
  table.nMask = ids.size();
  if (table.lo > table.K) return table;

  std::vector<const uint8_t*> masks;
  for (int id : ids) {
    masks.push_back(id == 0 ? nullptr : data.missingMasks[id].data());
  }
  const int nMask = masks.size();
  const MaskClasses mc = build_mask_classes(parent, child, data.nTip, masks);
  const NumericVector useRates = rates.size() ? rates : NumericVector(1, 1.0);

  const int tailLo = std::max(table.lo, kJointExactBelow);
  const int nTail = table.K - tailLo + 1;
  const bool direct = forceDirect || nTail <= kChebLevels[0];
  table.exactLo = table.lo;
  table.exactHi = direct ? table.K
                         : std::min(table.K, kJointExactBelow - 1);
  std::vector<double> kNew;
  for (int k = table.exactLo; k <= table.exactHi; ++k) kNew.push_back(k);
  if (direct) {
    append_variable_probs(mc, data, parent, child, scaledEdge, useRates,
                          betaScale, kNew, nMask, table);
    return table;
  }
  for (int n : kChebLevels) {
    for (double u : cheb_nodes_u(n)) kNew.push_back(1.0 + 1.0 / u);
    append_variable_probs(mc, data, parent, child, scaledEdge, useRates,
                          betaScale, kNew, nMask, table);
    kNew.clear();
    if (cheb_converged(table, n)) {
      // Keep the exact block and this level's nodes only.
      const int nKeep = table.exactHi >= table.exactLo
        ? table.exactHi - table.exactLo + 1 : 0;
      const int nK = table.kVals.size();
      std::vector<double> g(nMask * (nKeep + n));
      for (int m = 0; m < nMask; ++m) {
        std::copy_n(&table.g[m * nK], nKeep, &g[m * (nKeep + n)]);
        std::copy_n(&table.g[m * nK + nK - n], n, &g[m * (nKeep + n) + nKeep]);
      }
      std::vector<double> kv(table.kVals.begin(),
                             table.kVals.begin() + nKeep);
      kv.insert(kv.end(), table.kVals.end() - n, table.kVals.end());
      table.kVals = std::move(kv);
      table.g = std::move(g);
      table.nNode = n;
      return table;
    }
  }
  // No level converged: sum the tail directly.
  const int nKeep = table.exactHi >= table.exactLo
    ? table.exactHi - table.exactLo + 1 : 0;
  const int nK = table.kVals.size();
  std::vector<double> g(nMask * nKeep);
  for (int m = 0; m < nMask; ++m) {
    std::copy_n(&table.g[m * nK], nKeep, &g[m * nKeep]);
  }
  table.kVals.resize(nKeep);
  table.g = std::move(g);
  for (int k = tailLo; k <= table.K; ++k) kNew.push_back(k);
  append_variable_probs(mc, data, parent, child, scaledEdge, useRates,
                        betaScale, kNew, nMask, table);
  table.exactLo = table.lo;
  table.exactHi = table.K;
  return table;
}


double joint_asc_log(const McmcData& data, int partIdx,
                     const JointAscTable& table,
                     const std::vector<double>& logPrior, int lo) {
  if (lo > table.K || logPrior.empty()) return R_PosInf;
  const int nK = table.kVals.size();
  std::vector<double> w(nK, 0.0);
  const int nExact = table.exactHi >= table.exactLo
    ? table.exactHi - table.exactLo + 1 : 0;
  for (int i = 0; i < nExact; ++i) {
    const int k = table.exactLo + i;
    if (k >= lo) w[i] = std::exp(logPrior[k - lo]);
  }
  if (table.nNode > 0) {
    const int n = table.nNode;
    const std::vector<double> u = cheb_nodes_u(n);
    std::vector<double> bw(n);
    for (int j = 0; j < n; ++j) {
      bw[j] = ((j % 2) ? -1.0 : 1.0) *
        std::sin((2.0 * j + 1.0) * M_PI / (2.0 * n));
    }
    double* W = &w[nExact];
    std::vector<double> q(n);
    for (int k = std::max(lo, kJointExactBelow); k <= table.K; ++k) {
      const double pk = std::exp(logPrior[k - lo]);
      if (pk == 0.0) continue;
      const double uk = 1.0 / (k - 1.0);
      double s = 0.0;
      int hit = -1;
      for (int j = 0; j < n; ++j) {
        const double d = uk - u[j];
        if (d == 0.0) { hit = j; break; }
        q[j] = bw[j] / d;
        s += q[j];
      }
      if (hit >= 0) {
        W[hit] += pk;
        continue;
      }
      for (int j = 0; j < n; ++j) W[j] += pk * q[j] / s;
    }
  }
  const PartInfo& part = data.parts[partIdx];
  const MaskTally tally = tally_masks(data, part.globalCharIdx,
                                      part.tipStates.ncol());
  double total = 0.0;
  for (size_t i = 0; i < tally.ids.size(); ++i) {
    const double* gm = &table.g[table.row[tally.ids[i]] * nK];
    double D = 0.0;
    for (int j = 0; j < nK; ++j) D += w[j] * gm[j];
    if (!(D > 0.0)) return R_PosInf;
    total += tally.counts[i] * std::log(D);
  }
  return total;
}


const JointAscTable& JointAscCache::table(const McmcData& data, int partIdx,
                                          const IntegerVector& parent,
                                          const IntegerVector& child,
                                          const NumericVector& scaledEdge,
                                          const NumericVector& rates,
                                          double betaScale) {
  const int cls = data.parts[partIdx].classIdx;
  if ((int)slots_.size() <= cls) {
    slots_.resize(cls + 1);
    recent_.resize(cls + 1, 0);
  }
  const NumericVector useRates = rates.size() ? rates : NumericVector(1, 1.0);
  const double bs = data.qHeterogeneity ? betaScale : 0.0;
  auto same = [](const auto& a, const auto& b) {
    return std::equal(a.begin(), a.end(), b.begin(), b.end());
  };
  std::array<Slot, 2>& pair = slots_[cls];
  const int last = recent_[cls];
  for (int i : {last, 1 - last}) {
    const Slot& s = pair[i];
    if (s.data == &data && s.betaScale == bs && same(parent, s.parent) &&
        same(child, s.child) && same(scaledEdge, s.edge) &&
        same(useRates, s.rates)) {
      recent_[cls] = i;
      return s.table;
    }
  }
  Slot& s = pair[1 - last];
  s.data = &data;
  s.parent.assign(parent.begin(), parent.end());
  s.child.assign(child.begin(), child.end());
  s.edge.assign(scaledEdge.begin(), scaledEdge.end());
  s.rates.assign(useRates.begin(), useRates.end());
  s.betaScale = bs;
  s.table = joint_asc_table(data, partIdx, parent, child, scaledEdge,
                            useRates, bs);
  recent_[cls] = 1 - last;
  // Return:
  return s.table;
}


double trans_joint_asc_log(const McmcData& data, int partIdx,
                           const IntegerVector& parent,
                           const IntegerVector& child,
                           const NumericVector& scaledEdge,
                           const NumericVector& rates,
                           double betaScale, const KprimeHyper& hyper,
                           JointAscCache* cache) {
  const PartInfo& part = data.parts[partIdx];
  int lo;
  const std::vector<double> logPrior = kprime_log_prior_table(
    data, data.kObs[part.globalCharIdx[0]], hyper, &lo);
  if (cache) {
    return joint_asc_log(data, partIdx,
                         cache->table(data, partIdx, parent, child,
                                      scaledEdge, rates, betaScale),
                         logPrior, lo);
  }
  // Return:
  return joint_asc_log(
    data, partIdx,
    joint_asc_table(data, partIdx, parent, child, scaledEdge, rates,
                    betaScale),
    logPrior, lo);
}


// The kernels fold each node's children in reverse edge order, so a node
// must be reached from the root before its own edges are listed.
static void check_preorder(const IntegerVector& parent,
                           const IntegerVector& child, int nTip) {
  const int nEdge = parent.size();
  if (child.size() != nEdge) stop("parent and child differ in length.");
  int maxNode = nTip + 1;
  for (int e = 0; e < nEdge; ++e) {
    maxNode = std::max({maxNode, (int)parent[e], (int)child[e]});
    if (parent[e] <= nTip || child[e] < 1) {
      stop("Edges must join internal nodes to their children.");
    }
  }
  std::vector<uint8_t> reached(maxNode + 1, 0);
  reached[nTip + 1] = 1;
  for (int e = 0; e < nEdge; ++e) {
    if (!reached[parent[e]] || reached[child[e]]) {
      stop("Edges must be in preorder; a stale `order` attribute can hide "
           "this from `Preorder()`.");
    }
    reached[child[e]] = 1;
  }
}

// Pr(all observed tips in one state | k) for each mask (columns of `masks`,
// 1 = missing) at each k in kVals, under JC(k), or under the
// qHeterogeneity mixture when betaScale > 0. For tests of the kernels.
// [[Rcpp::export]]
NumericMatrix lumped_const_probs(IntegerVector parent, IntegerVector child,
                                 NumericVector edgeLen, int nTip,
                                 NumericVector rates, IntegerMatrix masks,
                                 NumericVector kVals, double betaScale = 0.0,
                                 int nBetaCat = 4) {
  check_preorder(parent, child, nTip);
  const int nMask = masks.ncol();
  std::vector<std::vector<uint8_t>> store(nMask, std::vector<uint8_t>(nTip));
  std::vector<const uint8_t*> ptr(nMask);
  for (int m = 0; m < nMask; ++m) {
    for (int t = 0; t < nTip; ++t) store[m][t] = masks(t, m) != 0;
    ptr[m] = store[m].data();
  }
  const std::vector<double> kv(kVals.begin(), kVals.end());
  std::vector<double> out(nMask * kv.size());
  if (betaScale > 0.0) {
    lumped_const_probs_het(parent, child, edgeLen, nTip, rates, betaScale,
                           nBetaCat, ptr, kv, out.data());
  } else {
    lumped_const_probs_jc(parent, child, edgeLen, nTip, rates, ptr, kv,
                          out.data());
  }
  NumericMatrix res(nMask, kv.size());
  for (int m = 0; m < nMask; ++m) {
    for (int i = 0; i < (int)kv.size(); ++i) res(m, i) = out[m * kv.size() + i];
  }
  // Return:
  return res;
}

// The joint ascertainment term of a transformational partition (1-based
// partIdx) at edge lengths already carrying the partition scale. For tests:
// `direct` sums every k exactly, and the result carries the table's node
// count as attribute "nNode".
// [[Rcpp::export]]
NumericVector joint_asc_log_cpp(SEXP dataPtr, int partIdx,
                                IntegerVector parent, IntegerVector child,
                                NumericVector scaledEdge, NumericVector rates,
                                double betaScale, double p, double alpha,
                                double beta, bool direct = false) {
  const McmcData* data = XPtr<McmcData>(dataPtr).get();
  if (partIdx < 1 || partIdx > (int)data->parts.size() ||
      data->parts[partIdx - 1].type != 1) {
    stop("partIdx must name a transformational partition.");
  }
  check_preorder(parent, child, data->nTip);
  KprimeHyper hyper;
  hyper.p = p;
  hyper.alpha = alpha;
  hyper.beta = beta;
  const PartInfo& part = data->parts[partIdx - 1];
  int lo;
  const std::vector<double> logPrior = kprime_log_prior_table(
    *data, data->kObs[part.globalCharIdx[0]], hyper, &lo);
  const JointAscTable table = joint_asc_table(
    *data, partIdx - 1, parent, child, scaledEdge, rates, betaScale, direct);
  NumericVector res(1, joint_asc_log(*data, partIdx - 1, table, logPrior,
                                     lo));
  res.attr("nNode") = table.nNode;
  // Return:
  return res;
}

// log pi(k) over the support the joint normaliser sums, for a character
// observed in kObs states; attribute "lo" gives its first k. For tests.
// [[Rcpp::export]]
NumericVector kprime_log_prior_table_cpp(SEXP dataPtr, int kObs, double p,
                                         double alpha, double beta) {
  const McmcData* data = XPtr<McmcData>(dataPtr).get();
  KprimeHyper hyper;
  hyper.p = p;
  hyper.alpha = alpha;
  hyper.beta = beta;
  int lo;
  const std::vector<double> lp = kprime_log_prior_table(*data, kObs, hyper,
                                                        &lo);
  NumericVector res(lp.begin(), lp.end());
  res.attr("lo") = lo;
  // Return:
  return res;
}
