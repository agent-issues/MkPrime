#ifndef MKPRIME_PRUNING_H
#define MKPRIME_PRUNING_H

// Felsenstein pruning with a persistent conditional-likelihood (CL) cache.
//
// Pruner<Model> holds CLs for every node and recomputes only nodes marked
// dirty: a node is dirty when its children, the lengths of the edges to its
// children, or any descendant's CLs have changed. set_edge_length() and
// set_tree() mark the affected nodes; set_model() and set_rates() mark all.
// An external MCMC loop brackets a proposal with checkpoint() and then
// accept() or reject(); reject() restores the CLs, tree, model and rates in
// place without recomputation.
//
// Trees are edge lists with 1-based node numbers (the ape / TreeTools
// convention): tips 1..nTip, internal nodes nTip+1..nNode, edges in any
// order. The root is the one node that is never a child; it may have two
// children (rooted) or three (unrooted). Internal nodes may have any number
// of children. Children of a node are combined in decreasing edge-index
// order, which reproduces MkPrime's reverse-preorder traversal exactly.
//
// Among-character rate variation (ACRV) uses nCat rate categories: the
// branch length for category i is t * rates[i], and site likelihoods are
// averaged over categories with equal weight (see acrv.h).
//
// Partitions: use one Pruner per partition. Each holds its own CLs and dirty
// flags, so an edge-length change must be applied to every Pruner.
//
// A character's CLs at a node are rescaled by a power of two once they fall
// below kRescaleBelow (scaling.h, header API v3), so a character's
// likelihood may fall far below 1e-308. node_cl() and node_up() return the
// stored, rescaled CLs; node_shift() and node_up_shift() give the shifts.

#include "models.h"
#include "scaling.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace MkPrime {

namespace detail {

// Header API v2: a model may define transition_edge(edge, src, dst, nChar, t)
// (and transition_transpose_edge) to receive the index of the edge, as given
// to set_tree(); Pruner then calls it in place of transition().
template <class M>
class has_transition_edge {
  template <class T>
  static auto test(int) -> decltype(std::declval<const T&>().transition_edge(
    0, static_cast<const double*>(nullptr), static_cast<double*>(nullptr), 0,
    0.0), std::true_type());
  template <class>
  static std::false_type test(...);

 public:
  static const bool value = decltype(test<M>(0))::value;
};

template <class M>
class has_transition_transpose_edge {
  template <class T>
  static auto test(int) -> decltype(
    std::declval<const T&>().transition_transpose_edge(
      0, static_cast<const double*>(nullptr), static_cast<double*>(nullptr), 0,
      0.0), std::true_type());
  template <class>
  static std::false_type test(...);

 public:
  static const bool value = decltype(test<M>(0))::value;
};

template <class M>
inline void transition(const M& m, int edge, const double* src, double* dst,
                       int nChar, double t, std::true_type) {
  m.transition_edge(edge, src, dst, nChar, t);
}
template <class M>
inline void transition(const M& m, int, const double* src, double* dst,
                       int nChar, double t, std::false_type) {
  m.transition(src, dst, nChar, t);
}
template <class M>
inline void transition_transpose(const M& m, int edge, const double* src,
                                 double* dst, int nChar, double t,
                                 std::true_type) {
  m.transition_transpose_edge(edge, src, dst, nChar, t);
}
template <class M>
inline void transition_transpose(const M& m, int, const double* src,
                                 double* dst, int nChar, double t,
                                 std::false_type) {
  m.transition_transpose(src, dst, nChar, t);
}

}  // namespace detail

template <class Model>
class Pruner {
 public:
  // True when Model defines transition_edge() (header API v2).
  static const bool kEdgeAware = detail::has_transition_edge<Model>::value;

  // tipStates: nTip x nChar, column-major (an R integer matrix); entries
  // 0..n_states()-1 or negative for missing (all states equally possible).
  // rates: one multiplier per ACRV category; {1.0} for none.
  Pruner(const Model& model, int nTip, int nChar, const int* tipStates,
         const std::vector<double>& rates = std::vector<double>(1, 1.0))
    : model_(model), nTip_(nTip), nChar_(nChar), k_(model.n_states()),
      stride_(static_cast<std::size_t>(nChar) * model.n_states()),
      rates_(rates), tipStates_(tipStates,
                                tipStates + static_cast<std::size_t>(nTip) * nChar),
      noShift_(nChar, 0) {
    if (nTip < 2) throw std::invalid_argument("Pruner: need >= 2 tips");
    if (nChar < 1) throw std::invalid_argument("Pruner: need >= 1 character");
    if (rates_.empty()) throw std::invalid_argument("Pruner: `rates` is empty");
    for (std::size_t i = 0; i < tipStates_.size(); ++i) {
      if (tipStates_[i] >= k_) {
        throw std::invalid_argument("Pruner: tip state exceeds n_states() - 1");
      }
    }
  }

  // Set or replace the tree. Edge lengths are in expected changes per
  // character. Only nodes whose children or child edge lengths differ from
  // the previous tree (and their ancestors) are marked dirty.
  void set_tree(const int* parent, const int* child, const double* edgeLength,
                int nEdge) {
    int nNode = 0;
    for (int e = 0; e < nEdge; ++e) {
      if (parent[e] > nNode) nNode = parent[e];
      if (child[e] > nNode) nNode = child[e];
      if (parent[e] < 1 || child[e] < 1) {
        throw std::invalid_argument("Pruner: node numbers must be 1-based");
      }
    }
    if (nNode <= nTip_) throw std::invalid_argument("Pruner: no internal nodes");

    std::vector<int> parentOf(nNode + 1, 0);
    std::vector<std::vector<int> > kids(nNode + 1);
    std::vector<double> len(nNode + 1, 0.0);
    std::vector<int> edgeOf(nNode + 1, -1);
    for (int e = nEdge - 1; e >= 0; --e) {
      const int p = parent[e], c = child[e];
      if (parentOf[c]) throw std::invalid_argument("Pruner: node has two parents");
      if (p <= nTip_) throw std::invalid_argument("Pruner: a tip has children");
      parentOf[c] = p;
      kids[p].push_back(c);
      len[c] = edgeLength[e];
      edgeOf[c] = e;
    }
    int root = 0;
    for (int n = 1; n <= nNode; ++n) {
      if (n <= nTip_ && !parentOf[n]) {
        throw std::invalid_argument("Pruner: a tip is not in the tree");
      }
      if (n > nTip_ && !parentOf[n]) {
        if (root) throw std::invalid_argument("Pruner: more than one root");
        root = n;
      }
      if (n > nTip_ && kids[n].empty()) {
        throw std::invalid_argument("Pruner: internal node has no children");
      }
    }
    if (!root) throw std::invalid_argument("Pruner: tree has no root");

    const bool fresh = nNode != nNode_;
    if (fresh) {
      save_all_cls();
      nNode_ = nNode;
      allocate();
    }
    for (int n = nTip_ + 1; n <= nNode; ++n) {
      bool changed = fresh || kids[n] != kids_[n];
      for (std::size_t i = 0; !changed && i < kids[n].size(); ++i) {
        changed = len[kids[n][i]] != len_[kids[n][i]] ||
          (kEdgeAware && edgeOf[kids[n][i]] != edgeOf_[kids[n][i]]);
      }
      if (changed) dirty_[n] = 1;
    }
    parentOf_.swap(parentOf);
    kids_.swap(kids);
    len_.swap(len);
    edgeOf_.swap(edgeOf);
    edgeChild_.assign(child, child + nEdge);
    root_ = root;
    build_postorder();
    // A node left dirty under the old tree may have new, clean ancestors.
    for (std::size_t i = 0; i < post_.size(); ++i) {
      const int n = post_[i];
      for (std::size_t j = 0; j < kids_[n].size(); ++j) {
        if (dirty_[kids_[n][j]]) dirty_[n] = 1;
      }
    }
  }

  // Change one edge length; `edge` indexes the arrays given to set_tree().
  void set_edge_length(int edge, double length) {
    require_tree();
    if (edge < 0 || edge >= static_cast<int>(edgeChild_.size())) {
      throw std::out_of_range("Pruner: edge index out of range");
    }
    const int c = edgeChild_[edge];
    if (len_[c] == length) return;
    len_[c] = length;
    mark_dirty(parentOf_[c]);
  }

  double edge_length(int edge) const { return len_[edgeChild_.at(edge)]; }

  // Mark the CLs that depend on one edge's transition kernel stale, e.g.
  // after an edge-aware model's weights for that edge change.
  void invalidate_edge(int edge) {
    require_tree();
    if (edge < 0 || edge >= static_cast<int>(edgeChild_.size())) {
      throw std::out_of_range("Pruner: edge index out of range");
    }
    mark_dirty(parentOf_[edgeChild_[edge]]);
  }

  void set_model(const Model& model) {
    if (model.n_states() != k_) {
      throw std::invalid_argument("Pruner: model changes n_states()");
    }
    model_ = model;
    invalidate();
  }

  const Model& model() const { return model_; }

  void set_rates(const std::vector<double>& rates) {
    if (rates.empty()) throw std::invalid_argument("Pruner: `rates` is empty");
    if (rates.size() != rates_.size() && nNode_) {
      save_all_cls();
      // Tip CLs may have been set by set_tip_partials(); carry them over.
      const std::vector<double> tips(node_cl(0, 1), node_cl(0, nTip_ + 1));
      rates_ = rates;
      allocate();
      for (int cat = 0; cat < n_cat(); ++cat) {
        double* dst = node_cl_mut(cat, 1);
        for (std::size_t j = 0; j < tips.size(); ++j) dst[j] = tips[j];
      }
    } else {
      rates_ = rates;
    }
    invalidate();
  }

  const std::vector<double>& rates() const { return rates_; }

  // Replace a tip's CLs (nChar * n_states() values), e.g. for ambiguous or
  // uncertain codings. `tip` is 1-based.
  void set_tip_partials(int tip, const double* cl) {
    require_tree();
    if (tip < 1 || tip > nTip_) throw std::out_of_range("Pruner: bad tip");
    journal(tip);
    for (int cat = 0; cat < n_cat(); ++cat) {
      double* dst = node_cl_mut(cat, tip);
      for (std::size_t j = 0; j < stride_; ++j) dst[j] = cl[j];
    }
    mark_dirty(parentOf_[tip]);
  }

  // Mark every internal node dirty.
  void invalidate() {
    for (int n = nTip_ + 1; n <= nNode_; ++n) dirty_[n] = 1;
  }

  // Bring all CLs up to date; returns the number of nodes recomputed.
  int update() {
    require_tree();
    int nDone = 0;
    for (std::size_t i = 0; i < post_.size(); ++i) {
      const int n = post_[i];
      if (!dirty_[n]) continue;
      journal(n);
      recompute(n);
      dirty_[n] = 0;
      ++nDone;
    }
    nRecomputed_ = nDone;
    return nDone;
  }

  // Per-character likelihoods, averaged over rate categories; out[nChar].
  // A likelihood below about 1e-308 underflows: see log_site_likelihoods().
  void site_likelihoods(double* out) {
    scaled_site_likelihoods();
    for (int c = 0; c < nChar_; ++c) {
      out[c] = siteShift_[c] ? std::ldexp(siteLik_[c], -siteShift_[c])
                             : siteLik_[c];
    }
  }

  // Per-character log likelihoods, averaged over rate categories; out[nChar].
  void log_site_likelihoods(double* out) {
    scaled_site_likelihoods();
    for (int c = 0; c < nChar_; ++c) {
      out[c] = siteLik_[c] > 0.0
        ? log_scaled(siteLik_[c], siteShift_[c])
        : -std::numeric_limits<double>::infinity();
    }
  }

  // Sum of log site likelihoods; -Inf if any character has likelihood <= 0.
  double log_likelihood() {
    scaled_site_likelihoods();
    double logLik = 0.0;
    for (int c = 0; c < nChar_; ++c) {
      if (siteLik_[c] <= 0.0) return -std::numeric_limits<double>::infinity();
      logLik += log_scaled(siteLik_[c], siteShift_[c]);
    }
    return logLik;
  }

  // Sum of site likelihoods: with constant_patterns() (ascertainment.h) as
  // tip data, the probability that a character is constant.
  double sum_site_likelihoods() {
    scaled_site_likelihoods();
    double total = 0.0;
    for (int c = 0; c < nChar_; ++c) {
      total += siteShift_[c] ? std::ldexp(siteLik_[c], -siteShift_[c])
                             : siteLik_[c];
    }
    return total;
  }

  // Upward (pre-order) pass, header API v2. For every node v, rate category
  // and character, up(v)[s] is the likelihood of the data outside the clade
  // below v, jointly with state s at v; at the root it holds the root
  // frequencies, read from root_site_lik(). Then down(v) . up(v) is the
  // site likelihood for every v, and node_marginals() normalises the
  // products. Needs Model::transition_transpose(src, dst, nChar, t), dst =
  // P(t)' src (or transition_transpose_edge() for edge-aware models).
  // Reflects the tree, model and rates at the time of the call; call again
  // after any change, including reject().
  void update_marginals() {
    update();
    up_.assign(cl_.size(), 0.0);
    upShift_.assign(shift_.size(), 0);
    std::vector<double> rootFreq(k_), unit(k_, 0.0);
    for (int s = 0; s < k_; ++s) {
      unit[s] = 1.0;
      rootFreq[s] = model_.root_site_lik(unit.data());
      unit[s] = 0.0;
    }
    std::vector<double> below, outside(stride_);
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double rate = rates_[cat];
      double* upRoot = node_up_mut(cat, root_);
      for (int c = 0; c < nChar_; ++c) {
        for (int s = 0; s < k_; ++s) {
          upRoot[static_cast<std::size_t>(c) * k_ + s] = rootFreq[s];
        }
      }
      // post_ holds internal nodes, children first: walk it backwards.
      for (std::size_t i = post_.size(); i-- > 0;) {
        const int n = post_[i];
        const std::vector<int>& nk = kids_[n];
        below.resize(nk.size() * stride_);
        for (std::size_t j = 0; j < nk.size(); ++j) {
          apply_transition(nk[j], node_cl(cat, nk[j]), below.data() + j * stride_,
                           len_[nk[j]] * rate);
        }
        const double* upN = up_.data() +
          (static_cast<std::size_t>(cat) * (nNode_ + 1) + n) * stride_;
        for (std::size_t j = 0; j < nk.size(); ++j) {
          for (std::size_t x = 0; x < stride_; ++x) outside[x] = upN[x];
          for (std::size_t o = 0; o < nk.size(); ++o) {
            if (o == j) continue;
            const double* b = below.data() + o * stride_;
            for (std::size_t x = 0; x < stride_; ++x) outside[x] *= b[x];
          }
          double* upJ = node_up_mut(cat, nk[j]);
          apply_transpose(nk[j], outside.data(), upJ, len_[nk[j]] * rate);
          int* shJ = shift_slot(upShift_, cat, nk[j]);
          const int* shN = shift_slot(upShift_, cat, n);
          for (int c = 0; c < nChar_; ++c) shJ[c] = shN[c];
          for (std::size_t o = 0; o < nk.size(); ++o) {
            if (o == j) continue;
            const int* shO = node_shift(cat, nk[o]);
            for (int c = 0; c < nChar_; ++c) shJ[c] += shO[c];
          }
          rescale_node(upJ, shJ);
        }
      }
    }
  }

  // Upward-pass vectors of `node` in category `cat`, laid out as node_cl();
  // character c's block is scaled by 2^node_up_shift(cat, node)[c].
  const double* node_up(int cat, int node) const {
    require_marginals();
    return up_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
  }

  const int* node_up_shift(int cat, int node) const {
    require_marginals();
    return upShift_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * nChar_;
  }

  // Posterior probability of each state at `node` (1-based, tips included),
  // given all the data: out[c * n_states() + s], averaged over rate
  // categories with weights proportional to their site likelihoods.
  void node_marginals(int node, double* out) const {
    if (node < 1 || node > nNode_) throw std::out_of_range("Pruner: bad node");
    for (std::size_t x = 0; x < stride_; ++x) out[x] = 0.0;
    // A category's weight is 2^-shift times its stored product; scale each
    // character relative to its largest-weight (least-shifted) category.
    std::vector<int> minShift(nChar_, std::numeric_limits<int>::max());
    for (int cat = 0; cat < n_cat(); ++cat) {
      const int* sd = node_shift(cat, node);
      const int* su = node_up_shift(cat, node);
      for (int c = 0; c < nChar_; ++c) {
        if (sd[c] + su[c] < minShift[c]) minShift[c] = sd[c] + su[c];
      }
    }
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double* down = node_cl(cat, node);
      const double* up = node_up(cat, node);
      const int* sd = node_shift(cat, node);
      const int* su = node_up_shift(cat, node);
      for (int c = 0; c < nChar_; ++c) {
        const int rel = minShift[c] - sd[c] - su[c];
        const std::size_t off = static_cast<std::size_t>(c) * k_;
        for (int s = 0; s < k_; ++s) {
          out[off + s] += rel ? std::ldexp(down[off + s] * up[off + s], rel)
                              : down[off + s] * up[off + s];
        }
      }
    }
    for (int c = 0; c < nChar_; ++c) {
      double* o = out + static_cast<std::size_t>(c) * k_;
      double total = 0.0;
      for (int s = 0; s < k_; ++s) total += o[s];
      for (int s = 0; s < k_; ++s) {
        o[s] = total > 0.0 ? o[s] / total
                           : std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  // Start a proposal. Until accept() or reject(), every CL that is
  // overwritten is saved first.
  void checkpoint() {
    saved_.model.assign(1, model_);
    saved_.rates = rates_;
    saved_.nNode = nNode_;
    saved_.root = root_;
    saved_.parentOf = parentOf_;
    saved_.kids = kids_;
    saved_.len = len_;
    saved_.edgeOf = edgeOf_;
    saved_.edgeChild = edgeChild_;
    saved_.post = post_;
    saved_.dirty = dirty_;
    saved_.allCl.clear();
    saved_.allShift.clear();
    saved_.allScaled.clear();
    saved_.hasAllCl = false;
    journalNodes_.clear();
    journalCl_.clear();
    journalShift_.clear();
    journalScaled_.clear();
    journaled_.assign(nNode_ + 1, 0);
    inProposal_ = true;
  }

  void accept() {
    inProposal_ = false;
    journalNodes_.clear();
    journalCl_.clear();
    journalShift_.clear();
    journalScaled_.clear();
    saved_.allCl.clear();
    saved_.allShift.clear();
    saved_.allScaled.clear();
  }

  void reject() {
    if (!inProposal_) throw std::logic_error("Pruner: reject() without checkpoint()");
    model_ = saved_.model[0];
    rates_.swap(saved_.rates);
    nNode_ = saved_.nNode;
    root_ = saved_.root;
    parentOf_.swap(saved_.parentOf);
    kids_.swap(saved_.kids);
    len_.swap(saved_.len);
    edgeOf_.swap(saved_.edgeOf);
    edgeChild_.swap(saved_.edgeChild);
    post_.swap(saved_.post);
    dirty_.swap(saved_.dirty);
    if (saved_.hasAllCl) {
      cl_.swap(saved_.allCl);
      shift_.swap(saved_.allShift);
      scaled_.swap(saved_.allScaled);
    } else {
      const std::size_t block = n_cat() * stride_;
      const std::size_t shiftBlock = static_cast<std::size_t>(n_cat()) * nChar_;
      for (std::size_t i = 0; i < journalNodes_.size(); ++i) {
        const int n = journalNodes_[i];
        for (int cat = 0; cat < n_cat(); ++cat) {
          double* dst = node_cl_mut(cat, n);
          const double* src = journalCl_.data() + i * block + cat * stride_;
          for (std::size_t j = 0; j < stride_; ++j) dst[j] = src[j];
          int* shDst = shift_slot(shift_, cat, n);
          const int* shSrc = journalShift_.data() + i * shiftBlock + cat * nChar_;
          for (int c = 0; c < nChar_; ++c) shDst[c] = shSrc[c];
          scaled_[static_cast<std::size_t>(cat) * (nNode_ + 1) + n] =
            journalScaled_[i * n_cat() + cat];
        }
      }
    }
    accept();
  }

  // CLs of `node` in rate category `cat`: nChar blocks of n_states() values,
  // block c scaled by 2^node_shift(cat, node)[c]. Up to date only after
  // update() or a likelihood call.
  const double* node_cl(int cat, int node) const {
    return cl_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
  }

  const int* node_shift(int cat, int node) const {
    const std::size_t i = static_cast<std::size_t>(cat) * (nNode_ + 1) + node;
    return scaled_[i] ? shift_.data() + i * nChar_ : noShift_.data();
  }

  int n_tip() const { return nTip_; }
  int n_node() const { return nNode_; }
  int n_char() const { return nChar_; }
  int n_states() const { return k_; }
  int n_cat() const { return static_cast<int>(rates_.size()); }
  int root() const { return root_; }
  int parent_of(int node) const { return parentOf_.at(node); }
  const std::vector<int>& children(int node) const { return kids_.at(node); }
  bool is_dirty(int node) const { return dirty_.at(node) != 0; }
  int n_recomputed() const { return nRecomputed_; }

 private:
  struct Saved {
    std::vector<Model> model;  // vector: Model need not be default-constructible
    std::vector<double> rates;
    int nNode = 0;
    int root = 0;
    std::vector<int> parentOf;
    std::vector<std::vector<int> > kids;
    std::vector<double> len;
    std::vector<int> edgeOf;
    std::vector<int> edgeChild;
    std::vector<int> post;
    std::vector<std::uint8_t> dirty;
    std::vector<double> allCl;
    std::vector<int> allShift;
    std::vector<std::uint8_t> allScaled;
    bool hasAllCl = false;
  };

  double* node_cl_mut(int cat, int node) {
    return cl_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
  }

  int* shift_slot(std::vector<int>& shifts, int cat, int node) {
    return shifts.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * nChar_;
  }

  void require_tree() const {
    if (!nNode_) throw std::logic_error("Pruner: call set_tree() first");
  }

  void require_marginals() const {
    if (up_.size() != cl_.size()) {
      throw std::logic_error("Pruner: call update_marginals() first");  // # nocov
    }
  }

  // Rescales each character block of `cl` that has fallen below
  // kRescaleBelow, adding the shift to shift[c].
  void rescale_node(double* cl, int* shift) const {
    for (int c = 0; c < nChar_; ++c) {
      double* block = cl + static_cast<std::size_t>(c) * k_;
      double sum = 0.0;
      for (int s = 0; s < k_; ++s) sum += block[s];
      if (sum < kRescaleBelow) shift[c] += rescale_cl(block, k_);
    }
  }

  // dst *= src, then rescale as rescale_node(), setting `any` and zeroing
  // `shift` first if it was not already in use.
  void multiply_rescale(double* dst, const double* src, int* shift,
                        bool& any) const {
    for (int c = 0; c < nChar_; ++c) {
      double* block = dst + static_cast<std::size_t>(c) * k_;
      const double* f = src + static_cast<std::size_t>(c) * k_;
      double sum = 0.0;
      for (int s = 0; s < k_; ++s) {
        block[s] *= f[s];
        sum += block[s];
      }
      if (sum < kRescaleBelow) {
        if (!any) {
          for (int j = 0; j < nChar_; ++j) shift[j] = 0;
          any = true;
        }
        shift[c] += rescale_cl(block, k_);
      }
    }
  }

  // Root site likelihoods averaged over categories, into siteLik_ scaled by
  // 2^siteShift_.
  void scaled_site_likelihoods() {
    update();
    const double invNCat = 1.0 / n_cat();
    siteLik_.assign(nChar_, 0.0);
    siteShift_.assign(nChar_, 0);
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double* rootCl = node_cl(cat, root_);
      const int* rootShift = node_shift(cat, root_);
      for (int c = 0; c < nChar_; ++c) {
        add_scaled(siteLik_[c], siteShift_[c],
                   model_.root_site_lik(rootCl + static_cast<std::size_t>(c) * k_),
                   rootShift[c]);
      }
    }
    for (int c = 0; c < nChar_; ++c) siteLik_[c] *= invNCat;
  }

  void allocate() {
    cl_.assign(static_cast<std::size_t>(n_cat()) * (nNode_ + 1) * stride_, 0.0);
    shift_.assign(static_cast<std::size_t>(n_cat()) * (nNode_ + 1) * nChar_, 0);
    scaled_.assign(static_cast<std::size_t>(n_cat()) * (nNode_ + 1), 0);
    dirty_.assign(nNode_ + 1, 1);
    for (int t = 0; t <= nTip_; ++t) dirty_[t] = 0;
    if (static_cast<int>(kids_.size()) != nNode_ + 1) kids_.assign(nNode_ + 1, std::vector<int>());
    if (static_cast<int>(len_.size()) != nNode_ + 1) len_.assign(nNode_ + 1, 0.0);
    journaled_.assign(nNode_ + 1, 1);  // a reallocation is undone wholesale
    init_tips();
  }

  void init_tips() {
    for (int tip = 1; tip <= nTip_; ++tip) {
      double* cl = node_cl_mut(0, tip);
      for (std::size_t j = 0; j < stride_; ++j) cl[j] = 0.0;
      for (int c = 0; c < nChar_; ++c) {
        const int state = tipStates_[static_cast<std::size_t>(c) * nTip_ + (tip - 1)];
        const std::size_t off = static_cast<std::size_t>(c) * k_;
        if (state < 0) {
          for (int s = 0; s < k_; ++s) cl[off + s] = 1.0;
        } else {
          cl[off + state] = 1.0;
        }
      }
      for (int cat = 1; cat < n_cat(); ++cat) {
        double* dst = node_cl_mut(cat, tip);
        for (std::size_t j = 0; j < stride_; ++j) dst[j] = cl[j];
      }
    }
  }

  // Before a reallocation inside a proposal, keep the whole CL array.
  void save_all_cls() {
    if (inProposal_ && !saved_.hasAllCl) {
      saved_.allCl = cl_;
      saved_.allShift = shift_;
      saved_.allScaled = scaled_;
      saved_.hasAllCl = true;
    }
  }

  void journal(int node) {
    if (!inProposal_ || saved_.hasAllCl || journaled_[node]) return;
    journaled_[node] = 1;
    journalNodes_.push_back(node);
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double* src = node_cl(cat, node);
      journalCl_.insert(journalCl_.end(), src, src + stride_);
      const int* sh = node_shift(cat, node);
      journalShift_.insert(journalShift_.end(), sh, sh + nChar_);
      journalScaled_.push_back(
        scaled_[static_cast<std::size_t>(cat) * (nNode_ + 1) + node]);
    }
  }

  void mark_dirty(int node) {
    while (node && !dirty_[node]) {
      dirty_[node] = 1;
      node = parentOf_[node];
    }
  }

  void build_postorder() {
    post_.clear();
    std::vector<int> stack(1, root_);
    while (!stack.empty()) {
      const int n = stack.back();
      stack.pop_back();
      if (n <= nTip_) continue;
      post_.push_back(n);
      for (std::size_t i = 0; i < kids_[n].size(); ++i) stack.push_back(kids_[n][i]);
    }
    for (std::size_t i = 0, j = post_.size(); i + 1 < j; ++i, --j) {
      const int tmp = post_[i];
      post_[i] = post_[j - 1];
      post_[j - 1] = tmp;
    }
  }

  void recompute(int n) {
    scratch_.resize(stride_);
    const std::vector<int>& nk = kids_[n];
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double rate = rates_[cat];
      double* dst = node_cl_mut(cat, n);
      const std::size_t slot = static_cast<std::size_t>(cat) * (nNode_ + 1) + n;
      int* shift = shift_.data() + slot * nChar_;
      bool any = false;
      for (std::size_t i = 0; i < nk.size(); ++i) {
        if (!scaled_[static_cast<std::size_t>(cat) * (nNode_ + 1) + nk[i]]) {
          continue;
        }
        if (!any) for (int c = 0; c < nChar_; ++c) shift[c] = 0;
        any = true;
        const int* shI = node_shift(cat, nk[i]);
        for (int c = 0; c < nChar_; ++c) shift[c] += shI[c];
      }
      apply_transition(nk[0], node_cl(cat, nk[0]), dst, len_[nk[0]] * rate);
      for (std::size_t i = 1; i < nk.size(); ++i) {
        apply_transition(nk[i], node_cl(cat, nk[i]), scratch_.data(),
                         len_[nk[i]] * rate);
        // Check after the last factor, and after every factor from the third
        // on, so that a polytomy cannot underflow before it is checked.
        if (i + 1 < nk.size() && i < 2) {
          for (std::size_t j = 0; j < stride_; ++j) dst[j] *= scratch_[j];
        } else {
          multiply_rescale(dst, scratch_.data(), shift, any);
        }
      }
      if (nk.size() == 1) {
        if (!any) std::fill(shift, shift + nChar_, 0);
        rescale_node(dst, shift);
        any = any || std::count(shift, shift + nChar_, 0) < nChar_;
      }
      scaled_[slot] = any;
    }
  }

  // P(t) src for the edge above `childNode`.
  void apply_transition(int childNode, const double* src, double* dst,
                        double t) const {
    detail::transition(model_, edgeOf_[childNode], src, dst, nChar_, t,
                       std::integral_constant<bool, kEdgeAware>());
  }

  void apply_transpose(int childNode, const double* src, double* dst,
                       double t) const {
    detail::transition_transpose(
      model_, edgeOf_[childNode], src, dst, nChar_, t,
      std::integral_constant<bool,
        detail::has_transition_transpose_edge<Model>::value>());
  }

  double* node_up_mut(int cat, int node) {
    return up_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
  }

  Model model_;
  int nTip_;
  int nChar_;
  int k_;
  std::size_t stride_;
  std::vector<double> rates_;
  std::vector<int> tipStates_;

  int nNode_ = 0;
  int root_ = 0;
  std::vector<int> parentOf_;
  std::vector<std::vector<int> > kids_;
  std::vector<double> len_;
  std::vector<int> edgeOf_;
  std::vector<int> edgeChild_;
  std::vector<int> post_;
  std::vector<std::uint8_t> dirty_;
  std::vector<double> cl_;
  // Shifts per (cat, node, character), laid out as cl_; read only where
  // scaled_[(cat, node)] is set, elsewhere every shift is 0.
  std::vector<int> shift_;
  std::vector<std::uint8_t> scaled_;
  std::vector<int> noShift_;
  std::vector<double> scratch_;
  std::vector<double> siteLik_;
  std::vector<int> siteShift_;
  std::vector<double> up_;  // upward-pass vectors (update_marginals())
  std::vector<int> upShift_;
  int nRecomputed_ = 0;

  bool inProposal_ = false;
  Saved saved_;
  std::vector<int> journalNodes_;
  std::vector<double> journalCl_;
  std::vector<int> journalShift_;
  std::vector<std::uint8_t> journalScaled_;
  std::vector<std::uint8_t> journaled_;
};

// One-shot log-likelihood, without caching.
template <class Model>
double log_likelihood(const Model& model, const int* parent, const int* child,
                      const double* edgeLength, int nEdge, int nTip, int nChar,
                      const int* tipStates,
                      const std::vector<double>& rates = std::vector<double>(1, 1.0)) {
  Pruner<Model> pruner(model, nTip, nChar, tipStates, rates);
  pruner.set_tree(parent, child, edgeLength, nEdge);
  return pruner.log_likelihood();
}

}  // namespace MkPrime

#endif  // MKPRIME_PRUNING_H
