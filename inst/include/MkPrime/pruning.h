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
// No rescaling is applied: CLs underflow only when a single character's
// likelihood falls below about 1e-308.

#include "models.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace MkPrime {

template <class Model>
class Pruner {
 public:
  // tipStates: nTip x nChar, column-major (an R integer matrix); entries
  // 0..n_states()-1 or negative for missing (all states equally possible).
  // rates: one multiplier per ACRV category; {1.0} for none.
  Pruner(const Model& model, int nTip, int nChar, const int* tipStates,
         const std::vector<double>& rates = std::vector<double>(1, 1.0))
    : model_(model), nTip_(nTip), nChar_(nChar), k_(model.n_states()),
      stride_(static_cast<std::size_t>(nChar) * model.n_states()),
      rates_(rates), tipStates_(tipStates,
                                tipStates + static_cast<std::size_t>(nTip) * nChar) {
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
        changed = len[kids[n][i]] != len_[kids[n][i]];
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
  void site_likelihoods(double* out) {
    update();
    const double invNCat = 1.0 / n_cat();
    for (int c = 0; c < nChar_; ++c) out[c] = 0.0;
    for (int cat = 0; cat < n_cat(); ++cat) {
      const double* rootCl = node_cl(cat, root_);
      for (int c = 0; c < nChar_; ++c) {
        out[c] += model_.root_site_lik(rootCl + static_cast<std::size_t>(c) * k_);
      }
    }
    for (int c = 0; c < nChar_; ++c) out[c] *= invNCat;
  }

  // Sum of log site likelihoods; -Inf if any character has likelihood <= 0.
  double log_likelihood() {
    siteLik_.resize(nChar_);
    site_likelihoods(siteLik_.data());
    double logLik = 0.0;
    for (int c = 0; c < nChar_; ++c) {
      if (siteLik_[c] <= 0.0) return -std::numeric_limits<double>::infinity();
      logLik += std::log(siteLik_[c]);
    }
    return logLik;
  }

  // Sum of site likelihoods: with constant_patterns() (ascertainment.h) as
  // tip data, the probability that a character is constant.
  double sum_site_likelihoods() {
    siteLik_.resize(nChar_);
    site_likelihoods(siteLik_.data());
    double total = 0.0;
    for (int c = 0; c < nChar_; ++c) total += siteLik_[c];
    return total;
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
    saved_.hasAllCl = false;
    journalNodes_.clear();
    journalCl_.clear();
    journaled_.assign(nNode_ + 1, 0);
    inProposal_ = true;
  }

  void accept() {
    inProposal_ = false;
    journalNodes_.clear();
    journalCl_.clear();
    saved_.allCl.clear();
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
    } else {
      const std::size_t block = n_cat() * stride_;
      for (std::size_t i = 0; i < journalNodes_.size(); ++i) {
        const int n = journalNodes_[i];
        for (int cat = 0; cat < n_cat(); ++cat) {
          double* dst = node_cl_mut(cat, n);
          const double* src = journalCl_.data() + i * block + cat * stride_;
          for (std::size_t j = 0; j < stride_; ++j) dst[j] = src[j];
        }
      }
    }
    accept();
  }

  // CLs of `node` in rate category `cat`: nChar blocks of n_states() values.
  // Up to date only after update() or a likelihood call.
  const double* node_cl(int cat, int node) const {
    return cl_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
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
    bool hasAllCl = false;
  };

  double* node_cl_mut(int cat, int node) {
    return cl_.data() + (static_cast<std::size_t>(cat) * (nNode_ + 1) + node) * stride_;
  }

  void require_tree() const {
    if (!nNode_) throw std::logic_error("Pruner: call set_tree() first");
  }

  void allocate() {
    cl_.assign(static_cast<std::size_t>(n_cat()) * (nNode_ + 1) * stride_, 0.0);
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
      model_.transition(node_cl(cat, nk[0]), dst, nChar_, len_[nk[0]] * rate);
      for (std::size_t i = 1; i < nk.size(); ++i) {
        model_.transition(node_cl(cat, nk[i]), scratch_.data(), nChar_,
                          len_[nk[i]] * rate);
        for (std::size_t j = 0; j < stride_; ++j) dst[j] *= scratch_[j];
      }
    }
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
  std::vector<double> scratch_;
  std::vector<double> siteLik_;
  int nRecomputed_ = 0;

  bool inProposal_ = false;
  Saved saved_;
  std::vector<int> journalNodes_;
  std::vector<double> journalCl_;
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
