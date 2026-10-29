// Test hooks for the public headers in inst/include/MkPrime. The umbrella
// header is included before Rcpp.h so that the build checks that it is
// self-contained and needs no R headers.
#include "../inst/include/MkPrime/MkPrime.h"
#include <Rcpp.h>
#include <algorithm>
#include <string>
#include <vector>

// [[Rcpp::export]]
int header_api_version() {
  return MKPRIME_HEADER_API_VERSION;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix header_jc_probs(int k, Rcpp::NumericVector t) {
  const MkPrime::JCRates rates(k);
  Rcpp::NumericMatrix out(t.size(), 4);
  for (R_xlen_t i = 0; i < t.size(); ++i) {
    const MkPrime::JCProbs p = rates.probs(t[i]);
    out(i, 0) = p.pSame;
    out(i, 1) = p.pDiff;
    out(i, 2) = p.diffCoeff;
    out(i, 3) = p.expTerm;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix header_mkn_probs(double rate_loss, Rcpp::NumericVector t) {
  const MkPrime::MkNRates rates(rate_loss);
  Rcpp::NumericMatrix out(t.size(), 4);
  for (R_xlen_t i = 0; i < t.size(); ++i) {
    const MkPrime::MkNProbs p = rates.probs(t[i]);
    out(i, 0) = p.p00;
    out(i, 1) = p.p01;
    out(i, 2) = p.p10;
    out(i, 3) = p.p11;
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::NumericVector header_lognormal_rates(double sd, Rcpp::NumericVector z) {
  Rcpp::NumericVector out(z.size());
  MkPrime::lognormal_rates(sd, z.size(), z.begin(), out.begin());
  return out;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix header_constant_patterns(int nTip, int nPattern,
                                             Rcpp::LogicalVector missing) {
  std::vector<unsigned char> miss(missing.begin(), missing.end());
  const std::vector<int> states = MkPrime::constant_patterns(
    nTip, nPattern, missing.size() ? miss.data() : nullptr);
  Rcpp::IntegerMatrix out(nTip, nPattern);
  std::copy(states.begin(), states.end(), out.begin());
  return out;
}

namespace {

// Applies `steps` to one Pruner, each inside checkpoint() ... accept() or
// reject(), and records the cached log-likelihood after each step against
// a from-scratch evaluation, plus the value after a rejection.
template <class Model, class MakeModel>
Rcpp::List run_trace(const Model& model, MakeModel makeModel,
                     const Rcpp::IntegerMatrix& tipStates,
                     const std::vector<double>& rates0,
                     const Rcpp::List& tree0, const Rcpp::List& steps) {
  const int nTip = tipStates.nrow(), nChar = tipStates.ncol();
  MkPrime::Pruner<Model> pruner(model, nTip, nChar, tipStates.begin(), rates0);
  Rcpp::IntegerVector p0 = tree0["parent"], c0 = tree0["child"];
  Rcpp::NumericVector l0 = tree0["edge_length"];
  pruner.set_tree(p0.begin(), c0.begin(), l0.begin(), p0.size());

  Rcpp::IntegerVector curP = Rcpp::clone(p0), curC = Rcpp::clone(c0);
  Rcpp::NumericVector curL = Rcpp::clone(l0);
  std::vector<double> curRates = rates0;
  Model curModel = model;

  const R_xlen_t n = steps.size();
  Rcpp::NumericVector cached(n), fresh(n), afterReject(n), before(n);
  Rcpp::IntegerVector recomputed(n);
  Rcpp::NumericVector initial = Rcpp::NumericVector::create(pruner.log_likelihood());
  for (R_xlen_t i = 0; i < n; ++i) {
    Rcpp::List step = steps[i];
    before[i] = pruner.log_likelihood();
    pruner.checkpoint();
    Rcpp::IntegerVector newP = curP, newC = curC;
    Rcpp::NumericVector newL = Rcpp::clone(curL);
    std::vector<double> newRates = curRates;
    Model newModel = curModel;
    if (step.containsElementNamed("parent")) {
      newP = step["parent"];
      newC = step["child"];
      newL = Rcpp::clone(Rcpp::NumericVector(step["edge_length"]));
      pruner.set_tree(newP.begin(), newC.begin(), newL.begin(), newP.size());
    }
    if (step.containsElementNamed("edge")) {
      Rcpp::IntegerVector edge = step["edge"];
      Rcpp::NumericVector len = step["length"];
      for (R_xlen_t j = 0; j < edge.size(); ++j) {
        pruner.set_edge_length(edge[j] - 1, len[j]);
        newL[edge[j] - 1] = len[j];
      }
    }
    if (step.containsElementNamed("rates")) {
      Rcpp::NumericVector r = step["rates"];
      newRates.assign(r.begin(), r.end());
      pruner.set_rates(newRates);
    }
    if (step.containsElementNamed("param")) {
      newModel = makeModel(Rcpp::as<double>(step["param"]));
      pruner.set_model(newModel);
    }
    pruner.update();
    recomputed[i] = pruner.n_recomputed();
    cached[i] = pruner.log_likelihood();
    fresh[i] = MkPrime::log_likelihood(newModel, newP.begin(), newC.begin(),
                                       newL.begin(), newP.size(), nTip, nChar,
                                       tipStates.begin(), newRates);
    if (Rcpp::as<bool>(step["accept"])) {
      pruner.accept();
      curP = newP;
      curC = newC;
      curL = newL;
      curRates = newRates;
      curModel = newModel;
      afterReject[i] = NA_REAL;
    } else {
      pruner.reject();
      afterReject[i] = pruner.log_likelihood();
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("initial") = initial, Rcpp::Named("before") = before,
    Rcpp::Named("cached") = cached, Rcpp::Named("fresh") = fresh,
    Rcpp::Named("after_reject") = afterReject,
    Rcpp::Named("recomputed") = recomputed);
}

}  // namespace

// Drives MkPrime::Pruner through a sequence of proposals. `model` is one of
// "jc" (param: unused; k = kStates), "jc_lumped" (kStates = kFull, kObs),
// "mkn" (param = rate_loss) or "f81" (pi, param = mu). Each element of
// `steps` is a list with `accept` and any of: parent / child / edge_length
// (new tree); edge / length (1-based edge indices); rates; param.
// [[Rcpp::export]]
Rcpp::List header_pruner_trace(std::string model, Rcpp::IntegerMatrix tip_states,
                               int kStates, int kObs, double param,
                               Rcpp::NumericVector pi,
                               Rcpp::NumericVector rates, Rcpp::List tree,
                               Rcpp::List steps) {
  const std::vector<double> r(rates.begin(), rates.end());
  if (model == "jc") {
    return run_trace(MkPrime::JCModel(kStates),
                     [kStates](double) { return MkPrime::JCModel(kStates); },
                     tip_states, r, tree, steps);
  }
  if (model == "jc_lumped") {
    return run_trace(MkPrime::JCLumpedModel(kStates, kObs),
                     [kStates, kObs](double) {
                       return MkPrime::JCLumpedModel(kStates, kObs);
                     },
                     tip_states, r, tree, steps);
  }
  if (model == "mkn") {
    return run_trace(MkPrime::MkNModel(param),
                     [](double x) { return MkPrime::MkNModel(x); },
                     tip_states, r, tree, steps);
  }
  if (model == "f81") {
    const std::vector<double> p(pi.begin(), pi.end());
    return run_trace(MkPrime::F81Model(p.data(), p.size(), param),
                     [p](double x) {
                       return MkPrime::F81Model(p.data(), p.size(), x);
                     },
                     tip_states, r, tree, steps);
  }
  Rcpp::stop("Unknown model `%s`.", model);
}

// Per-character likelihoods from a one-off Pruner (for checks against
// independent R implementations).
// [[Rcpp::export]]
Rcpp::NumericVector header_site_likelihoods(std::string model,
                                            Rcpp::IntegerVector parent,
                                            Rcpp::IntegerVector child,
                                            Rcpp::NumericVector edge_length,
                                            Rcpp::IntegerMatrix tip_states,
                                            int kStates, double param,
                                            Rcpp::NumericVector pi,
                                            Rcpp::NumericVector rates) {
  const std::vector<double> r(rates.begin(), rates.end());
  Rcpp::NumericVector out(tip_states.ncol());
  const int nTip = tip_states.nrow(), nChar = tip_states.ncol();
  if (model == "jc") {
    MkPrime::Pruner<MkPrime::JCModel> pr(MkPrime::JCModel(kStates), nTip, nChar,
                                         tip_states.begin(), r);
    pr.set_tree(parent.begin(), child.begin(), edge_length.begin(), parent.size());
    pr.site_likelihoods(out.begin());
  } else if (model == "mkn") {
    MkPrime::Pruner<MkPrime::MkNModel> pr(MkPrime::MkNModel(param), nTip, nChar,
                                          tip_states.begin(), r);
    pr.set_tree(parent.begin(), child.begin(), edge_length.begin(), parent.size());
    pr.site_likelihoods(out.begin());
  } else if (model == "f81") {
    MkPrime::Pruner<MkPrime::F81Model> pr(
      MkPrime::F81Model(pi.begin(), pi.size(), param), nTip, nChar,
      tip_states.begin(), r);
    pr.set_tree(parent.begin(), child.begin(), edge_length.begin(), parent.size());
    pr.site_likelihoods(out.begin());
  } else {
    Rcpp::stop("Unknown model `%s`.", model);
  }
  return out;
}

namespace {

// Test-only edge-aware model (header API v2): `Base` with each edge's length
// multiplied by its own factor.
template <class Base>
class EdgeScaled {
 public:
  EdgeScaled(const Base& base, const std::vector<double>& scale)
    : base_(base), scale_(scale) {}
  int n_states() const { return base_.n_states(); }
  void transition(const double* src, double* dst, int nChar, double t) const {
    base_.transition(src, dst, nChar, t);
  }
  void transition_edge(int edge, const double* src, double* dst, int nChar,
                       double t) const {
    base_.transition(src, dst, nChar, t * scale_.at(edge));
  }
  void transition_transpose_edge(int edge, const double* src, double* dst,
                                 int nChar, double t) const {
    base_.transition_transpose(src, dst, nChar, t * scale_.at(edge));
  }
  double root_site_lik(const double* cl) const {
    return base_.root_site_lik(cl);
  }

 private:
  Base base_;
  std::vector<double> scale_;
};

template <class Model>
Rcpp::List marginals_of(const Model& model, const Rcpp::IntegerVector& parent,
                        const Rcpp::IntegerVector& child,
                        const Rcpp::NumericVector& edgeLength,
                        const Rcpp::IntegerMatrix& tipStates,
                        const std::vector<double>& rates) {
  const int nTip = tipStates.nrow(), nChar = tipStates.ncol();
  MkPrime::Pruner<Model> pr(model, nTip, nChar, tipStates.begin(), rates);
  pr.set_tree(parent.begin(), child.begin(), edgeLength.begin(), parent.size());
  Rcpp::NumericVector site(nChar);
  pr.site_likelihoods(site.begin());
  pr.update_marginals();
  const int nNode = pr.n_node(), k = pr.n_states(), nCat = pr.n_cat();
  Rcpp::NumericVector marg(static_cast<R_xlen_t>(k) * nChar * nNode);
  Rcpp::NumericMatrix nodeLik(nChar, nNode);
  for (int node = 1; node <= nNode; ++node) {
    pr.node_marginals(node, marg.begin() +
                      static_cast<R_xlen_t>(node - 1) * k * nChar);
    for (int cat = 0; cat < nCat; ++cat) {
      const double* down = pr.node_cl(cat, node);
      const double* up = pr.node_up(cat, node);
      for (int c = 0; c < nChar; ++c) {
        double s = 0.0;
        for (int x = 0; x < k; ++x) s += down[c * k + x] * up[c * k + x];
        nodeLik(c, node - 1) += s / nCat;
      }
    }
  }
  marg.attr("dim") = Rcpp::IntegerVector::create(k, nChar, nNode);
  return Rcpp::List::create(Rcpp::Named("marginals") = marg,
                            Rcpp::Named("node_lik") = nodeLik,
                            Rcpp::Named("site_lik") = site);
}

}  // namespace

// P(t) and the kernel applied by transition_transpose(), as k x k matrices
// built column by column from unit vectors; `model` as for
// header_site_likelihoods().
// [[Rcpp::export]]
Rcpp::List header_transition_matrices(std::string model, int kStates,
                                      double param, Rcpp::NumericVector pi,
                                      double t) {
  Rcpp::NumericMatrix P(kStates, kStates), Pt(kStates, kStates);
  std::vector<double> unit(kStates), out(kStates);
  for (int j = 0; j < kStates; ++j) {
    std::fill(unit.begin(), unit.end(), 0.0);
    unit[j] = 1.0;
    for (int pass = 0; pass < 2; ++pass) {
      if (model == "jc") {
        const MkPrime::JCModel m(kStates);
        if (pass) m.transition_transpose(unit.data(), out.data(), 1, t);
        else m.transition(unit.data(), out.data(), 1, t);
      } else if (model == "mkn") {
        const MkPrime::MkNModel m(param);
        if (pass) m.transition_transpose(unit.data(), out.data(), 1, t);
        else m.transition(unit.data(), out.data(), 1, t);
      } else if (model == "f81") {
        const MkPrime::F81Model m(pi.begin(), pi.size(), param);
        if (pass) m.transition_transpose(unit.data(), out.data(), 1, t);
        else m.transition(unit.data(), out.data(), 1, t);
      } else {
        Rcpp::stop("Unknown model `%s`.", model);
      }
      for (int i = 0; i < kStates; ++i) (pass ? Pt : P)(i, j) = out[i];
    }
  }
  return Rcpp::List::create(Rcpp::Named("P") = P, Rcpp::Named("Pt") = Pt);
}

// Node marginals from Pruner::update_marginals(): `marginals` is a
// k x nChar x nNode array, `node_lik` the nChar x nNode products
// down(v) . up(v), averaged over categories, `site_lik` the site likelihoods.
// A non-empty `edge_scale` wraps the model in an edge-aware kernel that
// multiplies edge e's length by edge_scale[e].
// [[Rcpp::export]]
Rcpp::List header_marginals(std::string model, Rcpp::IntegerVector parent,
                            Rcpp::IntegerVector child,
                            Rcpp::NumericVector edge_length,
                            Rcpp::IntegerMatrix tip_states, int kStates,
                            double param, Rcpp::NumericVector pi,
                            Rcpp::NumericVector rates,
                            Rcpp::NumericVector edge_scale) {
  const std::vector<double> r(rates.begin(), rates.end());
  const std::vector<double> sc(edge_scale.begin(), edge_scale.end());
  const bool edge = !sc.empty();
  if (model == "jc") {
    const MkPrime::JCModel m(kStates);
    return edge ? marginals_of(EdgeScaled<MkPrime::JCModel>(m, sc), parent,
                               child, edge_length, tip_states, r)
                : marginals_of(m, parent, child, edge_length, tip_states, r);
  }
  if (model == "mkn") {
    const MkPrime::MkNModel m(param);
    return edge ? marginals_of(EdgeScaled<MkPrime::MkNModel>(m, sc), parent,
                               child, edge_length, tip_states, r)
                : marginals_of(m, parent, child, edge_length, tip_states, r);
  }
  if (model == "f81") {
    const MkPrime::F81Model m(pi.begin(), pi.size(), param);
    return edge ? marginals_of(EdgeScaled<MkPrime::F81Model>(m, sc), parent,
                               child, edge_length, tip_states, r)
                : marginals_of(m, parent, child, edge_length, tip_states, r);
  }
  Rcpp::stop("Unknown model `%s`.", model);
}

// An edge-aware Pruner given the same tree with its edges renumbered by
// `perm` (1-based): its cached likelihood against a fresh evaluation.
// [[Rcpp::export]]
Rcpp::NumericVector header_edge_renumber(Rcpp::IntegerVector parent,
                                         Rcpp::IntegerVector child,
                                         Rcpp::NumericVector edge_length,
                                         Rcpp::IntegerMatrix tip_states,
                                         int kStates,
                                         Rcpp::NumericVector edge_scale,
                                         Rcpp::IntegerVector perm) {
  typedef EdgeScaled<MkPrime::JCModel> Model;
  const Model m(MkPrime::JCModel(kStates),
                std::vector<double>(edge_scale.begin(), edge_scale.end()));
  const int nTip = tip_states.nrow(), nChar = tip_states.ncol();
  const int nEdge = parent.size();
  MkPrime::Pruner<Model> pr(m, nTip, nChar, tip_states.begin());
  pr.set_tree(parent.begin(), child.begin(), edge_length.begin(), nEdge);
  const double before = pr.log_likelihood();
  Rcpp::IntegerVector p2(nEdge), c2(nEdge);
  Rcpp::NumericVector l2(nEdge);
  for (int e = 0; e < nEdge; ++e) {
    p2[e] = parent[perm[e] - 1];
    c2[e] = child[perm[e] - 1];
    l2[e] = edge_length[perm[e] - 1];
  }
  pr.set_tree(p2.begin(), c2.begin(), l2.begin(), nEdge);
  const double cached = pr.log_likelihood();
  const double fresh = MkPrime::log_likelihood(m, p2.begin(), c2.begin(),
                                               l2.begin(), nEdge, nTip, nChar,
                                               tip_states.begin());
  return Rcpp::NumericVector::create(before, cached, fresh);
}
