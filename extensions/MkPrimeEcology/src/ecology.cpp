// Ecology-Biased Equilibrium (EBE) model as an MkPrime engine extension
// (extension ABI 3, header API 2). Specification: dev/ecology/ebe-spec.md on
// the archived branch (tag ecology-aware-2026-05-29), sections 1-3 and 7.
//
// The extension claims MkPrime's neomorphic characters. For character c on
// an edge, the MkN equilibrium is tilted in each ecology s by
//   logit(pi1) = logit(pi1_base) + s_z log(phi_s),
// s_z = +1 (z = 1), -1 (z = 2), 0 (z = 0 or s = reference ecology), with
// pi1_base = 1 / (1 + rate_loss); the decay rate stays MkPrime's base MkN rate
// lambda(rate_loss) = (1 + r)^2 / (2 r), as the spec's "tilt the attractor, not
// the tempo" (its fixed 2 is lambda at r = 1, before MkPrime rescaled MkN).
// The edge's kernel mixes these over ecology with weights wEdge(e, s): the
// mean of the JC(K) marginals of tip ecology at the edge's two ends, on the
// tree's own edge lengths with a uniform root. The root takes the base
// equilibrium.
// Neomorphic edge lengths are MkPrime's: edge length x neo_scale x ACRV rate.
#include <Rcpp.h>
#include <MkPrime/MkPrime.h>
#include <MkPrime/extension.h>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

const double kNegInf = -std::numeric_limits<double>::infinity();

struct EcoData {
  int nTip = 0, kEco = 0, refEco = 0, nNonRef = 0, nNeo = 0;
  int perEcology = 0;  // 0: one phi; 1: one phi per non-reference ecology
  int nPhi = 1;
  std::vector<int> ecoTip;     // nTip; 0..kEco-1, -1 missing
  std::vector<int> nonRefIdx;  // kEco; index among non-reference, -1 for ref
  double pi0A = 360, pi0B = 120, thetaA = 2, thetaB = 2, sigmaPhi = 1.5;
  std::vector<std::string> paramStr, intStr, moveStr;
  std::vector<const char*> paramNames, intNames, moveNames;
  std::vector<double> paramInit, moveWeight;
  std::vector<int> intInit, intLevels;
  MkPrimeExtension ext;
};

// Parameter layout: phi (nPhi), pi0, theta (nNonRef).
inline double pi0_of(const EcoData& d, const double* p) { return p[d.nPhi]; }
inline double theta_of(const EcoData& d, const double* p, int j) {
  return p[d.nPhi + 1 + j];
}

double log_dbeta(double x, double a, double b) {
  return (a - 1.0) * std::log(x) + (b - 1.0) * std::log1p(-x) -
    (std::lgamma(a) + std::lgamma(b) - std::lgamma(a + b));
}

// log P(z | pi0, theta_j): none with pi0; encouraged / discouraged split the
// rest as theta_j / (1 - theta_j).
inline double log_pz(int z, double pi0, double theta) {
  return z == 0 ? std::log(pi0) :
    std::log1p(-pi0) + (z == 1 ? std::log(theta) : std::log1p(-theta));
}

// Tilted equilibrium P(present) for influence z in ecology state s.
inline double tilted_pi1(const EcoData& d, int z, int s, double pi1Base,
                         const double* phi) {
  const int j = d.nonRefIdx[s];
  if (z == 0 || j < 0) return pi1Base;
  const double phiS = d.perEcology ? phi[j] : phi[0];
  const double logit = std::log(pi1Base / (1.0 - pi1Base)) +
    (z == 1 ? 1.0 : -1.0) * std::log(phiS);
  if (logit >= 0.0) return 1.0 / (1.0 + std::exp(-logit));
  const double e = std::exp(logit);
  return e / (1.0 + e);
}

// Edge-aware EBE kernel (header API v2). Column c uses influence row
// colZ[c * nNonRef ...]. With weights summing to one, the mixture
// sum_s w(e, s) P_s equals ex I + (1 - ex) 1 pibar', where
// pibar = sum_s w(e, s) (pi0_s, pi1_s) and ex = exp(-lambda t).
class EbeModel {
 public:
  EbeModel(const EcoData& d, const double* phi, double rateLoss,
           const std::vector<double>* wEdge, const std::vector<int>* colZ)
    : d_(&d), wEdge_(wEdge), colZ_(colZ), pi1Base_(1.0 / (1.0 + rateLoss)),
      lambda_(MkPrime::MkNRates(rateLoss).lambda()),
      pi1_(3 * d.kEco) {
    for (int z = 0; z < 3; ++z) {
      for (int s = 0; s < d.kEco; ++s) {
        pi1_[z * d.kEco + s] = tilted_pi1(d, z, s, pi1Base_, phi);
      }
    }
  }

  int n_states() const { return 2; }

  void transition(const double*, double*, int, double) const {
    throw std::logic_error("EbeModel: edge-aware kernel called without edge");
  }

  void transition_edge(int edge, const double* src, double* dst, int nChar,
                       double t) const {
    const double oneMinus = -std::expm1(-lambda_ * t);
    const double ex = 1.0 - oneMinus;
    const int K = d_->kEco, nNR = d_->nNonRef;
    const double* w = wEdge_->data() + static_cast<std::size_t>(edge) * K;
    for (int c = 0; c < nChar; ++c) {
      const int* zRow = colZ_->data() + static_cast<std::size_t>(c) * nNR;
      double pibar1 = 0.0;
      for (int s = 0; s < K; ++s) {
        const int j = d_->nonRefIdx[s];
        const int z = j < 0 ? 0 : zRow[j];
        pibar1 += w[s] * pi1_[z * K + s];
      }
      const double s0 = src[2 * c], s1 = src[2 * c + 1];
      const double base = oneMinus * ((1.0 - pibar1) * s0 + pibar1 * s1);
      dst[2 * c] = ex * s0 + base;
      dst[2 * c + 1] = ex * s1 + base;
    }
  }

  double root_site_lik(const double* cl) const {
    return (1.0 - pi1Base_) * cl[0] + pi1Base_ * cl[1];
  }

 private:
  const EcoData* d_;
  const std::vector<double>* wEdge_;
  const std::vector<int>* colZ_;
  double pi1Base_;
  double lambda_;
  std::vector<double> pi1_;
};

// wEdge (nEdge x kEco, row-major) from JC(K) marginals of tip ecology.
bool edge_weights(const EcoData& d, const MkPrimeTreeView* v,
                  std::vector<double>& wEdge) {
  MkPrime::Pruner<MkPrime::JCModel> pr(MkPrime::JCModel(d.kEco), d.nTip, 1,
                                       d.ecoTip.data());
  pr.set_tree(v->parent, v->child, v->edge_length, v->n_edge);
  pr.update_marginals();
  const int nNode = pr.n_node();
  std::vector<double> marg(static_cast<std::size_t>(nNode + 1) * d.kEco);
  for (int n = 1; n <= nNode; ++n) {
    pr.node_marginals(n, marg.data() + static_cast<std::size_t>(n) * d.kEco);
  }
  wEdge.assign(static_cast<std::size_t>(v->n_edge) * d.kEco, 0.0);
  for (int e = 0; e < v->n_edge; ++e) {
    const double* mp = marg.data() + static_cast<std::size_t>(v->parent[e]) * d.kEco;
    const double* mc = marg.data() + static_cast<std::size_t>(v->child[e]) * d.kEco;
    for (int s = 0; s < d.kEco; ++s) {
      const double w = 0.5 * (mp[s] + mc[s]);
      if (!std::isfinite(w)) return false;
      wEdge[static_cast<std::size_t>(e) * d.kEco + s] = w;
    }
  }
  return true;
}

// Log-likelihood of claimed characters `chars` (indices into the claimed
// block); entry i is scored with influence row i of `zRows`. Under variable
// coding each is conditioned on varying among its observed tips. perChar,
// when given, receives the per-entry terms.
double ebe_log_lik(const EcoData& d, const MkPrimeTreeView* v,
                   const double* params, const std::vector<int>& chars,
                   const std::vector<int>& zRows,
                   const std::vector<double>& wEdge,
                   std::vector<double>* perChar = nullptr) {
  const int nC = static_cast<int>(chars.size());
  if (perChar) perChar->assign(nC, 0.0);
  if (nC == 0) return 0.0;
  if (v->coding == 2) return kNegInf;  // informative coding: not supported
  const bool asc = v->coding == 1;
  const int nRep = asc ? 3 : 1;
  const int nCol = nRep * nC;
  std::vector<int> states(static_cast<std::size_t>(d.nTip) * nCol);
  std::vector<int> colZ(static_cast<std::size_t>(nCol) * d.nNonRef);
  for (int i = 0; i < nC; ++i) {
    const int* src = v->claimed_states + static_cast<std::size_t>(chars[i]) * d.nTip;
    for (int t = 0; t < d.nTip; ++t) {
      states[static_cast<std::size_t>(i) * d.nTip + t] = src[t];
      if (asc) {
        // Constant over the observed tips; missing tips stay missing.
        states[static_cast<std::size_t>(nC + i) * d.nTip + t] = src[t] < 0 ? -1 : 0;
        states[static_cast<std::size_t>(2 * nC + i) * d.nTip + t] = src[t] < 0 ? -1 : 1;
      }
    }
    for (int rep = 0; rep < nRep; ++rep) {
      for (int j = 0; j < d.nNonRef; ++j) {
        colZ[static_cast<std::size_t>(rep * nC + i) * d.nNonRef + j] =
          zRows[static_cast<std::size_t>(i) * d.nNonRef + j];
      }
    }
  }
  const double rateLoss = v->core[MKPRIME_CORE_RATE_LOSS];
  if (!(rateLoss > 0.0)) return kNegInf;
  EbeModel model(d, params, rateLoss, &wEdge, &colZ);
  const std::vector<double> rates(v->rate_cat, v->rate_cat + v->n_rate_cat);
  MkPrime::Pruner<EbeModel> pr(model, d.nTip, nCol, states.data(), rates);
  std::vector<double> len(v->edge_length, v->edge_length + v->n_edge);
  for (double& l : len) l *= v->neo_scale;
  pr.set_tree(v->parent, v->child, len.data(), v->n_edge);
  std::vector<double> site(nCol);
  pr.site_likelihoods(site.data());
  double ll = 0.0;
  for (int i = 0; i < nC; ++i) {
    double li = std::log(site[i]);
    if (asc) li -= std::log1p(-(site[nC + i] + site[2 * nC + i]));
    if (!std::isfinite(li)) return kNegInf;
    if (perChar) (*perChar)[i] = li;
    ll += li;
  }
  return ll;
}

bool params_ok(const EcoData& d, const double* p) {
  for (int j = 0; j < d.nPhi; ++j) if (!(p[j] > 0.0)) return false;
  const double pi0 = pi0_of(d, p);
  if (!(pi0 > 0.0 && pi0 < 1.0)) return false;
  for (int j = 0; j < d.nNonRef; ++j) {
    const double th = theta_of(d, p, j);
    if (!(th > 0.0 && th < 1.0)) return false;
  }
  return true;
}

bool view_ok(const EcoData& d, const MkPrimeTreeView* v) {
  return v->n_tip == d.nTip && v->n_claimed == d.nNeo && v->n_core >= 5;
}

// phi ~ LogNormal(0, sigmaPhi); pi0 ~ Beta; theta_j ~ Beta; z | pi0, theta.
double eco_log_prior(const double* p, const int* z, void* user) {
  const EcoData& d = *static_cast<const EcoData*>(user);
  if (!params_ok(d, p)) return kNegInf;
  const double sd = d.sigmaPhi;
  double lp = 0.0;
  for (int j = 0; j < d.nPhi; ++j) {
    const double lx = std::log(p[j]);
    lp += -lx - std::log(sd) - 0.5 * std::log(2.0 * M_PI) -
      0.5 * lx * lx / (sd * sd);
  }
  const double pi0 = pi0_of(d, p);
  lp += log_dbeta(pi0, d.pi0A, d.pi0B);
  for (int j = 0; j < d.nNonRef; ++j) {
    lp += log_dbeta(theta_of(d, p, j), d.thetaA, d.thetaB);
  }
  for (int c = 0; c < d.nNeo; ++c) {
    for (int j = 0; j < d.nNonRef; ++j) {
      lp += log_pz(z[c * d.nNonRef + j], pi0, theta_of(d, p, j));
    }
  }
  return lp;
}

double eco_log_lik(const MkPrimeTreeView* v, const double* p, const int* z,
                   void* user) {
  try {
    const EcoData& d = *static_cast<const EcoData*>(user);
    if (!view_ok(d, v) || !params_ok(d, p)) return kNegInf;
    std::vector<double> wEdge;
    if (!edge_weights(d, v, wEdge)) return kNegInf;
    std::vector<int> chars(d.nNeo);
    for (int c = 0; c < d.nNeo; ++c) chars[c] = c;
    const std::vector<int> zRows(z, z + static_cast<std::size_t>(d.nNeo) * d.nNonRef);
    return ebe_log_lik(d, v, p, chars, zRows, wEdge);
  } catch (...) {
    return kNegInf;
  }
}

// Move 0, "gibbs_z": draw one character's influence row from its full
// conditional. All 3^(K-1) rows are scored in one pruning pass. The log
// Hastings ratio log p(old | .) - log p(new | .) cancels the cold target
// ratio, so the cold chain always accepts.
double eco_propose(int move, const MkPrimeTreeView* v, double* p, int* z,
                   void* user) {
  try {
    const EcoData& d = *static_cast<const EcoData*>(user);
    if (move != 0 || !view_ok(d, v) || !params_ok(d, p) || d.nNeo == 0) {
      return kNegInf;
    }
    int c = static_cast<int>(v->unif_rand() * d.nNeo);
    if (c >= d.nNeo) c = d.nNeo - 1;
    int nConfig = 1;
    for (int j = 0; j < d.nNonRef; ++j) nConfig *= 3;
    std::vector<int> chars(nConfig, c);
    std::vector<int> zRows(static_cast<std::size_t>(nConfig) * d.nNonRef);
    int current = 0;
    for (int j = 0, mult = 1; j < d.nNonRef; ++j, mult *= 3) {
      current += z[c * d.nNonRef + j] * mult;
    }
    for (int cfg = 0; cfg < nConfig; ++cfg) {
      int rest = cfg;
      for (int j = 0; j < d.nNonRef; ++j) {
        zRows[static_cast<std::size_t>(cfg) * d.nNonRef + j] = rest % 3;
        rest /= 3;
      }
    }
    std::vector<double> wEdge, ll;
    if (!edge_weights(d, v, wEdge)) return kNegInf;
    if (!std::isfinite(ebe_log_lik(d, v, p, chars, zRows, wEdge, &ll))) {
      return kNegInf;
    }
    const double pi0 = pi0_of(d, p);
    std::vector<double> logw(nConfig);
    double mx = kNegInf;
    for (int cfg = 0; cfg < nConfig; ++cfg) {
      double lw = ll[cfg];
      for (int j = 0; j < d.nNonRef; ++j) {
        lw += log_pz(zRows[static_cast<std::size_t>(cfg) * d.nNonRef + j], pi0,
                     theta_of(d, p, j));
      }
      logw[cfg] = lw;
      if (lw > mx) mx = lw;
    }
    double total = 0.0;
    for (double lw : logw) total += std::exp(lw - mx);
    const double u = v->unif_rand() * total;
    double acc = 0.0;
    int pick = nConfig - 1;
    for (int cfg = 0; cfg < nConfig; ++cfg) {
      acc += std::exp(logw[cfg] - mx);
      if (u < acc) { pick = cfg; break; }
    }
    for (int j = 0; j < d.nNonRef; ++j) {
      z[c * d.nNonRef + j] = zRows[static_cast<std::size_t>(pick) * d.nNonRef + j];
    }
    return logw[current] - logw[pick];
  } catch (...) {
    return kNegInf;
  }
}

void eco_finalise(SEXP x) {
  delete static_cast<EcoData*>(R_ExternalPtrAddr(x));
  R_ClearExternalPtr(x);
}

// A tree view assembled from R values, for tests and diagnostics.
struct RView {
  MkPrimeTreeView v;
  std::vector<double> core, rates;
};

void r_view(RView& rv, Rcpp::IntegerVector& parent, Rcpp::IntegerVector& child,
            Rcpp::NumericVector& edgeLength, Rcpp::IntegerMatrix& states,
            double rateLoss, double neoScale, const Rcpp::NumericVector& rates,
            int coding) {
  rv.core.assign({1.0, rateLoss, 0.0, 1.0, 1.0});
  rv.rates.assign(rates.begin(), rates.end());
  MkPrimeTreeView& v = rv.v;
  v = MkPrimeTreeView();
  v.n_tip = states.nrow();
  v.n_edge = parent.size();
  v.parent = parent.begin();
  v.child = child.begin();
  v.edge_length = edgeLength.begin();
  v.unif_rand = &::unif_rand;
  v.n_core = 5;
  v.core = rv.core.data();
  v.neo_scale = neoScale;
  v.n_rate_cat = static_cast<int>(rv.rates.size());
  v.rate_cat = rv.rates.data();
  v.coding = coding;
  v.n_claimed = states.ncol();
  v.claimed_states = states.begin();
  v.claimed_char = nullptr;
}

const MkPrimeExtension* ext_of(SEXP ptr) {
  if (TYPEOF(ptr) != EXTPTRSXP ||
      R_ExternalPtrTag(ptr) != Rf_install("MkPrimeExtension")) {
    Rcpp::stop("Not an MkPrimeExtension pointer.");
  }
  const MkPrimeExtension* e =
    static_cast<const MkPrimeExtension*>(R_ExternalPtrAddr(ptr));
  if (!e || std::string(e->name) != "ecology") {
    Rcpp::stop("Not an ecology extension.");
  }
  return e;
}

}  // namespace

// Build the extension. `ecology`: tip ecology states, 0-based, NA missing,
// in MkPrime's tip order; `nNeo`: number of neomorphic characters MkPrime
// will hand over. `priors`: pi0 Beta shapes, theta Beta shapes, sigmaPhi.
// The returned pointer owns its data.
// [[Rcpp::export(.EcologyExtensionPtr)]]
SEXP ecology_extension_ptr(Rcpp::IntegerVector ecology, int kEco, int refEco,
                           int nNeo, bool perEcology,
                           Rcpp::CharacterVector paramNames,
                           Rcpp::NumericVector paramInit,
                           Rcpp::CharacterVector zNames,
                           Rcpp::NumericVector priors, double zWeight) {
  const int nNonRef = kEco - 1;
  const int nPhi = perEcology ? nNonRef : 1;
  const int nParam = nPhi + 1 + nNonRef;
  if (kEco < 2 || refEco < 0 || refEco >= kEco || nNeo < 0 ||
      paramNames.size() != nParam || paramInit.size() != nParam ||
      zNames.size() != static_cast<R_xlen_t>(nNeo) * nNonRef ||
      priors.size() != 5 || !(zWeight > 0.0)) {
    Rcpp::stop("Inconsistent ecology extension settings.");
  }
  EcoData* d = new EcoData();
  d->nTip = ecology.size();
  d->kEco = kEco;
  d->refEco = refEco;
  d->nNonRef = nNonRef;
  d->nNeo = nNeo;
  d->perEcology = perEcology ? 1 : 0;
  d->nPhi = nPhi;
  d->ecoTip.resize(d->nTip);
  for (int t = 0; t < d->nTip; ++t) {
    d->ecoTip[t] = Rcpp::IntegerVector::is_na(ecology[t]) ? -1 : ecology[t];
  }
  d->nonRefIdx.assign(kEco, -1);
  for (int s = 0, j = 0; s < kEco; ++s) {
    if (s != refEco) d->nonRefIdx[s] = j++;
  }
  d->pi0A = priors[0];
  d->pi0B = priors[1];
  d->thetaA = priors[2];
  d->thetaB = priors[3];
  d->sigmaPhi = priors[4];
  for (R_xlen_t i = 0; i < paramNames.size(); ++i) {
    d->paramStr.push_back(Rcpp::as<std::string>(paramNames[i]));
  }
  for (R_xlen_t i = 0; i < zNames.size(); ++i) {
    d->intStr.push_back(Rcpp::as<std::string>(zNames[i]));
  }
  d->moveStr.push_back("gibbs_z");
  for (const std::string& s : d->paramStr) d->paramNames.push_back(s.c_str());
  for (const std::string& s : d->intStr) d->intNames.push_back(s.c_str());
  for (const std::string& s : d->moveStr) d->moveNames.push_back(s.c_str());
  d->paramInit.assign(paramInit.begin(), paramInit.end());
  d->intInit.assign(zNames.size(), 0);
  d->intLevels.assign(zNames.size(), 3);
  d->moveWeight.assign(1, zWeight);

  MkPrimeExtension& e = d->ext;
  e = MkPrimeExtension();
  e.abi_version = MKPRIME_EXTENSION_ABI_VERSION;
  e.name = "ecology";
  e.n_params = nParam;
  e.param_names = d->paramNames.data();
  e.param_init = d->paramInit.data();
  e.user = d;
  e.n_moves = nNeo > 0 ? 1 : 0;
  e.move_names = d->moveNames.data();
  e.move_weight = d->moveWeight.data();
  e.claims = MKPRIME_CLAIM_NEOMORPHIC;
  e.n_ints = static_cast<int>(d->intStr.size());
  e.int_names = d->intNames.data();
  e.int_init = d->intInit.data();
  e.int_levels = d->intLevels.data();
  e.log_prior_ints = eco_log_prior;
  e.log_lik_ints = eco_log_lik;
  e.propose_ints = eco_propose;

  SEXP owner = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
  R_RegisterCFinalizerEx(owner, eco_finalise, TRUE);
  SEXP out = MkPrimeExtensionPtrProtect(&d->ext, owner);
  UNPROTECT(1);
  return out;
}

// The extension's log-likelihood on a tree given from R. `states`:
// nTip x nNeo, 0/1, -1 missing.
// [[Rcpp::export(.EcologyLogLik)]]
double ecology_log_lik(SEXP ptr, Rcpp::IntegerVector parent,
                       Rcpp::IntegerVector child, Rcpp::NumericVector edgeLength,
                       Rcpp::IntegerMatrix states, double rateLoss,
                       double neoScale, Rcpp::NumericVector rates, int coding,
                       Rcpp::NumericVector params, Rcpp::IntegerVector z) {
  const MkPrimeExtension* e = ext_of(ptr);
  RView rv;
  r_view(rv, parent, child, edgeLength, states, rateLoss, neoScale, rates,
         coding);
  return e->log_lik_ints(&rv.v, params.begin(), z.begin(), e->user);
}

// [[Rcpp::export(.EcologyLogPrior)]]
double ecology_log_prior(SEXP ptr, Rcpp::NumericVector params,
                         Rcpp::IntegerVector z) {
  const MkPrimeExtension* e = ext_of(ptr);
  return e->log_prior_ints(params.begin(), z.begin(), e->user);
}

// One gibbs_z proposal: the new z and its log Hastings ratio.
// [[Rcpp::export(.EcologyProposeZ)]]
Rcpp::List ecology_propose_z(SEXP ptr, Rcpp::IntegerVector parent,
                             Rcpp::IntegerVector child,
                             Rcpp::NumericVector edgeLength,
                             Rcpp::IntegerMatrix states, double rateLoss,
                             double neoScale, Rcpp::NumericVector rates,
                             int coding, Rcpp::NumericVector params,
                             Rcpp::IntegerVector z) {
  const MkPrimeExtension* e = ext_of(ptr);
  RView rv;
  r_view(rv, parent, child, edgeLength, states, rateLoss, neoScale, rates,
         coding);
  Rcpp::IntegerVector zNew = Rcpp::clone(z);
  Rcpp::NumericVector p = Rcpp::clone(params);
  const double logH = e->propose_ints(0, &rv.v, p.begin(), zNew.begin(),
                                      e->user);
  return Rcpp::List::create(Rcpp::Named("z") = zNew,
                            Rcpp::Named("log_hastings") = logH);
}
