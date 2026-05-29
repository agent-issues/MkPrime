# Independent R reference implementation ("oracle") of the Ecology-Biased
# Equilibrium (EBE) neomorphic likelihood.
#
# DERIVED SOLELY FROM dev/ecology/ebe-spec.md (the math contract). The author
# of this file did NOT read src/mcmc_ecology.cpp or any C++ kernel source, so
# the cross-check in test-ebe-likelihood.R is a genuine independent check and
# not a tautology. (The R orchestration in R/likelihood.R was read for the
# entry-point SIGNATURES only; it currently implements the *old* v2 rate-
# scaling semantics, which are deliberately NOT reproduced here.)
#
# ---------------------------------------------------------------------------
# Spec recap (the ONLY new math, ebe-spec.md s2):
#
#   Neomorphic binary char, state 0 = absent, 1 = present.
#   pi1_base = 1 / (1 + rate_loss),  pi0_base = rate_loss / (1 + rate_loss).
#   For a character c on an edge of ecology s:
#       s_z = +1 if z=1, -1 if z=2, 0 if z=0 or s == refEcology
#       logit(pi1) = logit(pi1_base) + s_z * log(phi_s)
#       pi1_{z,s}  = expit(...) ;  pi0_{z,s} = 1 - pi1_{z,s}
#   Decay rate is FIXED at lambda = 2 (tilt the attractor, not the tempo):
#       ex  = exp(-2 * tEff)        tEff = edge_length * rate_neo * ACRV_rate
#       P00 = pi0 + pi1*ex   P01 = pi1*(1 - ex)
#       P10 = pi0*(1 - ex)   P11 = pi1 + pi0*ex
#   Mixture over edge ecology: P_mix(e) = sum_s wEdge(e, s) * P_{z(c,s), s}.
#   Root distribution: (pi0_base, pi1_base) -- unchanged base equilibrium.
#   ACRV: average the per-character site likelihood over rate categories.
#   Variable coding: subtract log(1 - pConst) per character, where
#       pConst = P(all tips 0) + P(all tips 1)
#   computed by the SAME mixed-P pruning with all-0 / all-1 pseudo-tips, the
#   SAME (real-ecology) wEdge, and the SAME root distribution.
#
# phi_s indexing (spec s3): global mode (magnitudeMode == 0) uses phi[1];
#   per_ecology (mode == 1) uses phi[ecoState + 1] (1-based R indexing).
#
# theta / pi0 are PRIOR parameters only; gammaE is forced to 1 (spec s4), so
# NONE of theta / pi0 / gammaE enter this likelihood.
# ---------------------------------------------------------------------------


# expit / logit (avoid a hard dependency on stats::plogis being attached).
.ebe_logit <- function(p) log(p) - log1p(-p)
.ebe_expit <- function(x) 1 / (1 + exp(-x))


# Tilted equilibrium pi1 for state z and ecology s (spec s2).
# z: integer in {0, 1, 2}; ecoState: 0-based ecology index; refEcology: 0-based.
# phi: numeric vector; mode: 0 = global, 1 = per_ecology.
.ebe_pi1 <- function(pi1_base, z, ecoState, refEcology, phi, mode) {
  if (ecoState == refEcology || z == 0L) return(pi1_base)
  phi_s <- if (mode == 0L) phi[1] else phi[ecoState + 1L]
  s_z <- if (z == 1L) 1 else -1            # z == 2 -> -1
  .ebe_expit(.ebe_logit(pi1_base) + s_z * log(phi_s))
}


# EBE 2x2 transition matrix for one (z, s) on an edge with effective time
# tEff, returned in the orientation P[from + 1, to + 1] (rows = ancestor
# state, cols = descendant state) -- matching the pruning convention
# P %*% childPartial used below. lambda is hard-fixed at 2.
.ebe_P <- function(pi1_base, z, ecoState, refEcology, phi, mode, tEff) {
  pi1 <- .ebe_pi1(pi1_base, z, ecoState, refEcology, phi, mode)
  pi0 <- 1 - pi1
  ex  <- exp(-2 * tEff)
  matrix(c(pi0 + pi1 * ex, pi1 * (1 - ex),      # row 0 (from state 0): to0, to1
           pi0 * (1 - ex), pi1 + pi0 * ex),     # row 1 (from state 1): to0, to1
         nrow = 2L, byrow = TRUE)
}


# z-lookup mapping ecology state s (0-based) to the zMat column for character
# c. Reference ecology has no column and is treated as z = 0. Non-reference
# ecologies map to 1-based column index (s < refEcology) ? s + 1 : s, matching
# the existing kEco-1-column zMat layout in the test harness.
.ebe_z_for <- function(zMat, c, s, refEcology) {
  if (s == refEcology) return(0L)
  zCol <- if (s < refEcology) s + 1L else s
  as.integer(zMat[c, zCol])
}


# Felsenstein postorder pruning of ONE binary character with a per-edge,
# per-ecology mixed transition matrix. Returns the site likelihood (NOT log)
# for a single ACRV rate category.
#
# parent/child: edge endpoints, 1-based node ids (ape convention: tips
#   1..nTip, root = nTip + 1). edges may be in any postorder-compatible
#   ordering; we process them in REVERSE input order (children before
#   parents) exactly as the existing reference pruner does.
# tipPartial: nTip x 2 matrix of tip partial likelihoods (one-hot for an
#   observed state; c(1, 1) for missing). Used directly so the same pruner
#   serves the main pass and the all-0 / all-1 const-site pseudo-passes.
# tEffEdge: numeric vector, effective branch time per edge (edge_length *
#   rate_neo * ACRV_rate), already scaled.
# zVec: integer vector indexed by 1-based zMat column (length kEco - 1).
# wEdge: nEdge x kEco matrix of ecology mixture weights.
.ebe_prune_char <- function(parent, child, tEffEdge, tipPartial,
                            zVec, wEdge, pi1_base, rootFreqs,
                            refEcology, phi, mode) {
  nTip    <- nrow(tipPartial)
  maxNode <- 2L * nTip - 1L
  root    <- nTip + 1L
  kEco    <- ncol(wEdge)

  cl   <- matrix(0, nrow = maxNode, ncol = 2L)
  init <- rep(FALSE, maxNode)
  cl[seq_len(nTip), ] <- tipPartial
  init[seq_len(nTip)] <- TRUE

  # zVec is the per-character z row, indexed by 1-based column. Reconstruct a
  # single-row matrix so .ebe_z_for can index it uniformly.
  zMat1 <- matrix(as.integer(zVec), nrow = 1L)

  for (e in rev(seq_along(parent))) {
    pa <- parent[e]
    ch <- child[e]
    tEff <- tEffEdge[e]

    Pmix <- matrix(0, 2L, 2L)
    for (s in seq_len(kEco) - 1L) {           # s is 0-based ecology state
      w <- wEdge[e, s + 1L]
      if (w == 0) next
      z <- .ebe_z_for(zMat1, 1L, s, refEcology)
      Pmix <- Pmix + w * .ebe_P(pi1_base, z, s, refEcology, phi, mode, tEff)
    }

    msg <- as.vector(Pmix %*% cl[ch, ])
    if (!init[pa]) cl[pa, ] <- msg else cl[pa, ] <- cl[pa, ] * msg
    init[pa] <- TRUE
  }

  sum(rootFreqs * cl[root, ])
}


#' Independent EBE neomorphic log-likelihood oracle
#'
#' Pure-R reference for the EBE neomorphic likelihood, implementing
#' `dev/ecology/ebe-spec.md` section 2 exactly. Used to cross-check the C++
#' kernel (test contract T2). Handles multiple binary characters, ACRV, the
#' ecology mixture, and the variable-coding constant-site correction.
#'
#' @param parent,child Integer edge endpoints (1-based ape node ids).
#' @param tipStates nTip x nChar integer matrix; 0 / 1 observed, -1 missing.
#' @param rate_loss Loss:gain ratio. `pi1_base = 1 / (1 + rate_loss)`.
#' @param phi Numeric magnitude. Global mode: scalar (or length-1). Per-ecology
#'   mode: length kEco, indexed `phi[ecoState + 1]`.
#' @param zMat nChar x (kEco - 1) integer matrix, entries in {0, 1, 2}.
#' @param wEdge Optional nEdge x kEco mixture-weight matrix. If NULL, it is
#'   reconstructed from `ecologyTip` via `.EcologyNodeMarginals` +
#'   `.EcologyEdgeWeights` (the documented reconstruction helpers).
#' @param ecologyTip Optional integer tip-ecology vector (0-based states, -1
#'   missing); required when `wEdge` is NULL.
#' @param kEco Number of ecology states (required when reconstructing wEdge).
#' @param refEcology 0-based reference ecology index (no z column; z = 0).
#' @param rate_neo Neomorphic partition rate scalar.
#' @param edgeLen Numeric branch lengths (UNSCALED; `rate_neo` and ACRV are
#'   applied internally).
#' @param rate_log_sd,nCat ACRV lognormal sd and number of categories.
#' @param mode 0 = global phi, 1 = per_ecology phi.
#' @param coding 0 = none, 1 = variable (apply `-log(1 - pConst)` per char).
#'
#' @return Scalar log-likelihood (sum over characters).
#' @keywords internal
ebe_loglik_R <- function(parent, child, edgeLen, tipStates,
                         rate_loss, phi, zMat,
                         wEdge = NULL, ecologyTip = NULL, kEco = NULL,
                         refEcology = 0L,
                         rate_neo = 1.0,
                         rate_log_sd = 0, nCat = 1L,
                         mode = 0L, coding = 0L) {
  storage.mode(tipStates) <- "integer"
  nTip  <- nrow(tipStates)
  nChar <- ncol(tipStates)
  refEcology <- as.integer(refEcology)
  mode <- as.integer(mode)
  coding <- as.integer(coding)

  # Reconstruct wEdge from tip ecology if not supplied (documented helpers).
  if (is.null(wEdge)) {
    if (is.null(ecologyTip) || is.null(kEco)) {
      stop("ebe_loglik_R: supply either wEdge, or both ecologyTip and kEco")
    }
    eco <- as.integer(ecologyTip)
    eco[is.na(eco)] <- -1L
    marg  <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLen,
                                             eco, as.integer(kEco))
    wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  }
  kEco <- ncol(wEdge)

  pi1_base  <- 1 / (1 + rate_loss)
  pi0_base  <- rate_loss / (1 + rate_loss)
  rootFreqs <- c(pi0_base, pi1_base)

  # ACRV rate multipliers (same routine the package uses).
  rates <- DiscreteLognormalRates(rate_log_sd, nCat)
  nCat  <- length(rates)

  # Effective per-edge time for a given ACRV multiplier: edge * rate_neo * mult.
  neoEl <- edgeLen * rate_neo

  total <- 0
  for (c in seq_len(nChar)) {
    zVec <- if ((kEco - 1L) >= 1L) zMat[c, ] else integer(0)

    # --- main pass: average site likelihood over ACRV categories ---
    tipPartial <- matrix(0, nrow = nTip, ncol = 2L)
    for (i in seq_len(nTip)) {
      st <- tipStates[i, c]
      if (st < 0L) tipPartial[i, ] <- 1 else tipPartial[i, st + 1L] <- 1
    }

    siteLik <- 0
    for (cat_idx in seq_len(nCat)) {
      tEffEdge <- neoEl * rates[cat_idx]
      siteLik <- siteLik +
        .ebe_prune_char(parent, child, tEffEdge, tipPartial,
                        zVec, wEdge, pi1_base, rootFreqs,
                        refEcology, phi, mode)
    }
    avg <- siteLik / nCat
    if (avg <= 0) return(-Inf)
    charLogLik <- log(avg)

    # --- variable-coding constant-site correction ---
    # pConst = P(all 0) + P(all 1), ACRV-averaged FIRST, then -log(1-pConst).
    if (coding == 1L) {
      tip0 <- matrix(rep(c(1, 0), each = nTip), nrow = nTip)  # all tips state 0
      tip1 <- matrix(rep(c(0, 1), each = nTip), nrow = nTip)  # all tips state 1
      p0 <- 0
      p1 <- 0
      for (cat_idx in seq_len(nCat)) {
        tEffEdge <- neoEl * rates[cat_idx]
        p0 <- p0 + .ebe_prune_char(parent, child, tEffEdge, tip0,
                                   zVec, wEdge, pi1_base, rootFreqs,
                                   refEcology, phi, mode)
        p1 <- p1 + .ebe_prune_char(parent, child, tEffEdge, tip1,
                                   zVec, wEdge, pi1_base, rootFreqs,
                                   refEcology, phi, mode)
      }
      pConst <- (p0 + p1) / nCat
      charLogLik <- charLogLik - log(1 - pConst)
    }

    total <- total + charLogLik
  }

  total
}


# ---------------------------------------------------------------------------
# Standalone, package-DLL-free sanity helpers (used by the oracle SELF-tests
# in test-ebe-likelihood.R). These let us validate the oracle independent of
# the C++ kernel, which is being rewritten concurrently.
# ---------------------------------------------------------------------------

# Plain-R MkN pruner with NO ecology effect: a single 2-state asymmetric-
# equilibrium CTMC with lambda = 2, equilibrium (pi0_base, pi1_base), root
# (pi0_base, pi1_base). This is the z = 0 limit of the EBE kernel and is
# written from the textbook two-state formula, independent of .ebe_prune_char.
.ebe_baseline_mkn_loglik <- function(parent, child, edgeLen, tipStates,
                                     rate_loss, rate_neo = 1.0,
                                     rate_log_sd = 0, nCat = 1L,
                                     coding = 0L) {
  storage.mode(tipStates) <- "integer"
  nTip <- nrow(tipStates)
  nChar <- ncol(tipStates)
  maxNode <- 2L * nTip - 1L
  root <- nTip + 1L

  pi1 <- 1 / (1 + rate_loss)
  pi0 <- rate_loss / (1 + rate_loss)
  rootFreqs <- c(pi0, pi1)
  rates <- DiscreteLognormalRates(rate_log_sd, nCat)
  nCat <- length(rates)
  neoEl <- edgeLen * rate_neo

  Pmat <- function(tEff) {
    ex <- exp(-2 * tEff)
    matrix(c(pi0 + pi1 * ex, pi1 * (1 - ex),
             pi0 * (1 - ex), pi1 + pi0 * ex), nrow = 2L, byrow = TRUE)
  }

  pruneOne <- function(tipPartial, tEffEdge) {
    cl <- matrix(0, maxNode, 2L)
    init <- rep(FALSE, maxNode)
    cl[seq_len(nTip), ] <- tipPartial
    init[seq_len(nTip)] <- TRUE
    for (e in rev(seq_along(parent))) {
      pa <- parent[e]; ch <- child[e]
      msg <- as.vector(Pmat(tEffEdge[e]) %*% cl[ch, ])
      if (!init[pa]) cl[pa, ] <- msg else cl[pa, ] <- cl[pa, ] * msg
      init[pa] <- TRUE
    }
    sum(rootFreqs * cl[root, ])
  }

  total <- 0
  for (c in seq_len(nChar)) {
    tp <- matrix(0, nTip, 2L)
    for (i in seq_len(nTip)) {
      st <- tipStates[i, c]
      if (st < 0L) tp[i, ] <- 1 else tp[i, st + 1L] <- 1
    }
    sl <- 0
    for (cat_idx in seq_len(nCat)) sl <- sl + pruneOne(tp, neoEl * rates[cat_idx])
    cll <- log(sl / nCat)
    if (coding == 1L) {
      t0 <- matrix(rep(c(1, 0), each = nTip), nTip)
      t1 <- matrix(rep(c(0, 1), each = nTip), nTip)
      p0 <- 0; p1 <- 0
      for (cat_idx in seq_len(nCat)) {
        p0 <- p0 + pruneOne(t0, neoEl * rates[cat_idx])
        p1 <- p1 + pruneOne(t1, neoEl * rates[cat_idx])
      }
      cll <- cll - log(1 - (p0 + p1) / nCat)
    }
    total <- total + cll
  }
  total
}
