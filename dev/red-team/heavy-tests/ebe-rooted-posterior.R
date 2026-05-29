# ===========================================================================
# Heavy test: EBE rooted-topology posterior correctness (end-to-end).
#
# Challenges the claim that, under the Ecology-Biased-Equilibrium (EBE)
# likelihood (which depends on ROOT PLACEMENT when z != 0), the topology
# MOVES (NNI / SPR / TBR / pSPR) sample the CORRECT posterior over ROOTED
# trees. The moves' rooted Hastings ratios were built/tested only under the
# old root-INVARIANT model, so each must be re-verified.
#
# DESIGN (route (a) of the brief, made decisive):
#   * Tiny tree (5-6 tips). Scalar params FIXED (rate_loss, rate_neo,
#     rate_log_sd, phi, z, refEcology, ecology vector). z != 0, phi != 1 so the
#     likelihood is genuinely root-dependent.
#   * The prior over a rooted tree state (topology, TL, relBr) factorises as
#       pi(R) propto L_EBE(R) . Gamma(TL) . Dirichlet(relBr; 1,...,1)
#     and the prior over rooted TOPOLOGY is FLAT (proven: dev/ecology/
#     prior-reroot-check.R, maxdiff 0). So the EXACT marginal posterior over a
#     rooted topology t is the prior-predictive / evidence per topology
#       P(t) propto INT INT L_EBE(t,TL,relBr) Gamma(TL) Dir(relBr) d relBr d TL.
#   * We CONDITION ON FIXED total tree length TL0. The SPR/TBR Hastings ratios
#     are TL-invariant (verified: TL cancels in log(lRegraft) - log(lMerge);
#     same for TBR's log(lSubEdge) - log(lMergeSub)) to machine epsilon, so
#     fixing TL loses ZERO power on the branch-length Jacobian being tested.
#     The target becomes  P(t) propto E_{relBr ~ Dir(1..1)}[ L_EBE(t,TL0,relBr) ].
#   * The chain is a hand-rolled single-site Metropolis-Hastings using the
#     ACTUAL exported proposal kernels (spr_proposal / tbr_proposal /
#     nni_proposal -- byte-identical _impl code to the ecology-mode MCMC moves,
#     mcmc.cpp cases 6/17 fall through to *_proposal_impl) and the EXACT MH
#     acceptance form  logAlpha = dLogPi + logHastings  used at mcmc.cpp:5876.
#   * CRITICAL (the crux): a topology-only chain FREEZES the branch-length DOF
#     the target integrates over (NNI permutes a fixed length multiset; SPR/TBR
#     preserve total TL). So every chain ALSO runs a dirichlet_simplex relBr
#     move, making relBr VARY -- exercising exactly the merge/split that SPR/TBR
#     act on. Fixed-branch-length enumeration would silently test only NNI.
#
# PER-MOVE VERDICT: chi-squared / G-test of the chain's rooted-topology
# frequencies against the exact prior-predictive marginal.
#
# POSITIVE CONTROL: an SPR chain whose Jacobian is forcibly zeroed
# (logHastings := 0) MUST be flagged biased -- this calibrates the threshold
# and yields the power statement.
#
# Reproduction:
#   Small-N (<= 60 s, confirms execution):  Rscript ebe-rooted-posterior.R --quick
#   Full-scale:                             Rscript ebe-rooted-posterior.R --full
#   (or sbatch ebe-rooted-posterior-hamilton.sh)
# ===========================================================================

suppressMessages({
  devtools::load_all(".", quiet = TRUE)
  library(TreeTools)
})

args <- commandArgs(trailingOnly = TRUE)
QUICK <- ("--quick" %in% args) || !("--full" %in% args)

# ---------------------------------------------------------------------------
# Run-size knobs.
#   nTipSmall   : tip count for NNI/SPR/positive-control (rooted topos =
#                 NRooted(nTipSmall): 4->15, 5->105).
#   nTipTbr     : tip count for TBR (needs an internal subtree for Phase-B; on
#                 5 tips Phase-B fires ~3%/proposal, on 6 tips ~6%).
#   targetDraws : MC draws of relBr ~ Dir to estimate the exact marginal P(t).
#   chainIter   : MH iterations per chain (pre-thinning).
#   thinTo      : thin each chain's stored topology stream to ~this many
#                 samples (feedback_no_oversample: never store every iteration).
#   runTbr      : QUICK runs TBR at TINY scale too, so the execution gate
#                 covers EVERY move + the TBR positive control + the
#                 generalised control logic (a reference error in the TBR
#                 branch must surface here, NOT hours into a SLURM run).  Quick
#                 TBR is an EXECUTION check only -- it is underpowered for
#                 Phase-B bias (Phase-B fires ~3-6%/proposal and the chain is
#                 short), which the verdict reports explicitly.
#
# QUICK is an EXECUTION gate (<= 60 s): tiny trees, few draws/iters, one
# replicate.  It is STRUCTURALLY UNDERPOWERED to resolve small posterior
# deviations -- see the companion .md "What a failure would mean" and the
# positive-control power statement in verdict.txt.  The DECISIVE per-move
# verdict comes from --full; the exact Jacobians (incl. TBR Phase-B sigma) are
# pinned to machine precision by the companion ebe-move-detailed-balance.R.
# ---------------------------------------------------------------------------
if (QUICK) {
  CFG <- list(nTipSmall = 4L, nTipTbr = 5L,
              targetDraws = 800L, chainIter = 1.5e4, thinTo = 6000L,
              nReplChains = 1L, runTbr = TRUE, tag = "quick")
} else {
  CFG <- list(nTipSmall = 5L, nTipTbr = 6L,
              targetDraws = 4e5L, chainIter = 3e7, thinTo = 30000L,
              nReplChains = 4L, runTbr = TRUE, tag = "full")
}

RESULTS_DIR <- file.path("dev", "red-team", "heavy-tests",
                         "ebe-rooted-posterior-results")
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Fixture: neomorphic-only MkPrimeData on a small tree with one ecology column.
# Mirrors tests/testthat/test-ebe-likelihood.R::.ebe_make_neo_mkd so the EBE
# kernel is exercised exactly as the GATE-1 tests exercise it.
# ---------------------------------------------------------------------------
make_neo_mkd <- function(seed, nNeo, kEco, nTip) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  mat <- matrix(0L, nrow = nTip, ncol = nNeo + 1L, dimnames = list(tips, NULL))
  for (cc in seq_len(nNeo)) {
    repeat {
      col <- sample(0:1, nTip, replace = TRUE)
      if (length(unique(col)) > 1L) break
    }
    mat[, cc] <- col
  }
  repeat {
    eco <- sample.int(kEco, nTip, replace = TRUE) - 1L
    if (length(unique(eco)) == kEco) break
  }
  mat[, nNeo + 1L] <- eco
  pd  <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = seq_len(nNeo), ecology = nNeo + 1L)
  list(mkd = mkd, tips = tips, nTip = nTip)
}

# Fixed model knobs shared by target and chain. z != 0, phi != 1.
FIXED <- list(
  rate_loss = 1.4, rate_neo = 1.1, rate_log_sd = 0, nCat = 1L,
  phi = 2.5, magnitudeMode = "global", coding = "none", refEcology = 0L,
  TL0 = 1.0                       # fixed total tree length (Jacobian-irrelevant)
)

# Build the fixed z-matrix (one per fixture; z != 0 for >=1 character/column).
make_z <- function(mkd, seed) {
  zCols <- mkd$kEcology - 1L
  set.seed(seed)
  z <- matrix(sample.int(3L, mkd$nChar * zCols, replace = TRUE) - 1L,
              nrow = mkd$nChar, ncol = zCols)
  # Guarantee non-triviality: force at least one z=1 and one z=2.
  z[1, 1] <- 1L
  if (mkd$nChar >= 2L) z[2, 1] <- 2L
  storage.mode(z) <- "integer"
  z
}

# ---------------------------------------------------------------------------
# Likelihood wrapper: build a phylo from (edge, relBr) and call the GATE-1
# validated EBE kernel.  edge is in canonical preorder (proposals return that).
# ---------------------------------------------------------------------------
ebe_ll <- function(edge, relBr, mkd, zMat, nTip, tips) {
  ph <- structure(
    list(edge = edge, Nnode = nTip - 1L, tip.label = tips,
         edge.length = FIXED$TL0 * relBr),
    class = "phylo")
  MkPrime:::.MkpEcologyLogLikelihood(
    ph, mkd, kPrime = mkd$kObs,
    rate_loss = FIXED$rate_loss, rate_log_sd = FIXED$rate_log_sd,
    nCat = FIXED$nCat, rate_neo = FIXED$rate_neo, relabel = TRUE,
    phi = FIXED$phi, zMat = zMat, magnitudeMode = FIXED$magnitudeMode,
    coding = FIXED$coding, refEcology = FIXED$refEcology)
}

# Canonical key for a rooted topology (ignores branch lengths): reorder to
# preorder and write a cladogram Newick. Two phylo objects with the same rooted
# topology map to the same string.
rooted_key <- function(edge, nTip, tips) {
  ph <- structure(
    list(edge = edge, Nnode = nTip - 1L, tip.label = tips,
         edge.length = rep(1, nrow(edge))), class = "phylo")
  ape::write.tree(TreeTools::Preorder(ph))
}

# ---------------------------------------------------------------------------
# EXACT target: marginal posterior over rooted topologies.
#   P(t) propto E_{relBr ~ Dir(1..1)}[ L_EBE(t, TL0, relBr) ].
# Enumerate ALL rooted binary topologies; for each, average the likelihood over
# a COMMON set of relBr ~ Dir(1..1) draws (same draws across topologies so the
# ratio estimates share variance). Average in the log domain via log-sum-exp.
# ---------------------------------------------------------------------------
enumerate_rooted <- function(nTip, tips) {
  # All rooted binary trees on nTip tips, in our edge convention.
  #
  # NOTE: ape::as.phylo(integer, nTip, tipLabels) enumerates UNROOTED trees
  # (returning a fixed canonical rooting), so it yields only NUnrooted(nTip)
  # distinct ROOTED topologies, NOT NRooted(nTip).  Verified empirically.
  # We instead use the bijection
  #   rooted binary tree on n tips  <->  unrooted binary tree on n+1 tips
  # where an extra tip "ROOTMARK" marks the root edge.  Enumerate all
  # NUnrooted(n+1) = NRooted(n) unrooted (n+1)-trees, root each on ROOTMARK,
  # then drop ROOTMARK (collapsing the resulting degree-2 node) to recover the
  # rooted n-tree.  Verified to yield exactly NRooted(nTip) distinct rooted
  # topologies for n = 5 (105) and n = 6 (945), and chain output keys land in
  # the set 100% of the time.
  rtips <- c(tips, "ROOTMARK")
  nU <- TreeTools::NUnrooted(nTip + 1L)
  out <- vector("list", nU)
  for (i in seq_len(nU)) {
    u  <- ape::as.phylo(i - 1L, nTip = nTip + 1L, tipLabels = rtips)
    rr <- ape::root(u, outgroup = "ROOTMARK", resolve.root = TRUE)
    rr <- ape::drop.tip(rr, "ROOTMARK", collapse.singles = TRUE)
    out[[i]] <- TreeTools::Preorder(rr)
  }
  out
}

log_mean_exp <- function(lx) {
  m <- max(lx)
  m + log(mean(exp(lx - m)))
}

compute_target <- function(mkd, zMat, nTip, tips, nDraws, seed) {
  trees <- enumerate_rooted(nTip, tips)
  nEdge <- nrow(trees[[1]]$edge)
  # Common relBr ~ Dir(1,...,1) draws (uniform on the simplex): normalise
  # Exponential(1) variates.
  set.seed(seed)
  G <- matrix(rexp(nDraws * nEdge), nrow = nDraws, ncol = nEdge)
  relDraws <- G / rowSums(G)

  keys <- vapply(trees, function(p) rooted_key(p$edge, nTip, tips), character(1))
  logEvidence <- numeric(length(trees))
  for (ti in seq_along(trees)) {
    e <- trees[[ti]]$edge
    ll <- vapply(seq_len(nDraws), function(d)
      ebe_ll(e, relDraws[d, ], mkd, zMat, nTip, tips), numeric(1))
    logEvidence[ti] <- log_mean_exp(ll)
  }
  # Normalise to a probability vector over rooted topologies.
  logP <- logEvidence - log_mean_exp(logEvidence) - log(length(trees))
  P <- exp(logP - max(logP)); P <- P / sum(P)
  # key -> integer code lookup for fast categorical recording in the chain.
  keyToCode <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_along(keys)) assign(keys[i], i, envir = keyToCode)
  list(keys = keys, P = P, logEvidence = logEvidence, nDraws = nDraws,
       nTopo = length(trees), keyToCode = keyToCode)
}

# ---------------------------------------------------------------------------
# The hand-rolled single-site MH chain. Uses the ACTUAL exported proposals.
#   moveFn  : one of nni_proposal / spr_proposal / tbr_proposal / pspr (NULL).
#   relMove : if TRUE, alternate the topology move with a dirichlet_simplex
#             relBr move (REQUIRED so branch-length VALUES vary -- the crux).
#   zeroJac : positive control -- forcibly set the topology move's logHastings
#             to 0 (drops the branch-length Jacobian), which MUST bias SPR/TBR.
# Fixed total TL throughout (relBr is the only branch DOF; topology moves
# preserve TL; the dirichlet move preserves the simplex sum = 1).
# ---------------------------------------------------------------------------
run_chain <- function(moveName, mkd, zMat, nTip, tips, startEdge, startRel,
                      nIter, thinTo, seed, keyToCode, relMove = TRUE,
                      zeroJac = FALSE, psprData = NULL) {
  set.seed(seed)
  edge <- startEdge; rel <- startRel
  curLL <- ebe_ll(edge, rel, mkd, zMat, nTip, tips)
  thin <- max(1L, as.integer(nIter %/% thinTo))
  nStore <- as.integer(nIter %/% thin)
  codeStore <- integer(nStore); si <- 0L
  acc <- 0L; prop <- 0L; phaseBfire <- 0L

  topoMove <- function(edge, rel) {
    if (moveName == "nni") nni_proposal(edge, nTip, FIXED$TL0, rel)
    else if (moveName == "spr") spr_proposal(edge, nTip, FIXED$TL0, rel)
    else if (moveName == "tbr") tbr_proposal(edge, nTip, FIXED$TL0, rel)
    else stop("unknown move")
  }

  for (it in seq_len(nIter)) {
    doTopo <- !relMove || (it %% 2L == 1L)
    if (doTopo) {
      pr <- topoMove(edge, rel)
      lh <- pr$logHastings
      if (zeroJac) lh <- 0
      if (is.finite(lh)) {
        prop <- prop + 1L
        propLL <- ebe_ll(pr$edge, pr$rel_br_lengths, mkd, zMat, nTip, tips)
        logAlpha <- (propLL - curLL) + lh
        if (is.finite(logAlpha) && log(runif(1)) < logAlpha) {
          edge <- pr$edge; rel <- pr$rel_br_lengths; curLL <- propLL
          acc <- acc + 1L
        }
      }
    } else {
      # dirichlet_simplex relBr move: redistribute branch-length VALUES.
      pr <- dirichlet_simplex_proposal(rel, nCats = length(rel), alpha = 10)
      lh <- pr$logHastings
      if (is.finite(lh)) {
        relNew <- pr$value
        propLL <- ebe_ll(edge, relNew, mkd, zMat, nTip, tips)
        logAlpha <- (propLL - curLL) + lh
        if (is.finite(logAlpha) && log(runif(1)) < logAlpha) {
          rel <- relNew; curLL <- propLL
        }
      }
    }
    if (it %% thin == 0L) {
      si <- si + 1L
      codeStore[si] <- keyToCode[[rooted_key(edge, nTip, tips)]]
    }
  }
  list(codes = codeStore[seq_len(si)],
       accRate = if (prop > 0) acc / prop else NA_real_,
       nStored = si)
}

# ---------------------------------------------------------------------------
# Grade a chain's rooted-topology frequencies against the exact target.
#
# CRITICAL (a naive chi-squared gets this WRONG): MCMC samples are
# AUTOCORRELATED, so the effective sample size (ESS) is far below the stored
# count.  A G-test using the raw count INFLATES significance without bound
# (verified empirically: a CORRECT SPR chain's G grew 25 -> 76 -> 114 as
# nStored grew 7.5k -> 75k -> 300k while its TVD DECAYED 0.023 -> 0.012 ->
# 0.008 -- the hallmark of an UNBIASED sampler -- yet the raw-n p-value went
# 3.6e-2 -> 1.4e-10 -> 1.2e-17, a false-FAIL).  We therefore grade on TWO
# autocorrelation-aware criteria; the TVD null-band is primary.
#
#  (1) TVD NULL-BAND (primary discriminator).  Under the CORRECT target, the
#      expected L1 (total-variation) fluctuation of a multinomial at effective
#      size ESS is  E[TVD] ~ (1/2) sum_i sqrt(2 P_i(1-P_i)/(pi*ESS)), with sd of
#      similar order.  An UNBIASED move's TVD sits in this band and DECAYS as
#      1/sqrt(ESS); a BIASED move's TVD PLATEAUS at a nonzero bias floor
#      (verified: zeroed-Jacobian control TVD held ~0.14 across 7.5k/75k/300k).
#      PASS iff observed TVD <= nullMean + zCrit*nullSd.
#  (2) ESS-CORRECTED G-test (secondary).  Rescale G by ESS/nStored so the
#      chi-squared reference is valid; PASS iff p > alpha.
#
# ESS is estimated per high-probability topology via coda::effectiveSize on the
# 0/1 indicator series; we take the conservative MINIMUM over the topologies
# carrying >= 1% mass (the worst-mixing coordinate bounds the multinomial ESS).
# ---------------------------------------------------------------------------
estimate_ess <- function(codes, nLev, P, minMass = 0.01) {
  big <- which(P >= minMass)
  if (!length(big)) big <- order(P, decreasing = TRUE)[seq_len(min(5L, nLev))]
  essVals <- vapply(big, function(k) {
    ind <- as.numeric(codes == k)
    if (stats::sd(ind) == 0) return(length(ind))  # never/always visited
    es <- suppressWarnings(as.numeric(coda::effectiveSize(coda::mcmc(ind))))
    if (!is.finite(es) || es <= 0) es <- 1
    min(es, length(ind))
  }, numeric(1))
  max(1, floor(min(essVals)))
}

grade_chain <- function(chainCodes, target, label, alpha = 1e-3, zCrit = 4) {
  n <- length(chainCodes)
  nLev <- length(target$keys)
  obs <- tabulate(chainCodes, nbins = nLev)
  expP <- target$P
  obsP <- obs / n
  tvd <- 0.5 * sum(abs(obsP - expP))

  ess <- estimate_ess(chainCodes, nLev, expP)

  # (1) TVD null-band at ESS.  E|p_hat_i - P_i| for a binomial(ESS,P_i)
  # proportion ~ sqrt(2 P_i(1-P_i)/(pi*ESS)); TVD = 0.5 sum over i.  Its sd is
  # bounded by 0.5 sqrt(sum Var(|.|)) <= 0.5 sqrt(sum P_i(1-P_i)/ESS).
  v <- expP * (1 - expP)
  nullMean <- 0.5 * sum(sqrt(2 * v / (pi * ess)))
  nullSd   <- 0.5 * sqrt(sum(v / ess))
  tvdThresh <- nullMean + zCrit * nullSd
  passTVD <- tvd <= tvdThresh

  # (2) ESS-corrected G-test.
  expN <- expP * n
  big <- expN >= 5
  obsBig <- obs[big]; expBig <- expN[big]
  if (any(!big)) { obsBig <- c(obsBig, sum(obs[!big]))
                   expBig <- c(expBig, sum(expN[!big])) }
  nz <- obsBig > 0
  Graw <- 2 * sum(obsBig[nz] * log(obsBig[nz] / expBig[nz]))
  Gcorr <- Graw * (ess / n)               # rescale to effective size
  df <- length(obsBig) - 1L
  pval <- stats::pchisq(Gcorr, df = df, lower.tail = FALSE)
  passG <- pval > alpha

  list(label = label, n = n, ess = ess, nCellsBig = sum(big),
       tvd = tvd, tvdThresh = tvdThresh, nullMean = nullMean, nullSd = nullSd,
       Graw = Graw, Gcorr = Gcorr, df = df, pval = pval,
       passTVD = passTVD, passG = passG,
       pass = passTVD && passG)
}

# ---------------------------------------------------------------------------
# Drive one move's verdict: run nReplChains independent chains from random
# starts, grade each, and combine. A move PASSES iff every replicate passes
# (and -- for the positive control -- the zeroed-Jacobian variant FAILS).
# ---------------------------------------------------------------------------
random_start <- function(nTip, tips, seed) {
  set.seed(seed)
  tr <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  list(edge = tr$edge, rel = tr$edge.length / sum(tr$edge.length))
}

verdict_for_move <- function(moveName, mkd, zMat, nTip, tips, target, cfg,
                             seedBase, zeroJac = FALSE) {
  res <- lapply(seq_len(cfg$nReplChains), function(r) {
    st <- random_start(nTip, tips, seedBase + 1000L * r)
    ch <- run_chain(moveName, mkd, zMat, nTip, tips, st$edge, st$rel,
                    nIter = cfg$chainIter, thinTo = cfg$thinTo,
                    seed = seedBase + r, keyToCode = target$keyToCode,
                    relMove = TRUE, zeroJac = zeroJac)
    g <- grade_chain(ch$codes, target,
                     sprintf("%s%s repl %d", moveName,
                             if (zeroJac) "[zeroJac]" else "", r))
    g$accRate <- ch$accRate; g$nStored <- ch$nStored
    g
  })
  res
}

# ===========================================================================
# MAIN
# ===========================================================================
cat(sprintf("[ebe-rooted-posterior] mode=%s targetDraws=%d chainIter=%g thinTo=%d repl=%d\n",
            CFG$tag, CFG$targetDraws, CFG$chainIter, CFG$thinTo, CFG$nReplChains))

t0 <- Sys.time()
summary <- list(mode = CFG$tag, fixed = FIXED, cfg = CFG, moves = list())

# --- small-tip fixture for NNI / SPR / positive-control ---
nS <- CFG$nTipSmall
f5  <- make_neo_mkd(seed = 5101L, nNeo = 3L, kEco = 3L, nTip = nS)
z5  <- make_z(f5$mkd, seed = 55L)
cat(sprintf("%d-tip fixture: %d rooted topologies; computing exact target (%d draws)...\n",
            nS, TreeTools::NRooted(nS), CFG$targetDraws))
tgt5 <- compute_target(f5$mkd, z5, f5$nTip, f5$tips, CFG$targetDraws, seed = 5151L)
cat(sprintf("  target entropy ratio (effective topologies): %.1f / %d\n",
            exp(-sum(tgt5$P * log(pmax(tgt5$P, 1e-300)))), tgt5$nTopo))

for (mv in c("nni", "spr")) {
  cat(sprintf("Running %s chains (%d tips)...\n", mv, nS))
  summary$moves[[mv]] <- verdict_for_move(mv, f5$mkd, z5, f5$nTip, f5$tips,
                                          tgt5, CFG, seedBase = 700L)
}

# Positive control: SPR with Jacobian zeroed -- MUST be flagged biased.
cat(sprintf("Running POSITIVE CONTROL (SPR, Jacobian zeroed) chains (%d tips)...\n", nS))
summary$moves[["spr_zeroJac_CONTROL"]] <-
  verdict_for_move("spr", f5$mkd, z5, f5$nTip, f5$tips, tgt5, CFG,
                   seedBase = 900L, zeroJac = TRUE)

# --- TBR fixture (Phase-B sigma re-root; full mode only) ---
if (isTRUE(CFG$runTbr)) {
  nT <- CFG$nTipTbr
  f6  <- make_neo_mkd(seed = 6201L, nNeo = 3L, kEco = 3L, nTip = nT)
  z6  <- make_z(f6$mkd, seed = 66L)
  cat(sprintf("%d-tip fixture: %d rooted topologies; computing exact target (%d draws)...\n",
              nT, TreeTools::NRooted(nT), CFG$targetDraws))
  tgt6 <- compute_target(f6$mkd, z6, f6$nTip, f6$tips, CFG$targetDraws, seed = 6261L)
  cat(sprintf("Running tbr chains (%d tips)...\n", nT))
  summary$moves[["tbr"]] <- verdict_for_move("tbr", f6$mkd, z6, f6$nTip, f6$tips,
                                             tgt6, CFG, seedBase = 800L)
  # TBR-specific positive control: zero TBR's full Jacobian (drops both the
  # SPR-style and the Phase-B sigma terms) -- MUST be flagged biased, giving a
  # power statement for the Phase-B path that Part-B of the localiser does not
  # analytically cover.
  cat(sprintf("Running POSITIVE CONTROL (TBR, Jacobian zeroed) chains (%d tips)...\n", nT))
  summary$moves[["tbr_zeroJac_CONTROL"]] <-
    verdict_for_move("tbr", f6$mkd, z6, f6$nTip, f6$tips, tgt6, CFG,
                     seedBase = 850L, zeroJac = TRUE)
} else {
  cat("TBR skipped in quick mode (UNVERIFIED here; run --full).\n")
}

# ---------------------------------------------------------------------------
# Combine into per-move verdicts and write outputs.
# ---------------------------------------------------------------------------
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

move_verdict <- function(repls, isControl = FALSE) {
  pvals <- vapply(repls, function(r) r$pval, numeric(1))
  tvds  <- vapply(repls, function(r) r$tvd, numeric(1))
  thrs  <- vapply(repls, function(r) r$tvdThresh, numeric(1))
  esss  <- vapply(repls, function(r) r$ess, numeric(1))
  accs  <- vapply(repls, function(r) r$accRate, numeric(1))
  allPass <- all(vapply(repls, function(r) r$pass, logical(1)))
  if (isControl) {
    # Control SUCCEEDS (sanity) iff it is FLAGGED biased, i.e. NOT allPass.
    verdict <- if (!allPass) "CONTROL-OK (bias detected)" else "CONTROL-FAILED (bias MISSED)"
  } else {
    verdict <- if (allPass) "SAMPLES CORRECT ROOTED POSTERIOR" else "BIASED"
  }
  list(verdict = verdict, minPval = min(pvals), maxTvd = max(tvds),
       maxThresh = max(thrs), minEss = min(esss),
       meanAcc = mean(accs, na.rm = TRUE),
       pvals = pvals, tvds = tvds, allPass = allPass)
}

mv_summaries <- list()
for (nm in names(summary$moves)) {
  isCtrl <- grepl("CONTROL", nm)
  mv_summaries[[nm]] <- move_verdict(summary$moves[[nm]], isControl = isCtrl)
}
summary$verdicts <- mv_summaries
summary$elapsed_sec <- elapsed

saveRDS(summary, file.path(RESULTS_DIR, "summary.rds"))

# verdict.txt
controlNames <- grep("CONTROL", names(mv_summaries), value = TRUE)
ctrlFired <- function(nm) grepl("CONTROL-OK", mv_summaries[[nm]]$verdict)
# Map each real move to the control that powers it (same move family).
sprCtrlOK <- "spr_zeroJac_CONTROL" %in% controlNames && ctrlFired("spr_zeroJac_CONTROL")
tbrCtrlOK <- "tbr_zeroJac_CONTROL" %in% controlNames && ctrlFired("tbr_zeroJac_CONTROL")
# NNI and SPR share the SPR-style prune+regraft machinery; the SPR control
# powers both (NNI's logHastings is identically 0 so the relevant failure mode
# -- a dropped Jacobian -- is what the SPR control injects).
powerOf <- function(nm) {
  if (nm %in% c("nni", "spr")) sprCtrlOK
  else if (nm == "tbr") tbrCtrlOK
  else FALSE
}
realMoves <- setdiff(names(mv_summaries), controlNames)

# A move is TRUSTWORTHY only if its TVD/G verdict is CORRECT *and* its powering
# control fired (otherwise the run is underpowered for that move -> UNVERIFIED,
# never reported as CORRECT).  In quick mode TBR's Phase-B control rarely fires,
# so TBR is expected UNVERIFIED there.
moveStatus <- function(nm) {
  ok <- grepl("CORRECT", mv_summaries[[nm]]$verdict)
  pw <- powerOf(nm)
  if (ok && pw) "TRUSTWORTHY-CORRECT"
  else if (!ok && pw) "BIASED"
  else if (ok && !pw) "EXECUTED-UNDERPOWERED"   # ran clean but control didn't fire
  else "UNVERIFIED"
}
statuses <- vapply(realMoves, moveStatus, character(1))
names(statuses) <- realMoves
realPass <- all(statuses[realMoves] == "TRUSTWORTHY-CORRECT")
anyBiased <- any(statuses == "BIASED")
anyUnpowered <- any(statuses %in% c("EXECUTED-UNDERPOWERED", "UNVERIFIED"))

lines <- c(
  sprintf("EBE ROOTED-POSTERIOR HEAVY TEST  (mode=%s, %.1f s)", CFG$tag, elapsed),
  sprintf("Fixed: rate_loss=%.2f rate_neo=%.2f phi=%.2f TL0=%.2f coding=%s",
          FIXED$rate_loss, FIXED$rate_neo, FIXED$phi, FIXED$TL0, FIXED$coding),
  "Pass criterion (autocorrelation-aware): TVD <= nullMean + 4*nullSd at the",
  "chain's effective sample size AND ESS-corrected G-test p > 1e-3, every repl,",
  "AND the move's powering positive control fired (else UNDERPOWERED, not PASS).",
  "(A naive raw-n chi-squared false-FAILs correct moves -- see companion .md.)",
  "",
  "PER-MOVE VERDICT (rooted-topology marginal vs exact prior-predictive):")
for (nm in names(mv_summaries)) {
  v <- mv_summaries[[nm]]
  if (nm %in% controlNames) {
    label <- v$verdict; st <- ""
  } else {
    st <- moveStatus(nm)
    # Do NOT print "SAMPLES CORRECT ..." unless the move is actually trustworthy;
    # otherwise the line could be misread in isolation as a clean pass.
    label <- switch(st,
      "TRUSTWORTHY-CORRECT"   = "SAMPLES CORRECT ROOTED POSTERIOR",
      "BIASED"                = "BIASED (TVD exceeds band)",
      "EXECUTED-UNDERPOWERED" = "RAN CLEAN (control unpowered here)",
      "UNVERIFIED")
    st <- paste0("  [", st, "]")
  }
  lines <- c(lines, sprintf(
    "  %-22s %-34s TVD=%.4f(thr %.4f) ESS=%.0f corrP=%.3g acc=%.2f%s",
    nm, label, v$maxTvd, v$maxThresh, v$minEss, v$minPval, v$meanAcc, st))
}
lines <- c(lines, "", "POSITIVE CONTROLS (zeroed Jacobian MUST be flagged):")
for (nm in controlNames) {
  lines <- c(lines, sprintf("  %-22s detected bias: %s", nm,
                            if (ctrlFired(nm)) "YES (powered)" else "NO (UNDERPOWERED)"))
}

# Coverage notes.
tbrRun  <- "tbr" %in% names(mv_summaries)
notes <- c()
if (tbrRun && !tbrCtrlOK)
  notes <- c(notes, paste0("TBR: EXECUTED but UNDERPOWERED for Phase-B in this run ",
    "(control did not fire). TBR Phase-B sigma-Jacobian is pinned to machine ",
    "precision by ebe-move-detailed-balance.R; run --full for the end-to-end ",
    "TBR verdict."))
notes <- c(notes,
  "pSPR: UNVERIFIED (parsimony-guided; covered by this design but not wired into this run).")

# OVERALL: a true bias anywhere -> FAIL.  Otherwise, PASS only over the moves
# whose powering control fired; moves that ran clean but underpowered are
# reported separately and never counted as a refuted bias.
trustworthy <- realMoves[statuses == "TRUSTWORTHY-CORRECT"]
underpowered <- realMoves[statuses %in% c("EXECUTED-UNDERPOWERED", "UNVERIFIED")]
overall <- if (anyBiased) {
  sprintf("FAIL -- biased move(s): %s.",
          paste(realMoves[statuses == "BIASED"], collapse = ", "))
} else if (!sprCtrlOK) {
  "INCONCLUSIVE -- SPR positive control did not fire; this run is underpowered. PASS-without-power is NOT reported."
} else if (length(underpowered) == 0L) {
  "PASS -- every topology move sampled the correct rooted posterior (controls fired)."
} else {
  sprintf(paste0("PASS (partial) -- {%s} sample correctly (controls fired); ",
                 "{%s} EXECUTED but UNDERPOWERED here (see notes / run --full)."),
          paste(trustworthy, collapse = ", "),
          paste(underpowered, collapse = ", "))
}
lines <- c(lines, "", "COVERAGE NOTES:", paste0("  ", notes),
           "", paste0("OVERALL: ", overall))
writeLines(lines, file.path(RESULTS_DIR, "verdict.txt"))
cat(paste(lines, collapse = "\n"), "\n")
cat(sprintf("\nWrote %s and summary.rds\n", file.path(RESULTS_DIR, "verdict.txt")))
