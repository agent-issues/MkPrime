# EBE (Ecology-Biased Equilibrium) neomorphic likelihood — spec test contracts.
#
# Two tiers:
#   (A) Oracle SELF-tests — pure R, no package DLL. These guard the independent
#       oracle in helper-ebe-oracle.R and run regardless of kernel state. They
#       cover the oracle side of T3 (re-root) and T4 (normalisation), and the
#       z=0 graceful-degradation reduction.
#   (B) C++ GATE tests — cross-check the EBE C++ kernel against the oracle /
#       baseline (T1, T2, T3, T4). These are RUN AT THE INTEGRATION GATE by the
#       orchestrator after the DLL is rebuilt to EBE semantics; they will error
#       at collection time if the kernel has not yet been rebuilt, which is the
#       intended signal.
#
# The oracle (ebe_loglik_R) is derived from dev/ecology/ebe-spec.md section 2
# ALONE; no C++ kernel source was read. T2 is therefore a genuine independent
# check, not a tautology.

library("TreeTools")


# ---------------------------------------------------------------------------
# Shared tiny fixtures.
# ---------------------------------------------------------------------------

# Rooted 4-tip tree ((t1,t2),(t3,t4)); ape ids: tips 1..4, root 5,
# internal 6 = parent(1,2), 7 = parent(3,4).
.ebe_treeA <- function() {
  list(
    parent  = c(5L, 5L, 6L, 6L, 7L, 7L),
    child   = c(6L, 7L, 1L, 2L, 3L, 4L),
    edgeLen = c(0.4, 0.5, 0.6, 0.3, 0.7, 0.2)
  )
}

# The SAME unrooted tree re-rooted at old internal node 6, renumbered so the
# root is node nTip+1 = 5. Every PHYSICAL edge (and thus its wEdge row) is
# preserved; only the root placement moves. Used for the re-root contract.
#   new 5 -> 1, 5 -> 2, 5 -> 6 ; 6 -> 7 ; 7 -> 3, 7 -> 4
# physical-edge map vs treeA edges (1=5->6,2=5->7,3=6->1,4=6->2,5=7->3,6=7->4):
#   5->1 = old 6->1 (row 3); 5->2 = old 6->2 (row 4); 5->6 = old 5->6 (row 1);
#   6->7 = old 5->7 (row 2); 7->3 = old 7->3 (row 5); 7->4 = old 7->4 (row 6)
.ebe_treeB <- function() {
  list(
    parent  = c(5L, 5L, 5L, 6L, 7L, 7L),
    child   = c(1L, 2L, 6L, 7L, 3L, 4L),
    edgeLen = c(0.6, 0.3, 0.4, 0.5, 0.7, 0.2),
    permute = c(3L, 4L, 1L, 2L, 5L, 6L)   # row reorder of treeA wEdge
  )
}

# Hand-built valid stochastic wEdge (kEco = 2): rows >= 0, sum to 1.
.ebe_wEdgeA <- function() {
  w <- matrix(c(0.8, 0.2,
                0.3, 0.7,
                0.5, 0.5,
                0.9, 0.1,
                0.2, 0.8,
                0.6, 0.4), ncol = 2L, byrow = TRUE)
  stopifnot(all(abs(rowSums(w) - 1) < 1e-12))
  w
}


# ===========================================================================
# (A) Oracle SELF-tests — pure R, independent of the C++ kernel.
# ===========================================================================

test_that("oracle: z=0 reduces to baseline MkN for any phi (graceful degradation)", {
  trA <- .ebe_treeA()
  wA  <- .ebe_wEdgeA()
  rate_loss <- 1.7
  tipStates <- matrix(c(0L, 1L, 0L, 1L,
                        1L, 1L, 0L, 0L,
                        0L, 0L, 0L, 1L), nrow = 4L, ncol = 3L)
  zZero <- matrix(0L, nrow = 3L, ncol = 1L)

  for (sd in c(0, 0.4)) {
    for (cod in c(0L, 1L)) {
      ll_ebe <- ebe_loglik_R(
        trA$parent, trA$child, trA$edgeLen, tipStates,
        rate_loss = rate_loss, phi = 2.5, zMat = zZero, wEdge = wA,
        refEcology = 0L, rate_neo = 1.3, rate_log_sd = sd, nCat = 4L,
        mode = 0L, coding = cod
      )
      ll_base <- .ebe_baseline_mkn_loglik(
        trA$parent, trA$child, trA$edgeLen, tipStates,
        rate_loss = rate_loss, rate_neo = 1.3,
        rate_log_sd = sd, nCat = 4L, coding = cod
      )
      expect_equal(ll_ebe, ll_base, tolerance = 1e-12)
    }
  }
})


test_that("oracle T4: sum over all 2^nTip patterns of P(pattern) == 1 (coding=0)", {
  trA <- .ebe_treeA()
  wA  <- .ebe_wEdgeA()
  nTip <- 4L
  pats <- expand.grid(rep(list(0:1), nTip))

  enum_total <- function(zMat, sd, nCat) {
    s <- 0
    for (r in seq_len(nrow(pats))) {
      ts <- matrix(as.integer(unlist(pats[r, ])), nrow = nTip, ncol = 1L)
      s <- s + exp(ebe_loglik_R(
        trA$parent, trA$child, trA$edgeLen, ts,
        rate_loss = 1.7, phi = 2.5, zMat = zMat, wEdge = wA,
        refEcology = 0L, rate_neo = 1.3, rate_log_sd = sd, nCat = nCat,
        mode = 0L, coding = 0L
      ))
    }
    s
  }

  zZero <- matrix(0L, nrow = 1L, ncol = 1L)
  zMix  <- matrix(1L, nrow = 1L, ncol = 1L)   # z = 1
  expect_equal(enum_total(zZero, 0, 1L),   1, tolerance = 1e-10)
  expect_equal(enum_total(zMix,  0, 1L),   1, tolerance = 1e-10)
  expect_equal(enum_total(zMix,  0.4, 4L), 1, tolerance = 1e-10)
})


test_that("oracle T3: re-rooting is invariant at z=0, changes at z!=0", {
  trA <- .ebe_treeA(); trB <- .ebe_treeB()
  wA  <- .ebe_wEdgeA(); wB <- wA[trB$permute, ]
  rate_loss <- 1.7
  tipStates <- matrix(c(0L, 1L, 0L, 1L,
                        1L, 1L, 0L, 0L,
                        0L, 0L, 0L, 1L), nrow = 4L, ncol = 3L)
  zZero <- matrix(0L, nrow = 3L, ncol = 1L)
  zMix  <- matrix(c(1L, 2L, 0L), nrow = 3L, ncol = 1L)

  llA0 <- ebe_loglik_R(trA$parent, trA$child, trA$edgeLen, tipStates,
                       rate_loss, 2.5, zZero, wEdge = wA, refEcology = 0L)
  llB0 <- ebe_loglik_R(trB$parent, trB$child, trB$edgeLen, tipStates,
                       rate_loss, 2.5, zZero, wEdge = wB, refEcology = 0L)
  expect_equal(llA0, llB0, tolerance = 1e-12)

  llAm <- ebe_loglik_R(trA$parent, trA$child, trA$edgeLen, tipStates,
                       rate_loss, 2.5, zMix, wEdge = wA, refEcology = 0L)
  llBm <- ebe_loglik_R(trB$parent, trB$child, trB$edgeLen, tipStates,
                       rate_loss, 2.5, zMix, wEdge = wB, refEcology = 0L)
  expect_false(isTRUE(all.equal(llAm, llBm, tolerance = 1e-6)))
})


# ===========================================================================
# (B) C++ GATE tests — cross-check the EBE kernel against oracle / baseline.
#
# These call the EBE entry points (.MkpEcologyLogLikelihood, .PruningMknEcology)
# whose semantics are being rewritten concurrently. They are executed at the
# integration gate after the DLL is rebuilt to EBE semantics.
# ===========================================================================

# Build a neomorphic-only MkPrimeData on a 6-tip tree with one ecology column,
# so the EBE neomorphic kernel is exercised in isolation (no JC partition).
.ebe_make_neo_mkd <- function(seed = 2026L, nNeo = 3L, kEco = 3L,
                              nTip = 6L) {
  set.seed(seed)
  tips <- paste0("t", seq_len(nTip))
  mat <- matrix(0L, nrow = nTip, ncol = nNeo + 1L,
                dimnames = list(tips, NULL))
  # Neomorphic binary chars — ensure each is variable (avoid all-constant).
  for (c in seq_len(nNeo)) {
    repeat {
      col <- sample(0:1, nTip, replace = TRUE)
      if (length(unique(col)) > 1L) break
    }
    mat[, c] <- col
  }
  # Ecology column with kEco states, all states present.
  repeat {
    eco <- sample.int(kEco, nTip, replace = TRUE) - 1L
    if (length(unique(eco)) == kEco) break
  }
  mat[, nNeo + 1L] <- eco
  pd <- TreeTools::MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd, neomorphic = seq_len(nNeo), ecology = nNeo + 1L)
  tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  list(mkd = mkd, tree = tree, nNeo = nNeo, kEco = kEco)
}


# Reconstruct (parent, child, neoEl, wEdge, tipStates) for the oracle the SAME
# way the orchestrator does, so a T2 diff isolates the tilt math.
.ebe_recon <- function(tree, mkd, rate_neo) {
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  edgeLength <- tree$edge.length
  ecologyTip <- as.integer(mkd$ecology)
  ecologyTip[is.na(ecologyTip)] <- -1L
  marg  <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLength,
                                           ecologyTip, mkd$kEcology)
  wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)
  # neomorphic partition tip states (assume single neomorphic partition).
  neoPart <- Filter(function(p) p$type == "neomorphic", mkd$partitions)[[1]]
  tipStates <- neoPart$tip_states
  tipStates[is.na(tipStates)] <- -1L
  storage.mode(tipStates) <- "integer"
  list(parent = parent, child = child, edgeLength = edgeLength,
       wEdge = wEdge, tipStates = tipStates,
       neoEl = edgeLength * rate_neo)
}


test_that("T1: all z=0 EBE == baseline non-ecology MkN (coding 0 and 1)", {
  f <- .ebe_make_neo_mkd(seed = 11L)
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)

  for (cod in c("none", "variable")) {
    ll_base <- MkPrime:::.MkpLogLikelihood(
      f$tree, f$mkd, kPrime = f$mkd$kObs,
      rate_loss = 1.3, rate_log_sd = 0, nCat = 1L,
      coding = cod, rate_neo = 1.0, relabel = TRUE
    )
    # phi != 1 deliberately: at z=0 the tilt must vanish regardless of phi.
    ll_ebe <- MkPrime:::.MkpEcologyLogLikelihood(
      f$tree, f$mkd, kPrime = f$mkd$kObs,
      rate_loss = 1.3, rate_log_sd = 0, nCat = 1L,
      rate_neo = 1.0, relabel = TRUE,
      phi = 2.5, zMat = zMat, magnitudeMode = "global",
      coding = cod, refEcology = 0L
    )
    expect_equal(ll_ebe, ll_base, tolerance = 1e-10)
  }
})


test_that("T1 (ACRV): all z=0 EBE == baseline MkN with rate variation", {
  f <- .ebe_make_neo_mkd(seed = 12L)
  zMat <- matrix(0L, nrow = f$mkd$nChar, ncol = f$mkd$kEcology - 1L)
  for (cod in c("none", "variable")) {
    ll_base <- MkPrime:::.MkpLogLikelihood(
      f$tree, f$mkd, kPrime = f$mkd$kObs,
      rate_loss = 0.8, rate_log_sd = 0.5, nCat = 4L,
      coding = cod, rate_neo = 1.2, relabel = TRUE
    )
    ll_ebe <- MkPrime:::.MkpEcologyLogLikelihood(
      f$tree, f$mkd, kPrime = f$mkd$kObs,
      rate_loss = 0.8, rate_log_sd = 0.5, nCat = 4L,
      rate_neo = 1.2, relabel = TRUE,
      phi = 3.0, zMat = zMat, magnitudeMode = "global",
      coding = cod, refEcology = 0L
    )
    expect_equal(ll_ebe, ll_base, tolerance = 1e-10)
  }
})


test_that("T2: mixed z, phi != 1 — C++ EBE == R oracle (coding 0/1, +/- ACRV)", {
  f <- .ebe_make_neo_mkd(seed = 21L)
  zCols <- f$mkd$kEcology - 1L
  set.seed(7)
  zMat <- matrix(sample.int(3L, f$mkd$nChar * zCols, replace = TRUE) - 1L,
                 nrow = f$mkd$nChar, ncol = zCols)
  storage.mode(zMat) <- "integer"
  phi <- 2.5
  rate_loss <- 1.4
  rate_neo  <- 1.1

  grid <- list(
    list(sd = 0,   nCat = 1L),
    list(sd = 0.5, nCat = 4L)
  )
  for (g in grid) {
    for (cod in c("none", "variable")) {
      ll_cpp <- MkPrime:::.MkpEcologyLogLikelihood(
        f$tree, f$mkd, kPrime = f$mkd$kObs,
        rate_loss = rate_loss, rate_log_sd = g$sd, nCat = g$nCat,
        rate_neo = rate_neo, relabel = TRUE,
        phi = phi, zMat = zMat, magnitudeMode = "global",
        coding = cod, refEcology = 0L
      )
      rec <- .ebe_recon(f$tree, f$mkd, rate_neo)
      ll_r <- ebe_loglik_R(
        rec$parent, rec$child, rec$edgeLength, rec$tipStates,
        rate_loss = rate_loss, phi = phi, zMat = zMat,
        wEdge = rec$wEdge, refEcology = 0L, rate_neo = rate_neo,
        rate_log_sd = g$sd, nCat = g$nCat, mode = 0L,
        coding = if (cod == "variable") 1L else 0L
      )
      expect_equal(ll_cpp, ll_r, tolerance = 1e-10,
                   info = sprintf("coding=%s sd=%g nCat=%d", cod, g$sd, g$nCat))
    }
  }
})


test_that("T2 (per_ecology phi): C++ EBE == R oracle", {
  f <- .ebe_make_neo_mkd(seed = 22L, kEco = 3L)
  zCols <- f$mkd$kEcology - 1L
  set.seed(8)
  zMat <- matrix(sample.int(3L, f$mkd$nChar * zCols, replace = TRUE) - 1L,
                 nrow = f$mkd$nChar, ncol = zCols)
  storage.mode(zMat) <- "integer"
  phi <- c(1.5, 2.0, 0.7)   # length kEco; phi[ecoState + 1]
  rate_loss <- 1.2; rate_neo <- 1.0

  for (cod in c("none", "variable")) {
    ll_cpp <- MkPrime:::.MkpEcologyLogLikelihood(
      f$tree, f$mkd, kPrime = f$mkd$kObs,
      rate_loss = rate_loss, rate_log_sd = 0, nCat = 1L,
      rate_neo = rate_neo, relabel = TRUE,
      phi = phi, zMat = zMat, magnitudeMode = "per_ecology",
      coding = cod, refEcology = 0L
    )
    rec <- .ebe_recon(f$tree, f$mkd, rate_neo)
    ll_r <- ebe_loglik_R(
      rec$parent, rec$child, rec$edgeLength, rec$tipStates,
      rate_loss = rate_loss, phi = phi, zMat = zMat,
      wEdge = rec$wEdge, refEcology = 0L, rate_neo = rate_neo,
      rate_log_sd = 0, nCat = 1L, mode = 1L,
      coding = if (cod == "variable") 1L else 0L
    )
    expect_equal(ll_cpp, ll_r, tolerance = 1e-10, info = cod)
  }
})


test_that("T2 (raw kernel localiser): .PruningMknEcology == oracle (coding 0, nCat 1)", {
  # Calls the bare neomorphic pruner so a T2 failure localises to the kernel
  # vs the orchestration. coding=0, nCat=1 only (the raw pruner has no
  # const-site / ACRV-averaging wrapper of its own here).
  f <- .ebe_make_neo_mkd(seed = 23L)
  rate_loss <- 1.4; rate_neo <- 1.0
  rec <- .ebe_recon(f$tree, f$mkd, rate_neo)
  zCols <- f$mkd$kEcology - 1L
  set.seed(9)
  # neomorphic partition z-rows
  neoIdx <- Filter(function(p) p$type == "neomorphic",
                   f$mkd$partitions)[[1]]$char_indices
  zMatFull <- matrix(sample.int(3L, f$mkd$nChar * zCols, replace = TRUE) - 1L,
                     nrow = f$mkd$nChar, ncol = zCols)
  storage.mode(zMatFull) <- "integer"
  zPart <- zMatFull[neoIdx, , drop = FALSE]
  storage.mode(zPart) <- "integer"

  rootFreqs <- as.numeric(mkn_stationary_freqs(rate_loss))
  ll_cpp <- MkPrime:::.PruningMknEcology(
    rec$parent, rec$child, rec$neoEl, rec$tipStates,
    rate_loss, rootFreqs, 1,
    rec$wEdge, zPart, phi = 2.5, mode = 0L, refEcology = 0L
  )
  ll_r <- ebe_loglik_R(
    rec$parent, rec$child, rec$edgeLength, rec$tipStates,
    rate_loss = rate_loss, phi = 2.5, zMat = zPart,
    wEdge = rec$wEdge, refEcology = 0L, rate_neo = rate_neo,
    rate_log_sd = 0, nCat = 1L, mode = 0L, coding = 0L
  )
  expect_equal(ll_cpp, ll_r, tolerance = 1e-10)
})


test_that("T3 (C++): re-rooting changes EBE logLik when z!=0, invariant when z=0", {
  f <- .ebe_make_neo_mkd(seed = 31L)
  zCols <- f$mkd$kEcology - 1L

  # Re-root the SAME topology on a different outgroup tip; branch lengths and
  # topology are otherwise preserved, so only root placement differs. Try each
  # tip as outgroup until the edge table actually changes (a no-op re-root —
  # when the chosen tip is already root-adjacent — would spuriously fail the
  # z!=0 "changes" assertion).
  tips <- f$tree$tip.label
  treeA <- f$tree
  treeB <- NULL
  for (og in tips) {
    cand <- TreeTools::Preorder(ape::root(f$tree, outgroup = og,
                                          resolve.root = TRUE))
    if (!isTRUE(all.equal(cand$edge, treeA$edge)) ||
        !isTRUE(all.equal(cand$edge.length, treeA$edge.length))) {
      treeB <- cand
      break
    }
  }
  expect_false(is.null(treeB))   # a distinct rooting must exist on >=4 tips

  zZero <- matrix(0L, nrow = f$mkd$nChar, ncol = zCols)
  set.seed(3)
  zMix <- matrix(sample.int(3L, f$mkd$nChar * zCols, replace = TRUE) - 1L,
                 nrow = f$mkd$nChar, ncol = zCols)
  storage.mode(zMix) <- "integer"

  args0 <- list(kPrime = f$mkd$kObs, rate_loss = 1.0, rate_log_sd = 0,
                nCat = 1L, rate_neo = 1.0, relabel = TRUE,
                magnitudeMode = "global", coding = "none", refEcology = 0L)
  llA_z0 <- do.call(MkPrime:::.MkpEcologyLogLikelihood,
                    c(list(treeA, f$mkd, phi = 2.5, zMat = zZero), args0))
  llB_z0 <- do.call(MkPrime:::.MkpEcologyLogLikelihood,
                    c(list(treeB, f$mkd, phi = 2.5, zMat = zZero), args0))
  expect_equal(llA_z0, llB_z0, tolerance = 1e-10)

  llA_zm <- do.call(MkPrime:::.MkpEcologyLogLikelihood,
                    c(list(treeA, f$mkd, phi = 2.5, zMat = zMix), args0))
  llB_zm <- do.call(MkPrime:::.MkpEcologyLogLikelihood,
                    c(list(treeB, f$mkd, phi = 2.5, zMat = zMix), args0))
  expect_false(isTRUE(all.equal(llA_zm, llB_zm, tolerance = 1e-6)))
})


test_that("T4 (C++): sum over all 2^nTip tip patterns of P(pattern) == 1", {
  # PLUMBING-only normalisation check on the RAW neomorphic kernel
  # (.PruningMknEcology), not the orchestrator: this keeps the transition
  # model byte-identical across all 2^nTip evaluations (only tip states vary)
  # and avoids any data-dependent MkPrimeData reconstruction on the all-0 /
  # all-1 constant patterns. Holds for ANY stochastic P (spec T4), so a single
  # z=1 character + non-trivial phi suffices. coding is implicitly "none".
  nTip <- 5L
  tips <- paste0("t", seq_len(nTip))
  set.seed(41)
  tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = tips))
  parent <- tree$edge[, 1]
  child  <- tree$edge[, 2]
  edgeLength <- tree$edge.length

  # Fixed ecology vector → reconstruct wEdge ONCE (kEco = 3, all states used).
  repeat {
    eco <- sample.int(3L, nTip, replace = TRUE) - 1L
    if (length(unique(eco)) == 3L) break
  }
  marg  <- MkPrime:::.EcologyNodeMarginals(parent, child, edgeLength, eco, 3L)
  wEdge <- MkPrime:::.EcologyEdgeWeights(marg, parent, child)

  rate_loss <- 1.3
  rate_neo  <- 1.0
  neoEl <- edgeLength * rate_neo
  rootFreqs <- as.numeric(mkn_stationary_freqs(rate_loss))
  zMat1 <- matrix(1L, nrow = 1L, ncol = 3L - 1L)   # one char, z = 1

  pats <- expand.grid(rep(list(0:1), nTip))
  total <- 0
  for (r in seq_len(nrow(pats))) {
    ts <- matrix(as.integer(unlist(pats[r, ])), ncol = 1L)
    total <- total + exp(MkPrime:::.PruningMknEcology(
      parent, child, neoEl, ts, rate_loss, rootFreqs, 1,
      wEdge, zMat1, phi = 2.0, mode = 0L, refEcology = 0L
    ))
  }
  expect_equal(total, 1, tolerance = 1e-9)
})
