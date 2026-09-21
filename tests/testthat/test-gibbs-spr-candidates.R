# Deterministic tests for the corrected gibbs_spr kernel (GSPR-001/GSPR-004).
#
# The beta = 0 regression gate (dev/red-team/heavy-tests/gibbs-spr-db.R) is
# structurally blind to the two things the fix adds: at beta = 0 every
# selection weight is exp(0) = 1, so the selection ratio is identically zero
# and the GSPR-004 candidate-set asymmetry is invisible.  These tests settle
# both deterministically, with no sampling:
#
#   1. Pruning the same subtree from x and from any proposal y enumerates the
#      SAME set of edges of the shared residual tree R, and every edge's
#      tau = 1/2 reference weight agrees between the two directions — in
#      particular x's own (merged-edge) weight, computed through the
#      mergedEdge special case, equals y's regular-path evaluation of that
#      same edge.  That is the normaliser cancellation the MH ratio relies on.
#   2. The commit path recomputes state logPrior (retires GSPR-002).
#   3. Accepted moves no longer land in the pi-null set of bit-equal split
#      halves, and drawing the merged edge is a live branch-fraction move.

library("ape", quietly = TRUE)
library("TreeTools", quietly = TRUE)

# ---------------------------------------------------------------------------
# Fixture helpers
# ---------------------------------------------------------------------------

.CandContext <- function(seed = 7L, nTip = 6L, nChar = 5L,
                         model = MkPrimeModel()) {
  set.seed(seed)
  tree <- Preorder(ape::rtree(nTip, rooted = FALSE))
  mat  <- matrix(sample(0:2, nTip * nChar, replace = TRUE), nrow = nTip,
                 dimnames = list(tree$tip.label, NULL))
  mkd   <- suppressWarnings(MkPrimeData(MatrixToPhyDat(mat)))
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  list(mkd = mkd, model = model,
       dataPtr = MkPrime:::.InitMcmcData(mkd, model),
       tree = tree, nTip = nTip)
}

.StateFor <- function(ctx, tree) {
  sp <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, ctx$mkd,
                                                      ctx$model))
  fill_partition_cache(ctx$dataPtr, sp)
  allocate_cl_workspace(ctx$dataPtr, sp)
  sp
}

# Recursive tip-set lookup: returns a function node -> sorted tip ids below it
.TipSetBelow <- function(edge, nTip) {
  kids <- split(edge[, 2L], edge[, 1L])
  rec <- function(node) {
    if (node <= nTip) return(node)
    sort(unlist(lapply(kids[[as.character(node)]], rec), use.names = FALSE))
  }
  rec
}

# R-edge identity of an enumerated edge: child-side tip set minus the pruned
# subtree's tips.  Robust to the node renumbering Preorder may apply.
.REdgeKeys <- function(edges, tipSetFn, pruneTips) {
  apply(edges, 1L, function(e)
    paste(setdiff(tipSetFn(e[[2L]]), pruneTips), collapse = ","))
}

# Apply the SPR that gibbs_spr commits (regraft (u -> v) onto candRow at tau),
# then canonicalize — the R-side mirror of gibbs_spr_finish's proposal build.
.ApplySpr <- function(tree, pruneRow, candRow, tau) {
  edge <- tree$edge
  len  <- tree$edge.length
  u <- edge[pruneRow, 1L]
  v <- edge[pruneRow, 2L]
  parentRow <- which(edge[, 2L] == u)
  sibRow    <- which(edge[, 1L] == u & edge[, 2L] != v)
  sibNode   <- edge[sibRow, 2L]
  lMerge    <- len[parentRow] + len[sibRow]
  b    <- edge[candRow, 2L]
  lReg <- len[candRow]
  edge[parentRow, 2L] <- sibNode
  len[parentRow]      <- lMerge
  edge[candRow, 2L]   <- u
  len[candRow]        <- tau * lReg
  edge[sibRow, ]      <- c(u, b)
  len[sibRow]         <- (1 - tau) * lReg
  tree$edge        <- edge
  tree$edge.length <- len
  Preorder(tree)
}

# ---------------------------------------------------------------------------
# 1. GSPR-004: candidate-set symmetry and weight cancellation
# ---------------------------------------------------------------------------

test_that("GSPR-004: x and y enumerate the same edges of R with equal weights", {
  ctx  <- .CandContext(seed = 7L)
  nTip <- ctx$nTip
  root <- nTip + 1L
  x    <- ctx$tree

  # Two prune edges: one whose merge point g is the trifurcating root, one
  # deeper — the merged-edge special case must hold in both regimes.
  eligible <- which(x$edge[, 1L] != root)
  gOf <- x$edge[match(x$edge[eligible, 1L], x$edge[, 2L]), 1L]
  pruneRows <- unique(c(eligible[which(gOf == root)[1L]],
                        eligible[which(gOf != root)[1L]]))
  pruneRows <- pruneRows[!is.na(pruneRows)]
  expect_gte(length(pruneRows), 2L)

  for (pruneRow in pruneRows) {
    spX  <- .StateFor(ctx, x)
    enX  <- gibbs_spr_enumerate_cpp(ctx$dataPtr, spX, pruneRow)
    recX <- .TipSetBelow(x$edge, nTip)
    S    <- recX(enX$v)
    keyX <- .REdgeKeys(enX$edges, recX, S)
    nX   <- length(keyX)

    # Enumerated set is exactly E(R): every edge of the pruned tree, once.
    nDesc <- sum(vapply(x$edge[, 2L],
                        function(ch) all(recX(ch) %in% S), logical(1L)))
    expect_equal(nX, nrow(x$edge) - nDesc - 1L)
    expect_false(anyDuplicated(keyX) > 0L)

    # Build y by regrafting onto the first candidate at an asymmetric tau
    j <- 1L
    candRow <- which(x$edge[, 1L] == enX$edges[j, 1L] &
                     x$edge[, 2L] == enX$edges[j, 2L])
    expect_length(candRow, 1L)
    y    <- .ApplySpr(x, pruneRow, candRow, tau = 0.37)
    spY  <- .StateFor(ctx, y)
    recY <- .TipSetBelow(y$edge, nTip)

    # Locate the same prune edge (u -> v) in y by v's tip set
    pruneRowY <- which(vapply(seq_len(nrow(y$edge)), function(i)
      identical(recY(y$edge[i, 2L]), S), logical(1L)))
    expect_length(pruneRowY, 1L)

    enY  <- gibbs_spr_enumerate_cpp(ctx$dataPtr, spY, pruneRowY)
    keyY <- .REdgeKeys(enY$edges, recY, S)

    # Same set of edges of R from both endpoints — the GSPR-004 fix
    expect_equal(length(keyY), nX)
    expect_setequal(keyX, keyY)

    # Same weight for every edge of R, whichever endpoint enumerates it.
    # This includes x's own position: the mergedEdge special case (last row
    # of enX) must equal y's regular-path evaluation of that edge, and vice
    # versa — the normaliser cancellation is exact, not assumed.
    oy <- match(keyX, keyY)
    expect_false(anyNA(oy))
    expect_equal(enX$logLik, enY$logLik[oy], tolerance = 1e-8)

    # The merged edge of x is a regular candidate of y, and y's merged edge
    # is the candidate x moved to
    expect_true(keyX[nX] %in% keyY[-length(keyY)])
    expect_identical(keyY[length(keyY)], keyX[j])

    # lMerge from y equals the regraft edge length x split (lReg), and
    # vice versa (the Jacobian's two lengths swap roles between directions)
    lRegX <- x$edge.length[candRow]
    expect_equal(enY$lMerge, lRegX, tolerance = 1e-12)
  }
})

test_that("gibbs_spr_enumerate_cpp rejects invalid prune rows and non-partial paths", {
  ctx <- .CandContext(seed = 11L)
  sp  <- .StateFor(ctx, ctx$tree)
  root <- ctx$nTip + 1L
  expect_error(gibbs_spr_enumerate_cpp(ctx$dataPtr, sp, 0L),
               "out of range")
  expect_error(gibbs_spr_enumerate_cpp(ctx$dataPtr, sp,
                                       nrow(ctx$tree$edge) + 1L),
               "out of range")
  rootRow <- which(ctx$tree$edge[, 1L] == root)[1L]
  expect_error(gibbs_spr_enumerate_cpp(ctx$dataPtr, sp, rootRow),
               "root")

  ctxHet <- .CandContext(seed = 11L,
                         model = MkPrimeModel(qHeterogeneity = TRUE,
                                              nBetaCat = 4L, nCat = 1L,
                                              coding = "none"))
  spHet <- .StateFor(ctxHet, ctxHet$tree)
  eligibleHet <- which(ctxHet$tree$edge[, 1L] != root)[1L]
  expect_error(gibbs_spr_enumerate_cpp(ctxHet$dataPtr, spHet, eligibleHet),
               "partial-CL")
})

# ---------------------------------------------------------------------------
# 2. GSPR-002: the commit path recomputes logPrior
# ---------------------------------------------------------------------------

test_that("GSPR-002: accepted gibbs_spr recomputes logPrior, not the stored value", {
  ctx    <- .CandContext(seed = 13L)
  state0 <- MkPrime:::.InitState(ctx$tree, ctx$mkd, ctx$model)
  correct <- state0$log_prior
  # Poison the stored prior: only a commit path that RECOMPUTES the prior
  # can restore the correct value (the branch prior is flat Dirichlet, so
  # the true prior is invariant across gibbs_spr moves).
  state0$log_prior <- correct - 123.45
  sp <- MkPrime:::.InitMcmcChain(state0)
  fill_partition_cache(ctx$dataPtr, sp)
  allocate_cl_workspace(ctx$dataPtr, sp)

  set.seed(2)
  accepted <- FALSE
  for (i in seq_len(200L)) {
    if (do_move_cpp(ctx$dataPtr, sp, 10L, 0L, 0.5, 0.5, 1L, 1.0)) {
      accepted <- TRUE
      break
    }
  }
  expect_true(accepted)
  expect_equal(get_mcmc_state(sp)$logPrior, correct, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# 3. Kernel behaviour: real MH step, live merged-edge moves, no pi-null mass
# ---------------------------------------------------------------------------

test_that("gibbs_spr has a real accept/reject and both move flavours occur", {
  ctx <- .CandContext(seed = 17L)
  sp  <- .StateFor(ctx, ctx$tree)

  set.seed(3)
  nMove <- 300L
  nAcc <- 0L
  fracOnly <- 0L    # accepted merged-edge draws: fractions move, topology not
  topoMove <- 0L    # accepted regrafts onto another edge
  for (i in seq_len(nMove)) {
    before <- get_mcmc_state(sp)
    # get_mcmc_state returns relBrLengths as a view of the live state; force
    # a copy so the pre-move snapshot survives the move
    beforeRel <- c(before$relBrLengths)
    acc <- do_move_cpp(ctx$dataPtr, sp, 10L, 0L, 0.5, 0.5, 1L, 0.0)
    after <- get_mcmc_state(sp)
    if (acc) {
      nAcc <- nAcc + 1L
      if (identical(before$edge, after$edge)) {
        fracOnly <- fracOnly + 1L
        expect_false(identical(beforeRel, c(after$relBrLengths)))
      } else {
        topoMove <- topoMove + 1L
      }
      # GSPR-001 fingerprint: the deterministic 0.5 * lReg split put
      # bit-equal edge pairs in every accepted state; a drawn tau never does
      expect_false(anyDuplicated(after$relBrLengths) > 0L)
    } else {
      expect_identical(before$edge, after$edge)
      expect_identical(beforeRel, c(after$relBrLengths))
    }
  }
  expect_gt(nAcc, 0L)
  expect_lt(nAcc, nMove)   # a real MH step rejects some proposals
  expect_gt(fracOnly, 0L)  # merged-edge draws are live moves, not no-ops
  expect_gt(topoMove, 0L)
})

# ---------------------------------------------------------------------------
# 4. The two non-default evaluation paths commit a canonical state
# ---------------------------------------------------------------------------

.ExpectCoherentCommits <- function(ctx, nMove = 150L, seed = 5L) {
  sp <- .StateFor(ctx, ctx$tree)
  set.seed(seed)
  nAcc <- 0L
  for (i in seq_len(nMove)) {
    if (do_move_cpp(ctx$dataPtr, sp, 10L, 0L, 0.5, 0.5, 1L, 1.0)) {
      nAcc <- nAcc + 1L
      s <- get_mcmc_state(sp)
      expect_equal(sum(s$relBrLengths), 1.0, tolerance = 1e-9)
      # Committed logLik is the canonical full evaluation at the drawn tau
      expect_equal(s$logLik, eval_full_loglik_cpp(ctx$dataPtr, sp),
                   tolerance = 1e-9)
    }
  }
  expect_gt(nAcc, 0L)
  expect_lt(nAcc, nMove)
}

test_that("gibbs_spr full-fallback path (coding = informative) is coherent", {
  .ExpectCoherentCommits(.CandContext(seed = 19L,
                                      model = MkPrimeModel(coding = "informative")))
})

test_that("gibbs_spr Q-heterogeneity path is coherent", {
  .ExpectCoherentCommits(.CandContext(seed = 23L,
                                      model = MkPrimeModel(qHeterogeneity = TRUE,
                                                           nBetaCat = 4L,
                                                           nCat = 1L,
                                                           coding = "none")))
})

test_that("gibbs_spr Q-heterogeneity path with ascertainment is coherent", {
  .ExpectCoherentCommits(.CandContext(seed = 29L,
                                      model = MkPrimeModel(qHeterogeneity = TRUE,
                                                           nBetaCat = 4L,
                                                           nCat = 1L,
                                                           coding = "variable")))
})
