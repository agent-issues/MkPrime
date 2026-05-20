# Tests for .BuildMovesPartitioned() and .ParamNamesPartitioned() — the
# partition-aware move spec and sample-matrix column naming added in
# Layer 1 plumbing. Trivial spec collapses to the legacy spec/columns;
# non-trivial spec appends per-class moves and column names.

library("TreeTools")


.dummy_mkd <- function(nChar = 8L, nTip = 6L, seed = 17L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  MkPrimeData(MatrixToPhyDat(mat))
}


# ---- trivial spec: behaves exactly as the legacy functions ----

test_that(".BuildMovesPartitioned with NULL partition is identical to .BuildMoves", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = NULL, unlink = character(0), nClasses = 1L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  expect_identical(part, legacy)
})


test_that(".BuildMovesPartitioned with nClasses == 1 is identical to .BuildMoves", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = rep(1L, 8L), unlink = character(0), nClasses = 1L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  expect_identical(part, legacy)
})


test_that(".ParamNamesPartitioned with NULL partition is identical to .ParamNames", {
  mkd  <- .dummy_mkd()
  spec <- list(partition = NULL, unlink = character(0), nClasses = 1L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec)
  expect_identical(part, legacy)
})


# ---- non-trivial spec: per-class moves and columns appended ----

test_that("'shape' unlink emits one scale move per class", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  # Should be legacy + 2 new per-class moves
  expect_identical(length(part), length(legacy) + 2L)
  newNames <- setdiff(vapply(part, `[[`, character(1), "name"),
                      vapply(legacy, `[[`, character(1), "name"))
  expect_identical(newNames, c("scale_class_rate_log_sd_1",
                               "scale_class_rate_log_sd_2"))
  # All new moves carry the classIdx field
  newMoves <- part[(length(legacy) + 1L):length(part)]
  expect_identical(vapply(newMoves, `[[`, integer(1), "classIdx"), 1:2)
  expect_identical(unique(vapply(newMoves, `[[`, character(1), "type")),
                   "scale_class_rate_log_sd")
})


test_that("'ratemultiplier' unlink emits a single Dirichlet-simplex move on class_w", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L),
               unlink    = "ratemultiplier",
               nClasses  = 3L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  expect_identical(length(part), length(legacy) + 1L)
  newMove <- part[[length(part)]]
  expect_identical(newMove$name, "dirichlet_simplex_class_w")
  expect_identical(newMove$type, "dirichlet_simplex_class_w")
  expect_identical(newMove$target, "class_w")
  expect_identical(newMove$dim, 3L)
})


test_that("both 'shape' and 'ratemultiplier' unlinked emits both move kinds", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = c("shape", "ratemultiplier"),
               nClasses  = 2L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  # Legacy + 2 shape moves + 1 ratemultiplier move = legacy + 3
  expect_identical(length(part), length(legacy) + 3L)
  newNames <- vapply(part, `[[`, character(1), "name")[
    (length(legacy) + 1L):length(part)]
  expect_identical(newNames, c("scale_class_rate_log_sd_1",
                               "scale_class_rate_log_sd_2",
                               "dirichlet_simplex_class_w"))
})


# ---- ParamNames extensions ----

test_that(".ParamNamesPartitioned appends class<c>_rate_log_sd columns when 'shape' unlinked", {
  mkd  <- .dummy_mkd(nChar = 8L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec)
  appended <- setdiff(part, legacy)
  expect_identical(appended, c("class1_rate_log_sd", "class2_rate_log_sd"))
  # Legacy columns appear first, in unchanged order
  expect_identical(part[seq_along(legacy)], legacy)
})


test_that(".ParamNamesPartitioned appends w_<c> columns when 'ratemultiplier' unlinked", {
  mkd  <- .dummy_mkd(nChar = 8L)
  spec <- list(partition = c(1L, 1L, 1L, 2L, 2L, 3L, 3L, 3L),
               unlink    = "ratemultiplier",
               nClasses  = 3L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec)
  appended <- setdiff(part, legacy)
  expect_identical(appended, c("w_1", "w_2", "w_3"))
})


test_that(".ParamNamesPartitioned appends both column families when both unlinked", {
  mkd  <- .dummy_mkd(nChar = 8L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = c("shape", "ratemultiplier"),
               nClasses  = 2L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec)
  appended <- setdiff(part, legacy)
  # shape columns appended first, then ratemultiplier columns (per the
  # function's order of token handling)
  expect_identical(appended, c("class1_rate_log_sd", "class2_rate_log_sd",
                               "w_1", "w_2"))
})
