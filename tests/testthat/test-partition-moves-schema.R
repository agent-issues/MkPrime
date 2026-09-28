# Tests for .BuildMovesPartitioned() and .ParamNamesPartitioned() — the
# partition-aware move spec and sample-matrix column naming added in
# Layer 1 plumbing. Trivial spec collapses to the legacy spec/columns;
# non-trivial spec appends per-class moves and column names.



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

test_that("'shape' unlink + hyperprior emits per-class scale moves + scale_hyper_tau", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec)
  # When "shape" is unlinked the partitioned move list drops the legacy
  # scalar rate_log_sd / slice_rate_log_sd / joint_tl_rls moves (they would
  # break the lockstep with classRateLogSd[0]) and appends per-class moves
  # plus a scale_hyper_tau move (pooled hyperprior is the default).
  partNames   <- vapply(part,   `[[`, character(1), "name")
  legacyNames <- vapply(legacy, `[[`, character(1), "name")
  dropped <- setdiff(legacyNames, partNames)
  expect_setequal(dropped,
                  c("rate_log_sd", "slice_rate_log_sd", "joint_tl_rls"))
  added <- setdiff(partNames, legacyNames)
  expect_identical(added, c("scale_class_rate_log_sd_1",
                            "scale_class_rate_log_sd_2",
                            "scale_hyper_tau"))
  # Per-class moves carry the classIdx field; scale_hyper_tau does not.
  perClass <- part[grep("^scale_class_rate_log_sd_", partNames)]
  expect_identical(vapply(perClass, `[[`, integer(1), "classIdx"), 1:2)
  expect_identical(unique(vapply(perClass, `[[`, character(1), "type")),
                   "scale_class_rate_log_sd")
  tauMove <- part[[which(partNames == "scale_hyper_tau")]]
  expect_identical(tauMove$type, "scale_hyper_tau")
  expect_identical(tauMove$target, "hyper_tau")
  expect_null(tauMove$classIdx)
})


test_that("'shape' unlink with gamma_independent prior emits no scale_hyper_tau", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .BuildMoves(nEdge = 9L, nTrans = 8L, hasNeo = FALSE, mcmc = mcmc)
  part   <- .BuildMovesPartitioned(nEdge = 9L, nTrans = 8L, hasNeo = FALSE,
                                   mcmc = mcmc, partitionSpec = spec,
                                   priorOnClassRateLogSd = "gamma_independent")
  partNames   <- vapply(part,   `[[`, character(1), "name")
  legacyNames <- vapply(legacy, `[[`, character(1), "name")
  added <- setdiff(partNames, legacyNames)
  expect_identical(added, c("scale_class_rate_log_sd_1",
                            "scale_class_rate_log_sd_2"))
  expect_false("scale_hyper_tau" %in% partNames)
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
  partNames   <- vapply(part,   `[[`, character(1), "name")
  legacyNames <- vapply(legacy, `[[`, character(1), "name")
  added <- setdiff(partNames, legacyNames)
  expect_identical(added, c("scale_class_rate_log_sd_1",
                            "scale_class_rate_log_sd_2",
                            "scale_hyper_tau",
                            "dirichlet_simplex_class_w"))
})


# ---- ParamNames extensions ----

test_that(".ParamNamesPartitioned appends shape + hyperprior columns by default", {
  mkd  <- .dummy_mkd(nChar = 8L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec)
  appended <- setdiff(part, legacy)
  expect_identical(appended,
                   c("class1_rate_log_sd", "class2_rate_log_sd",
                     "hyper_tau",
                     "class1_rate_log_sd_z", "class2_rate_log_sd_z"))
  # Legacy columns appear first, in unchanged order
  expect_identical(part[seq_along(legacy)], legacy)
})


test_that(".ParamNamesPartitioned skips hyperprior columns under gamma_independent", {
  mkd  <- .dummy_mkd(nChar = 8L)
  spec <- list(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
               unlink    = "shape",
               nClasses  = 2L)
  legacy <- .ParamNames(mkd, nEdge = 9L)
  part   <- .ParamNamesPartitioned(mkd, nEdge = 9L, partitionSpec = spec,
                                   priorOnClassRateLogSd = "gamma_independent")
  appended <- setdiff(part, legacy)
  expect_identical(appended, c("class1_rate_log_sd", "class2_rate_log_sd"))
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
  # shape columns appended first, then ratemultiplier columns, then
  # hyperprior columns (hyper_tau + per-class z).
  expect_identical(appended, c("class1_rate_log_sd", "class2_rate_log_sd",
                               "w_1", "w_2",
                               "hyper_tau",
                               "class1_rate_log_sd_z",
                               "class2_rate_log_sd_z"))
})


# ---- #77: per-class moves share the scalar floor ----

test_that("per-class scalar moves receive the scalar weight floor", {
  mcmc <- MkPrimeMCMC(nIter = 100L, autoTune = FALSE, nRuns = 1L,
                      minWarmup = 50L, maxWarmup = 50L)
  spec <- list(partition = rep(1:2, length.out = 10L), unlink = "shape",
               nClasses = 2L)
  # A large tree, so the floor exceeds the per-class moves' raw weight of 1.
  moves <- MkPrime:::.BuildMovesPartitioned(200L, 5L, TRUE, mcmc, spec)
  w <- vapply(moves, `[[`, numeric(1), "weight")
  names(w) <- vapply(moves, `[[`, character(1), "name")
  perClass <- c("scale_class_rate_log_sd_1", "scale_class_rate_log_sd_2",
                "scale_hyper_tau")
  expect_gt(w[["tree_length"]], 1)
  # tree_length's raw weight is also 1: all four sit at the same floor.
  expect_equal(unname(w[perClass]), rep(w[["tree_length"]], 3))
  expect_true(all(perClass %in% MkPrime:::.ScalarFloorMoves(moves)))
})
