# Tests for M-127: partial CL evaluation for Dirichlet branch proposals.

test_that("dirichlet_simplex_proposal returns modifiedEdges", {
  set.seed(3917)
  x <- rep(1/10, 10)
  result <- dirichlet_simplex_proposal(x, 3L, 1.0)
  expect_true("modifiedEdges" %in% names(result))
  # 0-indexed edge indices, length = nCats
  expect_length(result$modifiedEdges, 3L)
  expect_true(all(result$modifiedEdges >= 0L))
  expect_true(all(result$modifiedEdges < 10L))
  # All indices are distinct

  expect_equal(length(unique(result$modifiedEdges)), 3L)
})

test_that("modifiedEdges identifies which elements changed", {
  set.seed(6204)
  x <- c(0.1, 0.2, 0.3, 0.15, 0.25)
  result <- dirichlet_simplex_proposal(x, 2L, 5.0)
  mod <- result$modifiedEdges + 1L  # convert to 1-indexed
  # Modified elements should have changed
  changed <- which(result$value != x)
  expect_true(all(mod %in% changed))
})

# ---------------------------------------------------------------------------
# Local Dirichlet (select_neighborhood + BFS)
# ---------------------------------------------------------------------------

test_that("local_dirichlet_proposal returns connected edges", {
  set.seed(4419)
  tree <- ape::rtree(15)
  x <- rep(1 / nrow(tree$edge), nrow(tree$edge))
  result <- local_dirichlet_proposal(x, tree$edge[, 1], tree$edge[, 2],
                                       5L, 2.0)
  expect_true("modifiedEdges" %in% names(result))
  expect_length(result$modifiedEdges, 5L)
  # All distinct
  expect_equal(length(unique(result$modifiedEdges)), 5L)

  # Check connectivity: build edge adjacency, verify BFS path exists
  mod <- result$modifiedEdges + 1L
  nodes <- unique(c(tree$edge[mod, 1], tree$edge[mod, 2]))
  # For K connected edges, should involve at most K+1 nodes
  # (a path of K edges has K+1 nodes, a star has 1 + K nodes)
  expect_lte(length(nodes), 6L)  # K+1
})

test_that("local_dirichlet_proposal preserves simplex sum", {
  set.seed(6821)
  tree <- ape::rtree(10)
  x <- runif(nrow(tree$edge))
  x <- x / sum(x)
  result <- local_dirichlet_proposal(x, tree$edge[, 1], tree$edge[, 2],
                                       3L, 5.0)
  expect_equal(sum(result$value), 1.0, tolerance = 1e-12)
})

test_that("local_dirichlet_proposal only modifies selected edges", {
  set.seed(1023)
  tree <- ape::rtree(10)
  x <- runif(nrow(tree$edge))
  x <- x / sum(x)
  result <- local_dirichlet_proposal(x, tree$edge[, 1], tree$edge[, 2],
                                       3L, 5.0)
  mod <- result$modifiedEdges + 1L
  unchanged <- setdiff(seq_along(x), mod)
  expect_equal(result$value[unchanged], x[unchanged])
})

test_that(".BuildMoves includes local_dirichlet when enabled", {
  mc <- MkPrimeMCMC(localDirichlet = TRUE)
  moves <- MkPrime:::.BuildMoves(20L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_true("local_dirichlet" %in% names)
  ld <- moves[[which(names == "local_dirichlet")]]
  expect_equal(ld$nCats, 6L)  # min(20, 6)
  expect_equal(ld$target, "rel_br_lengths")
})

test_that(".BuildMoves excludes local_dirichlet when disabled", {
  mc <- MkPrimeMCMC(localDirichlet = FALSE)
  moves <- MkPrime:::.BuildMoves(20L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_false("local_dirichlet" %in% names)
})

test_that(".kMoveTypes maps local_dirichlet to 24", {
  mt <- MkPrime:::.kMoveTypes
  expect_identical(mt[["local_dirichlet"]], 24L)
})

test_that("localDirichletK overrides default K", {
  mc <- MkPrimeMCMC(localDirichlet = TRUE, localDirichletK = 8L)
  moves <- MkPrime:::.BuildMoves(20L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  ld <- moves[[which(names == "local_dirichlet")]]
  expect_equal(ld$nCats, 8L)
})

test_that("dirichletK parameter validates correctly", {
  expect_error(MkPrimeMCMC(dirichletK = 1L), "at least 2")
  expect_error(MkPrimeMCMC(dirichletK = 0L), "at least 2")
  mc <- MkPrimeMCMC(dirichletK = 3L)
  expect_equal(mc$dirichletK, 3L)
  mc_null <- MkPrimeMCMC()
  expect_null(mc_null$dirichletK)
})

test_that("partial eval Dirichlet gives consistent results in MCMC", {
  skip_if_not_installed("TreeSearch")
  dat <- TreeSearch::inapplicable.phyData[["Vinther2008"]]
  mkd <- suppressWarnings(MkPrimeData(dat))
  model <- MkPrimeModel("variable")
  tipLabels <- names(dat)
  set.seed(7723)
  tree <- ape::rtree(length(tipLabels), tip.label = tipLabels, br = NULL)
  tree$edge.length <- rep(1, nrow(tree$edge))

  # Short run with Dirichlet enabled — should produce valid samples
  mcmc <- MkPrimeMCMC(nIter = 500L, thin = 1L, maxWarmup = 100L,
                       minWarmup = 100L, dirichletBranch = TRUE,
                       dirichletK = 3L)
  post <- suppressWarnings(RunMkPrime(mkd, tree = tree, model = model, mcmc = mcmc,
                                       fixTopology = TRUE))
  expect_false(any(is.na(post$samples)))
  expect_true(all(is.finite(post$samples[, "log_likelihood"])))
  # Likelihood should not be stuck at one value (would indicate
  # partial eval is always rejecting or accepting incorrectly)
  ll <- post$samples[, "log_likelihood"]
  expect_gt(length(unique(ll)), 5L)
})
