# Tests for M-125: Block Dirichlet branch-length proposal
# Tests proposal correctness, .kMoveTypes mapping, move construction,
# and C++ engine integration.

test_that("ProposeDirichletSimplex preserves simplex sum", {
  set.seed(6417)
  x <- c(0.3, 0.2, 0.15, 0.35)
  for (i in 1:30) {
    result <- ProposeDirichletSimplex(x, nCats = 2, alpha = 10)
    expect_equal(sum(result$value), sum(x), tolerance = 1e-10)
    expect_true(all(result$value > 0))
    expect_true(is.finite(result$logHastings))
  }
})

test_that("ProposeDirichletSimplex works with nCats = n (full Dirichlet)", {
  set.seed(8194)
  x <- c(0.25, 0.25, 0.25, 0.25)
  for (i in 1:20) {
    result <- ProposeDirichletSimplex(x, nCats = 4, alpha = 10)
    expect_equal(sum(result$value), sum(x), tolerance = 1e-10)
    expect_true(all(result$value > 0))
    expect_true(is.finite(result$logHastings))
  }
})

test_that("ProposeDirichletSimplex handles large simplex", {
  set.seed(2741)
  n <- 44  # typical nEdge for 23-taxon tree
  x <- rep(1, n)  # uniform simplex (sums to n)
  for (i in 1:20) {
    result <- ProposeDirichletSimplex(x, nCats = 10, alpha = 10)
    expect_equal(sum(result$value), sum(x), tolerance = 1e-8)
    expect_true(all(result$value > 0))
    expect_true(is.finite(result$logHastings))
  }
})

test_that("ProposeDirichletSimplex with higher alpha is more conservative", {
  set.seed(3958)
  x <- c(0.3, 0.2, 0.15, 0.35)
  diffs_low <- numeric(100)
  diffs_high <- numeric(100)
  for (i in 1:100) {
    r1 <- ProposeDirichletSimplex(x, nCats = 3, alpha = 5)
    r2 <- ProposeDirichletSimplex(x, nCats = 3, alpha = 100)
    diffs_low[i] <- max(abs(r1$value - x))
    diffs_high[i] <- max(abs(r2$value - x))
  }
  # Higher alpha should produce smaller perturbations on average
  expect_gt(mean(diffs_low), mean(diffs_high))
})

test_that(".kMoveTypes maps dirichlet_branch to 23", {
  mt <- MkPrime:::.kMoveTypes
  expect_identical(mt[["dirichlet_branch"]], 23L)
})

test_that(".BuildMoves includes dirichlet_branch when enabled", {
  mc <- MkPrimeMCMC(dirichletBranch = TRUE)
  moves <- MkPrime:::.BuildMoves(12L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_true("dirichlet_branch" %in% names)
  db <- moves[[which(names == "dirichlet_branch")]]
  expect_equal(db$nCats, 5L)  # min(12, 5) — M-127: K=5 default
  expect_equal(db$dim, 5L)
  expect_equal(db$target, "rel_br_lengths")
})

test_that(".BuildMoves nCats capped at nEdge for small trees", {
  mc <- MkPrimeMCMC(dirichletBranch = TRUE)
  moves <- MkPrime:::.BuildMoves(6L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_true("dirichlet_branch" %in% names)
  db <- moves[[which(names == "dirichlet_branch")]]
  expect_equal(db$nCats, 5L)  # min(6, 5) = 5 — M-127: K=5 default
})

test_that(".BuildMoves excludes dirichlet_branch for very small trees", {
  mc <- MkPrimeMCMC(dirichletBranch = TRUE)
  # nEdge = 3 → too small (< 4)
  moves <- MkPrime:::.BuildMoves(3L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_false("dirichlet_branch" %in% names)
})

test_that("dirichletK overrides default K in .BuildMoves", {
  mc <- MkPrimeMCMC(dirichletBranch = TRUE, dirichletK = 8L)
  moves <- MkPrime:::.BuildMoves(20L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  db <- moves[[which(names == "dirichlet_branch")]]
  expect_equal(db$nCats, 8L)
  expect_equal(db$dim, 8L)
})

test_that(".BuildMoves excludes dirichlet_branch when disabled", {
  mc <- MkPrimeMCMC(dirichletBranch = FALSE)
  moves <- MkPrime:::.BuildMoves(12L, 0L, FALSE, mc, FALSE, "geometric")
  names <- vapply(moves, `[[`, character(1), "name")
  expect_false("dirichlet_branch" %in% names)
})

test_that("dirichlet_alpha tuning default is 0.1", {
  mc <- MkPrimeMCMC()
  expect_equal(mc$tuning$dirichlet_alpha, 0.1)
})

test_that(".AdaptTuning adjusts dirichlet_alpha", {
  moves <- list(
    list(name = "dirichlet_branch", type = "dirichlet_simplex",
         target = "rel_br_lengths", weight = 3, dim = 10L)
  )
  tuning <- list(dirichlet_alpha = 10)

  # High acceptance → should decrease alpha (bolder proposals to lower acceptance)
  accept <- c(dirichlet_branch = 50L)
  propose <- c(dirichlet_branch = 100L)
  new_tuning <- MkPrime:::.AdaptTuning(tuning, accept, propose, moves)
  expect_lt(new_tuning$dirichlet_alpha, 10)

  # Low acceptance → should increase alpha (more conservative to raise acceptance)
  accept2 <- c(dirichlet_branch = 5L)
  propose2 <- c(dirichlet_branch = 100L)
  new_tuning2 <- MkPrime:::.AdaptTuning(tuning, accept2, propose2, moves)
  expect_gt(new_tuning2$dirichlet_alpha, 10)
})

test_that("dirichlet_branch appears in moveWeights validNames", {
  # Should not error when specifying dirichlet_branch weight
  mc <- MkPrimeMCMC(moveWeights = c(dirichlet_branch = 0.1))
  expect_equal(mc$moveWeights[["dirichlet_branch"]], 0.1)
})
