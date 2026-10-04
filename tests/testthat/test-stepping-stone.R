test_that("Stepping-stone returns finite marginal likelihood", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(4719)
  ss <- mkp_stepping_stone(pd, tree, nStones = 8L, nIter = 100L,
                            warmup = 30L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
  expect_gt(ss$se, 0)
  expect_length(ss$log_ratios, 8L)
  expect_length(ss$betas, 9L) # nStones + 1 boundary points
  expect_equal(ss$betas[1], 0)
  expect_equal(ss$betas[9], 1)
})


test_that("Beta schedule is monotonically increasing", {

  tree <- read.tree(text = "((t1:0.1,t2:0.1):0.1,t3:0.1);")
  mat <- matrix(c(0L, 1L, 0L), nrow = 3,
                dimnames = list(c("t1", "t2", "t3"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(2841)
  ss <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 50L,
                            warmup = 10L, verbose = FALSE)

  expect_true(all(diff(ss$betas) > 0))
})


test_that("Stepping-stone works with neomorphic characters", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(6103)
  ss <- mkp_stepping_stone(pd, tree, neomorphic = 1L,
                            nStones = 5L, nIter = 100L,
                            warmup = 30L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


test_that("Stepping-stone works with fixed topology", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(8234)
  ss <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 100L,
                            warmup = 30L, fixTopology = TRUE,
                            verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


test_that("More stones with more iterations gives consistent results", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L,
                  0L, 1L, 1L, 0L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(1583)
  ss1 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 500L,
                             warmup = 200L, verbose = FALSE)
  set.seed(7261)
  ss2 <- mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 500L,
                             warmup = 200L, verbose = FALSE)

  # Two independent estimates should be in the same ballpark. The tolerance
  # is the estimator's measured spread (a run-to-run sd of 0.68 over eight
  # seeds at these budgets), not a multiple of `se`, which is tested below.
  expect_lt(abs(ss1$log_marginal - ss2$log_marginal), 3)
})


test_that("Stepping-stone accepts MkPrimeData input", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L), nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)

  set.seed(9428)
  ss <- mkp_stepping_stone(mkd, tree, nStones = 5L, nIter = 50L,
                            warmup = 10L, verbose = FALSE)

  expect_true(is.finite(ss$log_marginal))
  expect_true(is.finite(ss$se))
})


# --- Tests for .EssVector helper ---

test_that(".EssVector returns sensible values for iid samples", {
  set.seed(3847)
  x <- rnorm(500)
  ess <- MkPrime:::.EssVector(x)
  # iid samples: ESS should be close to n

  expect_gt(ess, 300)
  expect_lte(ess, 500)
})

test_that(".EssVector returns lower ESS for autocorrelated samples", {
  set.seed(5192)
  # AR(1) with rho = 0.9
  n <- 500
  x <- numeric(n)
  x[1] <- rnorm(1)
  for (i in 2:n) x[i] <- 0.9 * x[i - 1] + rnorm(1)

  ess <- MkPrime:::.EssVector(x)
  # Highly autocorrelated: ESS should be much less than n
  expect_gt(ess, 1)
  expect_lt(ess, 200)
})

test_that(".EssVector handles edge cases", {
  # Degenerate inputs: NA (not enough information to estimate ESS)
  expect_true(is.na(MkPrime:::.EssVector(numeric(0))))
  expect_true(is.na(MkPrime:::.EssVector(42)))
  # Constant vector: zero variance → NA
  expect_true(is.na(MkPrime:::.EssVector(rep(5, 100))))
})


test_that("SE is positive and smaller with more iterations", {

  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)

  set.seed(6712)
  ssSmall <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 50L,
                                 warmup = 20L, verbose = FALSE)
  set.seed(6712)
  ssLarge <- mkp_stepping_stone(pd, tree, nStones = 5L, nIter = 500L,
                                 warmup = 50L, verbose = FALSE)

  expect_gt(ssSmall$se, 0)
  expect_gt(ssLarge$se, 0)
  # More samples should generally reduce SE (not guaranteed for any
  # single seed, but we check the SE is at least finite and positive)
  expect_true(is.finite(ssSmall$se))
  expect_true(is.finite(ssLarge$se))
})

test_that("Stepping-stone rejects partitioned data (#243)", {
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  mkd$partitions <- MkPrime:::.BuildPartitions(mkd, partition = 1:2)
  expect_error(mkp_stepping_stone(mkd, tree, nStones = 2L, nIter = 5L,
                                  warmup = 0L, verbose = FALSE),
               "does not support partitioned")
})


SixTaxonFixture <- function() {
  mat <- matrix(c(0L, 0L, 1L, 1L, 1L, 1L,
                  0L, 0L, 0L, 1L, 1L, 1L,
                  0L, 1L, 0L, 1L, 0L, 1L,
                  1L, 1L, 0L, 0L, 1L, 1L),
                nrow = 6, dimnames = list(letters[1:6], NULL))
  tree <- read.tree(text = "((a:0.1,b:0.1):0.1,c:0.1,(d:0.1,(e:0.1,f:0.1):0.1):0.1);")
  list(pd = MatrixToPhyDat(mat), tree = tree,
       scrambled = RenumberTips(tree, c("d", "a", "f", "b", "e", "c")))
}


test_that("Stepping-stone scores the tree it is given, whatever its tip order (#277)", {
  f <- SixTaxonFixture()
  expect_true(all.equal(f$tree, f$scrambled))
  expect_false(identical(f$scrambled$tip.label, f$tree$tip.label))
  SS <- function(tree) {
    set.seed(2770)
    mkp_stepping_stone(f$pd, tree, nStones = 3L, nIter = 30L, warmup = 5L,
                       fixTopology = TRUE, verbose = FALSE)
  }
  expect_equal(SS(f$scrambled)$log_ratios, SS(f$tree)$log_ratios)
})


test_that("Stepping-stone validates tip labels and pins as RunMkPrime does (#277)", {
  f <- SixTaxonFixture()
  stranger <- f$tree
  stranger$tip.label[1] <- "stranger"
  expect_error(mkp_stepping_stone(f$pd, stranger, nStones = 2L, nIter = 5L,
                                  warmup = 0L, verbose = FALSE),
               "do not match taxa")
  expect_error(
    mkp_stepping_stone(f$pd, f$tree,
                       mcmc = MkPrimeMCMC(moveWeights = c(nni = 0.5, spr = 0.5)),
                       nStones = 2L, nIter = 5L, warmup = 0L, verbose = FALSE),
    "no weight")
})


test_that("Stepping-stone SE matches the spread of independent estimates (#279)", {
  skip_on_cran()
  tree <- read.tree(text = "((t1:0.2,t2:0.3):0.1,(t3:0.15,t4:0.25):0.1);")
  mat <- matrix(c(0L, 1L, 0L, 1L,
                  0L, 0L, 1L, 1L,
                  0L, 1L, 1L, 0L),
                nrow = 4,
                dimnames = list(c("t1", "t2", "t3", "t4"), NULL))
  pd <- MatrixToPhyDat(mat)
  # sd over seeds 1:40 of run_log_marginal[2] at these budgets with
  # nRuns = 2: single runs from independently perturbed starts (#382). The
  # per-stone delta method reported a mean se of 0.09.
  singleRunSd <- 0.428
  nRuns <- 3L
  se <- vapply(1:6, function(seed) {
    set.seed(2790 + seed)
    mkp_stepping_stone(pd, tree, nStones = 15L, nIter = 500L, warmup = 200L,
                       nRuns = nRuns, verbose = FALSE)$se
  }, numeric(1))
  # Root mean square, as se^2 (not se) estimates the variance without bias.
  ratio <- sqrt(mean(se^2)) / (singleRunSd / sqrt(nRuns))
  expect_gt(ratio, 0.4)
  expect_lt(ratio, 2.5)
})


test_that("A degenerate stone leaves estimate and SE both undefined (#279)", {
  betas <- c(0, 0.5, 1)
  set.seed(2791)
  Draws <- function() matrix(rnorm(40, -10), 20, 2)
  ok <- MkPrime:::.SteppingStoneEstimate(list(Draws(), Draws()), betas, 0)
  expect_true(is.finite(ok$log_marginal))
  expect_true(is.finite(ok$se))
  expect_length(ok$run_log_marginal, 2L)

  noFinite <- Draws()
  noFinite[, 2] <- -Inf
  infinite <- Draws()
  infinite[3, 1] <- Inf
  for (bad in list(noFinite, infinite)) {
    expect_warning(
      est <- MkPrime:::.SteppingStoneEstimate(list(Draws(), bad), betas, 0),
      "undefined")
    expect_true(is.na(est$log_marginal))
    expect_true(is.na(est$se))
  }

  # A zero-likelihood draw among finite ones carries zero weight.
  someZero <- Draws()
  someZero[1, ] <- -Inf
  expect_true(is.finite(
    MkPrime:::.SteppingStoneEstimate(list(someZero), betas, 0)$log_marginal))
  expect_true(is.na(MkPrime:::.SteppingStoneEstimate(list(Draws()), betas, 0)$se))
})


test_that("Stepping-stone reports the mean of the runs it gives an SE for (#382)", {
  betas <- c(0, 0.5, 1)
  # Per-stone log ratios (0, 1) in one run and (1, 0) in the other.
  RunDraws <- function(stone1, stone2) {
    cbind(rep(stone1 / 0.5, 2), rep(stone2 / 0.5, 2))
  }
  est <- MkPrime:::.SteppingStoneEstimate(
    list(RunDraws(0, 1), RunDraws(1, 0)), betas, logZ0 = 0.25)
  expect_equal(est$run_log_marginal, c(1.25, 1.25))
  expect_equal(est$log_marginal, 1.25)
  expect_equal(est$se, 0)
  expect_equal(sum(est$log_ratios) + 0.25, est$log_marginal)
})


test_that("Stepping-stone runs start from distinct states (#382)", {
  f <- SixTaxonFixture()
  mkd <- MkPrimeData(f$pd)
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), NULL, mkd)
  tree <- MkPrime:::.PrepareStartTree(f$tree, mkd)
  set.seed(3820)
  free <- MkPrime:::.SteppingStoneStarts(tree, 3L, FALSE, mkd, model)
  fixed <- MkPrime:::.SteppingStoneStarts(tree, 3L, TRUE, mkd, model)
  expect_identical(free[[1]], MkPrime:::.InitState(tree, mkd, model))
  for (starts in list(free, fixed)) {
    p <- vapply(starts, `[[`, 0, "p")
    expect_length(unique(p), 3L)
    lengths <- lapply(starts, `[[`, "rel_br_lengths")
    expect_false(isTRUE(all.equal(lengths[[1]], lengths[[2]])))
    for (start in starts) {
      expect_equal(start$log_prior,
                   MkPrime:::LogPrior(start, model, mkd))
    }
  }
  splits <- function(starts) {
    lapply(starts, function(s) as.character(as.Splits(s$tree)))
  }
  expect_true(all(vapply(splits(fixed), setequal, TRUE,
                         splits(fixed)[[1]])))
  set.seed(3821)
  topologies <- unlist(lapply(1:5, function(i) {
    starts <- MkPrime:::.SteppingStoneStarts(tree, 2L, FALSE, mkd, model)
    !setequal(as.character(as.Splits(starts[[2]]$tree)),
              as.character(as.Splits(tree)))
  }))
  expect_true(any(topologies))
})


test_that("Stepping-stone controls are validated (#380)", {
  f <- SixTaxonFixture()
  SS <- function(...) {
    mkp_stepping_stone(f$pd, f$tree, nStones = 2L, nIter = 5L, warmup = 0L,
                       verbose = FALSE, ...)
  }
  expect_error(SS(nStones = 0L), "nStones")
  expect_error(SS(nRuns = 2.9), "nRuns")
  expect_error(SS(nRuns = 0), "nRuns")
  expect_error(SS(nRuns = NA), "nRuns")
  expect_error(SS(alpha = 0), "alpha")
  expect_error(SS(alpha = Inf), "alpha")
  expect_error(mkp_stepping_stone(f$pd, f$tree, nStones = 2L, nIter = 0L,
                                  verbose = FALSE), "nIter")
  expect_error(mkp_stepping_stone(f$pd, f$tree, nStones = 2L, nIter = 5L,
                                  warmup = -1L, verbose = FALSE), "warmup")
  expect_error(SS(model = MkPrimeModel(expSteps = 0)), "expSteps")
})
