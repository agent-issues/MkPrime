# A user's moveWeights pins share the schedule with the auto-pinned
# always-accepting moves; the free moves must never be left with no weight,
# and C++ must never receive a weight vector that is not a distribution.

test_that("pins that leave the free moves no budget abort the run", {
  set.seed(7101)
  tree <- BalancedTree(8)
  tree$edge.length <- rep(0.1, nrow(tree$edge))
  mat <- matrix(sample(0:2, 8 * 6, replace = TRUE), nrow = 8,
                dimnames = list(tree$tip.label, NULL))
  # Valid alone (it sums to at most 1), but gibbs_kPrime and
  # slice_rate_log_sd are auto-pinned on top of it.
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L, maxWarmup = 100L,
                      thin = 10L, nRuns = 1L, maxTime = 30, autoTune = FALSE,
                      moveWeights = c(nni = 0.95))
  expect_error(
    RunMkPrime(data = MkPrimeData(MatrixToPhyDat(mat)), tree = tree,
               model = MkPrimeModel(), mcmc = mcmc),
    "no weight"
  )
})

test_that(".SchedulePins names the threshold the pins must stay below", {
  moveWeights <- c(nni = 0.3, spr = 0.2, gibbs_kPrime = 0.3,
                   slice_rate_log_sd = 0.1, tree_length = 0.1)
  moveTypes <- c("nni", "spr", "gibbs_kprime_sweep", "slice", "scale")
  expect_error(
    MkPrime:::.SchedulePins(moveWeights, moveTypes, c(nni = 0.7)),
    "0.6"
  )
  # The documented example fits, and the auto-pins keep their shares.
  pins <- MkPrime:::.SchedulePins(moveWeights, moveTypes,
                                  c(nni = 0.3, spr = 0.2))
  expect_equal(pins, c(gibbs_kPrime = 0.3, slice_rate_log_sd = 0.1,
                       nni = 0.3, spr = 0.2))
  # A user pin on an auto-pinned move replaces its share, not adds to it.
  expect_equal(
    MkPrime:::.SchedulePins(moveWeights, moveTypes,
                            c(gibbs_kPrime = 0.1, nni = 0.75))[["nni"]],
    0.75
  )
  # Pinning every move is a complete schedule, not a starved one.
  expect_equal(
    sum(MkPrime:::.SchedulePins(moveWeights, moveTypes,
                                c(nni = 0.3, spr = 0.2, tree_length = 0.1))),
    1
  )
  expect_null(MkPrime:::.SchedulePins(moveWeights[c(1, 2, 5)],
                                      moveTypes[c(1, 2, 5)], NULL))
})

test_that(".PerturbMoveWeights never proposes a negative weight", {
  set.seed(7102)
  cands <- MkPrime:::.PerturbMoveWeights(
    c(a = 0.25, b = 0.25, c = 0.25, d = 0.25),
    pinnedWeights = c(a = 0.7, b = 0.45),
    moveNames = c("a", "b", "c", "d"),
    nPerturbations = 5L
  )
  for (w in cands) expect_true(all(w >= 0))
  # No budget remains for the free moves, so there is nothing to perturb.
  expect_length(cands, 0L)
})

test_that(".CheckMoveWeights rejects anything but a distribution", {
  expect_error(MkPrime:::.CheckMoveWeights(c(a = 1.2, b = -0.2)), "negative")
  expect_error(MkPrime:::.CheckMoveWeights(c(a = NaN, b = 1)), "finite")
  expect_error(MkPrime:::.CheckMoveWeights(c(a = 0.5, b = 0.4)), "sum to 1")
  expect_silent(MkPrime:::.CheckMoveWeights(c(a = 0.5, b = 0.5)))
})
