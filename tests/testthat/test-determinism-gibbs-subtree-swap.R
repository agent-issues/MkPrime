# MkPrime does not offer deterministic run output. This file tests an
# internal invariant used to detect state leaking between calls: for a
# FIXED move schedule, RunMkPrime() samples are a deterministic function of
# set.seed(), including on a repeat call in the same R session.
#
# The schedule itself is NOT seed-determined: during warmup
# .AdaptMoveWeights() scores each move by acceptances per wall-clock second,
# so two calls freeze slightly different weights and their RNG streams
# diverge.  That is a deliberate efficiency/reproducibility trade-off, so the
# tests pin every move weight (taken from a pilot run) to take timing out of
# the loop.  autoTune = FALSE skips the tuning phase, which also perturbs the
# schedule.
#
# The pilot run makes each comparison the 2nd and 3rd call in the session,
# so C++ static state carried across calls (e.g. counters in do_move_impl())
# fails the test.
skip_slow_tests()

library("ape")
library("TreeTools")

.DetTree <- function() {
  Preorder(read.tree(text =
    "((t1:0.1,t2:0.2):0.15,(t3:0.1,(t4:0.1,t5:0.2):0.12):0.18,t6:0.3);"))
}

.DetPd <- function() {
  MatrixToPhyDat(matrix(
    c(0, 1, 0, 1, 1, 0,
      0, 0, 1, 1, 0, 1,
      0, 1, 2, 0, 1, 0),
    nrow = 6, ncol = 3,
    dimnames = list(paste0("t", 1:6), NULL)
  ))
}

.DetMcmc <- function(gibbsSubtreeSwap, moveWeights = NULL) {
  MkPrimeMCMC(
    nRuns         = 1L,
    nIter         = 547L,
    thin          = 5L,
    maxWarmup     = 273L,
    minWarmup     = 273L,
    autoTune      = FALSE,
    gibbsSubtreeSwap = gibbsSubtreeSwap,
    moveWeights   = moveWeights
  )
}

# Pilot run fixes the schedule; two seeded runs under it must match exactly.
.ExpectReproducible <- function(gibbsSubtreeSwap) {
  tree <- .DetTree()
  pd   <- .DetPd()

  set.seed(1L)
  pilot <- RunMkPrime(pd, tree, mcmc = .DetMcmc(gibbsSubtreeSwap))
  expect_identical("gibbs_subtree_swap" %in% names(pilot$moveWeights),
                   gibbsSubtreeSwap)
  mcmc <- .DetMcmc(gibbsSubtreeSwap, moveWeights = pilot$moveWeights)

  set.seed(424242L)
  r1 <- RunMkPrime(pd, tree, mcmc = mcmc)
  set.seed(424242L)
  r2 <- RunMkPrime(pd, tree, mcmc = mcmc)

  expect_identical(r1$moveWeights, r2$moveWeights)
  expect_identical(r1$samples, r2$samples)
  expect_identical(r1$trees, r2$trees)
}

test_that("gibbsSubtreeSwap=TRUE is seed-reproducible under a pinned schedule", {
  .ExpectReproducible(TRUE)
})

test_that("gibbsSubtreeSwap=FALSE is seed-reproducible under a pinned schedule", {
  .ExpectReproducible(FALSE)
})
