# Joint-move correlations are estimated at default warmup lengths (#301).

test_that(".AccumulateRhoSnapshot records rate_neo", {
  paramNames <- c("tree_length", "rate_loss", "rate_log_sd", "rate_neo")
  buf <- NULL
  set.seed(3011)
  for (i in 1:60) {
    state <- list(treeLength = rlnorm(1), rateLogSd = 0.8, rateLoss = 0.3,
                  rateNeo = rlnorm(1))
    buf <- MkPrime:::.AccumulateRhoSnapshot(buf, state, hasNeo = TRUE,
                                            paramNames)
  }
  expect_gt(min(buf[, "rate_neo"]), 0)
  expect_equal(unname(buf[60, "rate_neo"]), state$rateNeo)
  expect_true(is.finite(cor(log(buf[, "tree_length"]),
                            log(buf[, "rate_neo"]))))
  expect_equal(MkPrime:::.EstimateJointRhos(buf, hasNeo = TRUE)$rho_tl_rn,
               cor(log(buf[, "tree_length"]), log(buf[, "rate_neo"])))
})

test_that(".EstimateJointRhos keeps rhos it cannot re-estimate", {
  prior <- list(rho_tl_rls = 0.5, rho_tl_rl = -0.3, rho_tl_rn = 0.2)
  short <- cbind(tree_length = rlnorm(10), rate_log_sd = rlnorm(10))
  expect_equal(MkPrime:::.EstimateJointRhos(short, TRUE, prior = prior), prior)
  expect_equal(MkPrime:::.EstimateJointRhos(NULL, TRUE, prior = prior), prior)

  # Only the pair the samples can estimate moves.
  set.seed(3012)
  tl <- rlnorm(100)
  long <- cbind(tree_length = tl,
                rate_log_sd = exp(0.7 * log(tl) + rnorm(100, 0, 0.2)))
  rhos <- MkPrime:::.EstimateJointRhos(long, TRUE, prior = prior)
  expect_gt(rhos$rho_tl_rls, 0.5)
  expect_equal(rhos[c("rho_tl_rl", "rho_tl_rn")],
               prior[c("rho_tl_rl", "rho_tl_rn")])
})

test_that("a default-length run adapts every joint-move correlation", {
  skip_under_memcheck()
  skip_on_cran()
  set.seed(301)
  tree <- ape::rtree(20)
  tree$edge.length <- tree$edge.length / 5
  Simulate <- function(k) {
    as.character(ape::rTraitDisc(tree, k = k, states = as.character(seq_len(k) - 1L)))
  }
  Variable <- function(m) m[, apply(m, 2, function(x) length(unique(x)) > 1)]
  trans <- Variable(replicate(50, Simulate(sample(2:3, 1))))
  neo <- Variable(replicate(10, Simulate(2)))
  mat <- cbind(trans, neo)
  rownames(mat) <- tree$tip.label

  seen <- new.env()
  buildRhos <- MkPrime:::.BuildJointRhoMatrix
  local_mocked_bindings(
    .BuildJointRhoMatrix = function(chainRhos, moves, nChains) {
      seen$rhos <- chainRhos[[1]]
      buildRhos(chainRhos, moves, nChains)
    }
  )
  setTimeLimit(elapsed = 240, transient = TRUE)
  on.exit(setTimeLimit(elapsed = Inf, transient = TRUE), add = TRUE)
  result <- allow_warning(
    RunMkPrime(MatrixToPhyDat(mat), neomorphic = ncol(trans) + seq_len(ncol(neo)),
               mcmc = MkPrimeMCMC(nIter = 30000L, nRuns = 1L, nCore = 1L,
                                  maxTime = 180)),
    "without stabilisation|clamped"
  )
  expect_gt(nrow(result$samples), 0L)
  # Warmup ends in far fewer than the 50 batches its snapshots need, so the
  # estimate comes from Tuning; before #301 all three stayed at 0.
  expect_true(all(unlist(seen$rhos) != 0))
})
