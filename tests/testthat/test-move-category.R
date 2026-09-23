# Every move the scheduler can run must be reported under a named category in
# the frozen-weight summary, rather than falling through to "Other".

.OtherLine <- function(moveNames) {
  w <- rep(1 / length(moveNames), length(moveNames))
  plain <- MkPrime:::.FormatMoveWeightsPlain(w, moveNames)
  grep("^Other:", strsplit(plain, "\n", fixed = TRUE)[[1]], value = TRUE)
}

test_that("every move .BuildMoves can build has a category", {
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L)
  arms <- expand.grid(
    kPrimePrior = c("geometric", "empirical_geometric", "beta_geometric",
                    "logseries"),
    likelihoodMode = c("sampled_k", "marginal_k"),
    qHeterogeneity = c(FALSE, TRUE),
    fixTopology = c(FALSE, TRUE),
    stringsAsFactors = FALSE
  )
  for (i in seq_len(nrow(arms))) {
    arm <- arms[i, ]
    moves <- suppressMessages(MkPrime:::.BuildMoves(
      20L, 5L, TRUE, mcmc,
      fixTopology = arm$fixTopology,
      kPrimePrior = arm$kPrimePrior,
      qHeterogeneity = arm$qHeterogeneity,
      likelihoodMode = arm$likelihoodMode
    ))
    expect_identical(.OtherLine(vapply(moves, `[[`, character(1), "name")),
                     character(0), label = paste(unlist(arm), collapse = "/"))
  }
})

test_that("partitioned per-class moves have a category", {
  mcmc <- MkPrimeMCMC(nIter = 200L, minWarmup = 100L)
  spec <- list(partition = rep(1:2, length.out = 10L), nClasses = 2L,
               unlink = c("shape", "ratemultiplier"))
  moves <- MkPrime:::.BuildMovesPartitioned(20L, 5L, TRUE, mcmc, spec)
  moveNames <- vapply(moves, `[[`, character(1), "name")
  expect_true("scale_class_rate_log_sd_2" %in% moveNames)
  expect_identical(.OtherLine(moveNames), character(0))
})

test_that("every registered move has a category", {
  expect_identical(.OtherLine(names(MkPrime:::.kMoveTypes)), character(0))
})
