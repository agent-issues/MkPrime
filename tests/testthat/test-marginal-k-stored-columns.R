# Under marginal_k the kPrime_ columns are stripped from the stored trace, so
# every stored column must still carry its own parameter, and trees must be
# built from the br_ columns (#266).

.MkscFixture <- function() {
  set.seed(266)
  tips <- paste0("t", 1:8)
  trans <- matrix(sample(0:2, 8 * 5, TRUE), 8)
  trans[1:3, ] <- matrix(0:2, 3, 5)
  neo <- matrix(c(0, 0, 1, 1, 0, 1, 1, 0), 8)
  mat <- cbind(trans, neo)
  dimnames(mat) <- list(tips, NULL)
  tree <- read.tree(text = paste0(
    "(((t1:0.05,t2:0.07):0.04,(t3:0.06,t4:0.05):0.03):0.05,",
    "((t5:0.04,t6:0.06):0.05,(t7:0.05,t8:0.04):0.06):0.04);"
  ))
  list(pd = MatrixToPhyDat(mat), tree = tree, neomorphic = 6L,
       model = MkPrimeModel(kPrimePrior = "geometric",
                            likelihoodMode = "marginal_k", coding = "none",
                            relabel = FALSE))
}

.MkscMcmc <- function(nIter, ckp) {
  MkPrimeMCMC(nRuns = 1L, nChains = 1L, nIter = nIter, thin = 10L,
              maxWarmup = 100L, minWarmup = 100L, autoTune = FALSE,
              checkEvery = 100L, nCore = 1L, maxTime = 60,
              checkpointFile = ckp)
}

.ExpectStoredColumnsCoherent <- function(samples, trees, ckp) {
  expect_false(any(grepl("^kPrime_", colnames(samples))))
  br <- samples[, grepl("^br_", colnames(samples)), drop = FALSE]
  expect_equal(unname(rowSums(br)), rep(1, nrow(br)), tolerance = 1e-8)
  expect_true(all(samples[, "p"] > 0 & samples[, "p"] < 1))

  # The final stored row was drawn at the final state, which the checkpoint
  # holds; compare column by column against it.
  chain <- readRDS(ckp)$runs[[1]]$chains[[1]]
  last <- samples[nrow(samples), ]
  expect_equal(last[["log_likelihood"]], chain$log_lik)
  expect_equal(last[["tree_length"]], chain$tree_length)
  expect_equal(last[["rate_loss"]], chain$rate_loss)
  expect_equal(last[["rate_log_sd"]], chain$rate_log_sd)
  expect_equal(last[["p"]], chain$p)
  expect_equal(last[["rate_neo"]], chain$rate_neo)
  expect_equal(unname(last[grepl("^br_", names(last))]), chain$rel_br_lengths)

  trees <- utils::tail(trees, nrow(samples))
  expect_equal(length(trees), nrow(samples))
  for (i in seq_along(trees)) {
    expect_equal(sort(trees[[i]]$edge.length),
                 sort(unname(br[i, ]) * samples[i, "tree_length"]),
                 tolerance = 1e-8)
  }
}

test_that("marginal_k stores each parameter in its own column (#266)", {
  skip_on_cran()
  fx <- .MkscFixture()
  td <- tempfile("mkp_mk_cols_")
  dir.create(td)
  on.exit(unlink(td, recursive = TRUE), add = TRUE)
  ckp <- file.path(td, "a.ckp")

  set.seed(1)
  result <- allow_warning(
    RunMkPrime(fx$pd, fx$tree, neomorphic = fx$neomorphic, model = fx$model,
               mcmc = .MkscMcmc(200L, ckp)),
    "without stabilisation"
  )
  .ExpectStoredColumnsCoherent(result$samples, result$trees, ckp)

  resumed <- allow_warning(
    ResumeMkPrime(ckp, fx$pd, fx$tree, neomorphic = fx$neomorphic,
                  mcmc = list(nIter = 300L)),
    "without stabilisation"
  )
  expect_gt(readRDS(ckp)$iter, 200L)
  .ExpectStoredColumnsCoherent(as.matrix(ReadMkLog(resumed$logFile)),
                               resumed$trees, ckp)
})
