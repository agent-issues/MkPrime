# #224: below the R API, tip i is paired with data row i by position. Each
# entry point must put the tree's tips into data order, or it scores and
# writes a different tree from the one it was given.

TipOrderOracle <- function() {
  set.seed(2240)
  truth <- ape::rtree(12, tip.label = paste0("t", 1:12))
  truth$edge.length <- runif(nrow(truth$edge), 0.05, 0.3)
  sim <- phangorn::simSeq(truth, l = 80, type = "USER", levels = c("0", "1"))
  mat <- as.character(sim)
  # Invariant characters would be dropped by MkPrimeData() but not by pml().
  mat <- mat[, apply(mat, 2, function(x) length(unique(x)) > 1L)]
  pd <- phangorn::phyDat(mat, type = "USER", levels = c("0", "1"))
  # Same topology and branch lengths; tips numbered out of data order.
  scrambled <- TreeTools::RenumberTips(truth, rev(truth$tip.label))
  list(truth = truth, scrambled = scrambled, pd = pd,
       mkd = MkPrimeData(pd))
}

# The sampler's logged likelihood must be that of the tree it writes, scored
# with the tips paired to the data by label.
ScoresWrittenTree <- function(tree, row, o, model) {
  row <- unlist(row)
  kPrime <- as.integer(row[startsWith(names(row), "kPrime_")])
  Param <- function(name) if (name %in% names(row)) row[[name]] else 1
  inDataOrder <- TreeTools::RenumberTips(tree, rownames(o$mkd$matrix))
  ll <- MkpLogLikelihood(inDataOrder, o$mkd, kPrime = kPrime,
                         rate_loss = Param("rate_loss"),
                         rate_log_sd = row[["rate_log_sd"]],
                         nCat = model$nCat, coding = model$coding,
                         rate_neo = Param("rate_neo"),
                         relabel = model$relabel)
  isTRUE(all.equal(ll, row[["log_likelihood"]], tolerance = 1e-8))
}

OracleRun <- function(o, tree, ...) {
  RunMkPrime(o$pd, tree, fixTopology = TRUE, nRuns = 1L, thin = 50L,
             maxWarmup = 50L, minWarmup = 50L, autoTune = FALSE, ...)
}


test_that("MkpLogLikelihood() does not depend on tip numbering (#224)", {
  skip_if_not_installed("phangorn")
  o <- TipOrderOracle()
  expect_false(identical(o$scrambled$tip.label, rownames(o$mkd$matrix)))

  LL <- function(tr) {
    MkpLogLikelihood(tr, o$mkd, coding = "none", rate_log_sd = 0,
                     relabel = FALSE)
  }
  expected <- phangorn::pml(o$truth, o$pd, model = "ER")$logLik
  expect_equal(LL(o$truth), expected, tolerance = 1e-8)
  expect_equal(LL(o$scrambled), expected, tolerance = 1e-8)
})


test_that("MkpLogLikelihood() rejects a tree whose tips are not the taxa", {
  skip_if_not_installed("phangorn")
  o <- TipOrderOracle()
  wrong <- o$truth
  wrong$tip.label[1] <- "stranger"
  expect_error(MkpLogLikelihood(wrong, o$mkd), "do not match")
})


test_that("RunMkPrime() writes trees with their own labels (#224)", {
  skip_if_not_installed("phangorn")
  o <- TipOrderOracle()
  set.seed(2241)
  result <- allow_warning(OracleRun(o, o$scrambled, nIter = 200L),
                          "without stabilisation")
  expect_true(ScoresWrittenTree(result$trees[[length(result$trees)]],
                                result$samples[nrow(result$samples), ],
                                o, result$model))
})


test_that("Default start tree is labelled correctly (#224)", {
  skip_if_not_installed("phangorn")
  skip_if_not_installed("TreeSearch")
  o <- TipOrderOracle()
  set.seed(2242)
  result <- allow_warning(OracleRun(o, NULL, nIter = 200L),
                          "without stabilisation")
  expect_true(ScoresWrittenTree(result$trees[[length(result$trees)]],
                                result$samples[nrow(result$samples), ],
                                o, result$model))
})


test_that("Resume keeps the labelling and the prior of the run (#224)", {
  skip_if_not_installed("phangorn")
  o <- TipOrderOracle()
  ckpFile <- tempfile(fileext = ".ckp")
  logFile <- tempfile(fileext = ".log")
  treeFile <- sub("\\.log$", "_trees.nwk", logFile)
  on.exit(unlink(c(ckpFile, logFile, treeFile)), add = TRUE)

  # expSteps is derived from the start tree, which the resume must not
  # rebuild: start from a poor tree so that a rebuilt one would differ.
  start <- o$scrambled
  start$tip.label <- sample(start$tip.label)
  set.seed(2243)
  first <- allow_warning(
    OracleRun(o, start, nIter = 200L, logFile = logFile,
              checkpointFile = ckpFile),
    "without stabilisation")
  ck <- readRDS(ckpFile)
  expect_gt(ck$model$expSteps, phangorn::parsimony(o$truth, o$pd))
  ck$mcmc$nIter <- 400L
  saveRDS(ck, ckpFile)

  # A model passed on resume (as the RB harness does) used to trigger a
  # fresh random-order start tree and a re-derived expSteps.
  resumed <- ResumeMkPrime(ckpFile, o$pd, model = MkPrimeModel())
  expect_identical(resumed$model$expSteps, ck$model$expSteps)
  written <- ape::read.tree(treeFile)
  logged <- ReadMkLog(logFile)
  expect_identical(length(written), nrow(logged))
  expect_gt(nrow(logged), nrow(first$samples))
  for (i in c(1L, nrow(logged))) {
    expect_true(ScoresWrittenTree(written[[i]], logged[i, ], o,
                                  resumed$model))
  }
})
