# Metropolis-coupled MCMC must leave the cold chain's target unchanged (#306).

.Mc3Data <- function() {
  # Two well-supported splits, two conflicting pairs and four autapomorphies:
  # every one of the 15 trees keeps some posterior mass.
  rows <- c(t1 = "11111100000101000",
            t2 = "11111100000010100",
            t3 = "00000000000100010",
            t4 = "00000011111010001",
            t5 = "00000011111000000")
  mat <- do.call(rbind, strsplit(rows, ""))
  rownames(mat) <- names(rows)
  list(tree = ape::read.tree(
         text = "((t1:0.1,t2:0.1):0.1,(t3:0.1,(t4:0.1,t5:0.1):0.1):0.1);"
       ),
       pd = MatrixToPhyDat(mat))
}

.Mc3Mcmc <- function(...) {
  do.call(MkPrimeMCMC, utils::modifyList(list(
    nRuns = 1L, nCore = 1L, thin = 10L, minWarmup = 2000L, maxWarmup = 2000L,
    autoTune = FALSE, heat = 0.3, maxTime = 120
  ), list(...)))
}

# Mean of `x` and its Monte Carlo standard error by non-overlapping batch
# means, which absorbs the autocorrelation an iid standard error would ignore.
.BatchMean <- function(x, nBatch = 20L) {
  batches <- split(x, cut(seq_along(x), nBatch, labels = FALSE))
  c(mean = mean(x), se = sd(vapply(batches, mean, 0)) / sqrt(nBatch))
}

test_that("tempering leaves the cold chain's marginals unchanged (#306)", {
  skip_on_cran()
  skip_under_memcheck()
  data <- .Mc3Data()
  Sample <- function(nChains) {
    allow_warning(
      RunMkPrime(data$pd, data$tree, mcmc = .Mc3Mcmc(
        nIter = 42000L, nChains = nChains
      )),
      "stabilis"
    )
  }
  set.seed(3061)
  single <- Sample(1L)
  set.seed(3063)
  coupled <- Sample(3L)

  # The comparison means nothing unless the cold chain really was coupled.
  expect_lt(coupled$betas[3], 0.5)
  expect_gt(sum(coupled$samples[, "swap_cold"]), 100)

  a <- single$samples
  b <- coupled$samples
  topologies <- union(a[, "topo_hash"], b[, "topo_hash"])
  Columns <- function(s) {
    cbind(
      s[, c("log_likelihood", "tree_length", "rate_log_sd")],
      vapply(topologies, function(h) as.numeric(s[, "topo_hash"] == h),
             numeric(nrow(s)))
    )
  }
  sa <- apply(Columns(a), 2, .BatchMean)
  sb <- apply(Columns(b), 2, .BatchMean)
  z <- (sb["mean", ] - sa["mean", ]) / sqrt(sa["se", ]^2 + sb["se", ]^2)
  # Over 12 seed pairs on main, |z| never exceeded 3.4 across these 18
  # statistics, though batch means on 4000 draws understate the error and the
  # single chain undersamples the tree-length tail.  Scoring the swap on the
  # tempered log-likelihood moves both the cold log-likelihood and the
  # topology frequencies by |z| > 35, so a wide cap costs no power.
  expect_lt(max(abs(z)), 6)
})

test_that("a resumed run starts from the ladder it saved (#306)", {
  skip_on_cran()
  skip_under_memcheck()
  data <- .Mc3Data()
  initial <- MkPrime:::.BuildTemperatureLadder(3L, 0.3)
  adapt <- MkPrime:::.AdaptTemperatures
  runBatch <- MkPrime:::.RunBatch
  seen <- new.env()
  local_mocked_bindings(
    .AdaptTemperatures = function(...) {
      out <- adapt(...)
      seen$adapted <- as.numeric(out)
      out
    },
    .RunBatch = function(...) {
      if (is.null(seen$resumed)) seen$resumed <- as.numeric(..3)
      runBatch(...)
    }
  )
  # Stop mid-warmup, with the ladder still adapting, and in Sample, where the
  # ladder is frozen.
  for (stopAt in c(1500L, 3500L)) {
    ckpFile <- tempfile(fileext = ".ckp")
    on.exit(unlink(ckpFile), add = TRUE)
    seen$adapted <- NULL
    set.seed(stopAt)
    allow_warning(
      RunMkPrime(data$pd, data$tree, mcmc = .Mc3Mcmc(
        nIter = stopAt, nChains = 3L, minWarmup = 3000L, maxWarmup = 3000L,
        checkpointFile = ckpFile
      )),
      "stabilis|without drawing"
    )
    saved <- seen$adapted
    expect_false(isTRUE(all.equal(saved, initial)))

    seen$resumed <- NULL
    allow_warning(
      ResumeMkPrime(ckpFile, data$pd, mcmc = list(nIter = 4000L)),
      "stabilis"
    )
    expect_equal(seen$resumed, saved)
  }
})
