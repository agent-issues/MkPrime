# Prior-input resolution (#256) and validation (#257).

PriorInputFixture <- function() {
  tree <- Preorder(ape::read.tree(text = "((a:1,b:1):1,(c:1,(d:1,e:1):1):1);"))
  mat <- matrix(c(0, 1, 2, 0, 1,
                  0, 1, 2, 3, 0,
                  0, 0, 1, 1, 0), nrow = 5,
                dimnames = list(tree$tip.label, NULL))
  pd <- phangorn::phyDat(mat, type = "USER", levels = as.character(0:3))
  list(tree = tree, mkd = MkPrimeData(pd))
}

# R LogPrior and the C++ prior built by .InitMcmcData, at a state with k' above
# kObs and kObs > 2, where Model A and Model B differ.
PriorPair <- function(model, f, rawModel = model) {
  modf <- MkPrime:::.FinalizeModel(model, f$tree, f$mkd)
  state <- MkPrime:::.InitState(f$tree, f$mkd, modf)
  state$p <- 0.3
  state$kPrime <- as.integer(f$mkd$kObs + 2L)
  dp <- MkPrime:::.InitMcmcData(f$mkd, modf)
  sp <- MkPrime:::.InitMcmcChain(state)
  fill_partition_cache(dp, sp)
  c(rFinal = LogPrior(state, modf, f$mkd),
    rRaw = LogPrior(state, rawModel, f$mkd),
    cpp = eval_log_prior_cpp(dp, sp))
}

test_that("a priorVariant-less model gets the same prior in R and C++ (#256)", {
  skip_if_not_installed("phangorn")
  f <- PriorInputFixture()
  expect_true(any(f$mkd$kObs > 2L))
  for (arm in c("geometric", "empirical_geometric")) {
    model <- MkPrimeModel(kPrimePrior = arm, expSteps = 1)
    model$priorVariant <- NULL
    lp <- PriorPair(model, f)
    expect_equal(lp[["rFinal"]], lp[["cpp"]], tolerance = 1e-10, label = arm)
    expect_equal(lp[["rRaw"]], lp[["cpp"]], tolerance = 1e-10, label = arm)
    # Anti-vacuity: the conditional variant gives a different value here.
    cond <- PriorPair(MkPrimeModel(kPrimePrior = arm, expSteps = 1,
                                   priorVariant = "conditional"), f)
    expect_gt(abs(cond[["cpp"]] - lp[["cpp"]]), 1e-3)
  }
})

test_that("finalisation fills missing prior fields as .InitMcmcData does (#256)", {
  skip_if_not_installed("phangorn")
  f <- PriorInputFixture()
  model <- MkPrimeModel(kPrimePrior = "logseries", expSteps = 1)
  model$kprimeLogseriesC <- NULL
  modf <- MkPrime:::.FinalizeModel(model, f$tree, f$mkd)
  expect_equal(modf$kprimeLogseriesC, 0.7)
  lp <- PriorPair(model, f)
  expect_equal(lp[["rRaw"]], lp[["cpp"]], tolerance = 1e-10)

  eg <- MkPrimeModel(kPrimePrior = "empirical_geometric", expSteps = 1)
  expect_null(eg$empiricalNObs)
  expect_s3_class(MkPrime:::.FinalizeModel(eg, f$tree, f$mkd)$empiricalNObs,
                  "MkPrimeEmpiricalPrior")
})

test_that("a mutated logseries c errors at finalisation (#256)", {
  skip_if_not_installed("phangorn")
  f <- PriorInputFixture()
  model <- MkPrimeModel(kPrimePrior = "logseries", expSteps = 1)
  model$kprimeLogseriesC <- 2
  expect_error(MkPrime:::.FinalizeModel(model, f$tree, f$mkd), "\\(0, 1\\)")
  # Only the logseries arm reads c.
  geo <- suppressWarnings(MkPrimeModel(kPrimePrior = "geometric", expSteps = 1,
                                       kprimeLogseriesC = 2))
  expect_no_error(MkPrime:::.FinalizeModel(geo, f$tree, f$mkd))
})

test_that("k'-prior hyperparameters must be positive finite scalars (#257)", {
  hyperMsg <- "kprimeHyperA.*kprimeHyperB.*positive finite"
  for (arm in c("geometric", "empirical_geometric")) {
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperA = -1), hyperMsg)
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperA = 0), hyperMsg)
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperB = 0), hyperMsg)
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperB = Inf), hyperMsg)
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperA = NA), hyperMsg)
    expect_error(MkPrimeModel(kPrimePrior = arm, kprimeHyperA = c(1, 2)),
                 hyperMsg)
  }
  bgMsg <- "kprimeAlpha.*kprimeBeta.*positive finite"
  expect_error(MkPrimeModel(kPrimePrior = "beta_geometric", kprimeAlpha = Inf),
               bgMsg)
  expect_error(MkPrimeModel(kPrimePrior = "beta_geometric",
                            kprimeBeta = NA_real_), bgMsg)
  expect_error(MkPrimeModel(kPrimePrior = "beta_geometric",
                            kprimeAlpha = c(1, 2)), bgMsg)
})

test_that("mutated hyperparameters error at finalisation (#257)", {
  skip_if_not_installed("phangorn")
  f <- PriorInputFixture()
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 1)
  model$kprimeHyperA <- -1
  expect_error(MkPrime:::.FinalizeModel(model, f$tree, f$mkd),
               "positive finite")
})

test_that("MkPrimeEmpiricalPrior rejects non-finite and P(2) = 0 bodies (#257)", {
  expect_error(MkPrimeEmpiricalPrior(body = c(Inf, 1), tail_decay = 0.5),
               "finite")
  expect_error(MkPrimeEmpiricalPrior(body = c(NA, 1), tail_decay = 0.5),
               "finite")
  expect_error(MkPrimeEmpiricalPrior(body = c(0.5, 0.5), tail_decay = NA),
               "tail_decay")
  expect_error(MkPrimeEmpiricalPrior(body = c(0, 1)), "N_obs = 2")
  expect_error(MkPrimeEmpiricalPrior(body = c(0, 1), tail_decay = 0.5),
               "N_obs = 2")
})

test_that("LogPrior is -Inf, not an error, for a NULL p (#257)", {
  skip_if_not_installed("phangorn")
  f <- PriorInputFixture()
  for (arm in c("geometric", "empirical_geometric")) {
    modf <- MkPrime:::.FinalizeModel(
      MkPrimeModel(kPrimePrior = arm, expSteps = 1), f$tree, f$mkd)
    state <- MkPrime:::.InitState(f$tree, f$mkd, modf)
    state$p <- NULL
    expect_equal(LogPrior(state, modf, f$mkd), -Inf, label = arm)
  }
})
