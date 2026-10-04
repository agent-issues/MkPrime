# .FinalizeModel checks hand-edited prior fields exactly: LogPrior,
# cpp_log_prior and .LogZ0 compare them with identical(), so any other value
# silently selects a different prior, or a different one in R than in C++
# (#366).

.FinalizeFixture <- function() {
  tree <- read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 2, 2, 0, 1, 0, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  list(tree = tree, mkd = MkPrimeData(MatrixToPhyDat(mat)))
}

.Finalize <- function(model) {
  f <- .FinalizeFixture()
  MkPrime:::.FinalizeModel(model, f$tree, f$mkd)
}

test_that("finalisation rejects a kPrimePrior that is not an arm", {
  for (bad in list(NULL, "Geometric", "geom", c("geometric", "logseries"),
                   1)) {
    model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 1)
    model$kPrimePrior <- bad
    expect_error(.Finalize(model), "kPrimePrior")
  }
})

test_that("finalisation rejects a priorVariant other than the two values", {
  for (bad in list("Unconditional", "uncond", NA_character_)) {
    model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 1)
    model$priorVariant <- bad
    expect_error(.Finalize(model), "priorVariant")
  }
})

test_that("finalisation checks kprimeTruncK and empiricalNObs", {
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 1)
  model$kprimeTruncK <- 300L
  expect_error(.Finalize(model), "kprimeTruncK")
  model$kprimeTruncK <- 1L
  expect_error(.Finalize(model), "kprimeTruncK")
  model$kprimeTruncK <- 8
  expect_identical(.Finalize(model)$kprimeTruncK, 8L)

  model <- MkPrimeModel(kPrimePrior = "empirical_geometric", expSteps = 1)
  model$empiricalNObs <- list(pmf = c(0.5, 0.5))
  expect_error(.Finalize(model), "MkPrimeEmpiricalPrior")
})

test_that("valid models still finalise", {
  for (arm in c("empirical_geometric", "geometric", "beta_geometric",
                "logseries")) {
    model <- .Finalize(MkPrimeModel(kPrimePrior = arm, expSteps = 1))
    expect_identical(model$kPrimePrior, arm)
  }
  model <- MkPrimeModel(kPrimePrior = "geometric", expSteps = 1)
  model$priorVariant <- NULL
  expect_identical(.Finalize(model)$priorVariant, "unconditional")
})

test_that("the tree-length prior must be positive (#380)", {
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2), "2")) {
    expect_error(MkPrimeModel(treeLengthShape = bad), "treeLengthShape")
    expect_error(MkPrimeModel(treeLengthRate = bad), "treeLengthRate")
    expect_error(MkPrimeModel(expSteps = bad), "expSteps")
  }
  model <- MkPrimeModel(expSteps = 1)
  model$treeLengthShape <- -2
  expect_error(.Finalize(model), "treeLengthShape")
  model <- MkPrimeModel(expSteps = 1)
  model$treeLengthRate <- 0
  expect_error(.Finalize(model), "treeLengthRate")
  model <- MkPrimeModel()
  model$expSteps <- -1
  expect_error(.Finalize(model), "expSteps")
  expect_equal(.Finalize(MkPrimeModel(treeLengthShape = 3, expSteps = 2))$
                 treeLengthRate, 1.5)
})
