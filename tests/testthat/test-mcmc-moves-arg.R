toggleNames <- c(
  "gibbsSpr", "gibbsSubtreeSwap", "tbr", "pSpr", "joint2d",
  "weightedBranchScale", "weightedSpr", "weightedSubtreeSwap",
  "blockGibbsBranch", "gibbsPMarginal", "dirichletBranch", "localDirichlet"
)

test_that("toggleNames covers every move-toggle argument", {
  expect_setequal(toggleNames, .kMoveToggles)
})

test_that("moves = list() leaves every toggle at its default", {
  expect_identical(MkPrimeMCMC(moves = list()), MkPrimeMCMC())
  defaults <- unclass(MkPrimeMCMC())[toggleNames]
  expect_identical(
    unlist(defaults),
    c(gibbsSpr = TRUE, gibbsSubtreeSwap = TRUE, tbr = TRUE, pSpr = TRUE,
      joint2d = TRUE, weightedBranchScale = FALSE, weightedSpr = FALSE,
      weightedSubtreeSwap = FALSE, blockGibbsBranch = FALSE,
      gibbsPMarginal = FALSE, dirichletBranch = TRUE, localDirichlet = TRUE)
  )
})

test_that("moves matches the individual arguments field for field", {
  viaMoves <- MkPrimeMCMC(
    moves = list(tbr = FALSE, weightedSpr = TRUE, localDirichlet = FALSE)
  )
  viaArgs <- MkPrimeMCMC(tbr = FALSE, weightedSpr = TRUE,
                         localDirichlet = FALSE)
  expect_identical(viaMoves, viaArgs)
  expect_false(viaMoves$tbr)
  expect_true(viaMoves$weightedSpr)
})

test_that("every toggle can be set through moves", {
  flipped <- !unlist(unclass(MkPrimeMCMC())[toggleNames])
  mcmc <- MkPrimeMCMC(moves = as.list(flipped))
  expect_identical(unlist(unclass(mcmc)[toggleNames]), flipped)
})

test_that("moves rejects malformed input", {
  expect_error(MkPrimeMCMC(moves = list(TRUE)), "named")
  expect_error(MkPrimeMCMC(moves = list(tbr = TRUE, FALSE)), "named")
  expect_error(MkPrimeMCMC(moves = c(tbr = TRUE)), "list")
  expect_error(MkPrimeMCMC(moves = list(tbr = "yes")), "tbr")
  expect_error(MkPrimeMCMC(moves = list(tbr = NA)), "tbr")
  expect_error(MkPrimeMCMC(moves = list(tbr = c(TRUE, FALSE))), "tbr")
  expect_error(MkPrimeMCMC(moves = list(tbr = TRUE, tbr = FALSE)), "tbr")
})

test_that("moves rejects unknown names and suggests the nearest", {
  expect_error(MkPrimeMCMC(moves = list(tbrr = TRUE)), "tbr")
  expect_error(MkPrimeMCMC(moves = list(nni = TRUE)), "nni")
  expect_error(MkPrimeMCMC(moves = list(gibbsSpr = TRUE, bogus = FALSE)),
               "bogus")
})

test_that("a move given in moves and as an argument must agree", {
  expect_error(MkPrimeMCMC(tbr = TRUE, moves = list(tbr = FALSE)), "tbr")
  expect_identical(
    MkPrimeMCMC(tbr = FALSE, moves = list(tbr = FALSE)),
    MkPrimeMCMC(tbr = FALSE)
  )
  expect_error(MkPrimeMCMC(tbr = TRUE, moves = list(tbr = FALSE)),
               "moves[$]tbr")
  wrapper <- function(tbr = TRUE, ...) MkPrimeMCMC(tbr = tbr, ...)
  expect_error(wrapper(moves = list(tbr = FALSE)), "moves[$]tbr")
})

test_that("moves combines with arguments for other toggles", {
  expect_identical(
    MkPrimeMCMC(pSpr = FALSE, moves = list(tbr = FALSE)),
    MkPrimeMCMC(pSpr = FALSE, tbr = FALSE)
  )
})
