# #32: with only `minTreeEss` configured, the progress ticker showed no ETA.
#
# The call site passed `mcmc$minEss` as the target unconditionally. With only a
# tree-ESS criterion that is NULL, so .EstimateEta() bailed; and the branch that
# should have substituted the tree criterion was itself gated on `minEss` being
# set, so it could never fire. The selection now lives in .EtaCriterion().

test_that(".EstimateEta declines a NULL target rather than guessing", {
  # This is the exact call the old site made when only minTreeEss was set.
  expect_null(.EstimateEta(50, NULL, 10))
  expect_null(.EstimateEta(NULL, 200, 10))
})

test_that("only minEss configured: the scalar criterion is used", {
  crit <- .EtaCriterion(
    diagCheck = list(minEss = 50, treeEss = NA_real_),
    mcmc = list(minEss = 200)
  )
  expect_identical(crit$current, 50)
  expect_identical(crit$target, 200)
})

test_that("only minTreeEss configured: the tree criterion is used", {
  crit <- .EtaCriterion(
    diagCheck = list(minEss = 50, treeEss = 30),
    mcmc = list(minTreeEss = 200)
  )
  expect_identical(crit$current, 30)
  expect_identical(crit$target, 200)

  # The point of the issue: an ETA is now actually produced.
  expect_false(is.null(.EstimateEta(crit$current, crit$target, 10)))
})

test_that("both configured: the binding (worse-ratio) criterion wins", {
  # Scalar is at 50/200 = 0.25; tree at 180/200 = 0.9. Scalar binds.
  crit <- .EtaCriterion(
    diagCheck = list(minEss = 50, treeEss = 180),
    mcmc = list(minEss = 200, minTreeEss = 200)
  )
  expect_identical(crit$target, 200)
  expect_identical(crit$current, 50)

  # Tree at 20/200 = 0.1 is worse than scalar at 190/200 = 0.95. Tree binds.
  crit2 <- .EtaCriterion(
    diagCheck = list(minEss = 190, treeEss = 20),
    mcmc = list(minEss = 200, minTreeEss = 200)
  )
  expect_identical(crit2$current, 20)
})

test_that("an unassessable window yields no criterion, not a wrong one", {
  # minEss is NA when every monitored parameter is constant (CONV-001), and
  # treeEss is NA whenever tree ESS was skipped.
  crit <- .EtaCriterion(
    diagCheck = list(minEss = NA_real_, treeEss = NA_real_),
    mcmc = list(minEss = 200, minTreeEss = 200)
  )
  expect_null(crit$current)
  expect_null(crit$target)
  expect_null(.EstimateEta(crit$current, crit$target, 10))
})
