# Unit tests for .ConvergenceStreak()

test_that(".ConvergenceStreak needs consecutive checkpoints", {
  first <- MkPrime:::.ConvergenceStreak(0L, TRUE)
  expect_equal(first, list(streak = 1L, stop = FALSE))
  expect_equal(MkPrime:::.ConvergenceStreak(first$streak, TRUE),
               list(streak = 2L, stop = TRUE))
})

test_that(".ConvergenceStreak resets on a failing checkpoint", {
  expect_equal(MkPrime:::.ConvergenceStreak(1L, FALSE),
               list(streak = 0L, stop = FALSE))
  expect_equal(MkPrime:::.ConvergenceStreak(1L, NA), list(streak = 0L,
                                                          stop = FALSE))
})
