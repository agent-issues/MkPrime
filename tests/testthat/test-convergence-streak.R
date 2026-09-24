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

test_that(".ConvergenceStreak holds on a check that saw no new samples (#199)", {
  expect_equal(MkPrime:::.ConvergenceStreak(1L, TRUE, fresh = FALSE),
               list(streak = 1L, stop = FALSE))
  expect_equal(MkPrime:::.ConvergenceStreak(0L, FALSE, fresh = FALSE),
               list(streak = 0L, stop = FALSE))
})

test_that(".CheckConvergenceFromLogs reports each log's row count (#199)", {
  paramNames <- c("log_posterior", "tree_length")
  set.seed(1990)
  logs <- replicate(2L, tempfile(fileext = ".log"))
  on.exit(unlink(logs), add = TRUE)
  for (i in 1:2) {
    nRow <- 20L + i
    write.table(data.frame(Sample = seq_len(nRow), log_posterior = rnorm(nRow),
                           tree_length = rnorm(nRow)),
                logs[i], sep = "\t", row.names = FALSE, quote = FALSE)
  }
  res <- MkPrime:::.CheckConvergenceFromLogs(logs, paramNames,
                                             list(minEss = 5))
  expect_equal(res$nRows, c(21L, 22L))
})
