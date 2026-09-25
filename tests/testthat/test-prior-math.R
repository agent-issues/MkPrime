# Accuracy of the shared log(1 - exp(x)) helper behind the truncation
# normalisers, in R (.Log1mExp) and C++ (mkp::log1m_exp).

test_that("log(1 - exp(x)) is accurate at both tails", {
  K <- 200L
  # Tiny p: 1 - (1-p)^(K-1) ~ (K-1) p underflows under log1p(-exp(x)).
  xSmall <- (K - 1) * log1p(-1e-19)
  # p = 0.5: 1 - 2^-(K-1) rounds to 1 under log(-expm1(x)).
  xLarge <- (K - 1) * log1p(-0.5)
  exact <- c(log((K - 1) * 1e-19), -2^-(K - 1))
  for (f in list(MkPrime:::.Log1mExp, MkPrime:::log1m_exp_cpp)) {
    got <- f(c(xSmall, xLarge))
    expect_true(all(is.finite(got)))
    expect_lt(max(abs(got / exact - 1)), 1e-12)
  }
  expect_identical(MkPrime:::log1m_exp_cpp(c(xSmall, xLarge, -log(2), -0.1)),
                   MkPrime:::.Log1mExp(c(xSmall, xLarge, -log(2), -0.1)))
})


test_that("geometric k' prior stays finite at tiny p in R and C++", {
  tree <- Preorder(
    read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);"))
  mat <- matrix(c(0, 1, 0, 1,
                  0, 1, 2, 0), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  for (variant in c("unconditional", "conditional")) {
    model <- suppressMessages(MkPrimeModel(
      expSteps = 10, kPrimePrior = "geometric", priorVariant = variant,
      kprimeTruncK = 200L))
    model <- MkPrime:::.FinalizeModel(model, tree, mkd)
    dataPtr <- MkPrime:::.InitMcmcData(mkd, model)
    state <- list(
      tree = tree, tree_length = 0.5,
      rel_br_lengths = tree$edge.length / sum(tree$edge.length),
      rate_loss = 1.0, rate_log_sd = 0.2, rate_neo = 1.0,
      kPrime = c(2L, 3L), p = 1e-19, log_lik = 0.0, log_prior = 0.0
    )
    lpR <- MkPrime:::LogPrior(state, model, mkd)
    lpCpp <- eval_log_prior_cpp(dataPtr, MkPrime:::.InitMcmcChain(state))
    expect_true(is.finite(lpR), info = variant)
    expect_equal(lpCpp, lpR, tolerance = 1e-12, info = variant)
  }
})
