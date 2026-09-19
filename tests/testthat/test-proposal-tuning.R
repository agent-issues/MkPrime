# Regression tests for the proposal / adaptive-tuning layer.
#
#   #4  the R fallback for logit_scale_p must use the same kernel as the C++
#       path, so the two are comparable at matched `scale_logit_p`;
#   #5  the scalar p moves' acceptance target must suit the Bactrian kernel
#       they actually use, not the Gaussian 1-D optimum;
#   #6  .kMoveTypes and the moveWeights whitelist must be one surface;
#   #7  the generic branch of the step-size tuner must be bounded above.


# --- #4: the R fallback uses the C++ kernel ------------------------------

test_that("logit_scale_p draws the same Bactrian step as the C++ path", {
  x <- 0.4
  sigma <- 1.3
  # bactrian_draws() is the Rcpp export of the same bactrian_perturbation()
  # that src/mcmc.cpp case 30 calls, so a shared seed pins both kernel shape
  # and kernel scale. A uniform or normal step fails this outright.
  set.seed(910L)
  expected <- qlogis(x) + sigma * MkPrime:::bactrian_draws(1L)
  set.seed(910L)
  got <- qlogis(MkPrime:::.ProposeLogitScale(x, tuning = sigma)$value)
  expect_equal(got, expected)
})

test_that("logit_scale_p step scale matches the C++ path", {
  sigma <- 1.3
  set.seed(20260919)
  steps <- vapply(seq_len(2000L), function(i) {
    qlogis(MkPrime:::.ProposeLogitScale(0.4, tuning = sigma)$value) -
      qlogis(0.4)
  }, numeric(1))
  # bactrian_perturbation() is variance-matched to Uniform(-0.5, 0.5), so the
  # step has sd sigma / sqrt(12) -- not sigma, as a normal step would give.
  expect_equal(sd(steps), sigma / sqrt(12), tolerance = 0.05)
  expect_equal(mean(steps), 0, tolerance = 0.05)
})

test_that("logit_scale_p Hastings ratio is the logit Jacobian", {
  set.seed(4L)
  x <- 0.3
  prop <- MkPrime:::.ProposeLogitScale(x, tuning = 0.7)
  expect_equal(prop$logHastings,
               log(prop$value * (1 - prop$value)) - log(x * (1 - x)))
})

test_that(".DoMove routes logit_scale_p through the shared proposal", {
  skip_if_not_installed("TreeTools")
  tree <- ape::read.tree(text = "((t1:0.1,t2:0.2):0.15,(t3:0.1,t4:0.3):0.2);")
  mat <- matrix(c(0, 1, 0, 1, 0, 0, 1, 1), 4, 2,
                dimnames = list(paste0("t", 1:4), NULL))
  mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(), tree, mkd)
  state <- MkPrime:::.InitState(TreeTools::Preorder(tree), mkd, model)

  seen <- NULL
  local_mocked_bindings(
    .ProposeLogitScale = function(x, tuning = 1.0) {
      seen <<- tuning
      list(value = x, logHastings = 0)
    },
    .package = "MkPrime"
  )
  move <- list(name = "mh_logit_p", type = "logit_scale_p", target = "p",
               weight = 1, dim = 1L)
  MkPrime:::.DoMove(move, state, mkd, model, list(scale_logit_p = 0.75))
  expect_equal(seen, 0.75)
})


# --- #5: acceptance target suits the Bactrian kernel ---------------------

test_that("scalar p moves keep the Bactrian-appropriate 0.35 target", {
  moves <- list(list(name = "mh_logit_p", type = "logit_scale_p",
                     target = "p", weight = 1, dim = 1L))
  # Measured on a 1-D N(0, 1) target, the Bactrian kernel's ESS/iteration
  # peaks near 0.30 acceptance and 0.44 costs ~14%, so acceptance above 0.35
  # must widen the step, not narrow it.
  out <- MkPrime:::.AdaptTuning(
    list(scale_logit_p = 1.0),
    c(mh_logit_p = 40), c(mh_logit_p = 100), moves
  )
  expect_gt(out$scale_logit_p, 1.0)
})


# --- #6: one move-name surface -------------------------------------------

test_that("MkPrimeMCMC accepts every registered move name", {
  registered <- names(MkPrime:::.kMoveTypes)
  bad <- registered[vapply(registered, function(nm) {
    inherits(try(MkPrimeMCMC(nIter = 10000L,
                             moveWeights = stats::setNames(0.01, nm)),
                 silent = TRUE), "try-error")
  }, logical(1))]
  expect_equal(bad, character(0))
})

test_that("MkPrimeMCMC accepts every name .BuildMoves can actually build", {
  mcmc <- MkPrimeMCMC(nIter = 10000L)
  built <- character(0)
  for (prior in c("geometric", "empirical_geometric", "beta_geometric",
                  "logseries")) {
    for (fixTop in c(FALSE, TRUE)) {
      for (mode in c("sampled_k", "marginal_k")) {
        mv <- MkPrime:::.BuildMoves(20, 5, TRUE, mcmc, fixTopology = fixTop,
                                    kPrimePrior = prior,
                                    likelihoodMode = mode)
        built <- union(built, vapply(mv, `[[`, character(1), "name"))
      }
    }
  }
  spec <- list(partition = rep(1:2, length.out = 10), nClasses = 2L,
               unlink = c("shape", "ratemultiplier"))
  mv <- MkPrime:::.BuildMovesPartitioned(20, 5, TRUE, mcmc, spec)
  built <- union(built, vapply(mv, `[[`, character(1), "name"))
  # Includes the per-class instances scale_class_rate_log_sd_1 / _2.
  expect_true(any(grepl("_[0-9]+$", built)))

  bad <- built[vapply(built, function(nm) {
    inherits(try(MkPrimeMCMC(nIter = 10000L,
                             moveWeights = stats::setNames(0.01, nm)),
                 silent = TRUE), "try-error")
  }, logical(1))]
  expect_equal(bad, character(0))
})

test_that("MkPrimeMCMC still rejects names no move can carry", {
  expect_error(MkPrimeMCMC(nIter = 10000L, moveWeights = c(not_a_move = 0.1)),
               "unknown move name")
  # A per-class suffix on a type that has no per-class instances.
  expect_error(MkPrimeMCMC(nIter = 10000L, moveWeights = c(nni_2 = 0.1)),
               "unknown move name")
})

test_that("an inert pinned weight warns rather than vanishing", {
  expect_warning(
    MkPrime:::.ResolvePinnedWeights(c(scale_hyper_tau = 0.3), c("nni", "spr")),
    "not in move pool"
  )
})

test_that("the dead mh_p move is no longer registered", {
  expect_false("mh_p" %in% names(MkPrime:::.kMoveTypes))
  expect_false("mh_p" %in% MkPrime:::.ValidMoveNames())
})

test_that("every registered move code is dispatched by do_move_impl", {
  src <- test_path("..", "..", "src", "mcmc.cpp")
  skip_if_not(file.exists(src), "src/ absent from the installed package")
  txt <- readLines(src, warn = FALSE)
  # Restrict to do_move_impl's own switch: other switches in this file reuse
  # the same small integers and would make the check vacuous.
  start <- grep("^static bool do_move_impl", txt)
  expect_length(start, 1L)
  body <- txt[seq(start, length(txt))]
  sw <- grep("  switch (moveType) {", body, fixed = TRUE)[1]
  expect_false(is.na(sw))
  body <- body[seq(sw, length(body))]
  # The switch ends at the function's closing brace in column 1.
  body <- body[seq_len(grep("^}", body)[1])]
  caseLines <- grep("^    case [0-9]+:", body, value = TRUE)
  implemented <- unique(as.integer(sub(":.*", "",
                                       sub("^    case ", "", caseLines))))
  # The slice families are branched on before the switch, in the run loop.
  sliceCodes <- unique(as.integer(sub("moveType == ", "",
    regmatches(txt, regexpr("moveType == (19|29)", txt)))))
  expect_setequal(sliceCodes, c(19L, 29L))
  expect_equal(setdiff(MkPrime:::.kMoveTypes, c(implemented, sliceCodes)),
               integer(0))
})


# --- #7: the step-size tuner is bounded above ----------------------------

test_that("relentless acceptance cannot run the generic step to infinity", {
  moves <- list(list(name = "rate_loss", type = "scale",
                     target = "rate_loss", weight = 1, dim = 1L))
  tuning <- list(scale_rate_loss = 0.5)
  for (i in seq_len(2000L)) {
    tuning <- MkPrime:::.AdaptTuning(
      tuning, c(rate_loss = 100), c(rate_loss = 100), moves
    )
  }
  expect_true(is.finite(tuning$scale_rate_loss))
  expect_lte(tuning$scale_rate_loss, 10)
})

test_that("the generic step floor still holds", {
  moves <- list(list(name = "rate_loss", type = "scale",
                     target = "rate_loss", weight = 1, dim = 1L))
  tuning <- list(scale_rate_loss = 0.5)
  for (i in seq_len(2000L)) {
    tuning <- MkPrime:::.AdaptTuning(
      tuning, c(rate_loss = 0), c(rate_loss = 100), moves
    )
  }
  expect_gte(tuning$scale_rate_loss, 0.01)
})
