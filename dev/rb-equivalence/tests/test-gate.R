# Tests for dev/rb-equivalence/R/gate.R (agent-issues/MkPrime#216), on
# synthetic chains. Needs only testthat, posterior and coda.
#
# Usage:
#   Rscript -e 'testthat::test_file("dev/rb-equivalence/tests/test-gate.R")'

suppressPackageStartupMessages(library(testthat))

source(file.path("..", "R", "gate.R"))

# compare.R's gate before #216, kept to show what the new gate catches.
OldGate <- function(mkChains, rbChains, targetRhat = 1.025, targetEss = 128) {
  mat <- do.call(cbind, c(mkChains, rbChains))
  rhat <- posterior::rhat_basic(mat)
  ess <- as.numeric(coda::effectiveSize(unlist(c(mkChains, rbChains))))
  rhat < targetRhat && ess > targetEss
}

OldDropBurnin <- function(x, frac = 0.25) {
  n <- length(x)
  start <- max(2L, floor(n * frac))
  x[start:n]
}

IidRuns <- function(nRuns, n, mu = 0) lapply(seq_len(nRuns), function(i) rnorm(n, mu))

test_that("DropBurnin drops exactly the first quarter (RB-111)", {
  expect_length(OldDropBurnin(1:8), 7L)
  expect_identical(DropBurnin(1:8), 3:8)
  expect_identical(DropBurnin(1:100), 26:100)
  df <- data.frame(a = 1:8)
  expect_identical(DropBurnin(df)$a, 3:8)
  expect_identical(DropBurnin(1:3), 1:3)
})

test_that("a 0.4-SD cross-sampler shift fails the gate (RB-110)", {
  set.seed(216)
  verdicts <- oldPass <- logical(0)
  for (i in 1:20) {
    mk <- IidRuns(2, 500)
    rb <- IidRuns(2, 500, mu = 0.4)
    oldPass[i] <- OldGate(mk, rb)
    verdicts[i] <- GateParam(mk, rb, 1.025, 128, nTests = 5L)$verdict == "FAIL"
  }
  expect_gt(mean(oldPass), 0.5)
  expect_true(all(verdicts))
})

test_that("identical samplers pass", {
  set.seed(1)
  g <- GateParam(IidRuns(2, 2000), IidRuns(2, 2000), 1.025, 128, nTests = 5L)
  expect_identical(g$verdict, "PASS")
  expect_identical(g$reasons, "")
  expect_lt(g$mdd_sd, 0.25)
})

test_that("MkPrime runs that disagree with each other fail, whatever RB does", {
  # pid 950's CID: MkPrime-only R-hat 1.23 hidden by RB's large pooled ESS.
  set.seed(950)
  mk <- list(rnorm(2000, 0), rnorm(2000, 0.5))
  rb <- IidRuns(2, 20000, mu = 0.25)
  expect_true(OldGate(mk, rb, targetRhat = 1.1))
  g <- GateParam(mk, rb, 1.025, 128)
  expect_identical(g$verdict, "FAIL")
  expect_match(g$reasons, "mkprime_unconverged")
})

test_that("low per-sampler ESS fails even when pooled ESS clears the target", {
  set.seed(3)
  ar <- function(n) as.numeric(arima.sim(list(ar = 0.995), n)) * sqrt(1 - 0.995^2)
  mk <- list(ar(2000), ar(2000))
  rb <- IidRuns(2, 5000)
  g <- GateParam(mk, rb, 1.025, 128)
  expect_lt(g$mkp_ess, 128)
  expect_identical(g$verdict, "FAIL")
  expect_match(g$reasons, "mkprime_unconverged")
})

test_that("a test too weak to see the tolerance is UNDERPOWERED, not PASS", {
  set.seed(4)
  g <- GateParam(IidRuns(2, 100), IidRuns(2, 100), 1.1, 64, nTests = 5L)
  expect_gt(g$mdd_sd, 0.25)
  expect_identical(g$verdict, "UNDERPOWERED")
})

test_that("log scale sees a location shift in a skewed parameter", {
  set.seed(5)
  mk <- lapply(1:2, function(i) rlnorm(2000, 0, 2))
  rb <- lapply(1:2, function(i) rlnorm(2000, 0.5, 2))  # 0.25 SD on log scale
  expect_identical(
    GateParam(mk, rb, 1.025, 128, logScale = TRUE, nTests = 5L)$verdict, "FAIL"
  )
  expect_error(MeanShiftTest(list(c(1, 0, 2, 3, 4)), list(1:5), logScale = TRUE),
               "Non-positive")
})

test_that("a missing scalar is an error, not a smaller comparison (RB-109)", {
  full <- data.frame(tree_length = 1, rate_log_sd = 1, rate_loss = 1, rate_neo = 1)
  expect_true(AssertScalarsPresent(list(full, full), "rb"))
  expect_error(AssertScalarsPresent(list(full, full[, -4]), "rb"),
               "rb run 2 lacks expected scalar\\(s\\): rate_neo")
})

test_that("MkPrime runs must have converged at the gate's own targets (RB-226)", {
  ok <- list(stop_reason = "converged", mcmc = list(maxRhat = 1.025, minEss = 128))
  expect_true(AssertMkRunTargets(ok, 1.025, 128))
  timedOut <- modifyList(ok, list(stop_reason = "max_time"))
  expect_error(AssertMkRunTargets(timedOut, 1.025, 128), "max_time")
  loose <- modifyList(ok, list(mcmc = list(maxRhat = 1.1, minEss = 32)))
  expect_error(AssertMkRunTargets(loose, 1.025, 128), "looser")
  expect_error(AssertMkRunTargets(list(), 1.025, 128), "NA")
})

RenderTemplate <- function(name) {
  txt <- readLines(file.path("..", "templates", name))
  txt <- gsub("\\$\\{NCHAR_NEO\\}", "10", txt)
  gsub("\\$\\{NCHAR_TRANS\\}", "20", txt)
}

rbMatchedSpec <- list(
  coding = "variable", nCat = 6L,
  treeLengthShape = 2, treeLengthRate = 2,
  rateLogSdShape = 1, rateLogSdRate = 1,
  rateLossMeanlog = 0, rateLossSdlog = 2,
  rateNeoMeanlog = 0, rateNeoSdlog = 2
)

test_that("RevModelSpec reads the settings RevBayes actually runs (RB-109)", {
  for (tmpl in c("by_nt_9v.template.Rev", "by_nt_kv.template.Rev")) {
    expect_identical(RevModelSpec(RenderTemplate(tmpl)), rbMatchedSpec)
  }
})

test_that("RevModelSpec reports a changed or unparseable template", {
  rev <- RenderTemplate("by_nt_9v.template.Rev")
  changed <- sub("rate_neo ~ dnLognormal(mean = 0, sd = 2)",
                 "rate_neo ~ dnLognormal(mean = 0, sd = 1)", rev, fixed = TRUE)
  expect_identical(RevModelSpec(changed)$rateNeoSdlog, 1)
  expect_false(isTRUE(all.equal(RevModelSpec(changed), rbMatchedSpec)))
  expect_error(RevModelSpec(sub("dnGamma( 1, 1 )", "dnExponential(1)", rev,
                                fixed = TRUE)), "rate_log_sd")
  firstCoding <- grep('coding = "variable"', rev, fixed = TRUE)[[1]]
  rev[firstCoding] <- sub('"variable"', '"all"', rev[firstCoding], fixed = TRUE)
  expect_error(RevModelSpec(rev), "one coding")
})

test_that("MkModelSpec matches RevModelSpec for the harness's model", {
  expect_error(MkModelSpec(NULL), "MkPrimeModel")
  skip_if_not_installed("MkPrime")
  source(file.path("..", "R", "utils.R"))
  expect_equal(MkModelSpec(RBMatchedModel()), rbMatchedSpec)
})

test_that("AssertSameSetup refuses mismatched or unstamped rds", {
  mk <- list(model_spec = rbMatchedSpec, burnin_frac = 0.25)
  expect_true(AssertSameSetup(mk, mk))
  expect_error(AssertSameSetup(list(), mk), "model_spec missing")
  expect_error(AssertSameSetup(mk, modifyList(mk, list(burnin_frac = 0))),
               "burnin_frac")
  other <- modifyList(mk, list(model_spec = list(nCat = 4L)))
  expect_error(AssertSameSetup(mk, other), "Model settings differ")
})

test_that("AppendCsv tolerates new columns and replaces same-key rows", {
  path <- tempfile(fileext = ".csv")
  write.csv(data.frame(pid = c("1", "2"), model = "m", param = "a", rhat = 1),
            path, row.names = FALSE)
  AppendCsv(data.frame(pid = "1", model = "m", param = "a", rhat = 2,
                       verdict = "PASS"), path)
  out <- read.csv(path, stringsAsFactors = FALSE)
  expect_setequal(names(out), c("pid", "model", "param", "rhat", "verdict"))
  expect_identical(nrow(out), 2L)
  expect_equal(out$rhat[out$pid == 1], 2)
  expect_true(is.na(out$verdict[out$pid == 2]))
})
