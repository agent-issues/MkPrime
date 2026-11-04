#!/usr/bin/env Rscript
# Checks the #179 study's forward models and scorer against the model MkPrime
# infers under. The empirical forward is T-EG-C's, gated by
# test-t-eg-c-forward.R; this adds the geometric one, and checks the exact
# per-k' likelihood against MkpLogLikelihood().
#
# Run from the repository root (about 30 s):
#   Rscript dev/red-team/heavy-tests/empirical-benefit/test-forward.R

.Ok <- function(label, cond) {
  cat(sprintf("%-66s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
source("dev/red-team/heavy-tests/empirical-benefit/forward.R")
source("dev/red-team/heavy-tests/empirical-benefit/score.R")

emp <- MkPrime::empiricalNObs
K <- 30L
nDraw <- 2e5L

.InferencePmf <- function(arm, p) {
  model <- list(kPrimePrior = arm, kprimeTruncK = K, empiricalNObs = emp)
  # Return:
  exp(.Mkp(".KprimeLogPriorTable")(model, p, 1, 1))
}

# Chi-square p-value of k' draws against a pmf on 2..K, pooling cells that
# expect fewer than 5 draws.
.ChisqP <- function(k, pmf) {
  observed <- tabulate(k, nbins = K)[2:K]
  expected <- length(k) * pmf
  small <- expected < 5
  if (any(small)) {
    observed <- c(observed[!small], sum(observed[small]))
    expected <- c(expected[!small], sum(expected[small]))
  }
  stat <- sum((observed - expected)^2 / expected)
  stats::pchisq(stat, df = length(expected) - 1L, lower.tail = FALSE)
}

set.seed(17901)
for (p in c(0.05, 0.3, 0.9)) {
  k <- .DrawKGeo(nDraw, p, K)
  .Ok(sprintf("geometric, p = %.2f: forward pmf matches the prior", p),
      .ChisqP(k, .InferencePmf("geometric", p)) > 1e-3)
}
.Ok("negative control: untruncated pmin() clamp fails",
    .ChisqP(pmin(2L + stats::rgeom(nDraw, 0.05), K),
            .InferencePmf("geometric", 0.05)) < 1e-10)

# Saturated 2-tip tree: tip states are iid uniform, so Pr(V | k) = 1 - 1 / k
# and recorded k' follow pi(k) (1 - 1 / k) / D.
saturated <- .SimTree(2L, shape = 2, expSteps = 50)
saturated$edge.length[] <- 1e3
for (truth in c("geometric", "empirical")) {
  arm <- c(geometric = "geometric", empirical = "empirical_geometric")[[truth]]
  recorded <- .SimDataset(saturated, 2e4L, truth, 0.5, emp, K)$kTrue
  joint <- .InferencePmf(arm, 0.5) * (1 - 1 / (2:K))
  .Ok(sprintf("%s: recorded k' follow the joint Pr(K, y | V)", truth),
      .ChisqP(recorded, joint / sum(joint)) > 1e-3)
}

tree <- .SimTree(12L, shape = 2, expSteps = 4)
d <- .SimDataset(tree, 20L, "empirical", 0.5, emp, K)
mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(d$states))
.Ok("MkPrimeData reads the simulated kObs in column order",
    all(mkd$kObs == d$kObs))
logLik <- .PerCharLogLik(tree, d$states, d$kObs, K)
for (arm in c("geometric", "empirical_geometric")) {
  kp <- d$kObs + sample(0:4, 20L, replace = TRUE)
  model <- MkPrimeModel(kPrimePrior = arm, kprimeTruncK = K, nCat = 1L)
  ref <- MkpLogLikelihood(tree, mkd, kPrime = kp, model = model, p = 0.4,
                          nCat = 1L)
  logPrior <- log(.InferencePmf(arm, 0.4))
  logD <- log(sum(exp(logPrior + .LogPVariable(tree, K))))
  .Ok(sprintf("%s: per-k' likelihood matches MkpLogLikelihood()", arm),
      abs(sum(logLik[cbind(1:20, kp - 1L)]) - 20 * logD - ref) < 1e-8)
}
