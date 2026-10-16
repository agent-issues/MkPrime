#!/usr/bin/env Rscript
# Checks that the T-EG-C forward matches the model MkPrime infers under:
# k' from empirical_geometric truncated to [2, K] and renormalized (#57, #392),
# and each recorded character from the joint Pr(K = k, y | V) (57628cf).
# A pmin() clamp, and redrawing y with k' held fixed, must fail the gates.
#
# Run from the repository root (about 20 s):
#   Rscript dev/red-team/heavy-tests/test-t-eg-c-forward.R

.Ok <- function(label, cond) {
  cat(sprintf("%-66s %s\n", label, if (isTRUE(cond)) "PASS" else "*** FAIL"))
  if (!isTRUE(cond)) stop(label)
}

suppressPackageStartupMessages({
  source("dev/red-team/heavy-tests/load-mkprime.R")
  LoadMkPrime(".")
})
source("dev/red-team/heavy-tests/kprime-viability/T-EG-C-forward.R")

emp <- MkPrime::empiricalNObs
K <- 30L
nDraw <- 2e5L

.InferencePmf <- function(p) {
  tab <- MkPrime:::.LogEmpGeomCapped(emp, p, K)
  exp(tab$logPk[2:K] - tab$logZ[2])
}

# Chi-square p-value of k' draws against a pmf on 2..K; cells expecting
# fewer than 5 draws are pooled so the asymptotic null holds.
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

.Clamped <- function(n, p) {
  pmin(MkPrime:::.SampleNObsEmpirical(n, emp) + stats::rgeom(n, p), K)
}

set.seed(5701)
for (p in c(0.02, 0.1, 0.5)) {
  k <- .DrawKTrue(nDraw, p, emp, K)
  .Ok(sprintf("p = %.2f: support is [2, K]", p), min(k) >= 2L && max(k) <= K)
  .Ok(
    sprintf("p = %.2f: forward pmf matches the inference prior", p),
    .ChisqP(k, .InferencePmf(p)) > 1e-3
  )
}
.Ok(
  "negative control: pmin() clamp fails the prior gate",
  .ChisqP(.Clamped(nDraw, 0.02), .InferencePmf(0.02)) < 1e-10
)

# On a 2-tip tree with saturated edges the tip states are iid uniform, so
# Pr(V | k) = 1 - 1 / k exactly and recorded k' follow pi(k) (1 - 1 / k) / D.
saturated <- .SimTree(2L, shape = 2, expSteps = 50)
saturated$edge.length[] <- 1e3
p <- 0.5
nChar <- 2e4L
recorded <- vapply(seq_len(nChar), function(i) {
  .SimVariableChar(saturated, p, emp, K)$k
}, integer(1))
joint <- .InferencePmf(p) * (1 - 1 / (2:K))
joint <- joint / sum(joint)
.Ok("recorded k' follow the joint Pr(K, y | V)", .ChisqP(recorded, joint) > 1e-3)
.Ok(
  "negative control: k' held fixed (per-k' Lewis form) fails",
  .ChisqP(.DrawKTrue(nChar, p, emp, K), joint) < 1e-10
)
