#!/usr/bin/env Rscript
# Checks that the T-EG-C forward draws k' from the prior the inference uses
# under kprimeTruncK = K: empirical_geometric truncated to [2, K] and
# renormalised (#57, #392). A pmin() clamp must fail the same gate.
#
# Run from the repository root (about 10 s):
#   Rscript dev/red-team/heavy-tests/test-t-eg-c-forward.R

ok <- function(label, cond) {
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

.ChisqP <- function(k, p) {
  counts <- tabulate(k, nbins = K)[2:K]
  expected <- length(k) * .InferencePmf(p)
  keep <- expected > 0
  stat <- sum((counts[keep] - expected[keep])^2 / expected[keep])
  stats::pchisq(stat, df = sum(keep) - 1L, lower.tail = FALSE)
}

.Clamped <- function(n, p) {
  pmin(MkPrime:::.SampleNObsEmpirical(n, emp) + stats::rgeom(n, p), K)
}

set.seed(5701)
for (p in c(0.02, 0.1, 0.5)) {
  k <- .DrawKTrue(nDraw, p, emp, K)
  ok(sprintf("p = %.2f: support is [2, K]", p), min(k) >= 2L && max(k) <= K)
  ok(
    sprintf("p = %.2f: forward pmf matches the inference prior", p),
    .ChisqP(k, p) > 1e-3
  )
}

p <- 0.02
k <- .DrawKTrue(nDraw, p, emp, K)
ok(
  "p = 0.02: no point mass at K (within 4 SE of the inference prior)",
  abs(mean(k == K) - .InferencePmf(p)[K - 1L]) <
    4 * sqrt(.InferencePmf(p)[K - 1L] / nDraw)
)
ok(
  "negative control: pmin() clamp fails the gate",
  .ChisqP(.Clamped(nDraw, p), p) < 1e-10
)

tr <- .SimTree(8L, shape = 2, expSteps = 0.05)
kObs <- vapply(seq_len(200L), function(i) {
  attr(.SimVariableChar(tr, 5L), "kObs")
}, integer(1))
ok("variable redraw: every character varies on a short tree", all(kObs >= 2L))
ok("variable redraw: kObs never exceeds k", all(kObs <= 5L))
