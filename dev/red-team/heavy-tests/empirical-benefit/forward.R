# Forward models for the empirical-vs-geometric k' prior study (#179).
#
# The empirical truth is T-EG-C's validated forward; the geometric truth
# mirrors it. Both draw k' from the inference prior truncated to [2, K], and
# discard an invariant (k', y) pair whole, as coding = "variable" conditions
# on the joint Pr(K = k, y | V). test-forward.R checks both against the
# inference prior. Needs MkPrime loaded and the working directory at the
# repository root.

source("dev/red-team/heavy-tests/kprime-viability/T-EG-C-forward.R")

# k' ~ 2 + Geometric(p) truncated to [2, K], as kPrimePrior = "geometric".
.DrawKGeo <- function(n, p, K) {
  kTrue <- integer(0)
  while (length(kTrue) < n) {
    k <- 2L + stats::rgeom(n, p)
    kTrue <- c(kTrue, k[k <= K])
  }
  # Return:
  kTrue[seq_len(n)]
}

.DrawK <- function(n, truth, p, emp, K) {
  switch(truth,
    empirical = .DrawKTrue(n, p, emp, K),
    geometric = .DrawKGeo(n, p, K),
    stop("Unknown truth: ", truth)
  )
}

# nChar recorded characters on `tree`, rows in tip-number order. kTrue[j] and
# kObs[j] belong to column j: score against these, never a re-read order (#54).
.SimDataset <- function(tree, nChar, truth, p, emp, K, maxTries = 1e5L) {
  nTip <- length(tree$tip.label)
  states <- matrix(NA_integer_, nTip, nChar,
                   dimnames = list(tree$tip.label, NULL))
  kTrue <- kObs <- integer(nChar)
  for (j in seq_len(nChar)) {
    for (i in seq_len(maxTries)) {
      k <- .DrawK(1L, truth, p, emp, K)
      cv <- .Canon(.SimJCChar(tree, k))
      if (attr(cv, "kObs") >= 2L) break
    }
    if (attr(cv, "kObs") < 2L) stop("No variable character in ", maxTries)
    states[, j] <- cv
    kTrue[[j]] <- k
    kObs[[j]] <- attr(cv, "kObs")
  }
  # Return:
  list(states = states, kTrue = kTrue, kObs = kObs)
}
