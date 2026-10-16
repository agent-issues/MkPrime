# Forward model for T-EG-C-modelA-sbc.R, kept apart so that
# test-t-eg-c-forward.R can check it against the inference prior.
# Needs MkPrime loaded (for .SampleNObsEmpirical).

# Tree length ~ Gamma(shape, shape / expSteps), split Dirichlet(1, ..., 1).
.SimTree <- function(nTip, shape, expSteps) {
  tr <- ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)))
  tl <- stats::rgamma(1, shape = shape, rate = shape / expSteps)
  g <- stats::rexp(nrow(tr$edge))
  tr$edge.length <- tl * g / sum(g)
  TreeTools::Preorder(tr)
}

# k' ~ empirical_geometric truncated to [2, K], as inference with
# kprimeTruncK = K (#392). The (N_obs, N_unobs) pair is redrawn together:
# redrawing only N_unobs would weight each N_obs by 1 / P(N_unobs <= K - N_obs).
.DrawKTrue <- function(n, p, emp, K) {
  kTrue <- integer(0)
  while (length(kTrue) < n) {
    k <- MkPrime:::.SampleNObsEmpirical(n, emp) + stats::rgeom(n, p)
    kTrue <- c(kTrue, k[k <= K])
  }
  # Return:
  kTrue[seq_len(n)]
}

# JC(k) root-to-tip on a preordered tree; rate convention as
# inst/include/MkPrime/transition.h (exp(-k t / (k - 1))).
.SimJCchar <- function(tree, k) {
  nTip <- length(tree$tip.label)
  states <- integer(2L * nTip - 1L)
  states[nTip + 1L] <- sample.int(k, 1L) - 1L
  edges <- tree$edge
  pSame <- 1 / k + (1 - 1 / k) * exp(-k * tree$edge.length / (k - 1))
  for (e in seq_len(nrow(edges))) {
    pa <- edges[e, 1L]
    ch <- edges[e, 2L]
    states[ch] <- if (stats::runif(1L) < pSame[e]) {
      states[pa]
    } else {
      (states[pa] + sample.int(k - 1L, 1L)) %% k
    }
  }
  # Return:
  states[seq_len(nTip)]
}

.Canon <- function(v) {
  uvals <- sort(unique(v))
  out <- match(v, uvals) - 1L
  attr(out, "kObs") <- length(uvals)
  out
}

# coding = "variable" conditions the data, not k', on variability: the
# likelihood is P(y | k') / P(variable | k'). So k' stays fixed and only y is
# redrawn; dropping invariant characters would tilt k' by P(variable | k').
# NULL when maxTries draws are all invariant.
.SimVariableChar <- function(tree, k, maxTries = 1e5L) {
  for (i in seq_len(maxTries)) {
    cv <- .Canon(.SimJCchar(tree, k))
    if (attr(cv, "kObs") >= 2L) {
      # Return:
      return(cv)
    }
  }
  NULL
}
