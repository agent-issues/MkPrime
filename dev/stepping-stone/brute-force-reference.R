# Reference log marginal likelihood for the brute-force test in
# tests/testthat/test-stepping-stone.R: simple Monte Carlo from the prior,
# Z = E_prior[L], on four tips, three binary characters of known k = 2,
# unrooted (three topologies, five edges). Four seeds of 5e4 draws gave
# -8.6273, -8.6233, -8.6327, -8.6288 (each MC se 0.0035): -8.628 +- 0.002.
library("MkPrime")

mat <- matrix(c(0L, 1L, 0L, 1L,
                0L, 0L, 1L, 1L,
                0L, 1L, 1L, 1L),
              nrow = 4, dimnames = list(paste0("t", 1:4), NULL))
mkd <- MkPrimeData(TreeTools::MatrixToPhyDat(mat),
                   knownStates = c("1" = 2L, "2" = 2L, "3" = 2L))
model <- MkPrime:::.FinalizeModel(MkPrimeModel(), NULL, mkd)
topologies <- lapply(c("((t1,t2),(t3,t4));", "((t1,t3),(t2,t4));",
                       "((t1,t4),(t2,t3));"),
                     function(text) ape::unroot(ape::read.tree(text = text)))

LogLikDraw <- function() {
  tree <- topologies[[sample.int(3L, 1L)]]
  treeLength <- rgamma(1, model$treeLengthShape,
                       MkPrime:::.TreeLengthRate(model))
  fraction <- rexp(nrow(tree$edge))
  tree$edge.length <- treeLength * fraction / sum(fraction)
  rateLogSd <- rgamma(1, model$rateLogSdShape, model$rateLogSdRate)
  MkpLogLikelihood(tree, mkd, rate_log_sd = rateLogSd, model = model)
}

set.seed(101)
logLik <- vapply(seq_len(5e4), function(i) LogLikDraw(), numeric(1))
w <- exp(logLik - max(logLik))
cat("log Z", max(logLik) + log(mean(w)),
    "MC se", sd(w) / mean(w) / sqrt(length(w)), "\n")
