# The chain's log-likelihood scored afresh, sharing no cached partials or
# joint ascertainment table with the state.
.ColdLogLik <- function(dataPtr, statePtr) {
  st <- get_mcmc_state(statePtr)
  cpp_log_likelihood_xptr(
    dataPtr, st$edge[, 1], st$edge[, 2], st$treeLength * st$relBrLengths,
    st$kPrime, st$rateLoss, st$rateLogSd, st$rateNeo, st$p, st$kprimeAlpha,
    st$kprimeBeta, st$betaScale
  )
}
