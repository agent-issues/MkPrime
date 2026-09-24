# A two-class partition-API chain whose class state has left its initial
# point. There every evaluator agrees (class rates 1, one shared sigma); away
# from it, an evaluator that ignores the per-class shape or rate multiplier
# scores a different model from cpp_log_likelihood_partitioned.
#
# kPrime starts one above kObs so that block_kprime_shift can move either way.
# `hyperprior = TRUE` samples sigma_c = tau * z_c under the pooled hyperprior,
# starting from tau = 1.
.PartitionedChain <- function(tree, mkd, model = MkPrimeModel(),
                              classW = c(0.3, 0.7),
                              classRateLogSd = c(0.3, 0.9),
                              hyperprior = FALSE) {
  partition <- rep(1:2, length.out = mkd$nChar)
  mkd$partitions <- MkPrime:::.BuildPartitions(mkd, partition = partition)
  model <- MkPrime:::.FinalizeModel(model, tree, mkd)
  st <- MkPrime:::.InitState(tree, mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(mkd, model)

  nCharPerClass <- tabulate(partition, 2L)
  classRate <- classW * mkd$nChar / nCharPerClass
  NewState <- function(logPrior) {
    init_mcmc_state(
      st$tree$edge[, 1L], st$tree$edge[, 2L],
      st$rel_br_lengths, st$tree_length,
      st$rate_loss, classRateLogSd[[1L]],
      st$rate_neo %||% 1, st$p %||% 0.5,
      as.integer(st$kPrime) + 1L, 0, logPrior,
      st$beta_scale %||% 1,
      st$kprime_alpha %||% 1, st$kprime_beta %||% 1,
      classRateLogSd = classRateLogSd, classW = classW,
      classRate = classRate, nCharPerClass = nCharPerClass,
      useHyperpriorOnSigma = hyperprior, hyperTau = 1, classZ = classZ
    )
  }
  classZ <- if (hyperprior) classRateLogSd else numeric(0)
  statePtr <- NewState(eval_log_prior_partitioned_cpp(
    dataPtr, NewState(0), classRateLogSd, classW, 1,
    useHyperpriorOnSigma = hyperprior, hyperTau = 1, classZ = classZ))
  fill_partition_cache(dataPtr, statePtr)
  allocate_cl_workspace(dataPtr, statePtr)

  # Return:
  list(dataPtr = dataPtr, statePtr = statePtr, nChar = mkd$nChar)
}

# Cold partitioned log-likelihood of `edge` / `edgeLen` (default: the chain's
# current tree) under the chain's current parameters.
.PartitionedLogLik <- function(chain, edge = NULL, edgeLen = NULL) {
  st <- get_mcmc_state(chain$statePtr)
  if (is.null(edge)) {
    edge <- st$edge
    edgeLen <- st$treeLength * st$relBrLengths
  }
  # Return:
  cpp_log_likelihood_partitioned_xptr(
    chain$dataPtr, edge[, 1L], edge[, 2L], edgeLen, st$kPrime,
    st$rateLoss, st$classRateLogSd, st$classRate, st$etaNeo, st$betaScale)
}
