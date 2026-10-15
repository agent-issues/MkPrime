# #433: McmcData / McmcState must own plain storage.  `.InitMcmcData()`
# passes tip states through `storage.mode<-`, which returns an ALTREP
# wrapper for vectors of more than 64 elements; Rcpp caches a pointer into
# the wrapper's payload, which a later writable DATAPTR() can replace and
# the GC can then free (use-after-free in build_cl_groups() under valgrind).

test_that("MCMC data and state hold no ALTREP storage", {
  set.seed(3847L)
  tree <- Preorder(ape::rtree(8L, rooted = FALSE))
  mat <- matrix(sample(0:2, 8L * 12L, replace = TRUE), nrow = 8L,
                dimnames = list(tree$tip.label, paste0("c", 1:12)))
  mkd <- MkPrimeData(MatrixToPhyDat(mat))
  model <- MkPrime:::.FinalizeModel(MkPrimeModel(qHeterogeneity = TRUE),
                                    tree, mkd)

  # The R-side input really is a wrapper (R makes one above 64 elements),
  # so this test exercises the hazard.
  sizes <- vapply(mkd$partitions, function(p) length(p$tip_states), 1)
  ts <- mkd$partitions[[which.max(sizes)]]$tip_states
  expect_gt(length(ts), 64L)
  ts[is.na(ts)] <- -1L                 # as .InitMcmcData() does
  storage.mode(ts) <- "integer"
  expect_match(utils::capture.output(.Internal(inspect(ts)))[1], "wrapper")

  dataPtr  <- MkPrime:::.InitMcmcData(mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(MkPrime:::.InitState(tree, mkd, model))
  altrep <- mcmc_altrep_storage(dataPtr, statePtr)
  expect_false(any(altrep), label = paste(names(altrep)[altrep], collapse = ", "))
})
