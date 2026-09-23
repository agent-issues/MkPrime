# Mixed-partition RB conformance for marginal_k (issue #25).
#
# test-marginal-k-rb-free-topology.R proves the RB identity
#   lse_{k'} [logPrior_S(k') + logLik_S(k')] == logPrior_M + logLik_M
# across random topologies -- but on a transformational-only fixture, where
# compute_partition_scales() returns (1, 1) and the partition rate scale is
# invisible. Every other marginal_k gate on main shares that blind spot.
#
# This fixture carries one neomorphic and two transformational characters of
# different kObs (three partitions), so neoScale/transScale are both != 1 and
# the marginal evaluator must apply transScale to transformational edges exactly
# as cpp_partition_log_likelihood does under sampled_k. The scales depend on the
# nNeo:nTrans ratio as well as rate_neo, so they are non-unit even at
# rate_neo == 1; rate_neo is varied anyway to exercise the full range.
#
# No known-k character: marginal_k rejects those outright
# (test-marginal-k-unsupported.R).


.rbmp_lse <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

# Three-character fixture: c1 neomorphic (binary), c2 transformational kObs 2,
# c3 transformational kObs 3. Drawn until every character is variable and c3
# shows all three states, so MkPrimeData drops nothing and the indices hold.
.rbmp_data <- function(tipLabels) {
  n <- length(tipLabels)
  repeat {
    m <- cbind(sample(c("0", "1"), n, replace = TRUE),
               sample(c("0", "1"), n, replace = TRUE),
               sample(c("0", "1", "2"), n, replace = TRUE))
    nStates <- apply(m, 2, function(x) length(unique(x)))
    if (identical(nStates, c(2L, 2L, 3L))) break
  }
  dimnames(m) <- list(tipLabels, c("c1", "c2", "c3"))
  MkPrimeData(MatrixToPhyDat(m), neomorphic = 1L)
}

# One evaluator closure per (mode, tree, data): build model + data once, then
# vary only the state across calls.
.rbmp_eval_closure <- function(mode, trp, mkd, K) {
  mod <- suppressMessages(MkPrimeModel(
    coding = "variable", nCat = 4L, kPrimePrior = "geometric",
    likelihoodMode = mode, priorVariant = "unconditional",
    kprimeTruncK = K, kprimeHyperA = 1, kprimeHyperB = 1, expSteps = 1.4))
  modf <- MkPrime:::.FinalizeModel(mod, trp, mkd)
  dp <- MkPrime:::.InitMcmcData(mkd, modf)
  TL <- sum(trp$edge.length)
  RBL <- trp$edge.length / TL
  function(mutate) {
    st <- MkPrime:::.InitState(trp, mkd, modf)
    st$tree_length <- TL
    st$rel_br_lengths <- RBL
    st <- mutate(st)
    sp <- MkPrime:::.InitMcmcChain(st)
    fill_partition_cache(dp, sp)
    eval_log_prior_cpp(dp, sp) + eval_full_loglik_cpp(dp, sp)
  }
}

test_that("mixed neo/trans fixture has three partitions, non-unit scales", {
  set.seed(4242L)
  mkd <- .rbmp_data(paste0("t", 1:8))
  expect_equal(unname(mkd$type), c("neomorphic", "transformational",
                                   "transformational"))
  expect_equal(unname(mkd$kObs), c(2L, 2L, 3L))
  expect_equal(length(mkd$partitions), 3L)
  # transScale = 1/(1 + rate_neo) * nTotal/nTrans; != 1 for every rate_neo here
  # because nNeo (1) != nTrans (2). This is what makes the fixture diagnostic.
  transScale <- function(rn) 1 / (1 + rn) * 3 / 2
  expect_false(isTRUE(all.equal(transScale(1), 1)))
  expect_false(isTRUE(all.equal(transScale(2.5), 1)))
})

test_that("marginal-k RB identity holds on a mixed neo/trans partition (#25)", {
  set.seed(1801L)
  K <- 10L
  for (rep in seq_len(3L)) {
    rn <- c(0.4, 2.5, 5)[rep]
    ntip <- sample(7:9, 1L)
    tr <- ape::rtree(ntip, tip.label = paste0("t", seq_len(ntip)))
    tr$edge.length <- runif(nrow(tr$edge), 0.03, 0.5)
    trp <- Preorder(tr)
    mkd <- .rbmp_data(tr$tip.label)
    transIdx <- which(mkd$type == "transformational")
    ks <- lapply(transIdx, function(i) mkd$kObs[i]:K)

    evalM <- .rbmp_eval_closure("marginal_k", trp, mkd, K)
    evalS <- .rbmp_eval_closure("sampled_k", trp, mkd, K)

    for (p in c(0.15, 0.45)) for (rls in c(0, 0.6)) {
      M <- evalM(function(s) {
        s$p <- p
        s$rate_log_sd <- rls
        s$rate_neo <- rn
        s
      })
      grid <- expand.grid(ks[[1]], ks[[2]])
      J <- vapply(seq_len(nrow(grid)), function(g) evalS(function(s) {
        s$p <- p
        s$rate_log_sd <- rls
        s$rate_neo <- rn
        s$kPrime[transIdx] <- as.integer(unlist(grid[g, ]))
        s
      }), numeric(1))
      tag <- sprintf("rep%d ntip=%d rn=%.1f p=%.2f rls=%.1f",
                     rep, ntip, rn, p, rls)
      # Non-vacuity: the joint genuinely varies over the k' grid.
      expect_gt(stats::sd(J), 1e-6,
                label = paste(tag, ": joint varies over k'"))
      # MAIN: summing the sampled_k joint over the k' grid reproduces the
      # marginal_k joint. Fails by whole nats when the marginal evaluator prunes
      # at unscaled transformational edge lengths.
      expect_equal(.rbmp_lse(J), M, tolerance = 1e-7,
                   label = paste(tag, ": lse(sampled joint) == marginal joint"))
    }
  }
})
