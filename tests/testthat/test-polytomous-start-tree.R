# #336: every tree move assumes a binary tree. A polytomous start tree used to
# reach the moves and segfault, hang or corrupt the chain.

PolyData <- function(tree) {
  set.seed(336L)
  tips <- tree$tip.label
  mat <- matrix(sample(0:2, length(tips) * 12L, TRUE), length(tips), 12L,
                dimnames = list(tips, NULL))
  MkPrimeData(TreeTools::MatrixToPhyDat(mat))
}

ChildCounts <- function(tree) {
  root <- setdiff(tree$edge[, 1], tree$edge[, 2])
  nChild <- table(tree$edge[, 1])
  list(root = unname(nChild[as.character(root)]),
       other = unname(nChild[names(nChild) != as.character(root)]))
}

rootPoly <- ape::read.tree(text =
  "((a:0.1,b:0.2):0.1,(c:0.1,d:0.1):0.2,(e:0.1,f:0.1):0.1,(g:0.2,h:0.1):0.1);")
innerPoly <- ape::read.tree(text =
  "((a:0.1,b:0.2,c:0.15):0.1,(d:0.1,e:0.1):0.2,((f:0.1,g:0.1):0.1,(h:0.2,i:0.1):0.1):0.1);")

test_that(".PrepareStartTree resolves polytomies with one warning", {
  for (tree in list(rootPoly, innerPoly)) {
    nTip <- length(tree$tip.label)
    expect_warning(
      out <- .PrepareStartTree(tree, PolyData(tree)),
      "polytom"
    )
    expect_true(ape::is.binary(out))
    expect_equal(nrow(out$edge), 2L * nTip - 2L)
    expect_equal(ChildCounts(out)$root, 2L)
    expect_true(all(ChildCounts(out)$other == 2L))
    expect_true(all(out$edge.length > 0))
    # Inserted edges are near-zero; every original edge keeps its length.
    inserted <- out$edge.length == 1e-8
    expect_equal(sum(inserted), nrow(out$edge) - nrow(tree$edge))
    expect_equal(sort(out$edge.length[!inserted]), sort(tree$edge.length))
  }
})

test_that(".PrepareStartTree keeps both binary shapes unchanged", {
  rooted <- ape::read.tree(text = "(((a:1,b:1):1,c:1):1,(d:1,e:1):1);")
  unrooted <- ape::read.tree(text = "((a:1,b:1):1,c:1,(d:1,e:1):1);")
  expect_true(ape::is.binary(rooted))
  expect_true(ape::is.binary(unrooted))
  for (tree in list(rooted, unrooted)) {
    expect_no_warning(out <- .PrepareStartTree(tree, PolyData(tree)))
    expect_equal(nrow(out$edge), nrow(tree$edge))
  }
})

test_that(".PrepareStartTree rejects single-child nodes", {
  tree <- ape::read.tree(text = "(((a:1,b:1):1):1,(c:1,d:1):1);")
  expect_error(.PrepareStartTree(tree, PolyData(tree)), "single child")
})

test_that("init_mcmc_state rejects a non-binary edge matrix", {
  for (tree in list(rootPoly, innerPoly)) {
    tree <- TreeTools::Preorder(tree)
    expect_error(
      init_mcmc_state(tree$edge[, 1], tree$edge[, 2],
                      tree$edge.length / sum(tree$edge.length),
                      sum(tree$edge.length), rateLoss = 1, rateLogSd = 0.5,
                      rateNeo = 1, p = 0.5, kPrime = rep(3L, 12L),
                      logLik = 0, logPrior = 0),
      "binary"
    )
  }
})

test_that("RunMkPrime completes from a polytomous start tree", {
  tree <- innerPoly
  set.seed(1L)
  warnings <- capture_warnings(
    res <- RunMkPrime(PolyData(tree), tree, verbosity = 0,
                      mcmc = MkPrimeMCMC(nRuns = 1L, nChains = 1L,
                                         nIter = 200L, thin = 2L,
                                         minWarmup = 20L, maxWarmup = 20L,
                                         autoTune = FALSE, maxTime = 30))
  )
  expect_match(warnings, "polytom", all = FALSE)
  nChild <- unique(unlist(lapply(res$trees, function(t) table(t$edge[, 1]))))
  expect_setequal(nChild, 2L)
})

test_that("mkp_stepping_stone resolves a polytomous start tree", {
  tree <- rootPoly
  set.seed(2L)
  expect_warning(
    ss <- mkp_stepping_stone(PolyData(tree), tree, nStones = 3L, nIter = 30L,
                             warmup = 10L, nRuns = 1L, verbose = FALSE),
    "polytom"
  )
  expect_true(is.finite(ss$log_marginal))
})
