# Moving the likelihood kernels to inst/include/MkPrime must not change any
# result in the last bit. legacy/pruning-legacy.cpp is a frozen copy of the
# code they replaced; compiled with the same toolchain, it must agree
# exactly on a fixed dataset.

skip_on_cran()
# covr builds the package at -O0; sourceCpp() builds the oracle at R's default
# flags. On FMA hardware (arm64) only optimised builds contract a * b + c, so
# the two differ in the last bit for reasons unrelated to the source.
skip_on_covr()
legacy <- new.env()
Rcpp::sourceCpp(test_path("legacy", "pruning-legacy.cpp"), env = legacy,
                cacheDir = file.path(tempdir(), "mkp-legacy"))

.Sun2018States <- function() {
  nexFile <- system.file("datasets/Sun2018.nex", package = "TreeSearch")
  skip_if(nexFile == "", message = "TreeSearch not available")
  tokens <- PhyDatToMatrix(ReadAsPhyDat(nexFile))
  states <- suppressWarnings(matrix(as.integer(tokens), nrow(tokens),
                                    dimnames = dimnames(tokens)))
  states[is.na(states)] <- -1L
  states
}

test_that("refactored pruning is bit-identical to the legacy code", {
  states <- .Sun2018States()
  colMax <- apply(states, 2, max)
  trans <- states[, colMax >= 1 & colMax <= 2, drop = FALSE]
  binary <- states[, colMax == 1, drop = FALSE]
  expect_gt(ncol(trans), 20)
  expect_gt(ncol(binary), 10)

  set.seed(5014)
  unrooted <- Preorder(RandomTree(rownames(states), root = FALSE))
  unrooted <- Preorder(RenumberTips(unrooted, rownames(states)))
  unrooted$edge.length <- rexp(nrow(unrooted$edge), 8)
  rooted <- Preorder(RootTree(unrooted, 1))
  rooted$edge.length <- rexp(nrow(rooted$edge), 8)
  rates <- header_lognormal_rates(0.6, qnorm((seq_len(4) - 0.5) / 4))
  expect_identical(rates, legacy$legacy_acrv_rates(0.6, 4L,
                                                   qnorm((1:4 - 0.5) / 4)))

  for (tr in list(unrooted, rooted)) {
    p <- tr$edge[, 1]
    ch <- tr$edge[, 2]
    el <- tr$edge.length
    expect_identical(
      pruning_jc(p, ch, el, trans, 3L, rep(1 / 3, 3)),
      legacy$legacy_pruning_jc(p, ch, el, trans, 3L, rep(1 / 3, 3))
    )
    expect_identical(
      pruning_jc(p, ch, el, trans, 3L, c(0.5, 0.3, 0.2)),
      legacy$legacy_pruning_jc(p, ch, el, trans, 3L, c(0.5, 0.3, 0.2))
    )
    expect_identical(
      pruning_jc_collapsed(p, ch, el, trans, 7L, 3L),
      legacy$legacy_pruning_jc_collapsed(p, ch, el, trans, 7L, 3L)
    )
    expect_identical(
      pruning_mkn(p, ch, el, binary, 0.55, mkn_stationary_freqs(0.55)),
      legacy$legacy_pruning_mkn(p, ch, el, binary, 0.55,
                                mkn_stationary_freqs(0.55))
    )
    expect_identical(
      pruning_jc_acrv(p, ch, el, trans, 3L, rep(1 / 3, 3), rates),
      legacy$legacy_pruning_jc_acrv(p, ch, el, trans, 3L, rep(1 / 3, 3), rates)
    )
    expect_identical(
      pruning_jc_acrv_collapsed(p, ch, el, trans, 5L, 3L, rates),
      legacy$legacy_pruning_jc_acrv_collapsed(p, ch, el, trans, 5L, 3L, rates)
    )
    expect_identical(
      pruning_mkn_acrv(p, ch, el, binary, 1.7, mkn_stationary_freqs(1.7),
                       rates),
      legacy$legacy_pruning_mkn_acrv(p, ch, el, binary, 1.7,
                                     mkn_stationary_freqs(1.7), rates)
    )
  }
})

test_that("kernel structs reproduce every replaced inline expression", {
  tt <- c(0, 2 ^ -40, 1e-9, 1e-4, 0.013, 0.25, 1, 3.5, 60, 800)
  for (k in c(2L, 3L, 4L, 9L, 31L)) {
    expect_identical(header_jc_probs(k, tt), legacy$legacy_jc_probs(k, tt))
  }
  for (r in c(0.05, 0.6, 1, 2.3, 40)) {
    expect_identical(header_mkn_probs(r, tt), legacy$legacy_mkn_probs(r, tt))
  }
})
