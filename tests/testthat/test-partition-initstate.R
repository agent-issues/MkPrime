# Tests for .InitStatePartitioned() — the partition-aware MCMC state
# initializer added in Layer 1. Trivial spec must equal legacy .InitState
# to ~1e-10 on log_lik; non-trivial spec must produce per-class fields
# with the correct shape and the §7b numeric-equivalence property at
# initial state.



.setup_init <- function(partition = NULL, unlink = character(0),
                        nChar = 8L, nTip = 6L, seed = 17L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, nTip * nChar, replace = TRUE),
                nrow = nTip, ncol = nChar,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- NJTree(pd, edgeLengths = TRUE)
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- pmax(tree$edge.length %||% rep(0.1, nrow(tree$edge)), 1e-8)
  }
  tree <- Preorder(tree)

  # Apply partition to mkd$partitions for the partition-aware path
  if (!is.null(partition)) {
    mkd$partitions <- .BuildPartitions(mkd, partition = partition)
  }
  spec <- .ValidatePartitionArgs(partition, unlink, mkd)

  model <- MkPrimeModel()
  model <- .FinalizeModel(model, tree, mkd)

  list(tree = tree, mkd = mkd, model = model, spec = spec)
}


# ---- trivial spec collapses to legacy initial log_lik (§7b at t=0) ----

test_that("trivial partition: log_lik equals legacy .InitState to ~1e-10", {
  s <- .setup_init(partition = NULL)  # spec = trivial; spec$partition NULL
  legacy <- .InitState(s$tree, s$mkd, s$model)
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)

  expect_equal(partitioned$log_lik, legacy$log_lik, tolerance = 1e-10)
  # All-1 class structure on the trivial spec
  expect_identical(partitioned$nChar_c, as.integer(s$mkd$nChar))
  expect_identical(partitioned$class_w, 1.0)
  expect_identical(partitioned$class_rate, 1.0)
  expect_identical(partitioned$class_rate_log_sd, legacy$rate_log_sd)
  expect_identical(partitioned$eta_neo, 1.0)
})


test_that("partition = rep(1L, nChar) with no unlink collapses to legacy", {
  s <- .setup_init(partition = rep(1L, 8L), unlink = character(0))
  legacy <- .InitState(s$tree, s$mkd, s$model)
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)
  expect_equal(partitioned$log_lik, legacy$log_lik, tolerance = 1e-10)
})


# ---- non-trivial partition with class_rate ≡ 1 matches legacy ----

test_that("multi-class partition starts at class_rate ≡ 1 (§7b at t=0)", {
  s <- .setup_init(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                   unlink    = c("shape", "ratemultiplier"))
  legacy <- .InitState(s$tree, s$mkd, s$model)
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)

  # nChar_c: tabulation of the partition vector
  expect_identical(partitioned$nChar_c, c(4L, 4L))
  # class_w on the unit simplex, initialised to nChar_c / nChar
  expect_equal(sum(partitioned$class_w), 1)
  expect_equal(partitioned$class_w, c(0.5, 0.5), tolerance = 1e-12)
  # class_rate = 1 by construction at initialisation
  expect_equal(partitioned$class_rate, c(1, 1), tolerance = 1e-12)
  # class_rate_log_sd is length nClasses because "shape" is in unlink
  expect_identical(length(partitioned$class_rate_log_sd), 2L)
  expect_true(all(partitioned$class_rate_log_sd == legacy$rate_log_sd))

  # §7b numeric-equivalence: initial log_lik matches legacy
  expect_equal(partitioned$log_lik, legacy$log_lik, tolerance = 1e-10)
})


test_that("ratemultiplier linked: class_rate is still derivable but w is uniform", {
  # Only "shape" in unlink; ratemultiplier is linked. class_rate is initialised
  # to ≡ 1 regardless of unlink choice (the linked case just means subsequent
  # moves won't touch w).
  s <- .setup_init(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                   unlink    = "shape")
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)
  expect_equal(partitioned$class_rate, c(1, 1), tolerance = 1e-12)
  expect_identical(length(partitioned$class_rate_log_sd), 2L)
})


test_that("shape linked: class_rate_log_sd is length 1", {
  s <- .setup_init(partition = c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                   unlink    = "ratemultiplier")
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)
  expect_identical(length(partitioned$class_rate_log_sd), 1L)
  expect_identical(partitioned$class_rate_log_sd, .InitState(s$tree, s$mkd, s$model)$rate_log_sd)
})


# ---- partition with imbalanced class sizes ----

test_that("imbalanced classes: class_w follows nChar_c / nChar; class_rate ≡ 1", {
  # 6 chars in class 1, 2 chars in class 2
  s <- .setup_init(partition = c(rep(1L, 6L), rep(2L, 2L)),
                   unlink    = c("shape", "ratemultiplier"))
  partitioned <- .InitStatePartitioned(s$tree, s$mkd, s$model, s$spec)

  expect_identical(partitioned$nChar_c, c(6L, 2L))
  expect_equal(partitioned$class_w, c(6/8, 2/8), tolerance = 1e-12)
  expect_equal(partitioned$class_rate, c(1, 1), tolerance = 1e-12)
  # mean-1 constraint: sum(nChar_c * class_rate) / nChar == 1
  expect_equal(sum(partitioned$nChar_c * partitioned$class_rate) / s$mkd$nChar,
               1, tolerance = 1e-12)
})


# ---- helpers: w ↔ class_rate round-trip ----

test_that("w → class_rate → w round-trips", {
  nCharPC <- c(4L, 6L, 10L)
  w <- c(0.2, 0.3, 0.5)
  cr <- .PartitionWToClassRate(w, nCharPC)
  w_back <- .PartitionClassRateToW(cr, nCharPC)
  expect_equal(w_back, w, tolerance = 1e-15)
  # mean-1 constraint
  expect_equal(sum(nCharPC * cr) / sum(nCharPC), 1, tolerance = 1e-12)
})


test_that(".PartitionNCharPerClass errors on empty classes", {
  expect_error(.PartitionNCharPerClass(c(1L, 1L, 3L, 3L), nClasses = 3L),
               regexp = "Empty user class",
               fixed = FALSE)
})


test_that(".PartitionNCharPerClass returns correct counts", {
  expect_identical(.PartitionNCharPerClass(c(1L, 1L, 2L, 2L, 2L), nClasses = 2L),
                   c(2L, 3L))
})
