# Casali-like end-to-end smoke test for Layer 1 partition API.
#
# Exercises all five treatments from the Casali (2023) production pipeline:
#   T0  — unpartitioned (partition = NULL)
#   T1  — 2-class anatomical (unlink = c("shape", "ratemultiplier"))
#   T2a — AutoPart isolated: 3-class with an isolated singleton body class
#   T2b — AutoPart merged:   2-class, same as T1 skeleton but merged
#   T4  — random control: 3-class random assignment
#
# Additionally asserts T3 (unlink = "brlens") errors with the Layer-2 message.
#
# All chains: 200 iter, 50 warmup, nRuns = 1, nChains = 1, autoTune = FALSE.
# Runtime is ~10-20s total; skip_slow_tests() is NOT called here because
# individual chains are short.  Add the env guard if runtime grows.
#
# Invariants checked per treatment:
#   - chain completes without error
#   - all log_likelihood values are finite
#   - T0: no per-class columns (§7a bit-compat / legacy schema)
#   - T1/T2a/T2b/T4 with "shape":         class<c>_rate_log_sd > 0 every row
#   - T1/T2a/T2b/T4 with "ratemultiplier": w_<c> columns sum to 1 (tol 1e-6)

# ---- shared helpers ----------------------------------------------------------

.casali_data <- function(seed = 99L, nChar = 20L, nTip = 10L) {
  set.seed(seed)
  mat <- matrix(
    sample(0:1, nTip * nChar, replace = TRUE),
    nrow = nTip, ncol = nChar,
    dimnames = list(paste0("t", seq_len(nTip)), NULL)
  )
  # Ensure every character is variable (Casali parity: no invariant chars)
  for (j in seq_len(nChar)) {
    if (length(unique(mat[, j])) < 2L) mat[1L, j] <- 1L - mat[1L, j]
  }
  pd  <- MatrixToPhyDat(mat)
  mkd <- MkPrimeData(pd)
  tree <- Preorder(
    NJTree(pd, edgeLengths = TRUE)
  )
  if (is.null(tree$edge.length) || any(tree$edge.length <= 0)) {
    tree$edge.length <- pmax(
      tree$edge.length %||% rep(0.1, nrow(tree$edge)), 1e-8
    )
  }
  list(mkd = mkd, tree = tree, nChar = mkd$nChar)
}

.casali_mcmc <- function() {
  MkPrimeMCMC(
    nIter            = 200L,
    maxWarmup        = 50L,
    minWarmup        = 50L,
    nChains          = 1L,
    thin             = 1L,
    autoTune         = FALSE,
    gibbsSubtreeSwap = FALSE
  )
}

# Run one treatment and return result invisibly.  Asserts chain health inline.
.run_treatment <- function(d, part, unlink, label) {
  set.seed(123L)
  result <- allow_warning(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = .casali_mcmc(),
      partition = part,
      unlink    = unlink
    ),
    "without stabilisation"
  )

  samp      <- result$samples
  nClasses  <- if (is.null(part)) 0L else max(part)
  allUnlink <- unlink

  # (a) Chain completed, samples present
  expect_gt(result$nSamples, 0L,
    label = paste0(label, ": nSamples > 0"))

  # (b) All LLs finite
  expect_true(all(is.finite(samp[, "log_likelihood"])),
    label = paste0(label, ": all log_likelihood finite"))

  # T0: no per-class columns
  if (is.null(part)) {
    perClassCols <- grep("^class[0-9]+_rate_log_sd$|^w_[0-9]+$",
                         colnames(samp), value = TRUE)
    expect_length(perClassCols, 0L)
    return(invisible(result))
  }

  # (c) shape columns: present and positive
  if ("shape" %in% allUnlink) {
    for (c in seq_len(nClasses)) {
      col <- paste0("class", c, "_rate_log_sd")
      expect_true(col %in% colnames(samp),
        label = paste0(label, ": ", col, " present"))
      expect_true(all(samp[, col] > 0),
        label = paste0(label, ": ", col, " positive"))
    }
  }

  # (d) ratemultiplier columns: present, in (0,1), simplex (row-sum ~1)
  if ("ratemultiplier" %in% allUnlink) {
    wCols <- paste0("w_", seq_len(nClasses))
    for (col in wCols) {
      expect_true(col %in% colnames(samp),
        label = paste0(label, ": ", col, " present"))
      expect_true(all(samp[, col] > 0 & samp[, col] < 1),
        label = paste0(label, ": ", col, " in (0,1)"))
    }
    wMat <- samp[, wCols, drop = FALSE]
    rowSums_ <- rowSums(wMat)
    expect_true(all(abs(rowSums_ - 1.0) < 1e-6),
      label = paste0(label, ": w columns sum to 1 (simplex)"))
  }

  invisible(result)
}

# ---- fixture -----------------------------------------------------------------

d <- .casali_data()
n <- d$nChar

# Partition vectors
part_t1  <- c(rep(1L, floor(n / 2)), rep(2L, n - floor(n / 2)))
# T2a: isolated class 3 gets just the last character
part_t2a <- c(rep(1L, floor(n / 2)), rep(2L, n - floor(n / 2) - 1L), 3L)
part_t2b <- part_t1     # merged: same 2-class layout as T1
# T4: random 3-class assignment
set.seed(7L)
part_t4  <- sample(1:3, n, replace = TRUE)
# ensure all three classes are non-empty
while (length(unique(part_t4)) < 3L) {
  part_t4 <- sample(1:3, n, replace = TRUE)
}

unlink_both <- c("shape", "ratemultiplier")

# ---- T0: unpartitioned -------------------------------------------------------

test_that("T0 (partition=NULL) runs without per-class columns", {
  .run_treatment(d, part = NULL, unlink = character(0), label = "T0")
})

# ---- T1: anatomical 2-class --------------------------------------------------

test_that("T1 (anatomical 2-class) runs, per-class columns correct", {
  .run_treatment(d, part = part_t1, unlink = unlink_both, label = "T1")
})

# ---- T2a: AutoPart isolated 3-class ------------------------------------------

test_that("T2a (isolated 3-class) runs, per-class columns correct", {
  .run_treatment(d, part = part_t2a, unlink = unlink_both, label = "T2a")
})

# ---- T2b: AutoPart merged 2-class --------------------------------------------

test_that("T2b (merged 2-class) runs, per-class columns correct", {
  .run_treatment(d, part = part_t2b, unlink = unlink_both, label = "T2b")
})

# ---- T4: random 3-class control ----------------------------------------------

test_that("T4 (random 3-class) runs, per-class columns correct", {
  .run_treatment(d, part = part_t4, unlink = unlink_both, label = "T4")
})

# ---- T3: brlens deferred to Layer 2 ------------------------------------------

test_that("T3 (unlink=brlens) errors with Layer-2 message", {
  expect_error(
    RunMkPrime(
      data      = d$mkd,
      tree      = d$tree,
      mcmc      = .casali_mcmc(),
      partition = part_t1,
      unlink    = "brlens"
    ),
    regexp = "brlens|Layer 2",
    ignore.case = TRUE
  )
})
