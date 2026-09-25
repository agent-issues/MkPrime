# Unit tests for dev/rb-equivalence/R/utils.R, run standalone (not part of
# the package's own tests/testthat/ suite -- these are dev-harness tests on
# small synthetic inputs, no RevBayes or ../neotrans checkout required).
#
# Usage (from anywhere, with MkPrime + TreeTools installed -- testthat::test_file()
# changes the working directory to the test file's own directory while it runs):
#   Rscript -e '.libPaths(c("<lib>", .libPaths()));
#               testthat::test_file("dev/rb-equivalence/tests/test-utils.R")'

suppressPackageStartupMessages({
  library(testthat)
  library(MkPrime)
  library(TreeTools)
})

source(file.path("..", "R", "utils.R"))

strip_symbols <- function(mat) {
  attr(mat, "nexusSymbols") <- NULL
  unname(mat)
}

write_nex <- function(path, ntax, nchar, symbols, rows) {
  writeLines(c(
    "#NEXUS",
    "BEGIN DATA;",
    sprintf("  DIMENSIONS NTAX=%d NCHAR=%d;", ntax, nchar),
    sprintf(
      "  FORMAT DATATYPE=STANDARD MISSING=? GAP=- SYMBOLS=\"%s\" INTERLEAVE=NO;",
      symbols
    ),
    "  MATRIX",
    rows,
    "  ;",
    "END;"
  ), path)
}

# --- RbStateK: RevBayes' max-index+1 rule, not distinct-token count --------

test_that("RbStateK matches RevBayes' max-observed-index+1 rule", {
  # column coded {1,2}: kObs (distinct tokens) = 2, RB k = 3 (agent-issues/MkPrime#214)
  expect_identical(RbStateK(c("1", "2")), 3L)
  expect_identical(RbStateK(c("0", "1")), 2L)
  expect_identical(RbStateK(c("?", "-", NA)), 0L)
  expect_identical(RbStateK(c("0", "0", "0")), 1L)
})

test_that("RbStateK counts polymorphic bits individually", {
  # {01} contributes states 0 and 1, same as if they appeared unbundled
  expect_identical(RbStateK(c("2", "{01}")), 3L)
  expect_identical(RbStateK(c("0", "(12)")), 3L)
})

test_that("ComputeRbK is RbStateK applied per column", {
  mat <- matrix(c("1", "2", "1", "2",
                  "0", "1", "0", "1",
                  "0", "0", "0", "0"),
                nrow = 4, ncol = 3)
  expect_identical(ComputeRbK(mat), c(3L, 2L, 1L))
})

# --- PrepareCellInfo: shared k/nChar/taxa derivation ------------------------

test_that("PrepareCellInfo derives RB-matched k, common taxa and post-drop nChar", {
  transPath <- tempfile(fileext = ".nex")
  neoPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(transPath, neoPath)))

  # trans: col1 in {1,2} (kObs=2, kRb=3); col2 invariant after intersection
  # (taxa A-E all "0" once F is dropped) -- but F only appears in neo, so
  # this column here uses all of A:E = c(2,1,0,0,0), non-invariant.
  write_nex(transPath, ntax = 5, nchar = 2, symbols = "012", rows = c(
    "    A  12", "    B  21", "    C  10", "    D  20", "    E  ?0"
  ))
  write_nex(neoPath, ntax = 6, nchar = 1, symbols = "01", rows = c(
    "    A  0", "    B  1", "    C  0", "    D  1", "    E  0", "    F  1"
  ))

  transMat <- ReadSplit(transPath)
  neoMat <- ReadSplit(neoPath)

  info <- suppressWarnings(PrepareCellInfo(transMat, neoMat, "by_nt_kv"))

  expect_identical(info$commonTaxa, c("A", "B", "C", "D", "E"))
  expect_identical(info$excludeNeoTaxa, "F")
  expect_identical(info$excludeTransTaxa, character(0))
  expect_identical(info$nTrans, 2L)
  expect_identical(info$nNeo, 1L)
  # col1: 1,2,1,2,? -> max index 2 -> k=3. col2: 2,1,0,0,0 -> max index 2 -> k=3
  expect_identical(info$kRb, c(3L, 3L))
  expect_identical(info$nNeoFinal, 1L)
  expect_identical(info$nTransFinal, 2L)
})

test_that("PrepareCellInfo rejects a cell RevBayes would silently truncate (k > 10)", {
  transPath <- tempfile(fileext = ".nex")
  neoPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(transPath, neoPath)))

  write_nex(transPath, ntax = 2, nchar = 1, symbols = "0123456789A", rows = c(
    "    A  0", "    B  A"  # state index 10 -> k = 11
  ))
  write_nex(neoPath, ntax = 2, nchar = 1, symbols = "01", rows = c(
    "    A  0", "    B  1"
  ))
  transMat <- ReadSplit(transPath)
  neoMat <- ReadSplit(neoPath)

  expect_error(PrepareCellInfo(transMat, neoMat, "by_nt_kv"), "k > 10")
})

# --- WriteOrCheckCellInfo: cross-process equality assertion -----------------

test_that("WriteOrCheckCellInfo writes once then asserts equality", {
  path <- tempfile()
  on.exit(unlink(path))
  info <- list(a = 1L, b = "x")

  expect_false(file.exists(path))
  WriteOrCheckCellInfo(path, info)
  expect_true(file.exists(path))

  expect_silent(WriteOrCheckCellInfo(path, info))

  mismatched <- info
  mismatched$a <- 2L
  expect_error(WriteOrCheckCellInfo(path, mismatched), "does not match")
})

# --- WriteFilteredNexus: round-trips through TreeTools::ReadCharacters -----

test_that("WriteFilteredNexus round-trips tokens including polymorphisms/gaps/missing", {
  srcPath <- tempfile(fileext = ".nex")
  outPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(srcPath, outPath)))

  write_nex(srcPath, ntax = 4, nchar = 3, symbols = "012", rows = c(
    "    A  1{01}2", "    B  210", "    C  10?", "    D  2(12)-"
  ))
  orig <- ReadSplit(srcPath)
  WriteFilteredNexus(orig, outPath)
  back <- ReadSplit(outPath)

  expect_identical(strip_symbols(orig), strip_symbols(back))
  expect_identical(rownames(orig), rownames(back))
})

test_that("RbStateK indexes by position in the declared SYMBOLS string, not sorted order", {
  # A file declaring SYMBOLS="ABC" (A=0, B=1, C=2): column {A,C} -> k = 3
  expect_identical(RbStateK(c("A", "C"), alphabet = c("A", "B", "C")), 3L)
  # the same tokens under a differently-ordered declaration give a different k
  expect_identical(RbStateK(c("A", "C"), alphabet = c("A", "C", "B")), 2L)
})

test_that("PrepareCellInfo uses the source file's own SYMBOLS declaration for kRb", {
  transPath <- tempfile(fileext = ".nex")
  neoPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(transPath, neoPath)))

  # Non-default SYMBOLS order: state "1" is declared second (index 1), "2"
  # third (index 2) -- same as the default numeric reading here, but this
  # would differ from a naive digit-value or sorted-order assumption.
  write_nex(transPath, ntax = 2, nchar = 1, symbols = "021", rows = c(
    "    A  1", "    B  2"
  ))
  write_nex(neoPath, ntax = 2, nchar = 1, symbols = "01", rows = c(
    "    A  0", "    B  1"
  ))
  transMat <- ReadSplit(transPath)
  neoMat <- ReadSplit(neoPath)
  info <- PrepareCellInfo(transMat, neoMat, "by_nt_kv")
  # symbols "021": "1" at position 2 (k=3), so max observed index for a
  # column containing "1" and "2" (positions 2, 1) is 2 -> kRb = 3
  expect_identical(info$kRb, 3L)
})

test_that("PrepareCellInfo sanitises polymorphic tokens to '?' for both matrices", {
  transPath <- tempfile(fileext = ".nex")
  neoPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(transPath, neoPath)))

  write_nex(transPath, ntax = 3, nchar = 1, symbols = "012", rows = c(
    "    A  1", "    B  {12}", "    C  0"
  ))
  write_nex(neoPath, ntax = 3, nchar = 1, symbols = "01", rows = c(
    "    A  0", "    B  1", "    C  0"
  ))
  transMat <- ReadSplit(transPath)
  neoMat <- ReadSplit(neoPath)
  info <- PrepareCellInfo(transMat, neoMat, "by_nt_kv")

  expect_identical(unname(info$transMatCanonical[, 1]), c("1", "?", "0"))
})

test_that("WriteFilteredNexus preserves a taxon subset and row order", {
  srcPath <- tempfile(fileext = ".nex")
  outPath <- tempfile(fileext = ".nex")
  on.exit(unlink(c(srcPath, outPath)))

  write_nex(srcPath, ntax = 3, nchar = 1, symbols = "01", rows = c(
    "    A  0", "    B  1", "    C  0"
  ))
  orig <- ReadSplit(srcPath)
  subset <- orig[c("C", "A"), , drop = FALSE]
  WriteFilteredNexus(subset, outPath)
  back <- ReadSplit(outPath)

  expect_identical(rownames(back), c("C", "A"))
  expect_identical(strip_symbols(back), strip_symbols(subset))
})
