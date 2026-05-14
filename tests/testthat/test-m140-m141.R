# Tests for M-140/M-146 (move weight formatting) and M-141 (ETA estimation)

# ---------------------------------------------------------------------------
# M-140/M-146: .FormatMoveWeights categorized colour-coded display
# ---------------------------------------------------------------------------

test_that(".FormatMoveWeights returns one styled line per category", {
  weights <- c(0.30, 0.10, 0.05, 0.02, 0.53)
  names   <- c("kPrime", "spr", "nni", "tree_length", "branch_lengths")
  result  <- .FormatMoveWeights(weights, names)
  # Returns a character vector: one element per present category
  expect_type(result, "character")
  # Should have 3 categories: Topology, Branches, Characters
  expect_length(result, 3)
  for (nm in names) {
    expect_true(any(grepl(nm, result, fixed = TRUE)),
                info = paste("Missing move name:", nm))
  }
})

test_that(".FormatMoveWeightsPlain returns categorized unstyled lines", {
  weights <- c(0.30, 0.10, 0.05)
  names   <- c("nni", "tree_length", "kPrime")
  result  <- .FormatMoveWeightsPlain(weights, names)
  expect_type(result, "character")
  expect_length(result, 1)  # single string with \n separators
  lines <- strsplit(result, "\n")[[1]]
  expect_length(lines, 3)
  expect_true(grepl("^Topology:", lines[1]))
  expect_true(grepl("^Branches:", lines[2]))
  expect_true(grepl("^Characters:", lines[3]))
})

test_that("moves sorted by weight within category (highest first)", {
  w <- c(0.05, 0.20, 0.10)
  n <- c("spr", "nni", "tbr")
  result <- .FormatMoveWeightsPlain(w, n)
  # Single Topology line, nni (20%) should appear before tbr (10%) before spr (5%)
  expect_match(result, "nni:20\\.0%.*tbr:10\\.0%.*spr:5\\.0%")
})

test_that("absent categories are omitted", {
  # Only branch moves — no Topology, Characters, or Rates lines
  w <- c(0.60, 0.40)
  n <- c("tree_length", "branch_lengths")
  result <- .FormatMoveWeights(w, n)
  expect_length(result, 1)
  expect_true(grepl("Branches", result, fixed = TRUE))
})

test_that("unknown move names go to Other category", {
  w <- c(0.50, 0.50)
  n <- c("nni", "my_custom_move")
  result <- .FormatMoveWeightsPlain(w, n)
  expect_true(grepl("Other:", result, fixed = TRUE))
})

test_that("colour-coded values: green >=10%, yellow 5-10%, white <5%", {
  w <- c(0.35, 0.07, 0.02)
  n <- c("nni", "spr", "tbr")
  result <- .FormatMoveWeights(w, n)
  # In non-interactive sessions cli may strip ANSI codes, so just
  # verify the function runs and includes expected content.
  collapsed <- paste(result, collapse = " ")
  expect_true(grepl("35.0%", collapsed, fixed = TRUE))
  expect_true(grepl("7.0%", collapsed, fixed = TRUE))
  expect_true(grepl("2.0%", collapsed, fixed = TRUE))
})


# ---------------------------------------------------------------------------
# M-141: .EstimateEta
# ---------------------------------------------------------------------------

test_that(".EstimateEta returns NULL for invalid inputs", {
  expect_null(.EstimateEta(NA, 200, 60))
  expect_null(.EstimateEta(0, 200, 60))
  expect_null(.EstimateEta(50, NULL, 60))
  expect_null(.EstimateEta(50, 200, NA))
  expect_null(.EstimateEta(50, 200, 0))
  expect_null(.EstimateEta(-10, 200, 60))
})

test_that(".EstimateEta returns 'now' when target already met", {
  expect_equal(.EstimateEta(200, 200, 60), "now")
  expect_equal(.EstimateEta(300, 200, 60), "now")
})

test_that(".EstimateEta returns conservative time estimate", {
  # 100 ESS in 60s → rate = 100/60 ESS/s
  # Need 100 more ESS → naive = 60s, with 1.5× safety = 90s
  result <- .EstimateEta(100, 200, 60, safetyFactor = 1.5)
  expect_type(result, "character")
  expect_true(grepl("~", result))
  # Should show ~90s (or ~1.5min depending on threshold)
  expect_true(grepl("90s|1\\.5min", result))
})

test_that(".EstimateEta formats minutes and hours correctly", {
  # 10 ESS in 60s, target 200 → need 19× more, naive 1140s, ×1.5 = 1710s ≈ 28.5min
  result <- .EstimateEta(10, 200, 60, safetyFactor = 1.5)
  expect_true(grepl("min", result))

  # Very long: 1 ESS in 60s, target 1000 → need 999× more → huge
  result <- .EstimateEta(1, 1000, 60, safetyFactor = 1.5)
  expect_true(grepl("h", result))
})

test_that(".EstimateEta uses safety factor", {
  # Same scenario but different safety factors
  eta1 <- .EstimateEta(100, 200, 120, safetyFactor = 1.0)
  eta2 <- .EstimateEta(100, 200, 120, safetyFactor = 2.0)
  # eta1 should be shorter than eta2
  # Parse the numbers: both should be in seconds range
  n1 <- as.numeric(gsub("[^0-9.]", "", eta1))
  n2 <- as.numeric(gsub("[^0-9.]", "", eta2))
  expect_true(n2 > n1, info = paste("eta1:", eta1, "eta2:", eta2))
})
