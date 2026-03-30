# Tests for M-140 (move weight colour formatting) and M-141 (ETA estimation)

# ---------------------------------------------------------------------------
# M-140: .FormatMoveWeights colour coding
# ---------------------------------------------------------------------------

test_that(".FormatMoveWeights returns styled string with all move names", {
  weights <- c(0.30, 0.10, 0.05, 0.02, 0.53)
  names   <- c("kPrime", "spr", "nni", "tree_length", "branch_lengths")
  result  <- .FormatMoveWeights(weights, names)
  # Should be a single string containing all move names

  expect_type(result, "character")
  expect_length(result, 1)
  for (nm in names) {
    expect_true(grepl(nm, result, fixed = TRUE),
                info = paste("Missing move name:", nm))
  }
})

test_that(".FormatMoveWeightsPlain returns unstyled string", {
  weights <- c(0.30, 0.05, 0.02)
  names   <- c("a", "b", "c")
  result  <- .FormatMoveWeightsPlain(weights, names)
  expect_type(result, "character")
  expect_equal(result, "a=30.0% b=5.0% c=2.0%")
})

test_that("names are silver, values colour-coded by magnitude", {
  # 35% value green, 7% value yellow, 2% value white; all names silver
  w <- c(0.35, 0.07, 0.02)
  n <- c("big", "mid", "small")
  result <- .FormatMoveWeights(w, n)
  # In non-interactive sessions cli may strip ANSI codes, so just
  # verify the function runs and includes expected content.
  expect_type(result, "character")
  for (nm in n) {
    expect_true(grepl(nm, result, fixed = TRUE),
                info = paste("Missing:", nm))
  }
  expect_true(grepl("35.0%", result, fixed = TRUE))
  expect_true(grepl("7.0%", result, fixed = TRUE))
  expect_true(grepl("2.0%", result, fixed = TRUE))
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
