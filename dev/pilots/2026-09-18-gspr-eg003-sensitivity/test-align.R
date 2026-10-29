# Synthetic-fixture checks for the character-order and grid helpers.
# Run: Rscript dev/pilots/2026-09-18-gspr-eg003-sensitivity/test-align.R
source("dev/pilots/2026-09-18-gspr-eg003-sensitivity/align.R")

Check <- function(label, ok) {
  cat(sprintf("%s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", label))
  if (!isTRUE(ok)) quit(status = 1L)
}

files <- sprintf("chr%d.nex", c(1, 2, 3, 10, 11))
gt <- data.frame(char_idx = c(1, 2, 3, 10, 11), kObs = c(2, 3, 4, 5, 6),
                 u_true = c(0, 1, 2, 3, 4))

Check("CharIdxFromFiles parses numbers", identical(CharIdxFromFiles(files), c(1L, 2L, 3L, 10L, 11L)))

lex <- CharIdxFromFiles(sort(files))
Check("lexical order differs from numeric", !identical(lex, sort(lex)))
al <- AlignGroundTruth(gt, lex)
Check("aligned to lexical sampler order", identical(al$kObs, c(2, 5, 6, 3, 4)))
Check("numeric order is identity", identical(AlignGroundTruth(gt, gt$char_idx), gt))
Check("subset (invariant characters dropped) aligns",
      identical(AlignGroundTruth(gt, c(10L, 2L))$u_true, c(3, 1)))
Check("unknown char_idx errors", inherits(try(AlignGroundTruth(gt, 99L), silent = TRUE), "try-error"))

grid <- expand.grid(tree = 1:2, arm = c("on", "off"), seed = 1:3, stringsAsFactors = FALSE)
Check("complete grid returns seed count", CheckCompleteGrid(grid, 1:2) == 3L)
Check("dropped cell errors and names it",
      grepl("tree 2 off has 2",
            tryCatch(CheckCompleteGrid(grid[-which(grid$tree == 2 & grid$arm == "off")[1], ], 1:2),
                     error = conditionMessage)))
Check("wholly missing tree errors",
      inherits(try(CheckCompleteGrid(grid[grid$tree == 1, ], 1:2), silent = TRUE), "try-error"))

# End-to-end: 01 on a fabricated tree whose sampler used numeric order (post-#54,
# char_idx recorded) and with invariant characters dropped. Corrected pairing must
# be perfect: k_post = kObs + u_true exactly.
tmp <- tempfile(); dir.create(tmp)
gtRoot <- file.path(tmp, "gt"); sumDir <- file.path(tmp, "summary")
dir.create(sumDir)
idx <- c(1L, 2L, 3L, 10L, 11L)
gtFull <- data.frame(char_idx = idx, kObs = c(2L, 3L, 4L, 5L, 6L), u_true = c(0L, 1L, 2L, 3L, 4L))
for (t in 1:26) {
  d <- file.path(gtRoot, sprintf("tree_%02d/rep_01", t))
  dir.create(d, recursive = TRUE)
  write.csv(gtFull, file.path(d, "ground_truth.csv"), row.names = FALSE)
  file.create(file.path(d, sprintf("chr%d.nex", idx)))
}
keep <- c(3L, 1L, 10L)  # sampler order, two characters dropped as invariant
kPost <- gtFull$kObs[match(keep, idx)] + gtFull$u_true[match(keep, idx)]
saveRDS(list(tag = "t01_r01", kp_means = kPost, char_idx = keep),
        file.path(sumDir, "mkp_geo_t01_r01.rds"))
saveRDS(list(), file.path(tmp, "eg.rds"))
out <- suppressWarnings(system2("Rscript",
  "dev/pilots/2026-09-18-gspr-eg003-sensitivity/01-character-order-audit.R",
  stdout = TRUE, stderr = TRUE,
  env = c(paste0("AUDIT_SUMMARY_DIR=", sumDir), paste0("AUDIT_GT_ROOT=", gtRoot),
          paste0("AUDIT_EG_POST=", file.path(tmp, "eg.rds")))))
Check("01 runs on a fixture", is.null(attr(out, "status")))
Check("01 keeps the invariant-bearing task", any(grepl("prior: geo \\( 1 tasks", out)))
Check("01 counts characters actually present, not 50 * n", any(grepl("of 3 char-tasks", out)))
Check("01 corrected pairing is exact (rho = 1)", any(grepl("rho .*corrected \\+1\\.0000", out)))
Check("01 corrected k' < kObs violations are zero", any(grepl("corrected 0 ", out)))
