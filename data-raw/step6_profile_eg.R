# Step 6: wall-time profile of moves under empirical_geometric.
#
# Runs a short eg MCMC on a representative tree-inference replicate
# (mkprime/tree-inference/tree_01/rep_01 — 50 chars, n_tips per dataset),
# instrumented with R's Rprof at the R level and Linux `perf` at the
# C++ level if available.
#
# Output: data-raw/step6-profile.Rprof.out + summary table.

suppressPackageStartupMessages({
  devtools::load_all(quiet = TRUE)
  library(TreeTools)
})

set.seed(2026)

DATA_ROOT <- "C:/Users/pjjg18/GitHub/mkprime/tree-inference"
tree_idx  <- 1L
rep_idx   <- 1L

# Read all chr*.nex files for tree_01/rep_01.
rep_dir <- file.path(DATA_ROOT, sprintf("tree_%02d", tree_idx),
                     sprintf("rep_%02d", rep_idx))
chr_files <- list.files(rep_dir, pattern = "\\.nex$", full.names = TRUE)
cat(sprintf("Loading %d characters from %s\n", length(chr_files), rep_dir))

# Combine into a single phyDat.
mats <- lapply(chr_files, function(f) {
  ch <- TreeTools::ReadCharacters(f)
  if (is.matrix(ch)) ch else matrix(ch, ncol = 1L,
                                     dimnames = list(names(ch), NULL))
})
# Align by tip-label union (use first as template).
tips <- rownames(mats[[1]])
mat <- do.call(cbind, lapply(mats, function(m) m[tips, , drop = FALSE]))
mode(mat) <- "character"
cat(sprintf("Combined matrix: %d tips x %d chars\n", nrow(mat), ncol(mat)))

pd <- TreeTools::MatrixToPhyDat(mat)
mkd <- MkPrimeData(pd)
cat(sprintf("Variable chars: %d  |  kObs distribution: %s\n",
            ncol(mkd$charMatrix),
            paste(names(table(mkd$kObs)), table(mkd$kObs), sep = ":",
                  collapse = " ")))

# NJ start tree.
start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)

# Short MCMC run, eg arm.
nIter <- 4000L
profPath <- "data-raw/step6-profile.Rprof.out"

cat(sprintf("\nProfiling %d iterations under empirical_geometric...\n", nIter))
t0 <- Sys.time()
Rprof(profPath, interval = 0.02, line.profiling = TRUE,
      memory.profiling = FALSE)
res <- RunMkPrime(
  mkd, start_tree,
  model = MkPrimeModel(coding = "variable",
                        kPrimePrior = "empirical_geometric"),
  mcmc = MkPrimeMCMC(nIter = nIter, thin = 50L,
                     maxWarmup = 1000L, minWarmup = 1000L,
                     autoTune = FALSE,
                     nRuns = 1L, nChains = 2L,
                     progressFn = function(...) invisible())
)
Rprof(NULL)
cat(sprintf("Wall time: %s\n", format(Sys.time() - t0, digits = 3)))

# Summarise.
cat("\n--- Per-move acceptance ---\n")
print(res$acceptance)

cat("\n--- Rprof self time (top 20) ---\n")
prof <- summaryRprof(profPath)
print(head(prof$by.self, 20))

cat("\n--- Rprof total time by function (top 20) ---\n")
print(head(prof$by.total, 20))

# Group moves by name keyword for readability.
cat("\n--- Sampling time share by move family ---\n")
top <- prof$by.self
move_keys <- c(
  "gibbs_kprime"   = "Cpp_RunBatch|gibbs_kprime",
  "spr/topology"   = "Cpp_RunBatch.*spr|nni|tbr",
  "branch_lengths" = "branch_length|dirichlet|tree_length|scale_tl",
  "rates"          = "rate_log|rate_loss|rate_neo",
  "p"              = "scale_p|gibbs_p|logit_scale",
  "kPrime_walk"    = "int_walk",
  "block_kPrime"   = "block_kprime"
)
cat("(Rprof only sees R callers; native moves appear as Cpp_RunBatch self time)\n")
cat(sprintf("Cpp_RunBatch self time: %.2fs  (this is the main MCMC engine)\n",
            sum(top[grepl("Cpp_RunBatch|Mcmc", rownames(top)), "self.time"])))

saveRDS(list(res = res, profile = prof, wall_time = Sys.time() - t0),
        "data-raw/step6-profile-result.rds")
cat("\nSaved data-raw/step6-profile-result.rds\n")
