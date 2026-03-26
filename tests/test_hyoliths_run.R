pkgbuild::compile_dll(debug = FALSE, quiet = TRUE)
devtools::load_all(quiet = TRUE)

nex_file <- system.file("extdata", "hyoliths.nex", package = "MkPrime")
pd       <- TreeTools::ReadAsPhyDat(nex_file)
cat("Loaded:", sum(attr(pd, "weight")), "chars,", length(pd), "taxa\n")

# Parse character types from STATELABELS
nex_lines <- readLines(nex_file)
sl_start  <- grep("^\\s+STATELABELS", nex_lines) + 1L
mx_line   <- grep("^\\s+MATRIX", nex_lines)
sl_end    <- mx_line[mx_line > sl_start][1L] - 1L
sl_block  <- nex_lines[sl_start:sl_end]
hdrs      <- grep("^\\s+\\d+\\s*$", sl_block)
cnums     <- as.integer(trimws(sl_block[hdrs]))
ranges    <- Map(function(s, e) sl_block[s:e], hdrs + 1L, c(hdrs[-1L] - 1L, length(sl_block)))
is_trans  <- vapply(ranges, function(l) any(grepl("Transformational character", l, ignore.case = TRUE)), logical(1))
neo_chars <- cnums[!is_trans]
cat("neomorphic:", length(neo_chars), "| transformational:", sum(is_trans), "\n")

# Starting tree
start_tree <- TreeTools::NJTree(pd, edgeLengths = TRUE)
start_tree <- TreeTools::UnrootTree(start_tree)
cat("Tree:", ape::Ntip(start_tree), "tips, edge.length NULL:", is.null(start_tree$edge.length), "\n")

# Short MCMC run
set.seed(8203)
result <- RunMkPrime(
  data       = pd,
  tree       = start_tree,
  neomorphic = neo_chars,
  model      = MkPrimeModel(coding = "variable", nCat = 6L, relabel = TRUE),
  mcmc       = MkPrimeMCMC(nIter  = 2000L, thin = 10L, warmup = 1000L,
                            nRuns  = 2L,   nChains = 1L)
)
cat("Full-data test: OK -", nrow(result$samples), "samples\n")
print(summary(result))
