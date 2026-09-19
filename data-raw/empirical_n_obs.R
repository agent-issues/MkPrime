# Build empiricalNObs: empirical pmf for the number of observed states in
# transformational characters, tabulated from the neotrans project corpus.
#
# Output: data/empiricalNObs.rda (an "MkPrimeEmpiricalPrior" object)
#
# Re-run with:  source("data-raw/empirical_n_obs.R")

stopifnot(requireNamespace("MkPrime", quietly = TRUE))
stopifnot(requireNamespace("TreeTools", quietly = TRUE))
stopifnot(requireNamespace("usethis", quietly = TRUE))

neotransProjects <- "C:/Users/pjjg18/GitHub/neotrans/inst/projects"
if (!dir.exists(neotransProjects)) {
  stop("neotrans project directory not found: ", neotransProjects)
}

transFiles <- list.files(neotransProjects, pattern = "\\.trans\\.nex$",
                          full.names = TRUE)
message("Reading ", length(transFiles), " transformational matrices.")

# Count distinct observed states in one column.  Polymorphisms encoded as
# {01} or (0,1) etc.; ReadCharacters returns each cell as a single string
# such as "0", "1", "(0,1)", "?", or "-".  Split polymorphisms into their
# constituent state tokens before counting.
.KObs <- function(col) {
  tokens <- col[!is.na(col) & col != "?" & col != "-" & col != ""]
  if (length(tokens) == 0L) {
    # Return:
    return(0L)
  }
  # Split polymorphic cells like "(0,1)" or "{01}" into individual states.
  isPoly <- grepl("[(){},]", tokens)
  states <- character(0)
  if (any(!isPoly)) {
    states <- c(states, tokens[!isPoly])
  }
  if (any(isPoly)) {
    polyStates <- unlist(strsplit(
      gsub("[(){},]", " ", tokens[isPoly]),
      "\\s+"
    ))
    states <- c(states, polyStates[nzchar(polyStates)])
  }
  # Return:
  length(unique(states))
}

kObsCounts <- integer(0)
for (file in transFiles) {
  chars <- tryCatch(
    TreeTools::ReadCharacters(file),
    error = function(e) {
      warning("Failed to read ", basename(file), ": ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(chars)) next
  perChar <- apply(chars, 2, .KObs)
  kObsCounts <- c(kObsCounts, perChar)
}

# Drop invariant characters
kObsCounts <- kObsCounts[kObsCounts > 1L]
message("Total informative transformational characters: ", length(kObsCounts))

nMax <- max(kObsCounts)
rawCounts <- tabulate(kObsCounts, nbins = nMax)
names(rawCounts) <- seq_len(nMax)
message("Raw count table (kObs -> count):")
print(rawCounts[rawCounts > 0L])

# Fit a geometric tail to kObs >= tailThreshold.  Geometric MLE on the
# shifted variable (kObs - tailThreshold) gives decay rate q = mean / (1 + mean)
# under the Geom(q) = q (1-q)^x parameterisation, so 1 - q is the per-step
# decay multiplier we want when extending the pmf to k > nMax.
tailThreshold <- 4L
tailCounts <- kObsCounts[kObsCounts >= tailThreshold]
if (length(tailCounts) < 5L) {
  stop("Too few tail observations (", length(tailCounts),
       ") to fit geometric decay.")
}
tailMean <- mean(tailCounts - tailThreshold)
# Probability of incrementing (decay multiplier per step beyond nMax).
qDecay <- tailMean / (1 + tailMean)
message(sprintf("Tail (kObs >= %d): n = %d, decay = %.4f",
                tailThreshold, length(tailCounts), qDecay))

# Convert raw counts to a probability mass function supported on k >= 2.
# `MkPrimeEmpiricalPrior()` takes the body unnormalised: it anchors a
# geometric tail one past the last body entry and rescales both to sum to 1.
bodyCounts <- rawCounts[seq.int(2L, nMax)]
totalChar <- sum(rawCounts)
bodyProps <- bodyCounts / totalChar

empiricalNObs <- MkPrime::MkPrimeEmpiricalPrior(
  body = bodyProps,
  tail_decay = qDecay,
  nSource = totalChar
)

stopifnot(abs(sum(empiricalNObs$body) +
                empiricalNObs$tail_start_p / (1 - empiricalNObs$tail_decay) -
                1) < 1e-12)

message("Body probabilities (kObs -> P):")
print(round(empiricalNObs$body, 5))
message(sprintf("Tail starts at k = %d with mass %.5f and decay %.4f",
                empiricalNObs$tail_start_k,
                empiricalNObs$tail_start_p,
                empiricalNObs$tail_decay))

usethis::use_data(empiricalNObs, overwrite = TRUE)
