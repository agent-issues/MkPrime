# aware-mode-characterisation.R --------------------------------------------
# Characterise the two aware modes in v4 (job 17263327): each PT run sits in
# a different topology mode despite essentially equal log-likelihood. List
# clades unique to each run, then check whether they involve taxa with the
# same vs different ecology coding (a clue to whether the modes are a
# phi/z-symmetry artefact).

suppressPackageStartupMessages({
  library(ape)
  library(TreeTools)
})

outDir <- "inst/ecology/scripts/rodent-comparison"
trees1 <- ape::read.tree(file.path(outDir, "rodent-MkNT-aware_trees_1.nwk"))
trees2 <- ape::read.tree(file.path(outDir, "rodent-MkNT-aware_trees_2.nwk"))
blindTrees <- do.call(c, lapply(Sys.glob(file.path(outDir, "rodent-MkNT-blind_trees_*.nwk")),
                                ape::read.tree))
class(blindTrees) <- "multiPhylo"
class(trees1) <- "multiPhylo"
class(trees2) <- "multiPhylo"

discard <- function(L) L[seq.int(ceiling(length(L)/4) + 1L, length(L))]
t1 <- discard(trees1)
t2 <- discard(trees2)
tB <- discard(blindTrees)

cat("Sample counts (post-burn-in): aware run1 =", length(t1),
    "  aware run2 =", length(t2),
    "  blind =", length(tB), "\n\n")

# --- Per-run majority consensus splits ------------------------------------
tipLabels <- sort(trees1[[1]]$tip.label)
n_all <- length(tipLabels)

tree_splits <- function(trees) {
  # Return: named list of split-strings -> support fraction
  pp <- ape::prop.part(trees)
  labels <- attr(pp, "labels")
  counts <- attr(pp, "number")
  total <- length(trees)
  out <- vapply(seq_along(pp), function(i) {
    s <- sort(labels[pp[[i]]])
    comp <- sort(setdiff(tipLabels, s))
    if (length(s) <= 1 || length(s) >= n_all - 1) return(NA_character_)
    if (length(s) < length(comp)) paste(s, collapse="+")
    else if (length(s) > length(comp)) paste(comp, collapse="+")
    else if (s[1] < comp[1]) paste(s, collapse="+") else paste(comp, collapse="+")
  }, character(1))
  freq <- counts / total
  data.frame(split = out, support = freq, stringsAsFactors = FALSE) |>
    subset(!is.na(split))
}

s1 <- tree_splits(t1)
s2 <- tree_splits(t2)
sB <- tree_splits(tB)

# Majority splits per source
maj_r1 <- s1$split[s1$support >= 0.5]
maj_r2 <- s2$split[s2$support >= 0.5]
maj_blind <- sB$split[sB$support >= 0.5]

cat("Majority (>=0.5) clades:\n")
cat("  Aware run 1: ", length(maj_r1), "\n")
cat("  Aware run 2: ", length(maj_r2), "\n")
cat("  Blind     : ", length(maj_blind), "\n")
cat("  Shared aware run1 & run2: ", length(intersect(maj_r1, maj_r2)), "\n")
cat("  Shared aware run1 & blind: ", length(intersect(maj_r1, maj_blind)), "\n")
cat("  Shared aware run2 & blind: ", length(intersect(maj_r2, maj_blind)), "\n")
cat("  Shared all three: ", length(Reduce(intersect, list(maj_r1, maj_r2, maj_blind))), "\n\n")

# --- Clades present in one run but not the other --------------------------
only_r1 <- setdiff(maj_r1, maj_r2)
only_r2 <- setdiff(maj_r2, maj_r1)

# Look up each clade's support in the opposite run
support_lookup <- function(split, srctab) {
  hit <- srctab[srctab$split == split, "support"]
  if (length(hit) == 0L) 0 else hit[1]
}

cat("==== Clades supported in run 1 (>=0.5) but NOT run 2 ====\n")
for (sp in only_r1) {
  s_other <- support_lookup(sp, s2)
  cat(sprintf("  [run1=%.2f, run2=%.2f] %s\n",
              support_lookup(sp, s1), s_other, sp))
}
cat("\n==== Clades supported in run 2 (>=0.5) but NOT run 1 ====\n")
for (sp in only_r2) {
  s_other <- support_lookup(sp, s1)
  cat(sprintf("  [run1=%.2f, run2=%.2f] %s\n",
              s_other, support_lookup(sp, s2), sp))
}

# --- Ecology-coding check -------------------------------------------------
# Read the ecology vector used in the run
nexFile <- "inst/ecology/data/rodent-X24848.nex"
mat <- TreeTools::ReadCharacters(nexFile)
ecoVec  <- mat[, 220]
extant  <- mat[, 221]
keepTaxa <- which(extant == "0")
ecoKeep <- ecoVec[keepTaxa]
ecoNames <- rownames(mat)[keepTaxa]
# Recode polymorphic + drop missing
poly <- grepl("^\\(", ecoKeep)
ecoKeep[poly] <- substr(sub("^\\(", "", ecoKeep[poly]), 1, 1)
hasEco <- !is.na(ecoKeep) & ecoKeep != "?" & !is.na(suppressWarnings(as.integer(ecoKeep)))
ecoKeep <- ecoKeep[hasEco]
ecoNames <- ecoNames[hasEco]
ecoMap <- setNames(as.integer(ecoKeep), ecoNames)
cat("\n==== Ecology composition of run1-only vs run2-only clades ====\n")
cat("Ecology states tabulated per clade (state of each tip in the split):\n")
report_eco <- function(splits, lbl) {
  for (sp in splits) {
    tips <- strsplit(sp, "+", fixed = TRUE)[[1]]
    tab <- table(ecoMap[tips])
    n <- sum(tab)
    cat(sprintf("  %-6s [n=%2d] eco{%s}: %s\n",
                lbl, n,
                paste(names(tab), collapse=","),
                paste(sprintf("%s=%d", names(tab), tab), collapse=" ")))
  }
}
report_eco(only_r1, "run1")
report_eco(only_r2, "run2")
