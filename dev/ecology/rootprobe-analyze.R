# Robust analysis of the root-invariance probe trees (avoids ape extract.clade /
# tip-label compression by reading each Newick independently and walking edges).
suppressMessages(library(ape))

tf <- "dev/ecology/rootprobe_trees.nwk"
lines <- readLines(tf); lines <- lines[nzchar(trimws(lines))]
trees <- lapply(lines, function(L) ape::read.tree(text = L))
cat(sprintf("read %d trees\n", length(trees)))

tipsUnder <- function(t, node) {
  nT <- ape::Ntip(t)
  if (node <= nT) return(t$tip.label[node])
  out <- character(0); stack <- node
  while (length(stack)) {
    cur <- stack[length(stack)]; stack <- stack[-length(stack)]
    kids <- t$edge[t$edge[, 1] == cur, 2]
    for (k in kids) if (k <= nT) out <- c(out, t$tip.label[k]) else stack <- c(stack, k)
  }
  sort(out)
}
rootBipart <- function(t) {
  nT <- ape::Ntip(t); rt <- nT + 1L
  kids <- t$edge[t$edge[, 1] == rt, 2]
  s <- vapply(kids, function(k) paste(tipsUnder(t, k), collapse = ","), character(1))
  paste(sort(s), collapse = "  |  ")
}
t1IsBasal <- function(t) {
  nT <- ape::Ntip(t); rt <- nT + 1L
  kids <- t$edge[t$edge[, 1] == rt, 2]
  any(vapply(kids, function(k) identical(tipsUnder(t, k), "t1"), logical(1)))
}

rb  <- vapply(trees, rootBipart, character(1))
t1b <- vapply(trees, t1IsBasal, logical(1))
tree0 <- trees[[1]]
rf  <- vapply(trees, function(t) ape::dist.topo(ape::unroot(t), ape::unroot(tree0)), numeric(1))

cat(sprintf("\n=== ROOT PROBE (%d trees) ===\n", length(trees)))
cat(sprintf("distinct root bipartitions across chain : %d\n", length(unique(rb))))
cat(sprintf("t1 still lone basal outgroup            : %d / %d (%.0f%%)\n",
            sum(t1b), length(t1b), 100 * mean(t1b)))
cat(sprintf("unrooted RF from tree[1] : min=%g max=%g  (moves accepted iff max>0)\n",
            min(rf), max(rf)))
cat("\nmost frequent root bipartitions:\n")
print(head(sort(table(rb), decreasing = TRUE), 6))
cat(if (length(unique(rb)) > 1)
      "\nVERDICT: ROOT FLOATS -> EBE-LIVE (Phase 2 must anchor the root)\n"
    else "\nVERDICT: root stable across moves -> moot\n")
