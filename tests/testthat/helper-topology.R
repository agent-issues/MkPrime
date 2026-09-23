# Canonical topology key: sorted non-trivial splits identified by tip names.
# Invariant under root placement and node renumbering, so it compares a rooted
# and an unrooted representation of the same topology.
.TopologyKey <- function(edge, tips) {
  # Nnode from the edge count, never nTip - 1: the state keeps the input tree's
  # rooting, and a wrong Nnode sends TreeTools' rerooting into an endless loop.
  tree <- structure(list(edge = edge, tip.label = tips,
                         Nnode = nrow(edge) - length(tips) + 1L,
                         edge.length = rep(1, nrow(edge))),
                    class = "phylo")
  splits <- as.logical(as.Splits(tree))
  tipNames <- colnames(splits)
  rowKeys <- apply(splits, 1, function(r) {
    sideA <- sort(tipNames[r])
    sideB <- sort(tipNames[!r])
    if (sideA[1] < sideB[1]) paste(sideA, collapse = ",")
    else paste(sideB, collapse = ",")
  })
  # Return:
  paste(sort(rowKeys), collapse = "|")
}
