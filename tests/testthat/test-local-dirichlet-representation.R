# local_dirichlet's edge set must depend on (topology, root, start edge)
# only, never on row or sibling order (#335).

EdgeClades <- function(parent, child, nTip) {
  vapply(seq_along(child), function(i) {
    below <- child[i]
    frontier <- child[i]
    while (length(frontier)) {
      frontier <- child[parent %in% frontier]
      below <- c(below, frontier)
    }
    paste(sort(below[below <= nTip]), collapse = ",")
  }, character(1))
}

# Same tree, shuffled rows and renumbered internal nodes.
ScrambleTree <- function(edge, nTip) {
  edge <- edge[sample.int(nrow(edge)), ]
  internal <- sort(unique(edge[edge > nTip]))
  relabel <- seq_len(max(edge))
  relabel[internal] <- sample(internal)
  matrix(relabel[edge], ncol = 2)
}

SelectionsByStart <- function(edge, nTip, nCats, nDraw) {
  clades <- EdgeClades(edge[, 1], edge[, 2], nTip)
  x <- rep(1 / nrow(edge), nrow(edge))
  sets <- replicate(nDraw, {
    picked <- local_dirichlet_proposal(x, edge[, 1], edge[, 2],
                                       nCats, 2)$modifiedEdges + 1L
    c(clades[picked[1]], paste(sort(clades[picked]), collapse = "|"))
  })
  bySet <- lapply(split(sets[2, ], sets[1, ]), unique)
  bySet[order(names(bySet))]
}

ExpectRepresentationFree <- function(tree, nCats) {
  nTip <- length(tree$tip.label)
  canonical <- TreeTools::Preorder(tree)$edge
  scrambled <- ScrambleTree(canonical, nTip)
  nDraw <- 30 * nrow(canonical)
  a <- SelectionsByStart(canonical, nTip, nCats, nDraw)
  b <- SelectionsByStart(scrambled, nTip, nCats, nDraw)
  expect_length(a, nrow(canonical))
  expect_true(all(lengths(a) == 1L))
  expect_identical(unlist(b), unlist(a))
}

test_that("local_dirichlet selects the same edges in any representation", {
  set.seed(3350)
  ExpectRepresentationFree(ape::rtree(8, rooted = TRUE), 6L)
  ExpectRepresentationFree(ape::rtree(8, rooted = FALSE), 6L)
  ExpectRepresentationFree(ape::rtree(20, rooted = TRUE), 6L)
})
