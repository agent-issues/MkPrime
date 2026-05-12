# sim3-simulate.R -------------------------------------------------------------
#
# Forward simulator for the convergent-ecology study.  Implements a
# transformational-only k-state Mk model with per-edge ecology rate
# modifiers driven by the spike-and-slab z latent:
#
#   - z_{c, e} = 0  →  edge rate factor = 1     (baseline)
#   - z_{c, e} = 1  →  edge rate factor = phi   (encouraged)
#   - z_{c, e} = 2  →  edge rate factor = 1/phi (discouraged)
#
# Neomorphic characters are deferred to a later sweep; Sim 3 uses
# transformational characters only.
#
# The simulator takes:
#   * tree            — phylo (Preorder), edge.length set
#   * edgeEcology     — integer vector length nEdge, values in 0..(kEco-1)
#   * z               — integer matrix nChar x kEco, values in {0L, 1L, 2L}
#   * phi             — positive scalar
#   * baseRate        — positive scalar; pre-modifier per-character rate
#   * kStates         — integer vector length nChar (number of states per char)
#
# Returns a character matrix nTip x nChar (states as "0", "1", ...).


#' Assign edge ecology by Fitch-like parsimony from tip ecology.
#'
#' For the canonical Sim 3 tree this returns the expected "obvious"
#' assignment: edges descending from a clade with homogeneous tip
#' ecology inherit that ecology; the root edge (and other branches
#' where the parsimony score is tied) defaults to ecology 0.
#'
#' Tie-breaking uses the smallest state index in the Fitch set; that
#' is sufficient for the symmetric Sim 3 design but should be
#' revisited if reused elsewhere.
#'
#' @param tree A `phylo` object in Preorder.
#' @param tipEcology Integer vector of length `NTip(tree)`, named by
#'   tip label, with non-negative values (`NA` allowed at tips and
#'   treated as fully ambiguous).
#' @return Integer vector of length `nrow(tree$edge)` giving the
#'   parent-side ecology state on each edge.
#'
#' @importFrom TreeTools NTip
.AssignEdgeEcology <- function(tree, tipEcology) {
  edge <- tree$edge
  nTip <- TreeTools::NTip(tree)
  nNode <- tree$Nnode
  nAll <- nTip + nNode
  tipEcology <- tipEcology[tree$tip.label]
  states <- sort(unique(stats::na.omit(tipEcology)))
  if (length(states) == 0L) {
    cli::cli_abort("All tip ecology values are NA.")
  }
  # Fitch pass 1 (postorder): compute candidate-state set at each node.
  cset <- vector("list", nAll)
  for (i in seq_len(nTip)) {
    cset[[i]] <- if (is.na(tipEcology[i])) states else tipEcology[i]
  }
  postIdx <- rev(seq_len(nrow(edge)))
  for (i in postIdx) {
    parent <- edge[i, 1L]
    child  <- edge[i, 2L]
    if (is.null(cset[[parent]])) {
      cset[[parent]] <- cset[[child]]
    } else {
      inter <- intersect(cset[[parent]], cset[[child]])
      cset[[parent]] <- if (length(inter) > 0L) inter
                       else union(cset[[parent]], cset[[child]])
    }
  }
  # Fitch pass 2 (preorder): pick a single state per node.
  state <- integer(nAll)
  rootIdx <- edge[1L, 1L]
  state[rootIdx] <- min(cset[[rootIdx]])
  for (i in seq_len(nrow(edge))) {
    parent <- edge[i, 1L]
    child  <- edge[i, 2L]
    if (child <= nTip) {
      state[child] <- if (is.na(tipEcology[child])) state[parent]
                      else tipEcology[child]
    } else if (state[parent] %in% cset[[child]]) {
      state[child] <- state[parent]
    } else {
      state[child] <- min(cset[[child]])
    }
  }
  # Parent-side ecology on each edge.
  as.integer(state[edge[, 1L]])
}


#' Forward-simulate characters with ecology-modulated rates.
#'
#' Per-character substitution dispatches on `type`:
#'
#' * **transformational** — symmetric k-state Mk (JC-style):
#'     `r(c, e) = baseRate * mu(z_{c, e}, phi)`
#'     where `mu(0, phi) = 1`, `mu(1, phi) = phi`, `mu(2, phi) = 1/phi`.
#'     Off-diagonal `(1 - exp(-k·r·t/(k-1)))/k`.
#'
#' * **neomorphic** — two-state asymmetric CTMC (gain/loss model from
#'   the M2-NT family of [neotrans]).  Base rates from a single ratio
#'   `rateLoss` (loss/gain):
#'     `rate01 = 2 / (1 + rateLoss)`        — gain
#'     `rate10 = 2 * rateLoss / (1 + rateLoss)` — loss
#'   Under `z_{c, e}`:
#'     * `z = 1` (encouraged): `(rate01 * phi, rate10 / phi)`
#'     * `z = 2` (discouraged): `(rate01 / phi, rate10 * phi)`
#'   Root drawn from edge-stationary distribution under the rates at
#'   the root edge.  This is the asymmetric mechanism that produces
#'   genuine state-level convergence between unrelated ecology-1 clades.
#'
#' @param tree Phylo object (Preorder).
#' @param edgeEcology Output of [.AssignEdgeEcology()].
#' @param z Integer matrix nChar x kEco with values in `{0L, 1L, 2L}`.
#' @param phi Positive scalar.
#' @param baseRate Positive scalar; per-character per-edge baseline
#'   rate (transformational characters only).
#' @param kStates Integer vector length nChar; states per character
#'   (transformational only).  Ignored for neomorphic (always 2).
#' @param type Character vector length nChar; `"transformational"` or
#'   `"neomorphic"`.  Recycled to nChar if length 1.
#' @param rateLoss Positive scalar; loss/gain ratio for neomorphic
#'   characters.  `rateLoss = 1` gives symmetric gain/loss.
#' @return A character matrix nTip x nChar; rownames = tip labels;
#'   states encoded as "0".."(k-1)".
#'
#' @importFrom TreeTools NTip
.SimulateMkPrimeEcology <- function(tree, edgeEcology, z, phi,
                                    baseRate = 1, kStates = 2L,
                                    type = "transformational",
                                    rateLoss = 1) {
  stopifnot(phi > 0, baseRate > 0, rateLoss > 0)
  edge   <- tree$edge
  brLen  <- tree$edge.length
  nEdge  <- nrow(edge)
  nTip   <- TreeTools::NTip(tree)
  nNode  <- tree$Nnode
  nAll   <- nTip + nNode
  nChar  <- nrow(z)
  if (length(kStates) == 1L) kStates <- rep(as.integer(kStates), nChar)
  if (length(type) == 1L) type <- rep(type, nChar)
  stopifnot(length(kStates) == nChar)
  stopifnot(length(type) == nChar)
  stopifnot(all(type %in% c("transformational", "neomorphic")))
  stopifnot(length(edgeEcology) == nEdge)
  mult <- c(1, phi, 1 / phi)  # indexed by z + 1L
  # Asymmetric Q for neomorphic, parametrised as in the M2-NT model.
  rate01Base <- 2 / (1 + rateLoss)
  rate10Base <- 2 * rateLoss / (1 + rateLoss)

  # Root state per character. Transformational: uniform.  Neomorphic:
  # stationary under the root-edge rate pair (ignores z; for sim-3 the
  # root edge has eco 0, so z = 0 and base rates apply).
  rootIdx <- edge[1L, 1L]
  stateMat <- matrix(NA_integer_, nrow = nAll, ncol = nChar)
  pi1Root <- rate01Base / (rate01Base + rate10Base)
  for (c in seq_len(nChar)) {
    if (type[c] == "neomorphic") {
      stateMat[rootIdx, c] <- as.integer(stats::runif(1) < pi1Root)
    } else {
      stateMat[rootIdx, c] <- sample.int(kStates[c], 1L) - 1L
    }
  }
  # Preorder walk.
  for (i in seq_len(nEdge)) {
    parent <- edge[i, 1L]
    child  <- edge[i, 2L]
    t      <- brLen[i]
    e      <- edgeEcology[i] + 1L  # 1-based for indexing
    for (c in seq_len(nChar)) {
      zce <- z[c, e]
      parentSt <- stateMat[parent, c]
      if (type[c] == "neomorphic") {
        # Asymmetric two-state CTMC with per-edge (rate01, rate10).
        if (zce == 0L) {
          a <- rate01Base; b <- rate10Base
        } else if (zce == 1L) {
          a <- rate01Base * phi; b <- rate10Base / phi
        } else {
          a <- rate01Base / phi; b <- rate10Base * phi
        }
        denom <- a + b
        pi1 <- a / denom
        pi0 <- b / denom
        ee  <- exp(-denom * t)
        if (parentSt == 0L) {
          pTo1 <- pi1 * (1 - ee)
          stateMat[child, c] <- as.integer(stats::runif(1) < pTo1)
        } else {
          pTo0 <- pi0 * (1 - ee)
          stateMat[child, c] <- 1L - as.integer(stats::runif(1) < pTo0)
        }
      } else {
        k <- kStates[c]
        if (k < 2L) {
          stateMat[child, c] <- parentSt
          next
        }
        r <- baseRate * mult[zce + 1L]
        pOff  <- (1 - exp(-k * r * t / (k - 1))) / k
        pSame <- 1 - (k - 1) * pOff
        if (stats::runif(1) < pSame) {
          stateMat[child, c] <- parentSt
        } else {
          choices <- setdiff(0:(k - 1L), parentSt)
          stateMat[child, c] <- as.integer(
            choices[sample.int(length(choices), 1L)]
          )
        }
      }
    }
  }
  tipStates <- stateMat[seq_len(nTip), , drop = FALSE]
  rownames(tipStates) <- tree$tip.label
  storage.mode(tipStates) <- "character"
  tipStates
}
