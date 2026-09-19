#' Empirical distribution of observed state counts (`N_obs`)
#'
#' Tabulated from the transformational character matrices of the `neotrans`
#' project corpus (7,017 informative characters across 87 published
#' morphological matrices).  The pmf is supported on `k >= 2`: an explicit
#' body covers `k = 2, ..., nMax`, and a fitted geometric tail extends to
#' `k > nMax`.
#'
#' Used as the empirical component of the `"empirical_geometric"` `k'` prior
#' (see [MkPrimeModel()]), where the prior on the total number of states is
#' the convolution of this empirical pmf with a `Geometric(p)` prior on the
#' number of unobserved states.
#'
#' Users may supply their own object of class `"MkPrimeEmpiricalPrior"` via
#' [MkPrimeEmpiricalPrior()] and pass it to [MkPrimeModel()] through the
#' `empiricalNObs` argument.
#'
#' @format An object of class `"MkPrimeEmpiricalPrior"`: a list with
#' \describe{
#'   \item{body}{Named numeric vector.  `body[i]` is `P(N_obs = i + 1)` for
#'     `i = 1, ..., nMax - 1`; names are the corresponding `k` values
#'     as strings.}
#'   \item{tail_decay}{Geometric decay rate `q` governing the tail beyond the
#'     body, so `P(N_obs = k + 1) = q * P(N_obs = k)` for `k >= tail_start_k`.}
#'   \item{tail_start_k}{Integer: smallest `k` governed by the tail.}
#'   \item{tail_start_p}{Probability mass at `tail_start_k`.}
#'   \item{nSource}{Number of source characters used to build the body.}
#' }
#'
#' @source `data-raw/empirical_n_obs.R`, run against
#'   `https://github.com/ms609/neotrans`.
#'
#' @seealso [MkPrimeEmpiricalPrior()], [MkPrimeModel()]
"empiricalNObs"


#' Construct an empirical prior on observed state counts
#'
#' Build a `"MkPrimeEmpiricalPrior"` object for use as the `empiricalNObs`
#' argument of [MkPrimeModel()].
#'
#' @param body Numeric vector of probabilities `P(N_obs = 2)`,
#'   `P(N_obs = 3)`, ..., `P(N_obs = nMax)`.  Need not sum to 1; the
#'   remaining mass is assumed to lie in the geometric tail.
#' @param tail_decay Geometric decay rate `q` for the tail, with
#'   `0 < q < 1`.  Default `0` (no tail; pmf truncated at `nMax`).  When
#'   `0`, `body` must sum to 1.
#' @param tail_start_k Integer specifying the smallest `k` governed by the
#'   tail; must be `length(body) + 2L`, as the tail is anchored on the last
#'   `body` entry.
#' @param nSource Optional integer.  Number of characters from which the
#'   body was derived (for provenance).
#'
#' @details
#' The full pmf is
#' \deqn{P(N_{obs} = k) = body[k - 1] \quad (k = 2, \ldots, nMax)}
#' \deqn{P(N_{obs} = k) = c \cdot q^{k - tail\_start\_k} \quad (k \ge tail\_start\_k)}
#' with `c = tail_start_p` chosen automatically so the pmf normalises to 1.
#'
#' @return An object of class `"MkPrimeEmpiricalPrior"`.
#' @examples
#' # Uniform on k = 2..5, no tail
#' MkPrimeEmpiricalPrior(body = rep(0.25, 4))
#'
#' # Body plus geometric tail
#' MkPrimeEmpiricalPrior(body = c(0.6, 0.3), tail_decay = 0.4)
#' @export
MkPrimeEmpiricalPrior <- function(body, tail_decay = 0,
                                   tail_start_k = NULL,
                                   nSource = NA_integer_) {
  if (!is.numeric(body) || length(body) < 1L || any(body < 0)) {
    cli::cli_abort("{.arg body} must be a non-negative numeric vector.")
  }
  if (!is.numeric(tail_decay) || length(tail_decay) != 1L ||
      tail_decay < 0 || tail_decay >= 1) {
    cli::cli_abort("{.arg tail_decay} must be a scalar in [0, 1).")
  }
  nBody <- length(body)
  if (is.null(tail_start_k)) {
    tail_start_k <- nBody + 2L
  }
  tail_start_k <- as.integer(tail_start_k)
  # The anchor below is the mass the tail would carry at nBody + 2. Starting
  # it later still normalises to 1, so nothing downstream errors -- it just
  # translates that mass outward, giving a pmf neither argument describes.
  if (is.na(tail_start_k) || tail_start_k != nBody + 2L) {
    cli::cli_abort(c(
      "{.arg tail_start_k} must be {.val {nBody + 2L}}, one past {.arg body}.",
      x = "Got {.val {tail_start_k}}.",
      i = "The tail's mass is fixed by the last {.arg body} entry and
           {.arg tail_decay}. Starting it further out shifts that mass to
           higher {.var k} rather than rescaling it, so the pmf would not be
           the one {.arg body} and {.arg tail_decay} describe."
    ))
  }

  bodySum <- sum(body)
  if (tail_decay == 0) {
    if (abs(bodySum - 1) > 1e-8) {
      cli::cli_abort(
        "With {.arg tail_decay} = 0, {.arg body} must sum to 1 (got
        {round(bodySum, 6)})."
      )
    }
    tailStartP <- 0
  } else {
    # Normalise: body + geometric tail starting with mass tailStartP must
    # sum to 1.  Free parameter: tailStartP.  Anchor the tail so it
    # continues the body's last value (smooth join).
    anchor <- body[nBody] * tail_decay
    tailMass <- anchor / (1 - tail_decay)
    total <- bodySum + tailMass
    body <- body / total
    tailStartP <- anchor / total
  }

  names(body) <- as.character(seq.int(2L, nBody + 1L))

  structure(
    list(
      body = body,
      tail_decay = tail_decay,
      tail_start_k = tail_start_k,
      tail_start_p = tailStartP,
      nSource = nSource
    ),
    class = "MkPrimeEmpiricalPrior"
  )
}


#' @export
print.MkPrimeEmpiricalPrior <- function(x, ...) {
  cli::cli_h1("MkPrime empirical prior on N_obs")
  cli::cli_ul(c(
    "Body: {.val {length(x$body)}} values, kObs in {names(x$body)[1]}-{names(x$body)[length(x$body)]}",
    "Tail decay: {.val {signif(x$tail_decay, 4)}}",
    "Tail starts at k = {.val {x$tail_start_k}} (mass {.val {signif(x$tail_start_p, 4)}})",
    "Source characters: {.val {x$nSource}}"
  ))
  invisible(x)
}
