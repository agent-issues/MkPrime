#' Extend the model with parameters and data from another package
#'
#' An extension contributes parameters, and their prior, to an MkPrime
#' analysis. Each parameter is a positive real that is logged alongside
#' MkPrime's own, checkpointed, carried through tempered-chain swaps and
#' updated by a Bactrian scale move.
#'
#' An extension may also contribute a log-likelihood term for its own data
#' given the tree, and moves of its own on its parameters. Characters an
#' extension models are given to the extension, not to [MkPrimeData()]. Beside
#' an extension likelihood, moves that score candidates with MkPrime's
#' likelihood alone (Gibbs and weighted topology moves, the block Gibbs
#' branch move, and slice sampling of tree length) are not used, and
#' [mkp_stepping_stone()] is unavailable.
#'
#' Extensions may also carry integer parameters (such as latent categories),
#' logged and checkpointed alongside the real ones and updated by the
#' extension's own moves; read MkPrime's scalars (`tree_length`, `rate_loss`,
#' `rate_log_sd`, `rate_neo`, `beta_scale`); and claim MkPrime's neomorphic
#' characters. MkPrime then leaves claimed characters out of its own
#' likelihood and passes their tip states to the extension, which scores
#' them, including any ascertainment correction. They still count towards
#' the normalisation of `rate_neo`, which MkPrime keeps sampling, together
#' with `rate_loss`. Claims are not available with a user partition.
#'
#' Extensions are written in C or C++ against
#' `system.file("include/MkPrime/extension.h", package = "MkPrime")`
#' (`LinkingTo: MkPrime`): a registered `.Call` entry point returns
#' `MkPrimeExtensionPtr(&ext)` for a static `MkPrimeExtension` struct.
#'
#' @param pointer Function of no arguments returning the extension's external
#'   pointer, typically `function() .Call(entryPoint)` defined in the
#'   extension package. It is called again in each process that runs a chain.
#' @return An `MkPrimeExtension` object, to pass to
#'   `MkPrimeModel(extensions)`.
#' @seealso [MkPrimeModel()]
#' @export
MkPrimeExtension <- function(pointer) {
  if (!is.function(pointer)) {
    cli::cli_abort("{.arg pointer} must be a function returning the extension pointer.")
  }
  info <- extension_info_cpp(pointer())
  # Return:
  structure(
    list(pointer = pointer, name = info$name, params = info$params,
         init = stats::setNames(info$init, info$params),
         log_lik = isTRUE(info$log_lik),
         moves = stats::setNames(info$move_weights, info$moves),
         ints = stats::setNames(info$int_init %||% integer(0),
                                info$ints %||% character(0)),
         int_levels = stats::setNames(info$int_levels %||% integer(0),
                                      info$ints %||% character(0)),
         claims_neomorphic = isTRUE(info$claims_neomorphic)),
    class = "MkPrimeExtension"
  )
}

# Column names MkPrime itself logs; an extension parameter may not take one.
.kReservedParamNames <- c(
  "log_posterior", "log_likelihood", "tree_length", "rate_loss",
  "rate_log_sd", "kprime_alpha", "kprime_beta", "p", "rate_neo",
  "beta_scale", "swap_cold", "topo_hash", "hyper_tau"
)

.CheckExtensions <- function(extensions) {
  if (inherits(extensions, "MkPrimeExtension")) extensions <- list(extensions)
  if (!is.list(extensions) ||
      !all(vapply(extensions, inherits, logical(1), "MkPrimeExtension"))) {
    cli::cli_abort(
      "{.arg extensions} must be a list of {.fn MkPrimeExtension} objects."
    )
  }
  params <- unlist(lapply(extensions, function(ext) {
    c(ext$params, names(ext$ints))
  }))
  reserved <- params[params %in% .kReservedParamNames |
                       grepl("^(br|kPrime|w)_|^class[0-9]", params)]
  if (length(reserved)) {
    cli::cli_abort("Extension parameter{?s} {.val {reserved}} clash{?es/} with \\
                    MkPrime's own columns.")
  }
  dup <- unique(params[duplicated(params)])
  if (length(dup)) {
    cli::cli_abort("Extension parameter{?s} {.val {dup}} {?is/are} defined twice.")
  }
  claimants <- vapply(extensions, function(ext) isTRUE(ext$claims_neomorphic),
                      logical(1))
  if (sum(claimants) > 1L) {
    cli::cli_abort("Only one extension may claim the neomorphic characters.")
  }
  # Return:
  extensions
}

.ExtensionPointers <- function(model) {
  lapply(model$extensions, function(ext) ext$pointer())
}

.ExtensionParamNames <- function(model) {
  as.character(unlist(lapply(model$extensions, `[[`, "params")))
}

.ExtensionIntNames <- function(model) {
  as.character(unlist(lapply(model$extensions, function(ext) names(ext$ints))))
}

# Log columns, in the order the engine writes them: real parameters, then
# integer ones.
.ExtensionColumns <- function(model) {
  c(.ExtensionParamNames(model), .ExtensionIntNames(model))
}

.ExtensionInit <- function(model) {
  as.numeric(unlist(lapply(unname(model$extensions), `[[`, "init")))
}

.ExtensionInitInts <- function(model) {
  as.integer(unlist(lapply(unname(model$extensions), `[[`, "ints")))
}

.ExtensionLogPrior <- function(params, model, ints = integer(0)) {
  lp <- 0
  at <- 0L
  atInt <- 0L
  for (ext in model$extensions) {
    n <- length(ext$params)
    nInt <- length(ext$ints)
    lp <- lp + extension_log_prior_cpp(ext$pointer(), params[at + seq_len(n)],
                                       as.integer(ints[atInt + seq_len(nInt)]))
    at <- at + n
    atInt <- atInt + nInt
  }
  # Return:
  lp
}

# Move types that score candidates with the core likelihood alone, so cannot
# run beside an extension likelihood; slice moves likewise.
.kCoreLikOnlyMoveTypes <- c(
  "gibbs_spr", "gibbs_subtree_swap", "weighted_branch_scale",
  "weighted_spr", "weighted_subtree_swap", "block_gibbs_branch", "slice"
)

.ExtensionHasLik <- function(model) {
  any(vapply(model$extensions, function(ext) isTRUE(ext$log_lik), logical(1)))
}

# One scale move per real extension parameter, addressed past the core
# scalars; then each extension's own moves (code 41), addressed by their
# position across extensions.
.WithExtensionMoves <- function(moves, model) {
  if (!length(model$extensions)) return(moves)
  if (.ExtensionHasLik(model)) {
    blocked <- vapply(moves, function(m) {
      isTRUE(m$type %in% .kCoreLikOnlyMoveTypes)
    }, logical(1))
    if (any(blocked)) {
      dropped <- vapply(moves[blocked], `[[`, "", "name")
      .Inform("{length(dropped)} move{?s} not used beside an extension \\
               likelihood: {.val {dropped}}.")
      moves <- moves[!blocked]
    }
  }
  extNames <- .ExtensionParamNames(model)
  scaleMoves <- lapply(seq_along(extNames), function(j) {
    .ScaleParamMove(extNames[j], index = length(.kScalarParams) + j - 1L)
  })
  own <- unlist(lapply(model$extensions, function(ext) {
    lapply(names(ext$moves), function(mv) {
      list(name = paste0(ext$name, "_", mv), type = "extension_move",
           target = ext$name, weight = ext$moves[[mv]], dim = 1L,
           code = 41L)
    })
  }), recursive = FALSE)
  for (j in seq_along(own)) own[[j]]$sliceParamIdx <- j - 1L
  # Return:
  c(moves, scaleMoves, own)
}
