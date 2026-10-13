#' Extend the model with parameters from another package
#'
#' An extension contributes parameters, and their prior, to an MkPrime
#' analysis. Each parameter is a positive real that is logged alongside
#' MkPrime's own, checkpointed, carried through tempered-chain swaps and
#' updated by a Bactrian scale move. In this version extension parameters do not
#' enter the likelihood.
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
         init = stats::setNames(info$init, info$params)),
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
  params <- unlist(lapply(extensions, `[[`, "params"))
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
  # Return:
  extensions
}

.ExtensionPointers <- function(model) {
  lapply(model$extensions, function(ext) ext$pointer())
}

.ExtensionParamNames <- function(model) {
  as.character(unlist(lapply(model$extensions, `[[`, "params")))
}

.ExtensionInit <- function(model) {
  unlist(lapply(unname(model$extensions), `[[`, "init"))
}

.ExtensionLogPrior <- function(params, model) {
  lp <- 0
  at <- 0L
  for (ext in model$extensions) {
    n <- length(ext$params)
    lp <- lp + extension_log_prior_cpp(ext$pointer(), params[at + seq_len(n)])
    at <- at + n
  }
  # Return:
  lp
}

# One scale move per extension parameter, addressed past the core scalars.
.WithExtensionMoves <- function(moves, model) {
  extNames <- .ExtensionParamNames(model)
  if (!length(extNames)) return(moves)
  c(moves, lapply(seq_along(extNames), function(j) {
    .ScaleParamMove(extNames[j], index = length(.kScalarParams) + j - 1L)
  }))
}
