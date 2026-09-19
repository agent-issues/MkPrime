# User-class partition API for RunMkPrime — Layer 1 plumbing.
#
# This file owns the new `partition` and `unlink` arguments. It is a sibling
# to (not a replacement of) the existing partition machinery in R/partition.R
# (which groups characters by (type, kObs) for the C++ pruning loop). The
# user-class layer wraps that grouping with an additional `classIdx` axis.
#
# Layer 1 currently implements:
#   - argument validation (partition shape, unlink token resolution)
#   - silent coercion of `unlink` when nothing is partitioned across
# and routes EVERY partition != NULL case to a deferred-implementation
# cli_abort. That keeps the §7a bit-identity guarantee for partition = NULL
# in force while later commits add the user-class state and C++ surface.

# Valid `unlink` tokens for Layer 1. `brlens` is reserved for Layer 2
# (T3 treatment) and rejected here with a "not yet implemented" message —
# it must validate cleanly through the prefix matcher (so users can pass
# uniform `unlink` vectors across treatments) but error before any code
# tries to honour it.
.kPartitionApiUnlinkTokens <- c("shape", "ratemultiplier", "brlens")
.kPartitionApiUnlinkLayer1  <- c("shape", "ratemultiplier")


# Validate and normalise `partition` and `unlink`.
#
# Returns a list with components:
#   partition  - NULL or an integer vector of length mkd$nChar with values
#                that form a contiguous range 1:nClasses.
#   unlink     - character vector of canonical tokens (subset of
#                .kPartitionApiUnlinkLayer1) with no duplicates. Always
#                character(0) when partition is NULL or nClasses == 1.
#   nClasses   - integer scalar (1L when partition is NULL).
#
# Errors hard on: wrong type, length mismatch, NAs, non-contiguous class
# ids, ambiguous unlink prefix, unknown unlink token (with `agrep` "did
# you mean" suggestion). Warns on resolvable partial-prefix matches.
.ValidatePartitionArgs <- function(partition, unlink, mkd) {

  # --- Unlink token resolution -------------------------------------------
  if (!is.null(unlink) && !is.character(unlink)) {
    cli::cli_abort("{.arg unlink} must be a character vector.")
  }
  unlink <- as.character(unlink)
  unlinkIn <- unlink

  resolved <- character(length(unlink))
  for (i in seq_along(unlink)) {
    raw <- unlink[i]
    if (!nzchar(raw)) {
      cli::cli_abort("{.arg unlink} contains an empty string at position {i}.")
    }
    tok <- tolower(raw)
    if (tok %in% .kPartitionApiUnlinkTokens) {
      resolved[i] <- tok
    } else {
      hits <- .kPartitionApiUnlinkTokens[startsWith(.kPartitionApiUnlinkTokens, tok)]
      if (length(hits) == 1L) {
        cli::cli_warn(
          "{.arg unlink} token {.val {unlinkIn[i]}} matched via prefix to {.val {hits}}."
        )
        resolved[i] <- hits
      } else if (length(hits) > 1L) {
        cli::cli_abort(c(
          "{.arg unlink} token {.val {unlinkIn[i]}} is an ambiguous prefix.",
          "i" = "Matches: {.val {hits}}."
        ))
      } else {
        validTokens <- .kPartitionApiUnlinkTokens
        sug <- agrep(tok, validTokens, value = TRUE,
                     max.distance = 0.4, ignore.case = TRUE)
        msgs <- c(
          "{.arg unlink} token {.val {unlinkIn[i]}} is not recognised.",
          "i" = "Valid tokens: {.val {validTokens}}."
        )
        if (length(sug) > 0L) {
          msgs <- c(msgs, "i" = "Did you mean: {.val {sug}}?")
        }
        cli::cli_abort(msgs)
      }
    }
  }
  if (anyDuplicated(resolved)) {
    cli::cli_warn("{.arg unlink} has duplicate tokens; deduplicating.")
    resolved <- unique(resolved)
  }
  unlink <- resolved

  # --- Partition shape validation ---------------------------------------
  nClasses <- 1L
  if (!is.null(partition)) {
    if (!is.numeric(partition)) {
      cli::cli_abort("{.arg partition} must be an integer vector or {.val NULL}.")
    }
    if (anyNA(partition)) {
      cli::cli_abort("{.arg partition} contains {.val NA}; supply an integer vector with no missing values.")
    }
    if (any(partition != floor(partition))) {
      cli::cli_abort("{.arg partition} must contain whole-number class IDs only.")
    }
    partition <- as.integer(partition)
    if (length(partition) != mkd$nChar) {
      cli::cli_abort(c(
        "{.arg partition} has length {length(partition)} but the data has \\
        {mkd$nChar} character{?s} (after invariant-character drop).",
        "i" = "Supply a length-{mkd$nChar} integer vector whose values form a \\
              contiguous range {.val 1}:{.val nClasses}."
      ))
    }
    if (any(partition < 1L)) {
      cli::cli_abort("{.arg partition} must contain values {.val >= 1}; got values < 1.")
    }
    nClasses <- max(partition)
    observed <- sort(unique(partition))
    expected <- seq_len(nClasses)
    if (!identical(observed, expected)) {
      missingCls <- setdiff(expected, observed)
      cli::cli_abort(c(
        "{.arg partition} skips class id{?s} {.val {missingCls}}.",
        "i" = "Class labels must form a contiguous range {.val 1}:{.val nClasses}."
      ))
    }
  }

  # --- Silent coercion when partition is trivial ------------------------
  if ((is.null(partition) || nClasses == 1L) && length(unlink) > 0L) {
    why <- if (is.null(partition)) {
      "no {.arg partition} supplied"
    } else {
      "only one user class in {.arg partition}"
    }
    cli::cli_alert_info(
      "{.arg unlink} ignored ({why}; nothing to unlink across)."
    )
    unlink <- character(0)
  }

  # Return:
  list(partition = partition, unlink = unlink, nClasses = nClasses)
}


# Layer-1 implementation gate. Called from RunMkPrime once validation has
# produced a (possibly non-trivial) partition spec. This commit opens the
# gate for the trivial spec (nClasses == 1, unlink = character(0)) which
# routes through cpp_log_likelihood_partitioned with length-1 per-class
# vectors — the §7b numeric-equivalence regime. Real partitioning (nClasses
# > 1 or non-empty unlink) stays closed until the per-class moves land.
.RequirePartitionImplemented <- function(spec) {
  # Trivial spec: partition = NULL (legacy path) OR the nClasses == 1 /
  # unlink = character(0) case that the MCMC loop now handles via the
  # partitioned likelihood (§7b contract).
  if (is.null(spec$partition) ||
      (spec$nClasses == 1L && length(spec$unlink) == 0L)) {
    # Return:
    return(invisible(NULL))
  }

  # Layer 1 gate: multi-class with "shape" and/or "ratemultiplier" unlink is
  # now supported (per-class moves implemented). "brlens" stays deferred
  # to Layer 2.
  if ("brlens" %in% spec$unlink) {
    cli::cli_abort(c(
      "{.val brlens} unlink is not yet implemented (Layer 2).",
      "i" = "Supported tokens in Layer 1: {.val shape}, {.val ratemultiplier}."
    ))
  }
  if (spec$nClasses > 1L &&
      !all(spec$unlink %in% .kPartitionApiUnlinkLayer1)) {
    cli::cli_abort(c(
      "Unsupported {.arg unlink} token(s) for multi-class partition.",
      "i" = "Layer 1 supports: {.val {.kPartitionApiUnlinkLayer1}}."
    ))
  }
}


# Derive per-class character counts from a partition vector.
# Returns integer vector of length nClasses with the number of characters
# in each user class. Used to convert between the Dirichlet w-simplex
# (char-weighted simplex) and class-rate vectors (§5.1).
.PartitionNCharPerClass <- function(partition, nClasses) {
  out <- tabulate(partition, nbins = nClasses)
  if (any(out == 0L)) {
    cli::cli_abort(
      "Empty user class{?es} {.val {which(out == 0L)}} in partition vector."
    )
  }
  # Return:
  as.integer(out)
}


# Convert a w-simplex (char-weighted) to per-class rate.
#   class_rate[c] = w[c] * nChar / nChar_c
# Inverse: w[c] = class_rate[c] * nChar_c / nChar.
# The mean-1 constraint sum(nChar_c * class_rate[c]) / nChar = 1 (§5.1) is
# satisfied iff w is on the unit simplex (sum(w) = 1).
.PartitionWToClassRate <- function(w, nCharPerClass) {
  nChar <- sum(nCharPerClass)
  # Return:
  as.numeric(w) * nChar / as.numeric(nCharPerClass)
}

.PartitionClassRateToW <- function(classRate, nCharPerClass) {
  nChar <- sum(nCharPerClass)
  # Return:
  as.numeric(classRate) * as.numeric(nCharPerClass) / nChar
}


# Initial state for the partition-aware MCMC path.
#
# Wraps .InitState (which is unchanged — §7a contract) and adds per-class
# fields when partitionSpec is non-trivial. Returned state contains:
#
#   - Everything in legacy state (tree, tree_length, rel_br_lengths,
#     rate_loss, rate_log_sd, kPrime, optional p / kprime_alpha+beta /
#     rate_neo / beta_scale, log_lik / log_prior / log_post).
#   - state$nChar_c: integer vector length nClasses.
#   - state$class_rate_log_sd: numeric vector. Length nClasses when "shape"
#     is unlinked (each c starts at rate_log_sd = 0.5); length 1 when
#     linked (== rate_log_sd).
#   - state$class_w: numeric vector length nClasses on the unit simplex
#     (Dirichlet domain). Initialised to nChar_c / nChar — the
#     "weights proportional to character count" prior mean, which gives
#     class_rate ≡ 1 in every class (matches the legacy rate = 1).
#   - state$class_rate: numeric vector length nClasses derived from w
#     (initially identically 1.0 by construction).
#   - state$eta_neo: scalar, 1.0 by default. Only consumed when hasNeo.
#
# log_lik is recomputed via cpp_log_likelihood_partitioned (sibling) so
# the initial likelihood reflects the partition-aware path even when
# class_rate is identically 1 (the §7b numeric-equivalence regime).
.InitStatePartitioned <- function(tree, mkd, model, partitionSpec) {
  # Build base state via the existing legacy initializer (bit-identical
  # to the §7a-locked path).
  state <- .InitState(tree, mkd, model)

  # Trivial spec: nothing to add. Caller can treat as legacy.
  if (is.null(partitionSpec$partition) || partitionSpec$nClasses == 1L) {
    state$nChar_c          <- as.integer(mkd$nChar)
    state$class_w          <- 1.0
    state$class_rate       <- 1.0
    state$class_rate_log_sd <- state$rate_log_sd
    state$eta_neo          <- 1.0
    # Return:
    return(state)
  }

  nClasses <- partitionSpec$nClasses
  nChar    <- mkd$nChar
  nCharPC  <- .PartitionNCharPerClass(partitionSpec$partition, nClasses)
  state$nChar_c <- nCharPC

  # class_w initialised so class_rate == 1.0 for every class:
  #   class_rate[c] = w[c] * nChar / nChar_c == 1
  #   => w[c] = nChar_c / nChar  (a valid simplex point: sum = 1)
  state$class_w    <- nCharPC / nChar
  state$class_rate <- .PartitionWToClassRate(state$class_w, nCharPC)

  # class_rate_log_sd is per-class only when "shape" is unlinked. Each c
  # starts at the legacy value 0.5 so the partition-aware path collapses
  # to the legacy LL at initial state regardless of unlink-shape choice.
  shapeUnlinked <- "shape" %in% partitionSpec$unlink
  if (shapeUnlinked) {
    state$class_rate_log_sd <- rep(state$rate_log_sd, nClasses)
  } else {
    state$class_rate_log_sd <- state$rate_log_sd
  }

  # Half-normal hyperprior on σ_c (non-centred): σ_c = hyper_tau · z_c.
  # Active only when shape is unlinked AND the model selects the pooled
  # hyperprior AND nClasses >= 2 (the structure is degenerate at K = 1;
  # the legacy single-σ Gamma prior covers that case). With τ = 1 and
  # z_c = rate_log_sd (= 0.5 by default), the derived σ_c equals the
  # legacy initial value, preserving the initial-state likelihood
  # equivalence to the legacy single-σ path.
  state$use_hyperprior_on_sigma <- shapeUnlinked &&
    nClasses >= 2L &&
    identical(model$priorOnClassRateLogSd, "hyperprior_pooled")
  if (state$use_hyperprior_on_sigma) {
    state$hyper_tau           <- 1.0
    state$class_rate_log_sd_z <- rep(state$rate_log_sd, nClasses)
  }

  # eta_neo is the per-§5.2 asymmetry parameter; only sampled when hasNeo.
  # Layer 1 freezes it at 1.0 (the no-asymmetry geometric mean) because
  # the Q-matrix-asymmetry semantics need clarification before the
  # eta_neo != 1 path can be honoured by cpp_log_likelihood_partitioned.
  # For Casali (hasNeo == FALSE) this is irrelevant; eta_neo is never
  # consulted.
  state$eta_neo <- 1.0

  # Recompute log_lik via the partition-aware sibling so the initial value
  # is consistent with the path the chain will subsequently take. At the
  # trivial class_rate ≡ 1.0 starting point this MUST equal the legacy
  # log_lik to ~1e-10 (§7b contract).
  tree_st <- state$tree
  parent  <- tree_st$edge[, 1]
  child   <- tree_st$edge[, 2]
  edgeLen <- tree_st$edge.length

  # The chain's own builder, so the initial value tracks the model's k'-prior
  # flags. A second McmcData here would not.
  dataPtr <- .InitMcmcData(mkd, model)

  state$log_lik <- cpp_log_likelihood_partitioned_xptr(
    dataPtr, parent, child, edgeLen, as.integer(state$kPrime),
    rateLoss   = state$rate_loss,
    rateLogSd  = state$class_rate_log_sd,
    classRate  = state$class_rate,
    etaNeo     = state$eta_neo,
    betaScale  = state$beta_scale %||% 1.0
  )
  # Recompute log_prior with the extended LogPrior() that handles per-class
  # fields (class_w, class_rate_log_sd). At the initial point (class_w ==
  # nChar_c / nChar, class_rate_log_sd == rate_log_sd for all classes),
  # the partitioned prior equals the legacy prior to ~1e-10 (§7b analogue).
  state$log_prior <- LogPrior(state, model, mkd)
  state$log_post  <- state$log_lik + state$log_prior

  # Return:
  state
}


# Build the move spec for the partition-aware MCMC path.
#
# Wraps .BuildMoves (which is unchanged — §7a contract for the legacy path).
# When `partitionSpec` is trivial, returns the legacy spec unchanged. When
# non-trivial, adds per-class moves consistent with the unlinked components:
#
#   "shape" in spec$unlink           -> scale_class_rate_log_sd_c per class
#   "ratemultiplier" in spec$unlink  -> dirichlet_simplex on class_w
#
# Per-class move specs use new `type` strings that the C++ dispatcher must
# learn to handle in a follow-up commit; until then they are listed in the
# move schedule but the .RequirePartitionImplemented gate prevents execution
# from ever reaching the dispatcher.
.BuildMovesPartitioned <- function(nEdge, nTrans, hasNeo, mcmc,
                                   partitionSpec,
                                   fixTopology   = FALSE,
                                   kPrimePrior   = "geometric",
                                   qHeterogeneity = FALSE,
                                   joint2d       = TRUE,
                                   priorOnClassRateLogSd = "hyperprior_pooled",
                                   likelihoodMode = "sampled_k") {
  # Always build the legacy spec first; trivial partition is a no-op.
  moves <- .BuildMoves(nEdge, nTrans, hasNeo, mcmc,
                       fixTopology   = fixTopology,
                       kPrimePrior   = kPrimePrior,
                       qHeterogeneity = qHeterogeneity,
                       joint2d       = joint2d,
                       likelihoodMode = likelihoodMode)

  if (is.null(partitionSpec$partition) || partitionSpec$nClasses == 1L) {
    # Return:
    return(moves)
  }

  nClasses <- partitionSpec$nClasses

  # "shape" unlinked: drop the legacy scalar rate_log_sd moves (case 2,
  # case 19-slice on rate_log_sd, and the joint_tl_rls 2D move). The
  # scalar `state->rateLogSd` only stays meaningful as the lockstep
  # mirror of `classRateLogSd[0]`; an independent scalar move would
  # break that lockstep and the partial-CL fallbacks in the C++ hot
  # path that consume state->rateLogSd would return wrong (or NaN)
  # log-likelihoods. Per-class moves (cases 31, 33) target σ_0
  # explicitly via classIdx = 1.
  if ("shape" %in% partitionSpec$unlink) {
    keep <- vapply(moves, function(mv) {
      !(identical(mv$target, "rate_log_sd") ||
        identical(mv$name, "rate_log_sd") ||
        identical(mv$name, "slice_rate_log_sd") ||
        identical(mv$name, "joint_tl_rls"))
    }, logical(1))
    moves <- moves[keep]
  }

  # "shape" unlinked: one MH scale move per class on class_rate_log_sd[c].
  # Each move targets a single component of a length-nClasses state vector.
  # Under the pooled hyperprior the C++ dispatch internally rescales z_c
  # (not σ_c) and derives σ_c = τ·z_c; the same move spec drives both
  # arms.
  if ("shape" %in% partitionSpec$unlink) {
    for (c in seq_len(nClasses)) {
      moves <- c(moves, list(
        list(
          name      = paste0("scale_class_rate_log_sd_", c),
          type      = "scale_class_rate_log_sd",
          target    = "class_rate_log_sd",
          weight    = 1,
          dim       = 1L,
          classIdx  = as.integer(c)
        )
      ))
    }

    # Single global Bactrian scale on τ (only meaningful when the pooled
    # hyperprior is active; rescaling τ moves every σ_c together).
    if (identical(priorOnClassRateLogSd, "hyperprior_pooled")) {
      moves <- c(moves, list(
        list(
          name   = "scale_hyper_tau",
          type   = "scale_hyper_tau",
          target = "hyper_tau",
          weight = 1,
          dim    = 1L
        )
      ))
    }
  }

  # "ratemultiplier" unlinked: one Dirichlet-simplex move on class_w
  # (length nClasses). class_rate is derived from w via the char-weighted
  # mean-1 map.
  if ("ratemultiplier" %in% partitionSpec$unlink) {
    moves <- c(moves, list(
      list(
        name   = "dirichlet_simplex_class_w",
        type   = "dirichlet_simplex_class_w",
        target = "class_w",
        weight = max(1, nClasses),
        dim    = as.integer(nClasses)
      )
    ))
  }

  # Return:
  moves
}


# Parameter name vector for the partition-aware sample matrix.
#
# Extends .ParamNames (unchanged — §7a contract on the legacy schema) with
# per-class columns when partitionSpec is non-trivial. Column-naming
# convention is `class<c>_rate_log_sd`, `w_<c>` (mirroring plan v4 §8 row 6).
.ParamNamesPartitioned <- function(mkd, nEdge, partitionSpec,
                                    kPrimePrior   = "geometric",
                                    qHeterogeneity = FALSE,
                                    priorOnClassRateLogSd = "hyperprior_pooled") {
  nms <- .ParamNames(mkd, nEdge, kPrimePrior = kPrimePrior,
                     qHeterogeneity = qHeterogeneity)

  if (is.null(partitionSpec$partition) || partitionSpec$nClasses == 1L) {
    # Return:
    return(nms)
  }

  nClasses <- partitionSpec$nClasses
  extra <- character(0)

  if ("shape" %in% partitionSpec$unlink) {
    extra <- c(extra, paste0("class", seq_len(nClasses), "_rate_log_sd"))
  }
  if ("ratemultiplier" %in% partitionSpec$unlink) {
    extra <- c(extra, paste0("w_", seq_len(nClasses)))
  }

  # Hyperprior columns: hyper_tau + per-class z_c. Appended last to keep
  # the legacy column order stable for AutoPart-style readers that index
  # by name and rely on prior column positions.
  if ("shape" %in% partitionSpec$unlink &&
      nClasses >= 2L &&
      identical(priorOnClassRateLogSd, "hyperprior_pooled")) {
    extra <- c(extra,
               "hyper_tau",
               paste0("class", seq_len(nClasses), "_rate_log_sd_z"))
  }

  # Order: legacy columns first, then per-class additions appended.
  # Downstream readers (AutoPart's 03_diagnostics.R / 04_metrics.R) read
  # columns by name, so appending is safe — see plan v4 §8 row 6.

  # Return:
  c(nms, extra)
}
