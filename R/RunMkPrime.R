# Main MCMC entry point for MkPrime
#
# Phase 3: single chain, fixed topology, R-side loop.

#' Run Bayesian MCMC under the MkPrime model
#'
#' Metropolis-Hastings MCMC on a fixed tree topology, sampling branch
#' lengths, model parameters, and per-character k' (for transformational
#' characters).
#'
#' @param data A `phyDat` object or `MkPrimeData` object.
#' @param tree A `phylo` object (the fixed topology).
#' @param neomorphic,known_states Passed to [MkPrimeData()] if `data` is
#'   a `phyDat` object.
#' @param model An `MkPrimeModel` object, or `NULL` for defaults.
#' @param mcmc An `MkPrimeMCMC` object, or `NULL` for defaults.
#'
#' @return An `MkPosterior` object.
#' @export
RunMkPrime <- function(data, tree,
                       neomorphic = integer(0),
                       known_states = integer(0),
                       model = NULL,
                       mcmc = NULL) {

  # --- Input processing ---
  if (inherits(data, "MkPrimeData")) {
    mkd <- data
  } else {
    mkd <- MkPrimeData(data, neomorphic = neomorphic,
                       known_states = known_states)
  }

  if (!inherits(tree, "phylo")) {
    cli::cli_abort("{.arg tree} must be a {.cls phylo} object.")
  }

  if (is.null(model)) model <- MkPrimeModel()
  if (is.null(mcmc)) mcmc <- MkPrimeMCMC()

  model <- .finalize_model(model, tree, mkd)

  # Precompute tree structure (fixed topology)
  tree <- ape::reorder.phylo(tree, "postorder")
  nEdge <- nrow(tree$edge)

  has_neo <- any(mkd$type == "neomorphic")
  trans_idx <- which(mkd$type == "transformational")
  has_trans <- length(trans_idx) > 0
  nTrans <- length(trans_idx)

  # --- Initialize MCMC state ---
  state <- .init_state(tree, mkd, model)

  # --- Build move schedule ---
  moves <- .build_moves(nEdge, nTrans, has_neo, mcmc)

  # --- Pre-allocate sample storage ---
  nSaved <- as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  param_names <- .param_names(mkd, nEdge)
  samples <- matrix(NA_real_, nrow = nSaved, ncol = length(param_names),
                    dimnames = list(NULL, param_names))
  tree_samples <- vector("list", nSaved)
  saved_idx <- 0L

  # Acceptance tracking
  accept_count <- integer(length(moves))
  propose_count <- integer(length(moves))
  names(accept_count) <- names(propose_count) <- vapply(
    moves, `[[`, character(1), "name"
  )

  # Adaptation state
  tuning <- mcmc$tuning

  # --- MCMC loop ---
  cli::cli_progress_bar(
    "MCMC", total = mcmc$nIter,
    format = "{cli::pb_bar} {cli::pb_current}/{cli::pb_total} | logP: {format(round(state$log_post, 1), nsmall = 1)} | acc: {format(round(recent_acc, 2), nsmall = 2)}"
  )
  recent_acc <- 0
  recent_window <- 500L
  recent_accepts <- logical(recent_window)
  recent_pos <- 0L

  for (iter in seq_len(mcmc$nIter)) {
    # Select move
    move_idx <- sample.int(length(moves), 1L, prob = vapply(
      moves, `[[`, numeric(1), "weight"
    ))
    move <- moves[[move_idx]]
    propose_count[move_idx] <- propose_count[move_idx] + 1L

    # Generate proposal
    accepted <- .do_move(move, state, mkd, model, tree, tuning)

    if (accepted$accept) {
      state <- accepted$state
      accept_count[move_idx] <- accept_count[move_idx] + 1L
    }

    # Track recent acceptance
    recent_pos <- (recent_pos %% recent_window) + 1L
    recent_accepts[recent_pos] <- accepted$accept

    # Adaptation during warmup
    if (iter <= mcmc$warmup && iter %% 200L == 0L) {
      tuning <- .adapt_tuning(tuning, accept_count, propose_count, moves)
    }

    # Save samples post-warmup
    if (iter > mcmc$warmup && (iter - mcmc$warmup) %% mcmc$thin == 0L) {
      saved_idx <- saved_idx + 1L
      samples[saved_idx, ] <- .state_to_row(state, mkd, nEdge)
      tree_samples[[saved_idx]] <- .state_to_tree(state, tree)
    }

    if (iter %% 100L == 0L) {
      n_recent <- min(iter, recent_window)
      recent_acc <- sum(recent_accepts[seq_len(n_recent)]) / n_recent
      cli::cli_progress_update()
    }
  }
  cli::cli_progress_done()

  # --- Build result ---
  MkPosterior(
    samples = samples,
    trees = tree_samples,
    acceptance = accept_count / pmax(propose_count, 1L),
    model = model,
    data = mkd,
    mcmc = mcmc,
    warmup = mcmc$warmup,
    tuning = tuning
  )
}


# --- Internal helpers ---

#' Initialize MCMC state from tree and data
#' @keywords internal
.init_state <- function(tree, mkd, model) {
  tree_length <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tree_length

  # Default k' = kObs for all characters
  kPrime <- mkd$kObs
  known_idx <- which(mkd$type == "known")
  if (length(known_idx)) kPrime[known_idx] <- mkd$known_k[known_idx]

  state <- list(
    tree_length = tree_length,
    rel_br_lengths = rel_br,
    rate_loss = 1.0,
    rate_log_sd = 0.5,
    kPrime = as.integer(kPrime),
    p = 0.5
  )

  # Compute initial likelihood and prior
  edge_length <- state$tree_length * state$rel_br_lengths
  tmp_tree <- tree
  tmp_tree$edge.length <- edge_length
  state$log_lik <- mkp_loglikelihood(
    tmp_tree, mkd,
    kPrime = state$kPrime,
    rate_loss = state$rate_loss,
    rate_log_sd = state$rate_log_sd,
    nCat = model$nCat,
    coding = model$coding,
    relabel = model$relabel
  )
  state$log_prior <- log_prior(state, model, mkd)
  state$log_post <- state$log_lik + state$log_prior

  state
}


#' Build move schedule
#' @keywords internal
.build_moves <- function(nEdge, nTrans, has_neo, mcmc) {
  moves <- list(
    list(name = "tree_length", type = "scale", target = "tree_length",
         weight = 1),
    list(name = "branch_lengths", type = "beta_simplex",
         target = "rel_br_lengths", weight = max(1, nEdge / 3))
  )

  if (nTrans > 0) {
    moves <- c(moves, list(
      list(name = "kPrime", type = "int_walk", target = "kPrime",
           weight = max(1, 2 * nTrans)),
      list(name = "p", type = "scale", target = "p",
           weight = 1)
    ))
  }

  if (has_neo) {
    moves <- c(moves, list(
      list(name = "rate_loss", type = "scale", target = "rate_loss",
           weight = 1.5)
    ))
  }

  moves <- c(moves, list(
    list(name = "rate_log_sd", type = "scale", target = "rate_log_sd",
         weight = 1.5)
  ))

  moves
}


#' Execute a move and return accept/reject decision
#' @keywords internal
.do_move <- function(move, state, mkd, model, tree, tuning) {
  proposed <- state

  switch(move$type,
    scale = {
      prop <- propose_scale(state[[move$target]],
                            tuning = tuning[[paste0("scale_", move$target)]])
      proposed[[move$target]] <- prop$value
      log_hastings <- prop$log_hastings
    },
    beta_simplex = {
      prop <- propose_beta_simplex(state$rel_br_lengths,
                                    tuning = tuning$beta_simplex)
      proposed$rel_br_lengths <- prop$value
      log_hastings <- prop$log_hastings
    },
    int_walk = {
      # Pick a random transformational character
      trans_idx <- which(mkd$type == "transformational")
      char_i <- sample(trans_idx, 1L)
      prop <- propose_bounded_int_walk(
        state$kPrime[char_i],
        lower = mkd$kObs[char_i],
        window = tuning$int_walk_window
      )
      proposed$kPrime[char_i] <- prop$value
      log_hastings <- prop$log_hastings
    }
  )

  # Early rejection if Hastings ratio is -Inf
  if (!is.finite(log_hastings)) {
    return(list(accept = FALSE, state = state))
  }

  # Compute proposed prior
  proposed$log_prior <- log_prior(proposed, model, mkd)
  if (!is.finite(proposed$log_prior)) {
    return(list(accept = FALSE, state = state))
  }

  # Compute proposed likelihood
  edge_length <- proposed$tree_length * proposed$rel_br_lengths
  tmp_tree <- tree
  tmp_tree$edge.length <- edge_length
  proposed$log_lik <- mkp_loglikelihood(
    tmp_tree, mkd,
    kPrime = proposed$kPrime,
    rate_loss = proposed$rate_loss,
    rate_log_sd = proposed$rate_log_sd,
    nCat = model$nCat,
    coding = model$coding,
    relabel = model$relabel
  )
  proposed$log_post <- proposed$log_lik + proposed$log_prior

  # MH accept/reject
  log_alpha <- proposed$log_post - state$log_post + log_hastings
  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    list(accept = TRUE, state = proposed)
  } else {
    list(accept = FALSE, state = state)
  }
}


#' Parameter names for the sample matrix
#' @keywords internal
.param_names <- function(mkd, nEdge) {
  names <- c("log_posterior", "log_likelihood", "tree_length",
             "rate_loss", "rate_log_sd", "p")

  trans_idx <- which(mkd$type == "transformational")
  if (length(trans_idx)) {
    names <- c(names, paste0("kPrime_", trans_idx))
  }

  # Branch length proportions (abbreviated)
  names <- c(names, paste0("br_", seq_len(nEdge)))

  names
}


#' Extract state values to a row vector for storage
#' @keywords internal
.state_to_row <- function(state, mkd, nEdge) {
  trans_idx <- which(mkd$type == "transformational")
  kp <- if (length(trans_idx)) as.numeric(state$kPrime[trans_idx]) else numeric(0)

  c(state$log_post, state$log_lik, state$tree_length,
    state$rate_loss, state$rate_log_sd, state$p,
    kp,
    state$rel_br_lengths)
}


#' Reconstruct a phylo object from current state
#' @keywords internal
.state_to_tree <- function(state, tree) {
  t <- tree
  t$edge.length <- state$tree_length * state$rel_br_lengths
  t
}


#' Adapt tuning parameters based on acceptance rates
#' @keywords internal
.adapt_tuning <- function(tuning, accept_count, propose_count, moves) {
  targets <- c(
    tree_length = 0.35, branch_lengths = 0.23, kPrime = 0.35,
    p = 0.35, rate_loss = 0.35, rate_log_sd = 0.35
  )

  tuning_keys <- c(
    tree_length = "scale_tree_length",
    branch_lengths = "beta_simplex",
    kPrime = "int_walk_window",
    p = "scale_p",
    rate_loss = "scale_rate_loss",
    rate_log_sd = "scale_rate_log_sd"
  )

  for (move in moves) {
    nm <- move$name
    if (propose_count[nm] < 20) next
    rate <- accept_count[nm] / propose_count[nm]
    target <- targets[nm]
    tk <- tuning_keys[nm]

    if (!is.na(tk) && !is.null(tuning[[tk]])) {
      # Multiplicative adjustment
      adj <- exp(0.5 * (rate - target))
      if (nm == "kPrime") {
        # Integer window: adjust multiplicatively, round, clamp to >= 1
        tuning[[tk]] <- max(1L, as.integer(round(tuning[[tk]] * adj)))
      } else if (nm == "branch_lengths") {
        # BetaSimplex: higher tuning = more conservative; invert direction
        tuning[[tk]] <- tuning[[tk]] / adj
        tuning[[tk]] <- max(2, tuning[[tk]])
      } else {
        tuning[[tk]] <- tuning[[tk]] * adj
        tuning[[tk]] <- max(0.01, tuning[[tk]])
      }
    }
  }

  # Reset counters implicitly — caller should track
  tuning
}
