# Main MCMC entry point for MkPrime
#
# Phase 3: single chain, fixed topology, R-side loop.
# Phase 4: topology moves (NNI, SPR) via mutable state$tree.
# Phase 5: parallel tempering (multiple chains, temperature ladder).

#' Run Bayesian MCMC under the MkPrime model
#'
#' Metropolis-Hastings MCMC sampling tree topology, branch lengths, model
#' parameters, and per-character k' (for transformational characters).
#' Supports parallel tempering with a geometric temperature ladder.
#'
#' @param data A `phyDat` object or `MkPrimeData` object.
#' @param tree A `phylo` object (starting topology).
#' @param neomorphic,known_states Passed to [MkPrimeData()] if `data` is
#'   a `phyDat` object.
#' @param model An `MkPrimeModel` object, or `NULL` for defaults.
#' @param mcmc An `MkPrimeMCMC` object, or `NULL` for defaults.
#' @param fix_topology Logical. If `TRUE`, tree topology is fixed (Phase 3
#'   behaviour). Default `FALSE` enables NNI and SPR topology proposals.
#'
#' @return An `MkPosterior` object.
#' @export
RunMkPrime <- function(data, tree,
                       neomorphic = integer(0),
                       known_states = integer(0),
                       model = NULL,
                       mcmc = NULL,
                       fix_topology = FALSE) {

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

  # Precompute tree structure
  tree <- ape::reorder.phylo(tree, "postorder")
  nEdge <- nrow(tree$edge)

  has_neo <- any(mkd$type == "neomorphic")
  trans_idx <- which(mkd$type == "transformational")
  has_trans <- length(trans_idx) > 0
  nTrans <- length(trans_idx)

  nChains <- mcmc$nChains

  # --- Temperature ladder ---
  betas <- .build_temperature_ladder(nChains, mcmc$heat)

  # --- Initialize per-chain states ---
  chains <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chains[[ch]] <- .init_state(tree, mkd, model)
  }

  # --- Build move schedule (shared across chains) ---
  moves <- .build_moves(nEdge, nTrans, has_neo, mcmc,
                        fix_topology = fix_topology)

  # --- Per-chain acceptance tracking and tuning ---
  chain_accept <- vector("list", nChains)
  chain_propose <- vector("list", nChains)
  chain_tuning <- vector("list", nChains)
  move_names <- vapply(moves, `[[`, character(1), "name")
  for (ch in seq_len(nChains)) {
    chain_accept[[ch]] <- integer(length(moves))
    chain_propose[[ch]] <- integer(length(moves))
    names(chain_accept[[ch]]) <- names(chain_propose[[ch]]) <- move_names
    chain_tuning[[ch]] <- mcmc$tuning
  }

  # --- Pre-allocate sample storage (cold chain only) ---
  nSaved <- as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  param_names <- .param_names(mkd, nEdge)
  samples <- matrix(NA_real_, nrow = nSaved, ncol = length(param_names),
                    dimnames = list(NULL, param_names))
  tree_samples <- vector("list", nSaved)
  saved_idx <- 0L

  # Tree file logging
  tree_file <- mcmc$tree_file
  if (!is.null(tree_file)) {
    writeLines("", tree_file)
  }

  # --- Swap tracking ---
  if (nChains > 1L) {
    swap_accept <- integer(nChains - 1L)
    swap_propose <- integer(nChains - 1L)
  }

  # --- MCMC loop ---
  cli::cli_progress_bar(
    "MCMC", total = mcmc$nIter,
    format = paste0(
      "{cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
      " | logP: {format(round(cold_logpost, 1), nsmall = 1)}",
      " | acc: {format(round(recent_acc, 2), nsmall = 2)}",
      if (nChains > 1L) " | chains: {nChains}" else ""
    )
  )
  cold_logpost <- chains[[1]]$log_lik + chains[[1]]$log_prior
  recent_acc <- 0
  recent_window <- 500L
  recent_accepts <- logical(recent_window)
  recent_pos <- 0L

  for (iter in seq_len(mcmc$nIter)) {
    # --- Propose moves for each chain ---
    for (ch in seq_len(nChains)) {
      move_idx <- sample.int(length(moves), 1L, prob = vapply(
        moves, `[[`, numeric(1), "weight"
      ))
      move <- moves[[move_idx]]
      chain_propose[[ch]][move_idx] <- chain_propose[[ch]][move_idx] + 1L

      accepted <- .do_move(move, chains[[ch]], mkd, model,
                           chain_tuning[[ch]], beta = betas[ch])

      if (accepted$accept) {
        chains[[ch]] <- accepted$state
        chain_accept[[ch]][move_idx] <- chain_accept[[ch]][move_idx] + 1L
      }

      # Track recent acceptance for cold chain only
      if (ch == 1L) {
        recent_pos <- (recent_pos %% recent_window) + 1L
        recent_accepts[recent_pos] <- accepted$accept
      }
    }

    # --- Chain swap proposal (M-030 placeholder) ---
    if (nChains > 1L) {
      swap_result <- .propose_chain_swap(chains, betas)
      chains <- swap_result$chains
      if (!is.null(swap_result$pair)) {
        pair_idx <- swap_result$pair[1]
        swap_propose[pair_idx] <- swap_propose[pair_idx] + 1L
        if (swap_result$accepted) {
          swap_accept[pair_idx] <- swap_accept[pair_idx] + 1L
        }
      }
    }

    # --- Adaptation during warmup ---
    if (iter <= mcmc$warmup && iter %% 200L == 0L) {
      # Adapt proposal tuning (per-chain)
      for (ch in seq_len(nChains)) {
        chain_tuning[[ch]] <- .adapt_tuning(
          chain_tuning[[ch]], chain_accept[[ch]],
          chain_propose[[ch]], moves
        )
      }
      # Adapt temperature ladder
      if (nChains > 1L) {
        betas <- .adapt_temperatures(betas, swap_accept, swap_propose)
      }
    }

    # --- Save cold chain samples post-warmup ---
    if (iter > mcmc$warmup && (iter - mcmc$warmup) %% mcmc$thin == 0L) {
      saved_idx <- saved_idx + 1L
      samples[saved_idx, ] <- .state_to_row(chains[[1]], mkd, nEdge)
      cur_tree <- .state_to_tree(chains[[1]])
      tree_samples[[saved_idx]] <- cur_tree
      if (!is.null(tree_file)) {
        cat(ape::write.tree(cur_tree), "\n", file = tree_file, append = TRUE)
      }
    }

    if (iter %% 100L == 0L) {
      n_recent <- min(iter, recent_window)
      recent_acc <- sum(recent_accepts[seq_len(n_recent)]) / n_recent
      cold_logpost <- chains[[1]]$log_lik + chains[[1]]$log_prior
      cli::cli_progress_update()
    }
  }
  cli::cli_progress_done()

  # --- Build result ---
  cold_acceptance <- chain_accept[[1]] / pmax(chain_propose[[1]], 1L)

  result <- MkPosterior(
    samples = samples,
    trees = tree_samples,
    acceptance = cold_acceptance,
    model = model,
    data = mkd,
    mcmc = mcmc,
    warmup = mcmc$warmup,
    tuning = chain_tuning[[1]]
  )

  # Attach tempering info
  if (nChains > 1L) {
    result$betas <- betas
    result$swap_rates <- swap_accept / pmax(swap_propose, 1L)
    result$chain_acceptance <- lapply(seq_len(nChains), function(ch) {
      chain_accept[[ch]] / pmax(chain_propose[[ch]], 1L)
    })
  }

  result
}


# --- Internal helpers ---

#' Build geometric temperature ladder
#' @keywords internal
.build_temperature_ladder <- function(nChains, heat) {
  if (nChains == 1L) return(1.0)
  # beta_1 = 1 (cold), beta_nChains = heat
  # beta_i = heat^((i-1) / (nChains-1))
  heat^(seq(0, 1, length.out = nChains))
}


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
    tree = tree,
    tree_length = tree_length,
    rel_br_lengths = rel_br,
    rate_loss = 1.0,
    rate_log_sd = 0.5,
    kPrime = as.integer(kPrime),
    p = 0.5
  )

  # Compute initial likelihood and prior (unheated)
  state$log_lik <- mkp_loglikelihood(
    tree, mkd,
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
.build_moves <- function(nEdge, nTrans, has_neo, mcmc,
                         fix_topology = FALSE) {
  moves <- list(
    list(name = "tree_length", type = "scale", target = "tree_length",
         weight = 1),
    list(name = "branch_lengths", type = "beta_simplex",
         target = "rel_br_lengths", weight = max(1, nEdge / 3))
  )

  if (!fix_topology && nEdge >= 5L) {
    moves <- c(moves, list(
      list(name = "nni", type = "nni", target = NULL,
           weight = max(1, nEdge / 2)),
      list(name = "spr", type = "spr", target = NULL,
           weight = max(1, nEdge / 4))
    ))
  }

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
#'
#' @param beta Inverse temperature (1 = cold chain, < 1 = heated).
#'   The heated MH acceptance is:
#'   `log_alpha = beta * (logLik_new - logLik_old) +
#'                (logPrior_new - logPrior_old) + log_hastings`
#' @keywords internal
.do_move <- function(move, state, mkd, model, tuning, beta = 1.0) {
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
    nni = {
      prop <- propose_nni(state$tree, state$tree_length,
                          state$rel_br_lengths)
      proposed$tree <- prop$tree
      proposed$rel_br_lengths <- prop$rel_br_lengths
      log_hastings <- prop$log_hastings
    },
    spr = {
      prop <- propose_spr(state$tree, state$tree_length,
                          state$rel_br_lengths)
      proposed$tree <- prop$tree
      proposed$rel_br_lengths <- prop$rel_br_lengths
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

  # Compute proposed prior (unheated)
  proposed$log_prior <- log_prior(proposed, model, mkd)
  if (!is.finite(proposed$log_prior)) {
    return(list(accept = FALSE, state = state))
  }

  # Compute proposed likelihood (unheated)
  tmp_tree <- proposed$tree
  tmp_tree$edge.length <- proposed$tree_length * proposed$rel_br_lengths
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

  # Heated MH acceptance: only likelihood is tempered

  log_alpha <- beta * (proposed$log_lik - state$log_lik) +
               (proposed$log_prior - state$log_prior) +
               log_hastings

  if (is.finite(log_alpha) && log(runif(1)) < log_alpha) {
    list(accept = TRUE, state = proposed)
  } else {
    list(accept = FALSE, state = state)
  }
}


#' Propose a swap between two adjacent chains
#'
#' Selects a random adjacent pair and proposes exchanging their full states.
#' @return List with `chains` (possibly swapped), `pair` (integer vector of
#'   length 2 giving the lower index), and `accepted` (logical).
#' @keywords internal
.propose_chain_swap <- function(chains, betas) {
  nChains <- length(chains)
  if (nChains < 2L) {
    return(list(chains = chains, pair = NULL, accepted = FALSE))
  }

  # Pick random adjacent pair
  i <- sample.int(nChains - 1L, 1L)
  j <- i + 1L

  # Swap acceptance: exp((beta_i - beta_j) * (logLik_j - logLik_i))
  log_alpha <- (betas[i] - betas[j]) *
               (chains[[j]]$log_lik - chains[[i]]$log_lik)

  accepted <- is.finite(log_alpha) && log(runif(1)) < log_alpha
  if (accepted) {
    tmp <- chains[[i]]
    chains[[i]] <- chains[[j]]
    chains[[j]] <- tmp
  }

  list(chains = chains, pair = c(i, j), accepted = accepted)
}


#' Adapt temperature ladder based on swap acceptance rates
#'
#' Adjusts the `heat` parameter (temperature of the hottest chain) to
#' achieve ~25% swap acceptance between adjacent pairs. The ladder is
#' then reconstructed with geometric spacing.
#'
#' @param betas Current temperature ladder.
#' @param swap_accept Integer vector of swap acceptances per adjacent pair.
#' @param swap_propose Integer vector of swap proposals per adjacent pair.
#' @param target Target swap acceptance rate (default 0.25).
#' @return Updated temperature ladder.
#' @keywords internal
.adapt_temperatures <- function(betas, swap_accept, swap_propose,
                                target = 0.25) {
  nChains <- length(betas)
  if (nChains < 2L) return(betas)

  total_propose <- sum(swap_propose)
  if (total_propose < 20L) return(betas)

  overall_rate <- sum(swap_accept) / total_propose

  # Adjust heat: if swap rate too low, increase heat (bring temps closer);
  # if too high, decrease heat (spread temps further)
  heat <- betas[nChains]
  # heat^adj where adj < 1 → increases heat (closer temps, easier swaps)
  # and adj > 1 → decreases heat (wider spread, harder swaps)
  adj <- exp(0.5 * (overall_rate - target))
  heat_new <- heat^adj

  # Clamp to reasonable range
  heat_new <- max(0.01, min(0.95, heat_new))

  .build_temperature_ladder(nChains, heat_new)
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
.state_to_tree <- function(state) {
  t <- state$tree
  t$edge.length <- state$tree_length * state$rel_br_lengths
  t
}


#' Adapt tuning parameters based on acceptance rates
#' @keywords internal
.adapt_tuning <- function(tuning, accept_count, propose_count, moves) {
  targets <- c(
    tree_length = 0.35, branch_lengths = 0.23,
    nni = 0.23, spr = 0.10,
    kPrime = 0.35,
    p = 0.35, rate_loss = 0.35, rate_log_sd = 0.35
  )

  tuning_keys <- c(
    tree_length = "scale_tree_length",
    branch_lengths = "beta_simplex",
    nni = NA_character_,
    spr = NA_character_,
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

  tuning
}
