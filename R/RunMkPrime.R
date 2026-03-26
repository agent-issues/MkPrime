# Main MCMC entry point for MkPrime
#
# Phase 3: single chain, fixed topology, R-side loop.
# Phase 4: topology moves (NNI, SPR) via mutable state$tree.
# Phase 5: parallel tempering, independent runs, convergence, stopping.

#' Run Bayesian MCMC under the MkPrime model
#'
#' Metropolis-Hastings MCMC sampling tree topology, branch lengths, model
#' parameters, and per-character k' (for transformational characters).
#' Supports parallel tempering with a geometric temperature ladder,
#' multiple independent runs, convergence monitoring, and early stopping.
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

  tree <- ape::reorder.phylo(tree, "postorder")
  nEdge <- nrow(tree$edge)

  has_neo <- any(mkd$type == "neomorphic")
  trans_idx <- which(mkd$type == "transformational")
  nTrans <- length(trans_idx)

  moves <- .build_moves(nEdge, nTrans, has_neo, mcmc,
                        fix_topology = fix_topology)

  nRuns <- mcmc$nRuns

  # --- Initialize per-run state ---
  runs <- vector("list", nRuns)
  for (run in seq_len(nRuns)) {
    start_tree <- if (run == 1L) tree else .perturb_start(tree)
    runs[[run]] <- .init_run(start_tree, mkd, model, mcmc, moves)
  }

  # --- Interleaved MCMC loop ---
  nSaved_per_run <- as.integer((mcmc$nIter - mcmc$warmup) / mcmc$thin)
  param_names <- .param_names(mkd, nEdge)

  # Pre-allocate storage per run
  for (run in seq_len(nRuns)) {
    runs[[run]]$samples <- matrix(NA_real_, nrow = nSaved_per_run,
                                  ncol = length(param_names),
                                  dimnames = list(NULL, param_names))
    runs[[run]]$tree_samples <- vector("list", nSaved_per_run)
    runs[[run]]$saved_idx <- 0L
  }

  # Tree file
  tree_file <- mcmc$tree_file
  if (!is.null(tree_file)) writeLines("", tree_file)

  # Stopping state
  stop_reason <- "max_iter"
  start_time <- proc.time()["elapsed"]

  # Progress bar
  cli::cli_progress_bar(
    "MCMC", total = mcmc$nIter,
    format = paste0(
      "{cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
      " | logP: {format(round(cold_logpost, 1), nsmall = 1)}",
      " | acc: {format(round(recent_acc, 2), nsmall = 2)}",
      if (nRuns > 1L) " | runs: {nRuns}" else "",
      if (mcmc$nChains > 1L) " | chains: {mcmc$nChains}" else ""
    )
  )
  cold_logpost <- runs[[1]]$chains[[1]]$log_lik +
                  runs[[1]]$chains[[1]]$log_prior
  recent_acc <- 0
  recent_window <- 500L
  recent_accepts <- logical(recent_window)
  recent_pos <- 0L

  for (iter in seq_len(mcmc$nIter)) {
    # --- Advance all runs by one iteration ---
    for (run in seq_len(nRuns)) {
      r <- runs[[run]]
      nChains <- mcmc$nChains

      # Propose moves for each chain
      for (ch in seq_len(nChains)) {
        move_idx <- sample.int(length(moves), 1L, prob = vapply(
          moves, `[[`, numeric(1), "weight"
        ))
        move <- moves[[move_idx]]
        r$chain_propose[[ch]][move_idx] <-
          r$chain_propose[[ch]][move_idx] + 1L

        accepted <- .do_move(move, r$chains[[ch]], mkd, model,
                             r$chain_tuning[[ch]], beta = r$betas[ch])

        if (accepted$accept) {
          r$chains[[ch]] <- accepted$state
          r$chain_accept[[ch]][move_idx] <-
            r$chain_accept[[ch]][move_idx] + 1L
        }

        # Track recent acceptance for run 1, cold chain
        if (run == 1L && ch == 1L) {
          recent_pos <- (recent_pos %% recent_window) + 1L
          recent_accepts[recent_pos] <- accepted$accept
        }
      }

      # Chain swaps
      if (nChains > 1L) {
        swap_result <- .propose_chain_swap(r$chains, r$betas)
        r$chains <- swap_result$chains
        if (!is.null(swap_result$pair)) {
          pair_idx <- swap_result$pair[1]
          r$swap_propose[pair_idx] <- r$swap_propose[pair_idx] + 1L
          if (swap_result$accepted) {
            r$swap_accept[pair_idx] <- r$swap_accept[pair_idx] + 1L
          }
        }
      }

      # Adaptation during warmup
      if (iter <= mcmc$warmup && iter %% 200L == 0L) {
        for (ch in seq_len(nChains)) {
          r$chain_tuning[[ch]] <- .adapt_tuning(
            r$chain_tuning[[ch]], r$chain_accept[[ch]],
            r$chain_propose[[ch]], moves
          )
        }
        if (nChains > 1L) {
          r$betas <- .adapt_temperatures(r$betas, r$swap_accept,
                                         r$swap_propose)
        }
      }

      # Save cold chain samples post-warmup
      if (iter > mcmc$warmup && (iter - mcmc$warmup) %% mcmc$thin == 0L) {
        r$saved_idx <- r$saved_idx + 1L
        r$samples[r$saved_idx, ] <- .state_to_row(r$chains[[1]], mkd, nEdge)
        cur_tree <- .state_to_tree(r$chains[[1]])
        r$tree_samples[[r$saved_idx]] <- cur_tree
        if (!is.null(tree_file)) {
          cat(ape::write.tree(cur_tree), "\n", file = tree_file,
              append = TRUE)
        }
      }

      runs[[run]] <- r
    }

    # Progress
    if (iter %% 100L == 0L) {
      n_recent <- min(iter, recent_window)
      recent_acc <- sum(recent_accepts[seq_len(n_recent)]) / n_recent
      cold_logpost <- runs[[1]]$chains[[1]]$log_lik +
                      runs[[1]]$chains[[1]]$log_prior
      cli::cli_progress_update()
    }

    # --- Stopping rule checks ---
    if (!is.null(mcmc$max_time)) {
      elapsed <- proc.time()["elapsed"] - start_time
      if (elapsed >= mcmc$max_time) {
        stop_reason <- "max_time"
        break
      }
    }

    if (iter > mcmc$warmup && !is.null(mcmc$check_every) &&
        iter %% mcmc$check_every == 0L && nRuns >= 2L) {
      diag <- .check_convergence(runs, param_names, mcmc)
      if (!is.null(diag) && diag$converged) {
        stop_reason <- "converged"
        break
      }
    }
  }
  cli::cli_progress_done()

  actual_iter <- min(iter, mcmc$nIter)

  # --- Build result ---
  .build_result(runs, model, mkd, mcmc, actual_iter, stop_reason)
}


# --- Run initialization ---

#' Initialize state for a single run
#' @keywords internal
.init_run <- function(tree, mkd, model, mcmc, moves) {
  nChains <- mcmc$nChains
  betas <- .build_temperature_ladder(nChains, mcmc$heat)

  chains <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chains[[ch]] <- .init_state(tree, mkd, model)
  }

  move_names <- vapply(moves, `[[`, character(1), "name")
  chain_accept <- chain_propose <- chain_tuning <- vector("list", nChains)
  for (ch in seq_len(nChains)) {
    chain_accept[[ch]] <- integer(length(moves))
    chain_propose[[ch]] <- integer(length(moves))
    names(chain_accept[[ch]]) <- names(chain_propose[[ch]]) <- move_names
    chain_tuning[[ch]] <- mcmc$tuning
  }

  swap_accept <- swap_propose <- if (nChains > 1L) {
    integer(nChains - 1L)
  } else {
    integer(0)
  }

  list(
    chains = chains,
    betas = betas,
    chain_accept = chain_accept,
    chain_propose = chain_propose,
    chain_tuning = chain_tuning,
    swap_accept = swap_accept,
    swap_propose = swap_propose
  )
}


# --- Convergence check during MCMC ---

#' Check convergence criteria (called during the loop)
#' @keywords internal
.check_convergence <- function(runs, param_names, mcmc) {
  if (!requireNamespace("coda", quietly = TRUE)) return(NULL)

  key_cols <- .key_param_cols(
    matrix(0, 1, length(param_names), dimnames = list(NULL, param_names))
  )

  # Gather saved samples from each run
  per_run_samples <- lapply(runs, function(r) {
    idx <- r$saved_idx
    if (idx < 10L) return(NULL)
    r$samples[seq_len(idx), key_cols, drop = FALSE]
  })

  if (any(vapply(per_run_samples, is.null, logical(1)))) return(NULL)

  # Compute PSRF
  chain_list <- lapply(per_run_samples, function(s) coda::mcmc(s))
  mcmc_list <- coda::mcmc.list(chain_list)

  gd <- tryCatch(
    coda::gelman.diag(mcmc_list, multivariate = FALSE),
    error = function(e) NULL
  )
  if (is.null(gd)) return(NULL)

  max_psrf <- max(gd$psrf[, 1], na.rm = TRUE)

  # Compute min ESS across all runs combined
  combined <- do.call(rbind, per_run_samples)
  ess <- apply(combined, 2, function(col) {
    if (sd(col, na.rm = TRUE) == 0) return(NA_real_)
    coda::effectiveSize(coda::mcmc(col))
  })
  min_ess <- min(ess, na.rm = TRUE)

  converged <- TRUE
  if (!is.null(mcmc$min_ess) && min_ess < mcmc$min_ess) {
    converged <- FALSE
  }
  if (!is.null(mcmc$max_psrf) && max_psrf > mcmc$max_psrf) {
    converged <- FALSE
  }

  list(converged = converged, min_ess = min_ess, max_psrf = max_psrf)
}


# --- Build final result ---

#' Build MkPosterior from all runs
#' @keywords internal
.build_result <- function(runs, model, mkd, mcmc, actual_iter, stop_reason) {
  nRuns <- length(runs)

  # Trim samples to actual saved count
  for (run in seq_len(nRuns)) {
    idx <- runs[[run]]$saved_idx
    if (idx > 0L) {
      runs[[run]]$samples <- runs[[run]]$samples[seq_len(idx), , drop = FALSE]
      runs[[run]]$tree_samples <- runs[[run]]$tree_samples[seq_len(idx)]
    } else {
      runs[[run]]$samples <- runs[[run]]$samples[integer(0), , drop = FALSE]
      runs[[run]]$tree_samples <- list()
    }
  }

  # Per-run summaries
  per_run_summaries <- lapply(runs, function(r) {
    cold_acc <- r$chain_accept[[1]] / pmax(r$chain_propose[[1]], 1L)
    result <- list(
      samples = r$samples,
      trees = r$tree_samples,
      acceptance = cold_acc
    )
    if (mcmc$nChains > 1L) {
      result$betas <- r$betas
      result$swap_rates <- r$swap_accept / pmax(r$swap_propose, 1L)
    }
    result
  })

  if (nRuns == 1L) {
    r <- per_run_summaries[[1]]
    result <- MkPosterior(
      samples = r$samples,
      trees = r$trees,
      acceptance = r$acceptance,
      model = model,
      data = mkd,
      mcmc = mcmc,
      warmup = mcmc$warmup,
      tuning = runs[[1]]$chain_tuning[[1]]
    )
    if (!is.null(r$betas)) {
      result$betas <- r$betas
      result$swap_rates <- r$swap_rates
      result$chain_acceptance <- lapply(seq_len(mcmc$nChains), function(ch) {
        runs[[1]]$chain_accept[[ch]] /
          pmax(runs[[1]]$chain_propose[[ch]], 1L)
      })
    }
  } else {
    all_samples <- do.call(rbind, lapply(per_run_summaries, `[[`, "samples"))
    all_trees <- do.call(c, lapply(per_run_summaries, `[[`, "trees"))
    avg_acceptance <- Reduce(`+`, lapply(per_run_summaries, `[[`,
                                         "acceptance")) / nRuns

    result <- MkPosterior(
      samples = all_samples,
      trees = all_trees,
      acceptance = avg_acceptance,
      model = model,
      data = mkd,
      mcmc = mcmc,
      warmup = mcmc$warmup,
      tuning = runs[[1]]$chain_tuning[[1]]
    )

    result$nRuns <- nRuns
    result$per_run <- per_run_summaries

    if (mcmc$nChains > 1L) {
      result$betas <- runs[[1]]$betas
      result$swap_rates <- per_run_summaries[[1]]$swap_rates
      result$chain_acceptance <- lapply(seq_len(mcmc$nChains), function(ch) {
        runs[[1]]$chain_accept[[ch]] /
          pmax(runs[[1]]$chain_propose[[ch]], 1L)
      })
    }
  }

  result$stop_reason <- stop_reason
  result$actual_iter <- actual_iter
  result
}


# --- Start perturbation ---

#' Generate a perturbed starting tree for independent runs
#' @keywords internal
.perturb_start <- function(tree) {
  nTip <- length(tree$tip.label)
  if (nTip < 4L) return(tree)

  n_nni <- sample(2:5, 1)
  for (i in seq_len(n_nni)) {
    tree_length <- sum(tree$edge.length)
    rel_br <- tree$edge.length / tree_length
    prop <- propose_nni(tree, tree_length, rel_br)
    tree <- prop$tree
    tree$edge.length <- tree_length * prop$rel_br_lengths
  }

  tree$edge.length <- tree$edge.length *
    exp(rnorm(length(tree$edge.length), sd = 0.1))
  tree
}


# --- Internal helpers ---

#' Build geometric temperature ladder
#' @keywords internal
.build_temperature_ladder <- function(nChains, heat) {
  if (nChains == 1L) return(1.0)
  heat^(seq(0, 1, length.out = nChains))
}


#' Initialize MCMC state from tree and data
#' @keywords internal
.init_state <- function(tree, mkd, model) {
  tree_length <- sum(tree$edge.length)
  rel_br <- tree$edge.length / tree_length

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

  if (!is.finite(log_hastings)) {
    return(list(accept = FALSE, state = state))
  }

  proposed$log_prior <- log_prior(proposed, model, mkd)
  if (!is.finite(proposed$log_prior)) {
    return(list(accept = FALSE, state = state))
  }

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
#' @keywords internal
.propose_chain_swap <- function(chains, betas) {
  nChains <- length(chains)
  if (nChains < 2L) {
    return(list(chains = chains, pair = NULL, accepted = FALSE))
  }

  i <- sample.int(nChains - 1L, 1L)
  j <- i + 1L

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
#' @keywords internal
.adapt_temperatures <- function(betas, swap_accept, swap_propose,
                                target = 0.25) {
  nChains <- length(betas)
  if (nChains < 2L) return(betas)

  total_propose <- sum(swap_propose)
  if (total_propose < 20L) return(betas)

  overall_rate <- sum(swap_accept) / total_propose

  heat <- betas[nChains]
  adj <- exp(0.5 * (overall_rate - target))
  heat_new <- heat^adj

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
      adj <- exp(0.5 * (rate - target))
      if (nm == "kPrime") {
        tuning[[tk]] <- max(1L, as.integer(round(tuning[[tk]] * adj)))
      } else if (nm == "branch_lengths") {
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
