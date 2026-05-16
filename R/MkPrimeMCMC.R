# MCMC configuration for MkPrime

#' Configure MCMC settings
#'
#' @param nIter Total iterations (including warmup). Default `Inf`, which
#'   relies on stopping criteria (`minEss`, `maxRhat`, `maxTime`) to
#'   terminate the run. Pass a finite integer to cap the number of iterations
#'   regardless of convergence.
#' @param thin Thinning interval (save every `thin`-th iteration).
#'   Default `"auto"`: initially set to the number of active move types
#'   (one "full cycle" per stored sample), then adapted at the first
#'   convergence check to `max(nMoves, round(maxACT * log(2)))` based on
#'   the observed autocorrelation time of the worst-mixing scalar
#'   parameter. This targets ~50\% correlation between consecutive
#'   stored samples. Pass a positive integer to fix the thinning interval
#'   and disable adaptation.
#' @param treeThin Thinning interval for tree samples (iterations).
#'   Default `NULL`, which stores a tree for every scalar sample
#'   (equivalent to `treeThin = thin`). Pass a positive integer that is
#'   a **multiple of `thin`** to store trees less frequently than scalar
#'   parameters. For example, `thin = 10, treeThin = 100` stores one tree
#'   for every 10 scalar samples. The multiple-of-thin constraint is
#'   validated at run time (after `thin = "auto"` is resolved).
#' @param warmup **Deprecated.** If supplied, treated as `maxWarmup`.
#'   Use `minWarmup` / `maxWarmup` instead. Retained for backward
#'   compatibility; a deprecation message is emitted when non-NULL.
#' @param minWarmup Minimum warmup iterations before the stabilisation
#'   detector can end warmup. Default 2000. Set lower for quick
#'   debugging runs.
#' @param maxWarmup Maximum warmup iterations. Warmup ends when the
#'   chain stabilises OR this ceiling is reached (with a warning).
#'   Default 50 000 when `nIter = Inf`; `nIter / 2` when finite.
#' @param autoTune Logical; enable the tuning phase after warmup
#'   (default `TRUE`). When `TRUE`, a short tuning phase optimises
#'   move weights by maximising min-ESS/s using a perturbation bandit.
#'   When `FALSE`, move weights are frozen at the end of warmup using
#'   the acceptance-rate heuristic (legacy behaviour). See section
#'   **Three-phase MCMC** below.
#' @param tuningBudget Maximum iterations in the tuning phase.
#'   Default 10 000. Only used when `autoTune = TRUE`.
#' @param tuningRounds Number of perturbation-evaluation rounds in
#'   the tuning phase. Default 5. Each round evaluates the current
#'   weights plus `nPerturbations` candidates. Only used when
#'   `autoTune = TRUE`.
#' @param nRuns Number of independent runs. Default 2. Each run has its
#'   own set of `nChains` chains. Convergence diagnostics (R-hat) require
#'   `nRuns >= 2`.
#' @param nChains Number of chains in the temperature ladder (parallel
#'   tempering). Default 1 (no tempering). Set to 4 for typical analyses.
#'   Chain 1 is the cold chain (beta = 1).
#' @param heat Temperature of the hottest chain (0 < heat < 1). Default
#'   0.2. The temperature ladder uses geometric spacing:
#'   `beta_i = heat^((i-1)/(nChains-1))` for i = 1, ..., nChains.
#'   Ignored when `nChains = 1`.
#' @param maxTime Maximum wall-clock time in seconds. `NULL` (default)
#'   means no time limit.
#' @param minEss Minimum effective sample size for early stopping.
#'   `NULL` (default) disables ESS-based stopping.
#' @param maxRhat Maximum R-hat (rank-normalized; Vehtari et al. 2021)
#'   for early stopping. `NULL` (default) disables R-hat-based stopping.
#'   The modern recommendation is 1.01 for reliable inference; 1.05 is
#'   a pragmatic threshold for phylogenetics where mixing is slower.
#'   Requires `nRuns >= 2`.
#'   Replaces the classical PSRF (Gelman-Rubin) statistic.
#' @param minTreeEss Minimum tree-topology ESS (median pseudo-ESS) for
#'   early stopping.  `NULL` (default) disables tree-ESS-based stopping.
#'   When set, tree ESS is computed adaptively during convergence checks:
#'   skipped when scalar ESS is far from `minEss`, coarse (500-tree
#'   subsample) when approaching, and fine (1000-tree subsample) when
#'   tree ESS is the binding constraint.  Requires **TreeDist**.
#'   Also used in the tuning-phase bandit: topology moves receive credit
#'   for improving tree ESS, preventing underallocation.
#' @param checkEvery Check convergence every this many iterations
#'   (default 1000). Only used when stopping criteria are set.
#' @param cancelFile Path to a cancel-signal file. `NULL` (default) disables
#'   cancel-file checking. When set, [RunMkPrime()] checks every 200
#'   iterations whether this file exists. If it does, the run flushes any
#'   buffered samples, saves a checkpoint (if `checkpointFile` is set), and
#'   exits with `stop_reason = "cancelled"`. Create the file to request a
#'   clean stop: `file.create(cancelFile)`. See also [MkCancelPath()].
#' @param checkpointFile Path to write checkpoint RDS files. `NULL`
#'   (default) auto-derives from `logFile` when set
#'   (e.g. `"run.log"` \u2192 `"run.ckp"`). Set to `FALSE` to
#'   disable checkpointing. Checkpoints are saved at each
#'   convergence check interval and on cancel.
#' @param treeFile Path to write sampled trees in Newick format.
#'   `NULL` (default) auto-derives from `logFile` when set
#'   (e.g. `"run.log"` → `"run_trees.nwk"`). Set to `FALSE` to
#'   disable tree file logging entirely. Trees are always stored
#'   in the returned `MkPosterior` object regardless.
#' @param logFile Path to write a tab-separated scalar parameter log
#'   (Tracer-compatible). `NULL` (default) keeps all samples in memory only.
#'   If the path has no file extension, `.log` is appended automatically.
#'   When set, samples are flushed to disk in batches of `bufferSize`, so the
#'   file can be opened in Tracer during the run.
#'   For multiple runs (`nRuns > 1`), separate files are created automatically
#'   by appending `_1`, `_2`, \dots before the file extension
#'   (e.g. `"run.log"` → `"run_1.log"`, `"run_2.log"`).
#'   See [ReadMkLog()] to load the log back into R after the run.
#' @param bufferSize Integer. Number of thinned samples to accumulate per run
#'   before flushing to `logFile`. Ignored when `logFile = NULL`. Default 500.
#'   Smaller values give more frequent flushes (lower Tracer latency, slightly
#'   more I/O); larger values amortise write overhead at the cost of slightly
#'   higher peak memory.
#' @param plotEvery Integer; invoke the progress callback every this many
#'   iterations. In interactive sessions, defaults to `checkEvery` so that
#'   live trace plots update alongside convergence checks. Set to `NULL`
#'   explicitly to disable. Typical values: 500--2000.
#' @param progressFn Progress callback function, the string `"default"`,
#'   or `NULL`. In interactive sessions, defaults to [MkpTracePlot()] for
#'   live base-R trace plots. In non-interactive sessions (e.g. `Rscript`),
#'   defaults to `NULL` (no plotting). Setting `plotEvery` without
#'   `progressFn` also implies `MkpTracePlot`. A custom function must
#'   accept a single argument: a named list with fields `iter`, `nIter`,
#'   `warmup`, `inWarmup`, `nRuns`, `nChains`, `runSamples`,
#'   `currentState`, `recentAcceptance`, and `elapsed`. See
#'   [MkpTracePlot()] for details.
#' @param gibbsSpr Logical; include the Gibbs SPR move (default `TRUE`).
#'   Cost: O(N) likelihood evaluations per proposal, where N is the number
#'   of candidate reattachment edges.
#' @param gibbsSubtreeSwap Logical; include the Gibbs subtree-swap move
#'   (default `TRUE`). Cost: O(N) likelihood evaluations.
#' @param weightedBranchScale Logical; include the weighted branch-length
#'   scale move (default `FALSE`). Cost: O(B) likelihood evaluations, where
#'   B = `nBranchBins`.
#' @param weightedSpr Logical; include the weighted SPR move (default
#'   `FALSE`). Cost: O(N * B) likelihood evaluations. Enable for improved
#'   mixing on difficult tree spaces at the cost of slower iterations.
#' @param weightedSubtreeSwap Logical; include the weighted subtree-swap
#'   move (default `FALSE`). Cost: O(N * B) likelihood evaluations.
#' @param tbr Logical; include the TBR (Tree Bisection and Reconnection)
#'   topology move (default `TRUE`). TBR is a superset of SPR: it additionally
#'   re-roots the pruned subtree at a random internal edge before regrafting,
#'   enabling larger jumps in tree space. Cost: O(nEdge) per proposal (same
#'   as SPR). The adaptive scheduler will downweight TBR on small trees where
#'   SPR is sufficient.
#' @param pSpr Logical; include the parsimony-guided SPR move (default
#'   `TRUE`). Like standard SPR but weights candidate regraft positions by
#'   Fitch parsimony score, so topologically better positions are proposed
#'   more often. Much cheaper than Gibbs SPR (no likelihood evaluation per
#'   candidate), but better guided than uniform SPR.
#'   See Yang & Rodríguez (2013); Ronquist et al. (2020).
#' @param joint2d Logical; include 2D joint Bactrian proposals for
#'   correlated parameter pairs (default `TRUE`). Proposes correlated
#'   updates to tree_length × rate_log_sd (and tree_length × rate_loss
#'   when neomorphic characters are present) using a bivariate Bactrian
#'   kernel. The correlation is learned adaptively during warmup from
#'   posterior sample correlations. No effect when parameters are
#'   uncorrelated (degenerates to independent proposals).
#' @param blockGibbsBranch Logical; include the block Gibbs branch-length
#'   sweep move (default `FALSE`). Each call sweeps over all edge pairs
#'   in random-permutation order, sampling each from an approximate
#'   conditional via `nBranchBins` bin evaluations and independent MH
#'   accept/reject. Cost: O(nEdge * B) likelihood evaluations per call.
#'   Useful when single-edge moves (BetaSimplex, weighted branch scale)
#'   mix slowly through the branch-length space. The adaptive scheduler
#'   accounts for the multi-dimensional nature of this move via its `dim`
#'   field.
#' @param dirichletBranch Logical; include the block Dirichlet simplex
#'   branch-length move (default `TRUE`). Each proposal selects K random
#'   edges, draws new fractions from a Dirichlet centered on current values,
#'   and rescales the remaining edges to maintain the simplex. Cost: one
#'   likelihood evaluation per proposal (same as BetaSimplex).
#' @param dirichletK Integer or `NULL`; number of edges to update per
#'   Dirichlet branch proposal. Default `NULL` uses `min(nEdge, 5)`.
#'   Smaller values give higher acceptance but smaller moves; larger
#'   values give bolder moves at the cost of lower acceptance.
#' @param localDirichlet Logical; include a localized Dirichlet branch
#'   proposal that selects K *connected* edges via BFS from a random
#'   starting edge (default `TRUE`). The connected selection produces
#'   compact dirty sets for efficient partial CL evaluation, and targets
#'   correlated local branch lengths.
#' @param localDirichletK Integer or `NULL`; number of connected edges
#'   for the localized Dirichlet proposal. Default `NULL` uses
#'   `min(nEdge, 6)`.
#' @param nBranchBins Integer; number of branch-fraction bins for weighted
#'   and block Gibbs moves (default `10L`). Used by `weightedBranchScale`,
#'   `weightedSpr`, `weightedSubtreeSwap`, and `blockGibbsBranch`.
#'   Ignored when all bin-based moves are disabled. Higher values increase
#'   accuracy of the Gibbs approximation but cost more likelihood
#'   evaluations.
#' @param moveWeights Named numeric vector of user-pinned move weights,
#'   or `NULL` (default). When non-NULL, each named entry fixes the
#'   probability of proposing that move type. Names must match valid move
#'   names (e.g., `"nni"`, `"spr"`, `"gibbs_spr"`, `"tree_length"`, etc.).
#'   Values must be positive and sum to at most 1. Remaining probability
#'   is distributed among un-pinned moves by the adaptive scheduler during
#'   warmup. Example: `moveWeights = c(nni = 0.3, spr = 0.2)` fixes NNI
#'   at 30% and SPR at 20%, with the remaining 50% allocated adaptively
#'   among other moves. To disable adaptive scheduling entirely, pin all
#'   moves (sum to 1). See section **Adaptive move scheduling** below.
#' @param nCore Integer. Number of parallel worker processes for independent
#'   runs. Default `getOption("mc.cores", 1L)`, matching the convention used
#'   by \pkg{TreeDist}. With `nCore = 1` (the default), runs execute serially.
#'   With `nCore > 1` and `nRuns > 1`, runs are dispatched as background R
#'   processes via [callr::r_bg()] and the parent process polls for
#'   convergence. If `nRuns == 1`, `nCore` is ignored (within-run
#'   parallelism is a separate facility). When `nRuns > nCore`, runs are
#'   dispatched from a rolling pool of `nCore` workers: as each finishes,
#'   the next pending run launches in its place. Cross-run convergence
#'   (`maxRhat`) cannot trigger early in this regime — it activates only
#'   once every run has produced samples — so prefer `nRuns <= nCore`
#'   when convergence-based early stopping matters.
#' @param pollInterval Integer. Seconds between convergence polls in parallel
#'   mode. Ignored when `nCore = 1`. Default `10L`.
#' @param cacheBonus Numeric; multiplier applied to partial-CL-eligible
#'   move weights (NNI, beta_simplex, Dirichlet, local_dirichlet) when the
#'   node CL cache is valid. Default 5. A value of 1 disables the boost.
#'   Higher values preferentially select fast partial-CL moves after
#'   another partial-CL move succeeds, exploiting the ~15\eqn{\times}
#'   speedup. Statistically valid: each component kernel individually
#'   satisfies detailed balance; only selection frequency changes.
#'   Ignored when Q-heterogeneity is enabled (partial CL not supported).
#' @param tuning Named list of initial tuning parameters for each move
#'   type. See Details.
#'
#' @details
#' Default tuning values:
#' - `scale_tree_length`: 0.5
#' - `beta_simplex`: 10
#' - `scale_rate_loss`: 0.5
#' - `scale_rate_log_sd`: 0.5
#' - `scale_p`: 0.5
#' - `int_walk_window`: 1
#'
#' ## Parallel tempering
#'
#' When `nChains > 1`, chains run at different temperatures. The cold
#' chain (beta = 1) targets the true posterior. Heated chains (beta < 1)
#' flatten the likelihood surface, aiding exploration. Chain swap proposals
#' (see M-030) exchange states between adjacent temperatures.
#'
#' ## Gibbs and weighted moves
#'
#' The Gibbs moves (`gibbsSpr`, `gibbsSubtreeSwap`) evaluate all candidate
#' topologies and sample proportional to their posterior weight. They cost
#' O(N) likelihood evaluations per proposal (where N is the number of
#' candidates, roughly the number of edges) but often achieve much better
#' mixing than standard NNI/SPR. On larger trees (> ~20 tips), the
#' per-invocation cost grows quadratically and may outweigh the mixing gain.
#' Initial Gibbs weights are therefore capped to limit their time share;
#' the adaptive scheduler refines further during warmup.
#'
#' The weighted moves (`weightedBranchScale`, `weightedSpr`,
#' `weightedSubtreeSwap`) additionally marginalise over a discrete grid
#' of branch-fraction placements. This produces proposals that are
#' approximately independent of the current state but costs
#' O(B) or O(N * B) likelihood evaluations. These are off by default
#' and recommended only when standard + Gibbs moves show poor mixing.
#'
#' The block Gibbs branch-length sweep (`blockGibbsBranch`) updates all
#' edge pairs in a single MCMC move via random-permutation-scan
#' MH-within-Gibbs. Each edge pair is sampled from an approximate
#' conditional (bin-based, same as `weightedBranchScale`) and
#' accepted/rejected independently. Cost: O(nEdge * B). Off by default.
#' The adaptive scheduler accounts for its multi-dimensional nature via
#' a `dim` field (= nEdge) in the score formula.
#'
#' ## Three-phase MCMC: Warmup / Tuning / Sample
#'
#' The MCMC engine runs in three phases, shown in the progress display:
#'
#' **Warmup.** Tuning parameters (scale, window), move weights
#' (acceptance-rate/cost heuristic), and temperatures (parallel
#' tempering) all adapt. The phase ends automatically when the cold
#' chain's log-posterior stabilises (Geweke z-score test), or when
#' `maxWarmup` is reached. No samples are saved.
#'
#' **Tuning.** The chain is approximately stationary. Move weights
#' are optimised to maximise min-ESS/s (minimum effective sample
#' size per second across parameters) using a perturbation bandit.
#' Each round evaluates the current and perturbed weight vectors
#' over short windows, adopting the best. Tuning samples are
#' discarded. Step-size tuning is frozen. The phase runs for
#' `tuningRounds` rounds or up to `tuningBudget` iterations.
#' Set `autoTune = FALSE` to skip this phase.
#'
#' **Sample.** Move weights are frozen and posterior samples are
#' collected. Convergence is monitored at `checkEvery` intervals.
#' The run terminates when `minEss`/`maxRhat` criteria are met,
#' `maxTime` is reached, or `nIter` iterations complete.
#'
#' Use `moveWeights` to pin specific move frequencies and exclude
#' them from adaptation in all phases.
#'
#' @return An S3 object of class `MkPrimeMCMC`.
#' @export
MkPrimeMCMC <- function(
    nIter = Inf,
    thin = "auto",
    treeThin = NULL,
    warmup = NULL,
    minWarmup = 2000L,
    maxWarmup = NULL,
    autoTune = TRUE,
    tuningBudget = 10000L,
    tuningRounds = 5L,
    nRuns = 2L,
    nChains = 1L,
    heat = 0.2,
    maxTime = NULL,
    minEss = NULL,
    maxRhat = NULL,
    minTreeEss = NULL,
    checkEvery = 1000L,
    cancelFile = NULL,
    checkpointFile = NULL,
    treeFile = NULL,
    logFile = NULL,
    bufferSize = 500L,
    plotEvery = NULL,
    progressFn = NULL,
    gibbsSpr = TRUE,
    gibbsSubtreeSwap = TRUE,
    tbr = TRUE,
    pSpr = TRUE,
    joint2d = TRUE,
    weightedBranchScale = FALSE,
    weightedSpr = FALSE,
    weightedSubtreeSwap = FALSE,
    blockGibbsBranch = FALSE,
    dirichletBranch = TRUE,
    dirichletK = NULL,
    localDirichlet = TRUE,
    localDirichletK = NULL,
    nBranchBins = 10L,
    moveWeights = NULL,
    cacheBonus = 5,
    tuning = list(),
    nCore = getOption("mc.cores", 1L),
    pollInterval = 10L,
    gibbsWarmupFactor = 1/3
) {
  nIter <- if (is.infinite(nIter)) Inf else as.integer(nIter)
  if (!identical(thin, "auto")) {
    thin <- as.integer(thin)
    if (is.na(thin) || thin < 1L) {
      cli::cli_abort("{.arg thin} must be {.val auto} or a positive integer.")
    }
  }
  if (!is.null(treeThin)) {
    treeThin <- as.integer(treeThin)
    if (is.na(treeThin) || treeThin < 1L) {
      cli::cli_abort("{.arg treeThin} must be a positive integer or NULL.")
    }
  }
  # --- Three-phase warmup / tuning / sampling parameters ---
  # Legacy `warmup` maps to `maxWarmup` with deprecation notice.
  # When warmup is explicitly passed and autoTune was not, disable
  # autoTune for backward compatibility (old two-phase behaviour).
  autoTuneExplicit <- "autoTune" %in% names(match.call())

  if (!is.null(warmup)) {
    cli::cli_warn(c(
      "{.arg warmup} is deprecated; use {.arg maxWarmup} instead.",
      "i" = "Setting {.arg maxWarmup} = {warmup}."
    ))
    if (is.null(maxWarmup)) maxWarmup <- as.integer(warmup)
    # Legacy compat: disable autoTune unless user explicitly asked for it
    if (!autoTuneExplicit) autoTune <- FALSE
  }
  if (is.null(maxWarmup)) {
    maxWarmup <- if (is.finite(nIter)) as.integer(nIter / 2L) else 50000L
  }
  maxWarmup  <- as.integer(maxWarmup)
  minWarmup  <- as.integer(minWarmup)
  autoTune   <- as.logical(autoTune)
  tuningBudget <- as.integer(tuningBudget)
  tuningRounds <- as.integer(tuningRounds)

  if (minWarmup > maxWarmup) {
    cli::cli_warn(c(
      "{.arg minWarmup} ({minWarmup}) exceeds {.arg maxWarmup} ({maxWarmup}).",
      "i" = "Setting {.arg minWarmup} = {.arg maxWarmup}."
    ))
    minWarmup <- maxWarmup
  }
  if (autoTune && tuningBudget < 1L) {
    cli::cli_abort("{.arg tuningBudget} must be a positive integer.")
  }
  if (autoTune && tuningRounds < 1L) {
    cli::cli_abort("{.arg tuningRounds} must be a positive integer.")
  }

  # Backward compat: store maxWarmup in the `warmup` slot so that

  # existing code paths (.RunMkPrimeSingleRun) can use it as the
  # hard ceiling. The state machine will handle auto-detection.
  warmup <- maxWarmup

  nRuns <- as.integer(nRuns)
  nChains <- as.integer(nChains)

  if (nRuns < 1L) {
    cli::cli_abort("{.arg nRuns} must be at least 1.")
  }
  if (nChains < 1L) {
    cli::cli_abort("{.arg nChains} must be at least 1.")
  }
  if (nChains > 1L) {
    if (heat <= 0 || heat >= 1) {
      cli::cli_abort("{.arg heat} must be in (0, 1), got {heat}.")
    }
  }
  if (!is.null(cancelFile)) {
    if (!is.character(cancelFile) || length(cancelFile) != 1L) {
      cli::cli_abort("{.arg cancelFile} must be a length-1 character string or NULL.")
    }
  }
  if (!is.null(logFile)) {
    if (!is.character(logFile) || length(logFile) != 1L) {
      cli::cli_abort("{.arg logFile} must be a length-1 character string or NULL.")
    }
    # Append .log if no extension
    if (!grepl("\\.", basename(logFile))) {
      logFile <- paste0(logFile, ".log")
    }
    # Auto-derive treeFile from logFile when not specified
    if (is.null(treeFile)) {
      treeFile <- sub("\\.[^.]+$", "_trees.nwk", logFile)
    }
    # Auto-derive checkpointFile from logFile when not specified
    if (is.null(checkpointFile)) {
      checkpointFile <- sub("\\.[^.]+$", ".ckp", logFile)
    }
  }
  # FALSE → disable (convert to NULL for downstream code)
  if (identical(treeFile, FALSE)) treeFile <- NULL
  if (identical(checkpointFile, FALSE)) checkpointFile <- NULL
  bufferSize <- as.integer(bufferSize)
  if (bufferSize < 1L) {
    cli::cli_abort("{.arg bufferSize} must be a positive integer.")
  }

  # Validate move toggles
  tbr <- as.logical(tbr)
  pSpr <- as.logical(pSpr)
  joint2d <- as.logical(joint2d)
  gibbsSpr <- as.logical(gibbsSpr)
  gibbsSubtreeSwap <- as.logical(gibbsSubtreeSwap)
  weightedBranchScale <- as.logical(weightedBranchScale)
  weightedSpr <- as.logical(weightedSpr)
  weightedSubtreeSwap <- as.logical(weightedSubtreeSwap)
  blockGibbsBranch <- as.logical(blockGibbsBranch)
  dirichletBranch <- as.logical(dirichletBranch)
  if (!is.null(dirichletK)) {
    dirichletK <- as.integer(dirichletK)
    if (is.na(dirichletK) || dirichletK < 2L) {
      cli::cli_abort("{.arg dirichletK} must be at least 2, got {dirichletK}.")
    }
  }
  localDirichlet <- as.logical(localDirichlet)
  if (!is.null(localDirichletK)) {
    localDirichletK <- as.integer(localDirichletK)
    if (is.na(localDirichletK) || localDirichletK < 2L) {
      cli::cli_abort("{.arg localDirichletK} must be at least 2, got {localDirichletK}.")
    }
  }
  nBranchBins <- as.integer(nBranchBins)
  if (nBranchBins < 2L) {
    cli::cli_abort("{.arg nBranchBins} must be at least 2, got {nBranchBins}.")
  }

  # Validate moveWeights (M-092: adaptive scheduler)
  if (!is.null(moveWeights)) {
    if (is.list(moveWeights)) moveWeights <- unlist(moveWeights)
    if (!is.numeric(moveWeights) || is.null(names(moveWeights))) {
      cli::cli_abort(
        "{.arg moveWeights} must be a named numeric vector or NULL."
      )
    }
    validNames <- c(
      # Continuous parameter moves
      "tree_length", "branch_lengths", "rate_loss", "rate_log_sd", "rate_neo",
      "beta_scale", "neo_joint",
      # Topology moves
      "nni", "spr", "tbr", "pspr",
      "gibbs_spr", "gibbs_subtree_swap",
      "weighted_branch_lengths", "weighted_spr", "weighted_subtree_swap",
      "block_gibbs_branch", "dirichlet_branch", "local_dirichlet",
      # k' moves
      "kPrime", "p", "gibbs_kPrime", "block_kPrime",
      # BG hyperparameter moves
      "kprime_alpha", "kprime_beta",
      "slice_kprime_alpha", "slice_kprime_beta",
      # Slice samplers
      "slice_rate_loss", "slice_rate_neo", "slice_rate_log_sd",
      "slice_beta_scale",
      # Joint 2D moves
      "joint_tl_rls", "joint_tl_rl"
    )
    bad <- setdiff(names(moveWeights), validNames)
    if (length(bad) > 0L) {
      cli::cli_abort(
        "{.arg moveWeights} contains unknown move name{?s}: {.val {bad}}."
      )
    }
    if (any(moveWeights <= 0)) {
      cli::cli_abort("All {.arg moveWeights} values must be positive.")
    }
    if (sum(moveWeights) > 1.0 + 1e-8) {
      cli::cli_abort(
        "{.arg moveWeights} sum to {sum(moveWeights)}, which exceeds 1.0."
      )
    }
  }

  defaults <- list(
    scale_tree_length = 0.5,
    beta_simplex = 10,
    scale_rate_loss = 0.5,
    scale_rate_log_sd = 0.5,
    scale_p = 0.5,
    scale_logit_p = 1.0,
    scale_rate_neo = 0.5,
    scale_neo_joint = 0.5,
    scale_beta_scale = 0.5,
    scale_kprime_alpha = 0.3,
    scale_kprime_beta = 0.5,
    scale_joint_tl_rls = 0.5,
    scale_joint_tl_rl = 0.5,
    dirichlet_alpha = 0.1,
    local_dirichlet_alpha = 0.1,
    int_walk_window = 1L,
    slice_width_rate_loss = 1.0,
    slice_width_rate_neo = 1.0,
    slice_width_rate_log_sd = 1.0,
    slice_width_tree_length = 1.0,
    slice_width_beta_scale = 1.0,
    slice_width_kprime_alpha = 0.3,
    slice_width_kprime_beta = 0.5
  )
  tuning <- modifyList(defaults, tuning)

  if (!is.null(checkEvery)) checkEvery <- as.integer(checkEvery)
  if (!is.null(plotEvery)) plotEvery <- as.integer(plotEvery)

  nCore <- as.integer(nCore)
  if (is.na(nCore) || nCore < 1L) {
    cli::cli_abort("{.arg nCore} must be a positive integer.")
  }
  if (nCore > parallel::detectCores(logical = FALSE)) {
    cli::cli_warn(c(
      "{.arg nCore} = {nCore} exceeds physical cores ({parallel::detectCores(logical = FALSE)}).",
      "i" = "Proceeding anyway; reduce if memory-bound."
    ))
  }
  pollInterval <- as.integer(pollInterval)
  if (pollInterval < 1L) {
    cli::cli_abort("{.arg pollInterval} must be a positive integer.")
  }
  cacheBonus <- as.numeric(cacheBonus)
  if (is.na(cacheBonus) || cacheBonus < 1) {
    cli::cli_abort("{.arg cacheBonus} must be >= 1, got {cacheBonus}.")
  }

  # Resolve progressFn / plotEvery defaults
  # Use match.call() to distinguish "user passed NULL" from "user omitted arg"
  mc <- match.call()
  plotEveryExplicit   <- "plotEvery"  %in% names(mc)
  progressFnExplicit <- "progressFn" %in% names(mc)

  # In interactive sessions, enable live trace plots unless explicitly disabled
  if (interactive() && !plotEveryExplicit && !progressFnExplicit) {
    plotEvery <- checkEvery
    progressFn <- MkpTracePlot
  }

  if (identical(progressFn, "default")) {
    progressFn <- MkpTracePlot
  }

  # If plotEvery is set but progressFn was not explicitly provided,
  # default to MkpTracePlot
  if (!is.null(plotEvery) && !progressFnExplicit && is.null(progressFn)) {
    progressFn <- MkpTracePlot
  }

  if (!is.null(progressFn) && !is.function(progressFn)) {
    cli::cli_abort(
      "{.arg progressFn} must be a function, {.val default}, or {.val NULL}."
    )
  }

  structure(
    list(nIter = nIter, thin = thin, treeThin = treeThin, warmup = warmup,
         minWarmup = minWarmup, maxWarmup = maxWarmup,
         autoTune = autoTune, tuningBudget = tuningBudget,
         tuningRounds = tuningRounds,
         nRuns = nRuns, nChains = nChains, heat = heat,
         maxTime = maxTime, minEss = minEss, maxRhat = maxRhat,
         minTreeEss = minTreeEss, checkEvery = checkEvery,
         cancelFile = cancelFile, checkpointFile = checkpointFile,
         treeFile = treeFile, logFile = logFile, bufferSize = bufferSize,
         plotEvery = plotEvery, progressFn = progressFn,
         tbr = tbr, pSpr = pSpr, joint2d = joint2d,
         gibbsSpr = gibbsSpr, gibbsSubtreeSwap = gibbsSubtreeSwap,
         weightedBranchScale = weightedBranchScale,
         weightedSpr = weightedSpr,
         weightedSubtreeSwap = weightedSubtreeSwap,
         blockGibbsBranch = blockGibbsBranch,
         dirichletBranch = dirichletBranch,
         dirichletK = dirichletK,
         localDirichlet = localDirichlet,
         localDirichletK = localDirichletK,
         nBranchBins = nBranchBins,
         moveWeights = moveWeights,
         cacheBonus = cacheBonus,
         tuning = tuning,
         nCore = nCore, pollInterval = pollInterval,
         gibbsWarmupFactor = gibbsWarmupFactor),
    class = "MkPrimeMCMC"
  )
}
