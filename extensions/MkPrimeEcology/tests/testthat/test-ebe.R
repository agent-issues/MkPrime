# The EBE extension against contracts in ebe-spec.md section 8, using an
# oracle written from the specification in plain R (not from src/).

# --- Independent R oracle -------------------------------------------------

.JCP <- function(t, K) {
  e <- exp(-K * t / (K - 1))
  matrix((1 - e) / K, K, K) + diag(e, K)
}

# Posterior ecology marginals at every node by enumerating internal states;
# JC(K), uniform root, raw edge lengths. Rows: nodes; columns: ecologies.
.OracleMarginals <- function(tree, eco, K) {
  edge <- tree$edge
  nTip <- length(eco)
  nNode <- max(edge)
  internal <- (nTip + 1):nNode
  root <- setdiff(edge[, 1], edge[, 2])
  Ps <- lapply(tree$edge.length, .JCP, K = K)
  configs <- as.matrix(expand.grid(rep(list(seq_len(K)), length(internal))))
  acc <- matrix(0, nNode, K)
  for (i in seq_len(nrow(configs))) {
    s <- integer(nNode)
    s[internal] <- configs[i, ]
    w <- 1 / K
    tipVec <- vector("list", nTip)
    for (e in seq_len(nrow(edge))) {
      p <- edge[e, 1]
      ch <- edge[e, 2]
      if (ch > nTip) {
        w <- w * Ps[[e]][s[p], s[ch]]
      } else {
        ind <- if (is.na(eco[ch])) rep(1, K) else as.numeric(seq_len(K) == eco[ch])
        v <- Ps[[e]][s[p], ] * ind
        tipVec[[ch]] <- v / sum(v)
        w <- w * sum(v)
      }
    }
    for (n in internal) acc[n, s[n]] <- acc[n, s[n]] + w
    for (tip in seq_len(nTip)) acc[tip, ] <- acc[tip, ] + w * tipVec[[tip]]
  }
  # Return:
  acc / rowSums(acc)
}

# Log-likelihood of neomorphic characters `x` (nTip x nChar, 0/1/NA) with
# influence z (nChar x nNonRef), by pruning with explicit mixed matrices.
.OracleLogLik <- function(tree, x, z, phi, rateLoss, eco, K, ref,
                          perEcology = FALSE, neoScale = 1, rates = 1,
                          coding = 0L) {
  marg <- .OracleMarginals(tree, eco, K)
  edge <- tree$edge
  wEdge <- (marg[edge[, 1], , drop = FALSE] + marg[edge[, 2], , drop = FALSE]) / 2
  nonRef <- setdiff(seq_len(K), ref)
  pi1b <- 1 / (1 + rateLoss)
  root <- setdiff(edge[, 1], edge[, 2])
  nTip <- nrow(x)
  post <- rev(seq_len(nrow(edge)))  # tree is in preorder
  CharLik <- function(col, zc) {
    Pi1 <- function(s) {
      j <- match(s, nonRef)
      if (is.na(j) || zc[j] == 0) return(pi1b)
      ph <- if (perEcology) phi[j] else phi[1]
      plogis(qlogis(pi1b) + (if (zc[j] == 1) 1 else -1) * log(ph))
    }
    mean(vapply(rates, function(rate) {
      cl <- matrix(1, max(edge), 2)
      for (tip in seq_len(nTip)) {
        if (!is.na(col[tip])) cl[tip, ] <- as.numeric(0:1 == col[tip])
      }
      for (e in post) {
        t <- tree$edge.length[e] * neoScale * rate
        ex <- exp(-2 * t)
        Pmix <- matrix(0, 2, 2)
        for (s in seq_len(K)) {
          p1 <- Pi1(s)
          p0 <- 1 - p1
          Ps <- rbind(c(p0 + p1 * ex, p1 * (1 - ex)),
                      c(p0 * (1 - ex), p1 + p0 * ex))
          Pmix <- Pmix + wEdge[e, s] * Ps
        }
        cl[edge[e, 1], ] <- cl[edge[e, 1], ] * as.vector(Pmix %*% cl[edge[e, 2], ])
      }
      sum(c(1 - pi1b, pi1b) * cl[root, ])
    }, numeric(1)))
  }
  ll <- 0
  for (ch in seq_len(ncol(x))) {
    zc <- z[ch, ]
    li <- log(CharLik(x[, ch], zc))
    if (coding == 1L) {
      obs <- !is.na(x[, ch])
      c0 <- ifelse(obs, 0, NA)
      c1 <- ifelse(obs, 1, NA)
      li <- li - log1p(-(CharLik(c0, zc) + CharLik(c1, zc)))
    }
    ll <- ll + li
  }
  # Return:
  ll
}

# --- Fixtures --------------------------------------------------------------

.EbeFixture <- function(nTip = 5L, nChar = 4L, K = 3L, seed = 8101L) {
  set.seed(seed)
  tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = paste0("t", seq_len(nTip)),
                                         br = function(n) runif(n, 0.05, 0.6)))
  x <- matrix(sample(0:1, nTip * nChar, replace = TRUE), nTip)
  x[, 1] <- c(0L, 1L, rep(NA, nTip - 2L))[seq_len(nTip)]
  if (nChar >= 2L) x[3, 2] <- NA
  eco <- sample.int(K, nTip, replace = TRUE)
  eco[1:K] <- seq_len(K)
  eco[nTip] <- NA
  list(tree = tree, x = x, eco = eco, K = K)
}

# An extension built directly from integer ecology (1-based) and nNeo.
.EbeExt <- function(fx, perEcology = FALSE, ref = 1L) {
  K <- fx$K
  nPhi <- if (perEcology) K - 1L else 1L
  MkPrimeEcology:::.EcologyExtensionPtr(
    as.integer(fx$eco) - 1L, K, ref - 1L, ncol(fx$x), perEcology,
    paste0("p", seq_len(nPhi + K)), rep(0.5, nPhi + K),
    paste0("z", seq_len(ncol(fx$x) * (K - 1L))), c(3, 2, 2, 2, 1.5), 1
  )
}

.EbeLL <- function(ptr, fx, tree, z, params, rateLoss, neoScale = 1,
                   rates = 1, coding = 0L) {
  x <- fx$x
  x[is.na(x)] <- -1L
  storage.mode(x) <- "integer"
  MkPrimeEcology:::.EcologyLogLik(
    ptr, tree$edge[, 1], tree$edge[, 2], tree$edge.length, x, rateLoss,
    neoScale, rates, coding, params, as.integer(t(z))
  )
}

# --- T2: correctness against the oracle -----------------------------------

test_that("EBE likelihood matches the independent oracle (T2)", {
  fx <- .EbeFixture()
  K <- fx$K
  z <- matrix(c(1L, 2L, 0L, 1L, 2L, 2L, 1L, 0L), ncol(fx$x), K - 1L)
  for (perEco in c(FALSE, TRUE)) {
    ptr <- .EbeExt(fx, perEcology = perEco, ref = 2L)
    phi <- if (perEco) c(2.5, 0.4) else 3
    params <- c(phi, 0.7, 0.3, 0.6)
    for (coding in 0:1) {
      for (rates in list(1, c(0.4, 1.6))) {
        got <- .EbeLL(ptr, fx, fx$tree, z, params, rateLoss = 0.8,
                      neoScale = 1.3, rates = rates, coding = coding)
        want <- .OracleLogLik(fx$tree, fx$x, z, phi, 0.8, fx$eco, K, ref = 2L,
                              perEcology = perEco, neoScale = 1.3,
                              rates = rates, coding = coding)
        expect_equal(got, want, tolerance = 1e-10,
                     label = paste(perEco, coding, length(rates)))
      }
    }
  }
})

# --- T3: non-stationarity ---------------------------------------------------

test_that("moving the root changes the likelihood only when z != 0 (T3)", {
  fx <- .EbeFixture(nTip = 6L, seed = 8102L)
  ptr <- .EbeExt(fx)
  tr <- fx$tree
  rootEdges <- which(tr$edge[, 1] == TreeTools::RootNode(tr))
  moved <- tr
  total <- sum(tr$edge.length[rootEdges])
  moved$edge.length[rootEdges] <- total * c(0.15, 0.85)
  params <- c(3, 0.7, 0.3, 0.6)
  z0 <- matrix(0L, ncol(fx$x), fx$K - 1L)
  z1 <- z0
  z1[, 1] <- 1L
  expect_equal(.EbeLL(ptr, fx, tr, z0, params, 0.8),
               .EbeLL(ptr, fx, moved, z0, params, 0.8), tolerance = 1e-10)
  expect_gt(abs(.EbeLL(ptr, fx, tr, z1, params, 0.8) -
                  .EbeLL(ptr, fx, moved, z1, params, 0.8)), 1e-6)
})

# --- T4: normalisation -------------------------------------------------------

test_that("pattern probabilities sum to one (T4)", {
  fx <- .EbeFixture(nTip = 5L, nChar = 1L, seed = 8103L)
  ptr <- .EbeExt(fx)
  z <- matrix(c(1L, 2L), 1L)
  patterns <- as.matrix(expand.grid(rep(list(0:1), 5L)))
  total <- sum(apply(patterns, 1, function(pat) {
    f <- fx
    f$x <- matrix(as.integer(pat), ncol = 1L)
    exp(.EbeLL(ptr, f, fx$tree, z, c(3, 0.7, 0.3, 0.6), 0.8,
               rates = c(0.5, 1.5)))
  }))
  expect_equal(total, 1, tolerance = 1e-12)
})

# --- Prior and the Gibbs move -------------------------------------------------

test_that("the prior is the documented one", {
  fx <- .EbeFixture()
  ptr <- .EbeExt(fx)
  z <- c(0L, 1L, 2L, 0L, 1L, 0L, 0L, 2L)
  params <- c(2, 0.7, 0.3, 0.6)
  pz <- function(z, th) ifelse(z == 0, 0.7, 0.3 * ifelse(z == 1, th, 1 - th))
  theta <- rep(c(0.3, 0.6), 4)
  want <- dlnorm(2, 0, 1.5, log = TRUE) + dbeta(0.7, 3, 2, log = TRUE) +
    dbeta(0.3, 2, 2, log = TRUE) + dbeta(0.6, 2, 2, log = TRUE) +
    sum(log(pz(z, theta)))
  expect_equal(MkPrimeEcology:::.EcologyLogPrior(ptr, params, z), want,
               tolerance = 1e-12)
  expect_identical(MkPrimeEcology:::.EcologyLogPrior(ptr, c(2, 1.2, 0.3, 0.6),
                                                     z), -Inf)
})

test_that("gibbs_z draws from the exact full conditional", {
  fx <- .EbeFixture(seed = 8104L)
  ptr <- .EbeExt(fx)
  params <- c(4, 0.6, 0.3, 0.7)
  x <- fx$x
  x[is.na(x)] <- -1L
  storage.mode(x) <- "integer"
  Target <- function(zv) {
    MkPrimeEcology:::.EcologyLogPrior(ptr, params, zv) +
      MkPrimeEcology:::.EcologyLogLik(ptr, fx$tree$edge[, 1],
                                      fx$tree$edge[, 2], fx$tree$edge.length,
                                      x, 0.8, 1, 1, 1L, params, zv)
  }
  set.seed(8105L)
  z <- integer(ncol(x) * (fx$K - 1L))
  changed <- 0L
  for (i in 1:60) {
    prop <- MkPrimeEcology:::.EcologyProposeZ(
      ptr, fx$tree$edge[, 1], fx$tree$edge[, 2], fx$tree$edge.length, x, 0.8,
      1, 1, 1L, params, z)
    expect_true(all(prop$z %in% 0:2))
    expect_lte(sum(prop$z != z), fx$K - 1L)
    # Exact Gibbs: the Hastings ratio cancels the target ratio.
    expect_equal(prop$log_hastings + Target(prop$z) - Target(z), 0,
                 tolerance = 1e-9)
    changed <- changed + any(prop$z != z)
    z <- prop$z
  }
  expect_gt(changed, 0L)
})

# --- T1 and integration through MkPrime --------------------------------------

.EbeRunFixture <- function(seed = 8106L) {
  set.seed(seed)
  nTip <- 9L
  mat <- matrix(sample.int(3, nTip * 16, replace = TRUE) - 1L, nTip, 16,
                dimnames = list(paste0("t", seq_len(nTip)), NULL))
  mat[, 1:8] <- pmin(mat[, 1:8], 1L)
  mat <- matrix(as.character(mat), nrow(mat), dimnames = dimnames(mat))
  mat[2, 3] <- "?"
  pd <- TreeTools::MatrixToPhyDat(mat)
  tree <- TreeTools::Preorder(ape::rtree(nTip, tip.label = rownames(mat),
                                         br = function(n) runif(n, 0.05, 0.3)))
  ecology <- stats::setNames(rep(c("land", "water", "air"), length.out = nTip),
                             rownames(mat))
  ecology["t5"] <- NA
  list(pd = pd, tree = tree, ecology = ecology,
       mkd = MkPrime::MkPrimeData(pd, neomorphic = 1:8))
}

.EbePointers <- function(fx, model) {
  model <- MkPrime:::.FinalizeModel(model, fx$tree, fx$mkd)
  state <- MkPrime:::.InitState(fx$tree, fx$mkd, model)
  dataPtr <- MkPrime:::.InitMcmcData(fx$mkd, model)
  statePtr <- MkPrime:::.InitMcmcChain(state)
  MkPrime:::fill_partition_cache(dataPtr, statePtr)
  MkPrime:::allocate_cl_workspace(dataPtr, statePtr)
  list(dataPtr = dataPtr, statePtr = statePtr)
}

test_that("with every z = 0 the model is MkPrime's MkN (T1)", {
  fx <- .EbeRunFixture()
  for (coding in c("none", "variable")) {
    ext <- EcologyExtension(fx$mkd, fx$ecology, refEcology = "land")
    plain <- .EbePointers(fx, MkPrime::MkPrimeModel(coding = coding))
    eco <- .EbePointers(fx, MkPrime::MkPrimeModel(coding = coding,
                                                  extensions = ext))
    s <- MkPrime:::get_mcmc_state(eco$statePtr)
    expect_true(all(s$extInts == 0L))
    expect_equal(s$logLik + s$extLogLik,
                 MkPrime:::get_state_log_lik(plain$statePtr),
                 tolerance = 1e-9, label = coding)
    set.seed(8107L)
    for (i in 1:10) {
      for (type in c(0L, 1L, 2L, 3L, 5L)) {
        MkPrime:::do_move_cpp(eco$dataPtr, eco$statePtr, type, 0L, 0.6, 1.0,
                              1L, 1.0)
        s <- MkPrime:::get_mcmc_state(eco$statePtr)
        expect_equal(s$logLik + s$extLogLik,
                     MkPrime:::eval_full_loglik_cpp(plain$dataPtr,
                                                    eco$statePtr),
                     tolerance = 1e-9, label = paste(coding, type))
      }
    }
  }
})

test_that("RunMkPrime samples the ecology model", {
  skip_on_cran()
  fx <- .EbeRunFixture(8108L)
  ext <- EcologyExtension(fx$pd, fx$ecology, neomorphic = 1:8)
  model <- MkPrime::MkPrimeModel(extensions = list(ext))
  set.seed(8109L)
  res <- suppressWarnings(MkPrime::RunMkPrime(
    fx$pd, fx$tree, neomorphic = 1:8, model = model,
    mcmc = MkPrime::MkPrimeMCMC(nRuns = 1L, nIter = 3000L, thin = 5L,
                                minWarmup = 500L, maxWarmup = 500L,
                                autoTune = FALSE, maxTime = 120)
  ))
  s <- res$samples
  zCols <- grep("^eco_z_", colnames(s), value = TRUE)
  expect_length(zCols, 8L * 2L)
  expect_true(all(c("eco_phi", "eco_pi0", "eco_theta_water", "eco_theta_air")
                  %in% colnames(s)))
  expect_true(all(as.matrix(s[, zCols]) %in% 0:2))
  expect_true(all(is.finite(s[, "log_likelihood"])))
  expect_gt(res$acceptance[["ecology_gibbs_z"]], 0.5)
  folded <- CanonicaliseEcology(s, ext)
  expect_true(all(folded[, "eco_phi"] >= 1))
})

test_that("CanonicaliseEcology() folds phi < 1", {
  ext <- list(ecology = list(label = "b", magnitude_mode = "global",
                             characters = 1:2))
  s <- cbind(eco_phi = c(0.5, 2), eco_theta_b = c(0.2, 0.2),
             eco_z_1_b = c(1, 1), eco_z_2_b = c(2, 0))
  out <- CanonicaliseEcology(s, ext)
  expect_equal(out[, "eco_phi"], c(2, 2))
  expect_equal(out[, "eco_theta_b"], c(0.8, 0.2))
  expect_equal(out[, "eco_z_1_b"], c(2, 1))
  expect_equal(out[, "eco_z_2_b"], c(1, 0))
})
