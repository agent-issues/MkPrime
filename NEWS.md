# MkPrime (development version)

* `RunMkPrime(tree = NULL)` now starts from a greedy parsimony
  stepwise-addition tree (`TreeSearch::AdditionTree`) when `TreeSearch`
  is installed, falling back to `TreeTools::NJTree` otherwise.  The
  parsimony tree is a methodologically preferred MCMC start for
  morphological data, and `AdditionTree` is essentially instant
  (no iterative search).  Behaviour is unchanged when a `tree`
  argument is supplied.

* Streaming log writer (`.FlushBuffer()`) now appends one row per
  `writeLines()` call on an explicitly opened append connection.  The
  previous implementation built the whole flush block as a single
  concatenated string and passed it to `cat()`; on networked filesystems
  the underlying `write()` could be split across syscalls, occasionally
  leaving the log file with a torn row that broke `ReadMkLog()` and
  `RelabelEcology()` on resume.  Per-row writes keep every payload well
  below `PIPE_BUF`, so each line is a single atomic OS write.

* `.BuildResult()` now aligns `z_samples` with the run's authoritative
  `saved_idx` before returning, warning when they differ.  Surplus
  snapshots (typically carried over from a resume that truncated log
  rows) are trimmed instead of corrupting the `MkPosterior` invariant.

* `RelabelEcology()` aborts with a rich diagnostic message when
  `result$z_samples` length does not match `nrow(samples)`, replacing
  the previous bare length-mismatch error.  A new `trimZSamples =
  "tail"` / `"head"` argument allows deliberate trimming when the
  surplus direction is known (most common: a streaming run interrupted
  mid-flush left extra in-memory z snapshots not recorded in the log).
  A shortfall of z entries still aborts unconditionally -- relabelling
  needs a paired z snapshot per sample.

* Multiple MCMC runs now execute in parallel when `nCore > 1`
  (defaults to `getOption("mc.cores", 1L)`). Backend: `callr::r_bg()`.
  Mirrors the TreeDist option-driven pattern.
* **Breaking (pre-release):** `parallel` argument to `MkPrimeMCMC()`
  removed; use `nCore` instead. The `future` package is no longer used.

* Tree-length prior now scales with the data.  When `expSteps` is left at its
  default (`NULL`), `RunMkPrime()` / `.FinalizeModel()` set
  `expSteps = expStepsInflation * parsimony(startingTree, data) / nChar`
  and derive `treeLengthRate = 2 / expSteps`.  Dividing by `nChar` converts
  the total parsimony score (sum of changes across all characters) into the
  same per-character expected-substitutions units used by `tree_length`
  (sum of edge lengths).  The `expStepsInflation` argument (default 1.05)
  encodes the expectation that the true tree length sits a few percent
  above the parsimony minimum.  Numeric `expSteps` supplied by the user is
  still respected verbatim (no division or inflation).  **This changes the
  prior for every run that did not previously pass an explicit
  `expSteps`**; the resolved value is printed at chain start and shown by
  `print(model)`.

* Ecology-aware substitution model.  `MkPrimeModel(ecologyAware = TRUE)`
  enables per-character, per-ecology rate modification via a
  spike-and-slab latent `z` and global magnitude `phi`.  Tip ecology
  is passed via `MkPrimeData(..., ecology = <integer vector>)`.
  Designed to reduce the influence of convergent character-state
  similarity in ecologically similar but phylogenetically unrelated
  lineages.

  The v2 (current) parameterisation uses a **reference ecology**
  (auto-selected as the most frequent), **gamma-normalised** rate
  factors so the prior-expected per-cell rate is 1 (analogous to ACRV
  mean-rate normalisation), and an **asymmetric slab prior**
  controlled by per-ecology `theta_e ~ Beta(thetaAlpha, thetaBeta)`
  (defaults `(2, 2)`).  Together these close the rate-time
  identifiability ridge, eliminate the K-vs-K-1 redundancy in z, and
  remove the assumption that the prior expects equal numbers of
  encouraged vs depressed characters.  See `vignette("ecology-details")`
  for the math and `vignette("rodent-ecology")` for an empirical
  case study identifying per-character ecology associations.
* New `kPrimePrior = "empirical_geometric"` prior (now the default).  Decomposes
  `k' = N_obs + N_unobs`, with `N_obs` drawn from an empirical pmf tabulated
  from real morphological matrices (`empiricalNObs`) and `N_unobs` from a
  `Geometric(p)` with `Beta(a, b)` hyperprior on `p`.  The convolution of the
  two distributions defines the prior on `k'` and counters the previous
  tendency of inference to collapse the number of unobserved states to zero.
* New `MkPrimeEmpiricalPrior()` constructor and `empiricalNObs` package
  dataset, allowing users to override the empirical component.
* `gibbs_p` Gibbs draw on `p` is disabled under `"empirical_geometric"`
  (no Beta conjugacy under the convolution); replaced with a logit-scale
  Metropolis-Hastings move `mh_logit_p` (move type 30).  The multiplicative
  `mh_p` move (type 8) is rejected for most proposals when the posterior on
  `p` concentrates near 1 (typical under the empirical_geometric prior),
  causing the adaptive scheduler to crush its weight to the floor and `p`
  to stay stuck; proposing on the unbounded logit scale fixes this.
