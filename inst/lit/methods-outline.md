# Methods — section outline
# MkPrime ecology-aware paper
# Draft 2026-05-19
#
# Status: outline with pre-drafted key sentences; parameters pulled from code.
# [DECISION] markers flag items requiring a choice before prose drafting.
# Code references are absolute paths from repo root.

---

## 1. The ecology-aware MkPrime model

### 1.1 Baseline model

- Both blind and ecology-aware analyses use the Mk' model (Smith 2020a): the
  true number of character states k' is treated as an unknown, with k' ≥ k_obs,
  inferred jointly with the phylogeny under a hierarchical prior.
- Transformational characters (multistate equal-rates CTMCs) use the Mkv
  ascertainment correction (Lewis 2001) for variable-only coding.
- Neomorphic (presence–absence) characters are modelled as asymmetric two-state
  CTMCs following Pagel 1994, with a rate_loss parameter (gain/loss asymmetry).
- Among-character rate variation (ACRV): discretised LogNormal, nCat=6
  categories; `rate_log_sd ~ Gamma(1, 1)`.
- **One sentence only on Mk'**: both blind and aware models share
  `kPrimePrior = "geometric"`, `coding = "variable"`, `expSteps = 10`;
  the prior on k' is not the subject of this paper.

### 1.2 Ecology-aware extension

- For each pair (character c, ecology state e), a latent influence category
  z_{c,e} takes values `none` (0), `encouraged` (1), or `discouraged` (2).
- Spike-and-slab prior:
  - P(z_{c,e} = none)        = π₀
  - P(z_{c,e} = encouraged)  = (1 − π₀) · θ_e
  - P(z_{c,e} = discouraged) = (1 − π₀) · (1 − θ_e)
- Rate modifier: when z_{c,e} = encouraged, the per-edge substitution rate for
  character c in ecology e is multiplied by φ; when z_{c,e} = discouraged, by
  1/φ. The reference ecology (state 0) carries no z column and no modification.
- Single global φ shared across all characters and ecologies
  (`magnitudeMode = "global"`).
- γ normalisation: the per-edge rate is normalised so that the weighted mean
  modification factor across all ecology states on each edge equals 1, ensuring
  that the tree-length parameter retains its standard Mk interpretation as the
  expected number of state changes per unit branch length.
  [MATH: placeholder — one-line expression for γ; see vignette("ecology-details")]
- Edge ecology marginals: each edge's ecology-state mixture weights w_e(s) are
  taken from the marginal of a standard Mk(K) process at the parent node,
  reconstructed on-the-fly from the current tree and tip ecology vector.

### 1.3 Prior choices (simulation study — multirep-v3)

Hyperpriors in force when the eight-replicate simulation study was run
(committed and installed on Hamilton before 2026-05-14; pre-R5 defaults):

| Parameter | Prior | Note |
|---|---|---|
| π₀ | Beta(75, 25) — mode 0.75, ESS 100 | Pre-R5 default |
| θ_e | Beta(1, 1) uniform | Pre-R5 default |
| φ | LogNormal(0, σ_φ = 1.0) | Pre-R5 default |
| tree_length | Gamma(2, 2/expSteps) | expSteps=10 |
| rate_loss | LogNormal(0, 2) | |
| k' (geometric) | Geometric(p), p ~ Beta(1, 1) | |

These defaults were subsequently revised (commit `2fd4dd5`, 2026-05-16) as part
of red-team round R5 to counteract pi0/phi/theta feedback under low per-cell
likelihood contrast; the updated defaults are Beta(360, 120), Beta(2, 2), and
σ_φ = 1.5. The simulation results presented in this paper used the pre-R5
values.

[DECISION] Whether to re-run the 8 reps under the R5 defaults before
submission (task #18 PT array; results pending). If PT-extended reps are
available, that run will use the current (R5) defaults and should be described
separately or as an extension of the main result.

### 1.4 Prior choices (empirical analysis — rodent)

Hyperpriors used in the rodent analysis (`run_rodent_aware_v2.R`, committed
2026-05-17, post-R5):

| Parameter | Prior | Note |
|---|---|---|
| π₀ | Beta(7, 3) — mode 0.7 | Explicitly set |
| θ_e | Beta(2, 2) | R5 default |
| φ | LogNormal(0, σ_φ = 0.5) | Explicitly set; tighter than R5 default |
| tree_length | Gamma(2, 2/10) | expSteps=10 |
| rate_loss | LogNormal(0, 2) | |
| k' (geometric) | Geometric(p), p ~ Beta(1, 1) | |

[DECISION] Whether the same hyperprior set should be used for the final
simulation runs; currently the two analyses use different values. A brief
sensitivity section or a unified prior table would simplify the methods.

### 1.5 Identifiability and post-hoc relabelling

- The ecology-aware likelihood is invariant under the joint transformation
  (φ → 1/φ, θ_e → 1 − θ_e, swap z codes 1 ↔ 2 for all cells), creating a
  discrete identifiability degeneracy.
- Addressed post hoc via `RelabelEcology()`: for each posterior sample, if
  φ < 1 the full reflection is applied so that the reported posterior is in
  the φ ≥ 1 canonical form.
- The relabelling step is applied before any posterior summary or scoring.
  "We follow convention by reporting the φ ≥ 1 representative of each
  posterior sample; this does not affect topology or tree-length inference."

---

## 2. Simulation design (multirep-v3)

### 2.1 Tree topology and branch lengths

- 16-tip tree, four four-tip clades A, B, C, D in balanced clade-quartet
  topology: `((A1..A4, A2..A4), (C1..C4))` rooting: true topology
  `((A,C), (B,D))`.
- Newick construction via `.BuildConvergentTree()` in
  `inst/simulations/ecology/sim3-helpers.R`.
- Branch lengths (v3 config, set in `inst/hamilton/sim3-multirep-v3/run_rep.R`):
  - tipBr   = 0.50 (within-clade terminal branches)
  - stemBr  = 0.30 (clade stem branches)
  - rootBr  = 0.15 (internal branch separating the two cherry-of-clades)
  - Truth tree length = 13.50 (sum of all branch lengths)

### 2.2 Ecology assignment

- Ecology 1: tips A1–A4 and B1–B4 (the two phylogenetically distant clades
  sharing ecology — the convergent group).
- Ecology 0: tips C1–C4 and D1–D4 (reference ecology).
- The true phylogeny groups A with C and B with D; shared ecology creates
  spurious synapomorphies linking A with B.

### 2.3 Character simulation

- 480 characters total: 120 neomorphic + 360 transformational.
- The z-assignment matrix (nChar × 2) is set so that all 120 neomorphic
  characters have z = encouraged for ecology 1 (θ = 1, φ = 4); all 360
  transformational characters are ecology-neutral (z = none).
- Simulation function: `.SimulateMkPrimeEcology()` in
  `inst/simulations/ecology/sim3-simulate.R` (normalize = TRUE, baseRate = 1.0,
  pi0 = 0.75, refEcology = 0L, rateLoss = 1).
- "The neomorphic characters evolve under a 4× elevated gain rate in
  ecology-1 lineages, systematically generating parallel presence signals
  across the phylogenetically distant clades A and B."
- Eight independent replicates: RNG seed = 20260601 + repId (repId 1–8).

---

## 3. MCMC settings (simulation)

- 100 000 iterations per chain, single chain per model arm per replicate
  (`nChains = 1L, nRuns = 1L`; no parallel tempering in the main analysis).
- Start tree: parsimony optimal tree from `TreeSearch::MaximizeParsimony()`
  on the simulated character matrix; uniform initial branch lengths of 0.1.
  Start-tree seed fixed at 1 (shared across arms).
- Warmup: adaptive, minWarmup = 5 000, maxWarmup = 40 000 iterations.
- Thinning: `thin = treeThin = max(1L, nIter %/% 500L) = 200` — approximately
  500 scalar samples and 500 tree samples retained per chain.
- Streaming logs written to disk; no in-memory sample accumulation beyond the
  thinned set (RAM budget constraint on Hamilton /nobackup).
- Burn-in: first 25 % of retained (post-warmup) samples discarded before
  scoring, matching the inline `discard()` function in `run_rep.R`.

[DECISION] Whether to report results from the main single-chain runs only, or
to include the PT-extended runs (task #14/#18). If PT is included, describe
separately: nChains = 4, heat = 0.2, geometric temperature ladder; same nIter
and thinning.

---

## 4. Robustness simulations (v4 family)

### 4.1 Shared structure

- Same 16-tip balanced clade-quartet topology; character simulation via
  `.SimulateMkPrimeEcology()`.
- Parameters: nNeo = 100, nTrans = 200, φ = 4, π₀ = 0.75, θ = 1.0
  (for v4a/b/c); see `inst/simulations/ecology/sim3v4-params.R`.
- Model and MCMC specifications match the main simulation (geometric k' prior,
  100k iterations, nChains = 1 for non-PT arms).

### 4.2 v4 geometry vs v4-cross geometry

- **v4 geometry**: ecology-1 assigned to `{A1..A4, C1..C4}` — a within-
  clade confound that breaks monophyly of A and C without directly contradicting
  the (A,C)|(B,D) sister grouping.
- **v4-cross geometry**: ecology-1 assigned to `{A1..A4, B1..B4}` — the same
  convergent assignment as multirep-v3; the spurious eco signal directly
  competes with the true (A,C) bipartition.
  [Helper: `inst/simulations/ecology/sim3v4cross-helpers.R`]

### 4.3 Branch-length regimes

| Variant | stemBrEco | stemBrClade | rootBr | Purpose |
|---|---|---|---|---|
| v4a | 0.05 | 0.30 | 0.15 | Conservative: both chains should recover truth |
| v4b | 0.08 | 0.15 | 0.08 | Balanced: moderate eco/ancestry competition |
| v4c | 0.15 | 0.05 | 0.04 | Stress: near-saturation eco signal |
| v4-cross-b (pt5) | 0.08 | 0.15 | 0.08 | v4b params, cross geometry, 5-rep PT array |

### 4.4 Purpose in the paper

- v4 family establishes the null result: "At the lower-homoplasy v4 and
  v4-cross parameter regimes (max false-clade P = 0.006 across 22 rescored
  chains), both blind and aware chains recover the true topology near-perfectly,
  confirming that the ecology-aware model does not degrade inference when the
  ecology confound is weak."
- v5break (pending, job 17226503): eco-1 = cross geometry, nNeo = 30,
  nTrans = 50, φ = 6, π₀ = 0.45, stemBrClade = 0.03 — designed to flip
  the blind chain to the false clade; if blind P(falseInner) > 0.1, the
  Goldilocks parameter regime becomes reportable.

[DECISION] Whether v5break results arrive before submission cutoff; if not,
the v4 family is the null-result section and the Goldilocks characterisation
remains future work.

---

## 5. Posterior scoring

### 5.1 Bipartition presence (root-invariant)

- Bipartitions of interest:
  - `trueSister`  = {A1..A4, C1..C4}: present in the true topology
  - `falseSister` = {A1..A4, B1..B4}: the ecology-driven false clade
- Bipartition presence tested via `HasBipartSplits()` in
  `inst/simulations/ecology/sim3-scoring.R`:
  enumerate unrooted bipartitions with `TreeTools::as.Splits()`; a target
  split is present iff some row of the logical bipartition matrix equals the
  target vector or its bitwise complement.
- **Legacy scoring note**: an earlier implementation used `ape::prop.part()`
  (root-dependent; returns descending-tip sets under the current MCMC root
  placement). When roots fall inside a target tip set, the split appears only
  as its complement and is missed. This artefact inflated apparent P(true) and
  suppressed apparent P(false) in 3/8 multirep-v3 replicates (most affected:
  reps 02, 05, 06). All results presented here use the corrected scorer (commit
  `7d1076d` and following). The bug does not affect CID or tree-length
  summaries.

### 5.2 Topology distance

- Normalised Clustering Information Distance (CID; Smith 2020b) to the true
  tree: `TreeDist::ClusteringInfoDistance(trees, refTree, normalize = TRUE)`.
  This is the normalised mutual-information variant of the CID metric; it is
  root-invariant and ranges from 0 (identical) to 1 (maximally distant).
- Per-replicate mean CID computed over post-burnin tree samples.

[DECISION] Confirm that `normalize = TRUE` gives the mutual-information
normalisation (not the sqrt-information normalisation). Check TreeDist
documentation or Smith 2020b for the exact formula underlying `normalize = TRUE`
in version used.

### 5.3 Credibility sets

- 95% credibility set (CS95): smallest set of distinct canonical topologies
  (unrooted, stored via `TreeTools` canonical form) whose cumulative posterior
  probability reaches or exceeds 0.95.
- Used to characterise posterior spread (number of unique topologies in CS95)
  rather than as the primary accuracy metric.

### 5.4 Burn-in

- First 25 % of post-warmup retained samples discarded: `keepFrom = ceiling(n/4) + 1`
  (function `.DiscardBurnin()` in `inst/simulations/ecology/sim3-scoring.R`).

---

## 6. Empirical analysis (rodent morphological matrix)

### 6.1 Dataset

- MorphoBank matrix X24848 (source paper: [CITATION NEEDED — confirm before
  submission; see intro draft note on P5]).
- Raw matrix: 102 taxa × 221 characters (including ecology character col 220
  and extant/extinct flag col 221).
- Taxon filtering: 42 extinct tips excluded; 60 extant tips retained.
- Character classification: `AutoDetectNeomorphic()` applied to the 219
  variable characters (excluding the ecology column); identifies 160 neomorphic
  and 59 transformational characters.
- MkPrimeData input: 219 characters; 2 invariant columns (indices 19 and 217
  relative to pdForDetect) removed automatically → **nChar = 217** characters
  (160 neomorphic + 57 transformational after invariant removal).

[DECISION] The split 160 + 57 = 217 is what rodent-v2.stdout reports ("nChar
217" and "160 | 59"), but 160 + 59 = 219 ≠ 217. The log output says
"Neomorphic chars: 160  | transformational: 59" before MkPrimeData is built,
and "nChar: 217" after — meaning 2 invariants are dropped but the
category-specific counts before dropping are 160 + 59. Verify the exact neo/trans
split *after* invariant removal (i.e. what mkd$type contains) from a completed
run log before finalising the character count table.

- Four ecology states: 0 = terrestrial (reference, 36 tips), 1 = arboreal
  (12 tips), 2 = semiaquatic (5 tips), 3 = fossorial (7 tips). Four tips
  with polymorphic or missing ecology were recoded to their primary state.

### 6.2 Model and MCMC

- Model: `ecologyAware = TRUE`, `magnitudeMode = "global"`, `kPrimePrior =
  "geometric"`, `coding = "variable"`, `expSteps = 10`.
- Ecology-aware hyperpriors: π₀ ~ Beta(7, 3), θ ~ Beta(2, 2),
  σ_φ = 0.5 (`rho0Alpha = 7, rho0Beta = 3, thetaAlpha = 2, thetaBeta = 2,
  sigmaPhi = 0.5`).
- Blind model: same base settings, `ecologyAware = FALSE`.
- MCMC: `nIter = 1 000 000`, single chain (`nChains = 1, nRuns = 1`),
  `thin = treeThin = 1 000`, minWarmup = 10 000, maxWarmup = 50 000.
- Start tree: parsimony optimal via `TreeSearch::MaximizeParsimony()`
  on full character matrix (set.seed(20260512)).
- Both chains resumed from checkpoints; aware chain: 1M total iterations
  (resumed from 371k; 20.68h wall); blind chain: 1M total iterations (19.5 min
  wall, much faster without ecology-aware likelihood).

### 6.3 Convergence target and realised values

- Convergence criterion applied: target minESS ≥ 200 (minimum effective sample
  size over all monitored continuous parameters at the final iteration).
- Realised minESS at 1M iterations: **aware = 88, blind = 38** (both below the
  target; from `.er` log files). Results are therefore pilot-quality.
- Aware continuation run (job 17222514) was still running at the time of
  writing; final comparison will be regenerated once the continuation completes.

[DECISION] Whether to report the 1M-iteration results as preliminary and flag
the ESS shortfall explicitly, or to wait for the continuation run before
describing the empirical analysis in the methods.

### 6.4 Posterior comparison

- Post-burnin: first 25 % of retained tree samples discarded.
- Pairwise CID distance matrix computed on a subsampled set (target N = 200
  per chain) for MDS visualisation.
- MDS: classical metric MDS (`cmdscale()`, k = 2) on the pairwise CID matrix.
- MR consensus: `ape::consensus(post_trees, p = 0.5, rooted = FALSE)`.
  The consensus computation is root-invariant (ape's internal `SHORTwise`
  canonicalisation; verified empirically against a manual `TreeTools::as.Splits`
  approach: 45 aware / 38 blind splits reproduced exactly, shared count = 32;
  commit `956bc6c`).
- Unique-split comparison: splits unique to each consensus identified via
  canonical bipartition matching (`pp_to_canonical_splits()` in
  `inst/scripts/rodent-comparison/rodent-comparison.R`).

---

## 7. Software and reproducibility

- **MkPrime R package** version 0.0.0.9000 (development version; URL:
  `https://github.com/mk-prime/r`, branch `worktree-ecology-aware`).
- Key dependencies: TreeTools (canonical tree operations, `as.Splits`,
  `Preorder`), TreeDist (CID, `ClusteringInfoDistance`), ape (consensus, Newick
  I/O), TreeSearch (parsimony start trees, `MaximizeParsimony`), phangorn
  (phylogenetic data structures).
- All simulation and analysis code is in the MkPrime package repository under
  `inst/hamilton/` (HPC dispatch scripts) and `inst/simulations/ecology/`
  (simulation helpers and scoring).
- Production MCMC runs executed on Durham University Hamilton8 HPC cluster
  (SLURM scheduler; 8-hour walltime per job, resumable via checkpoint files).
- Per-replicate RNG seeds: simulation seed = 20260601 + repId; start-tree seed
  fixed at 1 across all replicates.
- All posterior scoring uses the corrected `HasBipartSplits()` function
  (root-invariant, commit `7d1076d`).

[DECISION] Add CRAN/version pin for TreeDist and TreeTools (the CID normalise
argument behaviour may vary across versions; pin the versions used on Hamilton).

---

*Key sentences pre-drafted (section 1–7 above) and flagged [DECISION] items:*
1. *Which hyperprior set to report for the simulation (pre-R5 vs post-R5 R5; see §1.3).*
2. *Whether to unify simulation and empirical hyperpriors in revision (§1.4).*
3. *Whether to include PT-extended reps in main simulation results (§3).*
4. *Whether v5break results will be available before submission (§4.4).*
5. *Exact TreeDist normalise flag behaviour (§5.2).*
6. *Confirm neo/trans split after invariant removal in MkPrimeData (§6.1).*
7. *ESS shortfall strategy for empirical analysis (§6.3).*
8. *Version-pin TreeDist/TreeTools for reproducibility statement (§7).*
