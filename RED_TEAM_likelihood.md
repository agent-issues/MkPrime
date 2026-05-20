# Red-Team Audit: Mk′ likelihood & prior — why does mk_k40 beat mk on CID?

**Audit date:** 2026-05-19
**Auditor focus:** explain the monotone CID gradient
`mk (kObs) ≫ mk_kp1 ≫ … ≫ mk_k9 ≫ mk_k15 ≫ mk_k24 ≫ mk_k40`
where bigger fixed-k beats the truth-matching baseline by ~0.04 CID.

---

## Top-line conclusion

I found **no off-by-one bug, no padding bug, no recoding bug, and no
identifiable code defect** that would explain the paradox.

Across `mk`, `mk_kp1`, `mk_kp2`, `mk_k9`, `mk_k15`, `mk_k24`, `mk_k40` all
characters travel down the **same `type == "known"` partition code path**
(`src/mcmc_likelihood.cpp:1879-1957`), differing only in the integer `k`
passed to `pruning_jc(...)` / `constant_site_prob_jc(...)`. There is no
extra correction term, no relabelling term, no padding step, and no extra
prior contribution that distinguishes them.

So the gradient is, with high confidence, a **real statistical property of
the JC(k) + variable-coding + Gamma(2, 2/Fitch) tree-length prior** stack,
not a code bug. The finding I keep returning to is that **the same Gamma
prior on tree-length T is being shared across arms whose JC(k) likelihoods
parameterise "substitutions per site" differently as a function of k**.
That mis-anchoring is enough to bias which trees survive MH.

Below I tag each suspect and rule it in/out with file:line.

---

## NOTE-1 (most likely driver of the paradox): tree-length prior is k-invariant, but the JC(k) likelihood is not

**Files:**
- `R/MkPrimeModel.R:241-246` `.FinalizeModel()` sets
  `expSteps = max(1, .FitchScore(tree, mkd))` and
  `treeLengthRate = 2 / expSteps`. The Fitch score is computed once
  from the **data** (it uses `mkd$kObs`, not `known_k` or `kPrime`),
  so for a given dataset every arm gets the **identical**
  Gamma(2, 2/Fitch) prior on the total tree length T.
- `src/likelihood.cpp:80-86` (and parallel sites in
  `constant_site_prob_jc` at `src/ascertainment.cpp:66-73`) parameterise
  the JC(k) transition probability as
  `exp(-k·t / (k−1))`. The branch-length unit therefore corresponds to a
  **different** expected number of substitutions per site as k grows:
  - k = 2 → rate exponent `2t / 1 = 2t`
  - k = 40 → rate exponent `40t / 39 ≈ 1.026t`

**Why this explains the gradient.** The Gamma prior anchors **T**, but T
means different things in the two arms. For a fixed posterior-favoured
T, mk_k40 effectively allows the chain to absorb the observed tip
disagreement with roughly half as many substitution events as mk
(because each unit-T branch produces ~½ the substitution probability at
k = 40 vs k = 2). When LBA tempts the chain into placing long branches
to explain conflicting characters, the high-k model resists that more
strongly — the same T sits in a less changeful region of the JC
probability simplex. The high-k model therefore **acts as a smoothing
regulariser** that suppresses LBA-induced topology errors, at the cost
of biasing tree-length estimates downward.

This matches the **monotone** pattern you observe and the **direction**
(better CID at higher k) — and it predicts that the gain saturates once k
is large enough that `k/(k−1) ≈ 1`, which the k=24 vs k=40 plateau in the
data would be consistent with.

**Suggested test (cheap).** Re-run mk and mk_k40 with the **tree-length
prior re-scaled per arm** so that `E[T] = expSteps · (k−1)/k` (or
equivalently re-fit `treeLengthRate = 2k / ((k−1)·expSteps)`). If the
CID gap collapses, this is the cause.

---

## NOTE-2: variable-coding ascertainment is k-dependent in a way that also favours high k

**Files:**
- `src/mcmc_likelihood.cpp:1947-1955` (known partition, JC path):
  `knownConstProb = constant_site_prob_jc(... kStates ...)` then
  `ll -= part.tipStates.ncol() * log(1.0 - knownConstProb)`.
- `src/ascertainment.cpp:46-100` computes P(constant) on `kStates`
  symmetry-folded pseudo-characters.

For a binary character treated under JC(40), the model has 40 equally
likely root states; the probability of ending up with **all** tips in the
**same** state is much smaller than under JC(2). So `log(1 − P_const)` is
closer to 0 under JC(40), i.e. the ascertainment subtraction is smaller,
i.e. raw log-likelihoods get an extra positive nudge at high k.

This is **not** a constant offset that cancels in MH ratios for tree
moves — `P_const` depends on the **branch lengths** of the current tree
through `constant_site_prob_jc`. So this term shapes the posterior
over trees, again pushing in the direction of "less changeful" tree
configurations under high k.

I do **not** think this is a bug — Lewis (2001) coding is correctly
implemented (`ll -= nChar · log(1−P_const)`, character-wise). But it is
a second mechanism by which high-k models behave differently from the
truth-matching mk arm despite being "wrongly specified".

---

## Ruled out: off-by-one at the k′ = kObs boundary (CRITICAL/WARNING checked, NOT FOUND)

- `R/MkPrimeData.R:137-145`: validates `knownStates ≥ kObs`. Boundary
  case `knownStates == kObs` is **allowed** and produces a `known`
  partition with `k = kObs`. This is exactly what the mk arm uses (line
  `data-raw/hamilton/run_one.R:216` passes `knownStates = kObs_for_mk`).
  No off-by-one.
- `R/partition.R:60-78`: known characters partition by `known_k`, with
  partition `k = kv` and `kObs = max(mkd$kObs[sel])`. Correct.
- `src/mcmc_likelihood.cpp:1881` `int kStates = part.k;` — this is the
  user-supplied k (40, 24, …, or kObs), passed directly to `pruning_jc`
  with `rootFreqs(kStates, 1.0/kStates)`. Stationary frequencies are
  `1/k` on each of the k states. No off-by-one.

If there were an off-by-one (e.g. k passed as k−1 or k+1 somewhere), I'd
expect a non-monotone or asymmetric pattern; the smooth monotone ramp
is a positive signal that the dispatch is correct and we are seeing a
genuine modelling property.

## Ruled out: relabelling correction leaking into the mk_k* arms

- `src/corrections.cpp:35-44` `mk_prime_relabel_log` is only invoked from
  `mcmc_likelihood.cpp:2047-2050` and `:2146-2148`, both inside the
  `type == "transformational"` branch (`pinfo.type == 1`). The
  `type == 2` known branch at `mcmc_likelihood.cpp:1879-1957` does **not**
  add a relabel term.
- `prepare_mcmc_data` at `mcmc_likelihood.cpp:2247-2249` only populates
  `transIdxGlobal` for `charTypes_r[i] == "transformational"`. Since
  every variable character in mk and mk_k* is `known`, there are zero
  transformational characters and the relabel block is never reached.

So the relabel correction (which would have been a smoking gun — it
grows like `lgamma(k+1) − lgamma(k−kObs+1)` and would favour high k
trivially) is **NOT** present in any of these arms. Both mk and mk_k40
get raw JC log-likelihood + Lewis ascertainment correction, no relabel.

## Ruled out: differing MCMC move sets between arms

- `R/RunMkPrime.R:3106-3157`: kPrime moves (`int_walk`,
  `gibbs_kPrime`, `block_kPrime`, plus the `p` Gibbs move or hyperparameter
  moves) are gated on `nTrans > 0`. In mk and mk_k* arms, **all variable
  characters are `known`**, so `nTrans = 0` and **no kPrime moves are
  added**. Identical move set across mk, mk_kp1, …, mk_k40.
- Topology/branch-length/rate-heterogeneity moves are constructed
  identically (line `R/RunMkPrime.R:3004-3010` and `3170-3173`).
- Acceptance ratios will of course differ because likelihoods differ,
  but the *menu* of moves is the same. No "mk is stuck in a mode that
  mk_k40 can escape" via a move-availability route.

## Ruled out: state-space padding/leakage

- `R/MkPrimeData.R:294-305` `.PhyDatToIntMatrix` remaps each character's
  states to **contiguous 0-based integers** based on observed values
  only. So a binary character coded {1, 2} in the input becomes {0, 1}
  in the tip-state matrix, and the matrix never carries a state index
  ≥ kObs.
- `src/likelihood.cpp:57-72` initialises tip CLs with `1.0` at exactly
  one column out of `kStates`. For mk_k40 with a binary char,
  states 0 and 1 see CL=1, the other 38 see CL=0. The 38 latent
  states then accrue probability mass only through internal-node
  transitions. **This is correct JC(40) on data restricted to 2 observed
  labels — not a leak, just standard latent-state marginalisation.**
- Ambiguous tips set all `kStates` columns to 1.0 (line 62-65). Also
  correct.

## Ruled out: tree prior differs across arms

- Both mk and mk_k40 call `MkPrimeModel(coding = "variable")` with
  no override of `treeLengthShape`/`treeLengthRate`/`expSteps`. The
  Fitch score depends only on `mkd$kObs[j]` (line `R/MkPrimeModel.R:283`),
  which is identical across arms (kObs comes from the data). Therefore
  `treeLengthRate = 2 / FitchScore` is identical across arms. Confirms
  NOTE-1: the prior is the **same**, but it means **different things** in
  units of "expected substitutions per character per unit time" because
  the rate normaliser `k/(k−1)` differs.
- `rateLogSdShape = rateLogSdRate = 1` Gamma prior on ACRV — same.
- The branch-length simplex (Dirichlet(1,…,1)) prior — same.

## Ruled out: knownStates handling errors

- `R/RunMkPrime.R:117-123`: when `data` is a `phyDat`, `RunMkPrime`
  forwards `knownStates` to `MkPrimeData()`. If `data` is already
  `MkPrimeData`, `knownStates` is ignored (line 118-119) — and the mk
  arms always pass a phyDat (`data-raw/hamilton/run_one.R:213`). No
  silent drop.
- `R/MkPrimeData.R:69-159`: validates names, parses indices, drops
  invariants and **remaps the knownStates indices** (lines 112-122),
  validates `knownStates ≥ kObs` (137-145). Sound.

---

## Other findings unrelated to the paradox but worth flagging

### WARNING-1: empirical_geometric prior — moot for this paradox but live for mkp_eg

The `LogPrior` function for `empirical_geometric` (`R/MkPrimeModel.R:476-497`)
does **not** include a per-character truncation normaliser — see also the
`mkp_geo` arm comment in `data-raw/hamilton/run_one.R:344-347`, which
flags this as "EG-001 suspect normaliser". Not relevant to the mk vs
mk_k40 question (neither arm uses k′ inference), but it is a known live
issue for the mkp_eg arm and the user already has a pilot to test it.

### WARNING-2: `rate_log_sd == 0` boundary handling

`R/MkPrimeModel.R:454-462`: when `rate_log_sd == 0` and `shape == 1`,
`dgamma(0, 1, rate, log=TRUE)` is `log(rate)` (finite); the code skips
the `dgamma` call entirely in this branch. With the default
`rateLogSdShape = 1`, this means the prior contribution at the boundary
is treated as 0 rather than `log(rate)`. Tiny constant — does not affect
the paradox — but is technically a prior-normalisation inconsistency.

### NOTE-3: thinning differs

`data-raw/hamilton/run_one.R`: `mk` and `mk_kp1` use `thin_iters = 10L`
(default), while `mk_kp2`, `mk_k9`, `mk_k15`, `mk_k24` (and `mkp_geo`)
use `thin_iters = 100L`. Identical wall-clock budget therefore yields
~10× more samples on the mk arm, which **should help mk if anything**.
This argues the gradient is real and not a statistical-power artefact.

---

## What I would do next

1. **Test NOTE-1 directly.** Re-run mk_k40 with
   `treeLengthRate = 2 · k / ((k−1) · FitchScore)` so that the prior puts
   the same mass on "expected substitutions per site" as the mk arm does.
   If the CID gap collapses, NOTE-1 is the cause and the result is
   essentially "high-k arms get a free, advantageous prior on T".

2. **Sanity check via a known-tree simulation under JC(2).** Generate
   characters under JC(2) on a known tree with **no rate heterogeneity
   and no LBA-prone branches**. Both mk and mk_k40 should now perform
   similarly. If they don't, the JC(k) numerics themselves are the
   culprit and the paradox is not about LBA regularisation.

3. **Inspect posterior tree lengths per arm.** If NOTE-1 is right, the
   posterior `T̂` should be a roughly `k/(k−1)`-scaled version of the
   mk arm's `T̂`. That is the direct fingerprint.

4. **Re-derive the Lewis correction algebraically for JC(k) with
   kObs < k.** Confirm that `constant_site_prob_jc` correctly sums over
   *all* k constant patterns (including the 38 latent ones), not just
   the kObs observed ones. Reading `src/ascertainment.cpp:90-100` it
   does (the M-170 symmetry argument multiplies by `kStates`, which is
   correct), but this is the easiest place to be off and would amplify
   exactly with k.

---

## Files and lines cited (absolute paths)

- `C:\Users\pjjg18\GitHub\mkp\R\MkPrimeModel.R:241-246, 268-305, 391-532`
- `C:\Users\pjjg18\GitHub\mkp\R\MkPrimeData.R:47-178, 257-308`
- `C:\Users\pjjg18\GitHub\mkp\R\partition.R:18-81`
- `C:\Users\pjjg18\GitHub\mkp\R\RunMkPrime.R:76-188, 2999-3175`
- `C:\Users\pjjg18\GitHub\mkp\src\corrections.cpp:35-72`
- `C:\Users\pjjg18\GitHub\mkp\src\likelihood.cpp:30-100`
- `C:\Users\pjjg18\GitHub\mkp\src\ascertainment.cpp:32-101`
- `C:\Users\pjjg18\GitHub\mkp\src\mcmc_likelihood.cpp:1879-1957, 1958-2153, 2185-2280`
- `C:\Users\pjjg18\GitHub\mkp\data-raw\hamilton\run_one.R:203-490`
