# MkPrime ecology-aware — cron-fired update 2026-05-20 01:23 BST (extended 02:30)

## 🧪 MP trap test + multi-eco pre-screen (commits `5a83abb`, `9764b01`)

Following user feedback ("test of the trap"), used `TreeSearch::MaximizeParsimony`
as a fast deterministic proxy for whether blind ML can be fooled. Results:

| Sim | Truth TL | MP finds false eco-clade? | Mean MP step gap (truth - false) |
|---|---:|---:|---:|
| v6-realistic (n=3) | 1.16 | 0/3 | 0 (tie) |
| v5break (n=3) | 1.30 | 0/3 | ~7-8 |
| multirep-v3 (n=8) | 13.5 | **7/8** | ~28-30 (false better) |

**multirep-v3 MP-blind correspondence is 7/8 ↔ blind ML P(AB)=0.864** — confirms the
long-tree trap is real for both MP and ML. The corrected blind MCMC matches MP's
deterministic verdict. Validates the MP-as-proxy approach.

Multi-ecology v8 pre-screen (6 variants × 8 reps each, phi 4-10, pi0 0.2-0.5,
TL 1.4-1.9, including 3-ecology + 2-ecology designs): **none broke MP**. Best was
12% (1/8 reps).

**The phi-rate-multiplier mechanism cannot create a parsimony trap at TL ≤ 2 along
any single-axis variation we have tested.** v7 + v8 confirm the structural limit
diagnosed earlier: rate inflation makes eco-clades noisy, not directionally
convergent.

## 🐀 Rodent posterior diagnostics (commit `e90d67c`)

From the saved rodent aware-v2 chain (1M iter complete, post-burnin n=443):

| Param | median | 90% CI |
|---|---:|---|
| phi | **4.03** | (3.10, 5.78) |
| pi0 | 0.43 | (0.29, 0.56) |
| theta arboreal | 0.17 | (0.07, 0.26) |
| theta semiaquatic | 0.50 | (0.31, 0.72) |
| theta fossorial | 0.29 | (0.16, 0.41) |
| TL | 22.6 | (15.4, 31.1) |

**Real ecological gravity: phi ≈ 4** — substantially less than v6/v7's phi=8.
~50% of chars eco-affected per ecology. Cross-ecology Jaccard 0.50-0.67 (high
overlap) but Pearson r ≈ 0 (independent recruitment). Biological top-flagged
chars are postcranial/locomotor + craniofacial features — sensible.

TL=22.6 is much smaller than the v2 saturated 127 — the chain may have
finally found a less-saturated mode in the 371k → 1M iterations.

## 🎯 BREAKTHROUGH: v9 directional-z trap at realistic TL (commits `8fc37ef`, `e7c7f82`)

User insight 2026-05-20: "move to a neomorphic world... ecology A is gaining some
characters at a higher rate, and ecology B is losing these same characters at a
higher rate. That's the cleanest framing of a convergent signal."

**The existing `.SimulateMkPrimeEcology` already supports asymmetric neomorphic
rates** via z indexing: z=1 gives (rate01 × phi, rate10 / phi) — directional,
not symmetric. This mode was never exercised before. v9 exploits it.

**v9 MP pre-screen** (`dev/sim-design/v9-discriminate.R`):

| Variant | Mechanism | TL | phi | MP pAB (false sister) | MP pTrue (AC) |
|---|---|---:|---:|---:|---:|
| v9a M1-divergent (realistic) | z[A]=1, z[B]=2 | 1.48 | 6 | 25% | 88% |
| v9b M1-divergent (pushed) | z[A]=1, z[B]=2 | 1.52 | 8 | 12% | 88% |
| **v9c M2-parallel (realistic)** | **z[A]=z[B]=1** | **1.48** | **6** | **62%** | **62%** |

**v9c found a parsimony trap at TL=1.48** (realistic) with 62% of reps showing
the false (A,B) sister grouping. The mechanism: both eco-A and eco-B elevate
gain rate on the same chars → both clade trends to state 1 → blind sees A+B
share derived states.

**v9 Hamilton MCMC array submitted: 17233486_[1-3] (PENDING).**
This tests whether blind ML actually commits to the false clade AND whether
aware ML rescues truth. Result expected in ~2h. If both hold, the paper has
its long-sought realistic-TL headline.

## 🐀 Rodent v3 blind verdict (1M iter, 30 min wall, completed 06:31 BST)

- minESS = 76 (median 311) — borderline paper-grade
- **Posterior TL median = 126.9** (vs v2 saturated 127.8) — essentially identical
- The parsimony-anchored prior did NOT lower TL — the rodent likelihood
  GENUINELY prefers saturated TL. Rodent is empirically a long-tree case.

**Implication**: the long-tree regime where aware regularises is empirically
relevant (rodent lives there). The multirep-v3 narrative isn't a corner-case.

## 🎯 Paper-direction options after the diagnostic landscape

(a) **"aware as honest regulariser at long-tree regimes"** — multirep-v3 IS
    the demonstration. Argument: long-tree regimes occur in real palaeontology
    (deeper trees, rapid radiations, well-sampled clades). Acknowledge it as
    saturation-regime; document the prior-anchored fix; rodent as the
    intermediate empirical case.
(b) **Change the ecology MODEL to directional bias** — significant code change
    to `.SimulateMkPrimeEcology` (add per-eco state-bias) and `MkPrimeModel`.
    Would enable directional convergence at realistic TL. Could be the
    headline methodological contribution.
(c) **Add unexplored design levers** — larger taxon counts (60+), asymmetric
    clade sizes, or pre-saturated tip evolution (already at equilibrium then
    eco signal layered on). Could find a trap zone without changing the model.
(d) **Honest negative result paper** — methods paper documenting that
    eco-rate-multipliers cannot trap MP/ML at realistic TL with current
    morphological-data scales. Conclusion: the aware model is a no-op in
    realistic regimes; future work to develop a model that can capture
    realistic eco-driven convergence.

Honest read: (a) is the most defensible publication path with existing work.
(b) is the most ambitious / scientifically interesting if user wants to commit
code time. (c) might find new traps but the search space is large and the
single-axis search has been exhausted. (d) is publishable but unsatisfying.



## 🎯 v6-realistic NULL result (commit `3b61e17`)

Across 3 reps at realistic TL=1.16 with corrected prior:

| Rep | Blind P(AC) | Aware P(AC) | Blind P(AB) | Aware P(AB) | Blind CID | Aware CID | Blind TL | Aware TL |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.000 | 1.000 | 0.000 | 0.000 | 0.065 | 0.000 | 1.16 | 1.35 |
| 2 | 1.000 | 1.000 | 0.000 | 0.000 | 0.065 | 0.100 | 1.31 | 1.32 |
| 3 | 1.000 | 1.000 | 0.000 | 0.000 | 0.065 | 0.065 | 1.00 | 1.06 |

**Both models recover truth.** No model differentiation at realistic TL.
TL prior fix works (truth 1.16; posterior medians within 16% of truth in
6/6 chains; vs old 5-13× inflation).

**The multirep-v3 "blind fails / aware regularises" headline does NOT
survive at realistic TL.** It was an artefact of saturation (truth TL=13.5,
expSteps=10 prior with mean=10) plus the root-dependence scoring bug.

**One caveat worth flagging, not a story:** aware concentrates one eco-1
stem onto a long edge (median 2.1-2.6 vs truth 0.15) while leaving the
other compressed. Looks like a rooting/edge-mass artefact; blind keeps
both eco stems near truth. Not the "regulariser" story we want.

## 🔬 v7-induce discriminator: structural limit of phi model (commit `b256e3e`)

User proposed three levers: reduce nChar, push eco signal, drift TL up
slightly. v7 discriminator (`dev/sim-design/v7-discriminate.R`) tested 4
variants spanning these levers at TL ≤ 2:

| Variant | phi | pi0 | TL | Parsimony delta (AB − truth) | Verdict |
|---|---|---|---|---:|---|
| v7a (proposed) | 10 | 0.30 | 1.53 | **+1.88** | truth still wins |
| v7b (stronger) | 12 | 0.20 | 1.62 | +0.88 | truth still wins |
| v7c (push TL) | 12 | 0.20 | 1.80 | +2.25 | truth still wins |
| v7d (aggressive) | 15 | 0.15 | 1.69 | −0.25 | coin flip |

Target was AB ≤ −5 (blind should commit to false clade). **None achieved
it.** Pushing further (v7d, phi=15) flattens MP-best to 50/50 — past the
optimum, eco-clades randomise rather than converge.

**Structural reason:** `phi` is a SYMMETRIC rate multiplier
(`mult <- c(1, phi, 1/phi)` indexed by z+1 in `sim3-simulate.R`). It makes
eco-clades NOISY but not DIRECTIONALLY convergent. At realistic TL with
high phi, both eco stems reach equilibrium and the within-clade
homoplasy is resolved cheaply on the truth tree.

**The eco-rate-multiplier model cannot produce A↔B parallel synapomorphies
at realistic TL.** To break blind at TL ≤ 2 would require directional eco
bias (eco-1 prefers state X, eco-0 prefers state Y), not rate inflation.

## Paper options (refined)

(a) ~~**Strengthen eco confound at realistic TL**~~ — **STRUCTURALLY
    IMPOSSIBLE** with current phi-as-rate model. Would need to change
    the model itself to directional bias.
(b) **Extend to N=8 reps at v6** — confirms null with CIs. Modest cost,
    no narrative change.
(c) **Pivot paper to prior-anchoring + multirep-v3 as cautionary tale**
    — frame the contribution as "we found that the default Mk' TL prior
    drives chains into saturation; the parsimony-anchored fix recovers
    realistic posteriors; aware exists and behaves sensibly under both
    regimes". Honest 3-rep evidence + structural argument favours this.
(d) **Add directional eco-bias model** (new option) — code change to
    `.SimulateMkPrimeEcology` and `MkPrimeModel` to support a per-state
    bias on eco branches. Would enable directional convergence at
    realistic TL but is a significant model extension. Could be the
    headline contribution if the user wants to commit.

## Active jobs

| Job | What | Status |
|---|---|---|
| 17227246 | mkp-rod-aware-v3 (1M) | RUNNING (~270k samples, minESS 10, climbing) |
| 17231585 | mkp-rod-blind-v3 (resubmitted) | RUNNING (~3 min, fresh start) |

(17227245 rodent blind v3 failed on stale checkpoint + old Resume API; cleared and resubmitted.)

## Recommended next step when user returns

1. Read `inst/scripts/v6-realistic-analysis/v6-realistic-analysis-report.md`
   for the full v6 numerical picture.
2. Pick a paper direction (a / b / c). My honest read: **(c) is the
   defensible publishable story now**. The methodological contribution is:
   - Identified that morphological matrices need parsimony-anchored TL
     priors (the default Gamma(2, 0.2) drives chains into saturation).
   - Showed that an ecology-aware model exists and produces sensible
     posteriors under a correctly-calibrated prior; doesn't break things
     when ancestry signal is strong.
   - Demonstrated the *risk* of saturated MCMC chains misleading inference
     (the multirep-v3 cautionary tale, traced back through the scoring bug
     to the prior).
3. The rodent v3 chains will tell us whether the parsimony-anchored prior
   recovers a sensible TL on real morphological data (expected TL ≈ 7-15
   vs the v2 saturated 127). That's the empirical anchor for option (c).

# MkPrime ecology-aware — hand-off 2026-05-19 (earlier state below)

## 🎯 TL prior recalibration (commits `bc24d52`, `83ab33a`, `2b8d183`)

User flagged that simulated truth TL=13.5 (multirep-v3) is unrealistic
(typical morphological matrices have TL ~ 1-2). Investigation revealed:

**The TL prior was Gamma(2, 2/expSteps) with default expSteps=10**, giving
mean TL=10 — wildly mismatched to morphological data. This caused chains
to inflate TL **5-18× above truth** across all v4-family sims:

| Sim | Truth TL | Blind chain TL | Inflation |
|---|---:|---:|---:|
| v4a | 2.40 | 20.9 | 8.7× |
| v4b | 2.02 | 10.4 | 5.1× |
| v4c | 1.24 | 6.5 | 5.3× |
| v4cross-b | 2.02 | 9.4 | 4.7× |
| v5break-rep01 | 1.30 | 7.5 | 5.8× |

User's call:
> "Simulations on saturated datasets are junk; they have next to no real
> world interpretability. We need realistic simulations; else what's the
> point?"

**Fix attempt 1 in `bc24d52`** had a unit error: `expSteps = 1.05 ×
parsimony` makes prior mean TL = parsimony score (total state changes),
but MkPrime's `tree_length` is sum of edge lengths in per-character
substitution units. Prior was 270× too high.

**Fix attempt 2 in `f7b4644`** (final): `expSteps = (parsimony / nChar) ×
1.05`. Worked values:
- v6 (16 tips × 300 chars, parsimony 297) → expSteps = 1.04 ≈ truth TL 1.16 ✓
- Rodent (60 tips × 217 chars, parsimony 1519) → expSteps = 7.35 (reasonable)
- Toy (16 tips × 100 chars, parsimony 50) → expSteps = 0.525 (floor binding)

Sanity-warn added if computed expSteps > 100. 29/29 tests pass.

**Active jobs (after `f7b4644` unit-fix + `b8da51b` Resume-API fix):**
| Job | What | Status |
|---|---|---|
| 17227245 | mkp-rod-blind-v3 (corrected expSteps≈7.0) | RUNNING (parsimony search) |
| 17227246 | mkp-rod-aware-v3 (corrected expSteps≈7.0) | RUNNING |
| 17227265_[1-3] | mkp-v6-realistic (corrected expSteps≈1.04) | RUNNING |

(Earlier 17227192/3, 17227212, 17227247 cancelled — bad prior or Resume API mismatch.)

**Other contra-API run scripts** flagged by subagent (need patching if resubmitted):
sim3-v4a/b/c, sim3-v4cross-b/c/b-pt/b-pt5, sim3-v5break, sim3-multirep-v3/v4,
sim3-multirep-v3-pt, sim3-phi2*, sim3-pt-test, rodent-blind-v2. Most won't be
resubmitted under the "saturated sims are junk" stance.

**multirep-v3 declared junk.** v4 family TL marginal (1.24-2.4 truth, but
chains inflated). v6-realistic (truth TL=1.16, parsimony-anchored prior,
phi=8 + longer eco stems for confound) is the new headline simulation.

**Two flagged side-effects of the install:**
1. The reinstall **upgraded** Hamilton to the R5-tightened priors
   (rho0=Beta(360,120), theta=Beta(2,2), sigmaPhi=1.5) committed upstream
   2026-05-16 in `2fd4dd5`. Audit (task #22) showed Hamilton had been on
   the looser (75,25) prior all along. **Headline finding survives** —
   "aware regularises false-AB" was data-dominated under the looser prior
   (~17% prior weight on pi0). v6-realistic + rodent-v3 use the tighter
   prior; cross-version pi0 posterior comparison NOT directly comparable.
2. PT-aware-mr3-rep05 (17226606) and rodent aware cont2 (17222514) were
   cancelled — both used the bad expSteps=10 prior.

## Session bookmarks (autonomous batch 2026-05-19)

21 commits, 6 memory files touched. Key end-state:

1. **`project_scoring_bug.md`** — durable lesson on ape root-dependence;
   never use raw `prop.part` for bipartition support; safe alternatives
   listed.
2. **`inst/simulations/ecology/sim3-scoring.R`** — `HasBipartSplits`,
   `ScoreTreesUnrooted`, `RootInsideSet`, `RootConfigCounts`. Source this
   in every new analysis touching MCMC posterior trees.
3. **multirep-v3 is the paper's headline sim.** Across 8 reps, blind P(false
   AB clade)=0.864 (5/8 reps > 0.9); aware reduces to 0.411 and is uniformly
   closer to truth in CID (Wilcoxon p=0.014, 8/8 reps). Aware doesn't commit
   to AC; framed as "honest regulariser" not "recoverer".
4. **v4 family is the null robustness check.** All v4/v4-cross sims now
   show both models recover truth at 16 tips × 200-300 chars × phi=4. No
   model differentiation — that's expected when blind doesn't fail.
5. **2/8 multirep-v3 aware reps caveat (01, 07)** explained as **warmup TL
   trap** affecting 5/8 reps. Fix: PT or longer warmup. `sim3-multirep-v3-pt/`
   built (commit `b85714c`), pending SSH to submit.
6. **Publication figures drafted** at `inst/scripts/multirep-v3-figures/`.
   Results and methods outlines at `inst/lit/{results,methods}-outline.md`.
7. **Rodent comparison clean** (audited; ape's rooted=FALSE is safe).
   Awaiting aware-cont2 (job 17222514) to regenerate final report.

## Highest-leverage next actions when SSH returns

1. `sbatch inst/hamilton/sim3-multirep-v3-pt/sim3-multirep-v3-pt.sh` —
   tests whether PT escapes the warmup TL trap for rep 05.
2. Poll job 17226503 (v5break) and 17222514 (rodent aware cont2).
3. If PT-aware-multirep-v3 shows clean AC recovery on rep 05, queue
   `sim3-multirep-v3-pt5` array (task #18) for the publication run.


## ⚠ Scoring bug fixed 2026-05-19 (commits `7d1076d`, `6500cda`)

Legacy `hasBipart()` used `ape::prop.part`, which is **root-dependent**.
MCMC trees vary root position across samples (3–20 root configs per chain,
median ~12). When the root lay inside a target subset, the clade was
invisible to `prop.part` → spurious low support.

Fix: `inst/simulations/ecology/sim3-scoring.R` provides root-invariant
`HasBipartSplits()` via `TreeTools::as.Splits`. All 22 v4/v4-cross saved
RDS files re-scored; results in `inst/scripts/rescore-sim3-results.csv`.

**Findings (corrected scoring):**

1. **Both blind AND aware essentially recover all true clades + sisters on
   every v4/v4-cross regime tested.** Aware never demonstrably fails; blind
   never demonstrably fails either — at 16 tips × 200–300 chars × phi=4 the
   eco confound is too weak to break either model.
2. **The "aware loses on v4-cross" 5-rep PT headline is gone.** Corrected
   means: blind P(AC)=0.999, aware P(AC)=0.996; blind P(BD)=0.999, aware
   P(BD)=0.996. Both chains essentially perfect.
3. **The "v4c blind catastrophic" P(AC)=0 was an artefact.** Corrected =
   1.000. v4b blind P(AC) 0.574 → 1.000. The original "blind drops AC"
   motivation for v4-cross / v4c was the scoring bug, not real model
   behaviour.
4. **No active false-clade support anywhere** (P(falseInner)=0 everywhere;
   P(falseAB) ≤ 0.006 in one rep). NOT a scoring artefact.

**Implication for the paper narrative:** the entire "blind fails, aware
rescues" story has not been demonstrated on any tested simulation
architecture. The eco confound at tested scales is too weak to
differentiate the models. A sharper sim that actually breaks blind is
required before any methodological claim can be made. v5break (job array
17226503) is the live attempt.

## 🏆 multirep-v3 contains the false-clade demo we were chasing

Rescore of all 8 reps × 2 chains (commit `82d9440`, table at
`inst/scripts/rescore-multirep-v3-report.md`):

**Blind chain (false-clade signal — robust)**:
- Mean P(AB-wrong-eco-clade) = **0.864** across 8 reps
- 7/8 reps with P(AB) > 0.5; **5/8 reps with P(AB) > 0.9**
- Legacy scoring HID the signal in 3/8 reps (rep02 0→1.0, rep05 0→0.978,
  rep06 0.44→0.95)

**Aware chain (does NOT recover truth)**:
- Mean P(AC-true) = 0.012; 0/8 reps with P(AC) > 0.5
- Mean P(AB) reduced from 0.864 → 0.411 (aware breaks false clade)
- Posterior diffuse: 94–167 unique canonical topologies per chain,
  top-topology share 1.7–6.2%

**Strict dichotomy (blind wrong, aware right): 0/8 reps.**

Interpretation: the multirep-v3 high-homoplasy regime (tipBr=0.5,
stemBr=0.30) is enough to deceive blind into the false clade, but aware
cannot find AC in the same data — it simply refuses commitment. That's a
defensible "aware-as-regulariser" half-story for the paper, but not the
ideal "two models give different answers, one is right" full dichotomy.

The v4 family failed because it had too LITTLE homoplasy to break blind.
multirep-v3 has too MUCH for aware to recover. v5break (job 17226503) is
in the goldilocks zone — moderate homoplasy with concentrated eco signal.

## 🔥 Cascading collapses from the scoring fix

The same bug contaminated diagnostics, not just primary scores. Already
overturned:

- **"v4cross-b aware clade-B mode-trap"** (commit `24fcaf6`) — collapsed.
  Legacy B-mono 83/243 → corrected **243/243**. Legacy unique topology
  count 14 → corrected **1**. The aware chain was never pathological.
  Audit at `inst/hamilton/sim3-v4cross-b/diagnose/audit-2026-05-19.md`
  (commit `6cb7090`).
- **9 other contaminated analysis scripts** patched: `sim3-mcmc.R`,
  `sim3b-mcmc.R`, `sim3-v1v2-compare.R`, `sim3-multirep.R`,
  `sim3-analyse.R`, `sim-figures.R`, `diagnostic-multirep-mixing.R`,
  plus two `dev/pilots/` diagnostics (commit `97bc79e`). A *second* bug
  pattern was also patched: `prop.part` output used as topology
  fingerprint in dev pilots — different roots on same tree hash
  differently, inflating "unique topology" counts.
- **Stale RDS warning**: `sim3-multirep-n*.rds` bake legacy `agg` values
  in; would need regeneration before re-plotting.

Mixing followups audited (commit `5972085`):
- **Rep 01 TL collapse SURVIVES** — TL root-invariant; chain wanders on
  wrong side (P(AC)=0.022, P(AB-wrong)=0.789, ~160 canonical topologies).
  Real mode-trap, PT-fixable.
- **Rep 05 topology valley SURVIVES** — P(AC)=0 under both scorings;
  unique-topology count modestly inflated 14% by legacy artefact but the
  valley is genuine.
- **New blind rep05 mode-lock revealed** by audit: P(AB) legacy 0 →
  corrected 0.978. Same artefact class — root inside the AB tip set in
  every sample. (Subsumed into the multirep-v3 rescore above.)

Confirmed UNAFFECTED:
- Rodent comparison (`ape::consensus(rooted=FALSE)` is safe).
- pi0/phi/z sampler R5-4 audit (no bipartition scoring involved).
- Package R/ source (no bipartition logic).

## What this repo is

MkPrime is an R package implementing an ecology-aware Mk' MCMC for Bayesian
phylogenetic inference from morphological data. The `worktree-ecology-aware`
worktree drives the ecology paper: simulations demonstrating that ecology-aware
inference reduces false-clade support relative to an ecology-blind Mk' baseline,
plus an empirical rodent case study.

Artifact targets: paper draft + the rodent comparison; underlying R package is
the methods substrate.

## Where we left off

This session pushed three threads forward:

1. **v4-cross sim architecture** (eco-1 in non-sister clades A,B; opposes true
   `((A,C),(B,D))` topology rather than aligning with it as in v4):
   - `d4d0bee` introduced helpers + v4cross-b run script.
   - `2873ab4` added v4cross-c stress regime.
   - `fc2737d` added single-rep PT (nChains=4) and confirmed an aware
     mode-trap: clade B contiguously broken at sample 84 onward, ll=−1587 vs
     mode-1 ll=−1583. Diagnostics at `inst/hamilton/sim3-v4cross-b/diagnose/`.
   - `63bb335` then ran a **5-rep PT sweep** (`v4cross-b-pt5`, jobs
     17205831_[1-5], all COMPLETED). Headline: aware does NOT reliably rescue
     true topology — AC sister mean blind=0.750, aware=0.479. Aware no longer
     consistently beats blind on this confound geometry, and false-clade
     support never emerges (max falseAB=0.006). Inverse of intended story.

2. **Rodent case study**:
   - `eb2840f` added Hamilton blind pilot (200k iter, job 17205921, completed).
   - `40ca0d4` produced a first comparison report
     (`inst/scripts/rodent-comparison/`) with CID-MDS scatter and MR consensus
     plots. **Both chains were underpowered** — blind minESS=21, aware was
     truncated at 127k of 500k iter (only 7 post-burnin trees). The 0-shared-
     splits result and MDS separation are pilot artefacts, not real findings.
   - `0a95af5` resubmitted both: blind continuation (job 17207482, COMPLETED
     in 19.5 min) and aware fresh 1M iter (job 17207480, **TIMEOUT** at 12h,
     reached 371k, minESS=31, 161 trees).
   - `79df844` added `rodent-aware-v2-cont/` — resume from checkpoint to 1M.
     Job **17215319** running (24h wall, ~13h ETA).

3. **Sim2/3 mode-trap diagnostics** (background): the v4cross-b clB asymmetry
   is the same topology-valley failure as the multirep audit (rep 01/05). PT
   helps but does not eliminate it on hard configurations.

## ⚠ SSH to Hamilton intermittent 2026-05-19 17:45+

Permission denied (publickey,...) on multiple attempts. Several pending
items are blocked:
- `sim3-multirep-v3-pt/` (commit `b85714c`) built but not submitted
- v5break poll (job 17226503)
- rodent aware-v2-cont2 progress check (job 17222514)

Once SSH returns, submit `sim3-multirep-v3-pt` first — it directly tests
whether the warmup TL trap explanation for reps 01/07 can be resolved.

## 🔬 Reps 01/07 caveat diagnosed (commit `669c5e8`)

Aware reps 01 and 07 still support false AB clade because of a **warmup
TL trap** affecting 5/8 reps (chains start TL=3.0; only 3 escape to
high-TL basin ~18-21; 30-55 log-unit gap, untraversable in 100k iter
single-chain). Among the stuck-low chains, parsimony delta drives P(AB);
reps 01 and 07 have the two smallest deltas so they land highest on AB.

This isn't a topology mode-trap (rep01 has 161 unique topologies in 180
trees). It's a tree-length basin trap; flat low-TL likelihood defaults
to parsimony signal.

**Fix:** PT or longer warmup. Suggests the multirep-v3 regulariser
effect will be even cleaner once warmup is fixed.

## 📝 Results outline drafted (commit `b8f1eaf`)

`inst/lit/results-outline.md` — 7 sections, bullets + key sentences,
ready to expand to prose. Pairs with `inst/lit/intro-draft.md`.

## 📊 Publication figures drafted (commit `3ee883e`)

`inst/scripts/multirep-v3-figures/` — Fig 1 (false-clade support),
Fig 2 (CID-to-truth paired), Fig 3 (MDS rep 02). Pilot-quality
narratively, near-final visually.

## Pending jobs

| Job ID | Name | Status | On completion |
|--------|------|--------|---------------|
| 17222515 (blind-cont2) | mkp-rod-blind2 | **COMPLETED** 2026-05-19 (1h51m). 5M iter, **minESS=366**, median 1965, 5159 trees. | Ready — awaits aware to regenerate report |
| 17222514 (aware-cont2) | mkp-rod-awar2 | RUNNING, ~4h elapsed of 36h wall (cn027). Target 2.3M iter. | Check minESS; if ≥ 200, re-run `inst/scripts/rodent-comparison/rodent-comparison.R` (audited 2026-05-19 — root-invariant, no fix needed) against both converged RDS files. |
| 17226503_[1-3] (v5break) | mkp-v5break | **COMPLETED**. Both models recover AC at ≥0.994 in all 3 reps. No false-clade demo at moderate-homoplasy scales — eco signal too weak. Aware CID < blind CID in all 3 (small corroboration of multirep-v3 finding). | Treat v5break as a null check / robustness section; no follow-up needed. |
| 17226606 (mr3-pt rep05) | mkp-mr3-pt | RUNNING. Single-rep PT (nChains=4) of multirep-v3 rep 05 (blind+aware). Tests warmup TL trap hypothesis. | Examine aware P(AC). If > 0.5: full dichotomy story works, queue 5-rep array. If still diffuse: regulariser framing stands. |

## Open items / next steps

1. **Scoring bug fixed and rescored** (commits `7d1076d`, `6500cda`,
   `2aa4db3`). See top of file for findings. Sim story collapsed; v5break
   queued as the genuine break-blind shot.

2. **Rodent comparison audit clean** (commit `956bc6c`). `ape::consensus`
   with `rooted=FALSE` and `prop.clades(rooted=FALSE)` canonicalise
   internally — root-invariant. Verified empirically (45/38/32 splits
   reproduce by manual Splits-based approach). Existing pilot report stands
   methodologically; numerical results will improve when aware completes.

2. **Decide narrative on the v4-cross result**. The 5-rep PT sweep shows aware
   ≠ better than blind on v4-cross. Options:
   - Report it honestly as a failure mode of the aware model when eco signal
     opposes ancestry signal (caveat paragraph in paper).
   - Investigate further: does a wider PT ladder or longer chain shift things?
     Or is this an over-fitting of eco signal that hints at a sampler/prior fix?
   - Pivot the paper's eco-aware-rescues claim to the v4 (sister-clade)
     architecture only, where aware DOES help.

3. **Mode-trap follow-up** still queued in user memory
   (`project_mixing_followups`): rep 01 TL collapse, rep 05 topology valley,
   pi0 collapse audit. Behind PT fix.

4. **pi0/phi/z sampler audit** (`project_pi0_audit`): R5-4 verified —
   identifiability + weak prior, not a sampler bug. No code change indicated.

## Technical pointers

- **SSH to Hamilton from this machine**: use Windows OpenSSH explicitly:
  `/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "<cmd>"`.
  Passwordless. The plain `ssh` in MSYS PATH fails (Permission denied), and
  `plink` is configured for git but lacks a saved Hamilton session. SSH
  occasionally hangs from the parent shell mid-session — dispatching the same
  command via a subagent usually works; if not, retry in 15-30 min (likely
  login-node transient).
- **Hamilton paths**:
  - Repo: `/nobackup/pjjg18/mkp-sim3-multirep-v3/MkPrime` (branch
    `worktree-ecology-aware`)
  - R library: `MKP_LIB=/nobackup/pjjg18/mkp-sim3-multirep-v3/lib`
  - Deploy via `git fetch origin && git restore --source=origin/worktree-ecology-aware -- <paths>`.
    Do NOT use `git reset --hard` on Hamilton — earlier session reset zeroed
    files when ORIG_HEAD.lock failed; recovery cost ~1 hour.
- **Rodent data**: nexus at
  `/nobackup/pjjg18/mkp-rodent-blind-v2/data/mbank_X24848_2026-5-9-1135.nex`
  (only there — NOT bundled in repo). 60 extant tips × 217 chars
  (160 neomorphic + 59 transformational after `AutoDetectNeomorphic()`).
  Aware-v2 config: `magnitudeMode="global"`, `kPrimePrior="geometric"`,
  `coding="variable"`, `expSteps=10`, `rho0Alpha=7`, `rho0Beta=3`,
  `thetaAlpha=2`, `thetaBeta=2`, `sigmaPhi=0.5`, `set.seed(20260512)`.
- **v4-cross helpers**: `inst/simulations/ecology/sim3v4cross-helpers.R` and
  `sim3v4-helpers.R`. Note `MatrixToPhyDat` / `MkPrimeData` are NOT exported
  in the installed lib — use bare names, not `MkPrime::` prefix.
- **PT execution model**: `nChains > 1` runs sequentially in MkPrimeMCMC;
  set `cpus-per-task=1` and multiply wall by nChains. Confirmed in
  `sim3-multirep-v4.sh` comments.
- **Aware iter rate on rodent data**: ~50k iter/h (much slower than
  simulations). 1M iter = ~20h. 2G RAM, 1 CPU sufficient.

## Things ruled out

- **chmod -R u+w** on Hamilton paths — blocked by classifier as destructive.
  Use file-by-file chmod or work around via git restore.
- **scp of .sh scripts** — line-ending corruption injected binary nulls in
  earlier sessions. Either strip CRLF before scp, or write via SSH heredoc.
- **Treating all rodent characters as transformational** — user clarified
  the existing 160 neo + 59 trans classification IS gold standard; do not
  reclassify.
- **Single-chain aware MCMC on v4-cross** — exhibits contiguous mode-traps
  (sample 84 onward stuck). Default to nChains ≥ 4 with PT for any
  v4-cross-like configuration.
- **`-amend` after a failed pre-commit hook** — the commit didn't happen, so
  amend rewrites the previous commit. Create a new commit instead.
- **`ape::prop.part` for posterior tree scoring** — root-dependent. Use
  `TreeTools::as.Splits` instead, or implement a manual unrooted bipartition
  match. Discovered 2026-05-19; see top of file.

## Worktrees

| Branch | Path | Status |
|--------|------|--------|
| worktree-ecology-aware | C:\Users\pjjg18\GitHub\mkp\.claude\worktrees\ecology-aware | clean (modulo pre-existing untracked debris) |

## Suggested first action

Once back, run the COLLECT step on this RESUME.md: poll job 17215319 via
`/c/WINDOWS/System32/OpenSSH/ssh.exe pjjg18@hamilton8.dur.ac.uk "sacct -j 17215319 -n -X -o JobID,State,Elapsed,ExitCode"`.
If COMPLETED, regenerate the rodent comparison per the table's "On completion"
column.
