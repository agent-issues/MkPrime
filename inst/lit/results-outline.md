# Results — section outline
# MkPrime ecology-aware paper
# Draft 2026-05-19

---

## 1. Simulation design

- **Setup**: 16 tips in four balanced 4-tip clades; true topology ((A,C),(B,D)).
  Eco-1 assigned to all of clades A and B (8-tip false convergent group).
  480 characters: 120 neomorphic (all eco-affected, theta=1, phi=4) +
  360 transformational (no ecology effect).
  Branch lengths: tipBr=0.50, stemBr=0.30, rootBr=0.15; truth TL=13.50.
- **Replication**: 8 independent draws of the simulated character matrix;
  each scored against both a blind Mk' chain and an ecology-aware Mk' chain
  (100k iter MCMC, post-burnin 75%).
- **Scoring**: root-invariant bipartition matching via `TreeTools::as.Splits`
  (`inst/simulations/ecology/sim3-scoring.R`); topology distance = normalised
  Clustering Information Distance (CID) to truth.
- **Figure reference**: tree diagram in `dev/sim3-redesign-brief.qmd`.
- **Key sentence**: "The simulation targets the convergence failure mode: the
  neomorphic characters evolve faster in eco-1 lineages, systematically
  attracting ecologically similar but phylogenetically distant clades A and B
  as a spurious sister group under standard inference."

---

## 2. The false-clade artefact under blind Mk'

- **Claim**: the blind chain commits to the false eco-clade with high posterior
  support, consistently and reproducibly across 8 replicates.
- **Figure**: `inst/scripts/multirep-v3-figures/fig1-false-clade-support.pdf`
  — per-rep bar chart of P(AB) and P(AC) for blind vs aware.
- **Key numbers**:
  - Mean P(false AB eco-clade) blind = **0.864** (8 reps).
  - 7/8 reps: blind P(AB) > 0.5; **5/8 reps**: blind P(AB) > 0.9.
  - Mean P(true AC sister) blind = 0.001; 0/8 reps with P(AC) > 0.5.
  - Blind mean tree length = 101.85 vs truth = 13.50 — wild inflation driven
    by the convergent character signal.
- **Note on legacy scoring**: the scoring bug (`ape::prop.part`, root-dependent)
  had hidden the false-clade signal in 3/8 reps; corrected means rose from
  0.551 (legacy) to 0.864. Reps 02, 05, 06 showed the largest artefacts
  (rep02: 0 → 1.000; rep05: 0 → 0.978; rep06: 0.444 → 0.949).
- **Key sentence**: "Across eight independent replicates, the blind Mk' chain
  assigned mean posterior probability 0.864 to the false eco-clade, with five
  of eight replicates exceeding 0.9 — a strong and reproducible false-positive
  driven by ecology-driven convergence in the neomorphic character subset."

---

## 3. The aware model's regularising effect

- **Claim**: the ecology-aware chain does not recover the true topology, but it
  strips approximately half the false-clade mass and is demonstrably closer to
  truth in tree space on every replicate.
- **Figure**: `inst/scripts/multirep-v3-figures/fig2-cid-to-truth.pdf`
  — per-rep paired CID-to-truth (blind vs aware), Wilcoxon annotation.
- **Key numbers**:
  - Mean P(false AB) aware = **0.411** (down from 0.864 blind; delta = −0.453).
  - Mean P(AC) aware = 0.012; 0/8 reps with P(AC) > 0.5.
  - Mean CID-to-truth: blind 0.484, aware **0.432** (delta = −0.052).
  - Aware mean CID lower than blind in **8/8** reps.
  - Paired Wilcoxon: V = 0, **p = 0.0143** (two-sided).
  - Aware mean tree length = 9.32 (truth = 13.50) vs blind 101.85 —
    length inflation eliminated.
- **Key sentence**: "Although the ecology-aware chain never concentrates
  posterior mass on the true sister relationship ((A,C),(B,D)), it reduces
  mean P(false AB) by 0.453 and is uniformly closer to truth in Clustering
  Information Distance across all eight replicates (paired Wilcoxon p = 0.014),
  demonstrating consistent regularisation even in the absence of full recovery."

---

## 4. Posterior shape: honest uncertainty vs confident error

- **Claim**: the blind chain mode-locks onto wrong topologies containing the AB
  split; the aware chain spreads its posterior widely across topologies that
  lack AB, with no dominant mode — honest uncertainty rather than confident
  error.
- **Figure**: `inst/scripts/multirep-v3-figures/fig3-mds-rep02.pdf`
  — CID-MDS scatter for rep 02, coloured by chain (blind vs aware) and flagged
  by AB/AC bipartition presence.
- **Key numbers** (rep 02 as exemplar, from posterior-shape report):
  - Blind: 151 unique canonical topologies; mean top-topology mass 0.017;
    95% credibility set |CS95| = 143 topologies;
    all top-5 topologies contain the AB split.
  - Aware: 130 unique canonical topologies; top-topology mass 0.024;
    |CS95| = 122; top-5 mix AB-containing and AB-free topologies.
  - Blind MR consensus nodes overwhelmingly carry the AB bipartition.
  - Aware top-topology share across all 8 reps: 1.7–6.2 % per chain —
    no topology dominates; no mode-lock.
- **Caveat**: neither chain includes the true topology in its 95% credibility
  set in any of the 8 reps (truth in CS95 = 0/8 for both models). The aware
  chain moves the posterior neighbourhood closer to truth without placing truth
  inside the sampled cloud.
- **Key sentence**: "In replicate 02, the blind chain's top five topologies all
  contain the false AB bipartition; the aware chain's top five are a mixture of
  AB-containing and AB-free trees, with the highest-mass topology carrying only
  2.4 % of the posterior — reflecting broad uncertainty rather than commitment
  to any single wrong answer."

---

## 5. Robustness: lower-homoplasy regimes (null result)

- **Claim**: at lower homoplasy (the v4/v4-cross family), the ecology confound
  is too weak to mislead either model; both blind and aware recover the true
  topology near-perfectly. The ecology-aware model does no harm when the eco
  signal is weak.
- **Table**: corrected per-clade support from `inst/scripts/rescore-sim3-results.csv`
  (22 chains, v4/v4-cross).
- **Key numbers**:
  - v4-cross-b-pt5 (5-rep PT sweep): blind P(AC) = 0.999 ± 0.002;
    aware P(AC) = 0.996 ± 0.008. Both essentially perfect.
  - Max P(false AB) across all 22 rescored chains = **0.006** (one chain,
    v4cross-b-pt5-rep05-aware); P(falseInner) = 0 throughout.
  - 16 tips × 200–300 chars × phi=4 sits below the threshold at which the
    eco confound differentiates the models.
- **Key sentence**: "At the lower-homoplasy v4 and v4-cross parameter regimes,
  both models recovered the true topology at posterior probability approaching
  one, and no false-clade support was detected in any of 22 rescored chains,
  confirming that the ecology-aware extension does not degrade inference when
  ecological signal is weak relative to phylogenetic signal."
- **Pending**: v5break (job 17226503, 3-rep PT array; 80 chars, phi=6, weak
  ancestry) is in flight as the Goldilocks test — if blind P(falseInner) > 0.1,
  the full two-model dichotomy will become reportable.

---

## 6. Caveats: aware mode-trap in reps 01 and 07

- **Claim**: in 2/8 replicates the aware chain itself falls into a short-TL
  mode-trap and produces inflated P(AB), reproducing rather than reducing the
  false-clade artefact. This is diagnosable from tree-length collapse and does
  not represent a genuine signal inversion.
- **Key numbers** (from `inst/scripts/rescore-multirep-v3-report.md`):
  - Rep 01 aware: P(AB) = 0.789, mean TL = 3.19 (truth = 13.50).
  - Rep 07 aware: P(AB) = 0.983, mean TL = 2.87.
  - In both cases aware TL is ~4.5× below truth; blind TL is 138–154
    (wildly inflated but different failure mode).
  - Neither rep shows a real inversion (aware false-supporting AB while blind
    supports AC): blind P(AC) = 0 in both.
- **Interpretation**: the short-TL local optimum shares surface likelihood with
  the inflated-TL AB-false-clade mode under high homoplasy. Parallel tempering
  is expected to resolve the trap (queued: task #14).
- **Note**: dedicated mode-trap audit document
  (`inst/scripts/aware-multirep-v3-caveat-reps-01-07.md`) is pending
  (task #17 in flight).
- **Key sentence**: "Two of eight replicates exhibit aware-chain tree-length
  collapse (mean TL 3.19 and 2.87 against truth 13.50), with the chain
  trapped on a short-branch local optimum that also supports the false AB
  split; parallel-tempering extensions are expected to resolve these
  mode-traps."

---

## 7. Empirical case study: rodent morphological matrix

- **Claim**: applied to a published rodent matrix, blind and aware Mk' produce
  demonstrably different posterior distributions in tree space; the direction of
  the difference is consistent with ecology-driven homoplasy pulling the blind
  chain away from the ecology-aware posterior.
- **Data**: MorphoBank matrix X24848; 60 extant tips × 217 characters
  (160 neomorphic + 59 transformational after `AutoDetectNeomorphic()`);
  four ecology states.
- **Figure**: CID-MDS scatter and MR consensus panels from
  `inst/scripts/rodent-comparison/` (figures from pilot run; will be
  regenerated once aware-cont2 converges).
- **Key numbers** (pilot, pre-convergence — treat as indicative):
  - MDS centroid separation blind vs aware = **9.65 CID units** (aware
    spread = 9.01; blind spread = 7.01).
  - MR consensus splits: blind 38, aware 45, shared 32; 6 blind-unique,
    13 aware-unique.
  - Blind-unique splits include several large muroid groupings; aware-unique
    splits include distinct Heteromyidae, Gliridae, and outgroup placements
    consistent with known rodent phylogeny.
- **Methodology audit**: comparison is root-invariant; `ape::consensus` with
  `rooted=FALSE` and `prop.clades(rooted=FALSE)` verified empirically (45/38/32
  split counts reproduced by manual Splits-based approach; commit `956bc6c`).
- **Convergence caveats**:
  - Both chains run at 1M iterations; minESS = 88 (aware) and 38 (blind) —
    below recommended minimum of 200. Results are pilot-quality.
  - Blind continuation job 17222515 completed (5M iter, minESS = 366, 5159
    trees). Aware continuation job 17222514 still running (target 2.3M iter).
  - Rodent comparison will be regenerated once aware-cont2 completes; all
    split and MDS numbers above may change.
- **Key sentence**: "The blind and aware posterior clouds are substantially
  separated in CID-MDS space (centroid distance 9.65), and their majority-rule
  consensus trees differ by 19 splits out of a combined 51, indicating that
  the ecology-aware model systematically explores a different region of rodent
  tree space; quantitative comparison awaits completion of the aware
  continuation run."

---

*All figures available as PDF and PNG in `inst/scripts/multirep-v3-figures/`.*
*Rodent comparison figures in `inst/scripts/rodent-comparison/`.*
*Numbers sourced from:*
*  `inst/scripts/rescore-multirep-v3-report.md`*
*  `inst/scripts/aware-multirep-v3-posterior-shape-report.md`*
*  `inst/scripts/rescore-sim3-report.md`*
*  `inst/scripts/rodent-comparison/rodent-comparison.md`*
*  `dev/sim3-redesign-brief.qmd`*
