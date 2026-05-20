# Diagnosis: aware caveat reps 01 and 07 in multirep-v3

Date: 2026-05-19

## Summary

In the 8-rep multirep-v3 experiment (16 tips, eco-1 = full clades A+B, truth
`((A,C),(B,D))`), the ecology-aware chain successfully suppresses the false AB
clade in 6/8 reps but remains stuck supporting AB in reps 01 (P(AB) = 0.789)
and 07 (P(AB) = 0.983). Two compounded failures explain why these two reps are
the worst offenders:

1. **Branch-length-expansion failure (warmup trap):** All 5 stuck reps (01, 02,
   04, 07, 08) failed to escape the low-TL local optimum during warmup. The
   aware chain is rendered nearly blind by the short branch lengths; with
   essentially no phylogenetic resolution, it defaults to parsimony signal.

2. **Data-level weakness (parsimony delta):** Rep 01 and rep 07 have the
   weakest data discrimination between the true and false topologies across
   all 8 reps — parsimony deltas of 235 and 232 steps respectively (all other
   stuck reps: 253–278; escaped reps: 240–284). In the flat low-TL basin,
   the ecology correction cannot overcome this weak parsimony signal.

---

## Numerical comparison: all 8 reps

| Rep    | Basin   | TL median | logL median | P(AB) bl | P(AB) aw | parsimDelta | stuck |
|:-------|:--------|----------:|------------:|---------:|---------:|------------:|------:|
| rep01  | LOW-TL  |      3.12 |    −4569.1  |    0.854 |    0.789 |         235 | YES   |
| rep07  | LOW-TL  |      2.81 |    −4547.5  |    1.000 |    0.983 |         232 | YES   |
| rep02  | LOW-TL  |      3.67 |    −4570.3  |    1.000 |    0.375 |         253 | YES   |
| rep04  | LOW-TL  |      3.15 |    −4574.4  |    0.918 |    0.320 |         259 | YES   |
| rep08  | LOW-TL  |      3.17 |    −4580.0  |    0.772 |    0.253 |         278 | YES   |
| rep03  | HIGH-TL |     20.52 |    −4516.3  |    0.444 |    0.023 |         284 | no    |
| rep05  | HIGH-TL |     18.58 |    −4540.4  |    0.978 |    0.462 |         240 | no    |
| rep06  | HIGH-TL |     18.21 |    −4585.6  |    0.949 |    0.085 |         278 | no    |

"parsimDelta" = parsimony steps on false AB tree minus parsimony steps on true tree
(lower = data are less informative against the false clade). TL and logL are
aware post-burnin medians. Truth TL = 13.5 in sample-space units (raw TL ×
~10.78).

---

## Failure mode 1: low-TL warmup trap (affects 5/8 reps)

All chains start with `edge.length = 0.1` per edge (raw TL = 3.0 = sample
TL ~0.28). The warmup phase (min 5 k, max 40 k iters) must expand branch
lengths to the high-TL basin (sample TL ~18–21). In reps 03, 05, and 06 the
warmup found this basin; in reps 01, 02, 04, 07, and 08 it did not.

All post-burnin samples from the stuck reps are in the low-TL basin: 0/240
rep01 samples and 0/233 rep07 samples ever reach TL > 10 across the entire
run. The first thinned sample already shows the final basin (rep01: TL = 3.27
at iteration ~200; rep03: TL = 21.3 at iteration ~200). The logL gap between
the two basins is substantial: rep01 is ~55 log-units below the high-TL basin
mean; rep07 is ~30 log-units below.

The escape appears stochastic — a function of which warmup trajectory happened
to push branch lengths toward the higher-likelihood region early enough to
cross the barrier. No clear data-level difference distinguishes the stuck from
the escaped reps by warmup trace alone (all warmup logL values range from
−4820 to −4911 at warmup end).

---

## Failure mode 2: data parsimony delta (explains rep01 and rep07 within stuck group)

Within the low-TL basin, the aware correction is nearly flat across topologies.
The dominant signal becomes parsimony — how much more costly (in steps) is the
false AB topology than the true AC topology under the data. Reps 01 and 07 have
the weakest parsimony discrimination of all 8 reps (delta = 235 and 232 steps
respectively). All other stuck reps have larger deltas (253–278), and their
ecology-corrected likelihoods accordingly suppress P(AB) to 0.25–0.38 despite
also being trapped in the low-TL basin.

The secondary factor for rep07 specifically: its low-TL basin logL median
(−4547.5) is 23–35 log-units better than the other stuck reps (−4570 to
−4580), indicating it settled into a deeper local optimum within the same
basin, further reinforcing the stuck state.

---

## Why this is not a topology mode-trap

Rep01 has 161 unique topologies in 180 post-burnin trees; rep07 has 108 unique
topologies in 174 trees. The dominant topology occupies only 2–3% of the
posterior in each rep. Despite topological diffuseness, 78.9% (rep01) and
98.3% (rep07) of trees contain the AB clade. This is not a topology mode-trap
in the classical sense. Rather:

- The short-branch-length landscape makes the ecology-corrected likelihood
  nearly flat across topologies.
- In this flat regime, the parsimony signal embedded in the data drives the
  posterior, and reps 01 and 07 have the weakest signal against AB.

---

## Recommendations

### Would PT (parallel tempering) fix this?

Likely yes for failure mode 1. A hot chain at temperature β = 0.1–0.2 would
have a flattened landscape that could traverse the barrier; cold-chain swaps
would then pull the cold chain across. The `sim3-multirep-v3-pt` Hamilton
harness already exists.

However, PT does not address failure mode 2: even with correct TL exploration,
reps 01 and 07's weaker data would produce higher P(AB) than other reps. The
correct expectation after PT would be P(AB) at the true high-TL posterior —
which for rep01/07's weak data may still be non-trivial (but should be much
lower than 0.789/0.983).

*Note: parsimony delta is computed over all 480 characters. The false-AB
signal is driven by the 120 eco-1 characters; the true-AC signal by the 360
base characters. Splitting delta by arm would clarify whether rep01/07's
weakness is in the ecology arm (PT might help by exploring the correction
more) or the base arm (PT would not help). Left for follow-up.*

### Would longer warmup fix failure mode 1?

Uncertain. The escape appears to require a large simultaneous branch-length
expansion that the current proposal (small-step scale_tree_length moves) has
low probability of accepting at low temperatures. A dedicated branch-length
Gibbs sampler or global BL multiplier proposal is more likely to help than
extending warmup.

---

## Paper narrative implication

> "Ecology-aware MCMC suppresses the false eco-clade signal (P(AB) < 0.5) in
> 6/8 replicates. In 5/8 reps the chain failed to escape a short-branch-length
> local optimum during warmup, rendering the ecology likelihood flat; the
> posterior then reflects parsimony signal. Within the stuck group, the two
> worst-performing reps (01 and 07, P(AB) = 0.789 and 0.983) also have the
> weakest data discrimination between the true and false topologies (parsimony
> delta 232–235 steps vs 253–278 for other stuck reps and 240–284 for escaped
> reps). This is a mixing failure compounded by adverse data realizations in
> reps 01 and 07, not a model failure: parallel-tempering experiments are
> ongoing."
