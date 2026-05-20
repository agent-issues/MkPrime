# Rodent aware-v2 posterior diagnostics
Source: rodent-aware-v2-result.rds  (actual_iter=1000000; nSamples=886; retained post-burnin scalar=443; retained z=591)
Model prior: rho0=Beta(7,3), theta=Beta(2,2), sigmaPhi=0.5, relabel=TRUE

## phi (ecological rate multiplier)
median = 4.028   90% CI = [3.097, 5.781]   95% CI = [2.980, 6.294]

## pi0 (slab probability)
median = 0.430   95% CI = [0.246, 0.577]

## Characters flagged P(z!=0) > 0.5 per ecology
  ecology    eco_name thr30 thr50 thr70 thr90
1    eco1    arboreal   197   131    56    16
2    eco2 semiaquatic   217   213     4     0
3    eco3   fossorial   216   149    28     3

## Cross-ecology overlap
   eco_a eco_b      name_a      name_b n_flag_a n_flag_b n_both   jaccard
1   eco1  eco2    arboreal semiaquatic      131      213    127 0.5852535
11  eco1  eco3    arboreal   fossorial      131      149     93 0.4973262
2   eco2  eco3 semiaquatic   fossorial      213      149    145 0.6682028
     cor_pearson
1  -0.0477404239
11 -0.0122009501
2  -0.0004733397

Mean inter-ecology cor(P(z!=0)) = -0.020
Mean Jaccard overlap            = 0.584
Median characters flagged/eco   = 149

## Ecological-gravity verdict
phi posterior median 4.03  (95% CI [2.98, 6.29])
This is the empirical signal strength on real rodent morphology.
Median of 149 characters per ecology show P(z!=0)>0.5 (semiaquatic
saturates at 213/217 because the 5 semiaquatic tips give weak per-character
evidence, so the slab eats up almost everything).

## Sign / state composition
For arboreal (eco1), the top-10 characters are essentially all P(z=2)
(asymmetric / "discouraged"). For fossorial (eco3), again dominated by z=2.
For semiaquatic (eco2) the top characters are a mix: a number are P(z=1)
(symmetric / encouraged) and a number P(z=2). This is consistent with
semiaquatic having too few tips for the model to commit to a sign.

## Top characters with biological labels
arboreal:   axial body suture, medial tarsal sesamoid, third trochanter,
            infraorbital canal position, greater tuberosity of humerus
fossorial:  sphenopalatine foramen, third trochanter, nasolacrimal canal,
            radial fossa, medial tarsal sesamoid, mastoid foramen, orbit
            orientation
semiaquatic (weaker, mixed sign): premaxilla size, greater tuberosity,
            calcaneoastragalar facet, jugal ventral process
The biological signal is plausible: arboreal/fossorial signatures concentrate
on locomotor / postcranial characters (humerus, femur trochanters, tarsals)
plus a few craniofacial features (infraorbital canal, sphenopalatine foramen
in fossorial diggers).

## Implication for v8 sim design
- Use phi ~ 4 as the realistic "ecological gravity" centre (NOT 6 or 8).
  v8 should sweep phi in {2, 4, 6} centred on 4, with 4 as the canonical
  realistic value. phi=8 would overstate the real-world signal.
- pi0 ~ 0.43: of the chars NOT flagged as eco-affected, about 57% are
  classified slab (z=0). Equivalently, ~43% of (char, eco) cells show
  posterior eco-affected mass. Use a similar pi0 truth in v8 if simulating
  to match real-world sparsity, but be aware the rho0=Beta(7,3) prior is
  fairly informative here (prior median 0.7 pulls toward 0.7 — the posterior
  has shifted markedly toward more slab than the prior, indicating real
  signal).
- Number of eco-affected characters per ecology in truth should be ~125-150
  (matching the P(z!=0)>0.5 counts) — NOT 20 or 30 as in earlier sims.
- Inter-ecology overlap is high (Jaccard ~0.58, mean Pearson ~0): different
  ecologies recruit largely overlapping but uncorrelated character sets.
  v8 should let truth-z be drawn independently per ecology (not shared),
  which gives moderate overlap by chance.
- Sign: real data is strongly dominated by z=2 (asymmetric); v8 truth
  should not be 50/50 z=1/z=2, but more like 70/30 in favour of z=2.

