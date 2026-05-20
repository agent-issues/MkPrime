# Parsimony trap test — does MP recover the false eco-clade?

**Date:** 2026-05-20
**Driver:** `inst/ecology/scripts/parsimony-trap-test/parsimony-trap.R`
**Results:** `parsimony-trap-results.csv`, `parsimony-trap-summary.csv`

## Hypothesis

If `TreeSearch::MaximizeParsimony` recovers the FALSE eco-clade rather than truth,
ML can plausibly be fooled too. If MP robustly recovers truth even when eco
signal is strong, no "trap" exists for likelihood either. This bounds the
design space for breaking blind at realistic TL.

## Per-dataset trap fractions

False-clade-found fraction (any MPT in the search hit set contains the false eco bipartition):

| Dataset | TL | nChar | EW (cv=Inf) | IW k=3 | IW k=6 | IW k=12 | TrueSister found (any cv) |
|---|---|---|---|---|---|---|---|
| v6-realistic (n=3) | 1.16 | ~194-300 | 0/3 | 0/3 | 0/3 | 0/3 | 3/3 |
| multirep-v3 (n=8)  | 13.5 | 480     | 7/8 | 8/8 | 8/8 | 7/8 | 0/8 |
| v5break (n=3)      | 1.30 | 80      | 0/3 | 0/3 | 0/3 | 0/3 | 3/3 |

Mean (EW Fitch-step) gap between truth tree and best MPT:

| Dataset | mean truthLen - MPTlen (steps) |
|---|---|
| v6-realistic | 0.0 (truth IS the MP tree) |
| multirep-v3  | ~28 - 30 (truth is ~30 steps longer than the false-eco MPT) |
| v5break      | ~7 - 8 |

## Diagnostic — does MP behaviour predict the MCMC blind chain we already saw?

Yes, tightly:

* **multirep-v3 (TL=13.5, nChar=480).** MP recovers the false eco-clade in 7-8/8 reps under every weighting. This matches the rescored blind chain (P(AB) ~0.85 in rep01, similarly high in others). MP and ML/MCMC see the same data; both are equally fooled by saturated eco-convergent signal. **This is the proof-of-concept that ML's blind failure at saturation is a fundamental information-content problem, not an MCMC artefact.**
* **v6-realistic (TL=1.16, nChar=~200-300).** MP recovers the truth tree exactly under all weightings (truth IS the unique MP solution, 0 step gap). This matches blind P(AC)=1.00 in rep01 corrected scoring. At realistic TL the eco signal does not generate enough homoplasy to overwhelm the clade-stem synapomorphies — exactly the regime where v7-induce was needed to force a break-blind by lifting phi/pi0.
* **v5break (TL=1.30, nChar=80).** Despite the small matrix, MP still recovers the true sister under all weightings (false-cherry hit fraction = 0). The truth-tree gap is only ~7 steps (about 9% of MPT length), so the false signal is present but not dominant. This explains why v5break's blind chains in the rescored runs did *not* break in the way the v5break design predicted: parsimony can already tell truth from the false inner cherry on this matrix, and so can ML.

The clean rank ordering `multirep-v3 >> v5break ~ v6-realistic` in trap strength mirrors the rank ordering of blind chain failure we observed in the corrected MCMC scoring. **MP is a fast, deterministic proxy for whether the data contain a likelihood trap.**

## Does implied weighting amplify the trap?

Not in any of these datasets. IW (k=3, 6, 12) gives the same false/true clade verdict as EW on 13/14 datasets x rep cells. The single discrepancy is multirep-v3 rep06: EW and IW k=12 miss the false eco clade (the MPT under those weightings includes a different non-truth bipartition), while IW k=3 and k=6 find it. Substantively all four weightings fail to recover truth — no rep has IW flipping a recovered-truth verdict to a false-eco verdict that wasn't already present.

Interpretation: at high TL the false eco clade is supported by so many step-saving characters that even down-weighting homoplasy cannot disentangle it. At low TL the true clade signal is clean enough that down-weighting doesn't pull MP into the eco basin. Homoplasy-downweighting does not unlock a "trap zone" that EW would have hidden.

## Implications for v8 design

To find a genuine likelihood trap at realistic TL we need a regime where MP picks up the false eco clade at TL well below 13.5 — and certainly below the v5break range (TL=1.30). Candidate dials, holding nTip=16 and the convergent 8-clade topology:

* **Push phi higher and pi0 lower.** Multirep-v3 used phi=4, pi0=0.95-style; v6 used phi=8, pi0=0.45. The trap appears between these — but only at saturating TL. v8 should target phi >= 8, pi0 <= 0.3 with TL ~ 3-6 (intermediate).
* **More characters.** At TL=1.16 with nChar=200, truth wins by zero steps. At nChar=480 with same TL it would likely still win (variance shrinks). nChar=2000 is required to make the false-eco basin sample-size-deep at low TL.
* **Stem-length contrast.** v5break shrank stemBrClade to 0.03; this was the strongest "break" parameter. Going further (stemBrClade=0.015) plus higher stemBrEco (0.20) should put MP in the trap regime at TL ~ 2 with moderate nChar.
* **Watch IW.** Don't expect IW to flip the verdict; if MP under EW recovers truth, the data are not a trap for ML.

Recommended v8 cell: (TL ~ 2.5, nChar = 600, phi = 10, pi0 = 0.30, stemBrClade = 0.015, stemBrEco = 0.20). Pre-screen with MP-only EW; only run MCMC on cells where MP recovers the false eco clade.

## Commit

To be reported after commit.

## Files

* `inst/ecology/scripts/parsimony-trap-test/parsimony-trap.R` — driver
* `inst/ecology/scripts/parsimony-trap-test/parsimony-trap-results.csv` — per (dataset, rep, weighting) row
* `inst/ecology/scripts/parsimony-trap-test/parsimony-trap-summary.csv` — per (dataset, weighting) aggregate
* `inst/ecology/scripts/parsimony-trap-test/parsimony-trap.log` — run log
* `inst/ecology/scripts/parsimony-trap-test/v5break-data/v5break-rep0{1,2,3}.rds` — local copies of Hamilton blind-result RDS (read-only)
