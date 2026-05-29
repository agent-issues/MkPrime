# Per-move detailed-balance localiser for the EBE topology moves

## What this tests
A **localiser** companion to `ebe-rooted-posterior.R`. The end-to-end posterior
test is the decisive verdict on whether the topology moves sample the correct
rooted posterior; this script does NOT re-decide that. Its job is: IF the
end-to-end test flags a move biased, tell us WHICH factor of that move's
`logHastings` is wrong. It decomposes the proposal log-density ratio into the
two factors a fixed-dimension MH move must get right and checks each exactly.

## Theoretical basis
For an MH topology move to leave `pi propto L_EBE . Gamma(TL) . Dir(relBr)`
invariant, the reported `logHastings` must equal the TRUE proposal log-density
ratio `log[q(R'->R)/q(R->R')]` of the kernel actually executed (this is the
*only* condition; by MH construction the single-step probability flux is
identically balanced for ANY `logHastings`, so a wrong ratio can ONLY be
detected as a wrong STATIONARY distribution — which is `ebe-rooted-posterior.R`'s
job). That ratio factorises as:

- **(i) discrete choice-multiplicity ratio** `log[m(R'->R)/m(R->R')]`, where
  `m = (#prune edges) x (#regraft candidates) x (#subtree edges, TBR only)` is
  the inverse product of the move's uniform pick probabilities; and
- **(ii) continuous branch-length Jacobian** `log|J|` of the realised map from
  the free U(0,1) variates (tau for SPR; sigma,tau for TBR Phase-B) to the new
  branch lengths.

The code (`src/proposals.cpp:172`, `src/tree_moves.cpp:504`) sets `logHastings`
to ONLY factor (ii) (SPR: `log(lRegraft)-log(lMerge)`; TBR adds
`log(lSubEdge)-log(lMergeSub)`), i.e. it ASSERTS factor (i) is identically 1
(the TBR comment at `tree_moves.cpp:500-503` states this explicitly for the
candidate count; SPR omits it silently).

- **Part A (exact, no Monte-Carlo)** enumerates `m(R->R')` vs `m(R'->R)` over
  many states/choices and confirms they are equal (factor (i) == 1). A nonzero
  asymmetry is a smoking gun.
- **Part B (re-derivation-free Jacobian)** reads the realised edge lengths of a
  move and reconstructs the Jacobian purely from the MULTISET of changed lengths
  (no assumption about WHICH edges, so it cannot share the code's indexing blind
  spot), comparing `log` of it to the reported `logHastings`. A length-preserving
  merge/split move has, per merge-split pair, one DISAPPEARED "split-source"
  length = sum of two APPEARED lengths and one APPEARED "merge-target" = sum of
  two DISAPPEARED; the Jacobian is `prod(split-sources)/prod(merge-targets)`.
  This covers BOTH moves: SPR has ONE such pair (`lRegraft/lMerge`); TBR, when
  its **Phase-B sigma re-root** fires, has TWO (adding `lSubEdge/lMergeSub`),
  giving `exp(logHastings) == (lRegraft*lSubEdge)/(lMerge*lMergeSub)`. The
  split-count (1 vs 2) classifies SPR-equiv vs Phase-B moves and so also
  measures how often Phase-B fires. NOTE: finite differences are NOT usable --
  the (lengths)->(lengths) map is singular (the DOF lives in tau/sigma) -- which
  is why the multiset sum-matching is the right tool; GENERIC distinct branch
  lengths make the matching unambiguous.

Reference: detailed balance for fixed-dimension Metropolis-Hastings (Hastings
1970; Green 1995 for the dimension-matched branch-length Jacobian convention).

## Pass criterion
- **Part A SPR (discrete count):** `m(R->R') == m(R'->R)` EXACTLY for every
  enumerated forward/reverse pair (zero asymmetric pairs). Exact integer
  equality — no tolerance. (Observed: 0 asymmetric over 10k–430k pairs.)
- **Part B SPR (Jacobian):** total tree length preserved to `< 1e-9` relative,
  and `|logHastings - log(lRegraft/lMerge)| < 1e-8` for every decoded move.
  (Observed ~1e-15.)
- **Part B TBR (Jacobian, the crux):** total length preserved; Phase-B MUST
  actually fire (`nPhaseB > 0` — else the run is underpowered, reported as
  FAIL not PASS); and `|logHastings - log((lRegraft*lSubEdge)/(lMerge*
  lMergeSub))| < 1e-8` for every decoded Phase-B move (and the SPR-equiv
  branch `< 1e-8` too). (Observed ~9e-16 over hundreds of Phase-B moves.)

**Rigor status of the two TBR detailed-balance factors (be precise):**
- TBR **continuous Jacobian** (incl. Phase-B sigma): **VERIFIED-MEASURED** to
  machine precision (Part B above; ~9e-16 over ~14k fired Phase-B moves).
- TBR **discrete count symmetry** (`nPrune . nSubEdge . nCand` forward vs
  reverse): **ARGUED, not directly measured here.** The argument: `nPrune` is
  invariant (always `nEdge-2`); `nSubEdge` is invariant because re-rooting
  preserves the pruned subtree's taxon set and hence its internal edge count;
  `nCand = nEdge - |desc-edges(v)| - 2` is invariant because re-rooting
  preserves `|desc(v)|`; and the prune+regraft selection reduces to the SPR
  argument MEASURED exactly in Part A (433k pairs, 0 asymmetric). We do NOT wire
  in an apply-and-recount measurement for TBR because, unlike SPR, a single TBR
  move shifts the local sibling/candidate context of SEVERAL subtrees (the moved
  one plus collateral subtrees near the suppressed node and the regraft point),
  and an R-side heuristic for "which subtree was the one the C++ moved" cannot
  cleanly isolate it without re-implementing TBR's exact node tracking — a
  buggy recount would risk a FALSE asymmetry (a worse error than an honest
  gap). A genuine count asymmetry on the moved subtree would bias the
  STATIONARY distribution and is therefore backstopped by
  `ebe-rooted-posterior.R --full` (6 tips, ~1.8M Phase-B events), whose
  TVD-vs-ESS discriminator + TBR positive control measure counts and Jacobian
  jointly. The smoke run already showed correct-TBR TVD WITHIN its null band
  (no plateau) at ESS ~18k, consistent with no net count asymmetry.

## How to run
- **Small-N (<= 60 s):** `Rscript ebe-move-detailed-balance.R --quick`
  (200 states for Part A; 3000 SPR Jacobian samples). Observed runtime ~1 s —
  this script is proposal-only (no likelihood, no MCMC), so it is fast and the
  quick/full distinction is only sample count.
- **Full-scale:** `Rscript ebe-move-detailed-balance.R --full` (5000 states,
  1e5 SPR Jacobian samples; tips 5-9). Runs in seconds-to-minutes on a laptop;
  no HPC needed.

## Output interpretation
`ebe-move-detailed-balance-results/verdict.txt`: three PART lines (SPR
discrete-count, TBR subtree-conservation, SPR Jacobian) each PASS/FAIL, then
`OVERALL`. `summary.rds` holds the counts (pairs examined, asymmetric pairs,
max |asymmetry|, decoded-move count, max Jacobian residual).

## What a failure would mean
- **Part A SPR FAIL** (asymmetric `m`): SPR drops a candidate/prune-count term
  from `logHastings` → biased rooted posterior. Direct cause.
- **Part A TBR FAIL** (invariant violated): the TBR move corrupts the tree
  (lost/duplicated taxon or wrong edge count) → a structural move bug upstream
  of the Hastings ratio.
- **Part B SPR FAIL** (Jacobian mismatch / length not preserved): the
  branch-length redistribution Jacobian is wrong → biased over branch lengths
  and rooting.

## Coverage
This script analytically clears, to machine precision: the **SPR** Jacobian and
discrete count, AND **TBR's Phase-B sigma re-root Jacobian** (the task's crux —
`log(lSubEdge)-log(lMergeSub)`, the distinctive sigma/tau merge-split), with
Phase-B confirmed to actually fire. It is the PRIMARY, exact verdict on the
move Jacobians; the end-to-end `ebe-rooted-posterior.R` additionally confirms
the *stationary distribution* is correct (catching any residual discrete-count
error, which a Jacobian check cannot see) with SPR and TBR positive controls.
**pSPR** (parsimony-guided, data-dependent proposal density) is out of scope
here — its density ratio is not a simple count×Jacobian — and is covered
end-to-end by `ebe-rooted-posterior.R` (the posterior-match test is agnostic to
how the proposal density is built); until wired in there it is honestly marked
UNVERIFIED.
