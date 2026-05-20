# Audit: rep01 TL collapse + rep05 topology valley vs root-invariant scoring

**Date:** 2026-05-19
**Branch:** `worktree-ecology-aware`
**Scope:** Re-audit the two queued mixing pathologies (`project_mixing_followups`) against root-invariant bipartition + canonical-splits topology hashing, after the discovery that legacy `ape::prop.part` scoring and the chain-emitted `topo_hash` are root-dependent.

## Data audited

- `inst/ecology/simulations/multirep-v3-results/rep01/{aware,blind}-result.rds`
- `inst/ecology/simulations/multirep-v3-results/rep05/{aware,blind}-result.rds`

Sample matrices are empty in the blind RDS files (streaming-only); tree lists (`res$trees`) are saved and were used directly for bipartition + topology counting.

Driver: `dev/red-team/audit-multirep-pathologies-2026-05-19.R`
Full output: `dev/red-team/audit-multirep-pathologies-2026-05-19.txt`
Numeric results: `dev/red-team/audit-multirep-pathologies-2026-05-19.rds`

## Results

### Rep 01 -- AWARE chain (n=240 saved, 180 post-burnin)

| Metric | Legacy | Corrected | Delta |
| --- | --- | --- | --- |
| TL median | 3.125 | 3.125 | (root-invariant) |
| TL last-100 mean | 3.040 | 3.040 | (root-invariant) |
| P(AC) | 0.0111 | 0.0222 | +0.0111 |
| P(AB wrong) | 0.7889 | 0.7889 | 0 |
| Unique topologies | 161 (legacy `topo_hash`, n=180) | 160 (canonical splits, n=180) | -1 |
| Top topology share | -- | 0.0222 | -- |
| Distinct root configs | -- | 31 | -- |
| Top root-config share | -- | 0.1722 | -- |

**Truth TL = 13.500.** Chain TL collapses to ~3, 4 to 5x undershoot.

### Rep 01 -- BLIND chain (n=210 saved, 157 post-burnin)

| Metric | Legacy | Corrected | Delta |
| --- | --- | --- | --- |
| P(AC) | 0.0000 | 0.0000 | 0 |
| P(AB wrong) | 0.8535 | 0.8535 | 0 |
| Unique canonical topologies | -- | 157 / 157 | (every tree distinct) |
| Distinct root configs | -- | 42 | -- |

### Rep 05 -- AWARE chain (n=243 saved, 182 post-burnin)

| Metric | Legacy | Corrected | Delta |
| --- | --- | --- | --- |
| TL median | 18.729 | 18.729 | (root-invariant) |
| TL last-100 mean | 20.706 | 20.706 | (root-invariant) |
| P(AC) | 0.0000 | 0.0000 | 0 |
| P(AB wrong) | 0.4615 | 0.4615 | 0 |
| Unique topologies | 126 (legacy `topo_hash`, n=182) | 111 (canonical splits, n=182) | -15 (-12%) |
| Top topology share | -- | 0.0604 | -- |
| Top 3 topology shares | -- | 0.060 / 0.060 / 0.028 | -- |
| Distinct root configs | -- | 16 | -- |
| Top root-config share | -- | 0.2363 | -- |

**Note:** Memory entry says "TL fine (~15-19)"; actual aware TL is somewhat inflated (median 18.7 vs truth 13.5; last-100 mean 20.7). Not a collapse but not "fine" either.

### Rep 05 -- BLIND chain (n=240 saved, 180 post-burnin)

| Metric | Legacy | Corrected | Delta |
| --- | --- | --- | --- |
| P(AC) | 0.0000 | 0.0000 | 0 |
| P(AB wrong) | **0.0000** | **0.9778** | **+0.9778** |
| Unique canonical topologies | -- | 168 / 180 | -- |
| Top topology share | -- | 0.0222 | -- |
| Distinct root configs | -- | 16 | -- |
| Top root-config share | -- | 0.2667 | -- |

**Massive legacy-scoring artefact in blind rep05.** Legacy P(AB) read 0 because the chain was rooted *inside* the AB tip set in essentially every sample, hiding the wrong-bipartition mode-lock entirely. Corrected scoring exposes it: blind rep05 is locked onto AB with 97.8% support. Same failure mode as v4cross-b (legacy 14 -> corrected 1 reported in the prior audit), only here the topology count is large because canonical topologies *do* vary on the AB-correct side -- the 97.8% refers to the bipartition, not a specific topology.

## Claim-by-claim survival

| Claim | Status under corrected scoring | Notes |
| --- | --- | --- |
| Rep 01 TL collapse (aware) | **SURVIVES** | TL is root-invariant; chain genuinely collapses to ~3 vs truth 13.5. AB-wrong support 0.79, AC 0.02, topology count high (~160 distinct). Picture is: branches collapse, characters saturate, low-information posterior, chain visits many trees but locks onto AB. Mode-trap diagnosis intact. |
| Rep 05 topology valley (aware) | **SURVIVES** (with caveats) | P(AC)=0 under both scoring rules; P(AB) 0.46 under both. Unique-topology count 161 -> 111 with corrected hashing -- some root-shuffling inflation, but not a 14 -> 1 collapse. The "valley" is real: chain visits >100 distinct unrooted topologies post-burnin with no topology exceeding 6% share, never touches AC. Memory's "TL fine" claim slightly off (actually inflated to ~19), but headline topology mixing failure is genuine. |

## New finding (not part of original brief)

**Blind chain rep05** showed a previously hidden mode-lock: P(AB) legacy 0.000 -> corrected 0.978. The blind diagnostic in the original mixing-followups would have read "neither AC nor AB" and been misclassified. Under corrected scoring it is a confident wrong-answer chain, similar in character to the v4cross-b artefact (root inside the target tip set obscured legacy scoring). This is consistent with the broader audit finding (commit `cb3f9ea` and predecessors) that blind chains in this regime mode-lock onto AB and the legacy diagnostics missed it.

## Implications for "PT is needed for mixing"

- Rep 01 (TL collapse) is **not a PT-fixable problem in the canonical sense**. It is a branch-length attractor: the collapsed-TL posterior with MAP nuisance beats the truth posterior by ~30 nats (per the original diagnostic). Heating the chain may help cross the TL valley but the underlying issue is in the TL prior / joint TL-rate proposal, not topology mixing per se. The corrected scoring does not alter this judgement.
- Rep 05 (aware topology valley) remains a **plausibly PT-fixable problem**: TL is in the right neighbourhood (if elevated), truth has +106 nats over chain MAP (per original diagnostic), and the chain visits >100 distinct unrooted topologies post-burnin without ever reaching AC. The valley is real, not a scoring artefact.
- The blind rep05 mode-lock on AB (newly visible) re-confirms that the bipartition-scoring contamination has been hiding mixing failures in the *blind* arm specifically. The "aware does better than blind" headline is unaffected in direction (blind reaches AB; aware reaches a mix that includes some AC-adjacent topologies), but blind quality is worse than the legacy diagnostic suggested.

## Files committed

- `dev/red-team/audit-multirep-pathologies-2026-05-19.R` (driver)
- `dev/red-team/audit-multirep-pathologies-2026-05-19.md` (this report)
- `dev/red-team/audit-multirep-pathologies-2026-05-19.txt` (per-rep output)
- `dev/red-team/audit-multirep-pathologies-2026-05-19.rds` (numeric results)
