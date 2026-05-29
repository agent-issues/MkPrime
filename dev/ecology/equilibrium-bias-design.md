# Ecology-aware redesign: bias the equilibrium, not the rate

Date: 2026-05-29. Context: the rate-multiplier `ecologyAware` model is mis-specified
(see memory `project_aware_rodent_misspecified`). This note records the redesign
that emerged from (a) my analysis, (b) an independent unprimed subagent, and
(c) Martin's neomorphic-pooling instinct — all three converged.

## Premise vs implementation (the pivot)

- **Implementation was broken.** Rate-modulation changes the *pace* of change, not
  its *destination*. It (1) couples to branch length (`t·φ` is confounded with `t`;
  Sim3 aware logL ~65 nats below blind, TL→7500), (2) is non-directional for
  transformational chars and sign-unidentifiable for neomorphic under relabel
  (sim3-simulate.R:92-110), and (3) had no point-mass-at-null, so it could not
  switch itself off → "goes to pot." A correct model's null case IS the baseline.
- **Topology-improvement premise is identifiability-limited in general**, but NOT by
  confounding on the rodent data. Separating convergence from homology needs ecology
  decorrelated from clade structure — and rodent ecology IS decorrelated: 17 changes
  on 60 tips, with 4-5 INDEPENDENT origins per non-reference state (state1=4, state2=5,
  state3=5; eco-origins-check.R). Earlier "clade-confounded (z=-4.4)" was a misread —
  z=-4.4 is only the significance of *mild* clustering (17 vs random 22.3, ~23% below),
  not an effect size. Ecology is dispersed = the FAVORABLE regime for identifiability.
  The real limit is per-origin DEPTH: state2's 5 origins are 5 singletons, so per-cell
  info is sub-nat even though dispersion is high → pooling across characters is the fix.
- **This sharpens "implementation not premise":** the data carried dispersed,
  identifiable ecology signal and the rate-mechanism still destroyed the tree. The
  signal was there to use; the mechanism squandered it.
- **The actionable-evolution premise is sound and is the right deliverable.**

## The convergent design: Ecology-Biased Equilibrium (EBE)

Each character evolves under F81 whose **stationary distribution depends on edge
ecology**: `π_{c,e} = softmax(η_c + β_{c,e}·s_c)`, where `η_c` = baseline log-freqs,
`s_c` = fixed direction contrast, `β_{c,e}` = signed ecology effect. Edge process is
the ecology-weighted mixture `P_c(t)=Σ_e w(e)·exp(t·Q(π_{c,e}))` (reuse current
mixture + the existing `qHeterogeneity` F81 builder; F81 exp is closed-form).

**Why it fixes the failure:**
- `β` changes destination, not pace → **largely decoupled from tree length** (no
  rate-time ridge; the main structural win the subagent identified). Not *perfectly*
  orthogonal — a strong bias mildly reduces variability — but far better than `t·φ`.
- **Sign-identifiable for neomorphic** (equilibrium freqs are not label-symmetric;
  the relabel symmetry that collapsed enc↔disc does not apply). Delivers "which
  direction." Requires NOT relabelling neomorphic chars.
- **Graceful degradation by construction**: `β=0` ⇒ `π_{c,e}=π_c` ⇒ exact baseline
  Mk, a prior point mass with zero likelihood cost. Cannot go to pot.

**Pooling (the small-data move, = Martin's instinct, correctly scoped):**
- Pool the **magnitude** (shared slab SD `τ` across cells) — safe, high-value,
  identified from all ~160 neomorphic chars.
- Keep the **sign sparse per-character** (`β_{c,e}` with spike-and-slab, π0 mode 0.9,
  ESS scaled to nChar). Do NOT pool the sign into one per-ecology value: a single
  per-ecology gain:loss ratio cannot represent claw-gain AND eye-loss in the same
  ecology (averages mixed-direction convergence to null). Optionally pool the
  *direction propensity* hierarchically per ecology → "is ecology e net-elaborating
  vs net-regressive."
- **Direction lives in neomorphic chars; multistate contributes magnitude only**
  (symmetric multistate sign is root/polarity-weak). Restrict multistate `s_c` to a
  fixed corpus-derived contrast or disable the ecology layer for them. Both Martin
  and the subagent reached "neomorphic for direction" independently.

**Ecology history:** observed tips, marginal-mixture over edges (keep current arch),
but upgrade the edge weight from parent-node marginal → expected edge **occupancy**
(stochastic-mapping dwell-time fraction). Output joint ancestral-ecology maps.

**Fitting:** F81 closed-form exp; a `β_{c,e}` move invalidates only char c's partials
(far more local than current all-partition invalidation → better mixing); RJ/Gibbs on
inclusion indicators + RW on β; `τ`,`π0`,`μ_eco` hyper-moves; no β-topology block
needed (orthogonality); init all β=0 = baseline. PT likely unnecessary.

## Assessment vs the two goals

**(i) Topology — now genuinely testable on rodents, modest and SAFE.** Removes the
rate-time pathology; won't hurt (point-mass null). Reduces convergence mis-grouping
where ecology has independent origins — which rodent ecology HAS (4-5 per state). So
whether it improves topology is now an empirical question the model can answer, not a
foregone failure. Caveat the other way: blind already recovers sensible biology, so
convergence may not be badly distorting it; the win, if any, is correcting specific
ecologically-convergent nodes, not a wholesale re-topology.

**(ii) Actionable evolution — the real payoff, and viable on rodents.** Signed β =
"which traits, which direction, which ecologies, how strongly," sign-identified for
neomorphic (the incumbent provably can't). Because ecology is dispersed (not
clade-confounded), the pooled per-ecology effect is decorrelated from descent and
estimable: per-ecology net direction = elaboration vs regression (fossorial eye/orbit
reduction vs digging-apparatus elaboration). Pooled τ/π0/inclusion ranks are
trustworthy; per-CHARACTER cell claims are still limited by shallow per-origin depth
(state2 = 5 singletons) — which is exactly why pooling across the 160 neomorphic
characters is the right move. Plus ecology-conditioned ASR.

## Main NEW cost of going directional

Equilibrium bias makes the process **non-reversible** → introduces root-placement
sensitivity that reversible Mk avoids. Manageable for neomorphic (gain/loss asymmetry
already informs the root) but real for symmetric multistate. This is the principal
trade the equilibrium approach makes that the (broken) rate approach did not.

## Secondary design (subagent, ranked #2): ecology-pinned covarion/Markov-modulated

Hidden rate class tied to observed ecology. Rejected as primary: still rate-based
(re-creates rate-time coupling), directionless (fails goal ii), and covarion is
non-identifiable for binary chars (Allman-Rhodes-Sullivant 2009: #classes <
#observable states) — fatal for a 160-neomorphic matrix. Right choice only if the
question were "does ecology change tempo" not "does ecology bias toward states."

## Recommendation

Build EBE as a **neomorphic-only directional equilibrium layer** with shared-τ
pooling and a spike-at-0. The rodent dataset IS a legitimate testbed (dispersed
ecology, 4-5 origins/state) — pair it with a simulation that matches its structure
(dispersed regimes, shallow per-origin depth) to validate recovery. Frame the paper
around the actionable coevolution/ASR result (which traits, which direction, per
ecology); treat topology improvement as a secondary, empirical question.

ARCHITECTURAL COST to weigh before building: edge-dependent equilibria make the
process NON-STATIONARY → Felsenstein's pulley principle is lost, root placement
becomes part of the inference, and the current unrooted/reversible likelihood
machinery needs real change (not a kernel tweak). qHeterogeneity's F81 builder is a
starting point but `ecologyAware + qHeterogeneity` is currently forbidden in code
(beta_scale not threaded) — so adaptation, not free reuse. Also: equilibrium bias
models "drift TOWARD a state given a sustained regime," not an instantaneous switch;
if the target concept is literal switching, threshold/liability is the closer match.
