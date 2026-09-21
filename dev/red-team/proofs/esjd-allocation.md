# ESJD/s proposal scheduling — allocation theory, invariance, and measurement

**Lane:** math-prover (Fable-tier checkpoint "after this plan, before S1",
plan §7).
**Reviews:** `dev/plans/2026-08-12-esjd-proposal-scheduling.md` (assessment)
and `dev/plans/2026-08-12-esjd-implementation-plan.md` (the plan). This
document is the gate on whether S1 proceeds as specified.
**Code base:** worktree `esjd-proposal-scheduling-13cf5b`, branch
`claude/esjd-proposal-scheduling-13cf5b`. All `file:line` citations are
against this tree.
**Relation to prior proofs:** builds on
`dev/red-team/proofs/hastings-tree-moves.md` (cited below as HTM; its
line numbers are stale, mine are current) and
`dev/red-team/proofs/marginal-k-gibbs-p.md`.
**Numerical checks:** both counterexamples in this document (Q1 §1.4,
Q4 §4.5) were verified by one-off exact computation (16-state and 4-state
transition matrices, stationary distributions and asymptotic variances via
the fundamental matrix). Numbers quoted are from those runs. No test files
were left behind.

---

## 0. Verdict table

| Q | Question | Verdict |
|---|---|---|
| Q1 | Is ESJD the right criterion under multimodality? | **The multimodal motivation is overreach — proven.** One-step ESJD/s is a within-mode efficiency criterion; a verified counterexample shows the ESJD-optimal mixture can be worse than a balanced one by an unbounded factor on the asymptotic variance of the mode indicator, and the ESJD optimum *actively minimises* the crossing term whenever within-mode diffusion is rich. ESJD still strictly dominates the acceptance criterion (it credits $a D^2$ where acceptance credits $a$), so adoption is justified — but as a within-mode criterion, with an explicit **irreducibility insurance floor** on large-jump topology kernels and S0's island-switch rate as the multimodal arbiter. §1. |
| Q2 | Sum vs mean vs bottleneck-weighted credit within a block? | **Sum is a throughput objective and misaligns with the package's bottleneck gates; the LP-dual-weighted sum is correct to first order.** Mean-per-touched-coordinate is wrong in both directions. Sum is acceptable for v1 *only if* S1 instruments per-(move, coordinate) jump sums so S2 can test whether dual weighting flips rankings. Caveat: the shipped stopping rule never sees `br_` or `kPrime_` coordinates at all (`R/Convergence.R:271–272,57`), which weakens the plan's "matches the stopping rules" claim in both directions. §2. |
| Q3 | Joint moves measured on their primary block only? | **The measure-block-coordinates-only rule is coherent under the Q2 bottleneck logic; the bias is a bounded understatement (factor $1 + \hat\sigma^2_{z_1}/\hat\sigma^2_{z_2}$, = 2 at equal whitened scales) relative to a throughput view.** But the assignment must be a hand-maintained map: the registry's `target` field says `"tree_length"` for all three joint moves (`R/RunMkPrime.R:3810,3815,3821`), which would mis-assign them to the branch block. Assign to **rates** (as the plan's §0 table already does, contradicting its own §0 text) and record both jump components in S1. §3. |
| Q4 | Is the fixed-weight invariance argument airtight? | **The algebra is right; the premise is false.** (i) `gibbs_spr` is **not π-invariant** — proven, not a caveat; this upgrades HTM §5's "watertight with caveats". It is default-on. (ii) `gibbs_subtree_swap` lacks the neighbourhood-normaliser correction (confirms and sharpens HTM §6's open question). (iii) `weighted_subtree_swap`'s Hastings ratio is approximate on two counts (confirms HTM §8). (iv) Both slice samplers use a per-side capped stepping-out that violates the interval-selection symmetry slice sampling requires — new finding, exactness fails on cap-binding events. (v) The M-159 cache boost (default `cacheBonus = 5`) makes the mixture weights **state-dependent**, which voids $\pi(\sum_i w_i K_i) = \pi$ even for today's sampler — verified counterexample. It also biases the S1 measurements. **Gate condition: fix or drop `gibbs_spr`/`gibbs_subtree_swap` before S3, because ESJD/s will *up-weight* them** (always-accept, large RF jumps). "No SBC needed" holds only for the weight-rule delta conditional on these fixes. §4. |
| Q5 | Whitening from the warmup snapshot buffer? | **The diagnostic-not-kernel distinction holds** (scales feed only the score; `jointRhos` by contrast are kernel parameters, and both freeze at the Sample boundary). Legitimacy requires: a ≥50-row gate mirroring `.EstimateJointRhos`, an explicit ε-definition of "no usable scale" (σ̂→0 sends whitened jumps to +∞, not 0 — the plan's wording has it backwards), robust scale estimation or winsorised jumps, and awareness that warmup-transient drift *inflates* σ̂ and so *deflates* the score of drifting coordinates. §5. |

**S1 additions this review requires** (all cheap; consolidated list in §6):
per-(move, coordinate) jump ledgers for the branch and k′ blocks; split
ledgers by the M-159 `useBoost` flag; a crossing-tail counter (accepted
RF ≥ θ) per topology move; edge-identity (not row-index) keying for branch
jumps; both jump components for joint moves.

---

## 1. Q1 — ESJD under multimodality

### Theorem (informal)

Maximising one-step ESJD/s over mixture weights is equivalent to
minimising lag-1 autocorrelation per second, which controls relaxation
time only in regimes where the whole autocorrelation spectrum is slaved
to its first coefficient. On a two-island target this fails: the
ESJD-optimal mixture can have asymptotic variance for the island
indicator larger than a balanced mixture's by an arbitrary factor, and
larger precisely *because* it maximises ESJD.

### Assumptions

1. All candidate kernels $K_m$ are individually π-invariant (see §4 for
   where this fails in the actual registry; here we grant it).
2. The chain is at stationarity when ESJD is measured (warmup
   stabilisation, `R/RunMkPrime.R:1270–1284`, is a proxy for this).
3. "Mixing" is quantified by the asymptotic variance
   $\sigma^2_{\mathrm{as}}(f) = \lim_T T\,\mathrm{Var}(\bar f_T)$ of
   functionals $f$ of interest, equivalently by integrated
   autocorrelation time (Geyer 1992).

### 1.1 What ESJD is

For a stationary chain with kernel $K$ and any square-integrable $f$,

$$\mathrm{ESJD}(f) \;=\; \mathbb{E}\,|f(X_1)-f(X_0)|^2 \;=\; 2\,\mathrm{Var}_\pi(f)\,\bigl(1-\rho_1(f)\bigr),$$

so maximising ESJD of $f$ minimises the lag-1 autocorrelation of $f$
(Pasarica & Gelman 2010, eq. 3). The plan's criterion uses ESJD of the
*state embedding* — whitened coordinates for continuous blocks, RF for
topology — i.e. the sum of ESJD over a coordinate basis, not of the
functional a user cares about.

ESJD is a **linear functional of the kernel**:
$\mathrm{ESJD}(K) = \int \pi(dx)\,K(x,dy)\,d^2(x,y)$, hence linear in
mixture weights $w$. Any one-step, kernel-linear statistic shares the
blindness proven below; this is structural, not fixable by a better
metric $d$.

### 1.2 When maximising ESJD/s does maximise mixing

For a reversible $K$ with spectral measure $E_f$ of $f$
(spectral theorem for self-adjoint Markov operators; e.g. Geyer 1992 §3):

$$\mathrm{ESJD}(f) = 2\!\int (1-\lambda)\, dE_f(\lambda), \qquad
\sigma^2_{\mathrm{as}}(f) = \int \frac{1+\lambda}{1-\lambda}\, dE_f(\lambda).$$

ESJD weights each eigencomponent by $(1-\lambda)$;
$\sigma^2_{\mathrm{as}}$ weights it by $(1+\lambda)/(1-\lambda)$. The
criterion therefore **down-weights exactly the components on which the
target quantity diverges** ($\lambda \to 1$). Maximising the first
controls the second only when:

- **(a) Eigenfunction regime.** $E_f$ is (nearly) a point mass — $f$ is
  close to an eigenfunction — so $\rho_k \approx \rho_1^k$ and both
  quantities are monotone in $\rho_1$. Within one well-behaved mode this
  is the usual working approximation.
- **(b) Time-rescaling family.** The candidate kernels differ only by a
  common time rescaling of one limiting dynamics, so *all*
  autocorrelation functions scale together. This is exactly the
  Roberts–Gelman–Gilks (1997) diffusion-limit regime in which the
  0.234/ESJD machinery has its optimality pedigree: product targets,
  $d\to\infty$, a random-walk family indexed by step size. Changing the
  step size only changes the speed measure of one Langevin diffusion,
  so maximising ESJD maximises speed and minimises $\tau_{\mathrm{int}}$
  for every functional simultaneously (see also Roberts & Rosenthal
  2001; Sherlock & Roberts 2009).

Neither condition holds for a mixture of *heterogeneous* kernels on a
*multimodal* target: different mixture weights change the shape of the
spectrum, not just its scale, and the mode indicator's spectral mass
sits at $\lambda^\ast \to 1$ where the $(1-\lambda)$ weight vanishes.

### 1.3 Decomposition: within-mode ESJD masks zero conductance

Let $\{A, A^c\}$ be a two-set decomposition of the state space. Split

$$\mathrm{ESJD}(K) = \underbrace{\int_{x,y \text{ same side}} \pi(dx) K(x,dy) d^2(x,y)}_{\text{within}} \;+\; \underbrace{\int_{x,y \text{ cross}} \pi(dx) K(x,dy) d^2(x,y)}_{\text{cross}}.$$

Only the cross term relates to the conductance
$\Phi = \int_A \pi(dx)K(x,A^c)/\pi(A)$ that bounds the spectral gap
(Lawler & Sokal 1988; Sinclair 1992):
$\mathrm{cross} \le d^2_{\max}\,\pi(A)\,\Phi$ and
$\mathrm{cross} \ge d^2_{\min,\mathrm{cross}}\,\pi(A)\,\Phi$. The within
term is unconstrained by $\Phi$. Hence a mixture can have arbitrarily
large total ESJD and arbitrarily small conductance: the sum is
maximised by pouring weight into within-mode movers whenever their
per-second $d^2$ yield beats $a\,D^2$ of the crossing kernel.

### 1.4 Sharp counterexample (verified numerically)

**Construction.** States $(s,i)$, island $s\in\{0,1\}$, within-island
coordinate $i\in\{1..n\}$; $\pi$ uniform on the $2n$ states. Embedding
$x(s,i)=(sD, i)$, $d^2$ = squared Euclidean. Two kernels, both exactly
π-invariant:

- $K_{\mathrm{loc}}$ ("nni-like"): resample $i \sim U\{1..n\}$ within the
  current island (a Gibbs update). Per-draw ESJD $= (n^2-1)/6$.
- $K_{\mathrm{cross}}$ ("tbr-like"): with probability $a$ propose the
  mirror state $(1-s, i)$ (accepted, $\pi$ uniform), else hold. Per-draw
  ESJD $= a D^2$. The hold branch models proposals burnt in the
  inter-island valley; $a$ is the effective crossing acceptance.

Mixture $K_w = w K_{\mathrm{cross}} + (1-w) K_{\mathrm{loc}}$, equal
per-draw costs. Then

$$\mathrm{ESJD}(w) = w\,aD^2 + (1-w)\,\tfrac{n^2-1}{6},$$

linear in $w$: **whenever $aD^2 < (n^2-1)/6$, the ESJD-optimal $w$ is
the floor.** For the island indicator $f = \mathbf 1\{s=1\}$, the
$s$-process flips with probability $\varepsilon = wa$ per step, so
$\rho_k(f) = (1-2wa)^k$ and

$$\sigma^2_{\mathrm{as}}(f) = \tfrac{1}{4}\cdot\frac{1-wa}{wa} \;\approx\; \frac{1}{4wa}.$$

**Numbers** ($n=8$, $D=4$, $a=0.05$; exact computation on the 16-state
chain, both kernels checked π-invariant to $10^{-10}$):

| $w$ | ESJD | $\sigma^2_{\mathrm{as}}(f)$ |
|---|---|---|
| 0.01 (floor) | **10.40** (max) | **499.8** (worst) |
| 0.10 | 9.53 | 49.8 |
| 0.50 | 5.65 | 9.75 |
| 0.90 | 1.77 | 5.31 |

The ESJD-optimal mixture is 51× worse than $w=1/2$ on the functional
that defines the benchmark (island occupancy), and the ratio
$\sigma^2_{\mathrm{as}}(w_{\mathrm{floor}})/\sigma^2_{\mathrm{as}}(1/2)
\approx 1/(2w_{\mathrm{floor}})$ is unbounded as the floor shrinks.
Increasing $n$ (richer within-island space) makes the ESJD ranking
*more* lopsided against the only crossing kernel while leaving the
crossing problem unchanged — exactly backwards.

**Interpretation for MkPrime.** Substitute nni for $K_{\mathrm{loc}}$
(accepted jump deterministic RF² = 4, acceptance driven to 0.23 by
`.AdaptTuning`, `R/RunMkPrime.R:4790`) and tbr's island-crossing tail for
$K_{\mathrm{cross}}$ ($a \sim 10^{-2}$, RF² up to $(2(n_{\mathrm{tip}}-3))^2$).
Realistic numbers put the two ESJD contributions within an order of
magnitude, so starvation is *possible, not automatic* — an empirical
question S2 legitimately answers. What the counterexample establishes is
that the criterion has **no safeguard**: nothing in ESJD/s prevents the
crossing kernel being floored in precisely the regime where it is the
only thing that matters.

### 1.5 Fairness: ESJD still dominates the status quo

The current score credits a move with acceptance $a$ by $a \cdot
\mathrm{dim}/\mathrm{cost}$ (`R/RunMkPrime.R:4468`); it is blind to jump
magnitude, so tbr earns $a$ where ESJD earns $a D^2$ — a factor $D^2$
less blind. `.DecayLowAcceptMoves` (`R/RunMkPrime.R:4556–4621`)
additionally *punishes* $a < 0.02$ multiplicatively; its retirement
(plan §S3) is unambiguously right. So the redesign improves the
multimodal situation; it just does not fix it, and the plan should not
claim it does.

### 1.6 What repairs it

1. **Lag-$k$ / multi-step ESJD** — $\mathbb E\, d^2(X_t, X_{t+k}) =
   2\int(1-\lambda^k)dE$ reweights the spectrum usefully, but a $k$-step
   displacement is a property of the *mixture path*, not attributable to
   a single move. It can serve as a block-level metastability detector
   (slow growth in $k$ ⇒ trapped), not as a per-move allocation signal.
   Not recommended for the score.
2. **Crossing-tail counter (recommended, S1).** Per topology move,
   accumulate the count (and RF² sum) of accepted jumps with
   RF ≥ θ (θ ≈ 4, i.e. more than one bipartition). This is
   kernel-linear, per-move attributable, and estimates the *flow* term
   of §1.3 directly rather than diluting it in within-mode diffusion.
   One extra accumulator at the same instrumentation point.
3. **Floors as insurance, not measurability.** The plan's draw-count
   floors `min(kMinDraws/intervalDraws, 1/nValid)` (plan §S3) make the
   *share* shrink as intervals lengthen — correct for keeping a move
   measurable, wrong for irreducibility insurance, where the quantity to
   floor is the expected number of large-jump *attempts per unit time*.
   Keep both, separately justified: draw-count floors for measurement on
   all moves; a fixed share floor on the large-jump topology kernels
   (`tbr`, `spr`, `pspr`) stated as insurance. Since `nChains = 1L` is
   the default (`R/MkPrimeMCMC.R:342`), there is **no tempering backstop
   in default runs** — the move mixture is the only crossing mechanism,
   which is what makes the insurance floor non-optional. Under
   `nChains > 1` the ladder does the crossing (swap kernel,
   `src/mcmc.cpp:6201–6212`) and ESJD-within-mode is the right thing to
   tune; that argument does not extend to the default configuration.
4. **S0's island coverage / switch rate** (plan §S0 metrics) is the
   correct multimodal arbiter and must remain the gate for the topology
   block — never tree-ESJD alone, and (per assessment §5) never
   per-chain ESS.

### Edge cases

- $a \to 0$ (islands unreachable at any weight): no allocation fixes
  this; only tempering or better kernels. The insurance floor's cost is
  bounded by its share; its benefit is unbounded.
- Unimodal target: cross term vanishes, §1.2(a) applies, ESJD/s is the
  right criterion — the redesign's home turf.
- $k'$ block: integer walks on a near-unimodal conditional — ESJD fine.
- Topology posterior concentrated on one tree: all topology ESJD → 0;
  weights fall to floors; harmless.

### Verdict (Q1)

**The plan's §0 within-block ESJD/s is adoptable as a within-mode
efficiency criterion. The multimodal motivation is overreach and should
be reworded in the plan.** Required compensations: crossing-tail counter
in S1; insurance share floor on large-jump topology kernels in S3
(justified as irreducibility insurance, explicitly not by measurement);
S0 switch-rate as the topology-block gate. With those, "ESJD within
blocks" is coherent and strictly better than the status quo.

---

## 2. Q2 — Credit assignment within a block

### Theorem (informal)

If the objective is the package's own bottleneck criterion (min-ESS /
max-R̂ over parameters), the within-block score that is optimal to first
order is a **dual-weighted sum** $\sum_j \lambda_j e_{jm}/c_m$ with
weights $\lambda$ concentrated on the worst-mixing coordinates. The
unweighted sum (plan §1) is the special case $\lambda$ = uniform and can
be pessimal by a factor equal to the block's mixing heterogeneity.

### Assumptions

1. Per-draw expected whitened squared jump of move $m$ on coordinate
   $j$, $e_{jm} \ge 0$, and per-draw cost $c_m$, are stable over the
   adaptation window.
2. Proxy: coordinate-$j$ mixing rate under mixture $t$ (seconds
   allocated per move) is $R_j(t) = \sum_m (t_m/c_m)\, e_{jm}$. This
   uses ESJD's exact additivity over mixtures; ESS is not additive, so
   $R_j$ is a first-order proxy at the same fidelity level as the plan's
   entire approach (ESJD ≈ monotone in ESS via
   $\mathrm{ESS} \approx N(1-\rho_1)/(1+\rho_1)$, valid in the §1.2(a)
   regime).

### 2.1 The LP and its dual

Maximise $\min_j R_j(t)$ s.t. $\sum_m t_m = T$, $t \ge 0$. Standard LP;
at the optimum there exist $\lambda_j \ge 0$, $\sum_j \lambda_j = 1$,
supported on the binding (bottleneck) coordinates, such that every
funded move maximises the score

$$\sigma_m = \sum_j \lambda_j\, \frac{e_{jm}}{c_m}.$$

So the *form* of the right numerator is a weighted sum over coordinates;
the only question is $\lambda$.

### 2.2 Failure modes of each candidate rule

Let two moves have equal cost, two coordinates:

- **Sum** ($\lambda$ uniform). $e_1 = (9,\, 0.1)$, $e_2 = (1,\, 1)$. Sum
  picks $m_1$ (9.1 > 2); bottleneck rate 0.1 vs 1 — **10× worse** on the
  binding coordinate. This is precisely the "moves 3 edges hard vs all
  105 modestly" pattern the brief names, one level below the
  `dim`-vs-bottleneck tension the plan §1 removes at block level.
- **Mean per touched coordinate.** Same example: $m_1$ touches both
  (mean 4.55 > 1) — still wrong. Reverse failure: specialist
  $e_1 = (5, 0, \dots, 0)$ (mean 5) vs sweep $e_2 = (1, \dots, 1)$
  (mean 1): mean picks the specialist and starves $d_B - 1$
  coordinates. Mean is not the dual of any objective here; it penalises
  breadth and rewards concentration simultaneously. Reject.
- **Bottleneck-weighted sum** ($\lambda_j \propto 1/(\delta + \hat R_j)$
  with $\hat R_j$ the coordinate's accumulated whitened jump rate over
  the window, $\delta$ a noise floor). Approximates the LP dual by
  multiplicative reweighting. Failure modes: (i) noise-chasing — duals
  are extreme points, and $\hat R_j$ at warmup window lengths is noisy;
  mitigate with $\delta$, smoothing across batches, and the same
  annealed softmax the score already passes through; (ii) if *no* move
  differentiates on the bottleneck coordinate ($e_{jm}$ ≈ equal ∀m), the
  weighting is harmless (adds a constant); the LP handles this by
  spreading duals to the next constraint.

### 2.3 Reality check against the shipped stopping rule

The plan's premise "stopping rules are bottleneck criteria over all
parameters" is only two-thirds true. `ConvergenceDiagnostics` computes
`minEss`/`maxRhat` over `.KeyParamCols` = `^(log_|tree_|rate_|p$|kPrime_)`
(`R/Convergence.R:271–272`) and then *excludes* `kPrime_` and
`log_likelihood` from the min (`R/Convergence.R:57,76`). Branch
coordinates (`br_`) never enter. The tuning bandit's objective likewise
excludes `^(kPrime_|br_|log_likelihood)`
(`R/RunMkPrime.R:4983–4988`). Consequences:

- Within the **branch block**, the bottleneck justification is indirect
  (branch mixing feeds `tree_length`, `log_posterior` and tree-ESS, all
  of which *are* gated) — still real, but the clean "matches the
  stopping rule" argument is unavailable. The 105-coordinate bottleneck
  question is a mixing-quality question, not a literal gate question.
- Across blocks, the repaired min-ESS/s bandit inherits the same
  blindness: it cannot see a branch or k′ bottleneck except through
  tree-ESS (folded in at ≥20 trees, `R/RunMkPrime.R:4999–5008`). The
  plan should state this scope explicitly rather than claim full
  bottleneck coverage.

### Verdict (Q2)

**Sum is acceptable for v1 only as an explicitly provisional choice**,
with S1 instrumenting per-(move, coordinate) jump sums for the branch
and k′ blocks (≈ `nMoves × nEdge` doubles — trivial) so that S2 can
compute both rankings (uniform-λ vs inverse-rate-λ) from the same run.
If they disagree on the branch block, S3 should ship the
bottleneck-weighted sum with $\delta$-smoothed, annealed weights. Mean
per touched coordinate is rejected outright. Note for the topology
block: RF² is a block-scalar so this question does not arise there —
though RF is itself a sum over bipartitions and the same phenomenon
exists one level down (a move can churn cherry bipartitions while never
moving deep ones); no action for v1 beyond the crossing-tail counter
(§1.6).

---

## 3. Q3 — Joint moves measured on their primary block only

### Statement

`joint_tl_rls`, `joint_tl_rl`, `joint_tl_rn` propose correlated 2D
log-scale Bactrian perturbations with a shared step size and estimated
correlation `jointRhos` (`src/mcmc.cpp:5117–5127, 5129–5137, 5384–5390`;
both axes use the same `scaleTuning`, `:5120–5121`; Hastings
$= \log m_1 + \log m_2$, `:5126`). The plan (§0) measures each joint
move's jump over one block only.

### 3.1 Bias quantification

Whitened jump components: $\Delta z_1^2 = (\text{scale}^2 z_1^2)/\hat\sigma^2_{\log TL}$,
$\Delta z_2^2 = (\text{scale}^2 z_2^2)/\hat\sigma^2_{\log X}$ with
$\mathbb E z_1^2 = \mathbb E z_2^2$ (symmetric 2D Bactrian) and both
coordinates accepted or rejected together. Hence the expected
understatement of total whitened ESJD from crediting only coordinate
$k$ is the factor

$$1 + \frac{\hat\sigma^{2}_{\text{other}}{}^{-1}}{\hat\sigma^{2}_{k}{}^{-1}}
\Big|_{\text{whitened}} = 1 + \frac{\hat\sigma^2_{\log \cdot k}}{\hat\sigma^2_{\log \cdot \text{other}}},$$

which is exactly **2 when the two whitened scales are equal**, and
bounded only by the scale ratio in general. A joint move genuinely 1.5×
better than `slice_rate_log_sd` on total whitened ESJD/s scores 0.75×
under single-coordinate credit and gets annealed toward the floor — a
ranking flip within the plausible range.

### 3.2 Why block-coordinates-only is nonetheless the right rule

Under the Q2 bottleneck logic, a block's contest should be scored on
*that block's coordinates*: crediting a rates-block competitor with
tree-length displacement would let a joint move dominate the rates
block while barely moving rates (option (b) of the brief), which is a
worse failure than understatement. The clean statement is:

> A move's score in block $B$ is its whitened jump summed over $B$'s
> coordinates only (0 for coordinates outside $B$).

This makes `joint_tl_rls` vs `slice_rate_log_sd` an apples-to-apples
comparison on `rate_log_sd` displacement per second — precisely "does
the ridge-aligned move actually mix the rate coordinate better than the
axis-aligned one", which is the question the joint moves exist to
answer (`dev/notes/2026-05-27-rate-neo-ridge-and-joint-moves.md`).
The residual bias is a positive externality (tree-length displacement
credited nowhere), acceptable in v1 because `tree_length` has dedicated
movers and *is* seen by the across-block bandit
(`tree_` ∈ `.KeyParamCols`).

### 3.3 The assignment bug waiting to happen

Plan §0 text says "assign each move to the block of its *primary*
target", and plan §2 derives the block table from the moves' declared
`target` fields. But the registry declares `target = "tree_length"` for
all three joint moves (`R/RunMkPrime.R:3809–3811, 3814–3816,
3820–3823`), which under the plan's own §0 table maps to the **branch
block** — while the same table lists `joint_tl_*` under **rates**, and
the display categoriser also files them under Rates
(`R/RunMkPrime.R:4700–4704`). Resolution:

1. The block table must be a **hand-maintained move→block map**, not
   derived from `target` (the plan §2 already concedes this for
   topology/branch mutation; extend it to joint moves).
2. Assign each `joint_tl_x` to the block of its **non-TL coordinate**
   (rates) — that is where the historically slow coordinates live
   (`rate_log_sd` was the marginal SBC parameter) and where the
   marginal movers cannot traverse the ridge, so under-crediting the
   joint move there is the costly mistake.
3. S1 must record **both** whitened components per joint-move draw
   (2 doubles) so the assignment and the understatement factor are
   auditable rather than assumed.

### Edge cases

- `jointRho = 0` (no estimated correlation, `R/RunMkPrime.R:3398–3443`
  returns 0 below 50 rows): the joint move degenerates to two
  independent scale moves; single-coordinate credit understates by
  exactly its independent other-axis displacement; harmless.
- $|\rho| \to 0.95$ (clamp, `:3411`): components remain equal in
  expectation under the shared step size; factor-2 understatement holds.

### Verdict (Q3)

**Plan's rule is coherent; keep it, with the assignment corrected to a
hand-maintained map (rates block for all three) and both components
instrumented in S1.** The systematic understatement is bounded and
directionally safe under bottleneck scoring; option (b) (full-sum
credit inside one block) is rejected.

---

## 4. Q4 — Invariance of the mixture, kernel by kernel

### Statement under review

Plan §5: "Every registered move is a complete π-invariant kernel, so
for any fixed weights $\pi(\sum_i w_i K_i) = \sum_i w_i (\pi K_i) =
\pi$. … MkPrime freezes weights at the warmup boundary, so the sampling
phase remains a fixed-kernel MCMC," hence no diminishing-adaptation
argument and no SBC campaign.

### 4.1 The algebra and the freeze boundary — correct, with one wording fix

The displayed identity is correct for state-independent weights and
individually invariant kernels; mixtures of invariant kernels are
invariant, and mixtures of reversible ones are reversible (Tierney
1994 §2.4). The freeze is real but sits at the **Tuning→Sample**
boundary, not the warmup boundary: weights are adapted during Warmup
(`.AdaptMoveWeights`/`.DecayLowAcceptMoves`/`.WarmupGibbsCap`,
`R/RunMkPrime.R:1225–1252`) and perturbed by the bandit during Tuning
(`:1404, 1419`), then fixed on entry to Sample (`:1419–1431`).
Kernel-shape parameters freeze earlier or simultaneously: scale/slice
tunings adapt only in the Warmup branch (`:1194–1201`), `jointRhos`
in Warmup and Tuning (`:1259–1265, 1372–1379`), temperature ladder in
Warmup (`:1204`). No adaptation call is reachable in the Sample branch,
and only Sample-phase draws are retained (warmup saves nothing,
`R/RunMkPrime.R:3376`; tuning fills a discarded buffer, `:1150–1167`).
So: sampling phase is fixed-kernel **provided the premise holds**.
Plan §5's "warmup boundary" should read "tuning→sample boundary", and
the assert-the-boundary test (plan §5 item 2) is endorsed.

The rest of this section is about the premise.

### 4.2 `gibbs_spr` is not π-invariant — proof

**Construction** (`gibbs_spr_impl_full`, `src/mcmc.cpp:1801–1941`;
the het variant shares the sampling/commit semantics, `:1494, 1748`).
Prune edge $(u,v)$ chosen uniformly among non-root edges
(`:1806–1814`); for each regraft candidate $r$, the vacated edge pair
is merged to $\ell_{\mathrm{merge}}$ and the candidate edge $\ell_r$
is split **exactly in half**: `absLen[rr] = 0.5*lReg; absLen[sibRow] =
0.5*lReg` (`:1880–1882`); one of {self, candidates} is drawn
$\propto e^{\beta LL}$ (`:1896–1913`) and, if a candidate, **committed
without any MH step**, halves included (`:1921–1934`, `state->logLik`
overwritten at `:1937`).

**Claim.** $\pi K_{\mathrm{gspr}} \ne \pi$ for any target $\pi$ whose
branch-length marginal is absolutely continuous (true here: flat
Dirichlet on `rel_br_lengths` — `src/mcmc.cpp:2712–2714` — times an AC
tree-length prior).

**Proof.** Let $E$ = the event that some internal node subtends two
edges of exactly equal length. $E$ is a finite union of hyperplanes in
the branch-length simplex, so $\pi(E) = 0$. Whenever
$K_{\mathrm{gspr}}$ commits a non-self candidate, the accepted state
lies in $E$ (the two halves of the split edge are equal by
construction). The non-self probability is positive with positive
π-probability (any state where some candidate has finite likelihood).
Hence $(\pi K_{\mathrm{gspr}})(E) > 0 = \pi(E)$, so
$\pi K_{\mathrm{gspr}} \ne \pi$. The same holds for the full mixture
$K = \sum_i w_i K_i$ with $w_{\mathrm{gspr}} > 0$: no other kernel can
cancel a one-step discrepancy. $\square$

Equivalently in reversibility terms: from an accepted state $y$,
re-pruning the same edge can never regenerate $x$ (the reverse
candidate at $x$'s site carries the 50/50 split of
$\ell_{\mathrm{merge}}$, not $x$'s original fractions), so
$K(y,\{x\}) = 0$ while $K(x,\{y\}) > 0$.

**Relation to the prior proof.** HTM §5 identified the fixed
$\tau = 1/2$ and graded it "Watertight with caveats", reasoning that
the move is a Gibbs draw "conditional on a hard-coded mid-edge regraft"
and that separate branch moves make this "unproblematic". That
resolution is a composition error: a kernel invariant for a *different*
distribution (mass on the mid-edge manifold) mixed with π-invariant
kernels does not preserve π — mixture invariance requires every
component to preserve the *same* π. The measure argument above closes
what HTM left as a caveat. **HTM §5's verdict should be upgraded to
Implementation bug.**

**Scope.** `gibbsSpr` defaults to `TRUE` (`R/MkPrimeMCMC.R:356`); the
move is registered whenever topology is free and `nEdge ≥ 5`
(`R/RunMkPrime.R:3570–3575`), except under `marginal_k` where it is
disabled for unrelated reasons (`:3555–3568`). Default `sampled_k`
posteriors therefore currently carry this bias. Its magnitude is an
empirical question (every accepted `gibbs_spr` equalises one edge pair;
subsequent branch moves re-diffuse, so the stationary distortion is a
smoothed pull of adjacent edge-length pairs toward equality, plus
whatever that couples to). Testable corollary: the Hamilton SBC's
marginal `tree_length` (AD = 0.33) and `rate_log_sd` (AD = 0.26)
results were attributed to mixing; rerunning that arm with
`gibbsSpr = FALSE, gibbsSubtreeSwap = FALSE` would separate
bias-from-this-kernel from slow mixing.

> Caveat (4.2): magnitude of the `gibbs_spr` stationary bias on
> realistic data. Closing this requires either an exact small-case
> stationary computation or an SBC/posterior A-B with the move
> disabled. Should be picked up by lane mcmc-diagnostician (SBC arm) /
> numerical-auditor (small-case exactness).

**Fix sketch (non-trivial, no patch attached):** draw the regraft
fraction $f \sim$ Beta centred on ½ (or the M-087 bin scheme) and
MH-correct — i.e. make `gibbs_spr` a special case of `weighted_spr`,
whose construction is exact (§4.4). Alternatively drop `gibbs_spr` from
the default schedule; `weighted_spr` (currently default `FALSE`,
`R/MkPrimeMCMC.R:362`) is its corrected version.

### 4.3 `gibbs_subtree_swap` — normaliser asymmetry, confirms HTM §6

Construction (`gibbs_subtree_swap_impl_full`, `src/mcmc.cpp:2470–2575`):
node $A$ uniform via edge row (`:2475–2477`), partners enumerated,
candidates = swap($A$,$p$) with branch lengths travelling with the
subtrees (`:2510–2513`, commit `:2554–2558`) — measure-preserving, so
no §4.2-style collapse. Selection $\propto e^{\beta LL}$ over
{self} ∪ candidates with **no MH correction** (`:2538–2546`, commit
`:2548–2571`).

Detailed balance requires $Z(x; A) = Z(y; \cdot)$ where $Z$ is the
neighbourhood normaliser; the neighbourhoods
$\{x\} \cup \{\mathrm{swap}(x;A,p)\}$ and
$\{y\} \cup \{\mathrm{swap}(y;A,p')\}$ are not equal sets (a second
swap from $y$ composes two swaps), so the normalisers differ
generically. Proportional-to-π selection over a state-dependent
neighbourhood that does not partition the space is not a Gibbs update
and needs the normaliser-ratio acceptance
$\min\{1, Z(x)/Z(y)\}$ (Metropolised Gibbs; cf. the "locally balanced"
framework of Zanella 2020 — always-accept requires the balancing
function's normaliser to cancel, which it does not here). HTM §6
reached the same equation and left it open; nothing in the current code
closes it, and there is no mechanism that would make the kernel
invariant without detailed balance. **Presumed not invariant;
default-on** (`gibbsSubtreeSwap = TRUE`, `R/MkPrimeMCMC.R:357`).

> Caveat (4.3): a rigorous non-invariance proof (or a partner-set /
> normaliser invariance theorem rescuing the move) is open, as in HTM
> §6. A 5-tip exact stationary computation would settle it in an
> afternoon. Lane: numerical-auditor.

### 4.4 The corrected kernels — what "exact" looks like here

- **`weighted_spr`** (`src/mcmc.cpp:2886–3145`): candidate and reverse
  proposal families are both defined on the *shared prune-residual
  tree with conserved section lengths* ($\ell_{\mathrm{merge}}$
  conserved by merging, $\ell_r$ conserved by fraction split), so the
  forward and reverse normalisers are sums over the **same** set of
  (topology, bin) evaluations and cancel exactly; the remaining
  Hastings is the within-bin Beta density ratio times the bin-weight
  ratio (`:3116–3119`), and a full MH accept with prior ratio follows
  (`:3129–3143`). This is the template `gibbs_spr` and
  `gibbs_subtree_swap` fail to follow. Exact (upgrades HTM §7's
  "informal" note: the cancellation argument is closed by the
  conserved-section observation).
- **`weighted_branch_scale`** (`:2596–2697`): same-pair forward and
  reverse bin weights are identical evaluations (pair total conserved,
  other edges unchanged), Hastings `:2688–2691`, generic MH accept.
  Exact.
- **`block_gibbs_branch`** (`:2718–2863`): random-permutation
  composition of per-pair MH kernels, each accepted with
  $\beta \Delta LL + \log H$ (`:2839–2852`). Composition of invariant
  kernels with state-independent scan order ⇒ invariant. **Assumption
  made explicit:** the sweep skips prior recomputation on the grounds
  that the `rel_br_lengths` prior is flat Dirichlet(1,…,1)
  (`:2712–2714`, acceptance at `:2839` has no prior term). Exact today;
  silently wrong if the branch prior ever changes. This belongs in any
  future prior-change checklist.
- **`gibbs_kPrime`** (`:4520–4634`): random-scan Gibbs from per-character
  conditionals with weights $\beta\,LL + \log P(k'|p)$ — tempered
  likelihood, untempered prior, matching $\pi_\beta \propto L^\beta \cdot
  \mathrm{prior}$ (`src/mcmc.cpp:4253`). Exact up to the candidate-list
  truncation at relative weight $e^{-25}$
  (`kKprimeLogCutoff`, `src/mcmc_state.h:404`, applied `:4266–4272`)
  — mass ~$10^{-11}$, negligible; and the TRUNC-001 cap for the plain
  geometric arm (`:4586–4596`), previously proven. Invariant (not
  reversible — sweep — which only matters for §1.2's spectral argument,
  not for correctness).
- **`gibbs_p_marginal`**: exact by the standing proof
  (`dev/red-team/proofs/marginal-k-gibbs-p.md`, referenced at
  `src/mcmc.cpp:5197–5204`).
- **`mh_logit_p`** (`:5176–5195`): logit-scale Bactrian with Jacobian
  $\log p'(1-p') - \log p(1-p)$ — correct.
- **`joint_tl_*`** (`:5117–5137, 5384–5390`): multiplicative 2D with
  $\log H = \log m_1 + \log m_2$ — correct.
- **`nni`/`spr`/`tbr`/`pspr`**: not re-derived here; HTM §§1–3, 9
  proved them watertight and the marginal-k campaign re-verified their
  evaluator coherence. Carried as an assumption of this document.

### 4.5 `weighted_subtree_swap` — approximate Hastings, confirms HTM §8

`src/mcmc.cpp:3163–3349`. Forward: self as a *point* weight
$e^{\beta LL(x)}$ (`:3245`), candidates bin-marginalised. The MH
correction (`:3320–3323`) uses `wOrig` — the likelihood at $x$'s
*exact* fraction — where the reverse proposal density requires the
reverse candidate's **bin-midpoint** weight, and it implicitly sets
$Z(y) = Z(x)$, which fails because the other partners' totals
$br_A + br_{B'}$ change when $br_A$ changes (unlike `weighted_spr`,
there is no shared residual structure to force cancellation). Two
approximation errors, both bounded by within-bin likelihood variation
and neighbourhood asymmetry, neither zero. HTM §8 flagged the missing
partner-set correction; the two diagnoses are the same object (their
$|\mathrm{partners}|$ term is the uniform-weight shadow of the $Z$
ratio). **Not exactly invariant; default `FALSE`**
(`R/MkPrimeMCMC.R:363`), so not in default schedules — but it is in the
plan's topology block table and must not be promoted by S3 while
uncorrected.

### 4.6 Slice samplers — capped stepping-out breaks exactness (new finding)

Both `slice_scalar_impl` (`src/mcmc.cpp:3447–3565`) and
`slice_kprime_hyper_impl` (`:3596–3682`) step out with a
**deterministic cap of `maxSteps = 10` on each side**
(`:3466–3478`, `:3653–3662`).

**Requirement (derived).** The shrinkage acceptance "any $x_1 \in I$
with $f(x_1) \ge z$" yields detailed balance only if the interval
selection satisfies $P(I \mid x_0, z) = P(I \mid x_1, z)$ for every
such $x_1$ (pair the forward and reverse shrinkage sequences, which
have equal length; the interval-selection densities must then match).
With unbounded stepping out, $I$ depends only on the grid alignment
$u$ and the slice's connected component, giving equality. With the cap
binding on (say) the left, $I = [x_0 - uw - 10w,\; \cdot\,]$: to
reproduce the same $I$ from $x_1$ requires $u' = u + (x_1 - x_0)/w \in
(0,1)$, impossible whenever $|x_1 - x_0| \ge w$ — yet such $x_1$ lie in
$I$ (width up to $21w$) and are accepted by shrinkage. So
$P(I|x_1) = 0 < P(I|x_0)$: detailed balance fails **on cap-binding
events**. The standard repair is to randomise the apportionment of a
*total* step budget between the two sides (Neal 2003 §4.1's
stepping-out procedure does exactly this so that the capped variant
remains exact).

**Scope.** Post-adaptation the widths target ~3 expansions per call
(`.AdaptSliceWidths`, `R/RunMkPrime.R:4884–4908`, `target = 3`), so a
10-per-side cap rarely binds; but widths freeze at the Warmup boundary
(`:1198–1201`) and heavy-tailed conditionals (e.g. `rate_log_sd` under
weak data, the BG hyper-ridge that motivated the (s,r)
reparameterisation, `src/mcmc.cpp:3567–3588`) can bind it in the
Sample phase. The 100-iteration shrinkage cap with
restore-and-return-false (`:3482–3564`) is *not* a problem (forward and
reverse shrinkage sequences pair at equal length, and failure → hold is
invariant). The fix is ~6 lines but changes the RNG stream, so per the
trivial-fix policy (behavioural, test-pinned) it is **reported, not
patched**.

> Caveat (4.6): quantify cap-binding frequency. One counter
> (`nExp == 2*maxSteps` per call) added in S1 settles whether this is
> ever exercised in practice. Lane: numerical-auditor.

### 4.7 M-159 cache boost: state-dependent weights void the algebra (new finding)

**Mechanism.** When the node-CL cache is valid, move selection uses
boosted cumulative weights: types {4, 5, 6, 23, 24} × `cacheBonus`
(`src/mcmc.cpp:6038–6053`), switched **per chain per iteration** on
`states[ch]->nodeCL.ready()` (`:6137–6144`). `cacheBonus` defaults to
**5** (`R/MkPrimeMCMC.R:372`) and is passed in every phase, Sample
included (`R/RunMkPrime.R:1081`). So the realised kernel is
$K(x,\cdot) = \sum_i w_i(C(x))\,K_i(x,\cdot)$ with $C$ = cache flag —
**the plan §5 identity does not apply**, because the weights are not
constants.

**When state-dependent mixing is still safe.** If the auxiliary $C$
evolved *autonomously* (its next value depending only on its current
value and independent randomness), then $\pi \otimes \mu_C$ would be
stationary for the augmented chain and the $X$-marginal would remain
exactly $\pi$: given $C=c$, the applied kernel $\sum_i w_i(c) K_i$
preserves $\pi$ for each fixed $c$. MkPrime's $C$ is **not**
autonomous: readiness is created by evaluations and destroyed by
*accepted* mutations (`nodeCL.invalidate_all()` on accept:
`src/mcmc.cpp:1939, 2860, 3141`; granular invalidation in the slice
accept branch `:3537–3555`), and acceptance is correlated with the
state. Conditioning on "an accepting move just ran" tilts $X$ away from
$\pi$, and the boost then applies different weights to that tilted
conditional.

**Counterexample (verified exactly).** $X \in \{0,1\}$,
$\pi = (2/3, 1/3)$. $K_A$ = Metropolis flip (invariant), $K_B$ =
identity (invariant). Auxiliary $C' = 1$ iff move A ran **and
accepted** (the invalidate-on-accept pattern). Weights
$P(\text{choose }A \mid C{=}1) = 0.9$, $P(\cdot \mid C{=}0) = 0.1$. The
augmented 4-state chain's stationary $X$-marginal is
$(0.857, 0.143)$ — a **+0.19 bias** on $P(X{=}0)$ against
$\pi = (0.667, 0.333)$. Control: the same construction with constant
weights returns $\pi$ exactly. Both kernels are individually
π-invariant throughout; the bias is produced entirely by
acceptance-coupled state-dependent selection.

**What this means here.**

1. **Pre-existing correctness gap, independent of this plan.** The plan
   §8 already lists M-159 as an open interaction, but as a
   measurement question; it is also an invariance question. The plan
   §5 argument must either (a) be conditioned on `cacheBonus = 1`,
   (b) restore state-independence (e.g. decide the boost from an
   autonomous schedule — iteration parity, or a per-batch frozen flag —
   rather than live cache state), or (c) add the selection-probability
   ratio $w_{m}(x')/w_{m}(x)$ to each MH acceptance (impossible for the
   always-accept kernels, so (b) is the realistic fix). Magnitude in
   production is unquantified — the toy shows the mechanism is real,
   not that it is large; cache-flag/state correlation in a 30-move
   mixture is plausibly weak.
2. **Measurement bias for S1/S3 (the brief's question).** Yes, the
   measurement is biased in a specific, nameable way: boosted types are
   preferentially *sampled* when the cache is ready, which is exactly
   when they are cheapest (that is M-159's purpose), so their measured
   `moveTimeNs` per draw — the ESJD/s denominator — is the
   boost-conditional cost, not the schedule-average cost. Comparing
   boosted types (4, 5, 6, 23, 24) against unboosted blockmates
   (`tbr` = 17, `pspr` = 20, weighted = 12/13/14, `block_gibbs` = 15)
   inside the topology and branch blocks conflates kernel quality with
   cache regime; and S3 would then multiply the resulting weights by
   the boost *again* at runtime — double-counting. **Fix: S1 keeps two
   ledgers per (chain, move) — jump and time accumulated separately
   under `useBoost` true/false** (the flag is already computed at
   `src/mcmc.cpp:6137`). This costs one extra matrix pair and makes
   the bias measurable instead of arguable.
3. Reading state for jump measurement cannot *create* an M-159
   interaction as long as it performs no evaluations (does not touch
   `nodeCL`) and draws no RNG. The plan's bit-identical-chains gate
   (plan §5 item 1) catches both failure modes, because move selection
   depends on cache state: perturbing the cache would change the next
   draw and desynchronise the stream. Endorsed as the primary gate.

> Caveat (4.7): quantify the M-159 stationary bias. Cheapest decisive
> measurement: a small dataset with an exactly enumerable posterior (or
> the existing fast SBC) run at `cacheBonus = 1` vs `5`, same seeds.
> Lane: numerical-auditor / mcmc-diagnostician.

### 4.8 The amplification interaction — why this gates S3

These defects are not independent of the scheduling redesign. The
current score treats `gibbs_spr` as an ordinary move (acceptance ≈ 1
counts `return true` commits); the assessment (§3) itself notes
always-accept kernels currently score as `dim/cost`. **ESJD/s will
credit `gibbs_spr` and `gibbs_subtree_swap` with large, frequent RF
displacements — they are always-movers by construction — and will
plausibly *increase* their weight share within the topology block.** A
scheduler upgrade that re-allocates budget toward the two kernels with
open or proven invariance defects converts a background bias into a
larger one. Hence the gate condition in §0: before S3 flips any weight
rule, `gibbs_spr` must be fixed (fraction-draw + MH, §4.2) or removed
from the schedule, and `gibbs_subtree_swap` must get the
normaliser-ratio acceptance or the 5-tip exactness check that clears
it. S1 (measurement only) is unaffected and can proceed.

### Verdict (Q4)

**The fixed-weight invariance algebra is sound and the freeze boundary
is real (modulo the "warmup"→"tuning" wording), but the "every
registered move is π-invariant" premise is false today:** one proven
violation (`gibbs_spr`, default-on), one presumed violation
(`gibbs_subtree_swap`, default-on, open since HTM §6), one approximate
Hastings (`weighted_subtree_swap`, default-off), one conditional
exactness failure (capped slice stepping-out, binding-rare), and one
architecture-level violation (M-159 state-dependent weights, default
`cacheBonus = 5`). "No SBC campaign for the weight rule change" remains
a valid conclusion *for the weight-rule delta in isolation*; it is not
a certificate for the sampler, and S3 must not proceed until the two
Gibbs topology kernels are fixed, excluded, or exactness-checked.

---

## 5. Q5 — Whitening from the warmup snapshot buffer

### Statement

Plan §2: per-coordinate whitening scales estimated in R from the warmup
snapshot buffer, pushed to C++ per batch, mirroring the
`chain_rhos`/`jointRhos` machinery (`.AccumulateRhoSnapshot`
`R/RunMkPrime.R:3380–3393`, `.EstimateJointRhos` `:3398–3443`,
`.BuildJointRhoMatrix` `:3356–3371`); "a coordinate with no usable
scale contributes 0".

### 5.1 The diagnostic/kernel distinction — holds, with a precise boundary

Two different objects flow through the same plumbing:

- `jointRhos` are **kernel parameters**: they shape the 2D Bactrian
  proposal (`src/mcmc.cpp:5119`). Data-dependent kernel adaptation is
  licensed during Warmup/Tuning by the freeze (§4.1) — this is ordinary
  finite-adaptation MCMC, and needs no diminishing-adaptation condition
  because adaptation ceases before any retained draw.
- Whitening scales are **score inputs only**: they enter
  `jumpSq` accumulation → ESJD/s → weights. They never touch a proposal
  density. A wrong scale can therefore only *misallocate* effort, never
  bias the target — provided two invariants hold: (i) the scales (and
  the ESJD-per-iteration internal mode, plan §4) are never fed back
  into proposal tuning; (ii) the jump computation itself draws no RNG
  and perturbs no cache (§4.7.3, gated by the bit-identical test).

With (i) and (ii), the plan's claimed distinction is confirmed. Note
the weights themselves are kernel parameters, so the scales affect the
kernel *indirectly during adaptation* — which is fine for exactly the
same freeze reason, but it means "diagnostic" should be read as "cannot
affect correctness of the sampling phase", not "cannot affect the
chain at all".

### 5.2 Failure modes at the estimator

The snapshot buffer holds one row per warmup batch (cold chain only,
`R/RunMkPrime.R:1259–1261`), rolling 500 rows (`:3391`). Three hazards:

1. **σ̂ → 0.** The plan's sentence has the failure direction backwards:
   a vanishing scale sends whitened jumps $\Delta^2/\hat\sigma^2$ to
   **+∞**, not 0. A near-degenerate coordinate (e.g. `rate_neo` pinned
   by strong data, or a `rel_br_lengths` edge stuck near the simplex
   floor with log-scale whitening amplifying it) would dominate the
   block sum and hand the block to whichever move jiggles it. Required:
   an explicit usability definition — coordinate $j$ is usable iff the
   buffer has ≥ 50 rows (mirror `:3400`) *and*
   $\hat\sigma_j > \varepsilon_{\mathrm{abs}} + \varepsilon_{\mathrm{rel}}\,|\hat\mu_j|$;
   unusable coordinates contribute 0 (conservative under the sum rule;
   under §2's bottleneck weighting they must also be excluded from the
   min, or the allocator chases an unmeasurable coordinate). Add
   winsorisation of per-draw whitened jumps (cap at ~10²) so one
   outlier cannot own an interval — the `jumpSqSq` accumulator's SE
   makes such events visible.
2. **Transient inflation.** Warmup snapshots are by definition drawn
   from the transient; monotone drift inflates $\hat\sigma$ (variance
   about the mean of a trending series), which *deflates* whitened
   jumps on drifting coordinates — the opposite bias to (1), and it
   penalises exactly the coordinates still moving. Since the score only
   ranks *within* a block, this matters when moves differ in which
   coordinates they serve. Preferred source: estimate scales from the
   post-stabilisation tail of the buffer, or from the Tuning buffer
   (`tuningBuf`, `:1150–1167`), which is post-stabilisation by
   construction; re-estimate once per Tuning round.
3. **Few rows / degenerate cases.** Below the row gate everything is
   unusable and all scores are 0: S1 must define the 0/0 → softmax
   behaviour explicitly (uniform over the block is the sane answer, and
   is what `.AdaptMoveWeights` already does for non-scoreable moves via
   the `minProposals` gate, `R/RunMkPrime.R:4456–4459`). Constant
   columns (e.g. `beta_scale` absent, `p` under `logseries`) must be
   structurally excluded, not estimated.

One more structural point: `rel_br_lengths` is a simplex, so its ~105
log-scale coordinates carry one linear constraint and the diagonal
whitening ignores that dependence. For a *diagnostic* this is
acceptable (it mis-states total displacement by a bounded factor); an
ILR transform would be exact and is overkill for v1. Record the choice.

### Verdict (Q5)

**Legitimate as designed, with four required specifics:** ≥50-row gate;
ε-floor definition of "usable" (fixing the plan's inverted σ̂→0
wording); scales taken post-stabilisation (tail or Tuning buffer)
rather than whole-warmup; winsorised jumps plus explicit 0/0 handling.
None of these threaten correctness (the freeze protects that); all four
protect the *rankings* S2 will be trusted to adjudicate.

---

## 6. Consolidated S1 instrumentation additions

All at the single instrumentation point `src/mcmc.cpp:6163–6197`,
mirroring `moveTimeNs`; none change weights; all are covered by the
bit-identical gate.

1. `jumpSq(ch, moveIdx)`, `jumpSqSq(ch, moveIdx)` — as planned.
2. **Split by boost:** duplicate jump *and* time ledgers under the
   `useBoost` flag (`:6137`) — quantifies §4.7.2's measurement bias.
3. **Per-(move, coordinate) sums** for the branch block (`nMoves ×
   nEdge` doubles) and k′ block — enables §2's dual-vs-uniform ranking
   comparison in S2. **Branch jumps must be keyed by edge identity
   (e.g. child-node label), not row index:** every accepted topology
   commit reorders rows via `preorder_weighted_impl`
   (`src/mcmc.cpp:1926–1934, 3097–3101, 3133–3137`), so row-indexed
   deltas fabricate large spurious branch jumps; SPR additionally
   merges/splits edges, so ~3 edges per accepted SPR have
   create/destroy events that need a stated convention (recommend:
   count matched edges only; the topology block's RF already prices the
   structural change).
4. **Crossing-tail counter** per topology move: accepted-jump count and
   RF² sum for RF ≥ θ (§1.6.2).
5. **Both whitened components** for each `joint_tl_*` draw (§3.3.3).
6. **Slice cap-binding counter**: increment when a stepping-out loop
   exhausts `maxSteps` (§4.6) — piggybacks on the existing
   `sliceExpansions` plumbing (`:6171, 6177`).
7. RF and CID both, as the plan already commits (§3 of the plan).

---

## 7. References

- Geyer, C.J. (1992). Practical Markov chain Monte Carlo. *Statistical
  Science* 7:473–483. (Asymptotic variance, spectral representation.)
- Lawler, G. & Sokal, A. (1988). Bounds on the L² spectrum for Markov
  chains and Markov processes. *Trans. AMS* 309:557–580. (Conductance.)
- Neal, R.M. (2003). Slice sampling. *Annals of Statistics*
  31:705–767. (§4.1 stepping-out; randomised apportionment of the step
  budget.)
- Pasarica, C. & Gelman, A. (2010). Adaptively scaling the Metropolis
  algorithm using expected squared jumped distance. *Statistica Sinica*
  20:343–364.
- Roberts, G.O., Gelman, A. & Gilks, W.R. (1997). Weak convergence and
  optimal scaling of random walk Metropolis algorithms. *Ann. Appl.
  Prob.* 7:110–120.
- Roberts, G.O. & Rosenthal, J.S. (2001). Optimal scaling for various
  Metropolis–Hastings algorithms. *Statistical Science* 16:351–367.
- Sherlock, C. & Roberts, G.O. (2009). Optimal scaling of the random
  walk Metropolis on elliptically symmetric unimodal targets.
  *Bernoulli* 15:774–798.
- Sinclair, A. (1992). Improved bounds for mixing rates of Markov
  chains and multicommodity flow. *Combin. Probab. Comput.* 1:351–370.
- Tierney, L. (1994). Markov chains for exploring posterior
  distributions. *Annals of Statistics* 22:1701–1728. (Mixtures and
  cycles of invariant kernels.)
- Vehtari, A. et al. (2021). Rank-normalization, folding, and
  localization: An improved R̂. *Bayesian Analysis* 16:667–718.
- Zanella, G. (2020). Informed proposals for local MCMC in discrete
  spaces. *JASA* 115:852–865. (Locally balanced / normaliser
  conditions for informed always-accept-style proposals.)
- Prior lane artefacts: `dev/red-team/proofs/hastings-tree-moves.md`
  (HTM); `dev/red-team/proofs/marginal-k-gibbs-p.md`.
