# sim3v4-params.R -------------------------------------------------------------
#
# Named parameter sets for the redesigned Sim 3 ("v4") convergent-ecology
# study. Three regimes spanning conservative -> balanced -> strong
# convergent signal.
#
# Signal heuristics (binary character, baseRate = 1, phi = 4):
#   Transformational (Mk symmetric k=2):
#     P(change | t) = (1 - exp(-2 r t)) / 2
#     eco-1 branch (r = 4): P(change | t=0.05) ~ 0.165,
#                           P(change | t=0.10) ~ 0.275,
#                           P(change | t=0.20) ~ 0.399 (near saturation 0.5).
#   Neomorphic (asym CTMC, rateLoss = 1 -> base rates 1, 1):
#     Under z = 1: rates become (phi, 1/phi) = (4, 0.25). Per-edge change
#     probability on an eco-1 branch dominated by the 4x gain rate; rough
#     change prob is comparable in scale to the trans arm but biased toward
#     gain, producing the asymmetric synapomorphy signal that drives
#     genuine convergence between the two ecology-1 cherries (A2,A1) and
#     (C2,C1).
#
# Sub-saturation rule of thumb: keep stemBrEco below ~0.10 so that the
# expected number of changes per ecology-affected character on the eco
# stem stays below ~0.4. Above that, the spurious signal is washed out
# by saturation and BOTH chains fail to distinguish anything.
#
# The convergent signal that pulls a blind chain to the spurious clade
# {A1,A2,C1,C2} accumulates on TWO branches per ecology-affected
# character: the eco stem inside clade A and the eco stem inside clade C.
# Roughly 2 * P(change | stemBrEco) shared parallel changes per
# ecology-affected character.
#
# The ancestry signal that keeps clade A intact accumulates on the clade
# A stem (stemBrClade) and the within-clade pectinate branches (tipBr).
# The signal supporting (A,C) as true sisters accumulates on rootBr.

SIM3V4_PARAMS <- list(

  # --- v4a: conservative -----------------------------------------------------
  # Short eco stem (weak convergent signal), long clade stem & root
  # (strong ancestry signal). Expect: both blind and aware recover the
  # true topology; this is the "easy" regime to verify the harness works.
  #
  # Per ecology-affected binary character on each eco stem:
  #   trans: P(change | 0.05) ~ 0.165 ; ~0.33 expected parallel changes
  #          across both eco stems combined.
  # Per non-eco character on clade-A stem (r=1, t=0.30):
  #   P(change) = (1 - exp(-0.6))/2 ~ 0.226 ; strong ancestry signal.
  v4a = list(
    name        = "v4a_conservative",
    tipBr       = 0.10,
    stemBrEco   = 0.05,
    stemBrClade = 0.30,
    rootBr      = 0.15,
    nNeo        = 100L,
    nTrans      = 200L,
    phi         = 4,
    pi0         = 0.75,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1,
    kStates     = 2L
  ),

  # --- v4b: balanced (target sweet spot) -------------------------------------
  # Moderate eco stem with strong ecology effect; moderate clade stem
  # and root. The convergent signal is real but the ancestry signal is
  # still substantial. Expect: blind chain pulled toward the spurious
  # {A1,A2,C1,C2} clade much of the time; aware chain recovers truth.
  #
  # Per ecology-affected binary character on each eco stem:
  #   trans: P(change | 0.08) ~ 0.247 ; ~0.49 parallel changes across both.
  # Per non-eco character on clade stem (r=1, t=0.15):
  #   P(change) = (1 - exp(-0.3))/2 ~ 0.130 ; moderate ancestry signal.
  v4b = list(
    name        = "v4b_balanced",
    tipBr       = 0.05,
    stemBrEco   = 0.08,
    stemBrClade = 0.15,
    rootBr      = 0.08,
    nNeo        = 100L,
    nTrans      = 200L,
    phi         = 4,
    pi0         = 0.75,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1,
    kStates     = 2L
  ),

  # --- v4c: strong convergent signal (stress regime) -------------------------
  # Long eco stem (near saturation: P(change|0.15) ~ 0.40 per char),
  # short clade stem and root. The blind chain should be strongly drawn
  # to the spurious clade; aware chain must downweight ecology-affected
  # characters aggressively to recover truth. Analogous in spirit to the
  # v3 challenge regime but with only 4 ecology-1 tips instead of 8.
  #
  # Per ecology-affected binary character on each eco stem:
  #   trans: P(change | 0.15) ~ 0.398 ; ~0.80 parallel changes across both.
  # Per non-eco character on clade stem (r=1, t=0.05):
  #   P(change) = (1 - exp(-0.1))/2 ~ 0.048 ; weak ancestry signal.
  v4c = list(
    name        = "v4c_strong",
    tipBr       = 0.03,
    stemBrEco   = 0.15,
    stemBrClade = 0.05,
    rootBr      = 0.04,
    nNeo        = 100L,
    nTrans      = 200L,
    phi         = 4,
    pi0         = 0.75,
    theta       = 1.0,
    baseRate    = 1.0,
    rateLoss    = 1,
    kStates     = 2L
  )
)
