# Pre-registered decision rules for the Phase-2 overnight runs

Written 2026-06-02 eve BEFORE the results land (advisor-driven), so the read is
not rationalised after seeing the numbers. Job: `bx3s4g4qz` (driver log
`T-OVL-overnight-driver.log`); logs `T-OVL-moveson.log`, `T-OVL-precheck-200k.log`.

## (1) Moves-on overlap (T-OVL-moveson.log) — read RELATIVE to the gated baseline
The four re-enabled weighted/block moves are ON in BOTH modes. They inherit the
SAME p-mixing limitation as the baseline, so do NOT score against an absolute PASS.

- PASS criterion: **no NEW failure on a parameter the gated baseline passed.**
  tree_length and rate_log_sd must stay PASS; p-means must agree (small d_mean).
- A **fresh FAIL on tree_length or rate_log_sd** is the ONLY signal that implicates
  the re-enabled moves (proposal bias) — that would block and require per-move
  isolation (MARGINAL_K_OVL_EXTRA=<one move>).
- A p-INCONCLUSIVE or borderline p-d_sd is the PRE-EXISTING artifact (see below),
  NOT a Phase-2 regression. Do not round it up to a Phase-2 failure.

## (2) p d_sd recheck (T-OVL-precheck-200k.log) — three pre-registered outcomes
Context: the RB identity (test-marginal-k-rb-free-topology.R) ALREADY proves the
marginal posterior on p is identical between modes. So sd(p) must agree in the
limit. marginal_k integrates out k', so there is NO conjugate Gibbs-p — p moves
ONLY via mh_logit_p (MH). That is WHY marginal-p-ESS (157-256) << sampled_k (1259).
This underpowering is INTRINSIC, not a bug.

Check R-hat / the p-trace, not just ESS (slow mixing masquerades as a spread diff).

- (a) **d_sd shrinks toward noise as marginal-p-ESS rises** → artifact CONFIRMED.
      Refine the test (raise the sd-test ESS floor and/or a multiple-comparison-aware
      bar) and document. Baseline correctness stands.
- (b) **d_sd GROWS with ESS** → a real mh_logit_p target/spread bias under marginal_k.
      This is a SHIPPED-marginal_k bug, independent of Phase 2, and BLOCKS
      "baseline correctness nailed." Investigate mh_logit_p.
- (c) **d_sd stays ~3 AND marginal-p-ESS stays low (<~500 even at 200k)** →
      INCONCLUSIVE BY UNDERPOWERING, *not* a pass. Fall back to the RB identity
      (which already settles p-target-equivalence) and STOP chasing p via MCMC
      overlap. Do NOT report this as PASS.

## Operational
`bx3s4g4qz` holds the .so DLL while running. A no-source-change `pkgload::load_all`
is lock-safe (read-only DLL load), but RECOMPILING while the job holds the DLL is a
Windows lock race. The binary is CURRENT as of commit c68f3dc. Any later wakeup must
verify `bx3s4g4qz` finished ("OVERNIGHT DONE" in the driver log) before recompiling,
and skip the rebuild entirely if `src/` is unchanged.
