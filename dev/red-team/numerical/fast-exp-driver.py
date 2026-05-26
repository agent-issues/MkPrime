#!/usr/bin/env python3
"""Lane N1 - fast_exp + closed-form P(t) stability stress driver.

Audits src/fast_exp.h::fast_neg_exp and the JC/MkN/F81 closed-form
transition-probability formulas used in src/rate_matrix.cpp,
src/likelihood.cpp, src/gibbs_partial_cl.h, src/mcmc_likelihood.cpp.

There is NO eigendecomposition in this codebase; all P(t) is analytic.
So the audit reduces to:
  (1) measure fast_neg_exp vs std::exp (math.exp) vs mpmath.exp;
  (2) measure the assembled P(t) entries vs mpmath reference;
  (3) flag any regime where relative error exceeds 1e-10.

Reference: mpmath at dps=50.

Outputs:
  dev/red-team/numerical/fast-exp-results/scalar_exp.csv
  dev/red-team/numerical/fast-exp-results/jc.csv
  dev/red-team/numerical/fast-exp-results/mkn.csv
  dev/red-team/numerical/fast-exp-results/f81.csv
  dev/red-team/numerical/fast-exp-results/summary.txt
"""

from __future__ import annotations
import argparse
import csv
import math
import os
import struct
import sys
import time
from pathlib import Path

import mpmath as mp

mp.mp.dps = 50

OUT_DIR = Path(__file__).resolve().parent / "fast-exp-results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Bit-exact Python replica of src/fast_exp.h::fast_neg_exp.
#
# We replicate the C function so that:
#   - we can probe it from Python (no compiler needed);
#   - we can confirm it matches std::exp within its claimed bound;
#   - we can feed it into the JC/MkN/F81 closed forms identically to the
#     production code.
#
# All ops are pure IEEE-754 double; reinterpretation uses struct.pack.
# -----------------------------------------------------------------------------
def fast_neg_exp(x: float) -> float:
    # Argument is always <= 0 in our use case.
    if x < -708.0:
        return 0.0

    # Argument reduction: x = n * ln(2) + r, |r| <= ln(2)/2.
    LOG2E = 1.4426950408889634074
    LN2_HI = 6.93147180369123816490e-01
    LN2_LO = 1.90821492927058500170e-10

    nd = x * LOG2E
    # The C code: int n = static_cast<int>(nd - 0.5).
    # For x <= 0, nd <= 0, so nd - 0.5 <= -0.5, and C's truncation toward zero
    # gives floor for negative arguments-ish.  Python's int() also truncates
    # toward zero, so int(nd - 0.5) matches.
    n = int(nd - 0.5)

    r = (x - n * LN2_HI) - n * LN2_LO

    # Degree-11 Taylor polynomial, Horner form
    c2 = 1.0 / 2.0
    c3 = 1.0 / 6.0
    c4 = 1.0 / 24.0
    c5 = 1.0 / 120.0
    c6 = 1.0 / 720.0
    c7 = 1.0 / 5040.0
    c8 = 1.0 / 40320.0
    c9 = 1.0 / 362880.0
    c10 = 1.0 / 3628800.0
    c11 = 1.0 / 39916800.0
    p = 1.0 + r * (1.0 + r * (c2 + r * (c3 + r * (c4 +
        r * (c5 + r * (c6 + r * (c7 + r * (c8 +
        r * (c9 + r * (c10 + r * c11))))))))))

    # exponent injection: 2^n via raw IEEE-754 bit pattern
    bits = (n + 1023) << 52
    scale = struct.unpack("<d", struct.pack("<q", bits))[0]
    return p * scale


# -----------------------------------------------------------------------------
# Closed-form P(t) builders, evaluated in (a) fast_neg_exp, (b) math.exp,
# (c) mpmath (reference).
# -----------------------------------------------------------------------------
def jc_P(k: int, t: float, exp_fn):
    """JC(k): src/rate_matrix.cpp::jc_transition_probs, src/likelihood.cpp line 84."""
    alpha = k / (k - 1.0)
    E = exp_fn(-alpha * t)
    inv_k = 1.0 / k
    p_same = inv_k + (1.0 - inv_k) * E
    p_diff = inv_k - inv_k * E
    return p_same, p_diff


def jc_P_ref(k: int, t: float):
    k_ = mp.mpf(k)
    t_ = mp.mpf(t)
    alpha = k_ / (k_ - 1)
    E = mp.exp(-alpha * t_)
    inv_k = 1 / k_
    p_same = inv_k + (1 - inv_k) * E
    p_diff = inv_k - inv_k * E
    return p_same, p_diff


def mkn_P(rate_loss: float, t: float, exp_fn):
    """MkN 2-state: src/rate_matrix.cpp::mkn_transition_probs."""
    sum_rl = 1.0 + rate_loss
    rate01 = 2.0 / sum_rl
    rate10 = 2.0 * rate_loss / sum_rl
    lam = rate01 + rate10  # always 2.0
    E = exp_fn(-lam * t)
    i01 = rate01 / lam
    i10 = rate10 / lam
    P00 = i10 + i01 * E
    P01 = i01 - i01 * E
    P10 = i10 - i10 * E
    P11 = i01 + i10 * E
    return P00, P01, P10, P11


def mkn_P_ref(rate_loss: float, t: float):
    rl = mp.mpf(rate_loss)
    t_ = mp.mpf(t)
    sum_rl = 1 + rl
    rate01 = 2 / sum_rl
    rate10 = 2 * rl / sum_rl
    lam = rate01 + rate10
    E = mp.exp(-lam * t_)
    i01 = rate01 / lam
    i10 = rate10 / lam
    return (i10 + i01 * E, i01 - i01 * E, i10 - i10 * E, i01 + i10 * E)


def f81_P(pi: list[float], t: float, exp_fn):
    """F81 k-state: src/gibbs_partial_cl.h::f81_transition, src/mcmc_likelihood.cpp:1867.

    P_ij(t) = pi_j * (1 - d) + delta_ij * d,    d = exp(-mu*t)
    mu = 1 / (1 - sum pi_j^2)
    Returns the full k x k matrix as a flat tuple of (i,j,P_ij).
    """
    k = len(pi)
    sumPiSq = sum(p * p for p in pi)
    if sumPiSq >= 1.0:
        # impossible for valid pi but guard anyway
        mu = float("inf")
    else:
        mu = 1.0 / (1.0 - sumPiSq)
    arg = -mu * t
    d = exp_fn(arg) if arg > -708.0 else 0.0
    one_m_d = 1.0 - d
    out = []
    for i in range(k):
        for j in range(k):
            Pij = pi[j] * one_m_d + (d if i == j else 0.0)
            out.append(Pij)
    return out, mu, d


def f81_P_ref(pi: list[float], t: float):
    k = len(pi)
    pi_ = [mp.mpf(p) for p in pi]
    t_ = mp.mpf(t)
    sumPiSq = sum(p * p for p in pi_)
    mu = 1 / (1 - sumPiSq)
    d = mp.exp(-mu * t_)
    one_m_d = 1 - d
    out = []
    for i in range(k):
        for j in range(k):
            Pij = pi_[j] * one_m_d + (d if i == j else mp.mpf(0))
            out.append(Pij)
    return out, mu, d


def relerr(approx: float, ref) -> float:
    """Relative error |approx - ref| / |ref|, with mpmath-safe ref."""
    if ref == 0:
        return 0.0 if approx == 0 else float("inf")
    return float(abs(mp.mpf(approx) - ref) / abs(ref))


def make_acrv_rates(K: int, sigma: float):
    """K rate categories for lognormal ACRV with given log-sd, normalised to mean 1.

    Mirrors the discretisation used implicitly elsewhere: quantile midpoints
    of LogNormal(mu, sigma^2), with mu chosen so E[r] = 1.
    """
    mu_ = -sigma * sigma / 2.0
    qs = [(i + 0.5) / K for i in range(K)]
    # inverse normal CDF via mpmath
    rates = []
    for q in qs:
        z = float(mp.sqrt(2) * mp.erfinv(2 * q - 1))
        r = math.exp(mu_ + sigma * z)
        rates.append(r)
    # renormalise to mean 1
    m = sum(rates) / K
    return [r / m for r in rates]


# -----------------------------------------------------------------------------
# Test grids
# -----------------------------------------------------------------------------
def grid_full():
    return {
        "k": [2, 4, 6, 10],
        "rt": [1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1.0, 10.0, 100.0, 1e4],
        "pi_asym_logratio": [0.0, math.log10(3), 2.0, 6.0],  # uniform, mild, strong, pathological
        "sigma_acrv": [0.1, 1.0, 3.0, 10.0],
        "K_acrv": [4],
        "rate_loss": [1.0, 0.1, 0.01, 100.0, 1e6],  # MkN asymmetry
    }


def grid_quick():
    return {
        "k": [2, 4, 6],
        "rt": [1e-12, 1e-6, 1.0, 100.0],
        "pi_asym_logratio": [0.0, 2.0],
        "sigma_acrv": [1.0, 3.0],
        "K_acrv": [4],
        "rate_loss": [1.0, 100.0],
    }


def make_pi(k: int, asym_logratio: float):
    """Construct pi vector with given log10(max/min) asymmetry.

    Geometric spacing: pi_i = c * r^i, sum to 1, with r = 10^(-asym_logratio / (k-1)).
    """
    if asym_logratio == 0.0:
        return [1.0 / k] * k
    r = 10.0 ** (-asym_logratio / (k - 1))
    raw = [r ** i for i in range(k)]
    s = sum(raw)
    return [x / s for x in raw]


# -----------------------------------------------------------------------------
# Test 1: scalar fast_neg_exp vs std::exp vs mpmath
# -----------------------------------------------------------------------------
def test_scalar_exp(quick: bool):
    rows = [("x", "fast_neg_exp", "math_exp", "mpmath_exp", "rel_fast_vs_mp", "rel_std_vs_mp")]
    xs = [-1e-15, -1e-12, -1e-10, -1e-8, -1e-6, -1e-4, -1e-2, -0.1, -0.347, -0.5,
          -1.0, -2.0, -5.0, -10.0, -50.0, -100.0, -300.0, -500.0, -700.0,
          -707.5, -708.0]
    # boundary tests
    xs += [-708.5, -709.0]
    for x in xs:
        fe = fast_neg_exp(x)
        se = math.exp(x) if x > -745 else 0.0
        ref = mp.exp(mp.mpf(x))
        rf = relerr(fe, ref)
        rs = relerr(se, ref)
        rows.append((x, fe, se, float(ref), rf, rs))

    # random sweep
    import random
    random.seed(42)
    N = 2_000 if quick else 50_000
    for _ in range(N):
        # log-uniform over [-708, -1e-12]
        exp10 = random.uniform(-12, math.log10(708))
        x = -(10.0 ** exp10)
        fe = fast_neg_exp(x)
        se = math.exp(x)
        ref = mp.exp(mp.mpf(x))
        rows.append((x, fe, se, float(ref), relerr(fe, ref), relerr(se, ref)))

    fp = OUT_DIR / "scalar_exp.csv"
    with fp.open("w", newline="") as f:
        csv.writer(f).writerows(rows)

    # For the headline: exclude the underflow boundary (x < -708) where
    # fast_neg_exp deliberately returns 0 by design (header guard). Those
    # values are <2.2e-308 and are also subnormal-territory for std::exp.
    rels_fast_normal = [r[4] for r in rows[1:]
                        if isinstance(r[4], float) and r[0] >= -708.0]
    rels_std_normal = [r[5] for r in rows[1:]
                       if isinstance(r[5], float) and r[0] >= -708.0]
    # Also report what happens below -708 (informational).
    below = [r for r in rows[1:] if isinstance(r[0], float) and r[0] < -708.0]
    return {
        "max_fast": max(rels_fast_normal),
        "max_std": max(rels_std_normal),
        "n": len(rels_fast_normal),
        "below_708": len(below),
        "fp": fp,
    }


# -----------------------------------------------------------------------------
# Test 2: JC(k) closed form vs mpmath
# -----------------------------------------------------------------------------
def test_jc(grid):
    rows = [("k", "rt", "p_same_fast", "p_diff_fast",
             "p_same_std", "p_diff_std",
             "rel_same_fast", "rel_diff_fast",
             "rel_same_std", "rel_diff_std",
             "p_same_ref", "p_diff_ref")]
    max_rel_diff_fast = 0.0
    worst = None
    for k in grid["k"]:
        for rt in grid["rt"]:
            ps_f, pd_f = jc_P(k, rt, fast_neg_exp)
            ps_s, pd_s = jc_P(k, rt, math.exp)
            ps_r, pd_r = jc_P_ref(k, rt)
            rs_f = relerr(ps_f, ps_r)
            rd_f = relerr(pd_f, pd_r)
            rs_s = relerr(ps_s, ps_r)
            rd_s = relerr(pd_s, pd_r)
            rows.append((k, rt, ps_f, pd_f, ps_s, pd_s, rs_f, rd_f, rs_s, rd_s,
                         float(ps_r), float(pd_r)))
            if rd_f > max_rel_diff_fast:
                max_rel_diff_fast = rd_f
                worst = (k, rt, rd_f, rd_s)
    fp = OUT_DIR / "jc.csv"
    with fp.open("w", newline="") as f:
        csv.writer(f).writerows(rows)
    return {"max_rel_diff_fast": max_rel_diff_fast, "worst": worst, "fp": fp}


# -----------------------------------------------------------------------------
# Test 3: MkN closed form vs mpmath (asymmetric 2-state)
# -----------------------------------------------------------------------------
def test_mkn(grid):
    rows = [("rate_loss", "t", "P00_f", "P01_f", "P10_f", "P11_f",
             "max_rel_fast", "max_rel_std")]
    max_rel_fast = 0.0
    worst = None
    for rl in grid["rate_loss"]:
        for t in grid["rt"]:
            P_f = mkn_P(rl, t, fast_neg_exp)
            P_s = mkn_P(rl, t, math.exp)
            P_r = mkn_P_ref(rl, t)
            rels_f = [relerr(P_f[i], P_r[i]) for i in range(4)]
            rels_s = [relerr(P_s[i], P_r[i]) for i in range(4)]
            mf = max(rels_f); ms = max(rels_s)
            rows.append((rl, t, *P_f, mf, ms))
            if mf > max_rel_fast:
                max_rel_fast = mf; worst = (rl, t, mf)
    fp = OUT_DIR / "mkn.csv"
    with fp.open("w", newline="") as f:
        csv.writer(f).writerows(rows)
    return {"max_rel_fast": max_rel_fast, "worst": worst, "fp": fp}


# -----------------------------------------------------------------------------
# Test 4: F81 closed form vs mpmath (k-state asymmetric, ACRV-multiplied)
# -----------------------------------------------------------------------------
def test_f81(grid):
    rows = [("k", "pi_log_ratio", "sigma_acrv", "cat", "rate_cat",
             "t", "mu", "d_fast", "d_std", "d_ref",
             "max_rel_P_fast", "max_rel_P_std")]
    max_rel_fast = 0.0
    worst = None
    base_ts = [1e-10, 1e-6, 1e-2, 1.0, 100.0]
    for k in grid["k"]:
        if k < 2:
            continue
        for asym in grid["pi_asym_logratio"]:
            pi = make_pi(k, asym)
            for sigma in grid["sigma_acrv"]:
                K = grid["K_acrv"][0]
                rates = make_acrv_rates(K, sigma)
                for cat_i, rate_cat in enumerate(rates):
                    for base_t in base_ts:
                        t = base_t * rate_cat
                        P_f, mu_f, d_f = f81_P(pi, t, fast_neg_exp)
                        P_s, mu_s, d_s = f81_P(pi, t, math.exp)
                        P_r, mu_r, d_r = f81_P_ref(pi, t)
                        rels_f = [relerr(P_f[i], P_r[i]) for i in range(len(P_f))]
                        rels_s = [relerr(P_s[i], P_r[i]) for i in range(len(P_s))]
                        # Skip Pij entries where ref is 0 and approx is 0 (vacuous)
                        rels_f = [r for r in rels_f if math.isfinite(r)]
                        rels_s = [r for r in rels_s if math.isfinite(r)]
                        mf = max(rels_f) if rels_f else 0.0
                        ms = max(rels_s) if rels_s else 0.0
                        rows.append((k, asym, sigma, cat_i, rate_cat, t,
                                     mu_f, d_f, d_s, float(d_r), mf, ms))
                        if mf > max_rel_fast:
                            max_rel_fast = mf
                            worst = (k, asym, sigma, cat_i, t, mu_f, mf)
    fp = OUT_DIR / "f81.csv"
    with fp.open("w", newline="") as f:
        csv.writer(f).writerows(rows)
    return {"max_rel_fast": max_rel_fast, "worst": worst, "fp": fp}


# -----------------------------------------------------------------------------
# Test 5: Specific cancellation check — p_diff for JC at extreme rt.
# This is the headline finding: 1 - E for E ~ 1 catastrophically cancels.
# We also check the candidate fix (using -expm1(arg)) for comparison.
# -----------------------------------------------------------------------------
def test_cancellation(quick: bool):
    rows = [("k", "rt", "p_diff_curr", "p_diff_expm1",
             "rel_curr", "rel_expm1", "abs_logL_per_site_curr")]
    for k in [2, 4, 6, 10]:
        for log10_rt in range(-16, 1):  # rt = 1e-16 ... 1
            rt = 10.0 ** log10_rt
            alpha = k / (k - 1.0)
            arg = -alpha * rt
            # current implementation
            E_curr = math.exp(arg)
            pd_curr = (1.0 / k) - (1.0 / k) * E_curr
            # candidate fix: use expm1
            # 1 - E = -(E - 1) = -expm1(arg); pd_curr = (1/k)*(-expm1(arg))
            pd_fix = (1.0 / k) * (-math.expm1(arg))
            # reference
            E_r = mp.exp(mp.mpf(arg))
            pd_r = mp.mpf(1) / k - (mp.mpf(1) / k) * E_r
            rc = relerr(pd_curr, pd_r)
            rf = relerr(pd_fix, pd_r)
            # per-site logL impact if observation forces p_diff (worst case):
            # site_logL = log(p_diff); absolute logL error = rel error of p_diff
            abs_logL = rc if pd_curr > 0 else float("inf")
            rows.append((k, rt, pd_curr, pd_fix, rc, rf, abs_logL))
    fp = OUT_DIR / "cancellation_jc_pdiff.csv"
    with fp.open("w", newline="") as f:
        csv.writer(f).writerows(rows)

    # F81 cancellation check on (1 - d)
    rows2 = [("k", "asym", "sigma", "rate_cat", "t",
              "one_m_d_curr", "one_m_d_expm1",
              "rel_curr", "rel_expm1")]
    for k in [4]:
        for asym in [0.0, 2.0]:
            pi = make_pi(k, asym)
            sumPiSq = sum(p * p for p in pi)
            mu = 1.0 / (1.0 - sumPiSq)
            for log10_t in range(-16, 1):
                t = 10.0 ** log10_t
                arg = -mu * t
                d_curr = math.exp(arg)
                one_m_d_curr = 1.0 - d_curr
                one_m_d_fix = -math.expm1(arg)
                # ref
                d_r = mp.exp(mp.mpf(arg))
                one_m_d_r = 1 - d_r
                rc = relerr(one_m_d_curr, one_m_d_r)
                rf = relerr(one_m_d_fix, one_m_d_r)
                rows2.append((k, asym, "-", "-", t, one_m_d_curr,
                              one_m_d_fix, rc, rf))
    fp2 = OUT_DIR / "cancellation_f81_oneMd.csv"
    with fp2.open("w", newline="") as f:
        csv.writer(f).writerows(rows2)

    # extract the worst current vs fixed comparison
    curr_max = max(r[4] for r in rows[1:] if isinstance(r[4], float))
    fix_max = max(r[5] for r in rows[1:] if isinstance(r[5], float))
    f81_curr_max = max(r[7] for r in rows2[1:] if isinstance(r[7], float))
    f81_fix_max = max(r[8] for r in rows2[1:] if isinstance(r[8], float))

    return {
        "jc_pdiff_curr_max": curr_max,
        "jc_pdiff_fix_max": fix_max,
        "f81_oneMd_curr_max": f81_curr_max,
        "f81_oneMd_fix_max": f81_fix_max,
        "fp": fp, "fp2": fp2,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--quick", action="store_true",
                        help="Smoke-test grid; runs in <2 min.")
    args = parser.parse_args()

    t0 = time.time()
    grid = grid_quick() if args.quick else grid_full()

    print("== Lane N1: fast_exp + closed-form P(t) stability ==")
    print(f"Grid mode: {'QUICK' if args.quick else 'FULL'}")
    print(f"Reference: mpmath dps=50")
    print()

    print("Test 1: scalar fast_neg_exp vs std::exp vs mpmath")
    s1 = test_scalar_exp(args.quick)
    print(f"  max rel err (fast vs mp): {s1['max_fast']:.3e}")
    print(f"  max rel err (std  vs mp): {s1['max_std']:.3e}")
    print(f"  n samples (x >= -708): {s1['n']}  (below -708: {s1['below_708']})")

    print("Test 2: JC(k) closed-form P(t)")
    s2 = test_jc(grid)
    print(f"  max rel err on p_diff (fast): {s2['max_rel_diff_fast']:.3e}")
    print(f"  worst: {s2['worst']}")

    print("Test 3: MkN(2) closed-form P(t)")
    s3 = test_mkn(grid)
    print(f"  max rel err (fast): {s3['max_rel_fast']:.3e}")
    print(f"  worst: {s3['worst']}")

    print("Test 4: F81(k) closed-form P(t) with ACRV-multiplied rates")
    s4 = test_f81(grid)
    print(f"  max rel err (fast): {s4['max_rel_fast']:.3e}")
    print(f"  worst (k, asym, sigma, cat, t, mu, rel): {s4['worst']}")

    print("Test 5: cancellation of 1 - E (p_diff / one_m_d)")
    s5 = test_cancellation(args.quick)
    print(f"  JC p_diff  current max rel err: {s5['jc_pdiff_curr_max']:.3e}")
    print(f"  JC p_diff  expm1-fix max rel err: {s5['jc_pdiff_fix_max']:.3e}")
    print(f"  F81 1-d    current max rel err: {s5['f81_oneMd_curr_max']:.3e}")
    print(f"  F81 1-d    expm1-fix max rel err: {s5['f81_oneMd_fix_max']:.3e}")

    elapsed = time.time() - t0
    print(f"\nTotal time: {elapsed:.1f}s")

    summary = OUT_DIR / "summary.txt"
    with summary.open("w") as f:
        f.write("Lane N1 - fast_exp + closed-form P(t) stability\n")
        f.write(f"Mode: {'QUICK' if args.quick else 'FULL'}\n")
        f.write(f"Elapsed: {elapsed:.1f}s\n\n")
        f.write(f"Scalar fast_neg_exp max rel err vs mpmath: {s1['max_fast']:.3e}\n")
        f.write(f"Scalar std::exp     max rel err vs mpmath: {s1['max_std']:.3e}\n\n")
        f.write(f"JC p_diff max rel err  (fast): {s2['max_rel_diff_fast']:.3e}\n")
        f.write(f"  worst regime: {s2['worst']}\n\n")
        f.write(f"MkN max rel err (fast): {s3['max_rel_fast']:.3e}\n")
        f.write(f"  worst regime: {s3['worst']}\n\n")
        f.write(f"F81 max rel err (fast): {s4['max_rel_fast']:.3e}\n")
        f.write(f"  worst regime: {s4['worst']}\n\n")
        f.write("Cancellation study:\n")
        f.write(f"  JC p_diff current: {s5['jc_pdiff_curr_max']:.3e}\n")
        f.write(f"  JC p_diff with expm1 fix: {s5['jc_pdiff_fix_max']:.3e}\n")
        f.write(f"  F81 1-d   current: {s5['f81_oneMd_curr_max']:.3e}\n")
        f.write(f"  F81 1-d   with expm1 fix: {s5['f81_oneMd_fix_max']:.3e}\n")

    print(f"\nResults written under: {OUT_DIR}")
    print(f"Summary: {summary}")


if __name__ == "__main__":
    main()
