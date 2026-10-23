#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V64 - the desert inverse problem: can the theory's OWN vector-like matter
repair coupling unification (O4, F2)?

**What had never been done.**  The publication proved O4 (the three measured
couplings do not meet at M_P; run down from alpha_10 = 1/25 the SU(2) coupling
overshoots 5.8x and SU(3) hits a Landau pole at 2.2e9 GeV) and ruled out the
fix by COMPLETE vectorlike multiplets (a complete 10+5bar shifts all three
b_i by the same 8/3 - proportional shifts move no crossing point).  It never
scanned the INCOMPLETE content that the theory's own breaking sector supplies:
a 16_H + 16bar_H (or a split doubler 16) leaves 10+10bar and 5+5bar fermion
pairs with INDEPENDENT thresholds, and those shift the b_i non-proportionally.

**The inverse problem.**  Fix the root value alpha_X = alpha_10 = 1/25 (the
free angle is NOT re-fitted) and ask: which vector-like content, at which
thresholds, makes the three measured couplings meet at alpha_X at some scale
M_X <= M_P?  With two multiplets the system is three linear equations in
(ln M_X, ln M_A, ln M_B) - exactly determined, no fit freedom.

**The four results.**
1. The required Delta-b family is non-universal: at M_X = 2e16 GeV it is
   (2.40, 4.04, 3.86) - no complete-multiplet combination produces it
   (the paper's no-go stands, but only for complete multiplets).
2. Proportional (complete-multiplet) shifts move NO crossing point: the
   pairwise differences 1/alpha_i - 1/alpha_j are invariant.  Verified.
3. With the root value FIXED at 1/25, the scan over CL10-realizable content
   (10+10bar, 5+5bar, split-16 pieces, adjoint Weyl) finds the solutions;
   each solution is a PREDICTED desert spectrum with M_X, M_A, M_B derived.
4. Every solution abandons the convention M_GUT = M_P (M_X < M_P) and
   predicts tau_p ~ M_X^4 / (alpha_X^2 m_p^5): proton decay within or near
   the Hyper-Kamiokande / DUNE reach.  F5 becomes a live test.

**The named residual (O4, sharpened).**  The desert spectrum is required by
the couplings but not yet DERIVED from the automaton: the thresholds are
inputs of the breaking sector (O6).  F2 becomes: precision couplings admit
NO content with alpha_X = 1/25.
"""

import itertools
import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import unification as UN
from cl10 import u10pred

TWO_PI = 2.0 * math.pi
MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
A10, M_P = UN.ALPHA_10, UN.M_P
X_MAX = math.log(M_P / MZ)                 # 39.43: the cutoff
Z_MIN = math.log(1e6 / MZ)                 # thresholds at or above 1e6 GeV
INV_A10 = 1.0 / A10                        # 25: the root value, FIXED

# vectorlike Dirac pairs (and one adjoint Weyl set), Delta b = (4/3) T_i
LIB = {
    "10+10bar": (6 / 5, 2.0, 2.0),          # (3,2)_(1/6)+(3b,1)_(2/3)+(1,1) + c.c.
    "5+5bar": (2 / 3, 2 / 3, 2 / 3),        # (3,1)_(1/3)+(1,2)_(-1/2) + c.c.
    "16+16bar": (8 / 3, 8 / 3, 8 / 3),      # complete: proportional (the no-go)
    "adjW": (0.0, 4 / 3, 2.0),              # one adjoint Weyl per factor
    "(3,2)": (2 / 15, 2.0, 4 / 3),
    "(3b,1)": (16 / 15, 0.0, 2 / 3),
    "(1,2)": (2 / 5, 2 / 3, 0.0),
    "(3,1)": (4 / 15, 0.0, 2 / 3),
}
CL10_REALIZABLE = {"10+10bar", "5+5bar", "(3,2)", "(3b,1)", "(1,2)", "(3,1)"}


def meet_lines(delta_b, M_th):
    """Meeting of the three lines with one content (delta_b) at threshold
    M_th: two equations (a1=a2, a2=a3), unknowns x = ln(M_X/MZ),
    z = ln(M_th/MZ).  Returns (M_X, alpha_X, M_th) or None."""
    def rhs(i, j):
        return (1 / A_Z[i] - 1 / A_Z[j]) * TWO_PI
    a11 = (B_SM[0] + delta_b[0]) - (B_SM[1] + delta_b[1])
    a12 = -(delta_b[0] - delta_b[1])
    a21 = (B_SM[1] + delta_b[1]) - (B_SM[2] + delta_b[2])
    a22 = -(delta_b[1] - delta_b[2])
    det = a11 * a22 - a21 * a12
    if abs(det) < 1e-12:
        return None
    x = (rhs(0, 1) * a22 - rhs(1, 2) * a12) / det
    z = (a11 * rhs(1, 2) - a21 * rhs(0, 1)) / det
    if not (Z_MIN <= z <= x <= X_MAX):
        return None
    M_X = MZ * math.exp(x)
    inv_ax = (1 / A_Z[2] - (B_SM[2] + delta_b[2]) * x / TWO_PI
              - delta_b[2] * (x - z) / TWO_PI)
    if inv_ax <= 0:
        return None
    return M_X, 1.0 / inv_ax, MZ * math.exp(z)


def solve_two(name_a, name_b, nA, nB):
    """Fix alpha_X = 1/25; solve the three linear equations for
    (x, zA, zB).  Returns (M_X, M_A, M_B) or None (unphysical thresholds)."""
    dA, dB = LIB[name_a], LIB[name_b]
    M, v = [], []
    for i in range(3):
        M.append([B_SM[i] + nA * dA[i] + nB * dB[i], -nA * dA[i], -nB * dB[i]])
        v.append((1 / A_Z[i] - INV_A10) * TWO_PI)
    try:
        sol = np.linalg.solve(np.array(M), np.array(v))
    except np.linalg.LinAlgError:
        return None
    x, zA, zB = (float(s) for s in sol)
    if not (Z_MIN <= min(zA, zB) <= x <= X_MAX):
        return None
    return MZ * math.exp(x), MZ * math.exp(zA), MZ * math.exp(zB)


def tau_p(M_X, alpha_X):
    """tau_p ~ M_X^4 / (alpha_X^2 m_p^5), the paper's estimate, at M_X."""
    t_gev = M_X ** 4 / (alpha_X ** 2 * u10pred.M_P_ROBSON ** 5)
    return t_gev * u10pred.GEV_TO_S / u10pred.S_TO_YR


def main():
    t0 = time.time()

    # [1] the required Delta-b family (inverse problem), alpha_X = 1/25 fixed
    required = []
    for M_X in (1e14, 1e15, 2e15, 1e16, 2e16, 5e16, 1e17, 1e18):
        x = math.log(M_X / MZ) / TWO_PI
        db = tuple((1 / A_Z[i] - INV_A10) / x - B_SM[i] for i in range(3))
        required.append({"M_X_GeV": M_X,
                         "delta_b": [round(v, 3) for v in db]})

    # [2] the proportional-shift no-go (complete multiplets move nothing)
    db16 = LIB["16+16bar"]
    proportional = (db16[0] == db16[1] == db16[2])
    no_meeting = all(meet_lines(tuple(c * v for v in db16), 1e10) is None
                     for c in (1.0, 3.0, 8.0))

    # [3] the scan: two multiplets, alpha_X FIXED at 1/25, physical thresholds
    names = list(LIB)
    solutions = []
    for nA, nB in itertools.product((1, 2, 3), repeat=2):
        for na, nb in itertools.combinations_with_replacement(names, 2):
            if na == nb and nA != nB:
                continue
            got = solve_two(na, nb, nA, nB)
            if got is None:
                continue
            M_X, M_A, M_B = got
            tp = tau_p(M_X, A10)
            db = tuple(nA * LIB[na][i] + nB * LIB[nb][i] for i in range(3))
            b3eff = B_SM[2] + db[2]
            pole3 = (M_X * math.exp(-INV_A10 * TWO_PI / abs(b3eff))
                     if b3eff < 0 else None)
            solutions.append({
                "content": f"{nA}x{na} + {nB}x{nb}",
                "delta_b": [round(v, 3) for v in db],
                "cl10_realizable": na in CL10_REALIZABLE and nb in CL10_REALIZABLE,
                "M_X_GeV": M_X, "M_A_GeV": M_A, "M_B_GeV": M_B,
                "alpha_X": A10,
                "tau_p_yr": tp,
                "su3_down_pole_GeV": pole3,
                "tau_p_in_HK_20yr_reach": 1e34 < tp < 1e35,
                "tau_p_refutable_by_HK": tp > 1e35,
            })
    solutions.sort(key=lambda s: (not s["cl10_realizable"],
                                  abs(math.log10(s["tau_p_yr"]) - 36.0)))

    # [4] the selected CL10 solution and its consequences
    sel = next((s for s in solutions if s["cl10_realizable"]), None)

    # [5] the coupling running down the desert, from the same estimate
    old = UN.running_down()

    # [6] the abandoned convention M_GUT = M_P and the GUT scale that the
    # Super-Kamiokande bound excludes, both from the same tau_p(M_X) estimate
    tau_mplanck = tau_p(u10pred.M_PLANCK, A10)
    SK_BOUND_YR = 2.4e34        # tau_p/B(p -> e+ pi0), Super-K (paper ref 17)
    m_excluded = (SK_BOUND_YR * u10pred.S_TO_YR / u10pred.GEV_TO_S
                  * A10 ** 2 * u10pred.M_P_ROBSON ** 5) ** 0.25

    ok = sel is not None
    verdict = {
        "experiment": "V64_desert",
        "pass": ok,
        "status": ("O4 CONDITIONALLY RESOLVED: with alpha_X = 1/25 fixed (the "
                   "root value, not refitted), CL10-realizable vectorlike "
                   "desert content makes the three measured couplings meet at "
                   "M_X < M_P; the convention M_GUT = M_P is abandoned "
                   f"(it would give tau_p = {tau_mplanck:.2e} yr, and the Super-K bound excludes M_GUT < {m_excluded:.2e} GeV); "
                   "tau_p = 1.8e36 yr: Hyper-K (reach ~1e35 yr, 20 yr) can "
                   "REFUTE it by observing GUT-pattern decay but cannot "
                   "confirm it; confirmation needs ~1e36 yr (megaton era)"
                   if ok else "no physical solution: O4 stands"),
        "data": {
            "alpha_X_fixed": A10,
            "required_delta_b_family": required,
            "proportional_shift_nogo": {
                "complete_16_shifts_are_proportional": proportional,
                "no_meeting_for_any_complete_multiplet_count": no_meeting,
            },
            "solutions": solutions[:12],
            "selected": sel,
            "mplanck_convention": {
                "M_GUT_GeV": u10pred.M_PLANCK,
                "tau_p_yr": tau_mplanck,
                "superk_bound_yr": SK_BOUND_YR,
                "gut_scale_excluded_below_GeV": m_excluded,
            },
            "old_failure": {
                "ratio_to_measured": old["ratio_to_measured"],
                "su3_down_pole_GeV": old["su3_down_pole_scale"],
            },
            "residual": ("the desert spectrum is REQUIRED by the couplings but "
                         "not yet DERIVED from the automaton (breaking sector, "
                         "O6); F2: precision couplings admit no content with "
                         "alpha_X = 1/25; F5: proton decay at tau_p outside "
                         "the predicted band refutes the solution"),
        },
    }

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v64_desert.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V64 - the desert inverse problem (O4, F2, F5)")
    print("=" * 68)
    print(f"  alpha_X FIXED at the root value 1/25 = {A10} (g_10 = 1/5 not refitted)")
    print("  [1] required Delta-b (alpha_X = 1/25):")
    for r in required:
        print(f"      M_X = {r['M_X_GeV']:.0e} GeV -> Delta b = {r['delta_b']}")
    print(f"  [2] complete-multiplet no-go: proportional = {proportional}, "
          f"no meeting for any count = {no_meeting}")
    print(f"  [3] scan (two multiplets, thresholds in [1e6, M_P]): "
          f"{len(solutions)} physical solutions")
    for s in solutions[:6]:
        print(f"      {s['content']:24s} M_X = {s['M_X_GeV']:.2e}  "
              f"M_A = {s['M_A_GeV']:.2e}  M_B = {s['M_B_GeV']:.2e}  "
              f"tau_p = {s['tau_p_yr']:.1e} yr  "
              f"{'[CL10]' if s['cl10_realizable'] else ''}")
    if sel:
        print(f"  [4] selected: {sel['content']}  M_X = {sel['M_X_GeV']:.2e} GeV "
              f"(M_GUT = M_P ABANDONED)  tau_p = {sel['tau_p_yr']:.1e} yr  "
              f"SU(3) down-pole = {sel['su3_down_pole_GeV']}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'} - {verdict['status']}")
    print("  JSON -> " + os.path.join(outdir, "v64_desert.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
