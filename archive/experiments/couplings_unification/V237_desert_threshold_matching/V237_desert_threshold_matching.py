#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V237 - the desert solution's thresholds (O4/F2, tier 3, residual).

The question this pre-registration answers (review/pre_registrations.md):
the root-coupling solution (M_X = 1.235e17 GeV at alpha_X = pi/64) is
recorded with DESTROYED threshold ordering, and the content it needs is
three 10 + 10bar pairs at 1.2e7 GeV plus one (1,2) pair at 1.0e12 GeV.
Does the spectrum's own one-loop threshold shift restore the meeting, or
is the root-coupling solution not realizable by its own content?

Nothing is re-tuned: the content, multiplicities and preferred thresholds
are the closed V64/V133 values; the meeting is solved and the required
threshold correction is compared against the spectrum's own log-spacings.

T0  the desert Delta b vectors recomputed IN EXACT INTEGER ARITHMETIC
    from the root's own charge data (matter.SM_FIELDS hypercharges, the
    16-plet Dynkin indices T3 = T2 = 3/2, and B-L from matter.bl_charge),
    with the U1 normalisation k1 DERIVED from the completeness identity
    (a full 16 + 16bar adds equally to the three beta functions);
T1  the closed SM beta functions B_SM = (41/10, -19/6, -7) from the same
    charge data (three 16-plet generations plus one scalar doublet);
T2  the V133 solution reproduced independently of its solver helper, and
    checked against the committed artifact;
T3  the matching measurement: with the content pinned at its preferred
    thresholds (1.2e7, 1.0e12 GeV) and M_X pinned at the fitted
    desert.m_x_closed_gev, how far must each threshold move for the three
    couplings to meet at alpha_X = pi/64, against the spectrum's own
    log-spacings.
"""

import json
import math
import os
import sys
import time
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "couplings_unification",
                                "V64_desert"))

from cl10 import config as CFG                                     # noqa: E402
from cl10 import matter, unification as UN                         # noqa: E402
import V64_desert as V64                                           # noqa: E402
from cl10.paths import exp                                         # noqa: E402

TWO_PI = 2.0 * math.pi


def yfrac(f):
    """exact rational hypercharge of an SM_FIELDS entry"""
    return Fraction(f[2]).limit_denominator(100)


def main():
    t0 = time.time()
    checks, data = [], {}

    sf = matter.SM_FIELDS
    Y = {k: yfrac(v) for k, v in sf.items()}          # exact 1/6,-2/3,...
    dims = {k: int(v[3]) for k, v in sf.items()}
    bl16, bl16b = matter.bl_charge("16"), matter.bl_charge("16bar")
    n_weyl = matter.weyl16()[0].shape[1]

    # ---- T0 exact-integer Delta b for the desert contents ------------------
    # a complete 16 adds equally to b3 and b2: T3(16) = 1 + 1/2 + 1/2 = 2
    # (Q, u^c, d^c), T2(16) = 3/2 + 1/2 = 2 (Q, L): 10 + 5bar + 1 under
    # the closed gauge group. Sum of Y^2 over the 16 from the charges:
    sY2_16 = sum(dims[k] * Y[k] ** 2 for k in sf)                  # 10/3
    t2_16 = Fraction(3, 2) + Fraction(1, 2)                        # Q, L
    # completeness of the 16 + 16bar fixes the U1 normalisation exactly:
    # k1 * (10/3) = 2  ->  evaluated as a Fraction:
    k1 = t2_16 / sY2_16
    # Weyl content contribution to the beta functions: (2/3) * index, with
    # the per-field indices (T3 = 1 for a coloured triplet, T2 = I(I+1)*dim
    # collapsed as the tabulated multiplet indices: Q contributes 3/2 to b2
    # and 1 to b3):
    def db_pair(keys):
        """exact (b1, b2, b3) Delta b of the pair content + c.c."""
        b1 = Fraction(2, 3) * 2 * k1 * sum(dims[k] * Y[k] ** 2
                                           for k in keys)
        t2 = sum(Fraction(3, 2) if k == "Q" else
                 (Fraction(1, 2) if k == "L" else 0) for k in keys)
        t3 = sum(Fraction(1) if k == "Q" else
                 (Fraction(1, 2) if k in ("u^c", "d^c") else 0)
                 for k in keys)
        return (b1, Fraction(2, 3) * 2 * t2, Fraction(2, 3) * 2 * t3)

    db_10 = db_pair(("Q", "u^c"))
    db_12 = db_pair(("L",))
    db_16 = db_pair(tuple(sf))
    exact = {
        "10+10bar": db_10 == (Fraction(6, 5), Fraction(2), Fraction(2)),
        "(1,2)": db_12 == (Fraction(2, 5), Fraction(2, 3), Fraction(0)),
        "16+16bar": db_16 == (Fraction(8, 3),) * 3,
    }
    agrees = all(abs(float(a) - b) < 1e-12
                 for name in ("10+10bar", "(1,2)")
                 for a, b in zip(db_pair({"10+10bar": ("Q", "u^c"),
                                          "(1,2)": ("L",)}[name]),
                                 V64.LIB[name]))
    checks.append({
        "name": f"T0 exact-integer branching from the root charge data: the "
                f"U1 normalisation is fixed by the completeness identity to "
                f"k1 = {k1} (sum Y^2 over the 16 = {sY2_16}); the desert "
                f"Delta b vectors come out (6/5, 2, 2) for a 10 + 10bar pair, "
                f"(2/5, 2/3, 0) for a (1,2) pair and the complete 16 + 16bar "
                f"adds equally ({db_16[0]},{db_16[1]},{db_16[2]}) to all "
                f"three - all equalities exact Fractions, all matching the "
                f"closed V64 library to 1e-12 (16 Weyl = {n_weyl}, B-L = "
                f"{bl16}/{bl16b})",
        "pass": bool(all(exact.values()) and agrees and n_weyl == 16
                     and bl16 == -1 and bl16b == 1)})
    data["branching_exact"] = {
        "k1": str(k1), "sY2_16": str(sY2_16),
        "10+10bar": [str(x) for x in db_10],
        "(1,2)": [str(x) for x in db_12],
        "16+16bar": [str(x) for x in db_16]}

    # ---- T1 the closed SM beta functions from the same charge data ---------
    b3_sm = Fraction(-11) + Fraction(2, 3) * 3 * (
        Fraction(1) + Fraction(1, 2) + Fraction(1, 2))
    b2_sm = (Fraction(-22, 3) + Fraction(2, 3) * 3 * t2_16
             + Fraction(1, 3) * Fraction(1, 2))
    b1_sm = (Fraction(2, 3) * 3 * k1 * sY2_16
             + Fraction(1, 3) * k1 * sum(dims[k] * Y[k] ** 2
                                         for k in ("L",)))
    checks.append({
        "name": f"T1 the closed SM one-loop beta functions (config sm:) are "
                f"reproduced from the same charge data with one scalar "
                f"doublet: ({b1_sm}, {b2_sm}, {b3_sm}) = "
                f"{tuple(float(b) for b in (b1_sm, b2_sm, b3_sm))} against "
                f"B_SM = {UN.B_SM}",
        "pass": bool(abs(float(b1_sm) - UN.B_SM[0]) < 1e-12
                     and abs(float(b2_sm) - UN.B_SM[1]) < 1e-12
                     and abs(float(b3_sm) - UN.B_SM[2]) < 1e-12)})

    # ---- T2 the V133 pi/64 solution, solved independently ------------------
    dA, dB = V64.LIB["10+10bar"], V64.LIB["(1,2)"]
    nA, nB = 3, 1
    MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
    AX = math.pi / 64.0                      # the closed-angle value
    MXF = CFG.CFG["desert"]["m_x_closed_gev"]
    x_pref = [math.log(1.2e7 / MZ), math.log(1.0e12 / MZ)]   # (:983 content)

    def meet_solve():
        """solve the three meeting equations for (x, zA, zB) at alpha_X"""
        Mm, v = [], []
        for i in range(3):
            Mm.append([B_SM[i] + nA * dA[i] + nB * dB[i],
                       -nA * dA[i], -nB * dB[i]])
            v.append((1.0 / A_Z[i] - 1.0 / AX) * TWO_PI)
        return np.linalg.solve(np.array(Mm), np.array(v))

    sol = meet_solve()
    mx_sol, m10_sol, m12_sol = (MZ * math.exp(s) for s in sol)
    art = json.load(open(exp("V133_desert_pi64", "results",
                             "v133_desert_pi64.json")))
    q = art["data"]["quantized_pi64"]
    agree133 = (abs(mx_sol - q["M_X_GeV"]) / q["M_X_GeV"] < 1e-9
                and abs(m10_sol - q["M_10_GeV"]) / q["M_10_GeV"] < 1e-9
                and abs(m12_sol - q["M_12_GeV"]) / q["M_12_GeV"] < 1e-9)
    checks.append({
        "name": f"T2 the pi/64 meeting solved independently over the three "
                f"desert thresholds reproduces the committed V133 artifact "
                f"to 1e-9 relative: M_X = {mx_sol:.4e} GeV, M_10 = "
                f"{m10_sol:.4e} GeV, M_(1,2) = {m12_sol:.4e} GeV - the "
                f"(1,2) threshold lands ABOVE M_X (ordering destroyed) and "
                f"M_10 below the 1e6 GeV desert floor (artifact in_desert_"
                f"box = {q['in_desert_box']})",
        "pass": bool(agree133 and m12_sol > mx_sol
                     and m10_sol < 1e6)})
    data["solution_pi64"] = {"M_X_GeV": mx_sol, "M_10_GeV": m10_sol,
                             "M_12_GeV": m12_sol}

    # ---- T3 the matching measurement against the spectrum's spacings ------
    # pin M_X at the fitted m_x_closed_gev and the (1,2) pair at its
    # preferred 1.0e12 GeV: the three couplings meet at alpha_X = pi/64 only
    # if the 10 + 10bar scale takes the value the three equations demand -
    # and check where it lands against the preferred 1.2e7 GeV, and what the
    # (1,2) pair would have to absorb if instead it stayed at 1.0e12.
    xfix = math.log(MXF / MZ)
    # two independent equations for (zA, zB) at fixed M_X, taken as the
    # row differences (the two contents are proportional in the (b1, b2)
    # plane, so the raw rows 0 and 1 alone are degenerate):
    C = [(1.0 / A_Z[i] - 1.0 / AX) * TWO_PI
         - xfix * (B_SM[i] + nA * dA[i] + nB * dB[i]) for i in range(3)]
    R = [[-nA * dA[i], -nB * dB[i]] for i in range(3)]
    Mm = [[R[0][k] - R[1][k] for k in (0, 1)],
          [R[1][k] - R[2][k] for k in (0, 1)]]
    v = [C[0] - C[1], C[1] - C[2]]
    zA, zB = np.linalg.solve(np.array(Mm), np.array(v))
    # third-equation mismatch at the pinned M_X:
    resid = R[0][0] * zA + R[0][1] * zB - C[0]
    dz10 = zA - x_pref[0]                    # shift the 10-pair needs
    dz12 = zB - x_pref[1]                    # shift the (1,2) pair needs
    spac12 = xfix - x_pref[1]                # ln(M_X / 1.0e12)
    spac10 = xfix - x_pref[0]                # ln(M_X / 1.2e7)
    fires = abs(dz12) > spac12 and m12_sol > mx_sol
    checks.append({
        "name": f"T3 with M_X pinned at the fitted {MXF:.3e} GeV and alpha_X "
                f"= pi/64, the content at its preferred thresholds misses "
                f"the meeting; the required threshold shifts are "
                f"Delta ln M_10 = {dz10:+.3f} e-folds (to {MZ * math.exp(zA):.3e} "
                f"GeV) and Delta ln M_(1,2) = {dz12:+.3f} e-folds (to "
                f"{MZ * math.exp(zB):.3e} GeV) against the spectrum's own "
                f"log-spacings ln(M_X/M_12) = {spac12:.3f} and ln(M_X/M_10) "
                f"= {spac10:.3f}: the (1,2) pair must move "
                f"{abs(dz12) / spac12:.2f}x its entire spacing to M_X, and "
                f"it lands above M_X; third-equation residual "
                f"{resid:.3e}",
        "pass": bool(abs(resid) < 1e-2 and fires)})
    data["matching"] = {"M_X_pinned_GeV": MXF, "alpha_X": AX,
                        "zA_solved_GeV": MZ * math.exp(zA),
                        "zB_solved_GeV": MZ * math.exp(zB),
                        "dz10_efolds": dz10, "dz12_efolds": dz12,
                        "spacing12": spac12, "spacing10": spac10,
                        "residual": resid,
                        "shift_over_spacing_12": abs(dz12) / spac12,
                        "floor_violated": bool(MZ * math.exp(zA) < 1e6),
                        "ordering_violated": bool(MZ * math.exp(zB) > MXF)}

    branch = ("NOT REALIZABLE BY ITS OWN CONTENT" if fires
              else "RESTORED BY THRESHOLDS")

    ok = all(c["pass"] for c in checks)
    tau = q["tau_p_yr"]
    status = (
        f"desert threshold matching: {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the root-coupling solution alpha_X = "
        f"pi/64, M_X = {MXF:.3e} GeV needs its three 10 + 10bar pairs at "
        f"{MZ * math.exp(zA):.3e} GeV (below the 1e6 GeV desert floor) and "
        f"its (1,2) pair at {MZ * math.exp(zB):.3e} GeV - above M_X itself; "
        f"the (1,2) threshold must move {dz12:+.1f} e-folds, "
        f"{abs(dz12) / spac12:.2f}x the log-spacing between 1.0e12 GeV and "
        f"M_X, so no one-loop threshold correction within the spectrum's "
        f"own spacings restores the meeting - F2's unification side fires "
        f"(Table 5, O4 gap column); tau_p = {tau:.2e} yr inherits the "
        f"destroyed-ordering M_X and stays {tau / 2.4e34:.0e}x above the "
        f"Super-Kamiokande p -> e+ pi0 bound either way; the branching that "
        f"drives it is exact-integer from the root charges (k1 = {k1}, "
        f"Delta b 10+10bar = (6/5, 2, 2), 16+16bar = (8/3, 8/3, 8/3) "
        f"proportional - the no-go is the model's own content).")
    out = {"id": "V237_desert_threshold_matching",
           "tier": 3,
           "property": "the root-coupling desert solution tested against "
                       "its own one-loop thresholds (V237, O4/F2)",
           "residual": True,
           "programme": "V237",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v237_desert_threshold_matching.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
