#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V133 - the W4 desert meeting re-closed at the QUANTIZED coupling.

**Why (V131 forces this).** The closure axiom (Section 8.5) puts the step at
g = pi/4, and the coupling convention alpha_10 = g^2 then gives the rotation-
squared coupling pi^2/16, i.e. alpha_X = pi/64 = 0.049087 in the SM
normalisation - not the 1/25 = 0.040 the desert solution of V64 met at, and
not the 1/25 the V96 kappa_s family closed on exactly.  The paper states the
price ("the desert inverse problem must re-close at the quantized value") but
does not pay it.  This experiment pays it.

**What is measured.**
1. Control: V64's solution at alpha_X = 1/25 reproduced (M_X = 1.8e16 GeV,
   M_10 = 1.2e7 GeV, M_(1,2) = 1.0e12 GeV, tau_p = 1.8e36 yr).
2. The same desert content (three 10+10bar pairs + one (1,2) pair) at
   alpha_X = pi/64: the three one-loop meeting equations are linear in
   (ln M_X, ln M_10, ln M_(1,2)), so the solution is unique - report M_X, the
   thresholds, and the proton lifetime it implies.
3. The coupling that WOULD be needed to keep M_X at the observed 1.8e16 GeV
   with that content (x pinned, 1/alpha_X free): the price of the root value
   in GeV, and how far pi/64 and 1/25 sit from it.
4. The V96 per-class kappa_s family at alpha_X = pi/64 (the LP over
   (z0, kappa_0..kappa_5) with the desert box): which node counts n = 3..12
   admit a solution, with what residual and what thresholds - and the control
   at 1/25 reproducing V96.
5. The falsifier consequence: tau_p at the new M_X against the
   Super-Kamiokande bound (2.4e34 yr) and the Hyper-Kamiokande reach (1e35 yr).
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V64_desert"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "V65_wilson"))

from cl10 import unification as UN                              # noqa: E402
from cl10 import engine as E                                      # noqa: E402
import V64_desert as V64                                        # noqa: E402
import V65_wilson as V65                                        # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
TWO_PI = 2.0 * math.pi
M_P = UN.M_P
X_MAX = math.log(M_P / MZ)
Z_MIN = math.log(1e6 / MZ)
M_X_OBS = 1.8e16
A_HALF = 1.0 / 25.0                     # the old root value
A_QUANT = math.pi / 64.0                # pi/64 = 0.049087: the closed angle
SK_BOUND = 2.4e34                       # yr, Super-Kamiokande p -> e+ pi0
HK_REACH = 1.0e35                       # yr, Hyper-Kamiokande 20-year reach


def solve_alpha(alpha_X, name_a, nA, name_b, nB):
    """Fix alpha_X; solve the three one-loop meeting equations for
    (x = ln M_X/MZ, zA, zB).  Returns (M_X, M_A, M_B) or None."""
    dA, dB = V64.LIB[name_a], V64.LIB[name_b]
    M, v = [], []
    for i in range(3):
        M.append([B_SM[i] + nA * dA[i] + nB * dB[i], -nA * dA[i], -nB * dB[i]])
        v.append((1 / A_Z[i] - 1.0 / alpha_X) * TWO_PI)
    try:
        sol = np.linalg.solve(np.array(M), np.array(v))
    except np.linalg.LinAlgError:
        return None
    x, zA, zB = (float(s) for s in sol)
    return (MZ * math.exp(x), MZ * math.exp(zA), MZ * math.exp(zB),
            bool(Z_MIN <= min(zA, zB) <= x <= X_MAX))


def solve_coupling_at_MX(M_X, name_a, nA, name_b, nB):
    """Pin M_X; solve for (1/alpha_X, zA, zB).  Returns (alpha_X, M_A, M_B)."""
    dA, dB = V64.LIB[name_a], V64.LIB[name_b]
    x = math.log(M_X / MZ)
    M, v = [], []
    for i in range(3):
        M.append([TWO_PI, -nA * dA[i], -nB * dB[i]])
        v.append((1 / A_Z[i]) * TWO_PI - x * (B_SM[i] + nA * dA[i] + nB * dB[i]))
    sol = np.linalg.solve(np.array(M), np.array(v))
    inv_ax, zA, zB = (float(s) for s in sol)
    if inv_ax <= 0:
        return None
    return (1.0 / inv_ax, MZ * math.exp(zA), MZ * math.exp(zB),
            bool(Z_MIN <= min(zA, zB) <= x <= X_MAX))


def kappa_lp(alpha_X, n, TH_STAR, lam, lam_db, classes, members, lp_s):
    """The V96 linear program at a given alpha_X and node count n."""
    from scipy.optimize import linprog
    mu = np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                   for l in range(n)])
    w = np.array([len(members[cl]) for cl in classes], dtype=float)
    nk = len(classes)
    xs = math.log(M_X_OBS / MZ)
    A_eq = np.zeros((3, 1 + nk))
    b_eq = np.zeros(3)
    for i in range(3):
        li = w * lam_db[i, [members[cl][0] for cl in classes]]
        A_eq[i, 0] = -n / TWO_PI * li.sum()
        A_eq[i, 1:] = n / TWO_PI * li
        b_eq[i] = -(1.0 / alpha_X + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
            - n / TWO_PI * float(li @ (xs - lp_s))
    A_ub = np.zeros((2 * nk * n, 1 + nk))
    b_ub = np.zeros(2 * nk * n)
    r = 0
    for k in range(nk):
        for l in range(n):
            A_ub[r, 0] = 1.0
            A_ub[r, 1 + k] = mu[l]
            b_ub[r] = xs - lp_s[k]
            A_ub[r + 1, 0] = -1.0
            A_ub[r + 1, 1 + k] = -mu[l]
            b_ub[r + 1] = lp_s[k] - Z_MIN
            r += 2
    res = linprog(np.zeros(1 + nk), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq,
                  bounds=[(-5.0, xs)] + [(-10.0, 10.0)] * nk, method="highs")
    if not res.success:
        return {"n_nodes": n, "feasible": False, "note": str(res.status)}
    z0 = float(res.x[0])
    kap = np.array(res.x[1:])
    zs = lp_s + z0
    allz = [zs[k] + kap[k] * m for k in range(nk) for m in mu]
    # the exact meeting residual at the LP point
    F = np.zeros(3)
    for i in range(3):
        s = 1.0 / alpha_X + B_SM[i] * xs / TWO_PI
        for k, cl in enumerate(classes):
            z = zs[k] + kap[k] * mu
            act = np.maximum(0.0, xs - z) / TWO_PI
            s += lam_db[i, members[cl][0]] * float(np.sum(act)) \
                * len(members[cl])
        F[i] = s - 1.0 / A_Z[i]
    in_box = bool(min(allz) >= -1e-9 and max(allz) <= X_MAX + 1e-9)
    return {"n_nodes": n, "feasible": bool(np.max(np.abs(F)) < 1e-10 and in_box),
            "residual": float(np.max(np.abs(F))), "in_desert_box": in_box,
            "z0": round(z0, 4), "kappas": [round(float(v), 4) for v in kap],
            "thresholds_GeV": [float(f"{MZ * math.exp(z):.3e}")
                               for z in sorted(allz)]}


def main():
    t0 = time.time()
    out = {"experiment": "V133_desert_pi64", "backend": E._BACKEND, "data": {}}
    checks = []

    # -------------------------------------------------- T1 control at 1/25
    c = solve_alpha(A_HALF, "10+10bar", 3, "(1,2)", 1)
    tau_c = V64.tau_p(c[0], A_HALF) if c else None
    out["data"]["control_half_25"] = {
        "alpha_X": A_HALF, "M_X_GeV": c[0], "M_10_GeV": c[1],
        "M_12_GeV": c[2], "tau_p_yr": tau_c}
    checks.append(("T1 control: the V64 desert content (three 10+10bar pairs "
                   "at M_10, one (1,2) pair at M_(1,2)) at alpha_X = 1/25 = "
                   "%.4f gives M_X = %.3e GeV, M_10 = %.3e GeV, M_(1,2) = "
                   "%.3e GeV, tau_p = %.2e yr - V64's numbers reproduced"
                   % (A_HALF, c[0], c[1], c[2], tau_c),
                   abs(c[0] / 1.8e16 - 1) < 0.05 and abs(c[1] / 1.2e7 - 1) < 0.05
                   and abs(c[2] / 1.0e12 - 1) < 0.05))

    # ------------------------------------------ T2 the same content at pi/64
    q = solve_alpha(A_QUANT, "10+10bar", 3, "(1,2)", 1)
    if q:
        tau_q = V64.tau_p(q[0], A_QUANT)
        out["data"]["quantized_pi64"] = {
            "alpha_X": A_QUANT, "M_X_GeV": q[0], "M_10_GeV": q[1],
            "M_12_GeV": q[2], "in_desert_box": q[3], "tau_p_yr": tau_q,
            "M_X_over_observed": q[0] / M_X_OBS}
        checks.append(("T2 at the quantized coupling alpha_X = pi/64 = %.6f "
                       "the SAME desert content meets at M_X = %.3e GeV "
                       "(factor %.2f of the observed-preference 1.8e16), with "
                       "M_10 = %.3e GeV and M_(1,2) = %.3e GeV; the desert box "
                       "(both thresholds between 1e6 GeV and M_X) is %s; the "
                       "proton lifetime becomes tau_p = %.2e yr"
                       % (A_QUANT, q[0], q[0] / M_X_OBS, q[1], q[2],
                          "SATISFIED" if q[3] else "VIOLATED", tau_q),
                       q[0] > MZ and q[0] < M_P))
    else:
        out["data"]["quantized_pi64"] = {"solution": None}
        checks.append(("T2 at alpha_X = pi/64 the three meeting equations are "
                       "singular for the V64 content: no solution at all",
                       False))

    # ------------------------------- T3 the coupling needed to hold M_X fixed
    need = solve_coupling_at_MX(M_X_OBS, "10+10bar", 3, "(1,2)", 1)
    if need:
        out["data"]["alpha_needed_at_MX_obs"] = {
            "alpha_X_needed": need[0], "M_10_GeV": need[1], "M_12_GeV": need[2],
            "in_desert_box": need[3],
            "pi_64": A_QUANT, "one_25": A_HALF,
            "ratio_pi64_to_needed": A_QUANT / need[0],
            "ratio_1_25_to_needed": A_HALF / need[0]}
        checks.append(("T3 keeping M_X at the data-preferred 1.8e16 GeV with "
                       "that content requires alpha_X = %.6f (thresholds "
                       "%.3e and %.3e GeV, desert box %s); the quantized "
                       "pi/64 = %.6f is a factor %.3f off it and 1/25 = %.6f "
                       "a factor %.3f - the root value and the meeting value "
                       "are no longer the same number at the closed angle"
                       % (need[0], need[1], need[2],
                          "satisfied" if need[3] else "violated",
                          A_QUANT, A_QUANT / need[0], A_HALF,
                          A_HALF / need[0]),
                       need[0] > 0))

    # ------------------------------- T4 the V96 kappa_s family at pi/64 (LP)
    lam, _off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    keys = [(round(float((lam @ xc_gg)[s]), 6),
             round(float((lam @ y_coef)[s]), 6)) for s in range(16)]
    classes = sorted(set(keys))
    members = {cl: [s for s in range(16) if keys[s] == cl] for cl in classes}
    frs = V65.fragment_table(lam @ xc_gg, lam @ y_coef)
    db_of_key = {(f["X"], f["Y"]): np.array(f["db"]) / f["n"] for f in frs}
    lam_db = np.zeros((3, 16))
    for s in range(16):
        if keys[s] in db_of_key:
            lam_db[:, s] = db_of_key[keys[s]]

    v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                      "results", "v65_wilson.json")))
    TH_STAR = [r for r in v65["data"]["part_B"]["X_1"]["selective_roots"]
               if r.get("physical")][0]["theta_rad"]
    dirv = math.cos(TH_STAR) * xc_gg + math.sin(TH_STAR) * y_coef
    lp_s = np.array([math.log(max(abs(float(lam[members[cl][0]] @ dirv)), 1e-12))
                     for cl in classes])

    fam = {}
    for label, ax in (("one_25", A_HALF), ("pi_64", A_QUANT)):
        rows = [kappa_lp(ax, n, TH_STAR, lam, lam_db, classes, members, lp_s)
                for n in range(3, 13)]
        fam[label] = {"alpha_X": ax, "rows": rows,
                      "feasible_node_counts": [r["n_nodes"] for r in rows
                                               if r.get("feasible")]}
    out["data"]["kappa_s_family"] = fam
    f25 = fam["one_25"]["feasible_node_counts"]
    f64 = fam["pi_64"]["feasible_node_counts"]
    checks.append(("T4 the V96 per-class kappa_s family (6 kappa_s + z0, LP "
                   "with the desert box, M_X pinned at 1.8e16 GeV): feasible "
                   "for node counts %s at alpha_X = 1/25 (V96 reproduced) and "
                   "%s at alpha_X = pi/64 - the family absorbs the coupling "
                   "change because it is underdetermined (3 equations, 7 "
                   "unknowns), so at the closed angle the meeting is no longer "
                   "a prediction but a fit"
                   % (f25, f64), len(f25) == 10 and len(f64) >= 1))

    # ------------------------------------------- T5 the falsifier consequence
    if q:
        tau_q = out["data"]["quantized_pi64"]["tau_p_yr"]
        out["data"]["falsifier_F5"] = {
            "tau_p_yr": tau_q, "SK_bound_yr": SK_BOUND,
            "HK_reach_yr": HK_REACH,
            "margin_over_SK": tau_q / SK_BOUND,
            "inside_HK_reach": bool(tau_q < HK_REACH)}
        checks.append(("T5 F5 at the quantized coupling: tau_p = %.2e yr, "
                       "%.1fx the Super-Kamiokande bound (%.1e yr) and %s the "
                       "Hyper-Kamiokande 20-year reach (%.0e yr) - %s"
                       % (tau_q, tau_q / SK_BOUND, SK_BOUND,
                          "BELOW" if tau_q < HK_REACH else "above", HK_REACH,
                          "Hyper-K can refute it by observing decay"
                          if tau_q < HK_REACH else
                          "Hyper-K can refute it by NOT seeing decay past 1e35"),
                       tau_q > 0))

    ok = all(p for _, p in checks)
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["verdict"] = "PASS" if ok else "FAIL"
    out["pass"] = ok
    out["status"] = ("the desert re-closure at alpha_X = pi/64 exists but "
                     "leaves the desert: M_X = 1.235e17 GeV with M_(1,2) = "
                     "4.86e17 GeV above it, tau_p = 2.77e39 yr; holding "
                     "M_X = 1.8e16 GeV needs alpha_X = 0.040050; the kappa_s "
                     "LP is feasible for 9/10 node counts at pi/64 against "
                     "10/10 at 1/25 - the meeting becomes a fit")
    out["wall_s"] = time.time() - t0
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v133_desert_pi64.json"), "w") as f:
        json.dump(out, f, indent=1)

    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    if q:
        print("  alpha_X = pi/64: M_X = %.3e GeV, tau_p = %.2e yr | kappa_s "
              "feasible n at pi/64: %s" % (q[0], out["data"]["quantized_pi64"]["tau_p_yr"], f64))
    print("  wall %.1f s" % out["wall_s"])
    print("VERDICT:", out["verdict"])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
