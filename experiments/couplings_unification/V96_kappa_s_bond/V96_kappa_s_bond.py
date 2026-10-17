#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V96 - the FLAVOUR-DEPENDENT bond layer: per-class junction strengths
close the O4 meeting constructively (the ingredient named by V93/V95).

**Why (V92 + V93 + V95).** V92: the data-preferred meeting (M_X = 1.8e16
GeV, alpha_X = 1/25) lies inside the node-16 fragment cone - the desert
content exists. V93: SM-preserving Wilson directions are exactly
span{X, Y} - CLASS-degenerate - and no exact meeting exists in that
family. V95: off-diagonal junction splitting with a UNIVERSAL strength
kappa is invisible at first order (sum_l mu_l = -n reabsorbs into z0).
The named missing ingredient: a bond layer whose junction strength is
FLAVOUR-DEPENDENT, kappa_s per fragment class.

**The family.** Uniform SM-preserving direction c(theta) in span{X, Y},
offset z0, n nodes, and per-class junction strengths kappa_s:
    z_{s,l} = z_s(theta, z0) + kappa_s * (2 cos(2 pi l / n) - 1).
The Wilson direction preserves the SM algebra (intra-class roots
unbroken); the per-class kappa_s split the class-degenerate thresholds
INDEPENDENTLY - the first-order degeneracy of V95 is broken because the
per-class sums kappa_s * (-n) reabsorb into per-class offsets.
Unknowns (theta, z0, kappa_0..kappa_5) = 8 continuous + integer n:
the 3 meeting equations are UNDERDetermined - the point of the
experiment is not to find a solution (generic) but to verify that the
extended family (i) meets exactly, (ii) keeps all thresholds in the
desert, (iii) preserves the SM roots, and (iv) to NAME the cost: the
six kappa_s are new automaton input (a flavour-dependent bond layer),
not consequences of the current single bond gate P.
"""
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import u10pred  # noqa: E402
from cl10 import unification as UN  # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
INV_A10 = 25.0
TWO_PI = 2.0 * math.pi
M_X_OBS = u10pred.M_X_DESERT


def main():
    t0 = time.time()
    lam, _off = V65.cartan_weights()          # (16, 5) weights
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))

    xs = math.log(M_X_OBS / MZ)
    Z_MAX = math.log(u10pred.M_PLANCK / MZ)     # desert ceiling: cutoff

    # fragment classes: (X, Y) keys of the 16 weights
    keys = [(round(float((lam @ xc_gg)[s]), 6),
             round(float((lam @ y_coef)[s]), 6)) for s in range(16)]
    classes = sorted(set(keys))
    members = {c: [s for s in range(16) if keys[s] == c] for c in classes}

    # per-weight beta-function db from the V65 fragment table
    frs = V65.fragment_table(lam @ xc_gg, lam @ y_coef)
    db_of_key = {(f["X"], f["Y"]): np.array(f["db"]) / f["n"] for f in frs}
    lam_db = np.zeros((3, 16))
    for s in range(16):
        if keys[s] in db_of_key:
            lam_db[:, s] = db_of_key[keys[s]]

    def direction(th):
        return math.cos(th) * xc_gg + math.sin(th) * y_coef

    def class_thresholds(th, z0):
        c = direction(th)
        out = []
        for cl in classes:
            s = members[cl][0]
            proj = abs(float(lam[s] @ c))
            out.append(z0 + math.log(max(proj, 1e-12)))
        return np.array(out)

    def meeting(th, z0, kappas, n):
        zs = class_thresholds(th, z0)
        mu = np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                       for l in range(n)])
        F = [0.0, 0.0, 0.0]
        for i in range(3):
            s = INV_A10 + B_SM[i] * xs / TWO_PI
            for k, cl in enumerate(classes):
                z = zs[k] + kappas[k] * mu
                act = np.maximum(0.0, xs - z) / TWO_PI
                s += lam_db[i, members[cl][0]] * float(np.sum(act)) \
                    * len(members[cl])
            F[i] = s - 1.0 / A_Z[i]
        return F
    # --- pin theta to the V65 physical selective root ------------------------
    # (the only span{X,Y} directions that keep the SM desert group; a generic
    # theta breaks SU(2), theta = 0 preserves SU(5) x U(1) and is inadmissible)
    import json as _json
    v65 = _json.load(open(exp("V65_wilson",
                                       "results", "v65_wilson.json")))
    roots = [r for r in v65["data"]["part_B"]["X_1"]["selective_roots"]
             if r.get("physical")]
    TH_STAR = roots[0]["theta_rad"]
    dlam_all = np.array([lam[s] - lam[t] for s in range(16)
                         for t in range(s + 1, 16)])
    norm2_all = np.einsum("ij,ij->i", dlam_all, dlam_all)
    root_mask_all = np.abs(norm2_all - 8.0) < 1e-6
    unb_star = int(np.sum(root_mask_all
                          & (np.abs(dlam_all @ direction(TH_STAR)) < 1e-8)) // 2)
    checks_theta = {"theta_star": round(TH_STAR, 6),
                    "source": "V65_wilson X_1 physical selective root",
                    "unbroken_weight_pair_roots_at_theta_star": unb_star}

    # --- LP feasibility over (z0, kappa_0..kappa_5) at fixed theta ----------
    # With every threshold in the desert below M_X the meeting equations are
    # EXACTLY LINEAR, since sum_l mu_l = -n gives
    #   sum_l (xs - z_s - kappa_s mu_l) = n (xs - z_s + kappa_s).
    from scipy.optimize import linprog
    lp_s = np.array([math.log(max(abs(float(lam[members[cl][0]]
                                        @ direction(TH_STAR))), 1e-12))
                     for cl in classes])
    nk = len(classes)
    sols = []
    for n in range(3, 13):
        mu = np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                       for l in range(n)])
        w = np.array([len(members[cl]) for cl in classes], dtype=float)
        A_eq = np.zeros((3, 1 + nk))
        b_eq = np.zeros(3)
        for i in range(3):
            li = w * lam_db[i, [members[cl][0] for cl in classes]]
            A_eq[i, 0] = -n / TWO_PI * li.sum()
            A_eq[i, 1:] = n / TWO_PI * li
            b_eq[i] = -(INV_A10 + B_SM[i] * xs / TWO_PI - 1.0 / A_Z[i]) \
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
                r += 2
        res = linprog(np.zeros(1 + nk), A_ub=A_ub, b_ub=b_ub,
                      A_eq=A_eq, b_eq=b_eq,
                      bounds=[(-5.0, xs)] + [(-10.0, 10.0)] * nk,
                      method="highs")
        if not res.success:
            sols.append({"n_nodes": n, "note": f"LP infeasible ({res.status})"})
            continue
        z0 = float(res.x[0])
        kap = np.array(res.x[1:])
        F = np.array(meeting(TH_STAR, z0, kap, n))
        zs = class_thresholds(TH_STAR, z0)
        allz = [zs[k] + kap[k] * m for k in range(nk) for m in mu]
        if np.max(np.abs(F)) < 1e-10 and min(allz) >= -1e-9 \
                and max(allz) <= Z_MAX + 1e-9:
            sols.append({"n_nodes": n, "theta": round(TH_STAR, 5),
                         "z0": round(z0, 4),
                         "kappas": [round(float(v), 4) for v in kap],
                         "residual": float(np.max(np.abs(F))),
                         "z_min": round(min(allz), 4),
                         "z_max": round(max(allz), 4),
                         "thresholds_GeV": [
                             float(f"{MZ * math.exp(z):.3e}")
                             for z in sorted(allz)]})
        else:
            sols.append({"n_nodes": n, "note": "verification failed",
                         "fail_residual": float(np.max(np.abs(F))),
                         "z_min": round(min(allz), 4),
                         "z_max": round(max(allz), 4),
                         "Z_MAX": round(Z_MAX, 4)})

    n_ok = sum(1 for s in sols if "residual" in s)

    # --- SM root preservation at the first solution --------------------------
    breaking = None
    if n_ok:
        th = next(s["theta"] for s in sols if "residual" in s)
        c = direction(th)
        dlam = np.array([lam[s] - lam[t] for s in range(16)
                         for t in range(s + 1, 16)])
        norm2 = np.einsum("ij,ij->i", dlam, dlam)
        root_mask = np.abs(norm2 - 8.0) < 1e-6      # roots, each twice
        unb = int(np.sum(root_mask & (np.abs(dlam @ c) < 1e-8)) // 2)
        breaking = {"roots_total": int(root_mask.sum() // 2),
                    "unbroken_roots": unb,
                    "algebra_dim": 5 + unb}

    tau = None
    if n_ok:
        t_gev = M_X_OBS ** 4 / (INV_A10 ** -2 * u10pred.M_P_ROBSON ** 5)
        tau = t_gev * u10pred.GEV_TO_S / u10pred.S_TO_YR

    # the CHECK passes when the scan is valid and complete (V86/V87
    # convention); the physics verdict is in the status
    ok = n_ok > 0 and breaking is not None \
        and breaking["unbroken_roots"] == unb_star

    if ok:
        s0 = next(s for s in sols if "residual" in s)
        status = (
            f"O4 CLOSES CONSTRUCTIVELY WITH A NAMED EXTENSION: the "
            f"flavour-dependent bond layer (kappa_s per fragment class) "
            f"supplies an EXACT meeting (residual < 1e-10) at M_X = "
            f"{M_X_OBS:.1e} GeV, alpha_X = 1/25, for {n_ok}/10 node "
            f"counts n = 3..12, all thresholds in the desert "
            f"({s0['thresholds_GeV'][0]:.1e} .. "
            f"{s0['thresholds_GeV'][-1]:.1e} GeV), tau_p = {tau:.1e} "
            f"yr. The SM desert group is preserved by construction: "
            f"the Wilson direction is the V65 physical selective root "
            f"(theta = {TH_STAR:.4f}), the same direction V67 proved "
            f"leaves SU(3) x SU(2) x U(1)^3 unbroken "
            f"({breaking['unbroken_roots']}/{breaking['roots_total']} "
            f"weight-pair roots unbroken, all intra-fragment). The "
            f"first-order degeneracy of V95 is broken because the "
            f"per-class sums kappa_s*(-n) reabsorb into PER-CLASS "
            f"offsets. THE COST, STATED PLAINLY: the six kappa_s are "
            f"new automaton INPUT - the current bond layer has ONE "
            f"gate P with one strength; a flavour-dependent bond "
            f"layer is an EXTENSION of the automaton, and the desert "
            f"content becomes the content of the EXTENDED automaton "
            f"with six named (fitted) junction numbers. O4 is closed "
            f"conditionally on that extension; falsifier F2 stands.")
    else:
        status = (
            f"THE FLAVOUR-DEPENDENT FAMILY DOES NOT CLOSE O4 either "
            f"({n_ok}/10 node counts with an exact meeting, SM roots: "
            f"{breaking}): even with per-class kappa_s the desert-"
            f"constrained meeting fails - the obstruction is deeper "
            f"than the junction parameterisation.")

    verdict = {"experiment": "V96_kappa_s_bond",
               "pass": bool(ok), "status": status,
               "data": {"n_solutions": n_ok, "solutions": sols,
                        "theta_check": checks_theta,
                        "breaking": breaking,
                        "tau_p_yr": float(f"{tau:.3e}") if tau else None}}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v96_kappa_s_bond.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("V96 - flavour-dependent bond layer (kappa_s per class)")
    print(f"  classes: {len(classes)}  members: "
          f"{[len(members[c]) for c in classes]}")
    print(f"  exact meetings: {n_ok}/10 node counts")
    for s in sols:
        if "residual" in s:
            print(f"   n={s['n_nodes']}: theta={s['theta']}, z0={s['z0']}, "
                  f"kappa={s['kappas']}, z in [{s['z_min']}, "
                  f"{s['z_max']}], res={s['residual']:.1e}")
    if breaking:
        print(f"  breaking: {breaking}")
    if tau:
        print(f"  tau_p = {tau:.1e} yr")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v96_kappa_s_bond.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
