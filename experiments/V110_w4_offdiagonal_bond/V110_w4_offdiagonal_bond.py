#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V110 - are V96's six kappa_s REALIZABLE by one off-diagonal bond operator?

**Why (V96 + V105).** V96 closes the W4 meeting exactly with per-class
junction strengths kappa_s (six named inputs). V105 proved the kappa_s are
NOT the image of any on-site Cartan element: SM preservation forces
span{X,Y} (class-degenerate) and every Cartan splitting the classes breaks
the desert roots. V105's verdict: the required flavour dependence is
OFF-DIAGONAL bond structure. This experiment constructs that structure and
tests it.

**The construction.** The bond dressing is a Hermitian off-diagonal operator
J = sum_{root pairs (s,t)} c_st (E_s,t + E_t,s) on the 16, built from the
root vectors of the weight graph (weight pairs with |lam_s - lam_t|^2 = 8).
Its first-order class effect is degenerate (V95: universal kappa reabsorbs
into z0). The FIRST class-dependent effect is second-order level repulsion:
the class channel k (energy z_k = z0 + lp_k, the V96 threshold) shifts by
    kappa_k = sum_{m != k} C_km c_km^2 / (z_k - z_m),
C_km = the number of root-vector channels between classes k and m (the
degeneracy factor of the bond graph). ONE bond gate, one set of couplings
c_km >= 0, one offset z0: the six kappa_s become functions of K + 1
automaton parameters, K = the number of class-pair channels.

**The test.** Minimise the V96 meeting residual over (z0, c) with the
desert inequalities, and separately test the SIGN STRUCTURE: level
repulsion forces sign(kappa_k) to follow the energy ordering (the lowest
class can only be pushed down). Report whether the meeting closes
(residual < 1e-10, desert ok) under the repulsion convention, and the
exact obstruction if not.

PASS = the computation is self-consistent (graph verified, meeting
residual verified at the reported point) and the verdict is reported with
exact numbers. A negative realizability result PASSES by reporting the
obstruction.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import u10pred  # noqa: E402
from cl10 import unification as UN  # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
INV_A10 = 25.0
TWO_PI = 2.0 * math.pi
M_X_OBS = 1.8e16


def main():
    lam, _off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))

    xs = math.log(M_X_OBS / MZ)
    Z_MAX = math.log(u10pred.M_PLANCK / MZ)

    keys = [(round(float((lam @ xc_gg)[s]), 6),
             round(float((lam @ y_coef)[s]), 6)) for s in range(16)]
    classes = sorted(set(keys))
    members = {c: [s for s in range(16) if keys[s] == c] for c in classes}
    nk = len(classes)
    w = np.array([len(members[cl]) for cl in classes], dtype=float)

    frs = V65.fragment_table(lam @ xc_gg, lam @ y_coef)
    db_of_key = {(f["X"], f["Y"]): np.array(f["db"]) / f["n"] for f in frs}
    lam_db = np.zeros((3, 16))
    for s in range(16):
        if keys[s] in db_of_key:
            lam_db[:, s] = db_of_key[keys[s]]
    reps = [members[cl][0] for cl in classes]
    li_all = np.array([w * lam_db[i, reps] for i in range(3)])  # (3,6)

    v65 = json.load(open(os.path.join(ROOT, "experiments", "V65_wilson",
                                      "results", "v65_wilson.json")))
    roots = [r for r in v65["data"]["part_B"]["X_1"]["selective_roots"]
             if r.get("physical")]
    TH_STAR = roots[0]["theta_rad"]

    def direction(th):
        return math.cos(th) * xc_gg + math.sin(th) * y_coef

    lp_s = np.array([math.log(max(abs(float(lam[s] @ direction(TH_STAR))),
                                  1e-12)) for s in reps])

    # ---- the root-vector bond graph between classes ------------------------
    # root pairs: weight pairs with |lam_s - lam_t|^2 = 8 (the 120 roots of
    # D5 restricted to the 16 weights). Channel (k,m): root pairs joining
    # class k to class m. C_km = channel degeneracy.
    cl_of = {s: next(k for k, cl in enumerate(classes)
                     if s in members[cl]) for s in range(16)}
    C = np.zeros((nk, nk))
    n_root_pairs = 0
    for s in range(16):
        for t in range(s + 1, 16):
            d = lam[s] - lam[t]
            if abs(float(d @ d) - 8.0) < 1e-6:
                n_root_pairs += 1
                k, m = cl_of[s], cl_of[t]
                if k != m:
                    C[k, m] += 1.0
                    C[m, k] += 1.0
    channels = [(k, m) for k in range(nk) for m in range(k + 1, nk)
                if C[k, m] > 0]
    K = len(channels)

    def kappas_from(c2, z0):
        """Second-order level repulsion of the class channels."""
        kap = np.zeros(nk)
        z = z0 + lp_s
        for (k, m), cc in zip(channels, c2):
            dz = z[k] - z[m]
            kap[k] += C[k, m] * cc / dz
            kap[m] += C[m, k] * cc / (-dz)
        return kap

    def mu_of(n):
        return np.array([2.0 * math.cos(2.0 * math.pi * l / n) - 1.0
                         for l in range(n)])

    def meeting_residual(z0, kappas, n):
        F = np.zeros(3)
        for i in range(3):
            li = li_all[i]
            s = INV_A10 + B_SM[i] * xs / TWO_PI
            s += float(n / TWO_PI * (li @ (xs - lp_s - z0))
                       + n / TWO_PI * li @ kappas)
            F[i] = s - 1.0 / A_Z[i]
        return F

    def desert_ok(z0, kappas, n):
        mu = mu_of(n)
        allz = [z0 + lp_s[k] + kappas[k] * m for k in range(nk) for m in mu]
        return min(allz) >= -1e-9 and max(allz) <= Z_MAX + 1e-9

    # ---- sign structure: what level repulsion can and cannot do ------------
    # kappa_k = sum_m C_km c^2 / (z_k - z_m): the sign of every term is
    # sign(z_k - z_m), z0-independent. The LOWEST class can only go down,
    # the HIGHEST only up. Map the forced signs.
    order = np.argsort(lp_s)
    forced = []
    for k in range(nk):
        above = sum(1 for m in range(nk) if m != k and lp_s[m] > lp_s[k])
        below = nk - 1 - above
        sgn = "negative only" if below == 0 else (
            "positive only" if above == 0 else "mixed")
        forced.append({"class": k, "lp": round(float(lp_s[k]), 4),
                       "rank": int(np.where(order == k)[0][0]),
                       "reachable_sign": sgn})

    # ---- optimise the meeting over (z0, c^2) for each n --------------------
    from scipy.optimize import minimize
    results = {}
    closed = []
    for n in range(3, 13):
        def obj(x):
            z0, c2 = x[0], np.abs(x[1:])
            kap = kappas_from(c2, z0)
            F = meeting_residual(z0, kap, n)
            pen = 0.0
            mu = mu_of(n)
            for k in range(nk):
                z = z0 + lp_s[k] + kap[k] * mu
                pen += 10.0 * (np.sum(np.maximum(0.0, -z))
                               + np.sum(np.maximum(0.0, z - Z_MAX)))
            return float(np.max(np.abs(F))) + pen

        best = None
        rng = np.random.default_rng(110 + n)
        for trial in range(40):
            x0 = np.concatenate([[rng.uniform(-2.0, 2.0)],
                                 rng.uniform(0.0, 3.0, K)])
            r = minimize(obj, x0, method="Nelder-Mead",
                         options={"maxiter": 4000, "xatol": 1e-10,
                                  "fatol": 1e-12})
            if best is None or r.fun < best.fun:
                best = r
        z0, c2 = best.x[0], np.abs(best.x[1:])
        kap = kappas_from(c2, z0)
        F = meeting_residual(z0, kap, n)
        ok = desert_ok(z0, kap, n)
        res = float(np.max(np.abs(F)))
        results[n] = {"residual": res, "desert_ok": bool(ok),
                      "z0": round(float(z0), 4),
                      "kappas": [round(float(v), 4) for v in kap]}
        if res < 1e-10 and ok:
            closed.append(n)

    # ---- verification of the reported point at the best n ------------------
    best_n = min(results, key=lambda n: results[n]["residual"])
    z0b = results[best_n]["z0"]
    kapb = np.array(results[best_n]["kappas"])
    Fb = meeting_residual(z0b, kapb, best_n)
    verify_residual = float(np.max(np.abs(Fb)))

    # ---- the V96 reference vector: is IT reachable? ------------------------
    v96 = json.load(open(os.path.join(ROOT, "experiments", "V96_kappa_s_bond",
                                      "results", "v96_kappa_s_bond.json")))
    v96_kappas = None
    for sol in v96["data"]["solutions"]:
        if "kappas" in sol:
            v96_kappas = np.array(sol["kappas"])
            break

    ok = (abs(verify_residual - results[best_n]["residual"]) < 1e-3
          and n_root_pairs > 0 and K >= 3)
    verdict = ("the meeting CLOSES under bond level repulsion"
               if closed else
               "the meeting does NOT close under bond level repulsion")
    status = (
        f"OFF-DIAGONAL BOND REALIZABILITY TEST: the bond graph has "
        f"{n_root_pairs} root-vector pairs and {K} class-pair channels "
        f"(degeneracies {int(C.sum() / 2)} total). Level repulsion forces "
        f"the sign of kappa_k by the energy ordering: "
        f"{sum(1 for f in forced if f['reachable_sign'] != 'mixed')} of "
        f"{nk} classes have a FORCED sign "
        f"({', '.join('%d:%s' % (f['class'], f['reachable_sign']) for f in forced)}). "
        f"{verdict}: best residual over n=3..12 is "
        f"{results[best_n]['residual']:.3e} at n={best_n} "
        f"(desert_ok={results[best_n]['desert_ok']}), verified to "
        f"{verify_residual:.3e}; the closed set is {closed}. "
        f"The six kappa_s of V96 remain a NAMED input of the bond layer: "
        f"one Hermitian bond operator with {K + 1} parameters either "
        f"reproduces the meeting (closed: {closed}) or the sign/size "
        f"structure obstructs it (floor "
        f"{results[best_n]['residual']:.3e}).")
    out = {"experiment": "V110_w4_offdiagonal_bond", "pass": bool(ok),
           "status": status,
           "data": {"theta_star": round(TH_STAR, 6),
                    "lp_classes": [round(float(v), 4) for v in lp_s],
                    "n_root_pairs": int(n_root_pairs),
                    "n_channels": int(K),
                    "channel_degeneracies": {f"{k}-{m}": int(C[k, m])
                                             for (k, m) in channels},
                    "forced_signs": forced,
                    "per_n": results,
                    "closed_n": closed,
                    "best_n": int(best_n),
                    "verify_residual": verify_residual,
                    "v96_kappas": (None if v96_kappas is None else
                                   [round(float(v), 4) for v in v96_kappas])}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v110_w4_offdiagonal_bond.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
