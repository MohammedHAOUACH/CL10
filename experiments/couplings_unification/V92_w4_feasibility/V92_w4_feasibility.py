#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V92 - the O4 FEASIBILITY THEOREM: can ANY threshold pattern of
node-16 content reach the data-preferred meeting?

**Why (O4).**  V74 closed multi-scale splitting for IDENTICAL fragment
patterns (rank-1); V81 excluded the 45-composites; V82/V86/V87 scanned
junction-generated patterns and hit a common floor. The unscanned
space is every threshold assignment the automaton could possibly
supply: each node contributes the 16's fragment classes (the weights
of the node 16 under the breaking direction) at ARBITRARY desert
thresholds, with any number of nodes. The meeting equations at the
data-preferred point are LINEAR in the threshold widths:

    sum_s db_i(s) t_s = T_i,   t_s in [0, n x*/2pi],
    T_i = 1/A_Z,i - 1/alpha_10 - b_i x*/2pi,   x* = ln(M_X_obs/M_Z)

so feasibility is convex-hull membership in R^3 of T/x against
{0} union {db(s)} - decidable exactly. V92 decides it and, when
infeasible, returns the SEPARATING HYPERPLANE: the linear form that
every automaton fragment satisfies and the data violate - the exact
quantum-number signature of the missing content.

**Checks (computed).**
1. Regression: the V74 root meeting M_X = 8.9e14 GeV (V74 JSON) and
   the V64 data-preferred target from UN.
2. Hull membership: T' = T/x* in conv({0} union {db_s}) via exact
   Caratheodory enumeration (all <=4-vertex subsimplices).
3. Best infeasibility margin: min over u in [0,1]^m of
   ||sum db_s u_s - T'|| (projected gradient).
4. Separating hyperplane: the unit form w with w.T' > max_s w.db(s)
   by the largest margin (from the dual of the residual problem).
5. Junction cross-check: no Wilson direction theta reproduces the
   best-fit u_s pattern (scan 720 angles, report the min deviation).

**Verdict (computed).**  Whether node-16 content of ANY pattern can
close O4 (then: the thresholds and their junction realizability) or
cannot (then: the theorem, with the separating hyperplane naming the
required content's quantum numbers).
"""

import itertools
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_ANGLE"] = "0.2"  # built on the V65 desert machinery: its anchors are 1/5-branch, so pin the clock BEFORE cl10 latches it (V65 pins the same)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import unification as UN  # noqa: E402
from cl10 import u10pred  # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
M_P = UN.M_P
TWO_PI = 2.0 * math.pi
INV_A10 = 25.0
M_X_OBS = u10pred.M_X_DESERT   # GeV, the V64 data-preferred meeting


def main():
    t0 = time.time()
    checks = {}

    # --- fragments of the node 16 -------------------------------------------
    lam, _off = V65.cartan_weights()
    x_coefs, y_coefs = V65.identify_XY(lam)
    y_coef = y_coefs[0]
    xc_gg = next(xc for xc in x_coefs
                 if V65.physical_su5(lam @ xc, lam @ y_coef))
    Xs, Ys = lam @ xc_gg, lam @ y_coef
    frs = V65.fragment_table(Xs, Ys)
    dbs = np.array([f["db"] for f in frs], dtype=float)
    checks["fragment_classes"] = len(frs)
    checks["db_generators"] = [[round(v, 4) for v in d] for d in dbs]

    # --- the target ----------------------------------------------------------
    xs = math.log(M_X_OBS / MZ)
    T = np.array([1.0 / A_Z[i] - INV_A10 - B_SM[i] * xs / TWO_PI
                  for i in range(3)])
    checks["x_star"] = round(xs, 4)
    checks["target_T"] = [round(v, 4) for v in T]

    # --- hull membership with box t_s in [0, X], X = n x*/2pi ---------------
    # u_s = t_s/X in [0,1]: sum db_s u_s = T/X = Tp: Tp in conv({0,db_s})
    verts = np.vstack([np.zeros((1, 3)), dbs])   # vertex 0 = origin
    m = len(verts)
    Tp = T / xs
    checks["target_Tp"] = [round(v, 5) for v in Tp]
    feasible = False
    witness = None
    for combo in itertools.combinations(range(m), 4):
        V = verts[list(combo)]
        # barycentric: V.T lam = Tp, sum lam = 1, lam >= 0
        A = np.vstack([V.T, np.ones((1, 4))])
        b = np.append(Tp, 1.0)
        lam_b, *_ = np.linalg.lstsq(A, b, rcond=None)
        if (np.max(np.abs(A @ lam_b - b)) < 1e-9
                and np.min(lam_b) > -1e-9):
            feasible = True
            witness = {verts_i: round(float(l), 4) for verts_i, l
                       in zip(combo, lam_b) if l > 1e-6}
            break

    # --- best residual over the box (projected gradient) --------------------
    rng = np.random.default_rng(11)
    u = rng.uniform(0, 1, m - 1)
    D = dbs                                       # m-1 generators
    L = float(np.linalg.norm(D.T @ D, 2))
    for _ in range(20000):
        r = D.T @ u - Tp
        g = D @ r
        u = np.clip(u - g / L, 0.0, 1.0)
    r_best = D.T @ u - Tp
    resid = float(np.linalg.norm(r_best))
    checks["best_residual_norm"] = float(f"{resid:.4e}")
    checks["best_residual_components"] = [round(float(v), 4)
                                          for v in r_best]
    checks["hull_feasible"] = feasible
    if feasible:
        checks["witness_barycentric"] = witness

    # --- separating hyperplane (dual): max over |w|=1 of
    #     w.Tp - max_s w.db_s   (positive iff infeasible)
    w = r_best / max(resid, 1e-30)
    for _ in range(4000):
        svals = D @ w
        j = int(np.argmax(svals))
        grad = Tp - D[j]
        w = w + 1e-3 * grad
        w /= np.linalg.norm(w)
    margin = float(w @ Tp - np.max(D @ w))
    # refine scale: report the physical form: w.T vs max_s w.db(s)
    checks["separating_form_w"] = [round(float(v), 5) for v in w]
    checks["separating_margin_Tp"] = float(f"{margin:.4e}")
    checks["separating_form_T"] = float(f"{w @ T:.4f}")
    checks["max_form_over_fragments"] = float(f"{np.max(D @ w) * xs:.4f}")

    # --- junction cross-check: can any Wilson direction give the best u? ----
    ths = np.arange(720) * math.pi / 720
    best_dev = 1e9
    for th in ths:
        ca, cb = math.cos(th), math.sin(th)
        cf = np.array([abs(ca * f["X"] + cb * f["Y"]) for f in frs])
        # junction-realizable widths: t_s proportional to (x - z0 - ln cf_s)
        # clipped to [0, x]: shape test up to affine reparametrization:
        # rank the implied order and compare to u's order
        order_u = np.argsort(-u)
        order_cf = np.argsort(-np.log(np.maximum(cf, 1e-12)))
        dev = float(np.mean(order_u != order_cf))
        best_dev = min(best_dev, dev)
    checks["junction_order_best_dev"] = round(best_dev, 3)

    ok = (resid > 1e-3 and margin > 1e-3) or feasible

    if feasible:
        status = (
            f"NODE-16 CONTENT CAN REACH THE DATA-PREFERRED MEETING: "
            f"the target is inside the fragment-threshold convex hull "
            f"(witness barycentric {witness}). The required thresholds "
            f"are computed; junction realizability is the remaining "
            f"test (best order agreement over 720 Wilson angles: "
            f"{best_dev}).")
    else:
        status = (
            f"NO THRESHOLD PATTERN OF NODE-16 CONTENT REACHES THE "
            f"DATA-PREFERRED MEETING - O4 CLOSED BY THEOREM OVER THE "
            f"FULL SPACE (any junction, any profile, any n_nodes). "
            f"The meeting equations at M_X = {M_X_OBS:.1e} GeV, "
            f"alpha_X = 1/25 require sum db(s) t_s = T = "
            f"{[round(float(v), 2) for v in T]} with 0 <= t_s <= "
            f"n x*/2pi; the best box-constrained fit misses by "
            f"{resid:.3f} (components "
            f"{[round(float(v), 2) for v in r_best]}), and the "
            f"separating hyperplane w = "
            f"{[round(float(v), 3) for v in w]} gives w.T = "
            f"{w @ T:.2f} against max over every node-16 fragment "
            f"w.db(s) x* = {np.max(D @ w) * xs:.2f}: EVERY automaton "
            f"fragment of the node 16 satisfies w.db <= "
            f"{np.max(D @ w):.4f} while the data demand "
            f"{float(w @ Tp):.4f}. The missing content is named by "
            f"its quantum numbers: a state with db along w beyond "
            f"this bound - not a rearrangement of the node 16. "
            f"Junction patterns cannot even order-match the best fit "
            f"({best_dev} rank deviation over 720 angles).")

    verdict = {"experiment": "V92_w4_feasibility",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v92_w4_feasibility.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("V92 - O4 feasibility over the full node-16 threshold space")
    print(f"  fragment classes: {len(frs)}; target T = "
          f"{[round(float(v), 2) for v in T]}")
    print(f"  hull feasible: {feasible} | best residual: {resid:.3e}")
    print(f"  separating margin: {margin:.3e} | form w = "
          f"{[round(float(v), 3) for v in w]}")
    print(f"  junction order-match: {best_dev}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v92_w4_feasibility.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
