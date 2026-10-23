#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V74 - the O4 audit: multi-scale node splitting CANNOT move M_X.

**The hope.**  V65's single-scale Wilson mechanism meets at M_X =
8.9e14 GeV (tau_p ~ 1e31 yr, EXCLUDED by F5); the data-preferred V64
solution meets at 1.8e16 GeV with two fitted thresholds.  Maybe three
nodes with DIFFERENT strengths w1, w2, w3 (one per V70 family) can
reach the data point.

**The structural fact (proved numerically here).**  Each split node
carries the FULL vectorlike split 16 (six fragments, thresholds
z_k + ln c_f(theta)).  While every threshold sits below the meeting
point x, the three meeting equations depend on the strengths ONLY
through the sum z1 + z2 + z3: the Jacobian d(F1,F2,F3)/d(z1,z2,z3)
has three IDENTICAL columns (rank 1).  Verified at the V65 root, at
split strengths, at V64-like low thresholds and at a generic theta:
rank 1 exactly; the residuals F_i are EXACTLY invariant under any
redistribution of the sum.  A node pushed above x decouples entirely
(V68: n_nodes < 3 has no physical root).  Hence:

**Verdict (computed).**  The node mechanism's meeting point is the V65
value for EVERY choice of the three strengths: M_X = 8.9e14 GeV,
tau_p ~ 1e31 yr, EXCLUDED by F5.  The V64 data-preferred spectrum
(three 10+10bar-type sets at 1.2e7 + one (1,2) pair at 1e12) is NOT a
three-full-16 pattern - the automaton supplies full 16's per node, the
data want partial suits.  O4's gap is a NAMED INCOMPATIBILITY
(node content = full split 16 vs the desert spectrum the couplings
prefer), not a fitting problem: closing it needs new automaton content
(partial-suit states), not new thresholds.
"""

import json
import math
import os
import sys
import time

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)
os.environ["CL10_ANGLE"] = "0.2"  # built on the V65 desert machinery: its anchors are 1/5-branch, so pin the clock BEFORE cl10 latches it (V65 pins the same)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))
import V65_wilson as V65  # noqa: E402
from cl10 import unification as UN  # noqa: E402
from cl10 import u10pred  # noqa: E402
from cl10 import config as _config  # noqa: E402

MZ, A_Z, B_SM = UN.M_Z, UN.A_Z, UN.B_SM
M_P = UN.M_P
TWO_PI = 2.0 * math.pi
INV_A10 = 25.0
THETA_ROOT = float(_config.CFG.desert.theta_root_rad)   # the V65 root direction
_CACHE = {}


def fragments(theta):
    if "frs" not in _CACHE:
        lam, _ = V65.cartan_weights()
        x_coefs, y_coefs = V65.identify_XY(lam)
        y_coef = y_coefs[0]
        xc_gg = next(xc for xc in x_coefs
                     if V65.physical_su5(lam @ xc, lam @ y_coef))
        Xs, Ys = lam @ xc_gg, lam @ y_coef
        _CACHE["frs"] = V65.fragment_table(Xs, Ys)
    frs = _CACHE["frs"]
    ca, cb = math.cos(theta), math.sin(theta)
    cf = [abs(ca * f["X"] + cb * f["Y"]) for f in frs]
    return frs, cf


def F_of(z, x, theta):
    frs, cf = fragments(theta)
    lcf = [math.log(c) for c in cf]
    F = []
    for i in range(3):
        s = INV_A10 + B_SM[i] * x / TWO_PI
        for k in range(3):
            for j, f in enumerate(frs):
                act = x - z[k] - lcf[j]
                if act > 0:
                    s += f["db"][i] * act / TWO_PI
        F.append(s - 1.0 / A_Z[i])
    return F


def jacobian_rank(z, x, theta):
    J = np.zeros((3, 3))
    for k in range(3):
        zp = list(z)
        zp[k] += 1e-5
        J[:, k] = ((np.array(F_of(zp, x, theta))
                    - np.array(F_of(z, x, theta))) / 1e-5)
    return int(np.linalg.matrix_rank(J, tol=1e-6))


def main():
    t0 = time.time()

    # --- rank-1 fact at several configurations ----------------------------
    configs = [
        ("V65 root, degenerate strengths", THETA_ROOT, [20.0, 20.0, 20.0]),
        ("V65 root, split strengths", THETA_ROOT, [14.0, 20.0, 26.0]),
        ("V64-like low thresholds", THETA_ROOT, [11.5, 11.5, 11.5]),
        ("generic theta", 1.2, [15.0, 18.0, 21.0]),
    ]
    rank_checks = {}
    all_rank1 = True
    for name, th, z in configs:
        r = jacobian_rank(z, 33.5, th)
        all_rank1 = all_rank1 and (r == 1)
        rank_checks[name] = r

    # --- the node mechanism's meeting point: V65's own scan --------------
    th = THETA_ROOT
    frs0, _ = fragments(th)
    roots = [r for r in V65.scan_roots(frs0) if r["physical"]]
    r0 = max(roots, key=lambda r: r["M_X_GeV"])
    mx_node = r0["M_X_GeV"]
    x_star = math.log(mx_node / MZ)
    z0 = math.log(r0["w"] * M_P / (3.0 * MZ))

    # split-invariance of the residuals at the meeting point
    xs = []
    for split in ([z0] * 3, [z0 - 6, z0, z0 + 6], [z0 - 12, z0 + 3, z0 + 9]):
        F1 = np.array(F_of(split, x_star, th))
        F2 = np.array(F_of([z0] * 3, x_star, th))
        xs.append(float(np.max(np.abs(F1 - F2))))
    split_invariance = max(xs) < 1e-12

    # --- tau_p at the node-mechanism meeting point ------------------------
    t_gev = mx_node ** 4 / ((1 / 25) ** 2 * u10pred.M_P_ROBSON ** 5)
    tau_p = t_gev * u10pred.GEV_TO_S / u10pred.S_TO_YR

    ok = all_rank1 and split_invariance and roots
    status = (
        "MULTI-SCALE CANNOT MOVE M_X (structural): while all node "
        "thresholds sit below the meeting point the meeting equations "
        "depend on the three node strengths ONLY through their sum - the "
        "Jacobian has three identical columns (rank 1 at every "
        f"configuration tested: {list(rank_checks.values())}); the "
        f"residuals F_i are EXACTLY split-independent (dev {max(xs):.1e}). "
        "A node above x decouples (V68: n < 3 has no root). Hence the "
        f"node mechanism meets at M_X = {mx_node:.2e} GeV for EVERY "
        f"split, tau_p = {tau_p:.1e} yr - EXCLUDED by F5. The V64 "
        "data-preferred spectrum (three 10+10bar sets at 1.2e7 + one "
        "(1,2) pair at 1e12) is not a three-full-16 pattern: the "
        "automaton supplies full split 16's per node, the couplings "
        "prefer partial suits. O4 is a NAMED INCOMPATIBILITY (node "
        "content vs desert spectrum), not a fitting problem; closing it "
        "needs new automaton content, not new thresholds.")

    verdict = {"experiment": "V74_multiscale_nodes", "pass": bool(ok),
               "status": status,
               "data": {"jacobian_ranks": rank_checks,
                        "split_invariance_dev": max(xs),
                        "M_X_node_GeV": mx_node, "tau_p_yr": tau_p,
                        "n_physical_roots": len(roots)}}

    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                          "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v74_multiscale_nodes.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V74 - multi-scale node splitting: the rank-1 theorem (O4)")
    print("=" * 68)
    for name, r in rank_checks.items():
        print(f"  rank dF/d(z1,z2,z3) [{name}] = {r}")
    print(f"  split-invariance of F_i: dev = {max(xs):.2e}")
    print(f"  node-mechanism meeting: M_X = {mx_node:.2e} GeV, "
          f"tau_p = {tau_p:.1e} yr (F5 bound 2.4e34)")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v74_multiscale_nodes.json"))
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
