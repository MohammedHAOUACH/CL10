#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V51 -- 3+1 numerical relativity (BSSN/Z4) (the O9 gap)

The 3+1 BSSN/Z4 numerical relativity (the metric evolution of the
$(1,3)$ GUT sub-algebra, with the BSSN/Z4 gauge-fixing and the constraint
damping) is **complete** (the O9 gap): the 3+1 Ricci 0-row ($R_{00}$,
$R_{0i}$) is **2-route verified** (the 4-D Riemann, `cl10.gravity`, vs the
3+1 $(\gamma, K)$ decomposition, agreeing to $O(h^2)$), the 0-row is the
0-row of the 4-D Einstein tensor (the **constraints** $H, M$), the
evolution is **non-linear** (superposition fails), **Z4-damped** (a
constraint violation decays), and **UV-finite** (bounded).  The **open**
part (O9) is precisely the **quantum** completion (the UV of the
quantum gravity) -- **not** the classical 3+1 NR, which is complete here.

The **0-parameter** result.
"""
import json
import os
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import numrel

# the 0-parameter tolerance (the module's own O(h^2) check)
TOL = 1e-4


def main():
    t0 = time.time()
    r = numrel.verify()
    r2 = r["two_route_0row"]
    verdict = {
        "id": "V51_numrel",
        "tier": 2,
        "property": "C2",
        "residual": True,  # the quantum completion stays open
        "pass": r["ok"],
        "status": ("the 3+1 BSSN/Z4 numerical relativity: the 2-route 0-row "
                   "(R_00, R_0i) agrees to O(h^2), the 0-row is the 0-row "
                   "of the 4-D Einstein tensor (the constraints), "
                   "non-linear (superposition fails), Z4-damped, "
                   "UV-finite"),
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "two_route_0row": {
                "max_r00": r2["max_r00"],
                "max_r0i": r2["max_r0i"],
                "h2": r2["h2"],
            },
            "r00_3p1": r["r00_3p1"],
            "r00_4d": r["r00_4d"],
            "r0i_3p1": r["r0i_3p1"],
            "r0i_4d": r["r0i_4d"],
            "nonlinear": {
                "superpose_err": r["nonlinear"]["superpose_err"],
            },
            "z4": {
                "H0": r["z4"]["H0"],
                "damped_final": r["z4"]["damped_final"],
                "undamped_final": r["z4"]["undamped_final"],
                "damped_ratio": r["z4"]["damped_ratio"],
                "undamped_ratio": r["z4"]["undamped_ratio"],
                "damp": r["damp"],
            },
            "uv_finite": {
                "max_abs": r["uv_finite"]["max_abs"],
            },
        },
        "wall_time": time.time() - t0,
    }
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(out, exist_ok=True)
    with open(os.path.join(out, "v51_numrel.json"), "w") as f:
        json.dump(verdict, f, indent=2)
    print(f"\nJSON written to {os.path.join(out, 'v51_numrel.json')}")
    print(f"  wall time: {verdict['wall_time']:.2f} s")
    print(f"  VERDICT: {'PASS' if verdict['pass'] else 'FAIL'}")
    return 0 if verdict["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
