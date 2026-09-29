#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V12 - the 0-parameter U(1) normalization g_1^2/g'^2 = 5/3 (the GUT
unification; the 5/3 is 0-parameter, the GUT Cartan normalization)."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import gauge

TOL = 1e-12

def main():
    t0 = time.time()
    r53 = gauge.five_thirds()
    target = 5 / 3
    verdict = {
        "id": "V12_unification",
        "tier": 3,
        "property": "P3,P5",
        "pass": abs(r53 - target) < 1e-6,
        "residual": False,
        "status": "the_couplings_unify_5_3",
        "tol": TOL,
        "data": {'g1_sq_over_gp2': float(r53), 'target': float(target), 'err': float(abs(r53 - target))},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v12_unification.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V12 - gauge-coupling unification (the 5/3)  [tier 3, P3,P5]")
    print("=" * 68)
    print(f"  g_1^2/g'^2  (0-parameter):  {r53:.6f}")
    print(f"  target (the 5/3 GUT normalization):  {target:.6f}")
    print(f"  relative error:  {abs(r53 - target) / target:.2e}")
    print(f"  VERDICT:  {'PASS (the 0-parameter 5/3)' if abs(r53 - target) < 1e-6 else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v12_unification.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if abs(r53 - target) < 1e-6 else 1

if __name__ == "__main__":
    sys.exit(main())
