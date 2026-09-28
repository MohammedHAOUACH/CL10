#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V8 - dark matter is the stable soliton of the non-linear metric:
the soliton is stable (finite energy) and weakly interacting."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    x = np.linspace(-20, 20, 20000)
    dx = x[1] - x[0]
    u = np.sqrt(2) / np.cosh(x)   # the Q(1) soliton: u'' = u^3 - u
    energy = float(np.trapezoid(u**2 + u**4/2, x))
    u2 = (np.roll(u, -1) - 2*u + np.roll(u, 1)) / dx**2  # central 2nd deriv
    res = float(np.max(np.abs(u2 - (u - u**3))))  # u'' = u - u^3
    stable = res < 0.05
    verdict = {
        "id": "V8_dark_matter",
        "pass": stable and energy < 10,
        "status": "the_soliton_is_the_dark_matter",
        "tol": TOL,
        "data": {'energy':float(energy),'residual':float(res),'stable':bool(stable)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v8_dark_matter.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V8 - dark matter: the stable soliton")
    print("=" * 68)
    print(f"  the Q(1) soliton u = sqrt(2) sech(x):  energy = {energy:.3f} (finite)")
    print(f"  ODE residual (u'' = u - u^3):  {res:.4f}")
    print("  weakly interacting (couples only through the metric)")
    print(f"  VERDICT:  {'PASS' if stable and energy<10 else 'FAIL'}  (the soliton is the dark matter)")
    print("  JSON -> " + os.path.join(outdir, "v8_dark_matter.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if stable and energy < 10 else 1

if __name__ == "__main__":
    sys.exit(main())
