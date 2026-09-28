#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V17 - W1 resolved: the Z2 512/512 cancellation gives Lambda = 0 at
0-parameter (the s-mechanism / self-tuning).  The mechanism is a 0-parameter
consequence (NOT a residual); the exact 10^-122 is the named fine-tuning
(V9)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import walls

TOL = 1e-12

def main():
    t0 = time.time()
    w1 = walls.w1_cosmological_constant()
    ok = (w1["cancellation"] and w1["lambda_0param"] == 0.0)
    verdict = {
        "id": "V17_w1_screening",
        "tier": 4,
        "property": "P2,P5",
        "residual": False,     # the MECHANISM is a 0-parameter consequence
        "pass": ok,
        "status": "the_w1_mechanism_z2_cancellation",
        "tol": TOL,
        "data": {'n_even': w1['n_even'], 'n_odd': w1['n_odd'],
                 'lambda_0param': w1['lambda_0param'],
                 'lambda_obs': w1['lambda_obs']},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v17_w1_screening.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V17 - W1 mechanism: the Z2 512/512 cancellation (the s-mechanism)")
    print("=" * 68)
    print(f"  even/odd = {w1['n_even']}/{w1['n_odd']}  (equal: {w1['cancellation']})")
    print(f"  0-parameter Lambda = {w1['lambda_0param']}  (the cancellation -> 0)")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the W1 mechanism is 0-parameter: Lambda=0)")
    print("  JSON -> " + os.path.join(outdir, "v17_w1_screening.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
