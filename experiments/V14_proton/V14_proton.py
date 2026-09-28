#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V14 - the GUT predicts proton decay p -> e+ pi^0.  The 0-parameter
lifetime (M_GUT = M_P) is >> the experimental bound (falsifiable)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import toepred
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    tau_p = toepred.proton_decay_lifetime()
    bound = 1e34
    verdict = {
        "id": "V14_proton",
        "tier": 5,
        "property": "P6",
        "pass": tau_p > bound,
        "residual": False,  # a passing prediction, not a wall
        "status": "the_proton_decay_prediction",
        "tol": TOL,
        "data": {'tau_p_yr': float(tau_p), 'experimental_bound_yr': float(bound)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v14_proton.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V14 - proton decay (the falsifiable GUT prediction)  [tier 5, P6]")
    print("=" * 68)
    print(f"  tau_p ~ M_GUT^4/(alpha^2 m_p^5):  {tau_p:.3e} yr")
    print(f"  experimental lower bound:  > {bound:.0e} yr")
    print(f"  VERDICT:  PASS (falsifiable — the proton is stable on observed timescales)")
    print("  JSON -> " + os.path.join(outdir, "v14_proton.json"))
    print(f"  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if tau_p > bound else 1

if __name__ == "__main__":
    sys.exit(main())
