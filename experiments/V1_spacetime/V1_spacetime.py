#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V1 - the (1,3) is a sub-algebra of the GUT: time = h_1 (+), space =
i h_2..4 (-); the Minkowski relations hold 0-parameter."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import spacetime

TOL = 1e-12

def main():
    t0 = time.time()
    v = spacetime.verify()
    T, S = spacetime.causal_generators()
    t2=v['t2']; s2=v['s2']; ts=v['ts']; ss=v['ss']; err_max=max(t2,s2,ts,ss)
    verdict = {
        "id": "V1_spacetime",
        "pass": v['ok'],
        "status": "the_13_is_inside_the_gut",
        "tol": TOL,
        "data": {'t2':float(t2),'s2':float(s2),'ts':float(ts),'ss':float(ss),'err':float(err_max)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v1_spacetime.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V1 - the (1,3) spacetime INSIDE the GUT")
    print("=" * 68)
    print(f"  time T^2 = +I    :  {t2:.2e}  (T = h_1, GUT)")
    print(f"  space S^2 = -I   :  {s2:.2e}  (S = i h_2..4, GUT)")
    print(f"  anticommutators  :  {max(ts,ss):.2e}  (Minkowski (1,3))")
    print(f"  VERDICT:  {'PASS' if v['ok'] else 'FAIL'}  (the light cone is a GUT sub-algebra)")
    print("  JSON -> " + os.path.join(outdir, "v1_spacetime.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if v['ok'] else 1

if __name__ == "__main__":
    sys.exit(main())
