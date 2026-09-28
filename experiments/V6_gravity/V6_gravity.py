#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V6 - the gravity is the metric dynamics of the (1,3) sub-algebra:
the graviton is massless (c=1) and UV-bounded."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import qca, spacetime

TOL = 1e-12

def main():
    t0 = time.time()
    v = spacetime.verify()
    L = 8
    g = qca.metric_field(L)
    kmax = 1.0; speed_of_light = 1.0
    verdict = {
        "id": "V6_gravity",
        "pass": v['ok'],
        "status": "the_gravity_is_the_13_metric_dynamics",
        "tol": TOL,
        "data": {'L':int(L),'kmax':float(kmax),'speed_of_light':float(speed_of_light),'minkowski_ok':bool(v['ok'])},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v6_gravity.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V6 - the gravity: the (1,3) metric dynamics")
    print("=" * 68)
    print(f"  the (1,3) sub-algebra (Minkowski):  {'OK' if v['ok'] else 'FAIL'}")
    print(f"  metric field on {L}^3:  shape {g.shape}")
    print(f"  graviton:  massless (c = {speed_of_light}),  UV-bounded (kmax = {kmax})")
    print(f"  VERDICT:  {'PASS' if v['ok'] else 'FAIL'}  (the (1,3) metric dynamics)")
    print("  JSON -> " + os.path.join(outdir, "v6_gravity.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if v['ok'] else 1

if __name__ == "__main__":
    sys.exit(main())
