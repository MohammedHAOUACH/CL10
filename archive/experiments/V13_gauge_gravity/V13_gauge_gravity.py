#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V13 - matter, gauge, and gravity are in the ONE Cl(10): the (1,3)
causal (metric/gravity) and the 45 (gauge) are sub-algebras of the same 10."""

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
    gg = gauge.gauge_gravity_overlap()
    n_c = gg['n_causal']; n_g = gg['n_gauge']; n_gen = gg['n_generators']
    one_obj = gg['one_object'] and n_gen == 10 and n_g == 45 and n_c == 4
    verdict = {
        "id": "V13_gauge_gravity",
        "tier": 3,
        "property": "P3",
        "pass": one_obj,
        "residual": False,
        "status": "matter_gauge_gravity_in_one_object",
        "tol": TOL,
        "data": {'n_causal': int(n_c), 'n_gauge': int(n_g), 'n_generators': int(n_gen), 'one_object': bool(one_obj)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v13_gauge_gravity.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V13 - gauge-gravity unification (one object)  [tier 3, P3]")
    print("=" * 68)
    print(f"  the (1,3) causal (metric/gravity):  {n_c} of the {n_gen} generators")
    print(f"  the 45 (gauge):  {n_g} products h_a h_b")
    print(f"  both sub-algebras of the ONE {n_gen}-generator Cl(10):  {one_obj}")
    print(f"  VERDICT:  {'PASS (matter + gauge + gravity in one object)' if one_obj else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v13_gauge_gravity.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if one_obj else 1

if __name__ == "__main__":
    sys.exit(main())
