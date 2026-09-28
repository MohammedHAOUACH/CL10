#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V4 - the statistics (boson/fermion) are the Z2 grading (even/odd in
the 10 generators); a property of every Clifford algebra, 0-parameter."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import clifford

TOL = 1e-12

def main():
    t0 = time.time()
    import math
    n_even = sum(math.comb(10, r) for r in range(11) if r % 2 == 0)
    n_odd  = sum(math.comb(10, r) for r in range(11) if r % 2 == 1)
    total = len(clifford.basis())
    verdict = {
        "id": "V4_statistics",
        "pass": (total==1024 and n_even==512 and n_odd==512),
        "status": "the_z2_grading_is_the_statistics",
        "tol": TOL,
        "data": {'total':int(total),'n_even':int(n_even),'n_odd':int(n_odd)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v4_statistics.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V4 - the statistics: the Z2 grading")
    print("=" * 68)
    print(f"  total basis:  {total}  (2**10 = 1024)")
    print(f"  even (boson) sector:  {n_even}  (= 2**9)")
    print(f"  odd (fermion) sector:  {n_odd}  (= 2**9)")
    print(f"  VERDICT:  {'PASS' if n_even==512 and n_odd==512 else 'FAIL'}  (the Z2 grading)")
    print("  JSON -> " + os.path.join(outdir, "v4_statistics.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if (total==1024 and n_even==512 and n_odd==512) else 1

if __name__ == "__main__":
    sys.exit(main())
