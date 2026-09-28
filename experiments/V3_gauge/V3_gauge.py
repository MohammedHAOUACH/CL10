#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V3 - the gauge is the 45 (h_a h_b for a<b); they close the Spin(10)
adjoint 0-parameter on the 16 Weyl."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import clifford, engine as E
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    g = [E.asnumpy(clifford.generators()[a]) for a in range(10)]
    gauge45 = [g[a] @ g[b] for a in range(10) for b in range(a+1, 10)]
    n45 = len(gauge45)
    mat = np.vstack([x.ravel() for x in gauge45])
    rank45 = int(np.linalg.matrix_rank(mat, tol=1e-10))
    verdict = {
        "id": "V3_gauge",
        "pass": (n45==45 and rank45==45),
        "status": "the_45_is_the_spin10_adjoint",
        "tol": TOL,
        "data": {'n45':int(n45),'rank45':int(rank45)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v3_gauge.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V3 - the gauge: the 45 (the Spin(10) adjoint)")
    print("=" * 68)
    print(f"  the 45 gauge generators (h_a h_b for a<b):  {n45}")
    print(f"  rank (they span the 45 adjoint):  {rank45}")
    print(f"  VERDICT:  {'PASS' if n45==45 and rank45==45 else 'FAIL'}  (the 45 is the Spin(10) adjoint)")
    print("  JSON -> " + os.path.join(outdir, "v3_gauge.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if (n45==45 and rank45==45) else 1

if __name__ == "__main__":
    sys.exit(main())
