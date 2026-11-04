#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V0 - the root is Cl(10): the 10 generators anticommute, square to I;
the 2**10=1024 basis; the 16/16-bar split."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import clifford

TOL = 1e-12

def main():
    t0 = time.time()
    v = clifford.verify()
    b = clifford.basis()
    p, m = clifford.weyl_spinor()
    n_gen = v['n']; n_basis = len(b); n_weyl = p.shape[1]; n_16bar = m.shape[1]
    err_max = max(v['square_to_I'], v['anticommute'])
    verdict = {
        "id": "V0_root",
        "tier": 0,
        "property": "P0",
        "residual": False,
        "pass": (v['ok'] and n_gen == 10 and n_basis == 1024 and n_weyl == 16 and n_16bar == 16),
        "status": "the_cl10_is_the_root",
        "tol": TOL,
        "data": {'n_gen': int(n_gen), 'n_basis': int(n_basis), 'n_weyl': int(n_weyl), 'n_16bar': int(n_16bar), 'err': float(err_max)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v0_root.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V0 - the root: Cl(10) (10 generators)")
    print("=" * 68)
    print(f"  10 generators:  {n_gen}  (err {err_max:.2e})")
    print(f"  algebra basis:  {n_basis}  (2**10 = 1024)")
    print(f"  Weyl 16 / 16-bar:  {n_weyl} / {n_16bar}")
    print(f"  VERDICT:  {'PASS' if v['ok'] and n_weyl==16 else 'FAIL'}  (the Cl(10) is the root)")
    print("  JSON -> " + os.path.join(outdir, "v0_root.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if (v['ok'] and n_gen == 10 and n_basis == 1024 and n_weyl == 16 and n_16bar == 16) else 1

if __name__ == "__main__":
    sys.exit(main())
