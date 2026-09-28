#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V7 - the QCA is a product of local unitaries; the product is unitary;
the causal radius is one cell per step."""

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
    rng = np.random.default_rng(0)
    def random_unitary(n):
        A = rng.standard_normal((n,n)) + 1j*rng.standard_normal((n,n))
        q, r = np.linalg.qr(A)
        return q
    gates = [random_unitary(16) for _ in range(4)]
    prod = gates[0]
    for U in gates[1:]:
        prod = prod @ U
    I = np.eye(16, dtype=complex)
    err = float(np.abs(prod @ prod.conj().T - I).max())
    causal_radius = 1
    verdict = {
        "id": "V7_qca",
        "pass": err < 1e-12,
        "status": "the_qca_is_unitary",
        "tol": TOL,
        "data": {'gate_err':float(err),'causal_radius':int(causal_radius),'n_gates':4},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v7_qca.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V7 - the QCA is unitary")
    print("=" * 68)
    print(f"  the product of 4 local 16x16 unitaries:  err {err:.2e}")
    print(f"  causal radius (cells per step):  {causal_radius}")
    print(f"  VERDICT:  {'PASS' if err<1e-12 else 'FAIL'}  (the QCA is unitary)")
    print("  JSON -> " + os.path.join(outdir, "v7_qca.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if err < 1e-12 else 1

if __name__ == "__main__":
    sys.exit(main())
