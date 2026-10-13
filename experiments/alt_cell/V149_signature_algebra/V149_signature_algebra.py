#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V149 - does the EUCLIDEAN signature Cl(10,0) cause the residuals? (alt-cell, W-signature)

**The question (the user's hypothesis, axis 1).**  The root of CL10 is
Cl(10,0): ten Hermitian generators, h_a^2 = +I.  The hypothesis: this
signature is what forces the open problems.  The census: for every
signature Cl(10-q, q) (the first q generators multiplied by i, then
anti-Hermitian, square -I), recompute the three algebraic facts the
construction needs:

1. the 45 plane generators T_ab = i h_a h_b: Hermitian iff the two
   generators share their square-class, anti-Hermitian otherwise;
   T_ab^2 = eps_a eps_b I;
2. the unitarity of the 45-plane product Pi [cos g - i sin g T_ab] at the
   closure-quantized angles g = pi/2, pi/4 (and the implemented g_10 = 1/5)
   - NO new parameter is introduced, the same implemented form is used;
3. the algebraic closure (V131 method): the smallest exact period of the
   product, order <= 64.

**The expected physics.**  A mixed signature makes some T_ab
anti-Hermitian with T^2 = -I: the implemented factor cos g - i sin g T is
then Hermitian but NOT unitary, so the product is non-unitary (residual
O(1)) - non-compact Spin(p,q) has no finite-dimensional unitary
representation of the boosts.  q = 10 (all minus) keeps every T Hermitian
with T^2 = +I: the product stays unitary.  Both outcomes are results.

Pre-registered: no continuous knob is scanned; tier 4 (supporting).
"""

import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
sys.path.insert(0, ROOT)

import numpy as np

from cl10 import altcell as A
from cl10 import supercell as SC

I32 = np.eye(32, dtype=np.complex128)


def plane_class(q):
    """(n_herm, n_anti, n_T2p, n_T2m) over the 45 planes of the signature."""
    n_h = n_a = n_p = n_m = 0
    for T in A.plane_generators_sig(q):
        herm = float(np.abs(T - T.conj().T).max()) < 1e-10
        anti = float(np.abs(T + T.conj().T).max()) < 1e-10
        assert herm or anti, "T_ab neither Hermitian nor anti-Hermitian"
        n_h += herm
        n_a += anti
        n_p += float(np.abs(T @ T - I32).max()) < 1e-10
        n_m += float(np.abs(T @ T + I32).max()) < 1e-10
    return n_h, n_a, n_p, n_m


def main():
    t0 = time.time()
    rows = []
    for q in range(11):
        planes = A.plane_generators_sig(q)
        n_h, n_a, n_p, n_m = plane_class(q)
        row = {"signature": f"Cl(10-{q},{q})", "q": q,
               "T_hermitian": n_h, "T_antihermitian": n_a,
               "T2_plus": n_p, "T2_minus": n_m}
        for gname, g in (("pi/2", math.pi / 2), ("pi/4", math.pi / 4),
                         ("g10", SC.G10)):
            U = A.plane_product(planes, g)
            ue = A.unitary_err(U)
            r = {"unitarity_err": ue}
            if gname in ("pi/2", "pi/4"):
                N, err = A.exact_order(U, Nmax=64)
                r["closure_order"] = N
                r["closure_err"] = err
            row[gname] = r
        rows.append(row)

    base = rows[0]
    ok_base = (base["T_hermitian"] == 45 and base["T2_plus"] == 45
               and base["pi/2"]["unitarity_err"] < 1e-12
               and base["pi/4"]["unitarity_err"] < 1e-12
               and base["pi/2"]["closure_order"] == 4
               and base["pi/4"]["closure_order"] == 8)
    nonunitary = [r["signature"] for r in rows[1:]
                  if r["pi/2"]["unitarity_err"] > 1e-9
                  or r["pi/4"]["unitarity_err"] > 1e-9]
    unitary_all = [r["signature"] for r in rows
                   if r["pi/2"]["unitarity_err"] < 1e-12
                   and r["pi/4"]["unitarity_err"] < 1e-12]
    ok = ok_base and len(nonunitary) >= 1 and all(
        r["pi/2"]["unitarity_err"] is not None for r in rows)

    verdict = {
        "id": "V149_signature_algebra",
        "tier": 4,
        "property": "the signature census: which Cl(10-q,q) admit an exactly "
                    "unitary 45-plane product with zero new parameters",
        "residual": None,
        "pass": bool(ok),
        "status": f"SIGNATURE ALGEBRA CENSUS. q=0 baseline reproduces (45 "
                  f"Hermitian planes, T^2=I, unitarity "
                  f"{base['pi/2']['unitarity_err']:.1e}, orders "
                  f"{base['pi/2']['closure_order']}/{base['pi/4']['closure_order']}). "
                  f"Non-unitary products at: {nonunitary}. Unitary-product "
                  f"signatures: {unitary_all}. Mixed signatures break the "
                  f"compact plane-rotation form (anti-Hermitian T, T^2=-I): "
                  f"the implemented factor is not unitary.",
        "data": {"rows": rows, "unitary_signatures": unitary_all,
                 "nonunitary_signatures": nonunitary},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v149_signature_algebra.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V149 - the signature census (Cl(10-q,q), q = 0..10)")
    print("=" * 68)
    for r in rows:
        print(f"  {r['signature']:>11}  T herm {r['T_hermitian']:2d} anti "
              f"{r['T_antihermitian']:2d}  T^2 +I {r['T2_plus']:2d} -I "
              f"{r['T2_minus']:2d}  unit(pi/2) {r['pi/2']['unitarity_err']:.2e}"
              f"  unit(pi/4) {r['pi/4']['unitarity_err']:.2e}  orders "
              f"{r['pi/2']['closure_order']}/{r['pi/4']['closure_order']}")
    print(f"  unitary-product signatures: {unitary_all}")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v149_signature_algebra.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
