#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V48 — the **3- and 4-particle** QCA amplitudes (the V48 completion, the
**UV-finite** property, the O9 completion).

**The 0-parameter statement.**  The **n-particle** S-matrix is the
**Kronecker** product of the **1-particle** (``U (x) ... (x) U``, n times),
so it is **unitary** **iff** the **1-particle** is **unitary**:
``(U (x) ... (x) U)(U (x) ... (x) U)^dagger = (U U^dagger) (x) ... (x)
(U U^dagger) = I``.  The **3- and 4-particle** amplitudes are thus **finite**
(the **UV-finite** property: the n-step is **bounded** for all n).

**The 0-parameter checks.**  The **1- and 2-particle** unitarity is
**computed** **numerically** on the GPU (the **32x32** and **1024x1024**);
the **3- and 4-particle** (``32^3 = 32768`` and ``32^4 = 10^6``, too large
to store) follow **by** the **Kronecker** product property (the
**0-parameter** consequence).  The **only** input is ``g_10 = 0.2`` (from
the root by V38).

**The named residual** (O9): the **matching** of these amplitudes to the
**measured** cross-sections (the QCA lattice RG vs the perturbative QFT) —
the open part of the 't Hooft program.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import smatrix

TOL = 1e-8  # the 0-parameter identity tolerance


def main():
    t0 = time.time()
    r = smatrix.n_particle()

    ok = bool(
        r["n1_unitary_residual"] < TOL
        and r["n2_unitary_residual"] < TOL
        and r["n3_unitary_residual"] < TOL
        and r["n4_unitary_residual"] < TOL
        and r["n1_bounded"]
        and r["uv_finite"]
    )

    verdict = {
        "id": "V48_3particle",
        "tier": 3,
        "property": "C6",
        "residual": True,   # the matching to measured cross-sections (O9, the open part)
        "pass": ok,
        "status": "the_3and4particle_amplitudes_are_uvfinite_the_kronecker_product_property_the_nstep_is_bounded_for_all_n",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "n1_unitary_residual": r["n1_unitary_residual"],
            "n2_unitary_residual": r["n2_unitary_residual"],
            "n3_unitary_residual": r["n3_unitary_residual"],
            "n4_unitary_residual": r["n4_unitary_residual"],
            "n1_bounded": r["n1_bounded"],
            "uv_finite": r["uv_finite"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v48_3particle.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V48 - the 3- and 4-particle QCA amplitudes (the UV-finite property, O9)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] 1-particle unitarity:  residual = {r['n1_unitary_residual']:.2e} (0-parameter, computed)")
    print(f"  [check 2] 2-particle unitarity:  residual = {r['n2_unitary_residual']:.2e} (0-parameter, computed)")
    print(f"  [check 3] 3-particle unitarity:  residual = {r['n3_unitary_residual']:.2e} (32^3, by the Kronecker property)")
    print(f"  [check 4] 4-particle unitarity:  residual = {r['n4_unitary_residual']:.2e} (32^4, by the Kronecker property)")
    print(f"  [check 5] bounded:  max|S| <= 1  ({'YES' if r['n1_bounded'] else 'NO'})")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the 3- and 4-particle amplitudes are UV-finite: unitary by the Kronecker product property, the n-step is bounded for all n)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v48_3particle.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
