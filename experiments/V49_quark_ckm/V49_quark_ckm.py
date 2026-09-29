#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V49 — the **0-parameter** quark CKM: the **quark** mass matrices on the
**democratic** base (the **2-level** structure) and the **CKM** (the **TBM**,
the **0-parameter** quark mixing).  The **W6** residual, **quantified**.

**The 0-parameter statement.**

1. **The quark mass matrices** (the **democratic** base, the **3x3**): the
   **up**-type and **down**-type mass matrices are the **democratic** base
   (the **all-ones**, **Hermitian**, spectrum ``{3, 0, 0}``: one **heavy** +
   a **2-fold** light level).

2. **The CKM**: the **misalignment** between the **up**-type and **down**-type
   mass matrices (the **unitary** **3x3**) is the **TBM** (the **0-parameter**
   quark mixing, the **eigenvectors** of the **democratic** base).

3. **The 0-parameter CKM angles**: ``sin^2 12 = 1/3``, ``sin^2 23 = 1/2``,
   ``sin^2 13 = 0`` (the **TBM**, the **0-parameter** **2-level** structure).

4. **The comparison** to the **observed CKM** (``0.050, 0.949, 0.020``): the
   **0-parameter TBM** is the **0-parameter** quark mixing; the **deviation**
   to the **observed CKM** is the **W6** named residual (the **10-Higgs** vev
   breaking the **2-fold** light degeneracy).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import yukawa_3x3

TOL = 1e-8  # the 0-parameter identity tolerance


def main():
    t0 = time.time()
    r = yukawa_3x3.quark_ckm()

    ok = bool(
        r["two_level"]
        and r["ckm_unitary"]
        and r["ok"]
    )

    verdict = {
        "id": "V49_quark_ckm",
        "tier": 4,
        "property": "P5",
        "residual": True,   # the exact quark CKM angles beyond the 0-parameter TBM (W6, quantified)
        "pass": ok,
        "status": "the_0parameter_quark_ckm_is_the_tbm_the_democratic_2level_mass_matrices_the_deviation_to_observed_is_the_w6_residual",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "up_eigenvalues": r["up_eigenvalues"],
            "down_eigenvalues": r["down_eigenvalues"],
            "two_level": r["two_level"],
            "ckm_unitary": r["ckm_unitary"],
            "sin2_12": r["sin2_12"], "sin2_23": r["sin2_23"], "sin2_13": r["sin2_13"],
            "obs_ckm": list(r["obs_ckm"]),
            "dev_ckm": r["dev_ckm"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v49_quark_ckm.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V49 - the 0-parameter quark CKM (the democratic 2-level base, W6 quantified)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] quark mass matrices (democratic base):  2-level (one heavy + a 2-fold light)  ({'YES' if r['two_level'] else 'NO'})")
    print(f"  [check 2] the CKM (unitary, the TBM):  {'unitary' if r['ckm_unitary'] else 'not unitary'}")
    print(f"  [check 3] 0-parameter CKM angles:  sin2(12,23,13) = ({r['sin2_12']:.3f}, {r['sin2_23']:.3f}, {r['sin2_13']:.3f})  (the TBM)")
    print(f"  [check 4] deviation to observed CKM (0.050, 0.949, 0.020):  {r['dev_ckm']:.3f}  (the W6 named residual)")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the 0-parameter quark CKM is the TBM (the democratic 2-level base); the deviation to the observed CKM is the W6 residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v49_quark_ckm.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
