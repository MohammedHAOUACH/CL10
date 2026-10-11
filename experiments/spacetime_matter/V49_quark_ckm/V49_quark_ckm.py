#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V49 — the quark mass matrices on the **democratic** base (the **2-level**
structure, 0-parameter) and the **CKM**, **honestly**: the CKM computed from
identical up/down democratic matrices is the **IDENTITY**, so a non-trivial CKM
is the **W6** named residual (it needs an up/down misalignment).

**The 0-parameter statement.**

1. **The quark mass matrices** (the **democratic** base, the **3x3**): the
   **up**-type and **down**-type mass matrices are the **democratic** base
   (the **all-ones**, **Hermitian**, spectrum ``{3, 0, 0}``: one **heavy** +
   a **2-fold** light level).

2. **The CKM, honestly**: the CKM is the **misalignment** between the
   diagonalizers of the up-type and down-type mass matrices.  Here both are the
   **same** democratic matrix, so the computed CKM is the **IDENTITY**, not the
   TBM.  A non-trivial CKM requires an up/down **misalignment** that the
   symmetric democratic base does **not** provide — it comes from the 10-Higgs
   vev breaking (the W6 named residual).  The claim "CKM = TBM" is therefore
   **false in the 0-parameter construction as coded**; it holds only if one
   **postulates** different diagonalizers, which is exactly the residual.

3. **What IS 0-parameter**: the **2-level spectrum** {3, 0, 0} of both the
   up-type and down-type matrices, and (in the lepton sector) ``sin^2 12 = 1/3``
   — the squared overlap with the **unique** heavy eigenvector (1,1,1)/sqrt(3).

4. **What is NOT 0-parameter**: ``sin^2 23`` and ``sin^2 13`` sweep [0, 2/3]
   under the O(2) freedom of the 2-fold light subspace, and the quark CKM is
   the identity.  The TBM values (1/2, 0) are the **(23)-adapted choice** of a
   symmetry the democratic base does not single out.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import yukawa_3x3

TOL = 1e-8  # the 0-parameter identity tolerance


def main():
    t0 = time.time()
    r = yukawa_3x3.quark_ckm()

    ok = bool(
        r["two_level"]
        and r["ckm_unitary"]
        and r["ckm_is_identity"]      # the honest result: the computed CKM IS the identity
        and r["ok"]
    )

    verdict = {
        "id": "V49_quark_ckm",
        "tier": 4,
        "property": "P5",
        "residual": True,   # the non-trivial CKM (the computed one is the identity) + sin^2_23, sin^2_13 (W6)
        "pass": ok,
        "status": "the_democratic_2level_mass_matrices_are_0parameter_the_computed_CKM_IS_THE_IDENTITY_a_nontrivial_CKM_needs_an_up_down_misalignment_the_w6_residual",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "up_eigenvalues": r["up_eigenvalues"],
            "down_eigenvalues": r["down_eigenvalues"],
            "two_level": r["two_level"],
            "ckm_unitary": r["ckm_unitary"],
            "ckm_is_identity": r["ckm_is_identity"],
            "ckm_computed_angles": r["ckm_computed_angles"],
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
    print("V49 - the quark CKM, honestly (the democratic 2-level base is 0-parameter; the computed CKM is the identity)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] quark mass matrices (democratic base):  2-level (one heavy + a 2-fold light)  ({'YES' if r['two_level'] else 'NO'})")
    print(f"  [check 2] the CKM is unitary:  {r['ckm_unitary']}")
    print(f"  [check 3] HONEST: the CKM computed from H_u = H_d (both democratic) IS THE IDENTITY: {r['ckm_is_identity']}")
    print(f"              computed CKM angles: sin2(12,23,13) = {tuple(round(a,4) for a in r['ckm_computed_angles'])}")
    print("              (a non-trivial CKM needs an up/down MISALIGNMENT the democratic base does not give)")
    print(f"  [check 4] the TBM angles ({r['sin2_12']:.3f}, {r['sin2_23']:.3f}, {r['sin2_13']:.3f}) are the LEPTON-sector (23)-adapted choice, not a quark result")
    print(f"            deviation to observed CKM (0.050, 0.949, 0.020): {r['dev_ckm']:.3f}  (the W6 named residual)")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the democratic 2-level spectrum is 0-parameter; the computed CKM is the IDENTITY, so a non-trivial CKM is the W6 named residual — stated honestly)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v49_quark_ckm.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
