#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V47 — the **2-loop** running: the **UV-finite** QCA vs the **UV-divergent**
perturbative QFT (the W4 completion, the 2-loop UV completion).

**The 0-parameter contrast.**

1. **The QCA 2-step (finite).**  The QCA step ``U`` (the 32x32) is **unitary**;
   the **2-step** ``U^2`` is **unitary** (``U^2 (U^2)^* = 1``) — the
   **2-loop** (the 2-step) is **finite** (the **UV-finite** property: the QCA
   is a **unitary** automaton, the n-step is **bounded** for all n).

2. **The QCA 2-particle (bounded).**  The 2-particle S-matrix (the **2-loop**
   analog, the 1024x1024) has **bounded** amplitudes (``max_offdiag <= 1``) —
   the **UV-finite** property.

3. **The perturbative 2-loop (divergent).**  The **2-loop** beta coefficient
   (``b2``) makes the **U(1)** Landau pole **closer**: the **1-loop** pole is
   **above** M_P (the V12 UV-bounded), but the **2-loop** pole is **below**
   M_P (the **UV divergence** is **reached** **before** the cutoff).

**The 0-parameter result.**  The **QCA** is **UV-finite** (the 2-step is
**unitary**, the 2-particle is **bounded**), **whereas** the **perturbative**
QFT is **UV-divergent** (the 2-loop Landau pole).  This is the **UV-finite**
property of the QCA (the **W4** completion).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import unification

TOL = 1e-8  # the 0-parameter identity tolerance


def main():
    t0 = time.time()
    r = unification.running_2loop()

    ok = bool(
        r["qca_2step_finite"]
        and r["qca_2particle_bounded"]
        and r["pert_2loop_divergent"]
        and r["uv_finite"]
        and r["ok"]
    )

    verdict = {
        "id": "V47_2loop_uvfinite",
        "tier": 1,
        "property": "P2",
        "residual": True,   # the 2-loop UV completion: the QCA is UV-finite, the perturbative is UV-divergent
        "pass": ok,
        "status": "the_qca_is_uvfinite_the_2step_is_unitary_and_the_2particle_is_bounded_the_perturbative_2loop_is_uvdivergent",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "qca_2step_unitary_residual": r["qca_2step_unitary_residual"],
            "qca_2particle_max_offdiag": r["qca_2particle_max_offdiag"],
            "pert_1loop_landau_pole": r["pert_1loop_landau_pole"],
            "pert_2loop_landau_pole": r["pert_2loop_landau_pole"],
            "pert_2loop_divergent": r["pert_2loop_divergent"],
            "uv_finite": r["uv_finite"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v47_2loop_uvfinite.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V47 - the 2-loop running: the UV-finite QCA vs the UV-divergent perturbative QFT")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] QCA 2-step (U^2):  unitary residual = {r['qca_2step_unitary_residual']:.2e} (finite: the 2-loop is bounded)")
    print(f"  [check 2] QCA 2-particle:  max_offdiag = {r['qca_2particle_max_offdiag']:.3f} (bounded, <= 1: UV-finite)")
    print(f"  [check 3] perturbative 2-loop:  1-loop pole = {r['pert_1loop_landau_pole']:.2e} (above M_P), 2-loop pole = {r['pert_2loop_landau_pole']:.2e} (below M_P: UV-divergent)")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the QCA is UV-finite: the 2-step is unitary, the 2-particle is bounded; the perturbative 2-loop is UV-divergent: the Landau pole is below M_P)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v47_2loop_uvfinite.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
