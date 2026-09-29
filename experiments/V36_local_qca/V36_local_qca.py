#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V36 — the local QCA: the causal radius **implemented** (the 't Hooft
program, C1).

**The 0-parameter result.**  The QCA step of CL10 is a **product of local
(2-cell) unitaries**:

- the **internal layer** (the 0-parameter 45 Spin(10) plane rotations at
  ``g_10``) on each cell;
- the **bond layer** (the x, y, z directions **cycle** with the time step —
  one direction per step, so a signal gets **one cell of budget per step**;
  the 't Hooft staggered even / odd 2-cell bonds): the **0-parameter
  swap-rotation** ``U_bond = cos(g_10) I - i sin(g_10) P`` (``P`` = the
  0-parameter swap of the two cell states, a permutation — no angle).

The **only** number of the step is ``g_10`` (the 0-parameter GUT coupling).
Therefore:

- the step is **unitary** (the norm is preserved);
- a **localized** excitation stays **inside** the **light cone** (the
  Manhattan ball of radius k after k steps) — **causality is a property of
  the dynamics, not an axiom**;
- the **front** of the support **reaches** the light cone **exactly**
  (the **front speed is c = 1** — the light cone is **realized**);
- the light cone is a **structural** property (the **same** cone is
  realized for different 0-parameter angles — it is carried by the
  **2-cell operator structure**, not by the value of the angle).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import local_qca

TOL = 1e-9


def main():
    t0 = time.time()
    r = local_qca.verify(tol=TOL)

    ok = (r["ok"]
          and r["unitary_err"] < TOL
          and r["causal_ok"]
          and abs(r["front_speed"] - 1.0) < 1e-9
          and r["structural_ok"]
          and r["nontrivial_ok"])

    verdict = {
        "id": "V36_local_qca",
        "tier": 2,
        "property": "P2,P4",
        "residual": False,          # the local QCA is a 0-parameter PASS (no new residual)
        "pass": bool(ok),
        "status": "the_causal_radius_1_is_implemented_the_light_cone_is_realized",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "unitary_err": r["unitary_err"],
            "causal_ok": r["causal_ok"],
            "front_speed": r["front_speed"],
            "structural_ok": r["structural_ok"],
            "nontrivial_ok": r["nontrivial_ok"],
            "support_growth": r["support_growth"],
            "internal_only_growth": r["internal_only_growth"],
            "residual": "none new (the QCA lattice RG remains the named open problem, as before)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v36_local_qca.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V36 - the local QCA: the causal radius IMPLEMENTED (C1)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter angle:  g_10 = {r['g10']}  (the 45 rotations + the swap-rotations)")
    print(f"  the unitary (the norm is preserved):  residual = {r['unitary_err']:.2e}")
    print(f"  the LIGHT CONE is REALIZED:  the front speed c = {r['front_speed']:.3f}  (1 cell / step)")
    print(f"    the support after k steps (the ball of radius k):  {r['support_growth']}")
    print(f"    the internal layer alone (no 2-cell coupling):  {r['internal_only_growth']}")
    print(f"  the light cone is STRUCTURAL (independent of g_10):  {r['structural_ok']}")
    print(f"  VERDICT:  {'PASS (the causal radius 1 is a property of the local dynamics, not an axiom)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v36_local_qca.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
