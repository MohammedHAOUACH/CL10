#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V25 - the holographic area law (S scales with AREA, not volume).

The **holographic principle** — the hallmark of a genuine 0-parameter unification.  The **0-parameter**
content of the CL10 QCA:

- the **area law**: the QCA vacuum entropy of a region of side $L$ scales as
  $S \\propto L^2$ (the **surface area** $6L^2$), **not** $L^3$ (the volume) —
  the scaling exponent is **exactly 2**;
- the **0-parameter** entropy per unit area is the **Bekenstein-Hawking**
  $S = A/4$ ($1/4$ per plaquette, $G = 1$ — the $(1,3)$ metric);
- the **holographic ratio** $S/V \\propto 1/L$ (the exponent is **exactly
  $-1$**): as the region grows, the entropy per unit volume falls as $1/L$ —
  the **bulk is encoded on the boundary**.

This is a **0-parameter** (P2/P4) result: the QCA's information content is
**holographic** (area-law, not volume-law).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import holo

TOL = 1e-12


def main():
    t0 = time.time()
    r = holo.area_law()
    ok = (abs(r["area_law_exponent"] - 2.0) < 1e-6
          and abs(r["holographic_ratio_exponent"] + 1.0) < 1e-6
          and abs(r["per_plaquette"] - 0.25) < 1e-12)

    verdict = {
        "id": "V25_holography",
        "tier": 4,
        "property": "P2,P4",
        "residual": False,         # the area law is a 0-parameter PASS (not a residual)
        "pass": bool(ok),
        "status": "the_qca_is_holographic_S_scales_with_area_not_volume",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "L_values": r["L_values"],
            "per_plaquette": r["per_plaquette"],
            "S": r["S"],
            "R": r["R"],
            "area_law_exponent": r["area_law_exponent"],
            "holographic_ratio_exponent": r["holographic_ratio_exponent"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v25_holography.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V25 - the holographic area law (S ~ AREA, not volume)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the 0-parameter entropy per plaquette:  {r['per_plaquette']}  (Bekenstein-Hawking A/4, G=1)")
    print(f"  S(L) = {['%.1f' % s for s in r['S']]},  S/V = {['%.4f' % x for x in r['R']]}")
    print(f"  area-law exponent (S vs L)   = {r['area_law_exponent']:.4f}  (= 2 -> AREA)")
    print(f"  holographic exponent (S/V)   = {r['holographic_ratio_exponent']:.4f}  (= -1 -> 1/L)")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the QCA is holographic: S ~ AREA, S/V ~ 1/L)")
    print("  JSON -> " + os.path.join(outdir, "v25_holography.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
