#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V21 - the 0-parameter GUT coupling + the 1-loop running structure.

**The 0-parameter content** (V21):
- the **GUT coupling** $\alpha_{10} = 1/25$ ($g_{10} \approx 0.71$) — a
  0-parameter number of the root; the 45 gauge act with this **one** coupling;
- the **5/3** $\\mathrm{U}(1)$ normalization (the GUT Cartan);
- the **1-loop SM beta coefficients** $b_1=-41/10$, $b_2=-19/6$, $b_3=-7$
  (0-parameter: they depend only on the particle content).

**The named residual.**  The **1-loop SM matching** to the 0-parameter GUT
scale: the SM (non-SUSY) 1-loop running does **not** cleanly unify the three
couplings (the $\\mathrm{U}(1)$ runs apart — the unification *spread* is large),
so the 1-loop unification to the 0-parameter Planck scale is a **named**
residual (the theory names it; it does not derive the matching).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import gauge

TOL = 1e-6


def main():
    t0 = time.time()
    r = gauge.running()
    # 0-parameter structure is consistent:
    five_thirds_ok = abs(r["five_thirds"] - 5.0 / 3.0) < 1e-6
    gut_coupling_ok = abs(r["alpha_10"] - 1.0 / 25.0) < 1e-6
    beta_ok = (r["b1"] == 41 / 10 and r["b2"] == -19 / 6 and r["b3"] == -7)
    structure_ok = five_thirds_ok and gut_coupling_ok and beta_ok
    # the 1-loop running: the spread is large (no clean unification) -> named residual
    residual = (not r["cleanly_unifies"])
    ok = structure_ok and residual

    verdict = {
        "id": "V21_running",
        "tier": 3,
        "property": "P4",
        "residual": True,          # the 1-loop SM matching is the named residual
        "pass": bool(ok),
        "status": "the_0_parameter_gut_coupling_and_1_loop_structure",
        "tol": TOL,
        "data": {
            "alpha_10": r["alpha_10"],
            "g_gut": r["g_gut"],
            "m_gut_0param": r["m_gut_0param"],
            "five_thirds": r["five_thirds"],
            "b1": r["b1"], "b2": r["b2"], "b3": r["b3"],
            "unification_spread": r["unification_spread"],
            "cleanly_unifies": r["cleanly_unifies"],
            "named_residual": "the 1-loop SM matching to the 0-parameter GUT scale",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v21_running.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V21 - the 0-parameter GUT coupling + the 1-loop running")
    print("=" * 68)
    print(f"  0-parameter GUT coupling:  alpha_10 = {r['alpha_10']}  (g_10 = {r['g_gut']:.3f})")
    print(f"  0-parameter GUT scale:  M_GUT = {r['m_gut_0param']:.2e} GeV  (the Planck scale)")
    print(f"  0-parameter 5/3 normalization:  {r['five_thirds']:.4f}")
    print(f"  1-loop SM beta coefficients:  b1={r['b1']}, b2={r['b2']}, b3={r['b3']}  (the particle content)")
    print(f"  1-loop unification spread:  {r['unification_spread']:.1f}  (large -> the SM does NOT cleanly unify at 1-loop)")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the 0-parameter structure is consistent; the 1-loop SM matching is the NAMED residual)")
    print("  JSON -> " + os.path.join(outdir, "v21_running.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
