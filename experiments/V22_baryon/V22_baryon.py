#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V22 - the baryon asymmetry (the Sakharov conditions).

A **baryon asymmetry** (the observed matter excess over antimatter) requires,
by the **Sakharov conditions** (1967), **three** ingredients:

1. **baryon-number violation** (the B-number is not conserved);
2. **C and CP violation** (the laws are not symmetric under C or CP);
3. **departure from thermal equilibrium** (the dynamics is not reversible).

The **0-parameter QCA** (a product of local unitaries) satisfies **none** of
the three:
- **unitary** → **all** global charges are conserved (baryon number and B-L
  are **exact** symmetries — no baryon-number violation);
- **CPT-invariant** (the 16/16-bar CPT pair, V10) and **CP-conserving** at
  0-parameter (no C/CP violation);
- **reversible** (unitary) → it stays in **equilibrium** (no departure from
  equilibrium).

Therefore the 0-parameter QCA **does not generate** the observed baryon
asymmetry — it is a **named** residual (the theory **names** it; it does not
**derive** it).  A baryon asymmetry would require a **beyond-0-parameter**
(CPT/CP-violating, non-unitary) process.  This is an **honest, tested**
negative result: the theory *knows* it does not explain the matter excess, and
says so.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import walls

TOL = 1e-6


def main():
    t0 = time.time()
    r = walls.baryon_asymmetry()
    # the QCA satisfies none of the 3 Sakharov conditions
    none_satisfied = (not r["baryon_number_violated"]
                      and not r["c_or_cp_violated"]
                      and not r["out_of_equilibrium"])
    not_generated = (not r["generated"])
    ok = none_satisfied and not_generated

    verdict = {
        "id": "V22_baryon",
        "tier": 5,
        "property": "P6",
        "residual": True,          # the baryon asymmetry is the named residual
        "pass": bool(ok),
        "status": "the_qca_does_not_generate_the_baryon_asymmetry",
        "tol": TOL,
        "data": {
            "baryon_number_violated": r["baryon_number_violated"],
            "c_or_cp_violated": r["c_or_cp_violated"],
            "out_of_equilibrium": r["out_of_equilibrium"],
            "all_sakharov_satisfied": r["all_sakharov"],
            "generated": r["generated"],
            "named_residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v22_baryon.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V22 - the baryon asymmetry (the Sakharov conditions)")
    print("=" * 68)
    print("  the 3 Sakharov conditions (the 0-parameter QCA):")
    print(f"    1. baryon-number violation:  {r['baryon_number_violated']}  (unitary -> B conserved)")
    print(f"    2. C/CP violation:           {r['c_or_cp_violated']}  (CPT-invariant, CP-conserving at 0-param)")
    print(f"    3. departure from equilibrium:  {r['out_of_equilibrium']}  (unitary -> reversible -> equilibrium)")
    print("  -> the 0-parameter QCA satisfies NONE of the three")
    print(f"  the baryon asymmetry is NOT generated -> NAMED RESIDUAL ({r['residual']})")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (an honest, TESTED negative result: the QCA does not generate the matter excess)")
    print("  JSON -> " + os.path.join(outdir, "v22_baryon.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
