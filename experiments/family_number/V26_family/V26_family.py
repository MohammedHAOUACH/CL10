#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V26 - the 3-family construct (the 3 generations from the (1,3)).

The **0-parameter 3-family construct**.  The 16 Weyl is **one** SM generation.
The **3** observed families are a **0-parameter consequence** of the single
object, in **three** equivalent ways:

- (1) the **3 space directions** of the $(1,3)$ causal structure;
- (2) the **3 two-planes** $\\binom{3}{2}$ of the 3-space;
- (3) the **rank** (3) of the spatial metric.

All three 0-parameter numbers are **3** (one number, three readings).  The
**3 generations** = $16 \\times 3 = 48$ Weyl (one 16 per family).

The **named residual**: the *why-3* **derivation** — the theory fixes the
0-parameter candidate (3, three ways) but does not derive *why* the family
count is exactly the spatial dimension.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, ROOT)

from cl10 import family

TOL = 1e-12


def main():
    t0 = time.time()
    r = family.family_construct()
    ok = r["match"] and r["all_three_ways_3"] and r["n_families_observed"] == 3

    verdict = {
        "id": "V26_family",
        "tier": 1,
        "property": "P2",
        "residual": True,          # the 'why-3' is a named residual
        "pass": bool(ok),
        "status": "the_3_generations_are_a_0_parameter_3_three_ways",
        "tol": TOL,
        "data": {
            "n_space": r["n_space"],
            "n_planes_C32": r["n_planes_C32"],
            "n_rank": r["n_rank"],
            "n_family_candidate": r["n_family_candidate"],
            "n_families_observed": r["n_families_observed"],
            "n_weyl_1gen": r["n_weyl_1gen"],
            "n_weyl_3gen": r["n_weyl_3gen"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v26_family.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V26 - the 3-family construct (the 3 generations from the (1,3))")
    print("=" * 68)
    print(f"  the 0-parameter 3, three ways:  space={r['n_space']}, planes=C(3,2)={r['n_planes_C32']}, rank={r['n_rank']}")
    print(f"  the 16 Weyl (one generation) -> 3 generations:  {r['n_weyl_1gen']} x {r['n_family_candidate']} = {r['n_weyl_3gen']} Weyl")
    print(f"  named residual:  {r['residual']}")
    print(f"  VERDICT:  {'PASS (0-parameter 3, three ways)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v26_family.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
