#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V20 - the 0-parameter fermion mass generation (the 10-Higgs).

**Claim.**  In the Cl(10) GUT the **10** (the grade-1 generators, the GUT vector)
is **off-diagonal** in the 16/16-bar (matter/antimatter) split: it maps the 16
to the 16-bar and the 16-bar to the 16 (it **anticommutes** with the
pseudoscalar).  Therefore a **10-Higgs vev** (a specific grade-1 element)
couples the 16 to its 16-bar conjugate and gives the **0-parameter Dirac mass
term**

$$16\\;10\\;16\\bar \\;\\to\\; 1 .$$

The **0-parameter** content (V20):
- the **diagonal** blocks (16->16, 16-bar->16-bar) are **0** (the 10 is an odd
  element);
- the **off-diagonal** blocks (16->16-bar, 16-bar->16) are **nonzero** (the mass
  matrices);
- the 16 also contains the **sterile N** (the 1 of the 16, $B\\!-\\!L=0$), the
  0-parameter source of the **seesaw** (the 2-step W2).

This is the 0-parameter **structure** of the fermion mass term.  The **exact**
mass values (the 3-generation hierarchy, the mixing angles) are **named**
residuals (the theory gives the *structure*, not the *values*).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import matter

TOL = 1e-6


def main():
    t0 = time.time()
    r = matter.mass_generation()
    offdiag = r["off_diagonal"]
    diag_ok = (r["diag_16"] < TOL and r["diag_16bar"] < TOL)
    offdiag_nonzero = (r["A_16_16bar"] > 1.0 and r["B_16bar_16"] > 1.0)
    ok = offdiag and diag_ok and offdiag_nonzero and r["sterile_N_in_16"]

    verdict = {
        "id": "V20_massgen",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the exact mass values are named (the structure is 0-parameter)
        "pass": bool(ok),
        "status": "the_10_higgs_gives_the_0_parameter_mass_term",
        "tol": TOL,
        "data": {
            "n16": r["n16"], "n16bar": r["n16bar"],
            "diag_16": r["diag_16"], "diag_16bar": r["diag_16bar"],
            "A_16_16bar": r["A_16_16bar"], "B_16bar_16": r["B_16bar_16"],
            "mass_term": r["mass_term"],
            "sterile_N_in_16": r["sterile_N_in_16"],
            "named_residuals": "the exact mass values (the 3-gen hierarchy, the mixing) are named",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v20_massgen.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V20 - the 0-parameter fermion mass generation (the 10-Higgs)")
    print("=" * 68)
    print(f"  16 (matter) / 16-bar (antimatter):  {r['n16']} / {r['n16bar']}")
    print(f"  the 10 (grade-1) is OFF-DIAGONAL:  diag(16->16)={r['diag_16']:.1e}, diag(16bar->16bar)={r['diag_16bar']:.1e}  (0)")
    print(f"  off-diagonal mass matrices:  A(16->16bar)={r['A_16_16bar']:.3f},  B(16bar->16)={r['B_16bar_16']:.3f}  (nonzero)")
    print("  -> the 10-Higgs vev gives the 0-parameter Dirac mass term  16 . 10 . 16bar -> 1")
    print("  the sterile N (the 1 of the 16) is in the 16 -> the 0-parameter seesaw (W2)")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the 0-parameter mass *structure*; the exact values are NAMED residuals)")
    print("  JSON -> " + os.path.join(outdir, "v20_massgen.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
