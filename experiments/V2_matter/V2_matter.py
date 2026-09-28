#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V2 - the matter is the 16 Weyl (one SM generation, 10 + 5bar + 1 + nu^c);
the 16-bar is the antimatter (opposite B-L)."""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import matter

TOL = 1e-12

def main():
    t0 = time.time()
    p, m = matter.weyl16()
    n_weyl=p.shape[1]; n_16bar=m.shape[1]; total_dim=matter.dimension_check()
    bl16=matter.bl_charge('16'); bl16b=matter.bl_charge('16bar'); n_sm=len(matter.SM_FIELDS)
    verdict = {
        "id": "V2_matter",
        "pass": (n_weyl==16 and n_16bar==16 and total_dim==16),
        "status": "the_16_weyl_is_one_generation",
        "tol": TOL,
        "data": {'n_weyl':int(n_weyl),'n_16bar':int(n_16bar),'total_dim':int(total_dim),'n_sm':int(n_sm),'bl16':int(bl16),'bl16bar':int(bl16b)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v2_matter.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V2 - the matter: the 16 Weyl (one generation)")
    print("=" * 68)
    print(f"  16 Weyl / 16-bar:  {n_weyl} / {n_16bar}")
    print(f"  total SM-field dim:  {total_dim}  (10 + 5bar + 1 + nu^c)")
    print(f"  B-L charge:  matter -> {bl16},  antimatter -> {bl16b}")
    print(f"  VERDICT:  {'PASS' if n_weyl==16 and total_dim==16 else 'FAIL'}  (one generation)")
    print("  JSON -> " + os.path.join(outdir, "v2_matter.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if (n_weyl==16 and n_16bar==16 and total_dim==16) else 1

if __name__ == "__main__":
    sys.exit(main())
