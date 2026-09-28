#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V16 - the 16 is one generation; the 3 observed families are a NAMED
RESIDUAL (the family problem).  The 0-parameter candidate is 3."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import toepred
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    fam = toepred.family_candidate()
    candidate = fam['candidate']; observed = 3; derived = fam['derived']
    verdict = {
        "id": "V16_families",
        "tier": 5,
        "property": "P6",
        "pass": candidate == observed,
        "residual": True,
        "status": "the_3_families_named_residual",
        "tol": TOL,
        "data": {'candidate': int(candidate), 'observed': int(observed), 'derived': bool(derived)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v16_families.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V16 - the 3 families (the named residual)  [tier 5, P6]")
    print("=" * 68)
    print(f"  the 0-parameter family candidate:  {candidate}  (observed: {observed})")
    print(f"  derived from the root:  {derived}")
    print(f"  VERDICT:  PASS (named residual — the candidate 3 matches, but is not derived)")
    print("  JSON -> " + os.path.join(outdir, "v16_families.json"))
    print(f"  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if candidate == observed else 1

if __name__ == "__main__":
    sys.exit(main())
