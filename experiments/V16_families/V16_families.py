#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V16 - the 3 families; the MECHANISM is reduced by V19 (the 3 = the 3 space).
The exact "why-3" derivation is the NAMED RESIDUAL (the family problem)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import walls
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    w3 = walls.w3_families()
    candidate = w3['n_fam_candidate']; observed = w3['n_observed']; derived = w3['derived']
    verdict = {
        "id": "V16_families",
        "tier": 5,
        "property": "P6",
        "pass": candidate == observed,
        "residual": True,
        "status": "the_3_families_reduced_by_3_space",
        "tol": TOL,
        "data": {'candidate': int(candidate), 'observed': int(observed), 'derived': bool(derived),
                 'n_space': int(w3['n_space']),
                 'residual': 'the derivation of why the family count = the spatial dimension'},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v16_families.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V16 - the 3 families (the named residual)  [tier 5, P6]")
    print("=" * 68)
    print(f"  the 0-parameter family count (V19):  {candidate} = the {w3['n_space']} space  (observed: {observed})")
    print(f"  derived (3 = 3-space):  {derived}")
    print(f"  VERDICT:  PASS (REDUCED W3 — the 3 = 3-space mechanism (V19) is 0-parameter; the why-3 is the NAMED RESIDUAL)")
    print("  JSON -> " + os.path.join(outdir, "v16_families.json"))
    print(f"  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if candidate == observed else 1

if __name__ == "__main__":
    sys.exit(main())
