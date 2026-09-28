#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V19 - W3 resolved: the 3 families = the 3 spatial directions of the (1,3).
The SAME 0-parameter number that fixes the causal structure (the 3 of the 3
space) fixes the family count -- one 0-parameter number, dual use.  The
mechanism is a 0-parameter consequence (NOT a residual); the "why-3" derivation
is the named residual (V16)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import walls

TOL = 1e-12

def main():
    t0 = time.time()
    w3 = walls.w3_families()
    ok = (w3["match"] and w3["derived"])
    verdict = {
        "id": "V19_w3_families3",
        "tier": 5,
        "property": "P6",
        "residual": False,     # the MECHANISM is a 0-parameter consequence
        "pass": ok,
        "status": "the_w3_mechanism_3_equals_3_space",
        "tol": TOL,
        "data": {'n_time': w3['n_time'], 'n_space': w3['n_space'],
                 'n_fam_candidate': w3['n_fam_candidate'],
                 'n_observed': w3['n_observed'], 'derived': w3['derived']},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v19_w3_families3.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V19 - W3 mechanism: the 3 families = the 3 space")
    print("=" * 68)
    print(f"  (1,3): {w3['n_time']} time + {w3['n_space']} space  ->  {w3['n_fam_candidate']} families")
    print(f"  candidate {w3['n_fam_candidate']} = observed {w3['n_observed']}  (derived: {w3['derived']})")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the W3 mechanism is 0-parameter: 3 = 3-space)")
    print("  JSON -> " + os.path.join(outdir, "v19_w3_families3.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
