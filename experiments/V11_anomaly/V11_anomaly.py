#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V11 - the 16 Weyl is anomaly-free: all 8 SM gauge anomalies cancel
(the hardest unification requirement — a chiral gauge theory is inconsistent otherwise)."""

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
    a = matter.anomalies()
    max_abs = max(abs(v) for v in a.values())
    all_zero = max_abs < 1e-12
    verdict = {
        "id": "V11_anomaly",
        "tier": 3,
        "property": "P4",
        "pass": all_zero,
        "residual": False,
        "status": "the_16_weyl_is_anomaly_free",
        "tol": TOL,
        "data": {'max_abs_anomaly': float(max_abs), **{k: float(v) for k, v in a.items()}},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v11_anomaly.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V11 - anomaly cancellation (the hard unification requirement)  [tier 3, P4]")
    print("=" * 68)
    print("  the 8 SM gauge-anomaly coefficients (all must be 0):")
    for k, v in a.items():
        print(f"    {k:12s}:  {v:+.3e}")
    print(f"  max |anomaly|:  {max_abs:.3e}")
    print(f"  VERDICT:  {'PASS (the 16 Weyl is anomaly-free)' if all_zero else 'FAIL (inconsistent)'}")
    print("  JSON -> " + os.path.join(outdir, "v11_anomaly.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if all_zero else 1

if __name__ == "__main__":
    sys.exit(main())
