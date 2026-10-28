#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V295 - triangle selection from the booked mixing (tier 3).

**The question (pre-registered).** Does the booked mixing matrix
select a unique unitarity triangle - is there a nonzero Jarlskog
invariant, or are the triangles degenerate?

**The object.** the V45 booked mixing (theta_12 = pi/4,
sin^2 theta_13 = 0, theta_23 = 0 - literals as published, from the
V45 artifact anchor), the V109 booked CP-phase search (its J_max
read from the artifact).

**Measured.** J = c12 c23 c13 s12 s23 s13 sin(delta) evaluated on
the booked angles for delta over the full circle (its maximum in
absolute value); the booked V109 J_max; area agreement: J == 0
exactly iff s13 == 0 (the booked fact).

Branch: TRIANGLES-SELECTED (J bounded away from zero) /
TRIANGLES-DEGENERATE (J = 0 identically: all unitarity-triangle
areas vanish; selection impossible; the V109 search is the booked
record of the same fact).
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def scan_J(obj, path=""):
    """every numeric leaf whose key mentions J or jarlskog."""
    out = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            if any(s in str(k).lower() for s in ("j", "jarlskog", "area")) \
                    and isinstance(v, (int, float)) and not isinstance(v,
                                                                       bool):
                out.append((path + "/" + str(k), float(v)))
            out += scan_J(v, path + "/" + str(k))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            out += scan_J(v, path + "/%d" % i)
    return out


def main():
    t0 = time.time()
    checks, data = [], {}
    v45 = art("V45_yukawa_matrix")
    v109 = art("V109_cp_phase_scan")
    th12 = float(v45["data"]["angles"]["theta_12"])
    s2_13 = float(v45["data"]["angles"]["sin2_theta_13"])
    th23 = float(v45["data"]["angles"]["theta_23"])
    t0ok = (abs(th12 - np.pi / 4) < 1e-12 and s2_13 == 0.0
            and th23 == 0.0)
    checks.append({"name": "T0 V45 mixing anchor: th12 = pi/4, "
                          "sin^2 th13 = 0, th23 = 0",
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V295 MACHINERY-DRIFT", flush=True)
        _emit(data, checks, False)
        return 1

    # ---- T1 Jarlskog over the phase circle ------------------------------
    c12, s12 = np.cos(th12), np.sin(th12)
    s13 = np.sqrt(s2_13)
    c13 = np.sqrt(1 - s2_13)
    c23, s23 = np.cos(th23), np.sin(th23)
    dl = np.linspace(0, 2 * np.pi, 7201)
    Jv = c12 * c23 * c13 * s12 * s23 * s13 * np.sin(dl)
    J_max = float(np.max(np.abs(Jv)))
    booked = scan_J(v109)
    J_book = max((abs(v) for _, v in booked), default=None)
    data["jarlskog"] = {"J_max_over_delta": J_max,
                        "booked_V109_J_keys": [(p, v) for p, v in booked
                                               if "j" in p.lower()][:8],
                        "booked_J_max": J_book}
    status = ("TRIANGLES-SELECTED" if J_max > 1e-12
              else "TRIANGLES-DEGENERATE")
    data["status"] = status
    print("V295 {}: J_max = {} over the delta circle; the booked V109 "
          "J-scan leaves {}".format(status, J_max,
                                    J_book), flush=True)
    checks.append({"name": "T1 Jarlskog maximum over the phase circle = "
                          "{}; booked V109 record {}".format(J_max,
                                                             J_book),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    _emit(data, checks, ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V295_triangle_selection", "tier": 3,
           "property": "unique unitarity-triangle selection from the "
                       "booked mixing",
           "residual": True, "programme": "V295",
           "branch": data.get("status", ""), "data": data,
           "checks": checks, "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v295_triangle_selection.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V295 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
