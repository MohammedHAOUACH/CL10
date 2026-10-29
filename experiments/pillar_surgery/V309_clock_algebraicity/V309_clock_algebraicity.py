#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V309 - the algebraicity of the data clock (tier 3, O8/O4,
pillar-surgery).

**The question (pre-registered).** is g_obs = 0.19999798835751165
(V206's booked data-clock read, artifact) an algebraic number of small
height, i.e. could a root rule generate it?

**Measured (pre-registered).**
1. `mpmath.findpoly(g_obs, d)` for every degree d <=
   pillar_surgery.pslq_max_degree at 60 dp, height cap
   pillar_surgery.pslq_max_coeff (the mpmath pattern of V276 T2); any
   candidate polynomial is re-verified by substitution at 50 dp before
   it is booked.
2. the rational census: every a/b with 1 <= b <= 5000 inside
   pillar_surgery.clock_band_abs of g_obs (count, and the closest few
   listed whether or not they are inside).
3. the booked clock gap: g_cross_rung10 - g_obs from the V206 artifact,
   reported to full precision as the named gap of order 1.35e-4 that
   V310's expression census must supply.

**Acceptance (pre-registered).** CLOCK-ALGEBRAIC iff an integer
polynomial of degree <= pslq_max_degree and height <= pslq_max_coeff
vanishes at g_obs well below the float64 noise floor of the booked
decimal (substituted at 50 dp); the polynomial is then booked and the
V314 ledger must use it.  CLOCK-NO-RELATION-DECLARED otherwise: the
data clock carries no small algebraic rule on the reserved class, and
the gap and the rational census are booked as the bounded no-go.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402

from cl10 import config as _config  # noqa: E402

CFG = _config.CFG
DEN_MAX = 5000                 # the booked rational census range


def art(name, key):
    """A sibling verdict JSON, read by experiment name: a published
    number has a source, it is not a decimal restated here (cl10.paths
    resolves the category)."""
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def poly_string(coeffs):
    """An integer-coefficient list [c_0..c_d] as a readable polynomial."""
    terms = []
    for d, c in enumerate(coeffs):
        if c == 0:
            continue
        x = "" if d == 0 else ("x" if d == 1 else "x^%d" % d)
        terms.append(("%+d%s" % (c, x)) if terms else ("%d%s" % (c, x)))
    return " ".join(terms)


def main():
    import mpmath

    dmax = int(CFG.pillar_surgery.pslq_max_degree)
    hcap = int(CFG.pillar_surgery.pslq_max_coeff)
    band = float(CFG.pillar_surgery.clock_band_abs)
    checks, data = [], {}

    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    dm15 = float(art("V206_clock_tension",
                     "data.data_clock.g_obs_minus_fifth"))

    # ---- T0 gate: the booked clock arithmetic is consistent ------------
    fifth = 1.0 / 5.0
    t0_def = abs((fifth + dm15) - g_obs)
    data["t0"] = {"g_obs": g_obs, "g_cross_rung10": g_cross,
                  "g_obs_minus_fifth_booked": dm15,
                  "fifth_plus_booked_diff_minus_g_obs": t0_def}
    checks.append({"name": "T0 booked clock arithmetic: 1/5 + booked "
                          "g_obs_minus_fifth reproduces g_obs (|d| "
                          "{:.1e} < 1e-15)".format(t0_def),
                   "pass": bool(t0_def < 1e-15)})
    if t0_def >= 1e-15:
        data["status"] = "MACHINERY-DRIFT"
        print("V309 MACHINERY-DRIFT: the booked V206 clock arithmetic is "
              "not internally consistent; nothing is interpreted",
              flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- 1. findpoly over degrees 1..dmax -------------------------------
    mpmath.mp.dps = 60
    x = mpmath.mpf(repr(g_obs))
    noise = mpmath.mpf("1e-14")     # substitution floor: a booked float
                                    # is exact only to ~16 digits, so a
                                    # relation must vanish well below 5e-17
                                    # x (height x derivative) ~ 1e-14
    found = None
    found_verified = None
    per_degree = {}
    for d in range(1, dmax + 1):
        rel = mpmath.findpoly(x, d, maxcoeff=hcap)
        if rel is None:
            per_degree[d] = None
            continue
        res = abs(sum(mpmath.mpf(c) * x ** k
                      for k, c in enumerate(rel)))
        height = max(abs(int(c)) for c in rel)
        per_degree[d] = {"coeffs": [int(c) for c in rel],
                         "poly": poly_string([int(c) for c in rel]),
                         "residual_60dp": float(res), "height": height}
        if res < noise:
            # re-verify at 50 dp independently (pre-registration duty)
            old = mpmath.mp.dps
            mpmath.mp.dps = 50
            x50 = mpmath.mpf(repr(g_obs))
            res50 = abs(sum(mpmath.mpf(c) * x50 ** k
                            for k, c in enumerate(rel)))
            mpmath.mp.dps = old
            if res50 < noise:
                found = [int(c) for c in rel]
                found_verified = {"degree": d, "height": height,
                                  "residual_50dp": float(res50)}
        # keep searching all degrees even after a raw hit: the full
        # degree table is the booked content
    data["findpoly"] = {"degrees": "1..%d" % dmax, "height_cap": hcap,
                        "substitution_noise_floor": float(noise),
                        "per_degree": {str(d): per_degree[d]
                                       for d in per_degree}}
    mpmath.mp.dps = 15
    checks.append({"name": "findpoly degrees 1..%d, height cap %d: %s"
                          % (dmax, hcap,
                             "no relation" if found is None else
                             "candidate {} (verified {})".format(
                                 data["findpoly"]["per_degree"],
                                 found_verified)),
                   "pass": True})

    # ---- 2. the rational census ----------------------------------------
    inside = []
    all_r = []
    for b in range(1, DEN_MAX + 1):
        a = round(g_obs * b)
        if a < 1:
            continue
        d = abs(a / b - g_obs)
        all_r.append((d, a, b))
        if d < band:
            inside.append((d, a, b))
    all_r.sort()
    closest = [{"a": a, "b": b, "residual_abs": d,
                "inside_band": bool(d < band)} for d, a, b in all_r[:5]]
    data["rational_census"] = {"den_max": DEN_MAX, "band_abs": band,
                               "count_inside": len(inside),
                               "inside": [{"a": a, "b": b,
                                           "residual_abs": d}
                                          for d, a, b in inside[:10]],
                               "closest_five": closest}
    checks.append({"name": "rational census b <= {:d}: {:d} rationals "
                           "inside {:.0e} of g_obs; closest is {:d}/{:d} "
                           "at {:.3e} (the fifth clock 1/5 sits at "
                           "{:.3e} outside the band)".format(
                               DEN_MAX, len(inside), band,
                               closest[0]["a"], closest[0]["b"],
                               closest[0]["residual_abs"], abs(dm15)),
                   "pass": True})

    # ---- 3. the named clock gap -----------------------------------------
    gap = g_cross - g_obs
    data["clock_gap"] = {"g_cross_rung10": g_cross, "g_obs": g_obs,
                         "gap_abs": gap,
                         "named": "the rung-10 crossing minus the "
                                  "data-clock read = {:.6e}".format(gap)}
    checks.append({"name": "booked clock gap g_cross_rung10 - g_obs = "
                           "{:.17g} (named 1.35e-4 for the V310 "
                           "census)".format(gap),
                   "pass": True})

    # ---- branch ----------------------------------------------------------
    branch = ("CLOCK-ALGEBRAIC" if found else "CLOCK-NO-RELATION-DECLARED")
    data["status"] = branch
    data["relation"] = found
    data["relation_verified"] = found_verified
    print("V309 {}: g_obs = {:.17g}; findpoly(1..{}, cap {}) -> {}; "
          "rationals b <= {} inside {:.0e}: {}; closest {}; clock gap "
          "{:.17g}.".format(branch, g_obs, dmax, hcap,
                            found if found else "none at any degree",
                            DEN_MAX, band, len(inside),
                            "{}/{} at {:.3e}".format(
                                closest[0]["a"], closest[0]["b"],
                                closest[0]["residual_abs"]),
                            gap),
          flush=True)
    _emit(data, checks, ok=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V309_clock_algebraicity", "tier": 3,
           "property": "is the data clock g_obs an algebraic number of "
                       "small height - could a root rule generate it "
                       "(tier 3, O8/O4)",
           "residual": True,
           "programme": "pillar-surgery V309",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v309_clock_algebraicity.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V309 done] {}".format(out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
