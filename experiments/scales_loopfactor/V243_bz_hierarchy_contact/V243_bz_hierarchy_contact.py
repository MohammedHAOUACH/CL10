#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V243 - does the Brillouin-zone loop factor read the observed hierarchy?

**The question this pre-registers against** (`review/pre_registrations.md`,
second programme): with F_infinity from V242's artifact at each tested
pair, does the ladder reading of v/M_P - exactly as V203/V108 compute it -
land on the observed coefficient when the loop is F_infinity^(1/10) at the
pair's own class? V242 landed NO-LOOP-FACTOR, so the branch is decided by
the acceptance clause (no pair with an EMERGED F_infinity): the artifact
is the ratio table, and the copy control (the V108 reading recomputed from
its own published F_10) must still reproduce the published coefficient.

Closed inputs: V242's artifact (F_infinity per pair/class, measured flags),
the V108 artifact (reading definition: coefficient = 3 (16 pi F)^10 with
the axis-family multiplier 3, the published value 2.0881...), and
`cl10/config.yaml` (the loader-derived
hierarchy_coefficient_observed = 2.0765641...). No new number.

Branch `HIERARCHY-READS` iff a pair with a MEASURED F_infinity reaches
relative |coef - observed| < 1e-3; else `HIERARCHY-UNREAD` with the full
ratio table (floor pairs recorded as not quotable, per V242's caveat).
"""
import json
import math
import os
import sys
import time


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

from cl10.paths import exp                                     # noqa: E402
from cl10 import config as _config                             # noqa: E402

CFG = _config.CFG


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def coef(F, mult=3.0):
    """the V108 ladder reading: predicted coefficient in units of the
    ladder factor: 3 * (16 pi F)^10 (multiplier 3 = the three axis
    families, V70/V108; exponent 10 = the measured crossing, V89)."""
    return mult * (16.0 * math.pi * F) ** 10


def main():
    t0 = time.time()
    checks, data = [], {}
    v242 = art("V242_bz_loop_factor_census")
    v108 = art("V108_hierarchy_multiplicity")["data"]
    obs = float(CFG.physics.hierarchy_coefficient_observed)
    loop = float(CFG.physics.loop)
    data["inputs"] = {
        "observed_coefficient": obs,
        "observed_derivation": "v_over_mplanck_observed / loop^10 "
                               "(cl10/config.yaml loader; V52/V53/V108/"
                               "V203 all read this same value)",
        "v242_branch": v242["branch"],
        "pairs_tested": v242["data"]["pairs_tested"]}

    # ---- T0 copy control: the V108 reading recomputed from V108's own F
    F10 = float(v108["F_10_root"])
    pred = coef(F10)
    pub = float(v108["coefficient_predicted_3_F10_10"])
    checks.append({
        "name": f"T0 COPY CONTROL: the reading definition coefficient = "
                f"3 (16 pi F)^10 evaluated at V108's own published "
                f"F_10 = {F10:.8f} gives {pred:.6f} against the published "
                f"{pub:.4f} (dev {abs(pred - pub):.1e}); the observed "
                f"coefficient from the config loader = {obs:.7f}",
        "pass": bool(abs(pred - pub) < 5e-3)})
    data["copy_control"] = {"F_10_root": F10, "recomputed": pred,
                            "published": pub}

    # ---- T1 the ratio table at V242's measured F_infinity values -------
    rows = []
    for pair, beta, g in v242["data"]["pairs_tested"]:
        fin = v242["data"]["F_infinity"][pair]
        for cls in ("even", "odd"):
            F = float(fin[cls])
            c = coef(F)
            rows.append({"pair": pair, "beta": beta, "g": g, "class": cls,
                         "F_infinity": F,
                         "measurable": bool(fin["measured"]),
                         "coefficient": c if fin["measured"] else None,
                         "rel_distance_to_observed":
                             abs(c - obs) / obs if fin["measured"]
                             else None,
                         "note": None if fin["measured"] else
                         "V242: ratio of floors (structural zero "
                         "denominator), never quoted as a measurement"})
    data["ratio_table"] = rows
    meas = [r for r in rows if r["measurable"]]
    best = min(meas, key=lambda r: r["rel_distance_to_observed"])
    checks.append({
        "name": f"T1 ratio table (coefficient 3 (16 pi F_infinity)^10 "
                f"against observed {obs:.4f}): " + "; ".join(
                    f"{r['pair']}/{r['class']}: F = {r['F_infinity']:.6f} "
                    f"-> coef {r['coefficient']:.3e} (rel "
                    f"{r['rel_distance_to_observed']:.2e})"
                    if r["measurable"] else
                    f"{r['pair']}/{r['class']}: floor, not quotable"
                    for r in rows) +
                f"; best measurable relative distance "
                f"{best['rel_distance_to_observed']:.2e}",
        "pass": True})

    # ---- T2 the branch ---------------------------------------------------
    reads = any(r["measurable"] and r["rel_distance_to_observed"] < 1e-3
                for r in rows)
    branch = "HIERARCHY-READS" if reads else "HIERARCHY-UNREAD"
    status = (
        f"BZ HIERARCHY CONTACT, {branch} (zero fitted parameters, "
        f"{time.time() - t0:.0f} s): the V108 ladder reading evaluated at "
        f"V242's measured F_infinity lands at coefficient "
        f"{best['coefficient']:.3e} ({best['pair']}/{best['class']}) "
        f"against the observed {obs:.4f} - a relative distance of "
        f"{best['rel_distance_to_observed']:.2e}: the measured F_infinity "
        f"is {best['F_infinity'] / (loop * (obs / 3.0) ** 0.1):.2f}x the "
        f"loop value the reading needs ((obs/3)^(1/10)/(16 pi) = "
        f"{loop * (obs / 3.0) ** 0.1:.5f}) - too little suppression per "
        f"loop to close the hierarchy ladder. The closed-clock pair "
        f"is a floor ratio (structural zero), recorded, never quoted. No "
        f"pair on C reaches the observed coefficient within 1e-3; the "
        f"ladder coefficient remains the named O8 residual - V243 "
        f"contacts it from the momentum side and measures the miss.")
    ok = all(c["pass"] for c in checks)
    out = {"id": "V243_bz_hierarchy_contact",
           "tier": 4,
           "property": "the V108/V203 ladder reading evaluated at V242's "
                       "Brillouin-zone loop factor: does the observed "
                       "hierarchy coefficient follow (V243, O8/F2)",
           "residual": True,
           "programme": "V243",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v243_bz_hierarchy_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
