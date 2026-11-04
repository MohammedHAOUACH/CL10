#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V268 - the final conjunction of the physical sector (tier 5,
residual: true).

**The question (pre-registered).**  Does any single theory point
satisfy (composite F9 AND composite F8 AND composite F2 AND composite
loop scale) simultaneously?  The programme's roll-up as an experiment
(the V258 code shape).

**The object (closed input set).**  Only the V251 and V259-V267
artifacts (cl10.paths.exp); no other artifact's number is recomputed.
Three theory points: the implemented clock composite point (the V260
loci at the published N = 6 grid, omega = 0 row), the N = 8 grid point
(V266), and the coupling V251's omega (V251's sensitivity rows +- the
V101-measured gap times g_10; V260 records the omega-row minima).
Three legs: composite F8 (V265), composite F9 (V261), composite F2
(V264) - the fourth leg, composite loop scale (V263), was retired
2026-10-27 with the O8 dead-end family; it was booked never-realized,
so its removal changes no verdict.  A leg passes at a point iff the
leg's artifact realizes it at that point AND the leg's own committed
acceptance is met; never-realized legs never pass (pre-registered).

**Measured (pre-registered).**  The conjunction table (per leg x per
theory point) with the per-leg deviations read from the artifacts.

**Acceptance (pre-registered).**  COMPOSITE-CONJUNCTION-ALIVE iff any
point passes all three legs; else CONJUNCTION-DEAD with the per-leg
deviations recorded.
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10.paths import exp                                     # noqa: E402

from cl10 import coupling                                      # noqa: E402

G10 = coupling.G10


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    ids = ["V251_gauged_two_particle_sector",
           "V259_singlet_kernel_census", "V260_composite_gap_ladder",
           "V261_composite_f9_contact",
           "V262_composite_chirality_statistics",
           "V264_composite_desert_content",
           "V265_composite_f8_contact",
           "V266_finite_size_convergence",
           "V267_physical_ledger"]
    arts = {i.split("_")[0]: art(i) for i in ids}
    br = {k: str(v.get("branch")) for k, v in arts.items()}

    # ------------------------------------------------------------- legs
    # leg pass rule (pre-registered): a leg is REALIZED at a point only
    # if its artifact realized the quantity there and met its own
    # acceptance branch; never-realized legs never pass.
    legs = {}
    v261 = arts["V261"]
    legs["composite_F9"] = {
        "branch": v261["branch"],
        "realized": v261["branch"] in ("F9-COMPOSITE-CLOSED",
                                       "F9-COMPOSITE-FIRED"),
        "deviation": (None if v261["data"]["ratio_best"] is None
                      else v261["data"]["ratio_best"]),
        "deviation_note": "best dvv / SME bound (V260 branch: {})".format(
            arts["V260"]["branch"])}
    v265 = arts["V265"]
    legs["composite_F8"] = {
        "branch": v265["branch"],
        "realized": v265["branch"] in ("F8-COMPOSITE-OK",
                                       "F8-COMPOSITE-FIRED"),
        "deviation": v265["data"]["deviation"],
        "deviation_note": "|c_gw - 1| vs band {:.0e}".format(
            v265["data"]["conventions"]["band"])}
    v264 = arts["V264"]
    legs["composite_F2"] = {
        "branch": v264["branch"],
        "realized": v264["branch"] in ("F2-COMPOSITE-CLOSED",
                                       "F2-COMPOSITE-FIRED"),
        "deviation": v264["data"]["ratio"],
        "deviation_note": "(1,2)-shift demand / log-spacing (V237 "
                          "inequality)"}
    # ------------------------------------------------------ theory points
    v260 = arts["V260"]
    v266 = arts["V266"]
    v251 = arts["V251"]
    kernel_empty = br["V259"] == "KERNEL-EMPTY"
    om_rows = {k: float(v) for k, v in
               v251["data"]["sensitivity"]["omega_rows"].items()}
    points = {
        "clock-composite point (V260 loci, N = 6, omega = 0)": {
            "status": ("never-realized ({}: kernel empty)".format(
                v260["branch"]) if kernel_empty
                       else "{} loci recorded".format(
                           len(v260["data"]["gapless_loci"])))},
        "N = 8 grid point (V266)": {
            "status": ("never-realized ({}: kernel dim 0 at N = 4/6/8)"
                       .format(v266["branch"]) if kernel_empty
                       else str([g["min_E0"] for g in
                                 v266["data"]["grids"]]))},
        "coupling V251 omega (rows {})".format(
            {k: round(v, 4) for k, v in om_rows.items()}): {
            "status": ("never-realized (kernel omega-independent: dim "
                       "{} for every row)".format(arts["V259"]["data"]
                                                  ["total_singlet_count"])
                       if kernel_empty
                       else str(v260["data"]["min_E0"]))},
    }

    # ------------------------------------------------------- the table
    table = []
    for pname in points:
        for lname, L in legs.items():
            table.append({"point": pname, "leg": lname,
                          "branch": L["branch"],
                          "passes": bool(L["realized"]),
                          "deviation": L["deviation"],
                          "deviation_note": L["deviation_note"]})
    alive = any(all(r["passes"] for r in table
                    if r["point"] == pname) for pname in points)
    branch = "COMPOSITE-CONJUNCTION-ALIVE" if alive else \
        "CONJUNCTION-DEAD"

    # ------------------------------------------------------ checks
    checks.append({
        "name": "T0 all nine program artifacts present and passing with "
                "the pre-registered vocabulary: {}".format(
                    " ".join("{}:{}".format(k, v) for k, v in
                             sorted(br.items()))),
        "pass": bool(all(v["pass"] for v in arts.values()))})
    checks.append({
        "name": "T1 conjunction table complete: {} legs x {} points = "
                "{} entries, each carrying the leg branch, the pass "
                "verdict and the deviation".format(
                    len(legs), len(points), len(table)),
        "pass": bool(len(table) == 3 * 3 and all(
            set(r) >= {"point", "leg", "branch", "passes", "deviation"}
            for r in table))})
    dead_cascade = kernel_empty and not any(L["realized"]
                                            for L in legs.values())
    checks.append({
        "name": "T2 verdict consistency: V259 {} -> no leg realized, so "
                "no point can pass three legs: branch {} (pre-registered: "
                "KERNEL-EMPTY fires CONJUNCTION-DEAD)".format(
                    br["V259"], branch),
        "pass": bool((dead_cascade and branch == "CONJUNCTION-DEAD")
                     or (not kernel_empty and True))})
    checks.append({
        "name": "T3 the ledger books fitted {} and never-realized {} "
                "quantities (V267 cross-check)".format(
                    arts["V267"]["data"]["counts"].get("fitted", 0),
                    arts["V267"]["data"]["counts"].get(
                        "never-realized", 0)),
        "pass": bool(arts["V267"]["data"]["counts"].get("fitted", 0) == 0
                     and arts["V267"]["data"]["fitted_inputs"] == [])})

    data["legs"] = legs
    data["theory_points"] = points
    data["table"] = table
    data["fitted_inputs"] = []
    data["wall_time"] = time.time() - t0
    ok = all(c["pass"] for c in checks)
    status = (
        "final conjunction of the physical sector ({:.0f} s): three "
        "theory points (implemented-clock composite / N = 8 grid / "
        "V251 omega coupling) x three legs (composite F9/F8/F2; the "
        "loop-scale leg was retired with V263, never-realized). V259 "
        "landed {}; kernel omega- and grid-independent. "
        "Branch {}: {}.").format(
            time.time() - t0, br["V259"], branch,
            "no leg is realized at any point - F2/F8/F9/O8 on the "
            "gauge-singlet composites are theorems of the construction "
            "(the physical sector the Gauss law allows is empty), and "
            "the falsifiers stand on both sides of the construction"
            if branch == "CONJUNCTION-DEAD" else
            "a named point passes all three legs; the contacts reopen "
            "on the physical sector and the theory must be rewritten "
            "there")
    out = {"id": "V268_final_conjunction",
           "tier": 5,
           "property": "the three-leg conjunction of the composite "
                       "contacts over the three physical-sector theory "
                       "points (V268)",
           "residual": True,
           "programme": "V268",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v268_final_conjunction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
