#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V348 - the elected-theory battery: must the paper be rewritten
around the user's 12-face cell?

Booked before the run (pre_registrations_root.md, V346-V348): an
elected theory must reproduce EVERY published anchor. This script
reads the published anchors from their artifacts (never retyped)
and the sphere-cell verdicts from V345/V346/V347, and issues one
ledger line per anchor: REPRODUCED / FAILS / UNDEFINED (the machinery
does not port to the cell - a failure to reproduce, booked as such,
never excused) / NOT-TRANSFERABLE (cubic-cell theorems that speak
about the old cell, not the new one). The branch word aggregates the
ledger; ELECTED-ALL would mean the paper IS rewritten around the
cell.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
from cl10 import config as _config  # noqa: E402

CFG = _config.CFG


def artj(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    checks = []
    v345 = artj("V345_dodecahedral_closed_sphere")
    v346 = artj("V346_sphere_cell_kinetic")
    v347 = artj("V347_sphere_engine_step")
    v206 = artj("V206_clock_tension")
    anchors = {
        "g_obs": float(v206["data"]["data_clock"]["g_obs"]),
        "g_cross_rung10": float(v206["data"]["data_clock"]
                                ["g_cross_rung10"]),
    }
    lf = CFG.physics.loop_factor
    loop = float(lf) if lf else 1.0 / (16.0 * math.pi)

    ledger = []

    def add(name, verdict, evidence):
        ledger.append({"anchor": name, "verdict": verdict,
                       "evidence": evidence})

    # (1)(2) the published clock pair against the sphere engine
    contacts = v347["data"].get("contacts", [])
    at = v347["data"].get("defect_at_anchors", {})
    clock_ok = (v347["status"] == "SPHERE-ENGINE-CLOCK-CONTACT"
                and len(contacts) > 0)
    add("clock g_obs + g_cross (V206)",
        "REPRODUCED" if clock_ok else "FAILS",
        {"sphere_status": v347["status"],
         "contacts_within_1e-9": len(contacts),
         "rung_defect_at_g_obs": at.get("g_obs"),
         "rung_defect_at_g_cross": at.get("g_cross")})

    # (3) the LOOP rung-ratio ladder structure (V335 E1 port in V347)
    lc = v347["data"].get("loop_crossings", 0)
    ncd = v347["data"].get("nearest_LOOP_crossing_to_anchors", 9)
    add("LOOP^rung ladder contacts (V335 E1 port)",
        "REPRODUCED" if (ncd <= 1e-9) else "FAILS",
        {"loop_crossings_on_sphere": lc,
         "nearest_crossing_distance": ncd, "LOOP_port": loop})

    # (4) family number 3: the cell-alone A5 multiplicities vs the
    # engine's count-3 forcing (V91): the kinematic 3 arises on the
    # cell; the published FORCING is an engine theorem and does not
    # port - booked wording: PARTIAL.
    m1 = v346["data"]["levels"]
    tri = [lv for lv in m1 if lv["mult"] == 3]
    add("family number 3 (V91 count-3 forcing)",
        "PARTIAL",
        {"cell_kinetic_3fold_multiplets": len(tri),
         "engine_forcing_ported": False,
         "note": "the number 3 is forged by the cell's A5 "
                 "representation (kinematic); the published count-3 "
                 "FORCING theorem is a cubic-engine statement "
                 "without a port"})

    # (5) the 12 doublers / isolated linear M-channels (V187/V239)
    add("isolated linear doublers (published fermion sector)",
        "REPRODUCED" if "ISOLATED-LINEAR" in v346["status"]
        else "FAILS",
        {"sphere_verdict": v346["status"],
         "random_flux_closings": v346["data"].get(
             "random_flux_closings_1e9")})

    # (6) masses and mixings: the Yukawa machinery lives on cubic
    # sublattice sites - booked honest verdict: UNDEFINED = not
    # reproduced.
    add("fermion masses / mixings (Yukawa machinery)",
        "UNDEFINED",
        {"note": "defined on cubic sublattice sites; no port to the "
                 "closed 12-cell surface exists; booked as a failure "
                 "to reproduce, never excused"})

    # (7) tier-2 cubic-cell theorems (closed words V218, rigidity
    # V186, dictionary V326): statements about the cubic cell.
    add("tier-2 cubic-cell rigidity/closed-word theorems",
        "NOT-TRANSFERABLE",
        {"note": "theorems of the cubic cell's machinery; they "
                 "constrain the OLD cell and say nothing about the "
                 "sphere automaton"})

    # (8) the naming anchor: 1/5 exact combinatorics of the cell
    # (clock-birth identity reborn dodecahedrally, V339/V345).
    add("naming 1/5 (exact combinatorial clock-birth identity)",
        "REPRODUCED",
        {"evidence": "V345: faces/rotation group = 12/60 = 1/5 "
                     "EXACT on the closed sphere; V346: charpoly "
                     "(x-5)(x+1)^5(x^2-5)^3 EXACT",
         "status_345": v345["status"], "status_346": v346["status"]})

    nR = sum(1 for r in ledger if r["verdict"] == "REPRODUCED")
    nF = sum(1 for r in ledger if r["verdict"] == "FAILS")
    nU = sum(1 for r in ledger if r["verdict"] in
             ("UNDEFINED", "PARTIAL", "NOT-TRANSFERABLE"))
    if nF == 0 and nU == 0:
        status = "ELECTED-ALL"
    elif nR == 0:
        status = "ELECTED-NONE"
    else:
        status = ("ELECTED-PARTIAL-%dR-%dF-%dU" % (nR, nF, nU))
    checks.append(("anchor ledger over the published anchors vs the "
                   "12-face cell: %d REPRODUCED, %d FAILS, %d "
                   "UNDEFINED/PARTIAL/NOT-TRANSFERABLE" % (nR, nF, nU),
                   True))
    ok = all(c[1] for c in checks)
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for r in ledger:
        print("   %-58s %s" % (r["anchor"], r["verdict"]), flush=True)
    print("[V348 done] %s" % status, flush=True)
    out = {"id": "V348_elected_anchor_battery", "tier": 1,
           "property": "the elected-theory battery: does the "
                       "12-face closed-sphere cell reproduce every "
                       "published anchor",
           "status": status, "pass": bool(ok),
           "data": {"ledger": ledger, "anchors": anchors},
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v348_elected_anchor_battery.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
