#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V256 - does the interaction loop factor preserve the ladder
arithmetic? (tier 4, residual: true.)

**The question (pre-registered).**  The ladder rule of V197/V59 with
the interaction loop factor: loop' = F_int^(1/10) instead of the
published LOOP = (16 pi F_10)^-1 - do the three rung scales (singlet
A rung 3, singlet B rung 10, bare neutrino rung 17, values read from
the V197 artifact) still sit on integer rungs of MP x loop'^n?

**The object (closed input set).**  The V252 artifact: F_int (N = 32
primary rows, 96 candidates); the V197 artifact: ladder_rows (values
and integer rungs) and its rung-arithmetic rule (x = ln(m/M_P)/ln
loop, imported u10pred.M_PLANCK); the band is V197's rung arithmetic
itself (the published rungs reproduce to ~1e-15 on the published
ruler - the joint residual band for a new ruler is 1e-9).

**Measured (pre-registered).**  For each F_int: loop', x' per rung,
residual |x' - rung|, joint max over the three rungs.

**Acceptance (pre-registered).**  Branch LADDER-CONSISTENT iff the
joint residual <= 1e-9 for some F_int; else LADDER-BROKEN with the
best residual - then the published ladder is reported as an artefact
of the formal functional and F2's repair needs more than the loop
factor (the plan's named consequence).
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10.paths import exp                                     # noqa: E402

from cl10 import u10pred                                       # noqa: E402

MP = u10pred.M_PLANCK


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}

    v252 = art("V252_dynamical_self_energy_loop")
    v197 = art("V197_rung_arithmetic")
    rung_rows = [r for r in v197["data"]["ladder_rows"]
                 if r["rung"] in (3, 10, 17)]
    band = 1e-9

    cands = []
    for cnode, sub in sorted(v252["data"]["F_int"].items()):
        for tag, F in sorted(sub.items()):
            if "xN32/" in tag and tag.endswith("wcontact") \
                    and float(F) > 0.0:
                cands.append({"key": cnode + "|" + tag, "F": float(F)})

    # control: published ruler reproduces V197 rungs
    ctrl = max(abs(math.log(r["value_gev"] / MP)
                   / math.log(u10pred.LOOP) - r["rung"])
               for r in rung_rows)

    best = None
    rows = []
    for c in cands:
        lp = c["F"] ** 0.1
        if abs(lp - 1.0) < 1e-12:
            continue
        xs = []
        for r in rung_rows:
            xp = math.log(r["value_gev"] / MP) / math.log(lp)
            xs.append({"rung": r["rung"], "x_new": xp,
                       "residual": abs(xp - r["rung"])})
        joint = max(x["residual"] for x in xs)
        rows.append({"F_key": c["key"], "F_int": c["F"], "loop_new": lp,
                     "rungs": xs, "joint": joint})
        if best is None or joint < best["joint"]:
            best = rows[-1]

    consistent = bool(best is not None and best["joint"] <= band)

    # ---------------- checks ---------------------------------------------
    checks.append({
        "name": "T0 V197 rung-arithmetic control: the published ruler "
                "reproduces rungs 3/10/17 of the V197 artifact values "
                "to max residual {:.2e}".format(ctrl),
        "pass": bool(ctrl < 1e-9)})
    nzero = len(cands) - len(rows)
    checks.append({
        "name": "T1 candidate set: {} positive F_int from V252's N = 32 "
                "contact rows (zero/sign exclusions {} dropped)".format(
                    len(rows), nzero),
        "pass": bool(len(rows) >= 24)})
    checks.append({
        "name": "T2 the ruler is a tenth root, not a free dial: loop' = "
                "F_int^(1/10) per V197's loop-exponent rule (rung "
                "arithmetic imported from u10pred, no rescaling)",
        "pass": bool(all(abs(r["loop_new"] ** 10 - r["F_int"]) < 1e-12
                         for r in rows))})

    branch = "LADDER-CONSISTENT" if consistent else "LADDER-BROKEN"
    data["best"] = best
    data["rows"] = sorted(rows, key=lambda r: r["joint"])[:12]
    data["published_control_residual"] = ctrl
    data["band"] = band
    data["rung_values_gev"] = {r["name"]: r["value_gev"]
                               for r in rung_rows}
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    status = (
        "ladder rung reconstruction with the interaction loop factor "
        "({:.0f} s): loop' = F_int^(1/10), {} candidates against the "
        "V197 rung scales 3/10/17 (published ruler control {:.1e}). "
        "Branch {}: best joint rung residual {:.4f} at F_int = {:.4g} "
        "(loop' = {:.4f}, published LOOP {:.4f}) vs band {:.0e}: {}.".format(
            time.time() - t0, len(rows), ctrl, branch,
            best["joint"] if best is not None else float("nan"),
            best["F_int"] if best is not None else float("nan"),
            best["loop_new"] if best is not None else float("nan"),
            u10pred.LOOP, band,
            "the ladder arithmetic survives the substitution"
            if consistent else
            "integer rungs are lost - the published ladder is an "
            "artefact of the formal functional and F2's repair needs "
            "more than the loop factor"))
    out = {"id": "V256_ladder_rung_reconstruction",
           "tier": 4,
           "property": "the rung arithmetic re-derived with the "
                       "dynamically measured loop factor (V256)",
           "residual": True,
           "programme": "V256",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v256_ladder_rung_reconstruction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
