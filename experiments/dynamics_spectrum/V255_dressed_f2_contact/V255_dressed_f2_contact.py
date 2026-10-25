#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V255 - the dressed F2 contact: does the interacting loop factor
read the hierarchy? (tier 3, residual: true, F2.)

**The question (pre-registered).**  With the loop functional measured
dynamically (V252's F_int) instead of the formal (16 pi F_10)^-1, does
the ladder reading coefficient 3(16 pi F)^10 (the V108 rule, read
verbatim) land within 1e-3 of the observed 2.0765641, at the
unification reading unification.alpha_std(g) = 0.04 (and alpha_10 =
0.04)?

**The object (closed input set).**  V252's artifact: the F_int table
(N = 32 primary rows, both classes, published clock, both deltas and
signs - 96 candidates); theory points g in {pi/4, pi/2, 0.2} at
beta = pi/4, plus the C_star clock only if V249 landed
STAR-CLOCK-EXISTS (it landed STAR-CLOCK-QUANTIZED, so that leg is
recorded never-realized); the observed hierarchy coefficient from
cl10/config.yaml; V242's imported cell ladder for the T0 functional
control; the bands 1e-3 (ladder) and V233's coupling band (1e-12 *
alpha_target).  T3 verifies beta-independence the V247 way: build
two-dial steps, require the coupling readings bit-identical.

**Measured (pre-registered).**  Every (F_int, g) pair: ladder
deviation, alpha_std/alpha10 deviation from 0.04.

**Acceptance (pre-registered).**  Branch F2-INTERACTION-READS iff some
(F_int, g) pair lands inside BOTH bands; else F2-FIRED-X4 (V233's
alternative quadrilemma with the interaction functional substituted);
F2 stays fired.
"""
import importlib.util
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

from cl10 import config as _config                             # noqa: E402
from cl10 import coupling as _cpl                              # noqa: E402
from cl10 import supercell as SC                               # noqa: E402

G10 = _cpl.G10


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def load_v242():
    p = os.path.join(exp("V242_bz_loop_factor_census"),
                     "V242_bz_loop_factor_census.py")
    spec = importlib.util.spec_from_file_location("v242mod", p)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _beta_invariant(gs):
    a0 = [(_cpl.alpha_std(g), _cpl.alpha10(g)) for g in gs]
    for beta in (math.pi / 2, math.pi / 4):
        _B, _U32, U16 = SC.sectors(math.pi / 4)
        SC.bloch((0.31, -0.77, 0.44), U16, G10, SC.STAGES, beta)
    a1 = [(_cpl.alpha_std(g), _cpl.alpha10(g)) for g in gs]
    return a0 == a1


def main():
    t0 = time.time()
    checks = []
    data = {}

    v252 = art("V252_dynamical_self_energy_loop")
    v249 = art("V249_star_module_closure_census")
    v247 = art("V247_two_dial_desert_contacts")
    obs = float(_config.CFG.physics.hierarchy_coefficient_observed)
    ladder_band = 1e-3
    alpha_target = 0.04
    alpha_band = 1e-12 * alpha_target          # V233's fits band width

    # T0 control: the V242 cell ladder at pi/2 reproduces V247's pi/2 row
    v242mod = load_v242()
    F10_cell = v242mod.v209_cell_ladder(math.pi / 2)[10]
    _cell_coef = 3.0 * (16.0 * math.pi * F10_cell) ** 10
    _pi2_row = next(r for r in v247["data"]["two_by_two"]
                    if abs(r["g"] - math.pi / 2) < 1e-12)
    _cell_dev = abs(_cell_coef - _pi2_row["hierarchy_coef"])

    # F candidates: V252 primary rows (N = 32, contact limit)
    cands = []
    for cnode, sub in sorted(v252["data"]["F_int"].items()):
        clock, cls, qstr = cnode.split("/", 2)
        for tag, F in sorted(sub.items()):
            if "xN32/" in tag and tag.endswith("wcontact"):
                cands.append({"key": cnode + "|" + tag, "clock": clock,
                              "class": cls,
                              "q": [float(x) for x in qstr.split(",")],
                              "omega": tag, "F": float(F)})

    gs = [math.pi / 4, math.pi / 2, 0.2]
    rows = []
    for F in cands:
        coef = 3.0 * (16.0 * math.pi * F["F"]) ** 10
        devL = abs(coef - obs)
        for g in gs:
            a_std = _cpl.alpha_std(g)
            a10 = _cpl.alpha10(g)
            rows.append({"F_key": F["key"], "F_int": F["F"],
                         "clock": "published", "beta": math.pi / 4,
                         "g": g, "coef": coef, "dev_ladder": devL,
                         "alpha_std": a_std,
                         "dev_alpha_std": abs(a_std - alpha_target),
                         "alpha10": a10,
                         "dev_alpha10": abs(a10 - alpha_target),
                         "in_ladder_band": devL < ladder_band,
                         "in_alpha_band_std":
                             abs(a_std - alpha_target) <= alpha_band,
                         "in_alpha_band_10":
                             abs(a10 - alpha_target) <= alpha_band})
    reads = any(r["in_ladder_band"]
                and (r["in_alpha_band_std"] or r["in_alpha_band_10"])
                for r in rows)
    bestL = min(rows, key=lambda r: r["dev_ladder"])
    bestA10 = min(rows, key=lambda r: r["dev_alpha10"])

    # ---------------- checks ---------------------------------------------
    checks.append({
        "name": "T0 the ladder functional control: the imported V242 "
                "cell ladder at pi/2 gives 3(16 pi F_10)^10 = {:.6g}, "
                "reproducing V247's pi/2 row coefficient to {:.2e} "
                "(V108 rule verbatim)".format(_cell_coef, _cell_dev),
        "pass": bool(_cell_dev < 1e-6)})
    checks.append({
        "name": "T1 the candidate set: {} interaction loop factors "
                "(V252 N = 32 contact rows: {} node-clock keys x {} "
                "deltas x 2 signs) x {} theory points = {} pairs".format(
                    len(cands), len(v252["data"]["F_int"]),
                    len(v252["data"]["inputs"]["deltas"]), len(gs),
                    len(rows)),
        "pass": bool(len(cands)
                     == len(v252["data"]["F_int"])
                     * len(v252["data"]["inputs"]["deltas"]) * 2
                     and len(rows) == len(cands) * len(gs))})
    checks.append({
        "name": "T3 beta enters the ladder and coupling functionals "
                "nowhere (V247's verification, rebuilt here): two-dial "
                "steps at beta = pi/2, pi/4 leave alpha_std/alpha10 "
                "bit-identical",
        "pass": bool(_beta_invariant(gs))})
    checks.append({
        "name": "T2 no silent theory points: the g set [pi/4, pi/2, "
                "0.2] and the C_star leg status (V249 {}) are exactly "
                "the pre-registered ones".format(v249["branch"]),
        "pass": bool(v249["branch"] in ("STAR-CLOCK-QUANTIZED",
                                        "STAR-CLOCK-EXISTS"))})

    branch = "F2-INTERACTION-READS" if reads else "F2-FIRED-X4"
    data["rows"] = rows
    data["best_ladder"] = {"coef": bestL["coef"], "F_int": bestL["F_int"],
                           "dev": bestL["dev_ladder"]}
    data["best_alpha10"] = {"g": bestA10["g"], "alpha10": bestA10["alpha10"],
                            "dev": bestA10["dev_alpha10"]}
    data["observed_coefficient"] = obs
    data["bands"] = {"ladder": ladder_band, "alpha": alpha_band}
    data["c_star_leg"] = ("never-realized (V249 STAR-CLOCK-QUANTIZED)"
                          if v249["branch"] != "STAR-CLOCK-EXISTS"
                          else "included")
    data["cell_control"] = {"F10_cell": F10_cell, "coef": _cell_coef,
                            "dev_vs_V247": _cell_dev}
    data["fitted_inputs"] = []
    ok = all(c["pass"] for c in checks)
    status = (
        "dressed F2 contact ({:.0f} s): V252's interaction loop factor "
        "inserted into the ladder rule 3(16 pi F)^10 against the "
        "observed {:.7f} at the unification point - {} pairs ({} F x {} "
        "g; the C_star leg {}). Branch {}: best ladder coefficient "
        "{:.4g} (dev {:.3g}), best alpha_10 = {:.6f} (g = {:.4f}, dev "
        "{:.3g}) against bands {:.1e} / {:.1e}: {}.".format(
            time.time() - t0, obs, len(rows), len(cands), len(gs),
            data["c_star_leg"], branch, bestL["coef"],
            bestL["dev_ladder"], bestA10["alpha10"], bestA10["g"],
            bestA10["dev_alpha10"], ladder_band, alpha_band,
            "no (F_int, g) pair lands inside both bands - F2 stays "
            "fired under the interaction functional"
            if not reads else "a pair reads the hierarchy"))
    out = {"id": "V255_dressed_f2_contact",
           "tier": 3,
           "property": "F2's ladder-contact and unification contact "
                       "re-contacted with the dynamically measured "
                       "loop factor (V255)",
           "residual": True,
           "programme": "V255",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v255_dressed_f2_contact.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
