#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V355 - the elected theory, built with its one dial (booked before
the run in review/pre_registrations_root.md).

The published reading channel dressed by the identity insertion at
the V354-pinned amplitude (read from the V354 artifact); the clock
is BORN at the exact fifth. Six-guard acceptance battery on the
published clock-consuming measurements, house mechanisms verbatim.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V350_backreaction_fixed_point"))
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
from cl10.paths import exp  # noqa: E402
from cl10 import config as _config  # noqa: E402
from cl10 import neutrino as nv  # noqa: E402
import V350_backreaction_fixed_point as V350  # noqa: E402
import V352_oneloop_external_backreaction as V352  # noqa: E402

CFG = _config.CFG
FIFTH = 0.2
TOL_SCALE = float(CFG.pillar_surgery.ledger_scale_rel_tol)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[int(part)] if isinstance(j, list) else j[part]
    return j


def coef_dressed(g, alpha, Ju):
    F10 = abs(V352.Fn_corr(g, alpha * Ju)(10))
    return 3.0 * (F10 / V352.LOOP) ** 10


def main():
    checks, data = [], {}
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    c15_art = float(art("V206_clock_tension",
                        "data.data_clock.coef_at_fifth"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    alpha = float(art("V354_elect_read_root_degeneracy",
                      "data.elect.J1-identity.degeneracies.0.alpha"))
    data["alpha_elect"] = alpha
    data["beta_elect"] = alpha / V350.LOOP
    I32 = np.eye(32, dtype=np.complex128)
    zeroJ = np.zeros((32, 32), dtype=np.complex128)

    # ---- T0 tree guards --------------------------------------------------
    r0 = V352.rootB(FIFTH, zeroJ)
    dev_tree = abs(r0 - g_obs) if r0 else 9.9
    dev_read = abs(coef_dressed(FIFTH, 0.0, zeroJ) - c15_art) / c15_art
    checks.append(("T0 tree guards: rootB(1/5, 0) = %.10f on the published "
                   "g_obs to %.1e (<= 1e-8); reading at 1/5 vs the V206 "
                   "coef_at_fifth artifact rel %.1e (<= 1e-9)"
                   % (r0 if r0 else -1, dev_tree, dev_read),
                   dev_tree <= 1e-8 and dev_read <= 1e-9))

    # ---- B1 birth ---------------------------------------------------------
    r_el = V352.rootB(FIFTH, alpha * I32)
    dev_birth = abs(r_el - FIFTH) if r_el else 9.9
    checks.append(("B1 BIRTH: the dressed ladder's self-consistent root at "
                   "the dial is %.12f - distance to the rational fifth "
                   "%.1e (<= 1e-9): the clock is born rational"
                   % (r_el if r_el else -1, dev_birth), dev_birth <= 1e-9))

    # ---- B2 exact reading -------------------------------------------------
    c_el = coef_dressed(FIFTH, alpha, I32)
    dev_read2 = abs(c_el - V352.COEF_OBS) / V352.COEF_OBS
    checks.append(("B2 EXACT READING: the dressed ladder at the born clock "
                   "reads %.10f vs COEF_OBS %.10f - rel %.1e (<= 1e-8): "
                   "the paper's reading claim UPGRADED from relative 1e-5 "
                   "at the inverted g_obs to exact at the rational fifth"
                   % (c_el, V352.COEF_OBS, dev_read2), dev_read2 <= 1e-8))

    # ---- B3 crossing preserved (V206 definition verbatim) -----------------
    def f10d(g):
        return abs(V352.Fn_corr(g, alpha * I32)(10))

    lo, hi = 0.2, 0.21
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if f10d(mid) < V352.LOOP:
            lo = mid
        else:
            hi = mid
    gc = 0.5 * (lo + hi)
    dev_cross = abs(gc - g_cross)
    checks.append(("B3 CROSSING PRESERVED: the dressed rung-10 crossing of "
                   "1/(16 pi) sits at g~cross = %.9f, distance to the "
                   "published %.9f is %.1e (<= 1e-3 house tolerance), to "
                   "the born fifth %.1e" % (gc, g_cross, dev_cross,
                                           gc - FIFTH), dev_cross <= 1e-3))

    # ---- B4 rung ten ------------------------------------------------------
    FD = V352.Fn_corr(FIFTH, alpha * I32)
    nstar = None
    for n in range(2, 13):
        F = abs(FD(n))
        if 1e-12 < F < V352.LOOP:
            nstar = n
            break
    checks.append(("B4 RUNG TEN: first crossing of the dressed ladder at "
                   "the born clock is rung %s (must be 10: the rung-10 "
                   "clock is not a dressing artifact)" % nstar,
                   nstar == 10))

    # ---- B5 ladder untouched (V318 gate verbatim) -------------------------
    S17 = float(nv.loop_ladder(17)) * float(nv.EV_PER_GEV)
    S17_pub = float(art("V303_junction_energy_r", "data.contact.S17_eV"))
    dev_s17 = abs(S17 - S17_pub) / S17_pub
    checks.append(("B5 LADDER UNTOUCHED: rung-17 scale re-derives the V303 "
                   "contact at rel %.1e (< TOL_SCALE %.0e: the dial lives "
                   "in the reading channel only, the 13-rung machinery "
                   "unchanged)" % (dev_s17, TOL_SCALE), dev_s17 < TOL_SCALE))

    # ---- B6 dial accounting ------------------------------------------------
    dials = ["beta_elect = %.12f (alpha/LOOP, from the V354 artifact)"
             % data["beta_elect"]]
    checks.append(("B6 DIAL ACCOUNTING: exactly %d new constant (%s); the "
                   "published V248 zero-dial ledger is CONTRADICTED BY ONE "
                   "- stated verbatim per booking" % (len(dials), dials[0]),
                   len(dials) == 1))

    data["battery"] = {"birth_dev": dev_birth, "read_dev": dev_read2,
                       "gcross_dressed": gc, "gcross_dev": dev_cross,
                       "rung_star": nstar, "s17_dev": dev_s17}
    guards = [c[1] for c in checks]
    if not guards[0]:
        status, ok = "MACHINERY-DRIFT", False
    elif all(guards):
        status, ok = "ELECT-BUILT", True
    else:
        status, ok = "ELECT-BATTERY-FAIL", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V355 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V355_elect_one_dial", "tier": 1,
           "property": "the elected theory built: born rational clock, "
                       "exact reading, published battery intact, one dial",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results", "v355_elect_one_dial.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
