#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V344 - turn-unit root-clock census of the PUBLISHED ladder (V209's
exact mp_ladder machinery, no new object).

Framing booked in the pre-registration: Lindemann-Weierstrass forbids
any trig-polynomial machinery in g-RADIANS from resonating at the
rational g = 1/5 (sin 0.2 is transcendental); a rational clock can be
born only in TURN units (phase e^{2 pi i g}). V209 measured the
ladder's root-clock set as {8g in 2 pi Z} = denominator 8. The desert
has denominator 5. Census: every rational turn-clock with denominator
q in {3,4,5,6,8,10,12} (trig algebraic there - the 72-degree clock
lives in Q(phi), NOT forbidden), plus the g = 1/5 radian control that
must NOT resonate. Divergence flag: |F_2(1e-5)| / |F_2(1e-1)| > 1e8
(the pi/4 control: 7.95e6 -> 1.15e52).
"""
import json
import math
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V209_root_clock_degeneracy"))

from cl10 import config as _config  # noqa: E402

CFG = _config.CFG


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    import mpmath as mp
    from V209_root_clock_degeneracy import mp_ladder, skeleton
    t0 = time.time()
    checks = []
    data = {}
    N, blocks, dens, D, cart, planes, TL = skeleton()
    NMAX = 10                     # booked: screen level (V209 leaves ran
    DPS = 55                      # nmax 12, 55 dps; screen is cheaper)

    def F2(angle_str):
        r = mp_ladder(angle_str, NMAX, DPS, planes, blocks, dens, cart)
        return r[2]["F"]

    # ---- T0: V209 detuned leaves replayed exactly ----------------------
    v209 = art("V209_root_clock_degeneracy",
               "data.branch_structure.detuned")
    rel = {}
    for key in ("1e-1", "1e-5"):
        booked = float(v209[key]["2"])
        got = abs(float(F2("mp.pi/4 + mp.mpf('%s')" % key)))
        rel[key] = abs(got - booked) / booked
    t0ok = all(v < 1e-12 for v in rel.values())
    ctrl_ratio = (abs(float(F2("mp.pi/4 + mp.mpf('1e-5')")))
                  / abs(float(F2("mp.pi/4 + mp.mpf('1e-1')"))))
    t0_name = ("T0 V209 detuned F_2 leaves replayed EXACTLY through this "
               "lens (rel %s); pi/4 control ratio |F_2(1e-5)|/|F_2(1e-1)| "
               "= %.3e" % ({k: "%.1e" % v for k, v in rel.items()},
                           ctrl_ratio))
    checks.append((t0_name, t0ok))
    data["t0_rel"] = rel

    # ---- the census ----------------------------------------------------
    angles = []
    for q in (3, 4, 5, 6, 8, 10, 12):
        for p in (1, 2, 3):
            if p >= q or math.gcd(p, q) != 1:
                continue
            if (q, p) in ((4, 2),):
                continue
            angles.append((q, p))
    rows = {}
    for (q, p) in angles:
        base = ("mp.mpf('2')*mp.pi*mp.mpf('%d')/mp.mpf('%d')" % (p, q))
        big = abs(float(F2(base + " + mp.mpf('1e-5')")))
        small = abs(float(F2(base + " + mp.mpf('1e-1')")))
        f0 = F2(base)
        at0 = None if f0 is None else abs(float(f0))
        ratio = big / small if small > 0 else float("inf")
        rows["2pi*%d/%d" % (p, q)] = {
            "F2_1e-1": small, "F2_1e-5": big,
            "F2_exact": ("EXACT-0/0" if at0 is None else at0),
            "ratio": ratio, "divergent": bool(ratio > 1e8)}
    # the radian control: g = 1/5 (algebraic nonzero: trig transcendental)
    ctrl = {}
    for d in ("1e-1", "1e-5"):
        ctrl[d] = abs(float(F2("mp.mpf('0.2') + mp.mpf('%s')" % d)))
    ctrl["ratio"] = ctrl["1e-5"] / ctrl["1e-1"]
    ctrl["divergent"] = bool(ctrl["ratio"] > 1e8)
    rows["RADIAN-CONTROL 0.2"] = ctrl

    diverging = [k for k, v in rows.items() if v["divergent"]]
    q5 = [k for k in diverging if any(("%d/%d" % (p, q)) in k
                                      for (q, p) in angles
                                      if q in (5, 10))]
    checks.append(("census %d turn-clocks + radian control: diverging set "
                   "= %s; q in {5,10} diverging: %s; denominator-8 family "
                   "all diverging: %s; radian 0.2 control diverging: %s"
                   % (len(rows), diverging, q5,
                      all(rows.get("2pi*%d/%d" % (p, q), {}).get(
                          "divergent", True)
                          for (q, p) in angles if q in (4, 8)
                          and "2pi*%d/%d" % (p, q) in rows),
                      ctrl["divergent"]), True))
    data["rows"] = rows
    data["wall_s"] = time.time() - t0

    five_family = bool(q5)
    status = ("TURN-CLOCK-5" if five_family else
              "TURN-CLOCK-VOID" if (rows["2pi*1/4"]["divergent"] and
                                    not ctrl["divergent"]) else
              "TURN-CLOCK-UNEXPECTED")
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for k, v in rows.items():
        print("   %-18s ratio %12.3e  F2(1e-5) %10.3e  exact %s  %s"
              % (k, v["ratio"], v["F2_1e-5"],
                 v["F2_exact"] if isinstance(v["F2_exact"], str)
                 else "%.3e" % v["F2_exact"],
                 "DIVERGES" if v["divergent"] else "-"), flush=True)
    print("[V344 done] %s  wall %.0f s" % (status, data["wall_s"]),
          flush=True)
    out = {"id": "V344_turn_unit_clock_census", "tier": 2,
           "property": "does the desert's denominator 5 occur as a root "
                       "clock of the published ladder in TURN units "
                       "(Lindemann-Weierstrass framing)",
           "status": status, "pass": bool(all(c[1] for c in checks)),
           "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    with open(os.path.join(HERE, "results",
                           "v344_turn_unit_clock_census.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
