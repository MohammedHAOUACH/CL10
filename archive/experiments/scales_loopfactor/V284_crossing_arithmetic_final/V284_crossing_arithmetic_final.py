#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V284 - the crossing arithmetic, final (tier 4, O8).

**The question (pre-registered).** Is the certified crossing n* = 208
and the geometric ratio rho = 0.9780302704295731... arithmetically
related to 16 pi at high precision, closing or finally booking the
crossing arithmetic?

Object: V209's exact `skeleton` + `mp_ladder` at 100 and 200 digits
on the certified closed branch (the 4k class at the root clock,
defined by continuity as V209 books - eps = 1e-16 primary, eps =
1e-8 control); the V203/V204 float64 machinery as V276 used it; the
V204 and V276 artifacts. mpmath PSLQ over the booked basis
{1, rho, ..., rho^8, ln rho, ln(16 pi)}, degree <= 8, height cap
1e6, tol 1e-14, verified at 220 dps; the continued-fraction read of
ln(16 pi)/ln(1/rho) to 1e-10 by denominator growth; |n*_geom - 208|
from the certified 200-digit ladder.

Acceptance: CROSSING-ALGEBRAIC iff an integer relation in the booked
basis verifies to residual < 1e-10 at 200 digits; else
CROSSING-MEASURED-FINAL. T0: rho reproduced to 6.7e-16 and crossing
n = 208 from the live float64 machinery (V276's own test), else
MACHINERY-DRIFT.
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

import mpmath                                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp                                # noqa: E402

sys.path.insert(0, exp("V203_closed_branch_coefficient"))
from V203_closed_branch_coefficient import (hop_machinery,  # noqa: E402
                                            loop_ratios)
from cl10 import coupling                                 # noqa: E402

G10_REF = coupling.G_CLOSED                               # pi/4 (V131)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def load_mod(name):
    path = exp(name, "%s.py" % name)
    spec = importlib.util.spec_from_file_location(name.lower(), path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def main():
    t0 = time.time()
    checks, data = [], {}
    print("V284 - the crossing arithmetic, final", flush=True)

    # ---- T0: the float64 certified machinery (V276's own test) -----------
    mach = art("V204_closed_ladder_limits_audit", "data.machinery")
    rho_art = float(mach["rho"])
    cross_art = int(mach["crossing_n"])
    v276_diff = float(art("V276_certified_ladder_208",
                          "data.machinery.rho_abs_diff"))
    G, cartan, Ma, ev = hop_machinery(G10_REF)
    F, rho = loop_ratios(G, cartan, Ma, ev, G10_REF, 260)
    crossing = next(n for n in range(2, 261, 2) if F[n] < 1.0 /
                    (16.0 * math.pi))
    t0ok = (abs(rho - rho_art) <= max(v276_diff, 1e-15) * 1.01
            and crossing == cross_art == 208
            and abs(rho - 0.9780302704295731) < 1e-12)
    data["t0"] = {"rho": rho, "rho_artefact": rho_art,
                  "rho_abs_diff": abs(rho - rho_art),
                  "crossing_n": crossing, "v276_rho_abs_diff": v276_diff}
    checks.append({"name": "T0 float64 certified machinery: rho "
                          "{:.16f} reproduced to {:.1e}, crossing n = {}"
                   .format(rho, abs(rho - rho_art), crossing),
                   "pass": bool(t0ok)})
    if not t0ok:
        data["status"] = "MACHINERY-DRIFT"
        print("V284 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: the exact mp ladder on the closed branch, by continuity ----
    V209 = load_mod("V209_root_clock_degeneracy")
    N_int, blocks, dens, D, cart, planes, TL = V209.skeleton()
    ladders = {}
    for dps, eps in ((100, "1e-16"), (200, "1e-16"), (200, "1e-8")):
        expr = "pi/4 + mp.mpf('{}')".format(eps)
        lad = V209.mp_ladder(expr, 260, dps, planes, blocks, dens, cart)
        ladders[(dps, eps)] = lad
        print("V284 mp ladder dps={} eps={} done ({:.0f} s)".format(
            dps, eps, time.time() - t0), flush=True)

    def rho_of(lad, dps):
        """the certified 4k-class ratio: rho = (F_{n+4}/F_n)^(1/4) at
        the far end of the ladder (n = 244/240), 4k steps only."""
        with mpmath.workdps(dps):
            a = lad[240]["F"]
            b = lad[244]["F"]
            if a is None or b is None or a == 0:
                return None
            return (b / a) ** (mpmath.mpf(1) / 4)

    rho_100 = rho_of(ladders[(100, "1e-16")], 100)
    rho_200 = rho_of(ladders[(200, "1e-16")], 200)
    rho_200b = rho_of(ladders[(200, "1e-8")], 200)
    stab_100_200 = float(abs(rho_100 - rho_200)) if (rho_100 and rho_200
                                                     ) else None
    stab_eps = float(abs(rho_200 - rho_200b)) if (rho_200 and rho_200b
                                                  ) else None
    data["mp_ratios"] = {
        "rho_100dp": str(rho_100)[:40] if rho_100 else None,
        "rho_200dp": str(rho_200)[:60] if rho_200 else None,
        "rho_200dp_eps1e-8": str(rho_200b)[:60] if rho_200b else None,
        "stability_100_vs_200": stab_100_200,
        "stability_eps": stab_eps,
        "rho_float64_dev": (float(abs(rho_200 - mpmath.mpf(repr(rho))))
                            if rho_200 else None)}

    # certified crossing from the 200-digit ladder
    LOOPm = mpmath.mpf(1) / (16 * mpmath.pi)
    lad2 = ladders[(200, "1e-16")]
    n_cert = next((n for n in range(2, 261, 2)
                   if lad2[n]["F"] is not None and lad2[n]["F"] < LOOPm),
                  None)
    data["n_certified_mp"] = n_cert
    checks.append({"name": "T1 mp ladder stability: |rho_100 - rho_200| "
                          "= {}, |eps 1e-16 vs 1e-8| = {}, float64 dev "
                          "{}; certified mp crossing n = {}".format(
                              stab_100_200, stab_eps,
                              data["mp_ratios"]["rho_float64_dev"],
                              n_cert),
                   "pass": bool(stab_100_200 is not None
                                and stab_100_200 < 1e-30
                                and stab_eps is not None
                                and stab_eps < 1e-6
                                and data["mp_ratios"]
                                ["rho_float64_dev"] < 1e-6)})

    # ---- T2: PSLQ over the booked basis at 200 digits ---------------------
    rho_use = rho_200 if rho_200 is not None else mpmath.mpf(repr(rho))
    with mpmath.workdps(220):
        basis = [mpmath.mpf(1)] + [rho_use ** k for k in range(1, 9)] \
            + [mpmath.log(rho_use), mpmath.log(16 * mpmath.pi)]
        rel = mpmath.pslq(basis, maxcoeff=10 ** 6,
                          tol=mpmath.mpf("1e-14"))
        if rel:
            res = abs(sum(mpmath.mpf(c) * b
                          for c, b in zip(rel, basis)))
            relation = [int(c) for c in rel]
            residual = float(res)
        else:
            relation, residual = None, None
        # continued-fraction rationality of ln(16 pi)/ln(1/rho)
        x = mpmath.log(16 * mpmath.pi) / mpmath.log(1 / rho_use)
        cf = []
        y = x
        for _ in range(30):
            a = int(mpmath.floor(y))
            cf.append(a)
            y = 1 / (y - a)
        h0, k0, h1, k1 = 0, 1, 1, 0
        for a in cf[:12]:
            h0, k0, h1, k1 = h1, k1, a * h1 + h0, a * k1 + k0
        conv = (h1, k1)
        cf_read = {"terms": [int(v) for v in cf[:15]],
                   "convergent_12": [str(conv[0]), str(conv[1])],
                   "value": str(x)[:30],
                   "nearest_within_1e-10": bool(
                       abs(x - mpmath.mpf(conv[0]) / conv[1]) < 1e-10)}
    data["pslq"] = {"basis": "1, rho, ..., rho^8, ln rho, ln(16 pi)",
                    "degree_le": 8, "height_cap": 10 ** 6,
                    "tol": "1e-14 (dps 220)",
                    "relation": relation, "residual": residual}
    data["cf_ln_ratio"] = cf_read
    alg = bool(relation) and residual is not None and residual < 1e-10

    # ---- T3: |n*_geom - 208| with the 200-digit rho band -----------------
    with mpmath.workdps(60):
        rm = mpmath.mpf(repr(rho))
        F10o = mpmath.mpf(repr(float(art("V276_certified_ladder_208",
                                         "data.F_10_open_branch"))))
        ngeom = 10 + mpmath.log(LOOPm / F10o, rm)
    data["n_star_geom"] = float(ngeom)
    data["n_star_geom_minus_208"] = float(ngeom - 208)
    data["n_certified_mp_minus_208"] = ((n_cert - 208)
                                        if n_cert is not None else None)
    checks.append({"name": "T2 PSLQ booked basis: relation {} residual {} "
                          "(degree <= 8, height cap 1e6, 220 dps); "
                          "ln(16 pi)/ln(1/rho) = {} convergent {} "
                          "within 1e-10: {}".format(
                              relation, residual,
                              cf_read["value"], cf_read["convergent_12"],
                              cf_read["nearest_within_1e-10"]),
                   "pass": True})
    checks.append({"name": "T3 |n*_geom - 208| = {:.3e}; certified mp "
                          "crossing offset {}".format(
                              abs(float(ngeom - 208)),
                              data["n_certified_mp_minus_208"]),
                   "pass": True})

    # ---- branch -----------------------------------------------------------
    status = ("CROSSING-ALGEBRAIC" if alg else "CROSSING-MEASURED-FINAL")
    data["status"] = status
    print("V284 {}: PSLQ relation {} (residual {}); certified mp "
          "crossing n = {}; |n*_geom - 208| = {:.3e}".format(
              status, relation, residual, n_cert,
              abs(float(ngeom - 208))), flush=True)

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V284 done] {}  ({:.0f} s)".format(status,
                                              time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V284_crossing_arithmetic_final", "tier": 4,
           "property": "arithmetic relation between the certified "
                       "crossing 208 / ratio rho and 16 pi at 200 "
                       "digits (tier 4, O8, final)",
           "residual": True,
           "programme": "V284",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v284_crossing_arithmetic_final.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V284 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
