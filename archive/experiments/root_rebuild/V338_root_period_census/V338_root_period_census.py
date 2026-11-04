#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V338 - the numbers tribunal for candidate N1 (root-rebuild; booked in
review/pre_registrations_root.md as V337_root_period_census's successor -
the census took the V338 number when V337 became the engine change).

Does the 2I/E8 root FORGE the tribunal numbers (the clock g_obs, the
desert crossing g_cross, the coefficient COEF_OBS, the ladder base LOOP)
as periods/eigenvalues - or only rhyme with them? Dictionary tiers and
guards booked BEFORE the targets were read; V335's lessons adopted
verbatim: exact 50-60 digit mpmath throughout, nonzero-numerator guard,
convergents never evaluated through float64 rounding at large magnitude,
hit-count density printed beside every additive census.
Host arithmetic (mpmath scalars - no GPU path; GPU 1 is engaged by V337).
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)


def main():
    import mpmath
    mpmath.mp.dps = 50
    checks = []
    data = {}

    # booked targets (V206 leaves / config; booked in the pre-registration)
    TARGETS = {
        "g_obs": mpmath.mpf("0.19999798835751165"),
        "g_cross": mpmath.mpf("0.20013316730041675"),
        "COEF_OBS": mpmath.mpf("2.076564111173832"),
        "LOOP": mpmath.mpf("0.019894367886486918"),
    }

    phi = (1 + mpmath.sqrt(5)) / 2
    PI = mpmath.pi

    # ---- D0: geometric tier (π allowed) - booked NON-DISCOVERY tier -----
    best0 = {}
    for name, t in TARGETS.items():
        best = None
        for k in range(-4, 5):
            pip = PI ** k
            for a in range(-4, 5):
                for b in range(-4, 5):
                    for c in range(-4, 5):
                        base = pip * (mpmath.mpf(2) ** a) * (mpmath.mpf(3)
                                                             ** b) \
                            * (mpmath.mpf(5) ** c)
                        for d in range(-2, 3):
                            for f in range(-6, 7):
                                x = base * (mpmath.mpf(7) ** d) * phi ** f
                                if x == 0:
                                    continue
                                m = abs(t / x - 1)
                                if best is None or m < best[0]:
                                    best = (m, k, a, b, c, d, f)
        best0[name] = best
        data["D0_best_%s" % name] = {
            "margin_rel": mpmath.nstr(best[0], 6),
            "form": "pi^%d 2^%d 3^%d 5^%d 7^%d phi^%d" % best[1:]}

    # ---- D1: pi-FREE structural tier (the actual test) -------------------
    best1 = {}
    named1 = []
    for name, t in TARGETS.items():
        best = None
        for a in range(-4, 5):
            for b in range(-4, 5):
                for c in range(-4, 5):
                    for d in range(-2, 3):
                        base = (mpmath.mpf(2) ** a) * (mpmath.mpf(3) ** b) \
                            * (mpmath.mpf(5) ** c) * (mpmath.mpf(7) ** d)
                        for f in range(-6, 7):
                            for j in (0, 1):
                                x = base * phi ** f \
                                    * (mpmath.sqrt(3) ** j)
                                if x == 0:
                                    continue
                                m = abs(t / x - 1)
                                if best is None or m < best[0]:
                                    best = (m, a, b, c, d, f, j)
        best1[name] = best
        data["D1_best_%s" % name] = {
            "margin_rel": mpmath.nstr(best[0], 6),
            "form": "2^%d 3^%d 5^%d 7^%d phi^%d sqrt3^%d" % best[1:]}
        if best[0] <= mpmath.mpf("1e-12"):
            named1.append(name)

    # continued-fraction tier (q <= 1e6, exact evaluation, V335 guards)
    cf = {}
    for name, t in TARGETS.items():
        fr = Fraction(mpmath.nstr(t, 17)).limit_denominator(10 ** 6)
        p = mpmath.mpf(fr.numerator)
        q = mpmath.mpf(fr.denominator)
        resid = abs(t - p / q) if fr.numerator != 0 else None
        cf[name] = (str(fr), mpmath.nstr(resid, 6) if resid is not None
                    else "guarded")
    data["convergents_q1e6"] = cf

    # PSLQ per target (degree <= 6, H <= 20, 60 dp + substitution check)
    mpmath.mp.dps = 60
    pslq_out = {}
    for name, t in TARGETS.items():
        tv = mpmath.mpf(str(TARGETS[name]))
        rel = mpmath.pslq([tv ** k for k in range(7)], maxcoeff=20)
        ok = False
        if rel and max(abs(x) for x in rel) <= 20 and any(rel):
            val = sum(mpmath.mpf(str(c)) * tv ** k
                      for k, c in enumerate(rel))
            ok = abs(val) < mpmath.mpf("1e-30")
        pslq_out[name] = ([int(x) for x in rel] if rel else None, ok)
        if ok:
            named1.append(name + "/pslq")
    mpmath.mp.dps = 50
    data["pslq"] = pslq_out

    # the booked Coxeter probe: is the clock 6/30 to structural precision?
    cox = abs(30 * TARGETS["g_obs"] - 6)
    data["coxeter_probe_abs_30gobs_minus_6"] = mpmath.nstr(cox, 6)

    # n* = 208 additive census over the structural integers, density guard
    ints = [240, 120, 60, 45, 32, 30, 16, 12, 11, 13, 17, 8, 6, 4]
    hits = 0
    example = None
    for i1 in range(len(ints)):
        for i2 in range(len(ints)):
            for i3 in range(len(ints)):
                for c1 in (-4, -3, -2, -1, 1, 2, 3, 4):
                    s2 = c1 * ints[i1]
                    if abs(s2) > 600:
                        continue
                    for c2 in (-4, -3, -2, -1, 1, 2, 3, 4):
                        s3 = s2 + c2 * ints[i2]
                        if abs(s3) > 600:
                            continue
                        for c3 in (-4, -3, -2, -1, 1, 2, 3, 4):
                            if s3 + c3 * ints[i3] == 208:
                                hits += 1
                                if example is None:
                                    example = (c1, ints[i1], c2, ints[i2],
                                               c3, ints[i3])
    data["nstar_208_additive_hits"] = hits
    data["nstar_208_example"] = example

    # the booked non-discovery: LOOP at D0 (16 pi = 2^4 pi, float64 noise)
    loop_trivial = best0["LOOP"][0] < mpmath.mpf("1e-12")

    checks.append(("D1 (pi-free structural tier): margins %s; NAMED at "
                   "1e-12: %s; convergents q<=1e6 residuals %s; PSLQ "
                   "hits %s" % ({k: mpmath.nstr(v[0], 4)
                                 for k, v in best1.items()},
                                named1 or "none", cf,
                                {k: v[1] for k, v in pslq_out.items()}),
                   True))
    checks.append(("booked Coxeter probe: |30 g_obs - 6| = %s (structural "
                   "gate 1e-12: %s); n* = 208 additive census: %d hits "
                   "(density guard - a hit among %d+ proves nothing), "
                   "example %s; LOOP named at D0 (16 pi geometry): %s "
                   "(booked non-discovery)"
                   % (data["coxeter_probe_abs_30gobs_minus_6"],
                      "NO" if cox > mpmath.mpf("1e-12") else "YES",
                      hits, hits, example, loop_trivial), True))
    status = ("ROOT-NUMBERS-NAMED" if named1
              else "ROOT-NUMBERS-UNNAMED")
    data["status"] = status
    for name, passed in checks:
        print(("  ok  " if passed else " FAIL ") + name, flush=True)
    out = {"id": "V338_root_period_census", "tier": 2,
           "property": "does the 2I/E8 root forge the tribunal numbers "
                       "(root-rebuild N1 census)",
           "residual": True, "programme": "root-rebuild V338",
           "branch": status, "data": data,
           "checks": [{"name": n, "pass": p} for n, p in checks],
           "status": status, "pass": True}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v338_root_period_census.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    print("[V338 done] %s" % status, flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
