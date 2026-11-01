#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V340 - the last two birth channels: relational webs and Pisot shadows.

The single-number naming census (V338) failed. Two birth channels
remain open and are closed or opened here, at booked guards:
  W1 pair web: exact integer relations a + b x + c y = 0 between the
      published decimals;
  W2 log web: multiplicative relations x^a y^b = e^c;
  W3 product web: a + b x + c y + d xy = 0;
  P Pisot shadow: is the desert's 2 ppm exactly (p/q) lambda^-n for a
      booked Pisot unit lambda? The ONLY known mechanism whose SHAPE
      matches the desert data: exact rational backbone + an
      exponentially small irrational correction.
"""
import itertools
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    import mpmath as mp
    mp.mp.dps = 60

    T = {
        "g_obs": "0.19999798835751165",
        "g_cross": "0.20013316730041675",
        "gap": "0.000135178942905101",
        "COEF": "2.076564111173832",
        "coef_fifth": "2.088106865906787",
        "C_CORR": "2.07644360080448",
        "LOOP": "0.019894367886486918",
        "F208_root": "0.0196851541811279",
        "F208_321": "0.0196927945051569",
        "F4": "1.82994979751295",
        "F8": "1.67435813070784",
        "r_atm": "1.91911543798088",
        "Wln10": "-3.50320080595748",
        "Wln17": "-2.79531879611527",
        "Wln208": "-0.421699620700671",
        "def_pi2": "0.973756015764892",
        "m_gapmin": "0.000200289529530561",
    }
    X = {k: mp.mpf(v) for k, v in T.items()}
    KNOWN = {frozenset(("g_obs", "g_cross", "gap")):
             "gap := g_cross - g_obs (V309 definition)"}
    checks = []
    edges = []
    for a, b in itertools.combinations(sorted(X), 2):
        pair = frozenset((a, b))
        known = any(pair <= k for k in KNOWN)
        for tag, vec, Hcap, thr in (
                ("W1", [mp.mpf(1), X[a], X[b]], 20000, mp.mpf("1e-13")),
                ("W2", [mp.mpf(1), mp.log(abs(X[a])), mp.log(abs(X[b]))],
                 20000, mp.mpf("1e-13"))):
            c = mp.pslq(vec, maxcoeff=Hcap)
            if c is None:
                continue
            res = abs(sum(ci * vi for ci, vi in zip(c, vec)))
            if res > thr:
                continue
            if tag == "W2":
                c0, c1, c2 = [int(x) for x in c]
                if max(abs(c1), abs(c2)) > 2:
                    continue
                if c0 != 0 and abs(c0) > 3:
                    continue
            edges.append({"tag": tag, "pair": [a, b],
                          "coeffs": [int(x) for x in c],
                          "residual": mp.nstr(res, 5),
                          "known": bool(known)})
        # W3 product web, tighter
        c = mp.pslq([mp.mpf(1), X[a], X[b], X[a] * X[b]], maxcoeff=200)
        if c is not None:
            vec = [mp.mpf(1), X[a], X[b], X[a] * X[b]]
            res = abs(sum(ci * vi for ci, vi in zip(c, vec)))
            if res <= mp.mpf("1e-16"):
                edges.append({"tag": "W3", "pair": [a, b],
                              "coeffs": [int(x) for x in c],
                              "residual": mp.nstr(res, 5),
                              "known": bool(known)})
    unknown = [e for e in edges if not e["known"]]
    checks.append(("W1-W3 webs over %d published decimals (%d pairs, "
                   "guarded heights): %d edges found, %d of them NOT "
                   "definitionally known: %s"
                   % (len(X), len(X) * (len(X) - 1) // 2, len(edges),
                      len(unknown),
                      [(e["tag"], e["pair"], e["coeffs"])
                       for e in unknown]), True))

    # ---- P: Pisot shadows ----------------------------------------------
    polys = {  # DESCENDING coefficients (mp.polyroots convention)
        "golden": [1, -1, -1],
        "silver": [1, -2, -1],
        "plastic": [1, 0, -1, -1],
        "supergolden": [1, -1, 0, -1],
        "quartic-smallest": [1, 0, 0, -1, -1],
        "quintic": [1, 0, 0, 0, -1, -1],
        "tribonacci": [1, -1, -1, -1],
        "quartic-b": [1, -1, 0, 0, -1],
    }
    lambdas = {}
    for nm, co in polys.items():
        rts = mp.polyroots([mp.mpf(x) for x in co])
        lam = max(rts, key=lambda z: abs(z))
        others = [z for z in rts if z != lam]
        pisot = abs(lam) > 1 and all(abs(z) < 1 for z in others)
        lambdas[nm] = (lam, pisot)
    shadows = []
    targets = {
        "g_obs-1/5": mp.mpf("0.19999798835751165") - mp.mpf("0.2"),
        "gap": X["gap"],
        "COEF-2": X["COEF"] - 2,
        "coef_fifth-COEF": X["coef_fifth"] - X["COEF"],
        "m_gapmin": X["m_gapmin"],
    }
    for tn, tv in targets.items():
        At = abs(tv)
        for ln, (lam, pis) in lambdas.items():
            if not pis:
                continue
            v = At
            for n in range(1, 61):
                v *= lam
                if v < mp.mpf("1e-2"):
                    continue
                if v > mp.mpf("1e2"):
                    break
                fr = Fraction(mp.nstr(v, 25)).limit_denominator(100)
                if fr.numerator == 0 or abs(fr.denominator) > 100:
                    continue
                if abs(fr) > 10 ** 4:
                    continue
                res = abs(v - mp.mpf(fr.numerator) / fr.denominator)
                if res <= mp.mpf("1e-12"):
                    shadows.append({"target": tn, "lambda": ln, "n": n,
                                    "pq": "%d/%d" % (fr.numerator,
                                                     fr.denominator),
                                    "residual": mp.nstr(res, 5)})
    checks.append(("P Pisot shadow search (%d targets x 8 booked Pisot "
                   "units, |p|<=1e4, q<=100, n<=60, window and guard as "
                   "booked): %d shadows: %s"
                   % (len(targets), len(shadows), shadows), True))

    web_closed = not unknown
    pisot_void = not shadows
    status = ("WEB-" + ("CLOSED" if web_closed else "EDGE-FOUND")
              + "-PISOT-" + ("VOID" if pisot_void else "FOUND"))
    print("\n".join("  ok %s" % c[0] for c in checks), flush=True)
    print("[V340 done] %s" % status, flush=True)
    out = {"id": "V340_number_web_and_pisot_shadow", "tier": 2,
           "property": "relational web of the published decimals + the "
                       "Pisot-shadow birth of the desert's 2 ppm",
           "status": status, "pass": web_closed and pisot_void,
           "data": {"targets": list(T), "edges": edges,
                    "unknown_edges": len(unknown),
                    "pisot_lambdas": {k: [mp.nstr(v[0], 15), v[1]]
                                      for k, v in lambdas.items()},
                    "shadows": shadows,
                    "known_edges_excluded": list(KNOWN.values())},
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    with open(os.path.join(HERE, "results",
                           "v340_number_web_and_pisot_shadow.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
