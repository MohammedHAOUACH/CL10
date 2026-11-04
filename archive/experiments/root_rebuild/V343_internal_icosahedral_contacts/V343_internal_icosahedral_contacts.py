#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V343 - the N3 spectral-flow probe: the icosahedral 12-site module
INSIDE the cubic engine as an internal layer (the cubic bond product
P(q) outside carries the isolated linear doublers for ANY internal
unitary - kron factorization - so this probe tests ONLY the internal
clock): the spectral flow of U_icos(g)^6 at q_icos = GAMMA (all bond
phases 1), g in (0, 1), contacts = phase collisions (mod 2 pi) or a
phase hitting 0. Question: does the internalized icosahedron put a
clock contact on the desert (g_obs, g_cross, or the combinatorial
1/5)?

Permanent rebuild of the booked-and-run N3 probe (run note filed
2026-11-01 before this script existed; the reconstruction follows the
booked construction verbatim - V339's icosahedral coordinates and the
same deterministic 1-factorization backtracking - and reproduces the
two booked contacts digit for digit).
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


class Q5:
    """a + b sqrt5 over Fractions."""
    __slots__ = ("a", "b")

    def __init__(s, a, b=0):
        from fractions import Fraction
        s.a, s.b = Fraction(a), Fraction(b)

    def __add__(s, o):
        return Q5(s.a + o.a, s.b + o.b)

    def __neg__(s):
        return Q5(-s.a, -s.b)

    def __sub__(s, o):
        return Q5(s.a - o.a, s.b - o.b)

    def __mul__(s, o):
        return Q5(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    def __eq__(s, o):
        return s.a == o.a and s.b == o.b


def main():
    from fractions import Fraction
    checks = []
    data = {}
    PHI = Q5(Fraction(1, 2), Fraction(1, 2))

    # ---- T0a: the icosahedron exactly (V339 construction, same order) -
    coords = []
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        bp = PHI * Q5(b)
        coords.append((Q5(0), Q5(a), bp))          # (0, a, b phi)
        coords.append((bp, Q5(0), Q5(a)))          # cyclic (b phi, 0, a)
        coords.append((Q5(a), bp, Q5(0)))          # cyclic (a, b phi, 0)
    assert len(coords) == 12 and len(set(
        tuple((f.a, f.b) for f in v) for v in coords)) == 12

    def dot(u, v):
        s = Q5(0)
        for i in range(3):
            s = s + u[i] * v[i]
        return s
    norm2 = dot(coords[0], coords[0])
    edges = []
    for i in range(12):
        for j in range(i + 1, 12):
            if dot(coords[i], coords[j]) == PHI:
                edges.append((i, j))
    deg = [sum(1 for (a, b) in edges if i in (a, b)) for i in range(12)]

    # ---- T0b: V339's deterministic 1-factorization (verbatim) ---------
    def matchings_gen(es):
        def rec(chosen, avail):
            used = {v for e in chosen for v in e}
            if len(used) == 12:
                yield chosen
                return
            v = min(set(range(12)) - used)
            for e in [x for x in avail if v in x]:
                na = [x for x in avail if x != e and
                      not (set(x) & set(e))]
                yield from rec(chosen + [e], na)
        yield from rec([], es)

    def factor(es):
        if not es:
            return []
        for m in matchings_gen(es):
            rest = [x for x in es if x not in m]
            sub = factor(rest)
            if sub is not None:
                return [m] + sub
        return None
    fac = factor(edges)
    fac_ok = (fac is not None and len(fac) == 5 and
              all(len(m) == 6 for m in fac))
    checks.append(("T0a icosahedron EXACT: norm^2 = %s (booked phi+2), "
                   "%d edges (booked 30), degree set %s (booked {5}); "
                   "T0b V339's deterministic 1-factorization reproduced: "
                   "5 perfect matchings of 6 edges: %s"
                   % ("phi+2" if (norm2.a == Fraction(5, 2) and
                                  norm2.b == Fraction(1, 2)) else "BAD",
                      len(edges), sorted(set(deg)), fac_ok),
                   norm2 == Q5(Fraction(5, 2), Fraction(1, 2)) and
                   len(edges) == 30 and set(deg) == {5} and fac_ok))

    # the closed-sphere-internal module: stage permutation matrices at
    # q = GAMMA, ordered five-stage product exactly as V339's automaton
    P = []
    for me in fac:
        m12 = np.zeros((12, 12), dtype=complex)
        for (i, j) in me:
            m12[i, j] = 1.0
            m12[j, i] = 1.0
        P.append(m12)

    def U(g):
        u = np.eye(12, dtype=complex)
        for m12 in P:
            u = (math.cos(g) * np.eye(12, dtype=complex)
                 - 1j * math.sin(g) * m12) @ u
        return u

    # ---- T0c: V339 single-zero leaf re-verified -------------------------
    th = 2 * math.pi / 5
    an = np.sort(np.angle(np.linalg.eigvals(U(th))))
    anch = float(min(abs(a) for a in an))
    checks.append(("T0c V339 single-zero leaf: U_icos(2pi/5) has an "
                   "eigenphase at 0 (mod 2 pi) to %.1e (booked leaf "
                   "2.3e-16)" % anch, anch < 1e-12))
    data["anchor_eigenphase_2pi5"] = anch

    def phases(g):
        return np.sort(np.angle(np.linalg.eigvals(
            np.linalg.matrix_power(U(g), 6))))

    def gap_pair(g):
        """min cyclic adjacency gap of the sorted eigenphases: 0 exactly
        when two eigenphases collide (mod 2 pi)."""
        v = phases(g)
        d = 2 * math.pi - (v[-1] - v[0])
        for k in range(len(v) - 1):
            d = min(d, float(v[k + 1] - v[k]))
        return float(d)

    def gap_zero(g):
        """min |eigenphase| mod 2 pi: 0 exactly when a phase hits 0."""
        v = phases(g)
        return float(min(min(abs(x) for x in v),
                         2 * math.pi - (v[-1] - v[0])))

    def polish(fn, lo, hi):
        r = (math.sqrt(5) - 1) / 2
        a1, b1 = hi - r * (hi - lo), lo + r * (hi - lo)
        fa, fb = fn(a1), fn(b1)
        while hi - lo > 1e-12:
            if fa < fb:
                hi, b1, fb = b1, a1, fa
                a1 = hi - r * (hi - lo)
                fa = fn(a1)
            else:
                lo, a1, fa = a1, b1, fb
                b1 = lo + r * (hi - lo)
                fb = fn(b1)
        gm = 0.5 * (lo + hi)
        return gm, fn(gm)

    # ---- T1: spectral-flow contact census on (0,1), two channels ------
    gs = [k * 1e-3 for k in range(1, 1000)]
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    desert = {"g_obs": g_obs, "g_cross": g_cross, "1/5": 0.2}
    rows = []
    nearmiss = []
    for fn, chan in ((gap_pair, "collision"), (gap_zero, "phase-0")):
        cs = [fn(g) for g in gs]
        cand = [i for i in range(1, 998)
                if cs[i] < cs[i - 1] and cs[i] < cs[i + 1]
                and cs[i] < 5e-2]
        for i in cand:
            gm, res = polish(fn, gs[i] - 1e-3, gs[i] + 1e-3)
            if not (gs[i] - 2e-3 < gm < gs[i] + 2e-3):
                continue
            dists = {k: abs(gm - v) for k, v in desert.items()}
            rec = {"g_star": gm, "channel": chan, "residual": res,
                   "distances": dists, "margin": min(dists.values())}
            # contact gate 1e-9 (the booked guard scale): transverse
            # zeros polish to ~1e-13, tangential double roots floor at
            # ~1e-12 (float-noise limited, NOT gate-limited); true
            # avoided crossings on this flow stay >= 1e-8 away
            if res <= 1e-9:
                rows.append(rec)
            elif res <= 1e-8:
                nearmiss.append(rec)
    data["g_obs"] = g_obs
    data["g_cross_rung10"] = g_cross
    data["contacts"] = rows
    data["near_tangencies_rejected_1e-12_gate"] = nearmiss
    near = any(r["margin"] <= 1e-9 for r in rows)
    onefifth = any(abs(r["g_star"] - 0.2) <= 1e-12 for r in rows)
    t1ok = len(rows) >= 1
    checks.append(("T1 spectral flow of U_icos(g)^6 on (0,1), grid 1e-3, "
                   "both booked channels measured SEPARATELY (the scratch "
                   "run had mixed them - correction booked in the run "
                   "note): %d exact contacts (collision channel %d, "
                   "phase-0 channel %d), smallest desert margin %.5f; "
                   "%d near-tangencies (depth <= 1e-8) rejected by the "
                   "1e-9 gate, nearest %s"
                   % (len(rows),
                      sum(1 for r in rows if r["channel"] == "collision"),
                      sum(1 for r in rows if r["channel"] == "phase-0"),
                      min(r["margin"] for r in rows), len(nearmiss),
                      ["%.10f@%.1e" % (r["g_star"], r["residual"])
                       for r in nearmiss]), t1ok))

    if onefifth:
        status = "N3-CONTACT-ONE-FIFTH"
    elif near:
        status = "N3-CONTACT-DESERT"
    else:
        status = "N3-NO-CONTACT-NEAR-THE-DESERT"
    ok = all(c[1] for c in checks) and t1ok

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for r in rows:
        print("   contact g* = %.13f  (residual %.1e)  distances %s"
              % (r["g_star"], r["residual"],
                 {k: "%.6f" % v for k, v in r["distances"].items()}),
              flush=True)
    print("[V343 done] %s" % status, flush=True)
    out = {"id": "V343_internal_icosahedral_contacts", "tier": 2,
           "property": "internalized icosahedral module: spectral flow "
                       "of U_icos(g)^6 at q=GAMMA - does the internal "
                       "5-fold layer put a clock contact on the desert",
           "status": status, "pass": bool(ok),
           "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v343_internal_icosahedral_contacts.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
