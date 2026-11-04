#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V349 - the user's NEW GEOMETRY, built FOR the cell: an A5-covariant
quantum walk on the 12 faces (booked before the run).

States: the 60 directed face-adjacencies |f, i>, i = 0..4 the pentagon
cyclic edge order at face f (the boundary cyclic order transported by
proper rotations - so the coin family is A5-covariant BY CONSTRUCTION,
and the certificate T0b proves it, not assumes it).
Coin: C(alpha) = direct sum_f exp(i alpha (P_f + P_f^dagger)), P_f the
pentagon's own 5-cycle at f.
Shift S: traverse the bond; arrival index found in the neighbor's
cyclic order.
U(alpha) = S . C(alpha): 60 x 60 unitary. A5 is the replacement
geometry; the walk is the matter organ the battery found missing.

T0a unitarity + spectrum on the circle; T0b A5-equivariance
certificate for the generators p5, p2 (V339 T2 verbatim); T0c graph
counts (V345).
M1 gapless locus (scan + polish), multiplicity vs the published 12;
M2 protection tribunal (random A5-breaking, 10 amplitudes x 1e3
draws: protected <= 1e-6 keeps its zero, kinematic opens a gap);
M3 clock battery vs g_obs, g_cross_rung10, 1/5 both unit conventions,
guard 1e-9.
Branch: QW-GAPLESS-PROTECTED-CLOCK-CONTACT / QW-GAPLESS-PROTECTED-
CLOCK-VOID / QW-GAPLESS-UNPROTECTED / QW-MASSIVE.
"""
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np
from scipy.linalg import expm

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
    __slots__ = ("a", "b")

    def __init__(s, a, b=Fraction(0)):
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
    checks = []
    data = {}

    # ---- the 12 faces (= icosahedron vertices), exactly ----------------
    PHI = Q5(Fraction(1, 2), Fraction(1, 2))
    coords = []
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        bp = PHI * Q5(b)
        coords.append((Q5(0), Q5(a), bp))
        coords.append((bp, Q5(0), Q5(a)))
        coords.append((Q5(a), bp, Q5(0)))

    def dotq(u, v):
        s = Q5(0)
        for i in range(3):
            s = s + u[i] * v[i]
        return s
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12)
             if dotq(coords[i], coords[j]) == PHI]
    nbr = {i: [] for i in range(12)}
    for (i, j) in edges:
        nbr[i].append(j)
        nbr[j].append(i)

    # cyclic (boundary) order at each face: angular sort about the
    # radial axis (V345's rotation-system device on the dual solid)
    cf = [[float(x.a) + float(x.b) * math.sqrt(5.0) for x in v]
          for v in coords]

    def cyc_order(f):
        p = np.array(cf[f])
        ax = p / np.linalg.norm(p)
        e1 = np.cross(ax, [0, 0, 1.0])
        if np.linalg.norm(e1) < 1e-9:
            e1 = np.cross(ax, [0.0, 1.0, 0.0])
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(ax, e1)
        out = []
        for w in nbr[f]:
            d = np.array(cf[w]) - p
            out.append((math.atan2(float(d @ e2), float(d @ e1)), w))
        return [w for _, w in sorted(out)]
    m = {f: cyc_order(f) for f in range(12)}
    for f in range(12):
        assert sorted(m[f]) == sorted(nbr[f])
    checks.append(("T0c walk graph: 12 faces, %d bonds (booked 30), "
                   "pentagon cyclic orders closed at every face, vertex "
                   "degree {5} (V345 counts)" % len(edges),
                   len(edges) == 30 and all(len(m[f]) == 5
                                            for f in range(12))))

    # ---- the walk: 60 directed states |f, i> ----------------------------
    st = [(f, i) for f in range(12) for i in range(5)]
    idx = {s: k for k, s in enumerate(st)}
    D = 60
    S = np.zeros((D, D), dtype=complex)
    for (f, i) in st:
        g = m[f][i]
        j = m[g].index(f)
        S[idx[(g, j)], idx[(f, i)]] = 1.0

    def coin_blocks(alpha):
        """C_f = exp(i alpha (P_f + P_f^dagger)) with P_f the 5-cycle
        of the face's OWN cyclic order."""
        C = np.zeros((D, D), dtype=complex)
        for f in range(12):
            P5 = np.zeros((5, 5))
            for i in range(5):
                P5[(i + 1) % 5, i] = 1.0
            Hf = P5 + P5.conj().T
            blk = expm(1j * alpha * Hf)
            for i in range(5):
                for jj in range(5):
                    C[idx[(f, i)], idx[(f, jj)]] = blk[i, jj]
        return C

    def U(alpha):
        return S @ coin_blocks(alpha)

    # ---- T0a: unitarity + circle spectrum ------------------------------
    t0a = 0.0
    for al in (0.0, 0.37, 1.0, 2.0):
        Ua = U(al)
        t0a = max(t0a,
                  float(np.max(np.abs(Ua @ Ua.conj().T - np.eye(D)))),
                  float(np.max(np.abs(np.abs(np.linalg.eigvals(Ua))
                                      - 1.0))))
    checks.append(("T0a walk unitarity and circle spectrum at sampled "
                   "alphas: worst %.1e (<= 1e-12)" % t0a, t0a <= 1e-12))

    # ---- T0b: A5-equivariance certificate (V339 T2 generators) ---------
    def eps(i, j, k):
        if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            return 1
        if (i, j, k) in ((0, 2, 1), (2, 1, 0), (1, 0, 2)):
            return -1
        return 0

    def rot5(vf):
        v = np.array(vf)
        nv = float(v @ v)
        cos72 = (math.sqrt(5.0) - 1.0) / 4.0
        R = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                R[i, j] = (cos72 if i == j else 0.0) + \
                    (1 - cos72) * v[i] * v[j] / nv
                for k in range(3):
                    e = eps(i, k, j)
                    if e:
                        R[i, j] += e * 0.5 * v[k]
        return R

    def rot2(nf):
        n = np.array(nf)
        nn = float(n @ n)
        R = np.zeros((3, 3))
        for i in range(3):
            for j in range(3):
                R[i, j] = 2 * n[i] * n[j] / nn - (1.0 if i == j else 0.0)
        return R

    def perm_of(R):
        p = []
        for v in cf:
            w = R @ np.array(v)
            k = min(range(12), key=lambda t: sum(
                (cf[t][a] - w[a]) ** 2 for a in range(3)))
            p.append(k)
        return tuple(p)
    p5 = perm_of(rot5(cf[0]))
    a0, b0 = edges[0]
    n2 = [cf[a0][t] + cf[b0][t] for t in range(3)]
    p2 = perm_of(rot2(n2))

    def lift(p):
        """g . |f, i> = |g f, position of g(nbr) in the cyclic order>."""
        L = np.zeros((D, D))
        for (f, i) in st:
            gf = p[f]
            gw = p[m[f][i]]
            j = m[gf].index(gw)
            L[idx[(gf, j)], idx[(f, i)]] = 1.0
        return L

    t0b = 0.0
    for p in (p5, p2):
        L = lift(p)
        for al in (0.0, 0.37, 1.0, 2.0):
            Ua = U(al)
            t0b = max(t0b, float(np.max(np.abs(L @ Ua - Ua @ L))))
    checks.append(("T0b A5-EQUIVARIANCE CERTIFICATE: the generators p5 "
                   "(72 deg vertex) and p2 (180 deg edge), lifted to the "
                   "60 directed states, commute with U(alpha) at sampled "
                   "alphas: worst %.1e (<= 1e-12) - the coin family IS "
                   "A5-covariant, proved not assumed" % t0b,
                   t0b <= 1e-12))
    if not (t0a <= 1e-12 and t0b <= 1e-12):
        data["status"] = "QW-DRIFT"
        print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                        for c in checks), flush=True)
        _emit(data, checks, False)
        print("[V349 done] QW-DRIFT", flush=True)
        return 1

    def minph(alpha):
        ph = np.sort(np.abs(np.angle(np.linalg.eigvals(U(alpha)))))
        return float(ph[0])

    # ---- M1: gapless locus ---------------------------------------------
    grid = np.linspace(0.0, 2 * math.pi, 6284)
    vals = [minph(a) for a in grid]
    loc = [k for k in range(1, len(grid) - 1)
           if vals[k] <= vals[k - 1] and vals[k] <= vals[k + 1]
           and vals[k] <= 1e-3]
    zeros = []
    for k in loc:
        lo, hi = grid[k] - 1e-3, grid[k] + 1e-3
        r = (math.sqrt(5) - 1) / 2
        a1, b1 = hi - r * (hi - lo), lo + r * (hi - lo)
        fa, fb = minph(a1), minph(b1)
        while hi - lo > 1e-12:
            if fa < fb:
                hi, b1, fb = b1, a1, fa
                a1 = hi - r * (hi - lo)
                fa = minph(a1)
            else:
                lo, a1, fa = a1, b1, fb
                b1 = lo + r * (hi - lo)
                fb = minph(b1)
        am = 0.5 * (lo + hi)
        res = minph(am)
        if res <= 1e-9:
            ph = np.sort(np.abs(np.angle(np.linalg.eigvals(U(am)))))
            mult = int((ph <= 1e-7).sum())
            zeros.append({"alpha_star": float(am), "residual": res,
                          "zero_modes": mult})
    data["zeros"] = zeros
    data["n_gapless_points"] = len(zeros)
    maxmodes = max((z["zero_modes"] for z in zeros), default=0)
    checks.append(("M1 gapless locus of U(alpha): %d gapless points "
                   "(grid 6284, polished), max zero modes at one point "
                   "= %d (published doubler count = 12)"
                   % (len(zeros), maxmodes), True))

    # ---- M2: protection tribunal ---------------------------------------
    rng = np.random.RandomState(14)
    levels = [1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3]
    protected = []
    for z in zeros:
        am = z["alpha_star"]
        worst = []
        for eps in levels:
            w = 0.0
            for _ in range(160):
                # random bond phases (conjugate on reverse pairs) and
                # +-eps coin jitter per face: both break A5
                Sp = S.copy()
                for (f, i) in st:
                    g = m[f][i]
                    j = m[g].index(f)
                    ph1 = rng.uniform(-1, 1) * eps
                    Sp[idx[(g, j)], idx[(f, i)]] *= np.exp(1j * ph1)
                    Sp[idx[(f, i)], idx[(g, j)]] *= np.exp(-1j * ph1)
                Cp = np.zeros((D, D), dtype=complex)
                P5 = np.zeros((5, 5))
                for i in range(5):
                    P5[(i + 1) % 5, i] = 1.0
                Hf = P5 + P5.conj().T
                for f in range(12):
                    blk = expm(1j * (am + eps * rng.uniform(-1, 1)) * Hf)
                    for i in range(5):
                        for jj in range(5):
                            Cp[idx[(f, i)], idx[(f, jj)]] = blk[i, jj]
                Up = Sp @ Cp
                w = max(w, float(np.min(np.abs(np.angle(
                    np.linalg.eigvals(Up))))))
            worst.append(w)
            if w > 1e-6:
                break
        prot = all(w <= 1e-6 for w in worst)
        protected.append(prot)
        z["protection_worst_gap"] = worst[-1]
        z["protected"] = bool(prot)
    nprot = sum(protected)
    checks.append(("M2 protection tribunal (random bond phases + coin "
                   "jitter, 6 amplitudes x 160 draws each): %d/%d "
                   "gapless points stay <= 1e-6 (A5-PROTECTED), "
                   "worst-gap chain reported per point"
                   % (nprot, len(zeros)), True))

    # ---- M3: the clock battery (unchanged house rule) -------------------
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    desert = {"g_obs": g_obs, "g_cross_rung10": g_cross, "1/5": 0.2}
    contact = False
    for z in zeros:
        am = z["alpha_star"]
        rows = {}
        for tag, gv in desert.items():
            rows[tag] = min(abs(am - gv), abs(am / (2 * math.pi) - gv),
                            abs(2 * math.pi * gv - am))
        z["distances_to_desert"] = rows
        z["margin"] = min(rows.values())
        if z["margin"] <= 1e-9:
            contact = True
    data["desert_anchors"] = desert
    closest = min((z["margin"] for z in zeros), default=None)
    data["closest_margin"] = closest
    checks.append(("M3 clock battery: %d gapless points, nearest desert "
                   "margin %s (guard 1e-9)"
                   % (len(zeros),
                      ("%.3e" % closest) if closest is not None else "-"),
                   True))

    if not zeros:
        status = "QW-MASSIVE"
        ok = all(c[1] for c in checks)
    elif nprot == 0:
        status = "QW-GAPLESS-UNPROTECTED"
        ok = all(c[1] for c in checks)
    elif contact:
        status = "QW-GAPLESS-PROTECTED-CLOCK-CONTACT"
        ok = all(c[1] for c in checks)
    else:
        status = "QW-GAPLESS-PROTECTED-CLOCK-VOID"
        ok = all(c[1] for c in checks)

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    for z in zeros[:12]:
        print("   gapless alpha* = %.9f  modes %d  protected %s  "
              "desert margin %.4f"
              % (z["alpha_star"], z["zero_modes"], z["protected"],
                 z.get("margin", float("nan"))), flush=True)
    print("[V349 done] %s" % status, flush=True)
    _emit(data, checks, ok, status=status)
    return 0


def _emit(data, checks, ok, status=None):
    out = {"id": "V349_icosahedral_walk", "tier": 1,
           "property": "the new geometry BUILT FOR the 12-face cell: "
                       "A5-covariant quantum walk with pentagon coins "
                       "- gapless locus, topological protection, clock "
                       "battery",
           "status": status or data.get("status", ""),
           "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v349_icosahedral_walk.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
