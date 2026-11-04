#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V341 - N2b: the icosahedral cell dressed by the theory's OWN internal
Spin(10) module. V339's cell alone carries one gapless point at GAMMA
with QUADRATIC dispersion (z = 2); the published engine's nodes are
LINEAR because of the internal layer. supercell.bond_product says the
step factorizes EXACTLY as kron(P(q;beta), U(g)**n_stages) - so the
dressed spectrum is cell phases + internal phases, no new gate design
invented: the internal layer is the theory's own, one application per
stage (five for the icosahedral cycle)."""
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
from cl10 import supercell as SC  # noqa: E402


def wrap(x):
    return (x + math.pi) % (2 * math.pi) - math.pi


class F5:
    """a + b*sqrt(5) over Fractions."""
    __slots__ = ("a", "b")

    def __init__(self, a, b=Fraction(0)):
        self.a, self.b = Fraction(a), Fraction(b)

    def __add__(s, o):
        return F5(s.a + o.a, s.b + o.b)

    def __mul__(s, o):
        return F5(s.a * o.a + 5 * s.b * o.b, s.a * o.b + s.b * o.a)

    def __eq__(s, o):
        return s.a == o.a and s.b == o.b

    def val(s):
        return float(s.a) + float(s.b) * math.sqrt(5)


def cell_build():
    PHI = F5(Fraction(1, 2), Fraction(1, 2))
    ZERO = F5(0)
    ONE = F5(1)
    coords = []
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        bp = F5(PHI.a * b, PHI.b * b)
        coords.append((ZERO, F5(a), bp))
        coords.append((bp, ZERO, F5(a)))
        coords.append((F5(a), bp, ZERO))

    def dot(u, v):
        s = ZERO
        for i in range(3):
            s = s + u[i] * v[i]
        return s
    edges = [(i, j) for i in range(12) for j in range(i + 1, 12)
             if dot(coords[i], coords[j]) == PHI]

    def matchings_gen(es):
        def rec(chosen, avail):
            used = {v for e in chosen for v in e}
            if len(used) == 12:
                yield chosen
                return
            v = min(set(range(12)) - used)
            for e in [x for x in avail if v in x]:
                na = [x for x in avail if x != e
                      and not (set(x) & set(e))]
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
    edge_dir = {}
    for (i, j) in edges:
        d = [coords[j][t].val() - coords[i][t].val() for t in range(3)]
        nz = [t for t in range(3) if abs(d[t]) > 1e-9]
        if len(nz) == 1:
            edge_dir[(i, j)] = (nz[0], 1 if d[nz[0]] > 0 else -1)
    return edges, fac, edge_dir


def main():
    checks = []
    edges, fac, edge_dir = cell_build()
    checks.append(("T0c icosahedral cell rebuilt (V339 machinery re-run "
                   "inside): %d edges, 1-factorization %d matchings, %d "
                   "axis-parallel edges"
                   % (len(edges), len(fac), len(edge_dir)),
                   len(edges) == 30 and len(fac) == 5
                   and len(edge_dir) == 6))

    # ---- T0a: published factorization anchor through this lens ---------
    B, U32, U16 = SC.sectors()
    q = (0.31, -0.77, 0.44)
    full = SC.bloch(q, U16)
    P = SC.bond_product(q)
    kronv = np.kron(P, np.linalg.matrix_power(U16, 6))
    ferr = float(np.max(np.abs(full - kronv)))
    gM = float(SC.gap(SC.M_POINTS[0], U16))
    checks.append(("T0a published engine anchor: bloch == kron(P, U^6) to "
                   "%.1e; published M-point gaplessness %.1e"
                   % (ferr, gM), ferr < 1e-10 and gM < 1e-9))

    # ---- the dressed N2b spectrum --------------------------------------
    def cell_U(qv, g, stage_pow=1):
        u = np.eye(12, dtype=complex)
        for me in fac:
            M = np.zeros((12, 12), dtype=complex)
            for (i, j) in me:
                if (i, j) in edge_dir:
                    a, sg = edge_dir[(i, j)]
                    ph = np.exp(1j * sg * qv[a])
                    M[i, j] = ph
                    M[j, i] = np.conj(ph)
                else:
                    M[i, j] = 1.0
                    M[j, i] = 1.0
            u = (np.cos(g) * np.eye(12, dtype=complex)
                 - 1j * np.sin(g) * M) @ u
        return u

    clocks = {"2pi/5": 2 * math.pi / 5, "4pi/5": 4 * math.pi / 5,
              "pi/4": math.pi / 4}
    results = {}
    for cname, g in clocks.items():
        Uint = np.linalg.matrix_power(U16, 5)          # one layer per stage
        bph = np.sort(np.angle(np.linalg.eigvals(Uint)))
        # cell-only anchor at GAMMA
        cg = min(abs(np.angle(np.linalg.eigvals(cell_U((0., 0., 0.), g)))))
        # dressed grid scan
        N = 24
        best = (1e9, None)
        zeros = 0
        for ai in range(N):
            qx = -math.pi + 2 * math.pi * ai / N
            for bi in range(N):
                qy = -math.pi + 2 * math.pi * bi / N
                for ci in range(N):
                    qz = -math.pi + 2 * math.pi * ci / N
                    aph = np.angle(np.linalg.eigvals(cell_U((qx, qy, qz), g)))
                    gp = min(abs(wrap(p + b)) for p in aph for b in bph)
                    if gp < best[0]:
                        best = (gp, (qx, qy, qz))
                    if gp <= 1e-7:
                        zeros += 1
        # refine + z exponent along 5 rays
        cur = list(best[1])
        curgap = best[0]
        st = 0.2
        while st > 1e-7:
            improved = False
            for d in range(3):
                for s in (+1, -1):
                    t = cur[:]
                    t[d] = wrap(t[d] + s * st)
                    aph = np.angle(np.linalg.eigvals(cell_U(tuple(t), g)))
                    gp = min(abs(wrap(p + b)) for p in aph for b in bph)
                    if gp < curgap:
                        cur, curgap, improved = t, gp, True
            if not improved:
                st /= 2
        rays = {}
        for (nm, dv) in (("x", (1, 0, 0)), ("y", (0, 1, 0)), ("z", (0, 0, 1)),
                         ("golden", (1, (1 + 5 ** 0.5) / 2, 0)),
                         ("body", (1, 1, 1))):
            nr = math.sqrt(sum(d * d for d in dv))
            vals = []
            for eps in (1e-3, 3e-2):
                qv = tuple(wrap(c + eps * d / nr)
                           for c, d in zip(cur, dv))
                aph = np.angle(np.linalg.eigvals(cell_U(qv, g)))
                vals.append(min(abs(wrap(p + b))
                                for p in aph for b in bph))
            rays[nm] = round(math.log((vals[1] + 1e-300)
                                      / (vals[0] + 1e-300))
                             / math.log(30.0), 3)
        results[cname] = {"cell_only_gap_at_GAMMA": cg,
                          "dressed_grid_min": best[0],
                          "dressed_at": [round(x, 4) for x in cur],
                          "dressed_refined": curgap,
                          "grid_zeros": zeros, "z_per_ray": rays}
    g25 = results["2pi/5"]
    linear = any(0.8 <= z <= 1.2 for z in g25["z_per_ray"].values()) \
        and g25["dressed_refined"] <= 1e-9 and g25["grid_zeros"] > 0
    quad = all(z > 1.5 for z in g25["z_per_ray"].values())
    gapless = g25["dressed_refined"] <= 1e-9 and g25["grid_zeros"] > 0
    checks.append(("T1 dressed N2b at the five-clocks and pi/4: %s"
                   % json.dumps(results, default=str), True))
    checks.append(("T2 cell-alone GAMMA anchor at 2pi/5 reproduced: %.1e"
                   % g25["cell_only_gap_at_GAMMA"],
                   g25["cell_only_gap_at_GAMMA"] < 1e-9))

    allgreen = all(c[1] for c in checks)
    status = ("N2B-LINEAR-CHIRAL-NODES" if allgreen and linear else
              "N2B-QUADRATIC-TOUCHINGS" if allgreen and gapless and quad
              else "N2B-GAPPED" if allgreen and not gapless else
              "N2B-MIXED" if allgreen else "N2B-MACHINERY-DRIFT")
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V341 done] %s" % status, flush=True)
    out = {"id": "V341_icosahedral_internal_dressing", "tier": 2,
           "property": "N2b: icosahedral cell x the theory's own internal "
                       "Spin(10) layer, one application per stage: is the "
                       "GAMMA touching lifted to a linear chiral node?",
           "status": status, "pass": bool(allgreen),
           "data": {"results": results,
                    "published_anchor": {"kron_err": ferr,
                                         "M_gap_published": gM},
                    "internal_exponent": 5,
                    "tribunal": "g_obs naming remains V338/V340-closed"},
           "checks": [{"name": c[0][:600], "pass": c[1]} for c in checks]}
    with open(os.path.join(HERE, "results",
                           "v341_icosahedral_internal_dressing.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
