#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V339 - candidate N2a: the icosahedral 12-site cell itself.

A5 has no transitive action on 8 sites (V337's booked reason the digit 5
cannot share the cubic cell) but acts transitively on the 12 vertices of
the icosahedron. This experiment BUILDS that cell and judges it: exact
Q(sqrt5)/Q(sqrt5, sqrt3, t) combinatorics (t^2 = 10 + 2 sqrt5), A5 as
exact vertex permutations, the 1-factorization of the 30-edge graph into
5 stage matchings, the clock-birth identity |V|/|A5| = 12/60 = 1/5, a
momentum automaton whose 12 axis-parallel edges carry Bloch phases, and
the fermion tribunal (are there exactly gapless momenta at the closure
clock 2 pi/5?).
"""
import json
import math
import os
import sys
from fractions import Fraction

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


class F3:
    """Element of Q(sqrt5, sqrt3, t), t^2 = 10 + 2 sqrt5, as coefficients
    on the basis (sqrt5)^i (sqrt3)^j t^k, i,j,k in {0,1}: Fractions."""
    __slots__ = ("c",)

    def __init__(self, c):
        self.c = tuple(c)

    @staticmethod
    def q(x):
        c = [Fraction(0)] * 8
        c[0] = Fraction(x)
        return F3(c)

    def idx(self, i, j, k):
        return i + 2 * j + 4 * k

    def __add__(self, o):
        return F3([a + b for a, b in zip(self.c, o.c)])

    __radd__ = __add__

    def __neg__(self):
        return F3([-a for a in self.c])

    def __sub__(self, o):
        return F3([a - b for a, b in zip(self.c, o.c)])

    def __mul__(self, o):
        c = [Fraction(0)] * 8
        for p in range(8):
            if self.c[p] == 0:
                continue
            i1, j1, k1 = p % 2, (p // 2) % 2, p // 4
            for m in range(8):
                if o.c[m] == 0:
                    continue
                i2, j2, k2 = m % 2, (m // 2) % 2, m // 4
                coef = self.c[p] * o.c[m]
                # t^2 = 10 + 2 sqrt5: expand (t^k1)(t^k2) first
                kk = k1 + k2
                tparts = [(Fraction(0), 0), (Fraction(1), 0)]
                if kk == 1:
                    tparts = [(Fraction(1), 1)]
                elif kk == 2:
                    tparts = [(Fraction(10), 0), (Fraction(2), 2)]  # (coef, i-add)
                for (tc, extra_i) in tparts:
                    ii = (i1 + i2 + extra_i) % 2
                    carry5 = ((i1 + i2 + extra_i) // 2) * 5
                    jj = (j1 + j2) % 2
                    carry3 = ((j1 + j2) // 2) * 3
                    # carry constants multiply coef (both <= 5, 3)
                    base = coef * tc * carry5 * carry3
                    c[F3(0).idx(ii, jj, kk % 2 if kk < 2 else 0)] += base \
                        if kk < 2 else base  # handled below explicitly
                # (done per-part above; skip fallback)
        return F3(c)

    def __repr__(self):
        return repr(self.c)


def mul(a, b):
    """Correct tower multiplication (basis 8, reduction sqrt5^2=5,
    sqrt3^2=3, t^2=10+2 sqrt5)."""
    out = [Fraction(0)] * 8

    def basis_mult(p, m):
        i1, j1, k1 = p % 2, (p // 2) % 2, p // 4
        i2, j2, k2 = m % 2, (m // 2) % 2, m // 4
        # t part
        ksum = k1 + k2
        if ksum == 0:
            tp = [(Fraction(1), 0, 0)]  # coef, i-extra, k-result
        elif ksum == 1:
            tp = [(Fraction(1), 0, 1)]
        else:
            tp = [(Fraction(10), 0, 0), (Fraction(2), 1, 0)]
        res = []
        for (tc, ei, kr) in tp:
            isum = i1 + i2 + ei
            cf = tc * (5 if isum >= 2 else 1)
            ii = isum % 2
            jsum = j1 + j2
            cf *= (3 if jsum >= 2 else 1)
            jj = jsum % 2
            res.append((cf, ii, jj, kr))
        return res

    for p in range(8):
        if a[p] == 0:
            continue
        for m in range(8):
            if b[m] == 0:
                continue
            coef = a[p] * b[m]
            for (cf, ii, jj, kr) in basis_mult(p, m):
                out[ii + 2 * jj + 4 * kr] += coef * cf
    return out


def Fadd(x, y):
    return F3([p + q for p, q in zip(x.c, y.c)])


def Fmul(x, y):
    return F3(mul(x.c, y.c))


def Fsub(x, y):
    return F3([p - q for p, q in zip(x.c, y.c)])


def Finv(x):
    """Solve x * y = 1 by linear algebra on the multiplication matrix."""
    A = []
    cols = [mul(x.c, [Fraction(1) if b == c else Fraction(0)
                      for c in range(8)]) for b in range(8)]
    for r in range(8):
        A.append([cols[b][r] for b in range(8)])
    # A is the matrix of mult by x in the basis; invert by row reduction
    M = [row[:] + [Fraction(1) if r == c else Fraction(0)
                   for c in range(8)] for r, row in enumerate(A)]
    for col in range(8):
        piv = next(r for r in range(col, 8) if M[r][col] != 0)
        M[col], M[piv] = M[piv], M[col]
        pv = M[col][col]
        M[col] = [v / pv for v in M[col]]
        for r in range(8):
            if r != col and M[r][col] != 0:
                f = M[r][col]
                M[r] = [a - f * b for a, b in zip(M[r], M[col])]
    y = [M[r][8] for r in range(8)]
    # y is the inverse's coefficients? mult-by-x applied row-wise; verify
    z = mul(x.c, y)
    assert z == [Fraction(1)] + [Fraction(0)] * 7, "inverse check"
    return F3(y)


def Fq(x):
    return F3.q(x)


SQ5 = F3([Fraction(0)] * 8)
_c = [Fraction(0)] * 8
_c[1] = Fraction(1)
SQ5 = F3(_c)
_c = [Fraction(0)] * 8
_c[2] = Fraction(1)
W3 = F3(_c)
_c = [Fraction(0)] * 8
_c[4] = Fraction(1)
TT = F3(_c)
PHI = Fadd(Fq(1), SQ5)
PHI = F3([x / Fraction(2) for x in PHI.c])  # (1+sqrt5)/2


def main():
    checks = []

    # ---- T1: the icosahedron exactly ----------------------------------
    # cyclic permutations of (0, +-1, +-phi)
    coords = []
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        bp = F3([x * b for x in PHI.c])
        coords.append((Fq(0), Fq(a), bp))        # (0, a, b phi)
        coords.append((bp, Fq(0), Fq(a)))        # cyclic: (b phi, 0, a)
        coords.append((Fq(a), bp, Fq(0)))        # cyclic: (a, b phi, 0)
    assert len(coords) == 12 and len(set(tuple(tuple(f.c) for f in v)
                                          for v in coords)) == 12
    def dot(u, v):
        s = Fq(0)
        for i in range(3):
            s = Fadd(s, Fmul(u[i], v[i]))
        return s
    norm2 = dot(coords[0], coords[0])
    edges = []
    for i in range(12):
        for j in range(i + 1, 12):
            if dot(coords[i], coords[j]).c == PHI.c:
                edges.append((i, j))
    deg = [sum(1 for (a, b) in edges if i in (a, b)) for i in range(12)]
    checks.append(("T1 icosahedron EXACT: norm^2 = phi+2 (%s), %d edges "
                   "(booked 30), degree set {%s} (booked {5})"
                   % ("OK" if norm2.c[0] == Fraction(5, 2) and
                      norm2.c[1] == Fraction(1, 2) else "BAD",
                      len(edges), sorted(set(deg))),
                   len(edges) == 30 and set(deg) == {5}))

    # ---- T2: A5 exactly as vertex permutations ------------------------
    def eps(i, j, k):
        if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
            return 1
        if (i, j, k) in ((0, 2, 1), (2, 1, 0), (1, 0, 2)):
            return -1
        return 0

    cos72 = F3([x / Fraction(4) for x in Fsub(SQ5, Fq(1)).c])

    def rot5(v):
        """72 deg about the VERTEX axis v (|v|^2 = phi+2). The crucial
        exact identity: sin 72 = |v| / 2, so sin72 * u_k = v_k / 2 and
        every entry lives in Q(sqrt5, t) with NO new radicals."""
        assert (Fmul(v, v) if False else dot(v, v)).c == Fadd(
            PHI, Fq(2)).c
        inv2 = Finv(dot(v, v))
        R = [[Fq(0) for _ in range(3)] for _ in range(3)]
        half = F3([Fraction(1, 2)] + [Fraction(0)] * 7)
        for i in range(3):
            for j in range(3):
                uu = Fmul(Fmul(v[i], v[j]), inv2)
                R[i][j] = Fadd(R[i][j],
                               Fmul(cos72, Fq(1)) if i == j else Fq(0))
                R[i][j] = Fadd(R[i][j], Fmul(Fsub(Fq(1), cos72), uu))
                for k in range(3):
                    e = eps(i, k, j)
                    if e:
                        R[i][j] = Fadd(
                            R[i][j],
                            Fmul(Fq(e) if e == 1 else Fsub(Fq(0), Fq(1)),
                                 Fmul(half, v[k])))
        return R

    def rot2(n):
        """180 deg about the edge-midpoint axis n: entries in Q(sqrt5)."""
        inv2 = Finv(dot(n, n))
        R = [[Fq(0) for _ in range(3)] for _ in range(3)]
        for i in range(3):
            for j in range(3):
                R[i][j] = Fadd(Fmul(Fq(2), Fmul(Fmul(n[i], n[j]), inv2)),
                               Fsub(Fq(0), Fq(1)) if i == j else Fq(0))
        return R

    def perm_of(R):
        p = []
        for v in coords:
            w = [sum((Fmul(R[i][j], v[j]) for j in range(3)), Fq(0))
                 for i in range(3)]
            k = next(t for t in range(12)
                     if all(coords[t][i].c == w[i].c for i in range(3)))
            p.append(k)
        return p

    p5 = tuple(perm_of(rot5(coords[0])))
    a, b = edges[0]
    n2axis = [Fadd(coords[a][t], coords[b][t]) for t in range(3)]
    p2 = tuple(perm_of(rot2(n2axis)))
    p3 = tuple(p5[p2[i]] for i in range(12))

    def perm_order(p):
        seen = [False] * 12
        lc = 1
        from math import gcd
        for i in range(12):
            if not seen[i]:
                L = 0
                j = i
                while not seen[j]:
                    seen[j] = True
                    j = p[j]
                    L += 1
                lc = lc * L // gcd(lc, L)
        return lc
    G = {tuple(range(12))}
    fr = [tuple(range(12))]
    while fr:
        x = fr.pop()
        for g in (p5, p2):
            y = tuple(x[g[i]] for i in range(12))
            if y not in G:
                G.add(y)
                fr.append(y)
    orbit = set()
    x = tuple(range(12))
    # orbit of vertex 0 under G (as images: track where 0 goes)
    for g in G:
        orbit.add(g[0])
    # p5 acts as permutation: apply to coords: image of vertex v is p[v]
    checks.append(("T2 A5 EXACT: R5 (72 deg vertex) order %d, R2 (180 deg "
                   "edge) order %d, product R5R2 order %d (the (2,3,5) "
                   "presentation), group %d (booked 60), orbit of a vertex "
                   "%d (transitive), stabilizer %d (booked 5)"
                   % (perm_order(p5), perm_order(p2), perm_order(p3),
                      len(G), len(orbit), len(G) // len(orbit)),
                   perm_order(p5) == 5 and perm_order(p2) == 2 and
                   perm_order(p3) == 3 and len(G) == 60 and
                   len(orbit) == 12))

    # ---- T3: clock-birth identities + 1-factorization ------------------
    birth = (Fraction(12) / Fraction(60)) == Fraction(1, 5)

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
    # stage order certificate: S(g)^k = cos kg - i sin kg P: at g = 2pi/5,
    # S^5 = I exactly
    th = 2 * math.pi / 5
    S5check = abs(math.cos(5 * th) - 1) + abs(math.sin(5 * th))
    checks.append(("T3 clock birth: |V|/|A5| = 12/60 = 1/5 EXACT: %s; "
                   "1-factorization into 5 perfect matchings of 6 edges: "
                   "%s; stage S(2pi/5)^5 = I to %.1e (closure clock 2pi/5, "
                   "Coxeter h(H3) = 10)"
                   % (birth, fac_ok, S5check), birth and fac_ok and
                   S5check < 1e-12))

    # ---- T4: the automaton ---------------------------------------------
    idx = {v: k for k, v in enumerate(coords)}
    edge_dir = {}
    for (i, j) in edges:
        d = [Fsub(coords[j][t], coords[i][t]) for t in range(3)]
        nz = [t for t in range(3) if any(x != 0 for x in d[t].c)]
        if len(nz) == 1:
            a = nz[0]
            fl = [float(x) for x in d[a].c]
            sg = 1 if fl[0] > 0 else -1
            edge_dir[(i, j)] = (a, sg)
    per_axis = [sum(1 for e in edge_dir if edge_dir[e][0] == a)
                for a in range(3)]
    # (the golden-rectangle long sides have length 2 phi and are NOT
    # icosahedral edges: 6 axis-parallel edges = 2 per axis is the right
    # geometry; the booking of 4 per axis was wrong and is corrected)
    g = 2 * math.pi / 5

    def stage_m(qv, me):
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
        return np.cos(g) * np.eye(12, dtype=complex) - 1j * np.sin(g) * M

    def U(qv):
        u = np.eye(12, dtype=complex)
        for me in fac:
            u = stage_m(qv, me) @ u
        return u

    rng = np.random.RandomState(7)
    uni = 0.0
    for _ in range(40):
        qv = tuple(rng.uniform(-math.pi, math.pi, 3))
        uq = U(qv)
        uni = max(uni, float(np.max(np.abs(uq.conj().T @ uq
                                           - np.eye(12)))))
    # (the booked q -> -q mirror identity is FALSE for this automaton:
    # conj flips -i sin g M to +i sin g M*, so S_m(-q) = S_m(q)^T and
    # U(-q) is the TRANSPOSED-ORDER product - a five-stage ordered cycle
    # need not be isospectral to its reverse. Measured, not gated.)
    mir = 0.0
    for _ in range(40):
        qv = tuple(rng.uniform(-math.pi, math.pi, 3))
        a = np.sort(np.angle(np.linalg.eigvals(U(qv))))
        b = np.sort(np.angle(np.linalg.eigvals(U(tuple(-x for x in qv)))))
        mir = max(mir, float(np.max(np.abs(a - b))))
    checks.append(("T4 automaton: axis-parallel edges %s (6 total, 2 per "
                   "axis: only the golden-rectangle short sides are edges "
                   "- booked 4 per axis corrected before the run), stage "
                   "unitarity exact (max U^H U deviation over 40 random "
                   "q: %.1e); q -> -q spectral mirror residual %.1e -> "
                   "ABSENT (the ordered five-stage icosahedral cycle is "
                   "chiral by composition; the on-site +- mirror was a "
                   "cubic-cell property) - measured data, not a gate"
                   % (per_axis, uni, mir),
                   per_axis == [2, 2, 2] and uni < 1e-13))

    # ---- T5: fermion tribunal -----------------------------------------
    best = (1e9, None)
    zeros = 0
    N = 32
    for a in range(N):
        qx = -math.pi + 2 * math.pi * a / N
        for b in range(N):
            qy = -math.pi + 2 * math.pi * b / N
            for c in range(N):
                qz = -math.pi + 2 * math.pi * c / N
                qv = (qx, qy, qz)
                gap = min(abs(np.angle(np.linalg.eigvals(U(qv)))))
                if gap < best[0]:
                    best = (gap, qv)
                if gap <= 1e-8:
                    zeros += 1
    # refine: is the best point a true node? (local descent)
    (gapb, qb) = best
    cur = list(qb)
    curgap = gapb
    st = 0.2
    while st > 1e-6:
        improved = False
        for i in range(3):
            for s in (+1, -1):
                trial = cur[:]
                trial[i] = (trial[i] + s * st + math.pi) % (2 * math.pi) - math.pi
                gp = min(abs(np.angle(np.linalg.eigvals(U(tuple(trial))))))
                if gp < curgap:
                    cur, curgap, improved = trial, gp, True
        if not improved:
            st /= 2
    gapfam = {}
    for gg in (math.pi / 5, 2 * math.pi / 5, 4 * math.pi / 5, math.pi / 4):
        g = gg
        mb = 1e9
        for a in range(17):
            for b in range(17):
                for c in range(17):
                    qv = (-math.pi + 2 * math.pi * a / 16,
                          -math.pi + 2 * math.pi * b / 16,
                          -math.pi + 2 * math.pi * c / 16)
                    mb = min(mb, min(abs(np.angle(
                        np.linalg.eigvals(U(qv))))))
        gapfam["%.4f pi" % (gg / math.pi)] = mb
    # spectral structure probes at the closure clock
    g = 2 * math.pi / 5
    probes = {}
    for (qname, qv) in (("GAMMA", (0.0, 0.0, 0.0)),
                        ("random1", (0.3, -1.2, 2.0)),
                        ("random2", (-2.9, 0.7, -0.4))):
        uq = U(qv)
        an = np.sort(np.angle(np.linalg.eigvals(uq)))
        err5 = float(np.max(np.abs(np.linalg.matrix_power(uq, 5)
                                   - np.eye(12))))
        n1 = int(np.sum(np.abs(an) < 1e-9))
        probes[qname] = {"angles": [round(float(x), 6) for x in an],
                         "U5_minus_I": err5, "n_eigenvalue_1": n1}
    # dispersion of the nontrivial sector: second-smallest |phase| spread
    secs = []
    for _ in range(200):
        qv = tuple(rng.uniform(-math.pi, math.pi, 3))
        an = np.sort(np.abs(np.angle(np.linalg.eigvals(U(qv)))))
        nz2 = [x for x in an if x > 1e-9]
        secs.append(min(nz2) if nz2 else math.pi)
    flat = zeros > 512
    gapless = (0 < zeros <= 512) and not flat
    # Dirac dispersion at the isolated node: gap / |q| along several rays
    slopes = {}
    for (nm, dv) in (("x", (1.0, 0.0, 0.0)), ("y", (0.0, 1.0, 0.0)),
                     ("z", (0.0, 0.0, 1.0)),
                     ("golden", (1.0, float(PHI.c[0] + PHI.c[1] *
                                            math.sqrt(5)), 0.0)),
                     ("body", (1.0, 1.0, 1.0))):
        nrm = math.sqrt(sum(d * d for d in dv))
        vals = []
        for eps in (1e-3, 3e-3, 1e-2, 3e-2):
            qv = tuple(eps * d / nrm for d in dv)
            gp = min(abs(np.angle(np.linalg.eigvals(U(qv)))))
            vals.append(gp)
        z = math.log(vals[3] / vals[0]) / math.log(30.0)  # gap ~ |q|^z
        slopes[nm] = {"gap_over_eps": [round(v / e, 5) for v, e in
                                       zip(vals, (1e-3, 3e-3, 1e-2, 3e-2))],
                      "z": round(z, 3)}
    checks.append(("T5 fermions at the closure clock g = 2pi/5: grid min "
                   "gap %.3e at %s, refined %.3e; grid zeros <= 1e-8: %d "
                   "of %d (flat trivial-zero sector if many); nontrivial "
                   "sector min phase over 200 random q: %.3e (dispersed: "
                   "min %.3e max %.3e); Dirac slopes gap/|q| at the node "
                   "%s; probes %s; min gap by clock %s"
                   % (gapb, ["%.2f" % x for x in qb], curgap, zeros,
                      N ** 3, min(secs), min(secs), max(secs), slopes,
                      {k: (v["angles"], v["U5_minus_I"],
                           v["n_eigenvalue_1"])
                       for k, v in probes.items()},
                      {k: "%.3e" % v for k, v in gapfam.items()}), True))

    # ---- T6: closure order at GAMMA -------------------------------------
    U0 = U((0.0, 0.0, 0.0))
    kord = None
    p = np.eye(12, dtype=complex)
    for k in range(1, 121):
        p = p @ U0
        if np.max(np.abs(p - np.eye(12))) < 1e-9 or \
                np.max(np.abs(p + np.eye(12))) < 1e-9:
            kord = k
            break
    checks.append(("T6 closure order of U(GAMMA) at g = 2pi/5: %s (least "
                   "k <= 120 with U^k = +-I; None = the 5-stage product "
                   "has infinite order, measured data not a gate)" % kord,
                   True))

    allgreen = all(c[1] for c in checks)
    status = ("N2A-CELL-BUILT-GAPLESS-NODES" if allgreen and gapless else
              "N2A-CELL-BUILT-FLAT-TRIVIAL-ZERO" if allgreen and flat else
              "N2A-CELL-BUILT-NO-GAPLESS-NODES" if allgreen else
              "N2A-MACHINERY-DRIFT")
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V339 done] %s" % status, flush=True)
    out = {"id": "V339_icosahedral_cell", "tier": 2,
           "property": "candidate N2a: the icosahedral 12-site cell "
                       "(A5 transitive, closed 5-stage cycle, clock birth "
                       "12/60 = 1/5, fermion tribunal)",
           "status": status, "pass": bool(allgreen),
           "data": {"birth_12_60": "1/5 exact",
                    "edges": len(edges), "factorization": fac_ok,
                    "group_order": len(G),
                    "gap_clock_2pi5_refined": curgap,
                    "gap_at": [float(x) for x in cur],
                    "grid_zeros": zeros, "flat_zero_sector": flat,
                    "gapless_nodes": gapless,
                    "nontrivial_min_phase": min(secs),
                    "dirac_slopes": slopes,
                    "mirror_residual": mir,
                    "chiral_cycle": True,
                    "spectral_probes": probes,
                    "min_gap_by_clock": gapfam,
                    "closure_order_gamma": kord,
                    "tribunal_target_g_obs": 0.19999798835751165,
                    "closed_clock_radians": 2 * math.pi / 5,
                    "note": "the cell closes stage-by-stage at 2pi/5 "
                            "radians; the desert reads 0.1999980 radians; "
                            "12/60 = 1/5 names the BACKBONE exactly but "
                            "g_obs = 1/5 - 1.006e-5 remains un-named: the "
                            "deviation is the open theory"},
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    with open(os.path.join(HERE, "results",
                           "v339_icosahedral_cell.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
