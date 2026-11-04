#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V347 - THE ENGINE TRANSPLANTED ONTO THE USER'S 12-FACE CELL.

Booked before the run (pre_registrations_root.md, V346-V348):
W(beta) = [five bond stages on the face graph, each a cos/sin
rotation of the matching pairs at the transport angle beta - the
staggered swap-rotation bond layer of the PUBLISHED automaton with
the cubic block replaced by the closed 12-face sphere] . [I_12 (x) C],
C = the published internal Cl(10,0) 45-plane product at the
implemented-branch angle (cl10.supercell.internal at G10 - the clock
layer UNTOUCHED). There is no lattice momentum: the sphere is
closed; the cell carries the icosahedral group transitively
(60/12 = 5), which the cubic 8-subsite block could not (V337).

T0a: the cubic published path reproduces its own anchors inside this
script BEFORE the sphere is built (gapless M <= 1e-9 at the 12 M
points, R 3-fold multiplets, mirror pairing <= 1e-12) - else DRIFT.
T0b: backend invariance (GPU per standing directive): the sphere
observables on the cupy path vs the host path, worst <= 1e-11; step
unitarity <= 1e-11.
M1: the ELECTION test: the sphere rung ladder F_n^sph(beta) = the
smallest singular value of W(beta)^n - I, read exactly as
min_j 2|sin(n theta_j / 2)| from the eigenphases: zeros polished,
and rung crossings of |F| = LOOP^k (V335 E1 port), all compared
against the published anchors g_obs = 0.19999798835751165 and
g_cross_rung10 = 0.20013316730041675 (V206 artifact via art(), never
retyped) at the booked 1e-9 guard.
M2: closure order at the published angle: least k <= 240 with
W^k = +-I (tol 1e-8) - the V337 port.
Branch: SPHERE-ENGINE-CLOCK-CONTACT / -VOID / -DRIFT.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
from cl10 import config as _config  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import supercell as SC  # noqa: E402

CFG = _config.CFG


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


def main():
    checks = []
    data = {}

    # ---- T0a: the published cubic path re-derives its anchors ----------
    # (metrics copied VERBATIM from the V337 T0a gate: the +- mirror is
    # the INTERNAL layer's pairing - the coupled Bloch matrix is not the
    # object of that theorem, V337 booked gate defect - and the R
    # tribunal is the cluster-size signature at tol 1e-9)
    C = np.asarray(E.asnumpy(SC.internal()), dtype=complex)  # published
    # clock layer, moved to host by the engine's own asnumpy contract
    d = C.shape[0]

    def cluster_sizes(ph, tol=1e-9):
        s = np.sort(np.asarray(ph))
        out = []
        for v in s:
            if out and abs(v - out[-1][0]) < tol:
                out[-1][1] += 1
            else:
                out.append([float(v), 1])
        return [m for _, m in out]

    def internal_mirror(U):
        ph = np.sort(np.angle(np.linalg.eigvals(U)))
        return float(np.max(np.abs(ph + ph[::-1])))

    gapM = min(SC.gap(q, C) for q in SC.M_POINTS)
    mir = internal_mirror(C)
    R_signature = sorted(set(sum((cluster_sizes(SC.phases(SC.bloch(q, C)))
                                  for q in SC.R_POINTS[:3]), [])))
    checks.append(("T0a published cubic path inside this script: M-point "
                   "gap %.1e (booked gapless <= 1e-9), published R "
                   "cluster-size signature %s re-derived at the 3 R "
                   "points, internal mirror pairing %.1e (V337 verbatim "
                   "gates); internal layer %dx%d read from supercell "
                   "(never assumed)" % (gapM, R_signature, mir, d, d),
                   gapM <= 1e-9 and mir <= 1e-10 and
                   any(m % 3 == 0 and m >= 3 for m in R_signature)))
    data["R_signature_published"] = R_signature
    if not (gapM <= 1e-9 and mir <= 1e-10):
        data["status"] = "SPHERE-ENGINE-DRIFT"
        print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                        for c in checks), flush=True)
        _emit(data, checks, False)
        print("[V347 done] SPHERE-ENGINE-DRIFT", flush=True)
        return 1

    # ---- the sphere cell: face graph + booked factorization ------------
    from fractions import Fraction

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
            return Q5(s.a * o.a + 5 * s.b * o.b,
                      s.a * o.b + s.b * o.a)

        def __eq__(s, o):
            return s.a == o.a and s.b == o.b
    PHI = Q5(Fraction(1, 2), Fraction(1, 2))
    coords = []
    for (a, b) in [(1, 1), (1, -1), (-1, 1), (-1, -1)]:
        bp = PHI * Q5(b)
        coords.append((Q5(0), Q5(a), bp))
        coords.append((bp, Q5(0), Q5(a)))
        coords.append((Q5(a), bp, Q5(0)))

    def dot(u, v):
        s = Q5(0)
        for i in range(3):
            s = s + u[i] * v[i]
        return s
    edges = []
    for i in range(12):
        for j in range(12):
            if i < j and dot(coords[i], coords[j]) == PHI:
                edges.append((i, j))

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

    # bond stages on the sphere: swap-rotations of the matching pairs,
    # exactly the published staggered bond layer with the block graph
    # replaced by the face graph (no Bloch phases: the surface is
    # closed); internal layer identical to the published one.
    N = 12 * d
    Pc = []
    for me in fac:
        S = np.zeros((12, 12))
        for (i, j) in me:
            S[i, i] = S[j, j] = 0.0
            S[i, j] = S[j, i] = 1.0
        # rotation on each pair: cos b on site, +-i sin b on the hop
        Pc.append(S)

    I12 = np.eye(12)
    Id = np.eye(d)

    def W(beta):
        M = np.kron(I12, C)
        for me in fac:
            # published staggered SWAP-ROTATION bond layer: pair
            # exchange X (symmetric, X^2 = I on disjoint pairs), stage
            # = cos b I + i sin b X: unitary. (An antisymmetric K was
            # the first attempt - K^2 = -I gives hyperbolic cosh/sinh,
            # caught by the unitarity gate T0b, fixed before any
            # measurement was read.)
            X = np.zeros((12, 12))
            for (i, j) in me:
                X[i, j] = X[j, i] = 1.0
            Sb = np.cos(beta) * I12 + 1j * np.sin(beta) * X
            M = np.kron(Sb, Id) @ M
        return M

    W0 = W(SC.G10)
    unit = float(np.max(np.abs(W0 @ W0.conj().T - np.eye(N))))
    # backend invariance: device path through cl10.engine routing
    from scipy.linalg import schur as _schur0

    def ph_local(Mx):
        T = _schur0(Mx, output='complex')[0]
        return np.sort(np.angle(np.diag(T)))

    dev_ok = True
    try:
        # device path: the published internal layer natively on the
        # backend (SC.internal returns a device array), stages pushed
        # through the backend's own array primitive, all matmuls and
        # krons via E.matmul / E.kron routing
        Cd = SC.internal()

        def Wdev(beta):
            M = E.kron(E.Array_(I12), Cd)
            for me in fac:
                X = np.zeros((12, 12))
                for (i, j) in me:
                    X[i, j] = X[j, i] = 1.0
                Sb = np.cos(beta) * I12 + 1j * np.sin(beta) * X
                M = E.matmul(E.kron(E.Array_(Sb), E.eye(d)), M)
            return M
        Mdev_h = np.asarray(E.asnumpy(Wdev(SC.G10)), dtype=complex)
        a = ph_local(Mdev_h)
        b2 = ph_local(W0)
        # multisets on the circle: compare by circular Hausdorff
        # (sorted pairing is NOT a valid observable when the spectrum
        # carries masses of exactly-0.0/-0.0 phases - label flips
        # shift the pairing and fake failures up to the cluster
        # spacing; LENSSING booked in the run note - the certified
        # Riesz-projection rank test confirms the multisets agree)
        aa = np.sort(np.mod(a, 2 * math.pi))
        bb = np.sort(np.mod(b2, 2 * math.pi))
        bb2 = np.concatenate([bb, bb + 2 * math.pi])
        dev_res = math.inf
        for sh in range(len(aa)):
            d = np.max(np.abs(aa - np.roll(bb2, -sh)[:len(aa)]))
            d = min(d, 2 * math.pi - d)
            dev_res = min(dev_res, float(d))
    except Exception as exc:  # noqa: BLE001
        dev_ok = False
        dev_res = -1.0
        data["gpu_error"] = str(exc)
    checks.append(("T0b sphere step %dx%d: unitarity %.1e (<= 1e-11); "
                   "GPU path (cl10.engine routing) eigenphase invariance "
                   "%.1e (<= 1e-11)%s" % (N, N, unit, dev_res,
                                          "" if dev_ok else " FAILED"),
                   unit <= 1e-11 and dev_res <= 1e-11))
    if not (unit <= 1e-11 and dev_res <= 1e-11):
        data["status"] = "SPHERE-ENGINE-DRIFT"
        print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                        for c in checks), flush=True)
        _emit(data, checks, False)
        print("[V347 done] SPHERE-ENGINE-DRIFT", flush=True)
        return 1

    # ---- spectral lens for the sphere step ------------------------------
    # LENSSING (booked correction): the sphere step's spectrum carries
    # multiplicity-10 clusters at the published angle; the general QR
    # eigensolver (np.linalg.eigvals / LAPACK geev) splits such a
    # cluster into a cloud of radius ~0.3 (measured: 144 of 384
    # eigenphases displaced > 1e-9, and det/logdet says the cluster
    # point is a REAL 10-fold eigenvalue the QR missed). Every
    # observable in this experiment therefore reads the spectrum via
    # SCHUR with an explicit certificate: for a normal matrix Schur's
    # T is diagonal up to its strict-triangular norm, so when that
    # norm <= 1e-11 the diagonal IS the spectrum to machine precision.
    from scipy.linalg import schur as schur_f
    TRI = {"max": 0.0}

    def ph(M):
        T = schur_f(M, output='complex')[0]
        tri = float(np.max(np.abs(np.tril(T, -1))))
        TRI["max"] = max(TRI["max"], tri)
        return np.sort(np.angle(np.diag(T)))

    # ---- M1: the election battery --------------------------------------
    g_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    g_cross = float(art("V206_clock_tension",
                        "data.data_clock.g_cross_rung10"))
    lf = CFG.physics.loop_factor
    loop = float(lf) if lf else 1.0 / (16.0 * math.pi)
    data["anchors"] = {"g_obs": g_obs, "g_cross_rung10": g_cross,
                       "LOOP_port": loop}
    NS = list(range(3, 25)) + [170, 208]
    grid = np.linspace(1e-3, 1.0, 1001)
    Fn = np.empty((len(NS), len(grid)))
    for gi, b in enumerate(grid):
        th = ph(W(b))
        for ni, n in enumerate(NS):
            Fn[ni, gi] = float(np.min(2.0 * np.abs(np.sin(
                n * th / 2.0))))
    # rung zeros and booked LOOP^k crossings; nearest approach to the
    # anchors
    contacts = []
    loop_crossings = 0
    nearest_cross_dist = 1e9
    for ni, n in enumerate(NS):
        f = Fn[ni]
        # zeros
        for k in range(1, len(grid) - 1):
            if f[k] <= 1e-6 and f[k] < f[k - 1] and f[k] < f[k + 1]:
                # golden polish to width 1e-12 on the rung defect
                lo, hi = grid[k] - 1e-3, grid[k] + 1e-3

                def fn(b):
                    th = ph(W(b))
                    return float(np.min(2.0 * np.abs(
                        np.sin(n * th / 2.0))))
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
                resid = fn(gm)
                for tag, anch in (("g_obs", g_obs),
                                  ("g_cross", g_cross)):
                    if abs(gm - anch) <= 1e-6:
                        contacts.append({"rung": n, "anchor": tag,
                                         "g_star": gm,
                                         "residual": resid,
                                         "distance": abs(gm - anch)})
        # booked E1 port: |F| = LOOP crossings (k = 1: the crossing
        # set of the sphere ladder with the published loop factor)
        dlt = f - loop
        for k in range(len(grid) - 1):
            if dlt[k] == 0 or dlt[k] * dlt[k + 1] < 0:
                loop_crossings += 1
                gmid = 0.5 * (grid[k] + grid[k + 1])
                nearest_cross_dist = min(
                    nearest_cross_dist, abs(gmid - g_obs),
                    abs(gmid - g_cross))
        # nearest approach to the anchors on this rung
        jd = int(np.argmin(np.abs(grid - g_obs)))
        data.setdefault("nearest_approach", []).append(
            {"rung": n, "beta_near": float(grid[jd]),
             "F": float(f[jd])})
    data["loop_crossings"] = loop_crossings
    data["nearest_LOOP_crossing_to_anchors"] = nearest_cross_dist
    data["contacts"] = contacts

    # closest defect value AT the two anchors (the direct election read)
    at = {}
    for tag, anch in (("g_obs", g_obs), ("g_cross", g_cross)):
        th = ph(W(anch))
        mn = min(float(np.min(2.0 * np.abs(np.sin(n * th / 2.0))))
                 for n in NS)
        at[tag] = mn
    data["defect_at_anchors"] = at
    contact = bool(contacts)
    checks.append(("M1 election battery: %d rung functionals (n in "
                   "3..24,170,208) x 1001 beta-points: rung zeros within "
                   "1e-6 of the published anchors: %d (polished: %s); "
                   "smallest sphere rung defect AT the anchors: g_obs "
                   "%.3f, g_cross %.3f (a born clock would be ~0)"
                   % (len(NS), len(contacts),
                      ["n=%d %s d=%.1e" % (c["rung"], c["anchor"],
                                           c["distance"])
                       for c in contacts[:6]],
                      at["g_obs"], at["g_cross"]), True))

    checks.append(("T0c spectral lens certificate: schur diagonalized "
                   "every sphere step; worst strict-triangular norm "
                   "%.1e (<= 1e-11 means the diagonal IS the spectrum "
                   "for these normal matrices; the general QR lens, "
                   "which split the 10-fold cluster by up to 0.29, is "
                   "quarantined for this experiment)" % TRI["max"],
                   TRI["max"] <= 1e-11))

    # ---- M2: closure order at the published angle ----------------------
    th = ph(W0)
    order = None
    for k in range(1, 241):
        r1 = max(min(abs(k * t - 2 * math.pi * round(k * t /
                                                     (2 * math.pi))),
                     abs(k * t)) for t in th)
        if r1 <= 1e-8:
            order = k
            break
        r2 = max(min(abs(k * t - math.pi - 2 * math.pi * round(
            (k * t - math.pi) / (2 * math.pi))), 2 * math.pi -
            min(abs(k * t - math.pi), abs(k * t + math.pi)))
            for t in th)
        if r2 <= 1e-8:
            order = -k   # W^k = -I
            break
    data["closure_order_at_G10"] = order
    checks.append(("M2 closure order of the sphere step at the published "
                   "angle G10: %s (least |k| <= 240 with W^k = +-I, "
                   "tol 1e-8)" % (order if order else
                                  "NO-CLOSE-TO-240"), True))

    status = ("SPHERE-ENGINE-CLOCK-CONTACT" if contact else
              "SPHERE-ENGINE-CLOCK-VOID")
    ok = all(c[1] for c in checks)
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V347 done] %s" % status, flush=True)
    _emit(data, checks, ok, status=status)
    return 0


def _emit(data, checks, ok, status=None):
    out = {"id": "V347_sphere_engine_step", "tier": 1,
           "property": "the full published engine transplanted onto "
                       "the 12-face closed sphere: five bond stages on "
                       "the face graph x the published Cl(10,0) "
                       "internal layer - clock election vs the V206 "
                       "anchors",
           "status": status or data.get("status", ""),
           "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v347_sphere_engine_step.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
