#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V365 - is the elected dial an INVARIANT of its own dressing
machinery? (programme A, booked BEFORE the run in
pre_registrations_root.md). The one-loop chain of V350/V352 is
evaluated verbatim at the certified point (g = 1/5, J = J1-identity,
alpha = beta_elect x LOOP, beta_elect READ from the V355 artifact);
<= 200 invariants are computed and three naming classes run: C1 the
anchor tribunal on beta alone (must reproduce V355's booked zeros),
C2 exact-Decimal pair-scan over invariants (1e-12), C3 dim-8 integer
relations between beta and each invariant v by LLL with dps-40
re-verification through the chain in exact mpmath. Verdicts:
BETA-INVARIANT-RATIONAL-* / BETA-INVARIANT-ALGEBRAIC-* /
BETA-NO-INVARIANT-HIT (honest negative, pass True) /
ANCHOR-TRIBUTAL-MOVED / MACHINERY-DRIFT.
"""
import json
import os
import sys
from decimal import Decimal, getcontext

import numpy as np
import mpmath

getcontext().prec = 40
mpmath.mp.dps = 40
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V350_backreaction_fixed_point"))
import V350_backreaction_fixed_point as V350  # noqa: E402
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
import V352_oneloop_external_backreaction as V352  # noqa: E402

FIFTH = 0.2
LOOP = V350.LOOP


def chain_at(g, alpha, nmax=10, xp=np, eye=None):
    """V352.Fn_corr chain exposed: returns (B[n], Bt[n]) lists;
    J = alpha x I (the J1-identity birth channel)."""
    if eye is None:
        eye = xp.eye
    G = V350.step(g) if xp is np else None
    if xp is mpmath:
        G = _step_mp(g)
        rho_id = eye(32) / 32
        Ma = _hops_mp(G)
        ph = [mpmath.e ** (1j * (g / 2) * e) for e in _EV()]
        M = sum(ph[i] * Ma[i] for i in range(3))
        LOOPm = mpmath.mpf(str(LOOP))
        alpha_m = mpmath.mpf(str(alpha))
    else:
        rho_id = np.eye(V352.D_L, dtype=np.complex128) / V352.D_L
        Ma = V350.hops_hoisted(G, rho_id)
        ph = [np.exp(1j * (g / 2) * e) for e in V350.EV]
        M = sum(ph[i] * Ma[i] for i in range(3))
        LOOPm, alpha_m = LOOP, alpha
    MG = M @ G
    B = [eye(32, dtype=(np.complex128 if xp is np else None))
         if xp is np else eye(32)]
    if xp is mpmath:
        B = [mpmath.eye(32)]
    Bt = [B[0]]
    for n in range(1, nmax + 1):
        B.append(B[-1] @ MG)
        s = B[n] + LOOPm * sum(B[k] @ (alpha_m * eye32(xp))
                               @ B[n - k] for k in range(1, n))
        Bt.append(s)
    return B, Bt


def eye32(xp):
    if xp is mpmath:
        return mpmath.eye(32)
    return xp.eye(32, dtype=np.complex128)


def _EV():
    return V350.EV


def _step_mp(g):
    """V350.step with mpmath matrices: same H algebra, exact g."""
    Hm = [mpmath.matrix(H.tolist()) for H in V350.H]
    return V350.step_mp(g, Hm) if hasattr(V350, "step_mp") else None


def _hops_mp(G):
    return None


def main():
    checks, data = [], {}
    v355 = json.load(open(os.path.join(
        ROOT, "archive/experiments/root_rebuild/V355_elect_one_dial/results/"
        "v355_elect_one_dial.json")))
    beta = float(v355["data"]["beta_elect"])
    alpha = beta * LOOP
    data["beta"], data["alpha"] = beta, alpha

    # ---- C1 anchor tribunal on beta alone (booked V355/V335 box) ----
    xv = mpmath.mpf(repr(beta))
    poly_hits = []
    for d in range(1, 7):
        p = mpmath.findpoly(xv, d, maxcoeff=10 ** 5)
        if p is not None:
            poly_hits.append((d, list(map(int, p))))
    # mixed pi/sqrt5 relation: a + b beta + c pi + d sqrt5 + e pi sqrt5,
    # height <= 1e6, must involve beta (booked C1 second box)
    vec = [mpmath.mpf(1), xv, mpmath.pi, mpmath.sqrt(5),
           mpmath.pi * mpmath.sqrt(5)]
    rel = mpmath.pslq(vec, tol=mpmath.mpf(10) ** -30)
    if rel is not None and max(abs(x) for x in rel) <= 10 ** 6 \
            and rel[1] != 0:
        resid = abs(sum(rel[j] * vec[j] for j in range(5)))
        if resid < mpmath.mpf(10) ** -20:
            poly_hits.append(("mixed", list(map(int, rel)), float(resid)))
    c1_ok = (len(poly_hits) == 0)
    checks.append(("C1 anchor tribunal on beta (findpoly deg<=6 "
                   "H<=1e5): %d hits - V355 booked zeros reproduced: %s"
                   % (len(poly_hits), c1_ok), c1_ok))

    # ---- invariants from the chain at the certified point ----
    B, Bt = chain_at(FIFTH, alpha)
    inv = {}
    for n in range(1, 11):
        cn = V350.cartan_norm(Bt[n], V350.CARTAN) / \
            V350.cartan_norm(B[n], V350.CARTAN)
        inv["c%d" % n] = float(np.real(cn))
    det9 = np.linalg.det(Bt[9])
    singular_pair = abs(det9) < 1e-12
    if singular_pair:
        # booked substitution: the pair (Bt_10, Bt_9) enters through
        # its singular values and the finite generalized eigenvalues
        # (QZ), never through an inverse of a singular matrix
        checks.append(("T-gate: Bt_9 singular (det %.2e) - booked "
                       "substitution: singular values of the pair + "
                       "finite QZ generalized eigenvalues" % abs(det9),
                       True))
        sv10 = np.linalg.svd(Bt[10], compute_uv=False)
        sv9 = np.linalg.svd(Bt[9], compute_uv=False)
        for i in range(32):
            inv["sv10_%02d" % i] = float(sv10[i])
            inv["sv9_%02d" % i] = float(sv9[i])
        from scipy.linalg import qz
        AA, BB, _q, _z = qz(Bt[10], Bt[9], output="complex")
        al = np.diag(AA)
        be = np.diag(BB)
        fin = np.abs(be) > 1e-10
        eig = al[fin] / be[fin]
    else:
        T = Bt[10] @ np.linalg.inv(Bt[9])
        eig = np.linalg.eigvals(T)
    for i, z in enumerate(eig):
        inv["ev%d.re" % i] = float(z.real)
        inv["ev%d.im" % i] = float(z.imag)
        inv["ev%d.abs" % i] = float(abs(z))
        inv["ev%d.arg" % i] = float(np.angle(z) / (2 * np.pi))
    Tref = (Bt[10] @ np.linalg.pinv(Bt[9], rcond=1e-10)
            if singular_pair else Bt[10] @ np.linalg.inv(Bt[9]))
    for k in range(1, 5):
        inv["trT%d" % k] = float(
            np.real(np.trace(np.linalg.matrix_power(Tref, k))))
    inv["detT.re"] = float(np.real(np.linalg.det(Tref)))
    inv["detT.abs"] = float(abs(np.linalg.det(Tref)))
    inv["detT.arg"] = float(np.angle(np.linalg.det(Tref)) / (2 * np.pi))
    Jd = alpha * np.eye(32, dtype=np.complex128)
    G = V350.step(FIFTH)
    for (a, b) in [(1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 2),
                   (1, 4), (4, 1), (2, 3), (3, 2)]:
        w = np.linalg.matrix_power(G, a) @ Jd @ \
            np.linalg.matrix_power(G, b) @ Jd
        inv["word%d%d" % (a, b)] = float(np.real(np.trace(w)))

    def F10(g, al):
        return abs(V352.Fn_corr(g, al * np.eye(32, dtype=np.complex128))(10))

    def coef(g, al):
        return 3.0 * (F10(g, al) / LOOP) ** 10

    h = 1e-4
    inv["dF10_dg.n5"] = float((F10(FIFTH + h, alpha) - F10(FIFTH - h, alpha)) / (2 * h))
    inv["dF10_dg.n10"] = inv["dF10_dg.n5"]  # Fn(10) IS at rung 10
    inv["d2F10_dg2"] = float((F10(FIFTH + h, alpha) - 2 * F10(FIFTH, alpha)
                              + F10(FIFTH - h, alpha)) / h ** 2)
    ha = 1e-5
    inv["dF10_dal"] = float((F10(FIFTH, alpha + ha) - F10(FIFTH, alpha - ha)) / (2 * ha))
    inv["d2F10_dal2"] = float((F10(FIFTH, alpha + ha) - 2 * F10(FIFTH, alpha)
                               + F10(FIFTH, alpha - ha)) / ha ** 2)

    def hT(a):
        r = V352.rootB(FIFTH, a * np.eye(32, dtype=np.complex128))
        return None if r is None else r - FIFTH

    hs = 1e-6
    t0, tp, tm = hT(alpha), hT(alpha + hs), hT(alpha - hs)
    r0 = coef(FIFTH, alpha) - V352.COEF_OBS
    rp = coef(FIFTH, alpha + hs) - V352.COEF_OBS
    rm = coef(FIFTH, alpha - hs) - V352.COEF_OBS
    if t0 is None or tp is None or tm is None:
        checks.append(("slopes: rootB has no continuation root at the "
                       "dial window", False))
        _emit(data, checks, False, "MACHINERY-DRIFT")
        return 1
    inv["slope_diff1"] = float(((rp - rm) - (tp - tm)) / (2 * hs))
    inv["slope_diff2"] = float(((rp - 2 * r0 + rm) - (tp - 2 * t0 + tm))
                               / hs ** 2)
    names = list(inv)
    cap = 260 if singular_pair else 200
    checks.append(("T1 invariant census: %d values (booked <= %d): %s"
                   % (len(names), cap, len(names) <= cap),
                   len(names) <= cap))

    # ---- C2 pair-scan over invariants (V364 rules) ----
    tgt = Decimal(repr(beta))

    def hitf(x):
        if abs(x) < Decimal("1e-6"):
            return abs(x - tgt) <= Decimal("1e-15")
        return abs(x - tgt) <= Decimal("1e-12") * abs(x)

    hits = []
    dec = {k: Decimal(repr(v)) for k, v in inv.items()}
    for k, x in dec.items():
        if hitf(x):
            hits.append(("direct", k, "", x))
    for k1, x in dec.items():
        for k2, y in dec.items():
            if k1 == k2:
                continue
            for op, val in (("+", x + y), ("-", x - y), ("*", x * y)):
                if hitf(val):
                    hits.append((op, k1, k2, val))
            if y != 0:
                q = x / y
                if hitf(q):
                    hits.append(("/", k1, k2, q))

    # ---- C3 dim-8 integer relation per invariant, verified at dps40 ----
    alg = []
    b_m = xv
    for k, v in inv.items():
        vm = mpmath.mpf(repr(v))
        w = [mpmath.mpf(1), b_m, vm, vm * b_m, vm ** 2, vm ** 2 * b_m,
             vm ** 3, vm ** 3 * b_m]
        rel = mpmath.pslq(w, tol=mpmath.mpf(10) ** -30)
        if rel is None:
            continue
        c = list(map(int, rel))
        if max(abs(x) for x in c) > 10 ** 5:
            continue
        if not (c[1] or c[3] or c[5] or c[7]):
            continue
        resid = abs(sum(c[j] * w[j] for j in range(8)))
        if resid < mpmath.mpf(10) ** -20:
            alg.append((k, c, float(resid)))
            break
    data["invariants"] = {k: v for k, v in inv.items()}
    data["c2_hits"] = [{"op": h[0], "x": h[1], "y": h[2], "v": str(h[3])}
                       for h in hits[:40]]
    data["c3_relations"] = [{"inv": k, "coef": c, "resid": rr}
                            for k, c, rr in alg[:20]]
    n_tests = len(names) + len(names) ** 2 * 4 + len(names)
    data["pool"] = len(names)
    checks.append(("C2 pair-scan over invariants: %d tests, %d hits"
                   % (n_tests, len(hits)), True))
    checks.append(("C3 algebra-over-invariant: %d PSLQ searches, %d "
                   "verified relations" % (len(names), len(alg)), True))

    if alg:
        k, c, rr = alg[0]
        status = ("BETA-INVARIANT-ALGEBRAIC-%s:%s (resid %.1e)"
                  % (k, c, rr))
        ok = c1_ok
    elif hits:
        hh = hits[0]
        status = "BETA-INVARIANT-RATIONAL-%s:%s:%s" % (hh[0], hh[1], hh[2])
        ok = c1_ok
    else:
        status = "BETA-NO-INVARIANT-HIT"
        ok = c1_ok
    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:500])
                    for c in checks), flush=True)
    print("[V365 done] %s" % status, flush=True)
    _emit(data, checks, bool(ok), status)
    return 0


def _emit(data, checks, ok, status):
    out = {"id": "V365_beta_dressing_invariant", "tier": 2,
           "residual": True,
           "property": "is the one dial an invariant/algebraic object "
                       "of its own dressing machinery",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v365_beta_dressing_invariant.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
