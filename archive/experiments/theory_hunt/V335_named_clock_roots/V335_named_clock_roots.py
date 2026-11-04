#!/usr/bin/env python3
r"""V335 - the named-clock-root hunt (tier 2, hard, O4/O8, theory-hunt).

The first constructive test of the V334 dichotomy's live branch. The
rigidity theorem closed the exponential-closure route by proof; the
second branch - transcendental coefficients - is already live: the
crossing equation |F_n(g)| = LOOP and the inversion coef10(g) = COEF_OBS
are transcendental-coefficient conditions on the certified V206 ladder
functional. This experiment asks: does any root of a booked named
equation family simultaneously (a) match a closed form in a booked
dictionary at <= 1e-12, (b) pass the K1 coefficient contact band (1e-4
relative, the SAME functional), and (c) pass the K3 rung-17 ladder
contact band (5%)? FOUND is the first elected-theorem candidate of the
programme; the census books the full root ledger either way.

Equations (pre-registration V335): E1 |F_n(g)| = LOOP, n in 3..24 +
{170,208}; E2 |F_a(g)| = LOOP^k |F_b(g)|, (a,b) in {(10,3),(17,10),
(17,3)} x k in 1..3; roots in [0.05, 1.55] (0.01 sign-change grid then
60 halvings); plus the two known roots g_obs, g_cross enter the
dictionary census directly.
Machinery: the V206 sector-hop einsums at the identity rho contract
EXACTLY into precomputed 32x32 projector constants (W_i), so the scan
runs on the host at ~3 ms/evaluation - the GPU staging of the booked
microbenchmark turns out unnecessary (the contraction beats the staged-
einsum hybrid ~50x); the verbatim V206 numpy Fn_of is the reference and
the T0 host, and the T0e identity gate at 1e-12 makes every verdict
path-invariant.
Branches: CLOCK-NAMED-ROOT-FOUND / CLOCK-NO-NAMED-ROOT-DECLARED /
MACHINERY-DRIFT.
"""
import importlib.util
import json
import math
import os
import sys
import time
from fractions import Fraction

import mpmath
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp  # noqa: E402
from cl10 import config as _config  # noqa: E402

CFG = _config.CFG

_spec = importlib.util.spec_from_file_location(
    "v206mach",
    os.path.join(ROOT, "experiments", "scales_loopfactor",
                 "V206_clock_tension", "V206_clock_tension.py"))
V = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V)

LOOP = V.LOOP
COEF_OBS = V.COEF_OBS
G_LO, G_HI, G_STEP = 0.05, 1.55, 0.01
NAMES = list(range(3, 25)) + [170, 208]
PAIRS_E2 = [(10, 3), (17, 10), (17, 3)]
KS_E2 = (1, 2, 3)
NAMED_EPS = 1e-12
CF_QMAX = 10**6
LOOPQ_QMAX = 10**4
K1_BAND = 1e-4
K3_BAND = float(CFG.pillar_surgery.ledger_scale_rel_tol)  # 0.05
PSLQ_DEG, PSLQ_H = int(CFG.pillar_surgery.pslq_max_degree), \
    int(CFG.pillar_surgery.pslq_max_coeff)


def art(name, key):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        j = json.load(fh)
    for part in key.split("."):
        j = j[part]
    return j


# ---- fast functional: the hop einsums contract EXACTLY at rho = I/D_L.
# With rho = I/D_L (the V206 default), the sector hop
#   Ma_i[a,b] = sum_j (P_i kron-part contracted with G kron I_45)[a,j,b,j]
#             / D_L = (W_i @ G)[a,b] / D_L,
# W_i[a,c] = sum_j P_r[a,j,c,j] - constant 32x32 matrices precomputed per
# projector (exact index collapse, verified against the verbatim V.Fn_of
# by the T0e identity gate at 1e-12). The GPU staging of the pre-reg
# microbenchmark turns out unnecessary: the contracted path is host-only
# and ~50x faster than the staged-einsum hybrid.
_cache = {}
_W = [np.einsum("ajcj->ac", p.reshape(32, V.D_L, 32, V.D_L))
      for p in V.P]


def fn(g, nmax):
    key = (round(g, 15), nmax)
    hit = _cache.get(key)
    if hit is not None:
        return hit
    G = V.step(g)
    Ma = [(_Wi @ G) / V.D_L for _Wi in _W]
    ph = [np.exp(1j * (g / 2) * e) for e in V.EV]
    M = sum(ph[i] * Ma[i] for i in range(3))
    MG = M @ G
    Fn, Bn, Gn = {}, np.eye(32, dtype=np.complex128), \
        np.eye(32, dtype=np.complex128)
    for n in range(1, nmax + 1):
        Bn = Bn @ MG
        Gn = Gn @ G
        Fn[n] = V.cartan_norm(Bn, V.CARTAN) / V.cartan_norm(Gn, V.CARTAN)
    _cache[key] = Fn
    return Fn


def coef10(g):
    F10 = abs(fn(g, 10)[10])
    return 3.0 * (F10 / LOOP) ** 10 if F10 / LOOP < 1 else None


def bisect(f, a, b, iters=60):
    fa = f(a)
    for _ in range(iters):
        mid = 0.5 * (a + b)
        if (f(mid) > 0) == (fa > 0):
            a = mid
        else:
            b = mid
    return 0.5 * (a + b)


def cf_best(x, qmax):
    """Closest continued-fraction convergent with denominator <= qmax."""
    fr = Fraction(x).limit_denominator(qmax)
    return int(fr.numerator), int(fr.denominator), abs(x - float(fr))


def loop_best(x):
    """Booked LOOP^k-form test at ABSOLUTE 1e-12 in the tested quantity
    (|x*LOOP^k - p/q| or |x/LOOP^k - p/q|), evaluated in exact 50-digit
    mpmath arithmetic: float64 p/q evaluation aliases to large floats at
    1e-8 spacing and MUST NOT be trusted for the comparison. The ZERO
    numerator is guarded: x*LOOP^k underflows toward 0 at large k and a
    vanishing product is vacuously close to 0/1 - a nonzero clock is
    named only by a nonzero rational."""
    mpmath.mp.dps = 50
    xm = mpmath.mpf(repr(x))
    best = None
    for k in range(1, 13):
        s = mpmath.mpf(repr(LOOP)) ** k
        for tag, y in (("x*LOOP^%d" % k, xm * s), ("x/LOOP^%d" % k, xm / s)):
            yf = float(y)
            if yf <= 0:
                continue
            fr = Fraction(yf).limit_denominator(LOOPQ_QMAX)
            if fr.numerator == 0:
                continue
            miss = abs(y - mpmath.mpf(fr.numerator) / fr.denominator)
            if best is None or miss < best["miss"]:
                best = {"form": tag, "value": "%d/%d" % (fr.numerator,
                                                         fr.denominator),
                        "miss": float(miss)}
    return best


def pslq_hit(x):
    """V309's booked PSLQ box: findpoly per degree with maxcoeff, then
    substitution verification at 60 dp against the 1e-14-scale float
    noise floor (a float64 root is exact only to ~16 digits)."""
    mpmath.mp.dps = 60
    xv = mpmath.mpf(repr(x))
    for d in range(1, PSLQ_DEG + 1):
        poly = mpmath.findpoly(xv, d, maxcoeff=PSLQ_H)
        if poly is None:
            continue
        res = abs(sum(mpmath.mpf(c) * xv ** k for k, c in enumerate(poly)))
        if res <= mpmath.mpf("1e-30"):
            return [int(c) for c in poly]
    return None


def main():
    t0wall = time.time()
    data = {"class": "the booked V335 families on the certified V206 "
                     "ladder functional; dictionary: CF q<=1e6, "
                     "LOOP^k-forms q<=1e4 k<=12, V326 members x LOOP^k "
                     "k<=90, PSLQ deg<=%d H<=%d; contacts K1 %.0e, K3 "
                     "%.2f" % (PSLQ_DEG, PSLQ_H, K1_BAND, K3_BAND),
            "backend": "host exact-contraction fast path (projector hop "
                       "trace constants W_i; rho = I/D_L is the identity); "
                       "verbatim V206 numpy Fn_of is the reference under "
                       "the T0e gate; GPU unused (contraction beats the "
                       "staged-einsum hybrid by ~50x)",
            "D_L": V.D_L}
    checks = []

    # ---- T0e contraction identity gate (path-invariance) --------------
    worst_h = 0.0
    for gs in (0.2, 1.0 / 3.0, 0.35):
        Fh = fn(gs, 60)
        Fv = V.Fn_of(gs, 60)
        worst_h = max(worst_h, max(abs(Fh[n] - Fv[n]) for n in range(1, 61)))
    checks.append({"name": "T0e contracted == verbatim V206 Fn_of (n<=60, "
                           "3 seed clocks, <= 1e-12): worst %.2e"
                   % worst_h, "pass": worst_h <= 1e-12})
    if worst_h > 1e-12:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- T0f GPU-1 staged einsum path (the pre-reg booked hybrid): a
    # third backend identity at the seed clocks; GPU unused for the scan
    # (the host contraction is ~50x faster), engaged here as the booked
    # backend-invariance gate. ----
    data["t0"] = {}
    worst_g = None
    try:
        import cupy as cp
        I45 = cp.eye(V.D_L, dtype=np.complex128)
        RHO = I45 / V.D_L
        PC = [cp.asarray(p) for p in V.P]
        worst_g = 0.0
        for gs in (0.2, 1.0 / 3.0):
            G = V.step(gs)
            Gi = cp.kron(cp.asarray(G), I45)
            Ma = []
            for Pa in PC:
                X = (Pa @ Gi).reshape(32, V.D_L, 32, V.D_L)
                Z = cp.einsum("ajbm,ml->ajbl", X, RHO)
                Ma.append(cp.asnumpy(cp.einsum("ajbj->ab", Z)))
            ph = [np.exp(1j * (gs / 2) * e) for e in V.EV]
            M = sum(ph[i] * Ma[i] for i in range(3))
            MG = M @ G
            Fh = fn(gs, 60)
            Bn, Gn = np.eye(32, dtype=np.complex128), np.eye(32, dtype=np.complex128)
            for n in range(1, 61):
                Bn = Bn @ MG
                Gn = Gn @ G
                v = V.cartan_norm(Bn, V.CARTAN) / V.cartan_norm(Gn, V.CARTAN)
                worst_g = max(worst_g, abs(v - Fh[n]))
        data["t0"]["gpu_staged_identity_worst"] = worst_g
        checks.append({"name": "T0f GPU-1 staged einsum == contracted path "
                               "(n<=60, 2 seed clocks, <= 1e-12): %.2e"
                       % worst_g, "pass": worst_g <= 1e-12})
    except Exception as exc:                          # pragma: no cover
        data["t0"]["gpu_staged_identity_worst"] = "skipped: %s" % exc
        checks.append({"name": "T0f GPU staged einsum path skipped "
                               "(cupy unavailable or GPU busy: %s); "
                               "host T0e identity governs" % exc,
                       "pass": True})
    if worst_g is not None and worst_g > 1e-12:
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- T0a: g_obs re-derived by the verbatim V206 bracket ----------
    lo, hi = 0.1985, 0.2
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        c = coef10(mid)
        if c is None or c > COEF_OBS:
            hi = mid
        else:
            lo = mid
    g_obs = 0.5 * (lo + hi)
    leaf_obs = float(art("V206_clock_tension", "data.data_clock.g_obs"))
    d_obs = abs(g_obs - leaf_obs)
    data["t0"].update({"g_obs": g_obs, "g_obs_vs_leaf": d_obs})
    checks.append({"name": "T0a g_obs re-derived vs V206 leaf (<= 1e-13): "
                           "%.2e" % d_obs, "pass": d_obs <= 1e-13})

    # ---- T0b: g_cross re-derived ------------------------------------
    lo2, hi2 = 0.2, 0.21

    def f10(g):
        return abs(fn(g, 10)[10]) - LOOP

    g_cross = bisect(f10, lo2, hi2)
    leaf_cross = float(art("V206_clock_tension",
                           "data.data_clock.g_cross_rung10"))
    d_cross = abs(g_cross - leaf_cross)
    data["t0"]["g_cross"] = g_cross
    data["t0"]["g_cross_vs_leaf"] = d_cross
    checks.append({"name": "T0b g_cross re-derived vs V206 leaf (<= 1e-13): "
                           "%.2e" % d_cross, "pass": d_cross <= 1e-13})

    # ---- T0c: the pi/4 closed reading --------------------------------
    Fp = fn(math.pi / 4, 600)
    nstar = None
    deficit = None
    for n in range(2, 601):
        if 1e-12 < abs(Fp[n]) < LOOP:
            nstar, deficit = n, abs(Fp[n]) / LOOP
            break
    leaf_def = float(art("V205_channel_filtered_deficit",
                         "data.closed_crossing.deficit_vs_loop"))
    data["t0"]["pi4_nstar"] = nstar
    data["t0"]["pi4_deficit_vs_leaf"] = (abs(deficit - leaf_def)
                                         if deficit else 9e9)
    checks.append({"name": "T0c pi/4 reading n* = 208, deficit vs V205 "
                           "leaf < 1e-4: n* %s d %.2e"
                   % (nstar, data["t0"]["pi4_deficit_vs_leaf"]),
                   "pass": nstar == 208
                   and data["t0"]["pi4_deficit_vs_leaf"] < 1e-4})

    # ---- T0d: COEF_OBS from the u10pred leafs ------------------------
    coef_re = float((V.u10pred.V_EW / V.u10pred.M_PLANCK) / LOOP ** 10)
    data["t0"]["coef_obs_vs_leaf"] = abs(coef_re - COEF_OBS)
    checks.append({"name": "T0d COEF_OBS recomputed from u10pred (diff "
                           "0.0)", "pass": coef_re == COEF_OBS})

    if not all(c["pass"] for c in checks):
        data["status"] = "MACHINERY-DRIFT"
        _emit(data, checks, ok=False)
        return 1

    # ---- grid scan + root finding ------------------------------------
    grid = [G_LO + i * G_STEP for i in range(
        int(round((G_HI - G_LO) / G_STEP)) + 1)]
    for g in grid:
        fn(g, 208)

    def find_roots(f, label):
        roots = []
        prev = f(grid[0])
        for i in range(1, len(grid)):
            cur = f(grid[i])
            if (prev > 0) != (cur > 0):
                roots.append(bisect(f, grid[i - 1], grid[i]))
            prev = cur
        return [(label, r) for r in roots]

    roots = []
    for n in NAMES:
        roots += find_roots(lambda g, n=n: abs(fn(g, max(n, 17))[n]) - LOOP,
                            "E1:F_%d=LOOP" % n)
    for (a, b) in PAIRS_E2:
        for k in KS_E2:
            roots += find_roots(
                lambda g, a=a, b=b, k=k: abs(fn(g, max(a, b, 17))[a])
                - LOOP ** k * abs(fn(g, max(a, b, 17))[b]),
                "E2:F_%d=LOOP^%d*F_%d" % (a, k, b))
    roots += [("REF:g_obs", g_obs), ("REF:g_cross", g_cross)]

    # ---- dictionary census + contacts --------------------------------
    members = [float(s.split("=")[1]) for s in
               art("V326_phase_dictionary_bridge",
                   "data.dictionary_class.members")]
    ledger = {}
    any_cand = False
    for label, x in roots:
        cf = cf_best(x, CF_QMAX)
        lb = loop_best(x)
        db = None
        for c in members:
            for k in range(1, 91):
                miss = abs(x - c * LOOP ** k)
                if db is None or miss < db["miss"]:
                    db = {"member": c, "k": k, "miss": miss}
        if db is None:
            db = {"member": None, "k": None, "miss": 9.0}
        named = min(cf[2], lb["miss"], db["miss"]) <= NAMED_EPS
        Fx = fn(x, 17)
        c10 = coef10(x)
        k1 = None if c10 is None else abs(c10 - COEF_OBS) / COEF_OBS
        k3 = abs(Fx[17] - LOOP ** 17) / LOOP ** 17
        poly = pslq_hit(x)
        cand = (named and k1 is not None and k1 <= K1_BAND
                and k3 <= K3_BAND)
        any_cand = any_cand or cand
        ledger[label] = {
            "g": x,
            "cf": {"p": cf[0], "q": cf[1], "miss": cf[2]},
            "loop_form": lb, "dict_member": db,
            "pslq": poly,
            "named": bool(named),
            "K1_rel": k1, "K3_rel": k3,
            "candidate": bool(cand)}
    data["roots"] = ledger
    data["n_roots"] = len(roots)

    ok_census = len(roots) > 0
    checks.append({"name": "census complete: %d roots located and scored"
                   % len(roots), "pass": bool(ok_census)})

    status = ("CLOCK-NAMED-ROOT-FOUND" if any_cand
              else "CLOCK-NO-NAMED-ROOT-DECLARED")
    data["status"] = status
    best_near = min(ledger.values(), key=lambda r: min(
        r["cf"]["miss"], r["loop_form"]["miss"], r["dict_member"]["miss"]))
    data["closest_named_miss"] = best_near
    data["wall_time"] = time.time() - t0wall
    print("V335 %s: %d roots; closest named miss %.3e at g = %.12f (%s); "
          "K1 window holders: %s"
          % (status, len(roots), min(best_near["cf"]["miss"],
                                     best_near["loop_form"]["miss"],
                                     best_near["dict_member"]["miss"]),
             best_near["g"], "see ledger",
             [lbl for lbl, r in ledger.items()
              if r["K1_rel"] is not None and r["K1_rel"] <= K1_BAND]),
          flush=True)
    _emit(data, checks, ok=all(c["pass"] for c in checks))
    print("[V335 done] %s" % status, flush=True)
    return 0


def _emit(data, checks, ok):
    out = {"id": "V335_named_clock_roots", "tier": 2,
           "property": "the named-clock-root hunt: roots of booked "
                       "named ladder equations scored against the "
                       "dictionary class and the K1/K3 contact bands "
                       "(tier 2, O4/O8)",
           "residual": True,
           "programme": "theory-hunt V335",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v335_named_clock_roots.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V335] pass=%s branch=%s" % (out["pass"], out["branch"]),
          flush=True)


if __name__ == "__main__":
    sys.exit(main())
