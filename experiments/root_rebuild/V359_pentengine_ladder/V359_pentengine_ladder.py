#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V359 - the traced-word ladder OF THE PENTAGON ENGINE (booked
before the run in review/pre_registrations_root.md; tier 1).

F_n^P(g) = tr U_hop(g)^n / 36 = sum_j exp(i g n lambda_j) / 36 - the
closed-form ladder on the certified 36-cell hop spectrum (V358),
basis-free, no new eigensolves. T0: the V357 artifact's closed-walk
leaves W_n re-derived bit-for-bit from the face adjacency (the ladder
AT g = 0 of the 36-cell is 1 by normalization - the frame factor
between the 36-cell m-th moment tr K^m and W_m is measured and
printed: m = 2 gives the booked 3 = 180/60). Crossings: bisection to
1e-9 on ||F_10^P(g)| - LOOP| over the booked 1e-3 grid on (0,1), all
roots listed; any root within +-1e-4 of 1/5 named by continued
fraction (the booked cheap class q <= 1e4) against the V335 booked
dictionary class, printed density kept. Reading: C_P =
3(|F_10^P(1/5)|/LOOP)^10 against V352.COEF_OBS, window = the V355
artifact read_dev leaf (no hard-coded decimal). GPU-1 directive: the
(g, n) scan batched on cupy with the house 5-probe device-vs-host
identity gate <= 1e-12 (V352 T0e pattern).
"""
import json
import os
import sys

os.environ["CL10_ANGLE"] = "0.2"    # house pin BEFORE importing cl10

import numpy as np                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
from cl10.paths import exp          # noqa: E402
import V352_oneloop_external_backreaction as V352  # noqa: E402

FIFTH = 0.2

_CP = None


def _get_cp():
    global _CP
    if _CP is None:
        try:
            import cupy
            cupy.cuda.Stream  # noqa: B018  (probe)
            _CP = cupy
        except Exception:      # noqa: BLE01
            _CP = False
    return _CP or None


def _cf_name(x, qmax=10 ** 4):
    """Best continued-fraction name p/q with q <= qmax (the booked
    cheap class of the V335 dictionary)."""
    from fractions import Fraction
    fr = Fraction(x).limit_denominator(qmax)
    return int(fr.numerator), int(fr.denominator), abs(x - float(fr))


def main():
    checks, data = [], {}
    from cl10 import pentengine as P

    lam = P.hop_spectrum()          # certified 36 real eigenvalues of K

    # ---- T0: the V357 walk ladder re-derived; frame factor measured ---
    with open(exp("V357_pentagon_walk_ladder", "results",
                  "v357_pentagon_walk_ladder.json")) as fh:
        v357 = json.load(fh)
    W_art = {int(k): int(v) for k, v in
             v357["data"]["walks_n_1_30"].items()}
    vs, edges = P.face_geometry()
    A = np.zeros((12, 12), dtype=object)      # exact Python ints:
    for (u, v) in edges:                      # W_28 overflows int64
        A[u, v] += 1
        A[v, u] += 1
    W = {}
    M = A.copy()
    W[1] = int(sum(M[i][i] for i in range(12)))
    for n in range(2, 31):
        M = A @ M
        W[n] = int(sum(M[i][i] for i in range(12)))
    same = all(W[n] == W_art.get(n) for n in range(1, 31))
    F0 = [np.mean(np.exp(0j * lam * m)) for m in range(2, 11)]
    f0_one = all(abs(f - 1.0) <= 1e-15 for f in F0)
    trK2 = float(np.sum(lam * lam))
    frame2 = trK2 / W[2]
    data["frame_factor_m2"] = frame2
    data["walks_rederived"] = {str(n): W[n] for n in range(1, 11)}
    checks.append(("T0 the V357 artifact closed-walk leaves W_1..W_30 "
                   "re-derived bit-for-bit from the face adjacency: %s; "
                   "the 36-cell ladder at g = 0 is F_m^P(0) = 1 exactly "
                   "(m = 2..10, dev %.1e) - the booked frame factor: "
                   "36 F_m^P(0) = 36 vs W_m, and the m = 2 moment ratio "
                   "tr K^2 / W_2 = %.6f (booked 3 = 180/60)"
                   % ("bit-identical" if same else "MISMATCH",
                      max(abs(f - 1.0) for f in F0), frame2),
                   same and f0_one and abs(frame2 - 3.0) <= 1e-9))

    # ---- GPU batched (g, n) scan + house identity gate ----
    gs = np.arange(1, 1000) * 1e-3                    # (0,1), 1e-3 grid
    ns = np.arange(2, 15)
    cp = _get_cp()
    devname = "CPU fallback"
    if cp is not None:
        import cupy
        devname = "cupy " + str(cupy.cuda.runtime.getDeviceProperties(
            cupy.cuda.Device().id)["name"])
    if cp is not None:
        G = cp.asarray(gs)
        L = cp.asarray(lam)
        NS = cp.asarray(ns.astype(np.float64))
        # F_n^P(g) = mean_j exp(i g n lambda_j): batched (n, g) grid
        Fg = cp.exp(1j * G[None, :, None] * NS[:, None, None]
                    * L[None, None, :]).mean(axis=2)   # (13, 999)
        Fabs = cp.asnumpy(cp.abs(Fg))
        idx10 = int(np.where(ns == 10)[0][0])
        F10 = Fabs[idx10]
    else:
        ph = np.exp(1j * np.outer(gs, lam))           # (999, 36)
        F10 = np.abs(np.array([np.mean(p ** 10) for p in ph]))
        Fabs = np.stack([np.abs(np.array([np.mean(p ** n) for p in ph]))
                         for n in ns])
    probe_i = [167, 199, 244, 302, 451]
    host = np.array([abs(np.mean(np.exp(1j * gs[i] * 10 * lam)))
                     for i in probe_i])
    gpu_id = float(np.max(np.abs(F10[probe_i] - host)))
    data["gpu_backend"] = {"backend": "cupy" if cp is not None
                           else "numpy", "device": devname}
    checks.append(("T0e GPU batched ladder scan vs host ladder (5 "
                   "probes, F_10 at the booked grid): max dev %.1e (<= "
                   "1e-12); device in use: %s" % (gpu_id, devname),
                   gpu_id <= 1e-12))

    # ---- crossings: |F_10^P(g)| = LOOP, bisection to 1e-9 ----
    LOOP = float(V352.LOOP)
    COEF = float(V352.COEF_OBS)

    def f10(g):
        return abs(np.mean(np.exp(1j * g * 10 * lam))) - LOOP

    roots = []
    for i in range(len(gs) - 1):
        a, b = gs[i], gs[i + 1]
        fa, fb = f10(a), f10(b)
        if fa == 0.0:
            roots.append(a)
        elif fa * fb < 0:
            lo, hi = a, b
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if f10(mid) * fa <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
    if f10(gs[-1]) == 0.0:
        roots.append(gs[-1])
    data["LOOP"] = LOOP
    data["n_roots"] = len(roots)
    data["roots_g"] = roots[:60]

    # ---- the fifth: named by CF (booked cheap class) ----
    near = [r for r in roots if abs(r - FIFTH) <= 1e-4]
    named = [(r,) + _cf_name(r) for r in near]
    data["roots_near_fifth"] = [{"g": r, "cf": "%d/%d" % (p, q),
                                 "cf_dev": d} for (r, p, q, d) in named]
    v335_class = json.load(open(exp(
        "V335_named_clock_roots", "results",
        "v335_named_clock_roots.json")))["data"]["class"]
    checks.append(("CROSSINGS: %d root(s) of |F_10^P(g)| = LOOP on the "
                   "booked 1e-3 grid (bisection 1e-9); within +-1e-4 of "
                   "the fifth 1/5: %d named by CF q <= 10000 (booked "
                   "cheap class of the V335 dictionary: %s)"
                   % (len(roots), len(named), v335_class), True))

    # ---- reading at the exact fifth vs COEF_OBS, window = V355 leaf ---
    F10f = abs(np.mean(np.exp(1j * FIFTH * 10 * lam)))
    C_P = 3.0 * (F10f / LOOP) ** 10
    rel = abs(C_P - COEF) / COEF
    with open(exp("V355_elect_one_dial", "results",
                  "v355_elect_one_dial.json")) as fh:
        window = float(json.load(fh)["data"]["battery"]["read_dev"])
    data["C_P_at_fifth"] = {"F10_abs": F10f, "C_P": C_P, "COEF_OBS": COEF,
                            "rel_dev": rel, "window_read_dev": window}
    checks.append(("READING: C_P = 3(|F_10^P(1/5)|/LOOP)^10 = %.6g vs "
                   "COEF_OBS = %.6g: rel dev %.3e (dial window = the "
                   "V355 read_dev leaf %.3e: %s)"
                   % (C_P, COEF, rel, window,
                      "INSIDE" if rel <= window else "outside"), True))

    crossing = bool(named) or rel <= window
    if not (checks[0][1] and checks[1][1]):
        status, ok = "MACHINERY-DRIFT", False
    elif crossing:
        status, ok = "PENTLADDER-CROSSING", True
    else:
        status, ok = "PENTLADDER-NO-CROSSING", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V359 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V359_pentengine_ladder", "tier": 1,
           "property": "the pentagon engine's own traced-word ladder: "
                       "rung-10 crossings at the fifth and the reading "
                       "against COEF_OBS inside the elected dial window",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v359_pentengine_ladder.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
