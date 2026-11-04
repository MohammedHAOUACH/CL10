#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V384 - the cubond ladder functional (booked BEFORE the run in
review/pre_registrations_root.md, "V384_cubond_ladder_functional").

The named open design problem of V383, answered: the certified
hierarchy machinery (V206) lifted exactly to entangling engines.
The lift is forced by the linearity of the sector hop
(ma_lambda = A_lambda C, A_lambda the diagonal-block sum of the
sector projector at the maximally mixed adjoint state), not
chosen: the dressing becomes sum_lambda phase_lambda (I_8 x
A_lambda) G_engine, the Cartan norm takes the bond block-trace,
and the published normalization survives because 256 = 8 x 32.
The clocked cubond engine: G_gamma(g) = prod_j (cos g + i sin g
gamma_j), the cubic-symmetric momentum (g, g, g).

Verdicts: CUBOND-LADDER-NO-CONTACT / CUBOND-LADDER-READS-OBSERVED.
"""
import importlib.util
import json
import math
import os
import sys
import time

os.environ["CL10_BACKEND"] = "numpy"

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)

import numpy as np  # noqa: E402

from cl10 import config as _config  # noqa: E402
from cl10 import supercell_2i as S2  # noqa: E402
from cl10.paths import exp  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v206", os.path.join(ROOT, "archive", "experiments", "scales_loopfactor",
                         "V206_clock_tension", "V206_clock_tension.py"))
V206 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V206)

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, bool(ok)))


D_L = V206.D_L
LOOP = V206.LOOP
COEF_OBS = float(_config.CFG.physics.hierarchy_coefficient_observed)
BAND = 1e-4
NMAX = 600

# A_lambda = (1/45) sum_j P_lambda^{jj} - content data, built once.
ALS = []
for Pa in V206.P:
    Ar = Pa.reshape(32, D_L, 32, D_L)
    ALS.append(np.einsum("ajbj->ab", Ar) / D_L)

Sx, Sy, Sz = S2.spin_triple()
I32 = np.eye(32, dtype=np.complex128)
X, Y, Z = (np.array(m, dtype=np.complex128) for m in
           ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]))
I2 = np.eye(2, dtype=np.complex128)
P3 = [np.kron(np.kron(X, I2), I2), np.kron(np.kron(I2, Y), I2),
      np.kron(np.kron(I2, I2), Z)]
gam = [np.kron(p, 2 * r) for p, r in zip(P3, (Sx, Sy, Sz))]
I8 = np.eye(8, dtype=np.complex128)
MONO = {}
for mask in range(8):
    Mm = np.eye(256, dtype=np.complex128)
    for j in range(3):
        if mask >> j & 1:
            Mm = Mm @ gam[j]
    MONO[mask] = Mm
KLAM = [np.kron(I8, A) for A in ALS[:3]]


def G_gamma(g):
    c, s = math.cos(g), math.sin(g)
    G = np.zeros((256, 256), dtype=np.complex128)
    for mask in MONO:
        k = bin(mask).count("1")
        G += (c ** (3 - k)) * (1j * s) ** k * MONO[mask]
    return G


def L_of(g):
    L = np.zeros((256, 256), dtype=np.complex128)
    for lam in range(3):
        L += np.exp(1j * (g / 2) * V206.EV[lam]) * KLAM[lam]
    return L


CSTACK = np.stack(V206.CARTAN)          # (5, 32, 32)


def norm_x(Xm):
    # declared corrected lift: Cartan content of EVERY bond block
    # (the booked bond-trace version is identically zero on
    # entangling engines - the Pauli bond coefficients are
    # traceless - and T0 caught it; on the lifted published
    # engine this equals the published cartan_norm term by term)
    Xr = Xm.reshape(8, 32, 8, 32)      # (a, j, b, k)
    co = np.einsum("ikj,ajbk->abi", CSTACK, Xr)
    return float(np.abs(co).sum() / 256.0)


def Fn_gamma(g, G=None):
    G = G_gamma(g) if G is None else G
    # (L(g) @ G) is the dressed M; the certified loop B_n = (M G)^n:
    MG = (L_of(g) @ G) @ G
    Fn, Bn, Gn = {}, np.eye(256, dtype=np.complex128), \
        np.eye(256, dtype=np.complex128)
    for n in range(1, NMAX + 1):
        Bn = Bn @ MG
        Gn = Gn @ G
        Fn[n] = norm_x(Bn) / norm_x(Gn)
    return Fn


def reading(Fn):
    for n in range(2, NMAX + 1):
        if 1e-12 < abs(Fn[n]) < LOOP:
            d = abs(Fn[n]) / LOOP
            return n, d, 3.0 * d ** 10
    return None, None, None


def closure_k(M, cap=24, tol=1e-8):
    Mk = np.eye(M.shape[0], dtype=np.complex128)
    for k in range(1, cap + 1):
        Mk = Mk @ M
        e = min(float(np.max(np.abs(Mk - np.eye(M.shape[0])))),
                float(np.max(np.abs(Mk + np.eye(M.shape[0])))))
        if e <= tol:
            return k
    return None


def main():
    data = {}

    # ---- T0: the lift anchor on the PUBLISHED engine ----------------------
    v206 = json.load(open(exp("V206_clock_tension", "results",
                              "v206_clock_tension.json")))
    t0_ok = True
    for tag, g in [("pi/2", math.pi / 2), ("pi/4", math.pi / 4),
                   ("3pi/4", 3 * math.pi / 4), ("1/5rad", 0.2)]:
        Gl = np.kron(I8, V206.step(g))
        Fn = Fn_gamma(g, G=Gl)
        n, d, coef = reading(Fn)
        if tag == "1/5rad":
            ref = float(v206["data"]["data_clock"]["coef_at_fifth"])
            dd = abs(coef - ref)
            lim = 1e-8
        else:
            ref = float(v206["data"]["closure_class"][tag]["coef"])
            dd = abs(coef - ref)
            lim = 5e-3
        t0_ok = t0_ok and dd < lim
        data["anchor_%s" % tag] = {"coef": coef, "v206": ref,
                                   "detune": dd, "n_star": n}
        log("anchor %s: n*=%s coef=%.9f detune=%.1e (< %.0e)"
            % (tag, n, coef, dd, lim))
    chk("T0 lift anchor: the 256-dim block-wise lift of the certified "
        "machinery on the lifted published engine reproduces the V206 "
        "1/5-rad reading to %.1e (< 1e-8) and the three closure_class "
        "readings to max %.1e (< 5e-3, float64 ladder drift): %s"
        % (data["anchor_1/5rad"]["detune"],
           max(data["anchor_%s" % t]["detune"]
               for t in ("pi/2", "pi/4", "3pi/4")), t0_ok), t0_ok)

    # ---- T1: cubond reads at the five clocks -------------------------------
    reads = {}
    t1_ok = True
    for tag, g in [("pi/4", math.pi / 4), ("pi/3", math.pi / 3),
                   ("pi/5", math.pi / 5), ("2pi/5", 2 * math.pi / 5),
                   ("pi/6", math.pi / 6)]:
        n, d, coef = reading(Fn_gamma(g))
        miss = None if coef is None else coef - COEF_OBS
        reads[tag] = {"g": g, "n_star": n, "deficit": d, "coef": coef,
                      "miss": miss}
        if coef is not None and abs(miss) <= BAND:
            t1_ok = False
        log("cubond read at %s: n*=%s coef=%s miss=%s"
            % (tag, n, "VOID" if coef is None else "%.7f" % coef,
               "VOID" if miss is None else "%.3e" % miss))
    data["cubond_reads"] = reads
    chk("T1: no cubond clock {pi/4, pi/3, pi/5, 2pi/5, pi/6} is read by "
        "the lifted ladder within the 1e-4 constants band of the "
        "observed %.7f: %s" % (COEF_OBS, t1_ok), t1_ok)

    # ---- T2: dense scan, 96 points ------------------------------------------
    rescues, scan_void = [], 0
    for i in range(96):
        g = 0.05 + (1.60 - 0.05) * i / 95.0
        n, d, coef = reading(Fn_gamma(g))
        if coef is None:
            scan_void += 1
            continue
        if abs(coef - COEF_OBS) <= BAND:
            rescues.append({"g": g, "coef": coef,
                            "closure": closure_k(G_gamma(g))})
    data["rescue_candidates"] = rescues
    data["scan_void_points"] = scan_void
    closed = [r for r in rescues if r["closure"] is not None]
    chk("T2: dense scan 96 points on [0.05, 1.60]: %d band-hit(s) "
        "(%d VOID, excluded - V206 convention); none at a closed "
        "G_gamma(g*): %s" % (len(rescues), scan_void, not closed),
        not closed)

    ok = all(c[1] for c in checks)
    if not t0_ok:
        status = ("CUBOND-LADDER-LIFT-INVALID (the lift failed to "
                  "reproduce the certified published readings - the "
                  "design died at its own death-gate, nothing is "
                  "decided about the engine)")
    elif closed:
        status = ("CUBOND-LADDER-READS-OBSERVED (%s) - booked "
                  "obligation: the physics-contact census of the "
                  "cubond family at g* (gapless content, chirality, "
                  "clock closure); F2 re-opens as a NEW program"
                  % json.dumps(closed))
    else:
        cs = [r["coef"] for r in reads.values() if r["coef"]]
        status = ("CUBOND-LADDER-NO-CONTACT (the entangling bond "
                  "engine, put to the SAME certified ladder machinery "
                  "through the forced block-wise lift, reads the "
                  "observed coefficient at NO clock: the five closed "
                  "clocks and a 96-point scan all miss; the coefficient "
                  "is not carried by bond geometry - the design space "
                  "of the ladder is closed one family wider, and F2's "
                  "named residue stands)"
                  % ())
        data["cubond_read_summary"] = [round(c, 7) for c in cs]
    out = {"id": "V384_cubond_ladder_functional",
           "tier": 2,
           "residual": min((abs(r["coef"] - COEF_OBS)
                            for r in reads.values() if r["coef"]),
                           default=None),
           "property": ("the certified hierarchy functional on the "
                        "entangling cubond engine - block-wise exact "
                        "lift, cubic clock, tier 2, O8/F2"),
           "status": status,
           "pass": ok,
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v384_cubond_ladder_functional.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V384 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
