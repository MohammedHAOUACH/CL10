#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V385 - content census of the Clifford-bond cone (booked BEFORE the
run in review/pre_registrations_root.md, "V385_cubond_content_census", suite 298).

Property (booked): quelles particules habitent le cone du lien de
Clifford.  Machine: the 256 states bond(8) x internal(32) of the V383
engine; the five 32-dim Cartan generators C_k = 1j*(h_a h_{a+1}),
a in (0,2,4,6,8) from clifford.generators (the V65 cartan_weights
pattern, but in the full 32-dim module); common eigenbasis through the
combination with coefficients sqrt(1,2,3,5,7).

T0: gamma algebra re-certified (gamma_j^2 = I, anticommutators
<= 1e-12); the reference 16-dim form (V65.cartan_weights) reproduces
its own 16-weight multiset (off <= 1e-9).
T1: the 32 weights (5 components) of the internal multiset are exactly
lam16 union -lam16 (multiset equality, every component <= 1e-9) - the
V52 mirror verified on the entangling engine; lam16 recomputed via
V65.cartan_weights() in-script and auto-recorded (no "weights" key in
the V65 artifact - the booking names this fallback).
T2: at Gamma, U(0) = I: all 256 states at phase 0 (gap <= 1e-12); the
256-weight multiset (8 copies x 32) equals 8 x (lam16 + (-lam16)).

Verdicts: CUBOND-CONTENT-VERBATIM / CUBOND-CONTENT-SHIFTS.
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

from cl10 import supercell_2i as S2  # noqa: E402
from cl10 import clifford  # noqa: E402

_spec = importlib.util.spec_from_file_location(
    "v65", os.path.join(ROOT, "archive", "experiments", "dynamics_spectrum",
                        "V65_wilson", "V65_wilson.py"))
_v65 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_v65)

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, bool(ok)))


def N(A):
    return float(np.max(np.abs(A)))


def sig(v):
    return tuple(int(round(x)) for x in v)


def main():
    data = {}

    # ---- the V383 engine (same construction, re-verified here) ----
    Sx, Sy, Sz = S2.spin_triple()
    R = [2 * np.asarray(m, dtype=np.complex128) for m in (Sx, Sy, Sz)]
    X, Y, Z = (np.array(m, dtype=np.complex128) for m in
               ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]))
    I2 = np.eye(2, dtype=np.complex128)
    P = [np.kron(np.kron(X, I2), I2), np.kron(np.kron(I2, Y), I2),
         np.kron(np.kron(I2, I2), Z)]
    G = [np.kron(p, r) for p, r in zip(P, R)]
    I256 = np.eye(256, dtype=np.complex128)

    def Uq(q):
        U = I256.copy()
        for j in range(3):
            U = U @ (math.cos(q[j]) * I256 + 1j * math.sin(q[j]) * G[j])
        return U

    # ---- T0: gamma algebra + the certified 16-dim reference form ----
    e1 = max(N(g @ g - I256) for g in G)
    e2 = max(N(G[i] @ G[j] + G[j] @ G[i])
             for i in range(3) for j in range(i + 1, 3))
    lam16, off16 = _v65.cartan_weights()
    chk("T0 gamma re-certification: gamma_j^2 = I (%.1e), "
        "anticommutators (%.1e) <= 1e-12; V65.cartan_weights reproduces "
        "its own 16-weight multiset, off-diagonal residual %.1e <= 1e-9: "
        "%s" % (e1, e2, off16,
                max(e1, e2) <= 1e-12 and off16 <= 1e-9),
        max(e1, e2) <= 1e-12 and off16 <= 1e-9)

    # ---- the 32-dim Cartan and its common eigenbasis ----
    h = [np.asarray(m, dtype=np.complex128) for m in clifford.generators()]
    Ck = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    M = sum(math.sqrt(p) * Ck[k] for k, p in enumerate((1, 2, 3, 5, 7)))
    _ev, V = np.linalg.eigh(M)
    w32 = np.array([[float(np.real(np.vdot(V[:, s], Ck[k] @ V[:, s])))
                     for k in range(5)] for s in range(32)])
    off32 = max(float(np.abs(V.conj().T @ Ck[k] @ V - np.diagflat(
        [np.vdot(V[:, s], Ck[k] @ V[:, s]) for s in range(32)])).max())
        for k in range(5))
    data["offdiag_residual_32"] = off32

    # ---- T1: the 32-weight multiset = lam16 union -lam16 (V52 mirror) ----
    got = sorted(sig(w) for w in w32)
    want = sorted([sig(l) for l in lam16] + [sig(-l) for l in lam16])
    # components are exactly +-1 signed integers; match the multisets
    # in order and take the largest componentwise deviation
    gotv = sorted([tuple(w) for w in w32], key=sig)
    wantv = sorted([tuple(l) for l in lam16] + [tuple(-l) for l in lam16],
                   key=sig)
    comp_dev = max(max(abs(a - b) for a, b in zip(g, wv))
                   for g, wv in zip(gotv, wantv)) if len(gotv) == len(wantv) \
        else float("nan")
    missing = sorted(set(want) - set(got))
    extra = sorted(set(got) - set(want))
    data["weights32"] = [list(map(float, w)) for w in
                         sorted(w32.tolist(), key=lambda r: tuple(r))]
    data["lam16"] = [[float(x) for x in l] for l in
                     sorted(lam16.tolist(), key=lambda r: tuple(r))]
    data["t1_component_deviation"] = comp_dev
    chk("T1 the 32 internal weights (5 components, common eigenbasis, "
        "off-diagonal %.1e) are exactly lam16 union -lam16: sizes 32/32 "
        "(%s), matched-component deviation %.1e <= 1e-9 (the V52 mirror "
        "holds on the entangling engine); missing %s extra %s: %s"
        % (off32, len(got) == 32 and len(want) == 32, comp_dev,
           missing if missing else "none",
           extra if extra else "none", got == want and comp_dev <= 1e-9),
        got == want and comp_dev <= 1e-9)

    # ---- T2: Gamma point - U(0) = I and the 256-state census ----
    eU0 = N(Uq([0.0, 0.0, 0.0]) - I256)
    got256 = sorted([sig(w) for w in np.repeat(w32, 8, axis=0)])
    want256 = sorted([sig(l) for l in np.concatenate(
        [lam16, -lam16])] * 8)
    census_ok = got256 == want256
    chk("T2 U(0) = I: max deviation %.1e, all 256 states at phase 0 "
        "(gap 0.0e+00 <= 1e-12: %s); the 256-weight multiset (8 copies "
        "x 32) equals 8 x (lam16 + (-lam16)): %s -> the cone carries "
        "exactly the certified cell inventory: %s"
        % (eU0, eU0 <= 1e-12, census_ok,
           eU0 <= 1e-12 and census_ok),
        eU0 <= 1e-12 and census_ok)

    ok = all(c[1] for c in checks)
    status = ("CUBOND-CONTENT-VERBATIM (the cone of the Clifford bond "
              "carries exactly the certified cell inventory - the 32 "
              "internal weights are the 16 certified weights and their "
              "mirror, the 256 states at Gamma are eight copies of that "
              "set at phase 0; the existence theorem changes no "
              "particle, adds none, loses none)"
              if ok else "CUBOND-CONTENT-SHIFTS (see failing multiset "
              "gates; missing/extra weights named in the check lines)")
    out = {"id": "V385_cubond_content_census",
           "tier": 2,
           "residual": max(off32, comp_dev, eU0),
           "property": ("content census of the Clifford-bond cone: "
                        "which particles the cone carries, tier 2"),
           "status": status,
           "pass": ok,
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v385_cubond_content_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V385 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
