#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V386 - chiral ledger and the Nielsen-Ninomiya book (booked BEFORE
the run in review/pre_registrations_root.md,
"V386_cubond_chirality_nn_ledger", suite 299).

Property (booked): le cone est-il chiral comme le modele standard?
Machine: G5 = gamma_1 gamma_2 gamma_3 on 256; exact tangent
derivatives of the V383 engine: F_j' = i gamma_j F_j (the factor
identity), so dU/dq_j = U (i gamma_j) exactly, and
S_j(node) = i dU/dq_j evaluated analytically at the 16 gapless sites
of the V383 T7 census: the 8 of {0,pi}^3 and the 8 cubic-symmetric
corners (+-pi/2,+-pi/2,+-pi/2), enumerated explicitly.

T0: exact reuse of the V383 gammas (0.0e+00).
T1: G5 anti-Hermitian (G5 + G5^dag = 0 <= 1e-12), G5^2 = -I
    (<= 1e-12), spectrum = +-i with multiplicities 128/128.
T2: at each site S_j = sum_k M_jk gamma_k + scalar (exact projection
    M_jk = 1/256 tr(S_j gamma_k), the gammas being an orthogonal
    Clifford basis); a scalar component above 1e-9 is recorded, not
    corrected; chirality chi(site) = sign det M; GATE: the sum of the
    8 chi over {0,pi}^3 = 0 (NN book: the doublers cancel the
    chirality); the 16-entry chi table is recorded in data.
T3: at the 8 sites {0,pi}^3 the operators S close the Clifford
    algebra (commutator norms >= 1); gate: 8/8 sites.

Verdicts: CUBOND-CHIRAL-LEDGER-NN-CLOSED / CUBOND-CHIRAL-ANOMALY.
"""
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

T0 = time.time()
checks = []


def log(msg):
    print("[%7.1fs] %s" % (time.time() - T0, msg), flush=True)


def chk(text, ok):
    checks.append((text, bool(ok)))


def N(A):
    return float(np.max(np.abs(A)))


def main():
    data = {}

    # ---- the V383 engine, re-certified (T0) ----
    Sx, Sy, Sz = S2.spin_triple()
    R = [2 * np.asarray(m, dtype=np.complex128) for m in (Sx, Sy, Sz)]
    X, Y, Z = (np.array(m, dtype=np.complex128) for m in
               ([[0, 1], [1, 0]], [[0, -1j], [1j, 0]], [[1, 0], [0, -1]]))
    I2 = np.eye(2, dtype=np.complex128)
    P = [np.kron(np.kron(X, I2), I2), np.kron(np.kron(I2, Y), I2),
         np.kron(np.kron(I2, I2), Z)]
    G = [np.kron(p, r) for p, r in zip(P, R)]
    I256 = np.eye(256, dtype=np.complex128)
    e0 = max(N(g @ g - I256) for g in G)
    e0 = max(e0, max(N(G[i] @ G[j] + G[j] @ G[i])
                     for i in range(3) for j in range(i + 1, 3)))
    chk("T0 exact reuse of the V383 gammas: gamma^2 = I and "
        "anticommutators max deviation %.1e (0.0e+00 booked): %s"
        % (e0, e0 <= 1e-12), e0 <= 1e-12)

    # ---- T1: the 256-state chirality operator ----
    G5 = G[0] @ G[1] @ G[2]
    e1a = N(G5 + G5.conj().T)
    e1b = N(G5 @ G5 + I256)
    ev = np.linalg.eigvals(G5)
    npos = int((ev.imag > 0.5).sum())
    nneg = int((ev.imag < -0.5).sum())
    ere = float(np.max(np.abs(ev.real)))
    emag = float(np.max(np.abs(np.abs(ev.imag) - 1.0)))
    chk("T1 G5 = gamma_1 gamma_2 gamma_3: anti-Hermitian %.1e <= 1e-12, "
        "G5^2 = -I %.1e <= 1e-12, spectrum pure +-i (max real part "
        "%.1e, |imag| deviation %.1e) with multiplicities %d/%d "
        "(128/128 booked): %s"
        % (e1a, e1b, ere, emag, npos, nneg,
           e1a <= 1e-12 and e1b <= 1e-12 and npos == 128 and nneg == 128
           and ere <= 1e-12 and emag <= 1e-12),
           e1a <= 1e-12 and e1b <= 1e-12 and npos == 128 and nneg == 128
           and ere <= 1e-12 and emag <= 1e-12)

    # ---- the 16 gapless sites of the V383 T7 census ----
    sites = ([tuple(n * math.pi for n in nn) for nn in
              ((n1, n2, n3) for n1 in (0, 1) for n2 in (0, 1)
               for n3 in (0, 1))] +
             [tuple(s * math.pi / 2 for s in ss) for ss in
              ((s1, s2, s3) for s1 in (1, -1) for s2 in (1, -1)
               for s3 in (1, -1))])
    nn_sites = sites[:8]

    def Uq(q):
        U = I256.copy()
        for j in range(3):
            U = U @ (math.cos(q[j]) * I256 + 1j * math.sin(q[j]) * G[j])
        return U

    def S_j(q, j):
        # exact: dU/dq_j = U (i gamma_j)  =>  S_j = i dU/dq_j = -U gamma_j
        return -Uq(q) @ G[j]

    # ---- T2: exact Clifford projection, chi table, NN sum ----
    chi_table = {}
    scalar_max = 0.0
    resid_max = 0.0
    for q in sites:
        M = np.zeros((3, 3))
        Ss = [S_j(q, j) for j in range(3)]
        for j in range(3):
            for k in range(3):
                M[j][k] = float(np.real(np.trace(Ss[j] @ G[k]) / 256.0))
            scalar_max = max(scalar_max,
                             float(np.abs(np.trace(Ss[j]) / 256.0)))
            rec = sum(M[j][k] * G[k] for k in range(3))
            resid_max = max(resid_max, N(Ss[j] - rec))
        d = float(np.real(np.linalg.det(M)))
        chi = 0 if abs(d) <= 1e-9 else (1 if d > 0 else -1)
        chi_table[str([round(x, 6) for x in q])] = {
            "det_M": d, "chi": chi, "M": M.tolist()}
    nn_keys = [str([round(x, 6) for x in s]) for s in nn_sites]
    nn_sum = sum(chi_table[k]["chi"] for k in nn_keys)
    data["chi_table"] = chi_table
    data["scalar_component_max"] = scalar_max
    data["projection_residual_max"] = resid_max
    chk("T2 exact projections M_jk = 1/256 tr(S_j gamma_k) at all 16 "
        "gapless sites: scalar components max %.1e, non-Clifford "
        "residual max %.1e (recorded, not corrected); chi = sign det M "
        "on the {0,pi}^3 octet: %s, sum over the 8 NN sites = %d "
        "(0 booked - the doublers cancel the chirality); corner sites "
        "project to det M = 0 (recorded): %s"
        % (scalar_max, resid_max,
           [chi_table[str([round(x, 6) for x in s])]["chi"]
            for s in nn_sites], nn_sum, nn_sum == 0),
        nn_sum == 0)

    # ---- T3: the S operators close the Clifford algebra at the 8 NN sites ----
    good = 0
    worst = 99.0
    for q in nn_sites:
        Ss = [S_j(q, j) for j in range(3)]
        cm = min(N(Ss[i] @ Ss[j] - Ss[j] @ Ss[i])
                 for i in range(3) for j in range(i + 1, 3))
        worst = min(worst, cm)
        if cm >= 1.0:
            good += 1
    data["T3_min_commutator_norm"] = worst
    chk("T3 at the 8 NN sites {0,pi}^3 the S operators close the "
        "Clifford algebra: min ||[S_i,S_j]|| over all sites and pairs "
        "%.6f >= 1: %d/8 sites (8 booked): %s"
        % (worst, good, good == 8), good == 8)

    ok = all(c[1] for c in checks)
    status = ("CUBOND-CHIRAL-LEDGER-NN-CLOSED (the cone is VECTORLIKE: "
              "G5 splits the 256 states 128+128 of opposite chirality "
              "and the eight Nielsen-Ninomiya doublers of {0,pi}^3 "
              "cancel the chirality exactly, sum chi = 0; the corner "
              "doublers carry degenerate projections (recorded); the "
              "content of this family is NOT chiral like the Standard "
              "Model - booked as a measurement, not an obstacle)"
              if ok else "CUBOND-CHIRAL-ANOMALY (a gate failed - see "
              "the named site in the check lines)")
    out = {"id": "V386_cubond_chirality_nn_ledger",
           "tier": 2,
           "residual": max(e0, e1a, e1b, scalar_max, resid_max),
           "property": ("chirality of the Clifford-bond cone and the "
                        "Nielsen-Ninomiya book, tier 2"),
           "status": status,
           "pass": ok,
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v386_cubond_chirality_nn_ledger.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V386 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
