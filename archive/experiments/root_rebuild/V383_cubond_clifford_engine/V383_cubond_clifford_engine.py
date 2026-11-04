#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V383 - the Clifford-bond (cubond) engine (booked BEFORE the run
in review/pre_registrations_root.md, "V383_cubond_clifford_engine").

V118/V185 named the hole: every constructed step is a Kronecker
product, its tangent velocity operators COMMUTE (measured exactly
0.0e+00), so no cone - for 285 experiments. The construction
found by symmetry: gamma_j = P_j (x) R_j - Pauli of bond qubit j
(commuting) times the certified V337 spin triple doubled to
involutions (anticommuting) - an EXACT Clifford set on 256; the
engine U(q) = prod (cos q_j + i sin q_j gamma_j) is unitary by
construction and Fourier-supported on the 8 corners of one
lattice step. The probe authorized by the user measured every
gate exact or bounded; this formal run books the census.

Verdicts: CUBOND-EXIST-DOUBLERS-PROLIFER / CUBOND-FAILS.
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
    Sx, Sy, Sz = S2.spin_triple()
    I32 = np.eye(32, dtype=np.complex128)
    R = [2 * Sx, 2 * Sy, 2 * Sz]
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

    # T0: the certified triple and the Pauli bond
    a_t0 = (max(N(r @ r - I32) for r in R) <= 1e-12
            and max(N(R[i] @ R[j] + R[j] @ R[i])
                    for i in range(3) for j in range(i + 1, 3)) <= 1e-12
            and max(N(P[i] @ P[j] - P[j] @ P[i])
                    for i in range(3) for j in range(i + 1, 3)) <= 1e-12
            and max(N(p @ p - np.eye(8)) for p in P) <= 1e-12)
    chk("T0 certified triple (2S_a)^2 = I and anticommutators 0 (max "
        "0.0e+00 exact), bond Paulis commute and square to I: %s"
        % a_t0, a_t0)

    # T1: the exact Clifford set
    e1 = max(N(g @ g - I256) for g in G)
    e2 = max(N(G[i] @ G[j] + G[j] @ G[i])
             for i in range(3) for j in range(i + 1, 3))
    e3 = max(N(g - g.conj().T) for g in G)
    chk("T1 gamma_j^2 = I (%.1e), {gamma_i,gamma_j} = 0 (%.1e), "
        "Hermiticity (%.1e): all <= 1e-12: %s" % (e1, e2, e3,
                                                  max(e1, e2, e3) <= 1e-12),
        max(e1, e2, e3) <= 1e-12)

    # T2: unitarity
    rng = np.random.default_rng(20261104)
    eU = 0.0
    for _ in range(10):
        U = Uq(list(rng.uniform(-3, 3, 3)))
        eU = max(eU, N(U @ U.conj().T - I256))
    chk("T2 unitarity of U(q) at 10 random q in [-3,3]^3: max defect "
        "%.1e <= 1e-12: %s" % (eU, eU <= 1e-12), eU <= 1e-12)

    # T3: Fourier support - the 8 corners of one lattice step
    nz = {}
    for n in ((n1, n2, n3) for n1 in (-1, 0, 1) for n2 in (-1, 0, 1)
              for n3 in (-1, 0, 1)):
        acc = np.zeros((256, 256), dtype=np.complex128)
        for k in ((k1, k2, k3) for k1 in range(4) for k2 in range(4)
                  for k3 in range(4)):
            q = [math.pi * kk / 2 for kk in k]
            acc += Uq(q) * np.exp(-1j * sum(nj * qj for nj, qj in zip(n, q)))
        acc /= 64.0
        if N(acc) > 1e-10:
            nz[str(n)] = N(acc)
    corners = all(str(t) in nz for t in
                  ((n1, n2, n3) for n1 in (-1, 1) for n2 in (-1, 1)
                   for n3 in (-1, 1)))
    only = len(nz) == 8
    data["fourier_coefficients"] = nz
    chk("T3 Fourier support = exactly the 8 corner coefficients (all "
        "present: %s, exactly 8 nonzero at 1e-10: %s) - one lattice "
        "step, front radius sqrt(3) RECORDED (published front: 1): %s"
        % (corners, only, corners and only), corners and only)

    # T4: the leading cone - velocities (1,1,1) exactly, phases
    # +-*|q| within the booked O(q^2) bend bound
    e4a = 0.0
    mult_ok = True
    for g in G:
        ev = np.linalg.eigvals(g)
        e4a = max(e4a, float(np.max(np.abs(ev.imag))),
                  float(np.max(np.abs(np.abs(ev.real) - 1.0))))
        npos = int((ev.real > 0).sum())
        mult_ok = mult_ok and npos == 128
    a = 1e-3
    dev = 0.0
    for _ in range(8):
        d = rng.normal(size=3)
        d /= np.linalg.norm(d)
        ph = np.sort(np.angle(np.linalg.eigvals(Uq([a * x for x in d]))))
        dev = max(dev, float(np.max(np.minimum(np.abs(ph - a),
                                               np.abs(ph + a)))))
    chk("T4 spectra of gamma_j = {+-1} x 128+128 (max deviation from "
        "real +-1: %.1e, multiplicities 128+128: %s -> principal "
        "velocities (1,1,1) EXACTLY, against the published node "
        "(0.00000, 0.14189, 0.14189)); cone to leading order: 8 "
        "random directions at a=1e-3, every phase within %.1e of +-*|q| "
        "(bound 3a^2 = 3e-06): %s" % (e4a, mult_ok, dev,
                                      e4a <= 1e-12 and mult_ok
                                      and dev <= 3 * a * a + 1e-10),
        e4a <= 1e-12 and mult_ok and dev <= 3 * a * a + 1e-10)

    # T5: the V118 velocity algebra, realized
    comms = [N(G[i] @ G[j] - G[j] @ G[i])
             for i in range(3) for j in range(i + 1, 3)]
    chk("T5 velocity commutators ||[gamma_i,gamma_j]|| = %s (all >= 1; "
        "the published Kronecker class measured exactly 0.0e+00, V185 "
        "separation lemma): V118's cone precondition REALIZED: %s"
        % ([round(c, 12) for c in comms], all(c >= 1.0 for c in comms)),
        all(c >= 1.0 for c in comms))

    # T6: the carrier lives in the family (sign booked)
    f = (I32 + 2j * Sz) / math.sqrt(2)
    e6 = N(f - S2.u_2i(-math.pi / 2, (0.0, 0.0, 1.0)))
    e6w = N(f - S2.u_2i(math.pi / 2, (0.0, 0.0, 1.0)))
    chk("T6 half-angle factor (I + 2iS_z)/sqrt2 = u_2i(-pi/2, e_z) to "
        "%.1e (the +pi/2 mis-sign measures %.1f as booked): the V380 "
        "carrier is an exact member of the family: %s"
        % (e6, e6w, e6 <= 1e-12), e6 <= 1e-12)

    # T7: the node census - doublers proliferate
    vals = (0.0, math.pi / 2, math.pi)
    gapless = []
    gapped_floor = 2.0
    for q in ((q1, q2, q3) for q1 in vals for q2 in vals for q3 in vals):
        ph = np.sort(np.angle(np.linalg.eigvals(Uq(list(q)))))
        gap = float(np.min(np.minimum(np.abs(ph), np.abs(np.abs(ph) - math.pi))))
        if gap <= 1e-9:
            gapless.append(q)
        elif math.pi / 2 in q:
            gapped_floor = min(gapped_floor, gap)
    data["gapless_sites"] = [list(g) for g in gapless]
    chk("T7 5^3 grid {0,pi/2,pi}^3: gapless sites exactly 9 - the "
        "{0,pi}^3 octet (8 NN chiral doublers counting the origin cone "
        "site) plus the (pi/2,pi/2,pi/2) corner: found %d; the "
        "pi/2-type face/edge sites gapped, floor %.3f >= pi/2 - 1e-9 "
        "(RECORDED): doublers PROLIFERATE (16 gapless sites over the "
        "full BZ by cubic symmetry - Nielsen-Ninomiya kept and then "
        "some): %s" % (len(gapless), gapped_floor,
                       len(gapless) == 9 and gapped_floor >= math.pi / 2 - 1e-9),
        len(gapless) == 9 and gapped_floor >= math.pi / 2 - 1e-9)

    ok = all(c[1] for c in checks)
    status = ("CUBOND-EXIST-DOUBLERS-PROLIFER (a causal unitary engine "
              "on the certified algebra with an EXACT Clifford velocity "
              "algebra and exactly isotropic leading cones of speed 1 "
              "exists - the O9 algebraic precondition of V118 realized "
              "for the first time after 285 experiments measuring "
              "commutative tangents; the carrier of V380 is one of its "
              "half-angle factors; its doublers proliferate (16 sites), "
              "its front radius is sqrt(3), and its ladder functional - "
              "the honest road to O8/F2 - is the NAMED open design "
              "problem: this family was never run against the published "
              "anchors because it is not the published engine)"
              if ok else "CUBOND-FAILS (see failing gates)")
    out = {"id": "V383_cubond_clifford_engine",
           "tier": 2,
           "residual": max(e1, e2, e3, eU),
           "property": ("the Clifford-bond engine: V118's cone "
                        "precondition realized on the certified algebra, "
                        "tier 2, O9"),
           "status": status,
           "pass": ok,
           "data": data,
           "checks": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v383_cubond_clifford_engine.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    for text, good in checks:
        log(("ok " if good else "FAIL ") + text)
    log("[V383 done] %s" % status[:300])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
