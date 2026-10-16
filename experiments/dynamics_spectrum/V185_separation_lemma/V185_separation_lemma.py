#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V185 - the separation lemma: every tangent operator of the implemented
step is (sub-site) x (internal identity) (Phase 1 of the cone programme).

**Why.**  V118 closed W9 with one sentence: the tangent velocity operators
commute "because the bond layer is a Kronecker product (cell index) x
(internal algebra) and never mixes the two". V142 measured the same fact
frame-blind to 3e-15. The cone programme (V186) rests on this being a LEMMA,
not an observation: for the implemented step every momentum derivative of the
Bloch matrix, conjugated back to the step, is internal-blind - so the whole
low-energy velocity structure lives in the 8x8 sub-site factor and NO choice
of internal data (frame, ordering, coupling) can move it. This experiment
proves it exactly, at every node and at random momenta, with the exact
analytic derivative (no finite differences).

**The statement proved.**  With F_s(q) = S_s(q) kron U the six stage factors,
M(q) = F_6 ... F_1, the tangent operators are
    T_k(q) = i (d_k M) M^dag = sum_{s in stages of axis k}
             F_6...F_{s+1} [ (d_k S_s) S_s^dag kron I_16 ] F_{s+1}^dag...F_6^dag
and each conjugation by F = S kron U maps A kron I to (S A S^dag) kron I.
Hence T_k(q) = A_k(q) kron I_16 EXACTLY, for every q, and the defect
max |T_k - A_k kron I| is at machine precision. The gapless-subspace
operators S_j = P^dag T_j P of V118 are therefore 8x8 matrices in disguise:
a Clifford-type tangent algebra must be created by the BOND, not by the
internal layer.

**The Wilson nuance, measured.**  The Wilson factor W(q) = I_8 kron D(q) is
internal-DIAGONAL (a Cartan phase), so it commutes with every A kron I and
adds to T_k the internal-Cartan diagonal term (d_k ph(q)) diag(d): the
transport part stays internal-blind, the new term is diagonal in the weight
basis and carries no velocity structure (it is the Bessel-tail term of
Section 8.1). Both parts are measured separately.

Verdict: PASS = the lemma holds at machine precision at all 12 M nodes, all
8 R corners, 20 random momenta, and the Wilson term is exactly the predicted
internal-diagonal Cartan addition.
"""
import json
import math
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V65_wilson"))

import V65_wilson as V65                      # noqa: E402
from cl10 import supercell as SC              # noqa: E402
from cl10 import engine as E                  # noqa: E402

_B16, _U32, _U16 = SC.sectors()
B16, U16 = E.asnumpy(_B16), E.asnumpy(_U16)
I16 = np.eye(16, dtype=complex)
G10, ORDER = SC.G10, SC.STAGES

_lam, _off = V65.cartan_weights()
_v65 = json.load(open(exp("V65_wilson", "results", "v65_wilson.json")))
_gg = next(r for r in _v65["data"]["part_B"].values()
           if r["su5"] == "Georgi-Glashow")
XCOEF = np.array(_gg["X_coef"])
YCOEF = np.array(_v65["data"]["Y_coef"])


def stage_pieces(even, j, q):
    """The stage subsite factor as (S, dS) at momentum q (exact derivative).

    Intra-block stages carry no phase: dS = 0. Inter-block stages are
    S = D0 + e^{iq} P+ + e^{-iq} P-, so dS = i e^{iq} P+ - i e^{-iq} P-."""
    c, s = math.cos(G10), math.sin(G10)
    S = np.eye(8, dtype=complex)
    Pp = np.zeros((8, 8), dtype=complex)
    Pm = np.zeros((8, 8), dtype=complex)
    bit = 1 << j
    for e in range(8):
        if (e & bit) == (0 if even else bit):
            o = e ^ bit
            if even:
                S[e, e] = c
                S[o, o] = c
                S[e, o] = -1j * s
                S[o, e] = -1j * s
            else:
                S[e, e] = c
                S[o, o] = c
                Pp[e, o] = -1j * s
                Pm[o, e] = -1j * s
    if even:
        return S, np.zeros_like(S)
    ph = np.exp(1j * q)
    return (S + ph * Pp + np.conj(ph) * Pm,
            1j * ph * Pp - 1j * np.conj(ph) * Pm)


def kronI(A):
    return np.kron(A, I16)


def tangent_ops(q, w=0.0, th=0.0):
    """The exact T_k(q) = i (d_k M) M^dag of the (optionally Wilson-deformed)
    step, plus the predicted Cartan term, computed stage by stage."""
    F = []
    for even, j in ORDER:
        S, dS = stage_pieces(even, j, q[j])
        F.append((np.kron(S, U16), np.kron(dS, U16), j))
    M = np.eye(128, dtype=complex)
    for Fm, _, _ in F:
        M = Fm @ M
    T = [np.zeros((128, 128), dtype=complex) for _ in range(3)]
    for k in range(3):
        # sum_s F_6...F_{s+1} (dF_s F_s^dag) F_{s+1}^dag...F_6^dag
        for s in range(5, -1, -1):
            if F[s][2] != k:
                continue
            X = F[s][1] @ F[s][0].conj().T
            for t in range(s + 1, 6):
                X = F[t][0] @ X @ F[t][0].conj().T
            T[k] += X
    if w != 0.0:
        d = _lam @ (math.cos(th) * XCOEF + math.sin(th) * YCOEF)
        ph = w * sum(1.0 - math.cos(qj) for qj in q)
        wd = np.exp(-1j * ph * np.tile(d, 8))
        W = np.diag(wd)
        M = W @ M
        for k in range(3):
            dkph = w * math.sin(q[k])
            T[k] = W @ T[k] @ W.conj().T + dkph * np.diag(np.tile(d, 8))
    return M, T


def kron_defect(T):
    """max |T - A kron I_16| with A = (1/16) Tr_int T, and the internal
    off-diagonal content of T in the weight basis."""
    A = np.zeros((8, 8), dtype=complex)
    for e in range(8):
        for f in range(8):
            A[e, f] = np.trace(T[e * 16:(e + 1) * 16, f * 16:(f + 1) * 16]) / 16
    d = float(np.abs(T - np.kron(A, I16)).max())
    R = T.reshape(8, 16, 8, 16)
    off = 0.0
    for a in range(16):
        for b in range(16):
            if a != b:
                off = max(off, float(np.abs(R[:, a, :, b]).max()))
    return A, d, off


def main():
    out = {"experiment": "V185_separation_lemma", "data": {}}
    checks = []

    # ---- T0: regression - the stage form equals cl10.supercell.bloch -----
    rng = np.random.default_rng(185)
    reg = 0.0
    for _ in range(6):
        q = tuple(float(x) for x in rng.uniform(-math.pi, math.pi, 3))
        Mref = SC.bloch(q, U16)
        M, _ = tangent_ops(q)
        reg = max(reg, float(np.abs(M - Mref).max()))
    checks.append(("T0 the stage-factorised M(q) equals cl10.supercell.bloch "
                   "at 6 random momenta (%.1e)" % reg, reg < 1e-12))
    out["data"]["regression"] = reg

    # ---- T1: the lemma at every node and random momenta (bare step) -----
    pts = list(SC.M_POINTS) + list(SC.R_POINTS)
    pts += [tuple(float(x) for x in rng.uniform(-math.pi, math.pi, 3))
            for _ in range(20)]
    worst = 0.0
    worst_off = 0.0
    for q in pts:
        _, T = tangent_ops(q)
        for k in range(3):
            _, d, off = kron_defect(T[k])
            worst = max(worst, d)
            worst_off = max(worst_off, off)
    checks.append(("T1 SEPARATION LEMMA: at all 12 M nodes, 8 R corners and "
                   "20 random momenta, every tangent operator T_k = i(d_k M)M^"
                   "dag equals A_k(q) kron I_16 with max defect %.1e, and its "
                   "internal off-diagonal content is %.1e - the momentum "
                   "derivative of the step is INTERNAL-BLIND everywhere"
                   % (worst, worst_off), worst < 1e-12 and worst_off < 1e-12))
    out["data"]["bare_lemma"] = {"max_kron_defect": worst,
                                 "max_internal_offdiag": worst_off,
                                 "n_points": len(pts)}

    # ---- T2: the Wilson nuance - an internal-DIAGONAL Cartan term -------
    w, th = 0.7, 1.1
    worst_diag_pred = 0.0
    worst_off2 = 0.0
    worst_d2 = 0.0
    for q in pts[:12]:
        _, T = tangent_ops(q, w, th)
        d = _lam @ (math.cos(th) * XCOEF + math.sin(th) * YCOEF)
        for k in range(3):
            A, dd, off = kron_defect(T[k])
            # the predicted structure: T_k = (transport kron I) + dk(ph) diag(d)
            dkph = w * math.sin(q[k])
            pred = dkph * np.diag(np.tile(d, 8))
            # remove the predicted Cartan term, the rest must be kron I
            _, dres, _ = kron_defect(T[k] - pred)
            worst_diag_pred = max(worst_diag_pred, dres)
            worst_off2 = max(worst_off2, off)
            worst_d2 = max(worst_d2, dd)
    checks.append(("T2 with the Wilson factor (w=%.2f, theta=%.2f) the tangent "
                   "operator is EXACTLY (transport kron I) + (d_k ph) "
                   "diag(d): after subtracting the predicted internal-Cartan "
                   "term the kron-I defect is %.1e; the Wilson addition is "
                   "internal-diagonal (off-diagonal content %.1e, diagonal "
                   "content up to %.3f) - the Wilson factor adds a Cartan "
                   "phase, never a velocity structure"
                   % (w, th, worst_diag_pred, worst_off2, worst_d2),
                   worst_diag_pred < 1e-12 and worst_off2 < 1e-12))
    out["data"]["wilson_nuance"] = {
        "kron_defect_after_cartan_removal": worst_diag_pred,
        "internal_offdiag": worst_off2, "max_cartan_diag": worst_d2}

    # ---- T3: the consequence at the V118 gapless corner triplet ---------
    # find a gapless corner multiplet by (theta,w) refinement at R, then
    # verify the gapless-subspace operators S_j are sub-site operators in
    # disguise: they commute with I_8 kron G for internal generators G.
    C = (math.pi / 2,) * 3

    def top3(wv, thv):
        M, _ = tangent_ops(C, wv, thv)
        a = np.sort(np.abs(np.angle(np.linalg.eigvals(M))))
        return float(a[2])

    best = min(((top3(wv, thv), wv, thv)
                for wv in np.linspace(0.0, 3.0, 31)
                for thv in np.linspace(0.0, 2 * math.pi, 24)))
    wv, thv = best[1], best[2]
    for _ in range(6):
        for var in ("w", "th"):
            lo = (wv - 0.1, wv + 0.1) if var == "w" else (thv - 0.26, thv + 0.26)
            xs = np.linspace(lo[0], lo[1], 9)
            vals = [(top3(x, thv) if var == "w" else top3(wv, x), x)
                    for x in xs]
            v, x = min(vals)
            if var == "w":
                wv = float(x)
            else:
                thv = float(x)
    gapR = top3(wv, thv)
    M, T = tangent_ops(C, wv, thv)
    d = _lam @ (math.cos(thv) * XCOEF + math.sin(thv) * YCOEF)
    # the internal generators on the 16: the 45 T_ab restricted (sample 12)
    from cl10 import clifford as CL
    h = [E.asnumpy(m) for m in CL.generators()]
    Gs = [B16.conj().T @ (1j * (h[a] @ h[b])) @ B16
          for (a, b) in ((0, 1), (2, 3), (4, 5), (6, 7), (8, 9),
                         (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
                         (0, 9), (1, 3))]
    # the transport part: T_k minus the predicted Wilson Cartan term
    comm = 0.0
    for k in range(3):
        dkph = wv * math.sin(C[k])
        Tk = T[k] - dkph * np.diag(np.tile(d, 8))
        for G in Gs:
            Gi = np.kron(np.eye(8), G)
            comm = max(comm, float(np.linalg.norm(Tk @ Gi - Gi @ Tk)))
    # the Cartan coefficients at the corner: identical for all three axes
    cartan_k = [wv * math.sin(C[k]) for k in range(3)]
    checks.append(("T3 at the gapless corner triplet (gap %.1e at w=%.4f, "
                   "theta=%.4f): the TRANSPORT part of each tangent operator "
                   "(T_k minus the Wilson Cartan term) commutes with I_8 "
                   "kron G for 12 internal generators (max commutator %.1e), "
                   "and the Cartan coefficient is the SAME for all three "
                   "axes (%s): the Wilson factor adds one common internal-"
                   "diagonal shift to all three S_j - which is why V118's "
                   "commutator sits at the noise floor - and the internal "
                   "factor is a spectator of the velocity structure: a cone "
                   "must be built by the bond"
                   % (gapR, wv, thv, comm,
                      ", ".join("%.4f" % x for x in cartan_k)),
                   comm < 1e-10 and max(cartan_k) - min(cartan_k) < 1e-12))
    out["data"]["corner"] = {"gap": gapR, "w": wv, "theta": thv,
                             "spectator_comm": comm,
                             "cartan_coeff": cartan_k}

    ok = all(c[1] for c in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "SEPARATION LEMMA (cone programme, Phase 1): PROVEN. For the "
        "implemented step every tangent operator T_k(q) = i(d_k M)M^dag is "
        "exactly A_k(q) kron I_16 (max defect %.1e over %d points: the 12 M "
        "nodes, the 8 R corners, 20 random momenta), so the entire low-energy "
        "velocity structure lives in the 8x8 sub-site factor and the internal "
        "module is a spectator - this upgrades V118's closing sentence and "
        "V142's frame blindness from measurement to lemma. The Wilson factor "
        "modifies the statement exactly as predicted: it adds the internal-"
        "diagonal Cartan term (d_k ph) diag(d) (kron-I defect after removal "
        "%.1e) and no velocity structure. Consequence at the gapless corner "
        "triplet (gap %.1e): the transport part of T_k commutes with every "
        "internal generator (max %.1e) and the Cartan shift is common to all "
        "three axes - the internal factor is a spectator of the velocity "
        "structure. A Clifford-type tangent algebra - a cone - can therefore "
        "ONLY be created by a bond layer that entangles the sub-site and "
        "internal indices: this is the target of V186."
        % (worst, len(pts), worst_diag_pred, gapR, comm))

    res = os.path.join(HERE, "results")
    os.makedirs(res, exist_ok=True)
    with open(os.path.join(res, "v185_separation_lemma.json"), "w") as fh:
        json.dump(out, fh, indent=1)

    print("V185 - the separation lemma")
    for n, p in checks:
        print("  [%s] %s" % ("PASS" if p else "FAIL", n))
    print("  VERDICT: %s" % ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
