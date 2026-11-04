#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V372 - torus singlet tower census (booked BEFORE the run in
review/pre_registrations_root.md).

Question: does the two-quanta + flux composite class of the N=2^3
torus carry GAUSS-SINGLET channels beyond the single flux line?
The V370 pure-phase consequence (booked, no run): with ONE flux line
(ell=1) the deformation Sum_b J_b.L_b = -64 on every torus singlet
(V370 verbatim on the one touching bond; every other bond carries
neither matter nor flux), so e^{i eta J.L} is a global phase there:
no threshold move at ell=1. Binding needs the ell=2 sector, where
interior Gauss (L_f1 + L_f2)|psi> = 0 makes the touching contraction
(J_A - J_B).L_f1 - not a phase.

Machine (booked):
 A1 ker(32 (x) 32) = 2 (dense eigvalsh, 1024).
 A2 ker(45 (x) 45) = 1 (dense eigvalsh, 2025).
 A3 ker(E4) = ker(32x32x45) = 2 (LinearOperator eigsh SM, V366).
 W  exact Weyl-character census: mult(1 in V) =
    (1/|W|) sum_w det(w) m_V(w rho - rho) over the weight sumset;
    validated by A1-A3 (all three must match) and by the Casimir
    trace gates tr(sum T_a^2) = 40 on both 32 and 45.
 M  ell=2: mult(1 in 32x32x45x45) - the open number.
 T  torus counts N=2^3: 8 sites, 12 bonds, 24 adjacent bond pairs,
    3 bond orbits under translations.
Verdicts: TORUS-SINGLET-TOWER-ALIVE / -LINE1-ONLY / EXT-GROUP-MISMATCH
/ MACHINERY-DRIFT.
"""
import itertools
import json
import os
import sys
from collections import Counter

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "gauged_layer",
                                "V79_gauged_bond_layer"))
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)
from cl10 import altcell as AC  # noqa: E402

KTOL = 1e-8


def main():
    checks, data = [], {}

    def log(m):
        print(m, flush=True)

    Ts = [np.asarray(t, dtype=np.complex128) for t in root_generators()]
    f, _k, _a = structure_constants(Ts)
    Ls = [np.asarray(l, dtype=np.complex128) for l in adjoint_link_operators(f)]

    # ---- Casimir trace gate (normalization-free): the D5 weight model
    # predicts tr(sum T_a^2) = 32*11.25 = 360u on the 32 and 45*16 =
    # 720u on the adjoint (roots of length^2 2), so the ratio
    # tr45/tr32 must be exactly 2 whatever the generators' unit.
    tr32 = sum(np.trace(t @ t) for t in Ts)
    tr45 = float(np.real(sum(np.trace(l @ l) for l in Ls)))
    g1 = (abs(np.imag(tr32)) < 1e-9 and abs(2.0 * np.real(tr32) - tr45)
          < 1e-8 * max(1.0, abs(tr45)))
    checks.append(("Casimir trace ratio gate tr(C2,45)/tr(C2,32) = "
                   "%.6f (D5 weight model: 2.000000): %s"
                   % (tr45 / np.real(tr32), g1), g1))
    log("trace gate: tr C2(32) = %s, tr C2(45) = %.6f"
        % (np.round(tr32, 9), tr45))

    # ---- A1: dense ker(32 (x) 32) ----
    I32 = np.eye(32)
    C = np.zeros((1024, 1024), dtype=np.complex128)
    for t in Ts:
        q = np.kron(t, I32) + np.kron(I32, t)
        C += q @ q
    evA1 = np.linalg.eigvalsh((C + C.conj().T) / 2.0)
    kerA1 = int(np.sum(np.abs(evA1) < KTOL))
    a1 = kerA1 == 2
    checks.append(("A1 dense ker(32x32) = %d (booked 2), min eigs %s: %s"
                   % (kerA1, ["%.2e" % z for z in evA1[:4]], a1), a1))

    # ---- A2: dense ker(45 (x) 45) ----
    I45 = np.eye(45)
    Cb = np.zeros((2025, 2025), dtype=np.complex128)
    for l in Ls:
        q = np.kron(l, I45) + np.kron(I45, l)
        Cb += q @ q
    evA2 = np.linalg.eigvalsh((Cb + Cb.conj().T) / 2.0)
    kerA2 = int(np.sum(np.abs(evA2) < KTOL))
    a2 = kerA2 == 1
    checks.append(("A2 dense ker(45x45) = %d (booked 1), min eigs %s: %s"
                   % (kerA2, ["%.2e" % z for z in evA2[:4]], a2), a2))

    # ---- A3: eigsh SM on E4 = 32x32x45 (46080), V366 pattern ----
    def e4_matvec(x):
        h = np.asarray(x).reshape(32, 32, 45)
        y = np.zeros_like(h)
        for a in range(45):
            q = (np.einsum("Aa,abf->Abf", Ts[a], h)
                 + np.einsum("Bb,abf->aBf", Ts[a], h)
                 + np.einsum("Ff,abf->abF", Ls[a], h))
            y += (np.einsum("Aa,abf->Abf", Ts[a], q)
                  + np.einsum("Bb,abf->aBf", Ts[a], q)
                  + np.einsum("Ff,abf->abF", Ls[a], q))
        return y.ravel()
    op3 = spla.LinearOperator((46080, 46080), matvec=e4_matvec,
                              dtype=np.complex128)
    w3 = spla.eigsh(op3, k=10, which="SM", return_eigenvectors=False)
    kerA3 = int(np.sum(np.abs(w3) < KTOL))
    a3 = kerA3 == 2
    checks.append(("A3 eigsh SM ker(32x32x45) = %d (booked 2), min |eigs|"
                   " %s: %s" % (kerA3, ["%.2e" % abs(z) for z in
                                        sorted(w3, key=abs)[:4]], a3), a3))
    log("A1 %d  A2 %d  A3 %d" % (kerA1, kerA2, kerA3))

    # ---- exact dense-Casimir cluster census (machine correction
    # booked before the run: cluster method, Weyl formula dropped).
    # Model cluster SHAPES (booked from the weight calculus):
    # 32x32: dims {2, 90, 420, 20, 240, 252} (= 1^2, 45^2, 210^2,
    # 120+120bar, 126+126bar); 45x45: dims {1, 45, 54, 210,
    # 770, 945} (sym 1+54+945 = 1001, antisym 45+210+770 = 1025).
    # The 16x16bar = 1+45+210 anchors the C2 values of
    # the trivial, 45 and 210 sectors; pairing at any other shared
    # C2 value needs a unique common irrep-dim factor to contribute,
    # else 0 (no such shared value occurs on the machine); ambiguity
    # -> EXT-GROUP-MISMATCH.
    L_DIMS = (10, 45, 54, 120, 126, 210, 770, 945)

    def clusters(C, tol=1e-6):
        ev = np.linalg.eigvalsh((C + C.conj().T) / 2.0)
        out = []
        k = 0
        scale = max(1.0, float(np.max(np.abs(ev))))
        while k < len(ev):
            m = k
            while m < len(ev) and abs(ev[m] - ev[k]) <= tol * scale:
                m += 1
            out.append((float(np.mean(ev[k:m])), m - k))
            k = m
        return out

    def pair_build(repA, repB, nA, nB):
        C = np.zeros((nA * nB, nA * nB), dtype=np.complex128)
        Ia, Ib = np.eye(nA), np.eye(nB)
        for ta, tb in zip(repA, repB):
            q = np.kron(ta, Ib) + np.kron(Ia, tb)
            C += q @ q
        return C

    cM = clusters(pair_build(Ts, Ts, 32, 32))          # 32x32
    cF = clusters(pair_build(Ls, Ls, 45, 45))          # 45x45
    T16 = [np.asarray(t, dtype=np.complex128)
           for t in AC.weyl_generators()[1]]
    T16c = [-t.T for t in T16]        # V366 verbatim conjugate
    cP = clusters(pair_build(T16, T16c, 16, 16))       # 16x16bar

    pP = sorted(d for _, d in cP)
    pM = sorted(d for _, d in cM)
    pF = sorted(d for _, d in cF)
    okshape = (pP == [1, 45, 210]
               and pM == sorted([2, 90, 420, 20, 240, 252])
               and pF == sorted([1, 45, 54, 210, 770, 945]))
    checks.append(("cluster-shape gates: 16x16bar %s (1/45/210), "
                   "32x32 %s (20/240/252/2/420/90), 45x45 %s "
                   "(1/45/54/210/770/945): %s" % (pP, pM, pF,
                                                  okshape), okshape))

    table = {d: v for v, d in cP if d in (1, 45, 210)}   # rdim -> C2

    def value_at(v, cs, tol=1e-6):
        for vv, d in cs:
            if abs(vv - v) <= tol * max(1.0, abs(vv)):
                return d
        return 0

    def shared_mult(clA, clB):
        tot, unknown = 0, []
        for v, dA in clA:
            dB = value_at(v, clB)
            if dB == 0:
                continue
            rA = next((r for r, rv in table.items()
                       if abs(rv - v) <= 1e-6), None)
            if rA is not None:
                if dA % rA or dB % rA:
                    unknown.append((round(v, 6), dA, dB, rA))
                else:
                    tot += (dA // rA) * (dB // rA)
                continue
            if abs(v) < 1e-9:
                continue                     # r=1 anchored at v=0 only
            cand = [r for r in L_DIMS if dA % r == 0 and dB % r == 0]
            if len(cand) == 1:
                r = cand[0]
                tot += (dA // r) * (dB // r)
            elif len(cand) == 0:
                continue                     # collision, no common r
            else:
                unknown.append((round(v, 6), dA, dB, cand))
        return tot, unknown

    kerW_M = value_at(0.0, cM)               # trivial copies in 32x32
    kerW_F = value_at(0.0, cF)               # trivial copies in 45x45
    mult_E4, unE4 = shared_mult(cM, [(table[45], 45)])
    mult_E2, unE2 = shared_mult(cP, [(table[45], 45)])
    kerW_E1 = value_at(0.0, cP)
    n_2, un2 = shared_mult(cM, cF)
    unknown_all = unE4 + unE2 + un2
    n_m2, n_ff, n_1 = (kerW_M, 0), (kerW_F, 0), (mult_E4, 0)
    okgm = (n_m2 == (kerA1, 0) and n_ff == (kerA2, 0)
            and n_1 == (kerA3, 0) and kerW_E1 == 1 and mult_E2 == 1
            and not unknown_all and okshape)
    checks.append(("cluster census matches A1-A3 + E1 + E2: "
                   "%d/%d/%d/%d/%d vs %d/%d/%d/1/1; unknown shared "
                   "clusters %s: %s"
                   % (kerW_M, kerW_F, mult_E4, kerW_E1, mult_E2,
                      kerA1, kerA2, kerA3, unknown_all, okgm), okgm))
    log("cluster census: ker(32x32) %d, ker(45x45) %d, E4 %d, E1 %d,"
        " E2 %d, ell=2 %d" % (kerW_M, kerW_F, mult_E4, kerW_E1,
                              mult_E2, n_2))

    # ---- torus counts, N=2^3 ----
    sites = list(itertools.product(*[[0, 1]] * 3))
    bonds = set()
    for x in sites:
        for d in range(3):
            y = list(x)
            y[d] = (y[d] + 1) % 2
            bonds.add(frozenset((x, tuple(y))))
    nbonds = len(bonds)
    adj = 0
    for b1, b2 in itertools.combinations(sorted(bonds, key=lambda b:
                                                tuple(sorted(b))), 2):
        if b1 & b2:
            adj += 1
    # bond orbits under translations (xor by t)
    bset = list(bonds)

    def find(a):
        par[a] = par[par[a]]
        return par[a]
    par = list(range(nbonds))

    def uni(a, b):
        ra, rb = find(a), find(b)
        par[ra] = rb
    for t in itertools.product(*[[0, 1]] * 3):
        for i, b in enumerate(bset):
            j = bset.index(frozenset(tuple((x[d] ^ t[d]) for d in range(3))
                                     for x in b))
            uni(i, j)
    norb = len({find(i) for i in range(nbonds)})
    tc = (nbonds == 12 and adj == 24 and norb == 3 and len(sites) == 8)
    checks.append(("torus N=2^3 counts: %d bonds, %d adjacent pairs, "
                   "%d bond orbits (booked 12/24/3): %s"
                   % (nbonds, adj, norb, tc), tc))

    ker_ell1, ker_ell2 = n_1[0], n_2
    data["anchors"] = {"ker_32x32": kerA1, "ker_45x45": kerA2,
                       "ker_E4": kerA3, "min_eigs_32x32":
                       [float(z) for z in evA1[:4]],
                       "min_eigs_45x45": [float(z) for z in evA2[:4]],
                       "min_abs_eigs_E4": sorted(float(abs(z))
                                                 for z in w3)[:4]}
    data["cluster_census"] = {"ker_ell1": ker_ell1, "ker_ell2": ker_ell2,
                                                      "dim_ell1_fiber": 46080,
                           "dim_ell2_fiber": 2073600}
    data["torus"] = {"sites": 8, "bonds": nbonds, "adjacent_pairs": adj,
                     "bond_orbits": norb}
    data["pure_phase_ell1"] = (
        "booked consequence of V370 (no run here): on the torus with a "
        "single flux line every other bond is empty of both matter and "
        "flux, so Sum_b J_b.L_b = -64 on the whole ell=1 singlet tower "
        "and the deformed step is e^{-64 i eta} x identity there - no "
        "threshold move at ell=1; the ell=2 sector (interior Gauss "
        "(L_f1+L_f2)|psi>=0, touching contraction (J_A-J_B).L_f1) is "
        "the first place the published interaction can act.")

    ok_machine = g1 and a1 and a2 and a3
    if not ok_machine:
        status, ok = "MACHINERY-DRIFT", False
    elif not okgm:
        status, ok = "EXT-GROUP-MISMATCH", False
    elif ker_ell2 >= 1:
        status = ("TORUS-SINGLET-TOWER-ALIVE (ell=1: 2 channels x 12 "
                  "bonds, ell=2: %d channels x %d adjacent pairs -> "
                  "the V371 deformed-channel cut is licensed on the "
                  "ell<=2 sector; ell=1 alone is pure phase and "
                  "cannot bind)" % (ker_ell2, adj))
        ok = True
    else:
        status = ("TORUS-SINGLET-TOWER-LINE1-ONLY (the published "
                  "interaction CANNOT open a threshold through any "
                  "ell<=2 singlet channel: V371 closes at ell<=2; the "
                  "next honest steps are ell>=3 or a content "
                  "extension)")
        ok = True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0][:400])
                    for c in checks), flush=True)
    print("[V372 done] %s" % status, flush=True)
    out = {"id": "V372_torus_singlet_tower_census", "tier": 2,
           "residual": True,
           "property": "gauss-singlet channels of the torus composite "
                       "tower beyond one flux line (V371 prerequisite)",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v372_torus_singlet_tower_census.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
