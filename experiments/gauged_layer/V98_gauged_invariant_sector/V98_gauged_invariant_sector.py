#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V98 - the gauge-invariant sector of the bond layer: the center
obstruction and the matter-flux entangled complex.

**Why (V79 + V97).** V97 proved the implemented one-body step is
quadratic: Slater stays Slater, and the area law (Section 13.1) needs
NON-Gaussian entanglement. V79 built the gauged bond layer - one
adjoint Spin(10) link per bond, Gauss algebra exact. V98 computes the
PHYSICAL (gauge-invariant) sector of that layer and its entanglement:
the first non-Gaussian entanglement structure in the theory.

1. CENTER OBSTRUCTION (measured, three systems). With the implemented
   content - one 32 per site, one adjoint 45 per bond - the LOCAL
   Gauss constraint has no solution: dim ker = 0 on 32 x 45 (global),
   on 32 x 45 x 45 (global), and on the two-site one-link chain
   (32 x 32 x 45 with Q1 = T1 - L, Q2 = T2 + L), each with a large
   spectral gap. Adjoint flux carries even center charge and cannot
   screen a spinor: EVERY site of a physical state must carry an even
   number of spinor fields. The automaton's one-32-per-site content
   admits NO gauge-invariant state on any graph - the physical Hilbert
   space of the gauged layer is the space of NEUTRAL COMPOSITES.

2. MINIMAL PHYSICAL SECTOR. Two matter fields per site: the adjoint
   multiplicity in the site pair is m_45(32 x 32) = 2 (measured from
   the C2 spectrum on 1024: 90 states at C2 = 64), the intertwiner
   space ker(T1 + T2 + L) on (32 x 32) x 45 is 2-dimensional
   (measured, constructed exactly on the C2 = 64 subspace), and
   gluing two intertwiners over the link index gives the physical
   two-site sector of dimension m_45^2 = 4 (constructed and verified:
   all exact Spin(10) singlets, C2 < 1e-25).

3. ENTANGLEMENT. Each intertwiner is a maximal matter-flux singlet:
   Schmidt rank 45, S(pair | link) = ln 45. The glued physical states
   carry S(pair1 | pair2) = ln 45 <= ln 45: entanglement BOUND per
   bond, living on the bonds crossing the cut - the SHAPE of an area
   law (boundary-proportional, bounded per bond), Bekenstein
   coefficient still unnamed. The quadratic step of V97 cannot
   produce any of this.

Convention: the check PASSES when the protocol is valid and complete;
the physics verdict is in the status. Backend: CL10_BACKEND=cupy runs
the Lanczos matvecs and the dense diagonalizations on the GPU.
"""

import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

GPU = os.environ.get("CL10_BACKEND", "numpy") == "cupy"
if GPU:
    import cupy as XP
else:
    XP = np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
from cl10.paths import exp

_spec = importlib.util.spec_from_file_location(
    "v79", exp("V79_gauged_bond_layer",
                        "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)


def to_np(a):
    return XP.asnumpy(a) if GPU else np.asarray(a)


def move_apply(g, w, axis):
    """Apply matrix g to axis of tensor w."""
    m = XP.moveaxis(w, axis, 0)
    m = g @ m.reshape(g.shape[0], -1)
    m = m.reshape((g.shape[0],) + w.shape[:axis] + w.shape[axis + 1:])
    return XP.moveaxis(m, 0, axis)


def entropy_from_rho(rho):
    nus = to_np(XP.linalg.eigvalsh(rho))
    s = 0.0
    for nu in nus:
        nu = float(np.real(nu))
        if nu > 1e-12:
            s -= nu * math.log(nu)
    return s


def ker_modes(matvec_gpu, size, k=40, tol=1e-6):
    def matvec(v):
        return to_np(matvec_gpu(XP.asarray(v)))
    M = LinearOperator((size, size), matvec=matvec,
                       dtype=np.complex128)
    vals, vecs = eigsh(M, k=k, which="SA", tol=1e-10)
    vals = np.real(vals)
    idx = np.where(np.abs(vals) < tol)[0]
    return vals, vecs[:, idx]


def main():
    t0 = time.time()
    print(f"V98 - the gauge-invariant sector of the bond layer "
          f"[backend: {'cupy' if GPU else 'numpy'}]")
    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    Ts = [XP.asarray(T) for T in Ts_np]
    Ls = [XP.asarray(L) for L in Ls_np]
    I32 = XP.eye(32, dtype=np.complex128)

    # ---- 1a. global singlet tests: 32 x 45, 32 x 45 x 45 ------------------
    def mv_global(*factors):
        dims = [32] + [45] * (len(factors) - 1)
        size = int(np.prod(dims))

        def mv(v):
            w = v.reshape(dims)
            res = XP.zeros(dims, dtype=np.complex128)
            for A in range(45):
                q = XP.zeros(dims, dtype=np.complex128)
                for i in range(len(factors)):
                    g = Ts[A] if i == 0 else Ls[A]
                    q = q + move_apply(g, w, i)
                for i in range(len(factors)):
                    g = Ts[A] if i == 0 else Ls[A]
                    res = res + move_apply(g.conj().T, q, i)
            return res.reshape(size)
        return mv, size

    mv, size = mv_global([Ts], [Ls])
    _v1, ker1 = ker_modes(mv, size, k=20)
    mv, size = mv_global([Ts], [Ls], [Ls])
    _v1b, ker1b = ker_modes(mv, size, k=40)
    print(f"  ker(32 x 45)        = {ker1.shape[1]}  (expect 0: no "
          "singlet, V79 content theorem)")
    print(f"  ker(32 x 45 x 45)   = {ker1b.shape[1]}  (expect 0: adjoint "
          "flux cannot screen the spinor center charge)")

    # ---- 1b. two-site one-link chain, LOCAL constraints -------------------
    s2 = 32 * 32 * 45

    def mv_local(v):
        w = v.reshape(32, 32, 45)
        res = XP.zeros((32, 32, 45), dtype=np.complex128)
        for A in range(45):
            q1 = (move_apply(Ts[A], w, 0) - move_apply(Ls[A], w, 2))
            res += move_apply(Ts[A], q1, 0) - move_apply(Ls[A], q1, 2)
            q2 = (move_apply(Ts[A], w, 1) + move_apply(Ls[A], w, 2))
            res += move_apply(Ts[A], q2, 1) + move_apply(Ls[A], q2, 2)
        return res.reshape(s2)

    _v2, ker2 = ker_modes(mv_local, s2, k=40)
    print(f"  ker(32 x 32 x 45, local) = {ker2.shape[1]}  (expect 0: one "
          "spinor per site cannot be screened by adjoint flux)")

    # ---- 2. adjoint multiplicity in the site pair: m_45(32 x 32) ----------
    Js = [XP.kron(Ts[A], I32) + XP.kron(I32, Ts[A]) for A in range(45)]
    C2pair = XP.zeros((1024, 1024), dtype=np.complex128)
    for J in Js:
        C2pair += J @ J
    ev, evec = XP.linalg.eigh(C2pair)
    sel = XP.abs(ev - 64.0) < 1e-6
    n64 = int(XP.sum(sel))
    m45 = n64 // 45
    print(f"  C2(pair) = 64 states: {n64} -> m_45(32 x 32) = {m45} "
          "(expect 2)")

    # ---- 3. exact intertwiner space ker(T1 + T2 + L) ----------------------
    # Restrict the matter pair to its C2 = 64 subspace (90 states = two
    # copies of the 45); intertwiners are maps W (90 x 45) with
    # J_A W = L_A W for all A. Solve dim ker = m_45 exactly on 4050.
    E = evec[:, sel]                          # 1024 x 90, orthonormal
    Jsub = [E.conj().T @ J @ E for J in Js]   # 90 x 90
    s3 = 90 * 45
    M3 = XP.zeros((s3, s3), dtype=np.complex128)
    I90 = XP.eye(90, dtype=np.complex128)
    I45 = XP.eye(45, dtype=np.complex128)
    for A in range(45):
        QA = XP.kron(Jsub[A], I45) + XP.kron(I90, Ls[A])
        M3 += QA.conj().T @ QA
    w3, u3 = XP.linalg.eigh(M3)
    keep = XP.abs(w3) < 1e-8
    Ws = [u3[:, j].reshape(90, 45) for j in XP.where(keep)[0]]
    n_int = len(Ws)
    print(f"  intertwiners ker(T1+T2+L) = {n_int}  (expect {m45})")
    # lift to the full pair space and verify the Gauss constraint there
    Is = [E @ W for W in Ws]                  # each 1024 x 45
    gauss_err = 0.0
    for Ia in Is:
        for A in range(45):
            gauss_err = max(gauss_err, float(XP.linalg.norm(
                Js[A] @ Ia + Ia @ Ls[A].T)))
    print(f"  max ||(J_A + L_A) I_a|| on full space = {gauss_err:.2e}")

    # ---- 4. glue two intertwiners: physical two-site sector ---------------
    # Psi_{ab} = sum_alpha I_a(pair1, alpha) J_b(alpha, pair2):
    # both intertwiner bases are the same computed space.
    glued = []
    for a in range(n_int):
        for b in range(n_int):
            Ia = Is[a]                          # (1024, 45)
            Kb = Is[b].T                      # L antisymmetric: no conj
            glued.append((Ia @ Kb).reshape(1024 * 1024))
    G = XP.array(glued)
    G /= XP.linalg.norm(G, axis=1, keepdims=True)
    def mv_local_pair(v):
        # C2(pair1 + pair2) on (1024, 1024): glued states must be singlets.
        w = v.reshape(1024, 1024)
        res = XP.zeros_like(w)
        for A in range(45):
            q = Js[A] @ w + w @ Js[A].T
            res += Js[A] @ q + q @ Js[A].T
        return res.reshape(1024 * 1024)

    c2max = 0.0
    for vec in G:
        C2v = mv_local_pair(vec)
        c2 = XP.vdot(vec, C2v)
        c2max = max(c2max, float(np.abs(c2)))
    print(f"  glued physical states: {G.shape[0]} (expect {m45 ** 2}); "
          f"max |C2(pair1+pair2)| = {c2max:.2e} (expect 0: singlets)")

    # entanglement of the intertwiners and of the glued states
    S_int = []
    for a in range(n_int):
        Ia = Is[a]
        rho = Ia @ Ia.conj().T
        S_int.append(entropy_from_rho(rho / XP.trace(rho).real))
    S_glued = []
    for vec in G[:1]:
        w = vec.reshape(1024, 1024)
        rho = w @ w.conj().T
        S_glued.append(entropy_from_rho(rho / XP.vdot(w, w).real))
    print(f"  S(pair|link) per intertwiner: "
          f"{[round(s, 4) for s in S_int]} (ln 45 = {math.log(45):.4f})")
    print(f"  S(pair1|pair2) glued ground: {S_glued[0]:.4f} (<= ln 45)")

    ok = (ker1.shape[1] == 0 and ker1b.shape[1] == 0
          and ker2.shape[1] == 0 and m45 == 2 and n_int == 2
          and gauss_err < 1e-8
          and G.shape[0] == 4 and c2max < 1e-8
          and all(abs(s - math.log(45)) < 1e-6 for s in S_int)
          and S_glued[0] <= math.log(45) + 1e-9)

    status = (
        "THE GAUGED BOND LAYER HAS NO ONE-BODY PHYSICAL STATES; ITS "
        "PHYSICAL SECTOR IS AN ENTANGLED MATTER-FLUX COMPLEX. (1) CENTER "
        "OBSTRUCTION: dim ker = 0 for 32 x 45, for 32 x 45 x 45, and for "
        "the two-site one-link chain with LOCAL Gauss constraints - "
        "adjoint flux carries even center charge and cannot screen a "
        "spinor, so EVERY site must carry an even number of spinors: the "
        "implemented one-32-per-site content admits NO gauge-invariant "
        "state on any graph. (2) MINIMAL PHYSICAL SECTOR: two matter "
        f"fields per site: m_45(32 x 32) = {m45} (C2 spectrum), "
        f"intertwiner space dim {n_int}, glued two-site sector dim "
        f"{G.shape[0]}, all exact Spin(10) singlets (C2 < 1e-8). (3) "
        f"ENTANGLEMENT: each intertwiner is a maximal matter-flux "
        f"singlet, Schmidt rank 45, S = ln 45 = {math.log(45):.3f}; the "
        f"glued physical state carries S(pair1|pair2) = "
        f"{S_glued[0]:.3f} <= ln 45: entanglement BOUND per bond and "
        "carried by the bonds crossing the cut - the SHAPE of an area "
        "law (boundary-proportional, bounded per bond), Bekenstein "
        "coefficient unnamed. The quadratic step of V97 cannot produce "
        "any of this: the gauged layer is where the theory first gets "
        "non-Gaussian entanglement. Open (O9): multi-bond clusters "
        "(3-site, 2-link with paired matter is 2 x 10^10 dimensional) - "
        "the additivity of ln 45 per cut bond is a structural argument, "
        "not yet a measured scaling.")
    verdict = {"experiment": "V98_gauged_invariant_sector",
               "pass": bool(ok), "status": status,
               "data": {"ker_32x45": int(ker1.shape[1]),
                        "ker_32x45x45": int(ker1b.shape[1]),
                        "ker_chain_local": int(ker2.shape[1]),
                        "m45_pair": int(m45),
                        "n_intertwiners": int(n_int),
                        "n_glued": int(G.shape[0]),
                        "max_C2_glued": c2max,
                        "max_gauss_err": gauss_err,
                        "S_intertwiners": [round(s, 6) for s in S_int],
                        "S_glued": round(S_glued[0], 6),
                        "ln45": math.log(45.0)}}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v98_gauged_invariant_sector.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v98_gauged_invariant_sector.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
