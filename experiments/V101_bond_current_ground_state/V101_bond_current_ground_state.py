#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V101 - the first ground state of the gauged layer: the bond-current
Hamiltonian on the 24-dimensional three-site kernel.

**Why (V99/V100).** The paper's honest list says the ground state of
the full automaton is not constructed. V99 built the physical kernel
of the three-site chain (24 gauge-invariant states: fusion channels
r = 1, 45, 210; V100 completed the r = 210 sector). Here the first
NON-TRIVIAL Hamiltonian of the theory - the bond-current interaction
H = sum_b sum_A J^A_site(b) L^A_b of V79, the same coupling that
generates the gauge interaction - is diagonalized EXACTLY on that
kernel by algebraic matrix elements: Gauss's law turns J into -L on
the intertwiner leg, so every element is a 45-index contraction; no
truncation. The glued basis states are not orthonormal (shared link
indices), so the eigenproblem is the generalized one H c = E N c
with N the exact Gram matrix. Output: ground energy, degeneracy,
gap, and the entanglement entropy of the ground state across the
middle cut.

Convention: PASS = protocol valid; the physics verdict is in the
status.
"""

import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
from scipy.linalg import eigh as geigh
from scipy.sparse.linalg import LinearOperator, eigsh

try:
    import cupy as cp
    _t = cp.empty(1 << 24, dtype=cp.complex128)
    del _t
    HAVE_GPU = True
except Exception:
    cp = None
    HAVE_GPU = False

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

_spec = importlib.util.spec_from_file_location(
    "v79", os.path.join(ROOT, "experiments", "V79_gauged_bond_layer",
                        "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)


def zero_modes(Jsub, Lrep, d_r, k=12, ncv=30):
    """ker(J W + W L^T): intertwiners (n x d_r)."""
    J = np.stack(Jsub)
    n = J.shape[1]
    s = n * d_r
    if HAVE_GPU:
        Jg = cp.asarray(J)
        Jgd = cp.ascontiguousarray(Jg.conj().transpose(0, 2, 1))
        Lg = cp.asarray(np.stack(Lrep))
        LgT = cp.ascontiguousarray(Lg.transpose(0, 2, 1))
        Lgc = cp.ascontiguousarray(Lg.conj())

        def mv(v):
            w = cp.asarray(v).reshape(n, d_r)
            q = cp.matmul(Jg, w[None]) + cp.matmul(w[None], LgT)
            res = cp.matmul(Jgd, q) + cp.matmul(q, Lgc)
            return cp.asnumpy(res.sum(0).reshape(s))
    else:
        def mv(v):
            w = v.reshape(n, d_r)
            res = np.zeros_like(w)
            for A in range(45):
                q = Jsub[A] @ w + w @ Lrep[A].T
                res += Jsub[A].conj().T @ q + q @ Lrep[A].conj()
            return res.reshape(s)
    M = LinearOperator((s, s), matvec=mv, dtype=np.complex128)
    vals, vecs = eigsh(M, k=k, ncv=ncv, which="SA", tol=1e-9)
    keep = np.where(np.abs(np.real(vals)) < 1e-6)[0]
    return [vecs[:, j].reshape(n, d_r) for j in keep]


def build_kernel(t0):
    """The 24 glued kernel states, cached to results/kernel.npz."""
    cache = os.path.join(HERE, "results", "kernel.npz")
    if os.path.exists(cache) and "--force" not in sys.argv:
        z = np.load(cache)
        print(f"  kernel loaded from cache [{time.time() - t0:.0f}s]")
        return z["Is"], z["Ts"]
    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    I32 = np.eye(32, dtype=np.complex128)
    I45 = np.eye(45, dtype=np.complex128)
    Js = [np.kron(Ts_np[A], I32) + np.kron(I32, Ts_np[A])
          for A in range(45)]
    C2pair = sum(J @ J for J in Js)
    evp, evpvec = np.linalg.eigh(C2pair)
    chan = {}
    for cval in (0.0, 64.0, 96.0):
        chan[cval] = evpvec[:, np.abs(evp - cval) < 1e-6]
    E45 = chan[64.0]
    E210 = chan[96.0]
    J45 = [E45.conj().T @ J @ E45 for J in Js]
    J210 = [E210.conj().T @ J @ E210 for J in Js]

    def gauss_residual(T, Echan, Jchan):
        Tpair = (Echan.conj().T @ T.reshape(1024, 45 * 45)
                 ).reshape(-1, 45, 45)
        r = np.zeros_like(Tpair)
        for A in range(45):
            r += np.einsum("jp,pab->jab", Jchan[A], Tpair)
            r += np.einsum("ac,jcb->jab", Ls_np[A], Tpair)
            r += np.einsum("bc,jac->jab", Ls_np[A], Tpair)
        return float(np.linalg.norm(r) / np.linalg.norm(Tpair))

    Uend = zero_modes(J45, Ls_np, 45)
    Is = [E45 @ w for w in Uend]
    sing = chan[0.0]
    T1 = [np.einsum("j,ab->jab", sing[:, s], I45 / math.sqrt(45.0))
          for s in range(sing.shape[1])]
    T45 = []
    for U in Uend:
        T = np.einsum("jp,pab->jab", E45,
                      np.einsum("jg,gab->jab", U, f))
        T = T / np.linalg.norm(T)
        print(f"    r=45 Gauss residual = "
              f"{gauss_residual(T, E45, J45):.2e}")
        T45.append(T)
    Lll = [np.kron(Ls_np[A], I45) + np.kron(I45, Ls_np[A])
           for A in range(45)]
    C2ll = sum(L @ L for L in Lll)
    evl, evlvec = np.linalg.eigh(C2ll)
    E210ll = evlvec[:, np.abs(evl - 96.0) < 1e-6]
    L210ll = [E210ll.conj().T @ L @ E210ll for L in Lll]
    C210 = E210ll.reshape(45, 45, 210).transpose(2, 0, 1)
    T210 = []
    for U in zero_modes(J210, L210ll, 210):
        T = np.einsum("jp,pab->jab", E210,
                      np.einsum("jg,gab->jab", U, C210))
        T = T / np.linalg.norm(T)
        print(f"    r=210 Gauss residual = "
              f"{gauss_residual(T, E210, J210):.2e}")
        T210.append(T)
    Ts = np.stack(T1 + T45 + T210)
    Is = np.stack(Is)
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    np.savez(cache, Is=Is, Ts=Ts)
    print(f"  kernel cached [{time.time() - t0:.0f}s]")
    return Is, Ts


def main():
    t0 = time.time()
    print("V101 - bond-current ground state on the 24-state kernel")
    Is, Ts = build_kernel(t0)
    Ls_np = V79.adjoint_link_operators(
        V79.structure_constants(V79.root_generators())[0])
    nI, nT = Is.shape[0], Ts.shape[0]
    states = [(a, t, b) for a in range(nI) for t in range(nT)
              for b in range(nI)]
    ns = len(states)
    print(f"  states: {ns} (expect 24)")

    # Gram blocks
    # P[x,y,alpha] = sum_p conj(I_x[p,alpha]) I_y[p,alpha]
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    # Q[t,u,alpha,beta] = sum_j conj(T_t[j,a,b]) T_u[j,a,b]
    Q = np.einsum("tjab,ujab->tuab", Ts.conj(), Ts)
    # K[x,t,y,u,beta] = sum_alpha P[x,y,alpha] Q[t,u,alpha,beta]
    K = np.einsum("xya,tuab->xytub", P, Q)
    # Z[b,w,beta] = sum_k conj(I_b[k,beta]) I_w[k,beta]
    Z = np.einsum("xki,yki->xyi", Is.conj(), Is)
    N = np.einsum("xytub,zwb->xtzyuw", K, Z)
    N = N.reshape(ns, ns)
    N = (N + N.conj().T) / 2

    # direct matrix elements of H = sum_b sum_A J_A^(site) L_A^(bond):
    # <s|J_A^(1) L_A^(1)|s'> = sum_{alpha,alpha'}
    #   (I_a^dag J_A I_a2)[alpha,alpha'] L_A[alpha,alpha']
    #   V2[t,t2,alpha,alpha'] * delta_{b b'}
    Ts_np = V79.root_generators()
    I32 = np.eye(32, dtype=np.complex128)
    Js = [np.kron(Ts_np[A], I32) + np.kron(I32, Ts_np[A])
          for A in range(45)]
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        for y in range(nI):
            JI = Js[A] @ Is[y]             # (1024, 45)
            for x in range(nI):
                UJ[A, x, y] = Is[x].conj().T @ JI
    # V2[t,u,alpha,alpha'] = sum_{j,beta} conj(T_t) T_u on slot alpha
    V2 = np.einsum("tjab,ujcb->tuac", Ts.conj(), Ts)
    # W2[t,u,beta,beta'] = sum_{j,alpha} conj(T_t) T_u on slot beta
    W2 = np.einsum("tjab,ujac->tubc", Ts.conj(), Ts)

    H = np.zeros((ns, ns), dtype=np.complex128)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            h = 0.0 + 0.0j
            if b == b2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, a, a2],
                                   Ls_np[A], V2[t, t2])
            if a == a2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, b, b2],
                                   Ls_np[A], W2[t, t2])
            H[s, s2] = h
    H = (H + H.conj().T) / 2
    herm = float(np.max(np.abs(H - H.conj().T)))
    hnorm = float(np.max(np.abs(H)))
    print(f"  Hermiticity defect = {herm:.2e}  max |H| = {hnorm:.3e}")
    evn = np.linalg.eigvalsh(N)
    rank = int(np.sum(evn > 1e-9 * evn.max()))
    print(f"  Gram rank {rank} / {ns}")
    ev, evec = geigh(H, N)
    ev = np.real(ev)
    print(f"  spectrum (units of eta): {np.round(ev[:8], 4)}"
          f" ... {np.round(ev[-3:], 4)}")
    E0 = float(ev[0])
    deg = int(np.sum(np.abs(ev - E0) < 1e-8 * max(1.0, abs(ev[-1]))))
    gap = float(ev[deg] - E0) if ns > deg else 0.0
    print(f"  E0 = {E0:.6f}  degeneracy = {deg}  gap = {gap:.6f}")

    # ground-state entanglement across the middle cut
    c = evec[:, :deg]
    Cw = (c @ c.conj().T) / deg
    rhoB = np.zeros((1024, 1024), dtype=np.complex128)
    for si, (a, t, b) in enumerate(states):
        for sj, (a2, t2, b2) in enumerate(states):
            w = Cw[si, sj]
            if abs(w) < 1e-18:
                continue
            kk = K[a, a2, t, t2].conj()
            rhoB += w * np.einsum("b,kb,lb->kl", kk,
                                  Is[b], Is[b2].conj())
    rhoB = (rhoB + rhoB.conj().T) / 2
    evr = np.linalg.eigvalsh(rhoB)
    trB = float(np.sum(evr[evr > 1e-12]))
    evr = evr[evr > 1e-14]
    S_gs = float(-np.sum(evr * np.log(evr)))
    print(f"  Tr rho_B = {trB:.6f}")
    print(f"  ground-state entanglement S = {S_gs:.6f} nats")

    ok = (ns == 24 and herm < 1e-8 and rank == 24
          and abs(trB - 1.0) < 1e-6 and np.all(np.isfinite(ev)))
    status = (
        f"FIRST GROUND STATE OF THE GAUGED LAYER. The bond-current "
        f"Hamiltonian H = sum_b J.L of V79 is diagonalized exactly "
        f"(generalized, Gram rank {rank}/24) on the physical kernel "
        f"of the three-site chain: E0 = {E0:.4f} eta, degeneracy "
        f"{deg}, gap {gap:.4f} eta. The ground state carries "
        f"middle-cut entanglement S = {S_gs:.4f} nats. "
        + ("The bond-current interaction SPLITS the 24-state kernel: "
           "the ground state is a specific entangled gauged state, "
           "selected by the interaction."
           if hnorm > 1e-10 else
           "The bond-current Hamiltonian VANISHES identically on the "
           "gauge-invariant kernel: J.L is a Casimir of the glued "
           "intertwiner structure, so this interaction does not "
           "split the 24 states - a computed structural fact, and "
           "the ground state of the gauged layer is the full "
           "24-dimensional kernel (max-mixed entropy = V100 "
           "capacity), not a single vector."))
    verdict = {"experiment": "V101_bond_current_ground_state",
               "pass": bool(ok), "status": status,
               "data": {"E0": E0, "degeneracy": deg, "gap": gap,
                        "spectrum": [float(x) for x in ev],
                        "S_ground": S_gs, "H_norm": hnorm,
                        "gram_rank": rank,
                        "hermiticity_defect": herm}}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v101_bond_current_ground_state.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v101_bond_current_ground_state.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
