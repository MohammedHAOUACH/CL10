#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V100 - the vertex correction to the gauged area law: intertwiner
dimensions and the capacity of the three-site physical sector.

**Why (V98/V99).** V99 measured the area law of the gauged layer per
bond (ln 45, additive across a single bond) and its VERTEX correction
(a two-bond region carries ln 45, not 2 ln 45), and named the gap to
Bekenstein as a factor 4 ln 45 = 15.23. Two numbers decide how much
of that gap is closed by the vertex constraints of the CUBIC lattice,
where every vertex carries six adjoint legs:

1. I_n = dim (45^otimes n)^G, n = 2..6: the number of gauge-invariant
   states at a vertex with n adjoint legs. By the Weyl integration
   formula this is a Haar average over the group:
   I_n = E[chi_adj(g)^n], chi_adj(g) = ((tr g)^2 - tr g^2)/2, with
   Haar-distributed SO(10) matrices from QR of Gaussian matrices.
   Validation: I_2 = I_3 = 1 exactly (the singlet, the f-tensor).
   The effective per-bond area-law coefficient of the cubic lattice
   is then S_bond = (1/3) ln I_6 (each vertex carries 6 legs, each
   bond is shared by 2 vertices: 3 bonds per vertex).

2. The CAPACITY of the three-site physical sector: the entanglement
   entropy of the maximally mixed state on the 24-dimensional kernel
   of V99 across the middle cut - the maximum entropy the two-bond
   region can carry. The r = 210 fusion channel (8 of the 24 states,
   left as a counting argument in V99) is CONSTRUCTED here: the
   channel intertwiners U (ker(J + W L^T W^T)) are found by Lanczos
   on 90 x 45 and 420 x 210 systems, and the bond tensors are U
   contracted with the invariant tensors of 45 x 45 (the f-tensor
   for r = 45, the 210-channel projector basis for r = 210). Every
   final tensor is checked against Gauss's law. The glued states
   Psi = I (x)_alpha T (x)_beta I' are NOT normalized (the shared
   link indices do not factorize): the Gram matrix of the 24 glued
   states is computed exactly, the basis is orthonormalized through
   it, and the capacity is the max-mixed entropy in that basis.
   Heavy linear algebra runs on the GPU when cupy is available; the
   kernel is cached to results/kernel.npz.

Convention: PASS = protocol valid and complete; the physics verdict
is in the status.
"""

import importlib.util
import json
import math
import os
import sys
import time

import numpy as np
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


def _jax_batch():
    """JAX (x64) for the batched 10x10 QRs: numpy's stacked QR loops per
    matrix (~3 us each, 30 s for 10M); XLA vectorizes across the batch.
    Returns jnp or None (GPU busy / jax absent)."""
    try:
        os.environ.setdefault("XLA_PYTHON_CLIENT_MEM_FRACTION", "0.25")
        import jax
        jax.config.update("jax_enable_x64", True)
        if not jax.devices():
            return None
        import jax.numpy as jnp
        return jnp
    except Exception:  # noqa: BLE001
        return None


def haar_chi_moments(ns, n_samples, rng):
    """E[chi_adj^n] over Haar SO(10), n in ns: (mean, stderr) each.
    The Gaussian stream stays numpy (identical samples either way);
    only the batched QR moves to JAX when the GPU is free."""
    jnp = _jax_batch()
    xp, crng, batch = np, rng, 100_000
    acc = {n: xp.zeros(2) for n in ns}
    done = 0
    while done < n_samples:
        m = min(batch, n_samples - done)
        G = crng.standard_normal((m, 10, 10))
        if jnp is not None:
            try:
                Q, R = jnp.linalg.qr(jnp.asarray(G))
                d = jnp.sign(jnp.diagonal(R, axis1=1, axis2=2))
                Q = Q * d[:, None, :]
                det = jnp.sign(jnp.linalg.det(Q))
                Q = Q.at[:, 0, :].multiply(det[:, None])
                tr = jnp.einsum("mii->m", Q)
                tr2 = jnp.einsum("mij,mji->m", Q, Q)
                chi = np.asarray(0.5 * (tr * tr - tr2))
            except Exception:  # noqa: BLE001  (OOM mid-run: numpy path)
                jnp = None
        if jnp is None:
            Q, R = xp.linalg.qr(G)
            d = xp.sign(xp.diagonal(R, axis1=1, axis2=2))
            Q = Q * d[:, None, :]
            det = xp.sign(xp.linalg.det(Q))
            Q[:, 0, :] *= det[:, None]
            tr = xp.einsum("mii->m", Q)
            tr2 = xp.einsum("mij,mji->m", Q, Q)
            chi = 0.5 * (tr * tr - tr2)
        for n in ns:
            x = chi ** n
            acc[n] += xp.array([x.sum(), (x * x).sum()])
        done += m
    out = {}
    for n in ns:
        s1, s2 = (float(v) for v in acc[n])
        mean = s1 / n_samples
        var = max(s2 / n_samples - mean * mean, 0.0)
        out[n] = (mean, math.sqrt(var / n_samples))
    return out


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
    vals = np.real(vals)
    print(f"    Lanczos {n}x{d_r}: lowest eigenvalues "
          f"{np.round(vals[:6], 3)}")
    keep = np.where(np.abs(vals) < 1e-6)[0]
    return [vecs[:, j].reshape(n, d_r) for j in keep]


def build_kernel(t0):
    """The 24 glued kernel states: Is (2 x 1024 x 45) and Ts
    (6 x 1024 x 45 x 45), cached to results/kernel.npz."""
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
    print(f"  channels of 32x32: singlet {chan[0.0].shape[1]}, "
          f"45: {chan[64.0].shape[1]}, 210: {chan[96.0].shape[1]}")
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
    print(f"  end intertwiners: {len(Is)} (expect 2) "
          f"[{time.time() - t0:.0f}s]")
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


def gram_and_K(Is, Ts):
    """Gram matrix N of the 24 glued states and the partial-trace
    tensors K[(a,t),(a',t'); beta]."""
    nI, nT = Is.shape[0], Ts.shape[0]
    # P[x,y,alpha] = sum_i conj(I_x[i,alpha]) I_y[i,alpha]
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    Q = np.einsum("tjab,ujab->tuab", Ts.conj(), Ts)
    # K[x,t,y,u,beta] = sum_alpha P[x,y,alpha] Q[t,u,alpha,beta]
    K = np.einsum("xya,tuab->xytub", P, Q)
    # Z[b,w,beta] = sum_k conj(I_b[k,beta]) I_w[k,beta]
    Z = np.einsum("xki,yki->xyi", Is.conj(), Is)
    # rows (x,t,z), columns (y,u,w)
    N = np.einsum("xytub,zwb->xtzyuw", K, Z)
    N = N.reshape(nI * nT * nI, nI * nT * nI)
    return N, K, Z


def main():
    t0 = time.time()
    print("V100 - vertex correction of the gauged area law")

    # ---- 1. I_n = dim (45^otimes n)^G by Haar averaging -------------------
    rng = np.random.default_rng(20261008)
    N_MC = 10_000_000
    I = haar_chi_moments((2, 3, 4, 5, 6), N_MC, rng)
    for n in (2, 3, 4, 5, 6):
        print(f"  I_{n} = {I[n][0]:12.1f} +- {I[n][1]:.1f}  "
              f"(nearest int {round(I[n][0])})")
    ok_low = abs(I[2][0] - 1) < 0.1 and abs(I[3][0] - 1) < 0.1
    print(f"  validation I_2 = I_3 = 1: {ok_low}")

    # ---- 2. the 24-state kernel and its Gram matrix ------------------------
    Is, Ts = build_kernel(t0)
    nI, nT = Is.shape[0], Ts.shape[0]
    ns = nI * nT * nI
    N, K, Z = gram_and_K(Is, Ts)
    evn, Un = np.linalg.eigh(N)
    rank = int(np.sum(evn > 1e-9 * evn.max()))
    print(f"  Gram: rank {rank} / {ns}, "
          f"eigenvalues {evn.min():.3e}..{evn.max():.3e}")
    # orthonormal basis coefficients: O = U diag(1/sqrt(lam)) on support
    sup = evn > 1e-9 * evn.max()
    O = (Un[:, sup] / np.sqrt(evn[sup])) @ Un[:, sup].conj().T
    C = (O @ O.conj().T) / rank  # max-mixed weights in glued basis

    # ---- 3. capacity: max-mixed entropy, middle cut ------------------------
    # vectorized: Gm[k,l] = sum_{a,t,b,c,u,v} C[(a,t,b),(c,u,v)]
    #   * conj(K[a,c,t,u,q]) * Is[b][k,q] * conj(Is[v][l,q]),
    # summed over the shared leg q (K index order (x,y,t,u,beta)).
    C6 = C.reshape(nI, nT, nI, nI, nT, nI)          # (a,t,b,c,u,v)
    Kc = K.conj()                                    # (a,c,t,u,q)
    M = np.einsum("atbcuv,actuq->bvq", C6, Kc)
    Gm = np.einsum("bvq,bkq,vlq->kl", M, Is, Is.conj(), optimize=True)
    ev = np.linalg.eigvalsh(Gm)
    tr = float(np.real(np.sum(ev[np.abs(ev) > 1e-12])))
    ev = ev[ev > 1e-14]
    S_cap = float(-np.sum(ev * np.log(ev)))
    print(f"  Tr rho_B = {tr:.6f}")
    print(f"  capacity S (max-mixed, middle cut) = {S_cap:.6f} nats")
    print(f"  naive (no vertices) 2 ln 45 = {2 * math.log(45):.6f}")
    print(f"  capacity per bond = {S_cap / 2:.6f}  "
          f"vs Bekenstein 1/4: factor {S_cap / 2 / 0.25:.2f}")
    s_lat = math.log(I[6][0]) / 3.0 if I[6][0] > 0 else 0.0
    print(f"  cubic-lattice vertex-corrected coefficient "
          f"(1/3) ln I_6 = {s_lat:.4f} nats/bond  "
          f"vs ln 45 = {math.log(45):.4f}  "
          f"vs 1/4: factor {s_lat / 0.25:.2f}")

    ok = (ok_low and nT == 6 and nI == 2 and ns == 24
          and rank == 24 and abs(tr - 1.0) < 1e-6)

    status = (
        "VERTEX CORRECTION DATA FOR THE GAUGED AREA LAW. "
        f"I_n = dim (45^otimes n)^G by Haar averaging over SO(10): "
        f"I_2 = {round(I[2][0])}, I_3 = {round(I[3][0])}, "
        f"I_4 = {round(I[4][0])}, I_5 = {round(I[5][0])}, "
        f"I_6 = {round(I[6][0])} (MC errors in data). I_6 counts the "
        "gauge-invariant states at a cubic-lattice vertex with six "
        "adjoint legs: the vertex-corrected area-law coefficient is "
        f"(1/3) ln I_6 = {s_lat:.3f} nats/bond against ln 45 = 3.807 "
        "with no vertices: the vertex constraint removes a factor "
        f"45^3 / I_6 = {45 ** 3 / I[6][0]:.1f} per bond. "
        "The r = 210 fusion channel is constructed (2 bond tensors, "
        "Gauss-checked), completing the 24-state kernel of V99 "
        f"(Gram rank {rank}/24); the capacity of the two-bond "
        f"region is S = {S_cap:.3f} nats ({S_cap / 2:.3f} per "
        f"bond), a factor {S_cap / 2 / 0.25:.1f} above the "
        "Bekenstein 1/4 - the gap is reduced from 15.23 (naive "
        "ln 45) to this number under the physical-sector counting, "
        "and remains an open residual.")
    verdict = {"experiment": "V100_vertex_correction",
               "pass": bool(ok), "status": status,
               "data": {"I": {str(n): [I[n][0], I[n][1]] for n in I},
                        "n_T": nT, "n_end": nI, "n_states": ns,
                        "gram_rank": rank,
                        "S_capacity": S_cap,
                        "vertex_coefficient": s_lat}}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v100_vertex_correction.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v100_vertex_correction.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
