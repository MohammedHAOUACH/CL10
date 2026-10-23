#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""CL10 alt-interactions GPU route (numba CUDA, exploration-only).

The heavy loop of the interaction census (V172/V173) is the batched
eigendecomposition of the 128x128 Bloch step U(q) over the momentum
grid (13824 / 32768 points).  numba.cuda has no LAPACK, so this module
implements the batched solver as a CUDA kernel in the gtc2017-numba
style: one thread block per matrix, cyclic JACOBI sweeps for NORMAL
matrices (Fasi/HSVD variant): each 2x2 core is diagonalized by the
UNITARY polar factor G = X (X^dagger X)^{-1/2} of its eigenvector
matrix X (closed form for 2x2 Hermitian positive definite), applied as
the similarity U <- G^dagger U G; V accumulates the eigenvectors.
Converges quadratically for normal (here unitary) matrices: CPU
simulation on random 128x128 unitaries reaches off-diagonal 3e-10 in
59 sweeps, eigenphase error 8e-15 vs numpy eig.

The module is a VERIFICATION route, not a verdict: ``python3 -m
cl10.altgpu`` recomputes the V172 Fock self-energy on the GPU and
compares with a numpy pass on the SAME Bloch matrices.  The automaton's
anti-unitary symmetry pairs the spectrum in +-theta: at points where
the two smallest |theta| coincide the single-mode projector is an
arbitrary basis choice, so the comparison uses the convention-
invariant SUBSPACE projector there (single mode when isolated).
Agreement < 1e-6 confirms the GPU route.  It never touches the master
modules.
"""
import math
import sys
import time

import numpy as np

from . import altcell as A
from . import supercell as SC
from . import config as _config
_TOL = _config.CFG.tolerances   # the shared bands: cl10/config.yaml

N = 128   # the 16-sector Bloch dimension (8 sub-sites x 16 internal)


def build_kernels():
    from numba import cuda

    @cuda.jit(device=True)
    def _cabs2(z):
        return z.real * z.real + z.imag * z.imag

    @cuda.jit(device=True)
    def _cconj(z):
        return complex(z.real, -z.imag)

    @cuda.jit(device=True)
    def _csqrt(z):
        r = math.sqrt(z.real * z.real + z.imag * z.imag)
        if r == 0.0:
            return 0.0 + 0.0j
        t = math.sqrt(0.5 * (r + z.real))
        u = 0.5 * z.imag / t if t > 1e-300 else 0.0
        return t + u * 1j

    @cuda.jit
    def jacobi_normal(U, V, n, sweeps, tol):
        """Batched cyclic Jacobi for normal matrices: U (m, n, n) unitary
        -> diagonal by unitary polar-factor similarities; V accumulates
        the eigenvectors.  Launch one block of n threads PER MATRIX:
        m = blockIdx.x."""
        m = cuda.blockIdx.x
        Um = U[m]
        Vm = V[m]
        sh = cuda.shared.array(1, dtype=np.float64)
        for sweep in range(sweeps):
            for p in range(n - 1):
                for q in range(p + 1, n):
                    a = Um[p, p]
                    b = Um[p, q]
                    c0 = Um[q, p]
                    d = Um[q, q]
                    if _cabs2(b) + _cabs2(c0) < 1e-28:
                        continue
                    tr = a + d
                    dif = a - d
                    disc = dif * dif + 4.0 * b * c0
                    sr = _csqrt(disc)
                    l1 = 0.5 * (tr - sr)
                    l2 = 0.5 * (tr + sr)
                    # eigenvector matrix X of the core (columns)
                    if _cabs2(b) >= _cabs2(c0):
                        x00, x10 = b, l1 - a
                        x01, x11 = b, l2 - a
                    else:
                        x00, x10 = l1 - d, c0
                        x01, x11 = l2 - d, c0
                    # H = X^dagger X (2x2 Hermitian positive definite)
                    h00 = _cabs2(x00) + _cabs2(x10)
                    h11 = _cabs2(x01) + _cabs2(x11)
                    h01 = _cconj(x00) * x01 + _cconj(x10) * x11
                    detH = h00 * h11 - _cabs2(h01)
                    if detH < 1e-30:
                        continue
                    s = math.sqrt(detH)
                    t = math.sqrt(h00 + h11 + 2.0 * s)
                    st = s * t
                    # H^{-1/2} = (1/(s t)) [[h11+s, -h01], [-conj(h01), h00+s]]
                    i00 = (h11 + s) / st
                    i11 = (h00 + s) / st
                    i01 = -h01 / st
                    i10 = -_cconj(h01) / st
                    # G = X H^{-1/2} (unitary)
                    g00 = x00 * i00 + x01 * i10
                    g01 = x00 * i01 + x01 * i11
                    g10 = x10 * i00 + x11 * i10
                    g11 = x10 * i01 + x11 * i11
                    cuda.syncthreads()
                    # rows p,q <- G^dagger rows (thread t owns col t)
                    t0 = cuda.threadIdx.x
                    r0 = Um[p, t0]
                    r1 = Um[q, t0]
                    Um[p, t0] = _cconj(g00) * r0 + _cconj(g10) * r1
                    Um[q, t0] = _cconj(g01) * r0 + _cconj(g11) * r1
                    cuda.syncthreads()
                    # cols p,q <- cols G (thread t owns row t); V <- V G
                    u0 = Um[t0, p]
                    u1 = Um[t0, q]
                    Um[t0, p] = g00 * u0 + g10 * u1
                    Um[t0, q] = g01 * u0 + g11 * u1
                    v0 = Vm[t0, p]
                    v1 = Vm[t0, q]
                    Vm[t0, p] = g00 * v0 + g10 * v1
                    Vm[t0, q] = g01 * v0 + g11 * v1
                    cuda.syncthreads()
            if cuda.threadIdx.x == 0:
                off = 0.0
                for i in range(n):
                    for j in range(n):
                        if i != j:
                            off += _cabs2(Um[i, j])
                sh[0] = math.sqrt(off)
            cuda.syncthreads()
            if sh[0] < tol:
                break

    return jacobi_normal


def gpu_modes(bloch_fn, qs, sweeps=80, tol=1e-10):
    """Batched GPU eigendecomposition of the Bloch step at the points
    ``qs``: returns (thetas, us) where thetas[i] are the TWO smallest
    |theta| eigenphases at q_i (sorted by |theta|) and us[i] the
    corresponding gauge-fixed eigenvectors as (N, 2) columns."""
    from numba import cuda
    jacobi_normal = build_kernels()
    M = len(qs)
    thetas, us = [], []
    CHUNK = 64   # persistent device buffers: 2 x 128x128 complex128 each
    U = np.empty((CHUNK, N, N), dtype=np.complex128)
    EYE = np.array([np.eye(N, dtype=np.complex128)] * CHUNK)
    dU = cuda.to_device(U)
    dV = cuda.to_device(EYE)
    for lo in range(0, M, CHUNK):
        hi = min(lo + CHUNK, M)
        Mc = hi - lo
        for i in range(lo, hi):
            U[i - lo] = bloch_fn(qs[i])
        cuda.to_device(EYE, to=dV)
        cuda.to_device(U, to=dU)
        jacobi_normal[Mc, 128](dU, dV, N, sweeps, tol)
        Ug = dU.copy_to_host()
        Vg = dV.copy_to_host()
        for i in range(Mc):
            th = np.angle(np.diag(Ug[i]))
            order = np.argsort(np.abs(th))[:2]
            U2 = Vg[i][:, order]
            for j in range(2):
                p = int(np.argmax(np.abs(U2[:, j])))
                U2[:, j] = U2[:, j] * np.exp(-1j * np.angle(U2[p, j]))
            thetas.append(np.abs(th[order]))
            us.append(U2)
    return thetas, us


def _sigma_from_modes(th, U2, Uk):
    """One q-point contribution to the Fock sum with the CONVENTION-
    INVARIANT projector: the single smallest-|theta| mode when it is
    isolated, the 2D subspace of the exactly degenerate +-theta pair
    otherwise (the automaton's anti-unitary symmetry pairs the whole
    spectrum; inside a degenerate pair the single-mode choice is an
    arbitrary basis, identical in both routes only as a subspace)."""
    if th[1] - th[0] < _TOL.loose:
        P = U2 @ U2.conj().T
    else:
        P = np.outer(U2[:, 0], U2[:, 0].conj())
    ov = Uk.conj().transpose(0, 2, 1) @ P @ Uk
    return (-1.0 / th[0]) * (ov @ ov.conj().transpose(0, 2, 1))


def main():
    t0 = time.time()
    from numba import cuda
    if not cuda.is_available():
        print("no CUDA device visible - free the GPU and rerun run_altgpu.sh")
        return 1
    C = A.cell((2, 2, 2))
    ST = A.stages_baseline()

    def bf(q):
        return A.bloch_cell(C, SC.sectors()[2], SC.G10, ST, q)

    QSTAR = (math.pi / 2, math.pi / 2, 0.0)

    # 1. kernel verification: random unitaries vs numpy eig
    rng = np.random.default_rng(7)
    M0 = 8
    Aq = rng.standard_normal((M0, N, N)) + 1j * rng.standard_normal((M0, N, N))
    Us = np.array([np.linalg.qr(m)[0] for m in Aq])
    jacobi_normal = build_kernels()
    dU = cuda.to_device(Us.copy())
    dV = cuda.to_device(np.array([np.eye(N, dtype=np.complex128)] * M0))
    jacobi_normal[M0, 128](dU, dV, N, 80, 1e-11)
    Ug = dU.copy_to_host()
    Vg = dV.copy_to_host()
    err_ev, err_pr = 0.0, 0.0
    for i in range(M0):
        th_gpu = np.sort(np.angle(np.diag(Ug[i])))
        th_np = np.sort(np.angle(np.linalg.eigvals(Us[i])))
        err_ev = max(err_ev, float(np.abs(th_gpu - th_np).max()))
        j = int(np.argmin(np.abs(np.angle(np.diag(Ug[i])))))
        Pg = np.outer(Vg[i][:, j], Vg[i][:, j].conj())
        z, W = np.linalg.eig(Us[i])
        k = int(np.argmin(np.abs(np.angle(z))))
        Pn = np.outer(W[:, k], W[:, k].conj())
        err_pr = max(err_pr, float(np.abs(Pg - Pn).max()))
    print(f"altgpu anchor: Jacobi-normal vs numpy eig: eigenphase err "
          f"{err_ev:.2e}, projector err {err_pr:.2e}")

    # 2. the V172 grid pass on the GPU vs numpy on the ACTUAL Bloch
    #    matrices (grid 16, 12 kappas), subspace-projector convention
    dirs = A.sphere_dirs(3)
    kappas = [tuple(1e-3 * x for x in d) for d in dirs[:12]]
    grid_n = 16
    g = np.arange(grid_n) * (2 * math.pi / grid_n) - math.pi
    qpts = []
    for qx in g:
        for qy in g:
            for qz in g:
                q = (float(qx), float(qy), float(qz))
                if math.sqrt(sum((a - b) ** 2
                                 for a, b in zip(q, QSTAR))) >= 0.15:
                    qpts.append(q)
    th_g, us_g = gpu_modes(bf, qpts)
    Uk = np.stack([A.node_modes(bf, tuple(a + b for a, b in zip(QSTAR, kap)),
                                n=1)[0] for kap in kappas])
    S_g = np.zeros((len(kappas), 1, 1), complex)
    S_n = np.zeros_like(S_g)
    amb = 0
    for q, th, U2 in zip(qpts, th_g, us_g):
        if th[1] - th[0] < _TOL.loose:
            amb += 1
        S_g += _sigma_from_modes(th, U2, Uk)
        z, V = np.linalg.eig(bf(q))
        tha = np.abs(np.angle(z))
        order = np.argsort(tha)[:2]
        U2n = V[:, order]
        S_n += _sigma_from_modes(tha[order], U2n, Uk)
    S_g *= -1.0 / len(qpts)
    S_n *= -1.0 / len(qpts)
    dev = max(float(np.abs(sg - sn).max()) for sg, sn in zip(S_g, S_n))
    print(f"altgpu anchor: Fock Sigma GPU vs numpy (grid {grid_n}, "
          f"subspace convention): max dev {dev:.2e}  "
          f"[{amb}/{len(qpts)} points with exact +-theta pair]")
    ok = err_ev < 1e-8 and err_pr < 1e-7 and dev < _TOL.loose
    print(f"GPU route: {'CONFIRMED' if ok else 'MISMATCH'}   "
          f"wall {time.time() - t0:.1f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
