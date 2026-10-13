#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""CL10 alt-interactions GPU route (numba CUDA, exploration-only).

The heavy loop of the interaction census (V172/V173) is the batched
eigendecomposition of the 128x128 Bloch step U(q) over the momentum
grid (13824 / 32768 points).  numba.cuda has no LAPACK, so this module
implements the batched solver as a CUDA kernel in the gtc2017-numba
style: one thread block per matrix, cyclic JACOBI-SCHUR sweeps - a
sequence of 2x2 unitary SIMILARITIES that diagonalize a normal
(unitary) matrix; the eigenphases are the angles of the diagonal, the
eigenvectors accumulate in V.

The module is a VERIFICATION route, not a verdict: ``python3 -m
cl10.altgpu`` recomputes the V172 Fock self-energy on the GPU and
compares with the CPU (numpy) computation; agreement < 1e-6 confirms
the GPU route.  It never touches the master modules.
"""

import math
import sys
import time

import numpy as np

from . import altcell as A
from . import supercell as SC

N = 128   # the 16-sector Bloch dimension (8 sub-sites x 16 internal)


def build_kernels():
    from numba import cuda

    @cuda.jit(device=True)
    def _cabs(z):
        return math.sqrt(z.real * z.real + z.imag * z.imag)

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
    def jacobi_schur(U, V, n, sweeps, tol):
        """Batched cyclic Jacobi-Schur: U (m, n, n) unitary -> diagonal by
        2x2 unitary similarities; V accumulates the eigenvectors."""
        m = cuda.grid(1)
        Um = U[m]
        Vm = V[m]
        for sweep in range(sweeps):
            for p in range(n - 1):
                for q in range(p + 1, n):
                    a = Um[p, p]
                    b = Um[p, q]
                    c = Um[q, p]
                    d = Um[q, q]
                    tr = a + d
                    dif = a - d
                    disc = dif * dif + 4.0 * b * c
                    sr = _csqrt(disc)
                    l1 = 0.5 * (tr - sr)
                    l2 = 0.5 * (tr + sr)
                    if _cabs(b) >= _cabs(c):
                        g00, g10 = b, l1 - a
                        g01, g11 = b, l2 - a
                    else:
                        g00, g10 = l1 - d, c
                        g01, g11 = l2 - d, c
                    n1 = math.sqrt(_cabs2(g00) + _cabs2(g10))
                    n2 = math.sqrt(_cabs2(g01) + _cabs2(g11))
                    if n1 < 1e-300 or n2 < 1e-300:
                        g00, g10, g01, g11 = 1.0 + 0j, 0j, 0j, 1.0 + 0j
                    else:
                        g00 /= n1
                        g10 /= n1
                        g01 /= n2
                        g11 /= n2
                        ip = _cconj(g00) * g01 + _cconj(g10) * g11
                        g01 -= ip * g00
                        g11 -= ip * g10
                        n2 = math.sqrt(_cabs2(g01) + _cabs2(g11))
                        if n2 < 1e-300:
                            g00, g10, g01, g11 = 1.0 + 0j, 0j, 0j, 1.0 + 0j
                        else:
                            g01 /= n2
                            g11 /= n2
                    cuda.syncthreads()
                    t = cuda.threadIdx.x
                    u0 = Um[t, p]
                    u1 = Um[t, q]
                    Um[t, p] = u0 * g00 + u1 * g10
                    Um[t, q] = u0 * g01 + u1 * g11
                    v0 = Vm[t, p]
                    v1 = Vm[t, q]
                    Vm[t, p] = v0 * g00 + v1 * g10
                    Vm[t, q] = v0 * g01 + v1 * g11
                    cuda.syncthreads()
                    r0 = Um[p, t]
                    r1 = Um[q, t]
                    Um[p, t] = _cconj(g00) * r0 + _cconj(g10) * r1
                    Um[q, t] = _cconj(g01) * r0 + _cconj(g11) * r1
                    cuda.syncthreads()
            off = 0.0
            for i in range(n):
                for j in range(n):
                    if i != j:
                        off += _cabs2(Um[i, j])
            if math.sqrt(off) < tol:
                break

    return jacobi_schur


def gpu_modes(bloch_fn, qs, sweeps=60, tol=1e-10):
    """Batched GPU eigendecomposition of the Bloch step at the points
    ``qs``: returns (theta, u) for the SMALLEST-|theta| mode at each
    point (theta in (-pi, pi], u the gauge-fixed eigenvector)."""
    from numba import cuda
    jacobi_schur = build_kernels()
    M = len(qs)
    U = np.empty((M, N, N), dtype=np.complex128)
    for i, q in enumerate(qs):
        U[i] = bloch_fn(q)
    V = np.array([np.eye(N, dtype=np.complex128)] * M)
    dU = cuda.to_device(U)
    dV = cuda.to_device(V)
    jacobi_schur[M, 128](dU, dV, N, sweeps, tol)
    U = dU.copy_to_host()
    V = dV.copy_to_host()
    thetas, us = [], []
    for i in range(M):
        th = np.angle(np.diag(U[i]))
        j = int(np.argmin(np.abs(th)))
        u = V[i][:, j]
        p = int(np.argmax(np.abs(u)))
        u = u * np.exp(-1j * np.angle(u[p]))
        thetas.append(float(th[j]))
        us.append(u)
    return thetas, us


def main():
    t0 = time.time()
    from numba import cuda
    if not cuda.is_available():
        print("no CUDA device visible - free the GPU and rerun run_gpu.sh")
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
    jacobi_schur = build_kernels()
    dU = cuda.to_device(Us.copy())
    dV = cuda.to_device(np.array([np.eye(N, dtype=np.complex128)] * M0))
    jacobi_schur[M0, 128](dU, dV, N, 80, 1e-11)
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
    print(f"altgpu anchor: Jacobi-Schur vs numpy eig: eigenphase err "
          f"{err_ev:.2e}, projector err {err_pr:.2e}")

    # 2. the V172 grid pass on the GPU vs CPU (grid 16, 12 kappas)
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
    S = np.zeros((len(kappas), 1, 1), complex)
    for u, w in zip(us_g, th_g):
        Uq = u[:, None]
        ov = Uk.conj().transpose(0, 2, 1) @ Uq
        S += (-1.0 / abs(w)) * (ov @ ov.conj().transpose(0, 2, 1))
    S *= -1.0 / len(qpts)
    S_cpu = A.fock_self_energy(bf, QSTAR, 1.0, grid_n, kappas)
    dev = max(float(np.abs(sg - sc).max())
              for sg, sc in zip(S, S_cpu))
    print(f"altgpu anchor: Fock Sigma GPU vs CPU (grid {grid_n}): "
          f"max dev {dev:.2e}")
    ok = err_ev < 1e-8 and err_pr < 1e-7 and dev < 1e-6
    print(f"GPU route: {'CONFIRMED' if ok else 'MISMATCH'}   "
          f"wall {time.time() - t0:.1f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
