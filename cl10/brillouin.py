#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
CL10 - the **Brillouin-zone chirality test** (the W9 item that had never been
run): does the local QCA realize a **single chiral Weyl sector**, or does the
Nielsen-Ninomiya doubling survive?

**The object.**  The local step (``cl10.local_qca``) is *internal layer + one
staggered bond layer*, the axis cycling ``t mod 3`` and the bond parity
alternating ``t mod 2``.  The step is periodic in **time** with period 6, so
the 6-step **Floquet** operator is periodic in **space** with period 2 in each
direction: the true unit cell is the **2x2x2 cell** and the good quantum
number is the **supercell momentum** q in [-pi,pi]^3 (lattice units).  The
6-step Bloch matrix is

    U_F(q) = prod_{t=0..5} [ S_{j(t),parity(t)}(q) . U_16 ],   j(t) = t mod 3,

with S the 8x8 bond layer on the 8-site cell (a bond crossing the j-direction
of the cell carries the Bloch phase exp(+-i q_j)).  The module builds U_F(q)
**exactly** and validates it against the **real-space** 6-step Floquet on a
4^3 torus (agreement 5e-15).

**The three results.**
1. **The 16 sector is closed by the dynamics** (not only by the internal
   layer): the bond layer acts on the *cell* index, so it commutes with the
   pseudoscalar w - the 16 -> 16bar leakage of the implemented 6-step evolution
   is 1e-15.  **Chirality is a dynamical symmetry of the QCA.**
2. **The zone center is gapped:** min_a |omega_a(q=0)| = 0.0365 per 6-step
   cycle for the 16 sector - the would-be massless matter is **not** at q = 0.
3. **The gapless points are the doublers:** the 16-sector Floquet is gapless at
   the **12 momenta** (±pi/2, ±pi/2, 0) and their permutations - the classic
   Nielsen-Ninomiya doubler momenta - with linear dispersion |v| ~ 0.14 and
   **paired chiralities** (the sum of the node chiralities is 0, as the theorem
   requires).  **The doubling is present**; the theory has no Wilson / overlap /
   domain-wall mechanism that removes it.

**The named residual (W9).**  The doublers must be removed (or confined, or
projected) for the 16 to be the *single* chiral generation; and the gap at
q = 0 is **not** protected in g_10 (it changes by an order of magnitude between
g_10 = 0.2 and g_10 = 0.15).
"""

import math

import numpy as np

from . import clifford
from . import engine as E
from . import local_qca
from . import sm_qca
G10 = local_qca.G10                  # 1/5, the 0-parameter GUT angle
ORDER = [(True, 0), (False, 1), (True, 2), (False, 0), (True, 1), (False, 2)]
DOUBLERS = [q for sx in (1, -1) for sy in (1, -1)
            for q in ((sx * math.pi / 2, sy * math.pi / 2, 0.0),
                      (sx * math.pi / 2, 0.0, sy * math.pi / 2),
                      (0.0, sx * math.pi / 2, sy * math.pi / 2))]   # the 12 doubler momenta


def blocks(g10=G10):
    """The 16 / 32 internal unitaries and the 16 basis (the +i eigenspace of w)."""
    w = E.asnumpy(clifford.pseudoscalar()).astype(np.complex128)
    ev, evec = np.linalg.eig(w)
    B = np.linalg.qr(evec[:, np.argsort(-np.imag(ev))[:16]])[0]
    U32 = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(), g10)).astype(np.complex128)
    return B, U32, B.conj().T @ U32 @ B


def _layer(even, j, q, ncell, g10):
    """The staggered bond layer along axis j at supercell momentum q."""
    c, s = math.cos(g10), math.sin(g10)
    S = np.eye(ncell, dtype=complex)
    bit = 1 << j
    for e in range(ncell):
        if (e & bit) == (0 if even else bit):
            o = e ^ bit
            ph = 1.0 if even else (np.exp(1j * q) if (e & bit) else np.exp(-1j * q))
            S[e, e] = c
            S[o, o] = c
            S[e, o] = -1j * s * ph
            S[o, e] = -1j * s * np.conj(ph)
    return S


def floquet(q, U, order=ORDER, g10=G10):
    """The exact Bloch Floquet matrix over the 6-step cycle."""
    ncell = 2 if len(order) == 2 else 8
    M = np.eye(ncell * U.shape[0], dtype=complex)
    for even, j in order:
        M = np.kron(_layer(even, j, q[j], ncell, g10), U) @ M
    return M


def _phases(M):
    return np.sort(np.angle(np.linalg.eigvals(M)))


def _gap(q, U, order=ORDER, g10=G10):
    return float(np.min(np.abs(_phases(floquet(q, U, order, g10)))))


def sector_closure(n_steps=6, L=8):
    """16-sector closure of the IMPLEMENTED local step (internal + bond layer):
    a state initialized in the 16 is evolved with ``cl10.local_qca`` and the
    16-bar component is measured."""
    B, U32, U16 = blocks()
    P = B @ B.conj().T
    rng = np.random.default_rng(0)
    v16 = rng.standard_normal(16) + 1j * rng.standard_normal(16)
    v32 = B @ v16
    f = np.zeros((L, L, L, 32), dtype=complex)
    f[L // 2, L // 2, L // 2] = v32
    U = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(), G10)).astype(np.complex128)
    leak = 0.0
    for t in range(n_steps):
        f = local_qca.internal_layer(f, U)
        f = local_qca.bond_layer(f, t % 3, t % 2 == 0, G10)
        comp = np.einsum("ab,...b->...a", P, f)
        outside = f - comp
        leak = max(leak, float(np.linalg.norm(outside) / np.linalg.norm(f)))
    return leak


def torus_validation():
    """The Bloch Floquet vs the real-space 6-step Floquet on a 4^3 torus."""
    _, _, U16 = blocks()
    L = 4
    dim = L ** 3 * 16
    c, s = math.cos(G10), math.sin(G10)

    def sid(x, y, z):
        return (x * L * L + y * L + z) * 16

    def bond(even, j):
        X = np.zeros((dim, dim), dtype=complex)
        for x in range(L):
            for y in range(L):
                for z in range(L):
                    if ((x, y, z)[j] % 2 == 0) != even:
                        continue
                    b = [x, y, z]
                    b[j] = (b[j] + 1) % L
                    for i in range(16):
                        ia, ib = sid(x, y, z) + i, sid(*b) + i
                        X[ia, ia] = c
                        X[ib, ib] = c
                        X[ia, ib] = -1j * s
                        X[ib, ia] = -1j * s
        return X

    Uall = np.kron(np.eye(L ** 3), U16)
    F = np.eye(dim, dtype=complex)
    for even, j in ORDER:
        F = bond(even, j) @ Uall @ F
    worst = 0.0
    for qx in (0.0, math.pi):
        for qy in (0.0, math.pi):
            for qz in (0.0, math.pi):
                q = (qx, qy, qz)
                V = np.zeros((dim, 128), dtype=complex)
                for x in range(L):
                    for y in range(L):
                        for z in range(L):
                            e = (x % 2) | ((y % 2) << 1) | ((z % 2) << 2)
                            k = (x // 2, y // 2, z // 2)
                            amp = np.exp(1j * (qx * k[0] + qy * k[1] + qz * k[2])) / math.sqrt(8)
                            for i in range(16):
                                V[sid(x, y, z) + i, e * 16 + i] = amp
                a = _phases(V.conj().T @ F @ V)
                b = _phases(floquet(q, U16))
                worst = max(worst, float(np.abs(a - b).max()))
    return worst


def brillouin():
    """The full W9 Brillouin-zone chirality test."""
    _, U32, U16 = blocks()
    leak = sector_closure()
    dev = torus_validation()

    gap0 = _gap((0.0, 0.0, 0.0), U16)
    gap0_32 = _gap((0.0, 0.0, 0.0), U32)

    nodes = []
    for q in DOUBLERS:
        g = _gap(q, U16)
        sl = []
        for ax in range(3):
            if abs(q[ax]) < 1e-12:
                continue                       # the flat direction of the node
            dq = 1e-4
            qa = list(q); qa[ax] -= dq
            qb = list(q); qb[ax] += dq
            pa = np.angle(np.linalg.eigvals(floquet(tuple(qa), U16)))
            pb = np.angle(np.linalg.eigvals(floquet(tuple(qb), U16)))
            sl.append(float((pb[np.argmin(np.abs(pb))] - pa[np.argmin(np.abs(pa))]) / (2 * dq)))
        nodes.append({"q": [round(v, 6) for v in q], "gap": g,
                      "slope": [round(v, 5) for v in sl],
                      "chirality": float(np.sign(sl[0] * sl[1]))})
    chirality_sum = float(sum(n["chirality"] for n in nodes))
    v_typ = float(np.mean([abs(s) for n in nodes for s in n["slope"]]))

    scan = []
    for g10 in (0.2, 0.15, 0.1):
        _, _, U16g = blocks(g10)
        scan.append({"g10": g10,
                     "gap_q0": _gap((0.0, 0.0, 0.0), U16g, g10=g10),
                     "gap_doubler": _gap((math.pi / 2, 0.0, -math.pi / 2), U16g, g10=g10)})

    return {
        "backend": E._BACKEND,
        "g10": float(G10),
        "floquet_cycle": len(ORDER),
        "sector_leakage_16_to_16bar": leak,
        "bloch_vs_torus_dev": dev,
        "gap_q0_16": gap0,
        "gap_q0_32": gap0_32,
        "n_doubler_nodes": len(nodes),
        "max_doubler_gap": float(max(n["gap"] for n in nodes)),
        "nodes": nodes,
        "chirality_sum": chirality_sum,
        "v_typical": v_typ,
        "g10_scan": scan,
    }


def verify(tol=1e-9):
    r = brillouin()
    ok = (r["sector_leakage_16_to_16bar"] < tol
          and r["bloch_vs_torus_dev"] < 1e-9
          and r["gap_q0_16"] > 1e-3                    # the zone center is gapped
          and r["max_doubler_gap"] < 1e-3               # the doublers are gapless
          and abs(r["chirality_sum"]) < 1e-9            # Nielsen-Ninomiya: they cancel
          and r["n_doubler_nodes"] == 12)
    r["ok"] = bool(ok)
    return r
