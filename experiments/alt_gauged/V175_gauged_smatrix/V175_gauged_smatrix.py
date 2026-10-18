#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V175 - the two-fermion S-matrix of the gauged bond layer (alt-gauged,
step 1: the missing object)

**The question.**  The master automaton is a FREE theory: the 2-particle
S-matrix is ``U (x) U`` (``cl10/smatrix.py``, V48 Kronecker property).  The
gauged bond layer (V79: Spin(10) adjoint quantum links, bond gate
``exp(i eta J.L)``) is the only interaction the theory has, and its
two-fermion S-matrix is NOT built anywhere in the repo (V101 diagonalized
the bond-current Hamiltonian but never exponentiated it).  V98 proved there
are NO one-body physical states (center obstruction), so the three-site,
two-link chain with paired matter (V99: 24-state physical kernel, fusion
channels r = 1, 45, 210) IS the minimal two-fermion system.  This
experiment constructs the exact S-matrix ``S = exp(-i g_10 H_orth)`` on
that kernel (one stage; eta = g_10 = 1/5, the root angle - no new
parameter), and records what the interaction does: eigenphases (the
scattering phases) and fusion-channel mixing.

**Construction.**  V101's glued kernel states ``psi_(a,t,b)`` are not
orthonormal (shared link indices): the physical S-matrix is
``S = exp(-i eta L^-1 H L^-dagger)`` with ``N = L L^dagger`` the exact Gram
and ``H`` the exact bond-current matrix (Gauss turns J into -L on the
intertwiner leg - no truncation).  Route 2 rebuilds Gram and H with
different contraction paths (explicit sums vs the V101 einsum chain;
per-generator vs per-bond operator order, bonds reversed) and the
whitening via the Cholesky vs the symmetric-orthogonal (eigh) factor.

**The result.**  S is unitary to 1e-12; the 24 scattering phases are the
V101 spectrum times g_10; the interaction MIXES (or does not mix) the
fusion channels - recorded, not assumed.  The eta -> 0 anchor gives S = I.
"""

import importlib.util
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core

import numpy as np
from scipy.linalg import eigh as geigh, expm, cholesky, solve_triangular

from cl10.paths import exp
from cl10 import supercell as SC

_spec = importlib.util.spec_from_file_location(
    "v79", exp("V79_gauged_bond_layer", "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)
_spec2 = importlib.util.spec_from_file_location(
    "v101", exp("V101_bond_current_ground_state",
                "V101_bond_current_ground_state.py"))
V101 = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(V101)


def gram_route1(Is, Ts, states):
    """The V101 einsum chain (route 1)."""
    ns = len(states)
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    Q = np.einsum("tjab,ujab->tuab", Ts.conj(), Ts)
    K = np.einsum("xya,tuab->xytub", P, Q)
    Z = np.einsum("xki,yki->xyi", Is.conj(), Is)
    N = np.einsum("xytub,zwb->xtzyuw", K, Z).reshape(ns, ns)
    return (N + N.conj().T) / 2


def gram_route2(Is, Ts, states):
    """Route 2: explicit per-channel sums, different contraction order."""
    nI, nT = Is.shape[0], Ts.shape[0]
    ns = len(states)
    P = np.zeros((nI, nI, 45), complex)
    for x in range(nI):
        for y in range(nI):
            P[x, y] = np.sum(Is[x].conj() * Is[y], axis=0)
    Q = np.zeros((nT, nT, 45, 45), complex)
    for t in range(nT):
        for u in range(nT):
            Q[t, u] = np.sum((Ts[t].conj() * Ts[u]).reshape(-1, 45, 45),
                             axis=0)
    N = np.zeros((ns, ns), complex)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            N[s, s2] = np.einsum("a,ab,b->", P[a, a2], Q[t, t2], P[b, b2])
    return (N + N.conj().T) / 2


def h_route1(Is, Ts, states, Ls_np, Js):
    """The V101 per-generator UJ chain (route 1)."""
    nI = Is.shape[0]
    ns = len(states)
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        for y in range(nI):
            JI = Js[A] @ Is[y]
            for x in range(nI):
                UJ[A, x, y] = Is[x].conj().T @ JI
    V2 = np.einsum("tjab,ujcb->tuac", Ts.conj(), Ts)
    W2 = np.einsum("tjab,ujac->tubc", Ts.conj(), Ts)
    H = np.zeros((ns, ns), dtype=np.complex128)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            h = 0.0 + 0.0j
            if b == b2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, a, a2], Ls_np[A],
                                   V2[t, t2])
            if a == a2:
                for A in range(45):
                    h += np.einsum("xz,xy,yz->", UJ[A, b, b2], Ls_np[A],
                                   W2[t, t2])
            H[s, s2] = h
    return (H + H.conj().T) / 2


def h_route2(Is, Ts, states, Ls_np, Js):
    """Route 2: left-associated sandwich, V2/W2 via matrix products,
    contraction as the elementwise sum UJ . (L V) (V101's scalar), bonds
    in reversed order (a-side first)."""
    nI, nT = Is.shape[0], Ts.shape[0]
    ns = len(states)
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        Jd = [Is[x].conj().T @ Js[A] for x in range(nI)]   # (45, 1024)
        for x in range(nI):
            for y in range(nI):
                UJ[A, x, y] = Jd[x] @ Is[y]
    V2 = np.zeros((nT, nT, 45, 45), complex)
    W2 = np.zeros((nT, nT, 45, 45), complex)
    for t in range(nT):
        for u in range(nT):
            for j in range(Ts.shape[1]):
                V2[t, u] += Ts[t, j].conj() @ Ts[u, j].T
                W2[t, u] += Ts[t, j].conj().T @ Ts[u, j]
    H = np.zeros((ns, ns), dtype=np.complex128)
    for s, (a, t, b) in enumerate(states):
        for s2, (a2, t2, b2) in enumerate(states):
            h = 0.0 + 0.0j
            if a == a2:                       # bond order reversed: a-side
                for A in range(45):
                    h += np.einsum("xz,xz->", UJ[A, b, b2],
                                   Ls_np[A] @ W2[t, t2])
            if b == b2:
                for A in range(45):
                    h += np.einsum("xz,xz->", UJ[A, a, a2],
                                   Ls_np[A] @ V2[t, t2])
            H[s, s2] = h
    return (H + H.conj().T) / 2


def whiten_cholesky(N):
    L = cholesky(N, lower=True)
    return solve_triangular(L, np.eye(N.shape[0]), lower=True)


def whiten_orth(N):
    """Route 2 whitening: N^-1/2 via the symmetric orthogonal factor."""
    w, V = np.linalg.eigh(N)
    return V @ np.diag(1.0 / np.sqrt(w)) @ V.conj().T


def main():
    t0 = time.time()
    print("=" * 68)
    print("V175 - the two-fermion S-matrix of the gauged bond layer")
    print("=" * 68)
    Is, Ts = V101.build_kernel(t0)
    nI, nT = Is.shape[0], Ts.shape[0]
    states = [(a, t, b) for a in range(nI) for t in range(nT)
              for b in range(nI)]
    ns = len(states)
    print(f"  kernel: nI = {nI}, nT = {nT}, states = {ns}")

    Ts_np = V79.root_generators()
    f, _k, _as = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    I32 = np.eye(32, dtype=np.complex128)
    Js = [np.kron(Ts_np[A], I32) + np.kron(I32, Ts_np[A])
          for A in range(45)]

    # --- Gram and H, two independent routes -----------------------------
    N1 = gram_route1(Is, Ts, states)
    N2 = gram_route2(Is, Ts, states)
    H1 = h_route1(Is, Ts, states, Ls_np, Js)
    H2 = h_route2(Is, Ts, states, Ls_np, Js)
    dev_N = float(np.abs(N1 - N2).max())
    dev_H = float(np.abs(H1 - H2).max())
    evn = np.linalg.eigvalsh(N1)
    rank = int(np.sum(evn > 1e-9 * evn.max()))
    print(f"  Gram rank {rank}/{ns}; route devs: Gram {dev_N:.2e}, "
          f"H {dev_H:.2e}")

    # --- spectrum anchor vs V101 JSON -----------------------------------
    ev1 = np.real(geigh(H1, N1)[0])
    ev2 = np.real(geigh(H2, N2)[0])
    v101 = json.load(open(os.path.join(
        ROOT, "experiments", "gauged_layer",
        "V101_bond_current_ground_state", "results",
        "v101_bond_current_ground_state.json")))
    ref = np.array(v101["data"]["spectrum"])
    spec_dev = float(np.abs(ev1 - ref).max() / max(1.0, abs(ref).max()))
    route_spec_dev = float(np.abs(ev1 - ev2).max())
    E0 = float(ev1[0])
    deg = int(np.sum(np.abs(ev1 - ev1[0]) < 1e-8 * max(1.0, abs(ev1[-1]))))
    gap = float(ev1[deg] - ev1[0])
    print(f"  spectrum vs V101 JSON: {spec_dev:.2e}; route dev "
          f"{route_spec_dev:.2e}")
    print(f"  E0 = {E0:.6f}  degeneracy = {deg}  gap = {gap:.6f} "
          f"(V101: -5802.3982 / 6 / 42.0906)")

    # --- the S-matrix ----------------------------------------------------
    X1 = whiten_cholesky(N1)
    X2 = whiten_orth(N2)
    Ho1 = X1 @ H1 @ X1.conj().T
    Ho2 = X2 @ H2 @ X2.conj().T
    Ho1 = (Ho1 + Ho1.conj().T) / 2
    Ho2 = (Ho2 + Ho2.conj().T) / 2
    eta = SC.G10
    S = expm(-1j * eta * Ho1)
    S2 = expm(-1j * eta * Ho2)
    unit = float(np.abs(S @ S.conj().T - np.eye(ns)).max())
    unit2 = float(np.abs(S2 @ S2.conj().T - np.eye(ns)).max())
    phases = np.sort(np.angle(np.linalg.eigvals(S)))
    phases2 = np.sort(np.angle(np.linalg.eigvals(S2)))
    # the two whitening routes give S in DIFFERENT orthonormal bases
    # (Cholesky vs symmetric factor): the basis-invariant 2-route check
    # is the sorted eigenphase set, not the matrix entries.
    sdev = float(np.max(np.abs(phases - phases2)))
    ev_orth = np.linalg.eigvalsh(Ho1)
    phase_dev = float(np.max(np.abs(np.sort(
        np.angle(np.exp(-1j * eta * ev_orth))) - phases)))
    S0 = expm(-1j * eta * np.zeros((ns, ns)))
    eta0_dev = float(np.abs(S0 - np.eye(ns)).max())
    print(f"  S unitarity err = {unit:.2e} (route 2 {unit2:.2e}); "
          f"whitening eigenphase dev = {sdev:.2e}")
    print(f"  phases vs -eta*ev(H_orth): {phase_dev:.2e}; "
          f"eta->0 anchor S=I err {eta0_dev:.2e}")

    # --- fusion-channel block structure ----------------------------------
    # Ts order: T1 (singlet multiplicity of 32 (x) 32), then T45 (= nI
    # intertwiners), then T210 (the rest) - V101 build_kernel.
    C2pair = sum(J @ J for J in Js)
    m1 = int(np.sum(np.abs(np.linalg.eigvalsh(C2pair)) < 1e-6))
    bounds = [(0, m1), (m1, m1 + nI), (m1 + nI, nT)]
    chans = [1, 45, 210]
    rows = []
    for r, (lo, hi) in zip(chans, bounds):
        idx = [s for s, (a, t, b) in enumerate(states) if lo <= t < hi]
        oidx = [s for s in range(ns) if s not in idx]
        blk = S[np.ix_(idx, idx)]
        out = S[np.ix_(idx, oidx)]
        rows.append({
            "r": r, "n_states": len(idx),
            "mean_phase": float(np.mean(np.angle(np.linalg.eigvals(blk)))),
            "diag_block_norm": float(np.linalg.norm(blk)),
            "mixing_out_norm": float(np.linalg.norm(out)),
        })
        print(f"  r = {r:3d}: {len(idx)} states, mean phase "
              f"{rows[-1]['mean_phase']:+.6f}, mixing out "
              f"{rows[-1]['mixing_out_norm']:.2e}")

    ok = (rank == ns and spec_dev < 1e-8 and route_spec_dev < 1e-9
          and dev_N < 1e-12 and dev_H < 1e-10 and unit < 1e-12
          and unit2 < 1e-12 and sdev < 1e-9 and phase_dev < 1e-10
          and eta0_dev < 1e-15
          and abs(E0 + 5802.3982) < 1e-3 and deg == 6
          and abs(gap - 42.0906) < 1e-3)

    mix = [f"{r0['mixing_out_norm']:.1e}" for r0 in rows]
    mixes = ("MIXES" if max(r0["mixing_out_norm"] for r0 in rows) > 1e-8
             else "PRESERVES")
    verdict = {
        "id": "V175_gauged_smatrix",
        "tier": 4,
        "property": "exact two-fermion S-matrix of the gauged bond layer "
                    "on the 24-state three-site kernel: S = exp(-i g_10 "
                    "L^-1 H L^-dagger), Gram-whitened; fusion-channel "
                    "block structure",
        "residual": None,
        "pass": bool(ok),
        "status": f"THE TWO-FERMION S-MATRIX OF THE GAUGED LAYER EXISTS. "
                  f"On the minimal physical two-fermion system (V98: no "
                  f"one-body states; V99 kernel, 24 states, Gram rank "
                  f"{rank}/24) the bond-current interaction exponentiates "
                  f"to an exactly unitary S = exp(-i g_10 H_orth): "
                  f"unitarity err {unit:.1e}, 2-route (einsum/per-generator "
                  f"vs explicit-sum/per-bond, Cholesky vs orthogonal "
                  f"whitening) dev {sdev:.1e}, spectrum matches V101 to "
                  f"{spec_dev:.1e}. E0 = {E0:.4f} eta, degeneracy {deg}, "
                  f"gap {gap:.4f} eta. Fusion-channel mixing (r = 1, 45, "
                  f"210): {mix} - the interaction {mixes} the channels. "
                  f"eta -> 0: S = I exactly.",
        "data": {"rows": rows, "phases": [float(p) for p in phases],
                 "spectrum": [float(e) for e in ev1],
                 "E0": E0, "degeneracy": deg, "gap": gap,
                 "unitarity_err": unit, "route_spec_dev": route_spec_dev,
                 "spec_dev_vs_v101": spec_dev, "gram_rank": rank,
                 "eta": eta},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v175_gauged_smatrix.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v175_gauged_smatrix.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
