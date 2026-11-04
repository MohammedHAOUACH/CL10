#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V362 - the complete rhombic engine and the elected battery
(user directive 2026-11-01: "tu n'as pas cree un engine complet
et test si Le Dodecaedre Rhombique (12 faces) est l'elue - fait
le"; booked before execution in review/pre_registrations_root.md;
tier 1).

The published engine transplanted VERBATIM in structure onto the
rhombic dodecahedral cell: cell = 12 faces x 3 frames = 36 with
the canonical edge-axis transport (the +60-degree rotation about
the shared-edge axis (d_u x d_v)/sqrt3: rational entries of
denominator 6 certified in exact Fractions, O_h-equivariant
including the improper elements, R_vu = R_uv^{-1}; the V361
table-order transport is booked as a non-equivariant rule
artifact) x the published internal 32 VERBATIM
(cl10.supercell, the V358 pattern): U(g, gin) = (I36 (x) Uin .
U_hop(g) (x) I32)^6 on 1152 states, U_hop = W e^{i g lambda}
W^dagger by the symmetric eigensolve. The twisted-hop gap of the
EQUIVARIANT transport is measured (the probe booked opens a gap
- 0.276-class, no target); the exact zero modes of the closed
cell live in the UNTWISTED hopping, re-derived here from the
V361 certified spectrum. The engine's traced-word ladder
F_n^Rh(g) = sum_j e^{i g n lambda_j}/36 then faces the elected
theory's own acceptance: rung-10 crossings of LOOP on (0,1)
bisected to 1e-9, names by CF q <= 10000 for fifth-near roots,
distance against the published dressed crossing V355.gcross_-
dressed with its window gcross_dev; and the rung-10 reading at
the exact fifth C_R = 3(|F_10^Rh(1/5)|/LOOP)^10 against
V352.COEF_OBS inside the elected dial window V355.read_dev - all
windows read as leaves from the artifacts, no hard-coded
decimals. The BORN fixed-point map (V352.rootB: hoisted hop
lifts, Cartan contractions of the 45-link ladder) is published-
engine specific and booked NON-TRANSFERRED; the clock acceptance
on this engine is the T2/T3 battery. GPU-1 directive: the 48
O_h covariances batched on cupy with the house 3-probe device-
vs-host identity gate.
"""
import json
import os
import sys
from fractions import Fraction
from itertools import product, combinations

os.environ["CL10_ANGLE"] = "0.2"    # house pin BEFORE importing cl10

import numpy as np                  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(HERE))))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "archive", "experiments", "root_rebuild",
                                "V352_oneloop_external_backreaction"))
from cl10.paths import exp          # noqa: E402
from cl10 import pentengine as P    # noqa: E402
import V352_oneloop_external_backreaction as V352  # noqa: E402

DIRS = sorted({d for d in product((-1, 0, 1), repeat=3)
               if sum(map(abs, d)) == 2})
DI = {d: i for i, d in enumerate(DIRS)}
FIFTH = 0.2


def _oh():
    gens = [np.array([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
            np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0]]),
            np.array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]]),
            -np.eye(3, dtype=int)]
    I3 = np.eye(3, dtype=int)
    G, seen, fr = [I3], {tuple(I3.ravel())}, [I3]
    while fr:
        nf = []
        for M in fr:
            for H in gens:
                MH = M @ H
                k = tuple(MH.ravel())
                if k not in seen:
                    seen.add(k)
                    G.append(MH)
                    nf.append(MH)
        fr = nf
    return sorted(G, key=lambda M: tuple(M.ravel()))


def _order(M):
    R = np.eye(3, dtype=int)
    for k in range(1, 13):
        R = R @ M
        if np.array_equal(R, np.eye(3, dtype=int)):
            return k
    return -1


def _transport(u, v):
    """+60 deg rotation about (d_u x d_v)/sqrt3, in exact
    Fractions: R = I/2 + E(e)/2 + e e^T/6 (n = |e|^2 = 3 always
    for adjacent rhombic faces, sin(60)/sqrt3 = 1/2 exactly)."""
    du, dv = DIRS[u], DIRS[v]
    e = (du[1] * dv[2] - du[2] * dv[1],
         du[2] * dv[0] - du[0] * dv[2],
         du[0] * dv[1] - du[1] * dv[0])
    H = [[Fraction(0) for i in range(3)] for j in range(3)]
    for i in range(3):
        H[i][i] += Fraction(1, 2)
    Ee = [[0, -Fraction(e[2], 2), Fraction(e[1], 2)],
          [Fraction(e[2], 2), 0, -Fraction(e[0], 2)],
          [-Fraction(e[1], 2), Fraction(e[0], 2), 0]]
    for i in range(3):
        for j in range(3):
            H[i][j] += Ee[i][j] + Fraction(e[i] * e[j], 6)
    return H, e


def main():
    checks, data = [], {}
    G = _oh()
    census = {}
    for M in G:
        census[_order(M)] = census.get(_order(M), 0) + 1
    perms = [tuple(int(DI[tuple(M @ np.array(d))]) for d in DIRS)
             for M in G]
    D = np.array(DIRS)
    adj = [(i, j) for i, j in combinations(range(12), 2)
           if int(D[i] @ D[j]) == 1]

    # ---- T0a: the group exactly ----
    checks.append(("T0a O_h rebuilt exactly: order %d (booked 48), "
                   "census %s (booked {1:1,2:19,3:8,4:12,6:8}), 12-face "
                   "action %d distinct kernel %d"
                   % (len(G), {k: census[k] for k in sorted(census)},
                      len(set(perms)),
                      sum(1 for p in perms if p == tuple(range(12)))),
                   len(G) == 48 and len(set(perms)) == 48
                   and sum(1 for p in perms if p == tuple(range(12)))
                   == 1))

    # ---- T0b: transport + hop in exact Fractions ----
    ok_trans = True
    Kex = {}
    trK2 = Fraction(0)
    for u, v in adj:
        H, e = _transport(u, v)
        # transport maps d_u -> d_v exactly
        mapped = [sum(H[a][b] * DIRS[u][b] for b in range(3))
                  for a in range(3)]
        if tuple(mapped) != DIRS[v]:
            ok_trans = False
        # orthogonality exactly in Fractions
        for a in range(3):
            for b in range(3):
                s = sum(H[a][c] * H[b][c] for c in range(3))
                if s != (Fraction(1) if a == b else Fraction(0)):
                    ok_trans = False
        Kex[(u, v)] = H
        Kex[(v, u)] = [[H[b][a] for b in range(3)] for a in range(3)]
    # tr K^2 over every directed adjacent pair, exact Fractions
    trK2 = sum(Kex[(u, v)][a][c] * Kex[(v, u)][c][a]
               for (u, v) in list(adj) + [(b, a) for a, b in adj]
               for a in range(3) for c in range(3))
    # equivariance EXACTLY (Fractions): h K h^T == K for all 48
    ok_equiv = True
    for h in G:
        for (u, v) in adj:
            hu, hv = DI[tuple(h @ np.array(DIRS[u]))], \
                DI[tuple(h @ np.array(DIRS[v]))]
            L = Kex[(hu, hv)]
            R = Kex[(u, v)]
            # (h R h^T)[a][b]
            for a in range(3):
                for b in range(3):
                    s = sum(h[a][c] * R[c][d] * h[b][d]
                            for c in range(3) for d in range(3))
                    if s != L[a][b]:
                        ok_equiv = False
    Kf = np.zeros((36, 36))
    for (u, v), H in Kex.items():
        Kf[3 * u:3 * u + 3, 3 * v:3 * v + 3] = \
            np.array([[float(x) for x in row] for row in H])
    lam, W = np.linalg.eigh(Kf)
    hopgap = float(np.min(np.abs(lam)))
    data["hop_equivariant"] = {"trK2_exact": int(trK2),
                               "transport_rule": "+-60 deg about "
                               "(d_u x d_v)/sqrt3, Fractions d=6",
                               "equivariant_exact_Fractions": ok_equiv,
                               "RHBIC_HOPGAP_equivariant": hopgap,
                               "spectrum": [float(x) for x in lam]}
    # untwisted closed-cell zero modes re-derived (V361 T1 leaf)
    A = np.zeros((12, 12))
    for i, j in adj:
        A[i, j] = A[j, i] = 1
    wA = np.linalg.eigvalsh(A)
    nzA = int(np.sum(np.abs(wA) < 1e-9))
    checks.append(("T0b the canonical edge-axis transport in exact "
                   "Fractions: maps every adjacent pair (d_u -> d_v), "
                   "orthogonality exact, O_h-equivariance EXACT in "
                   "Fractions for all 48 elements: %s; tr K^2 = %s "
                   "(booked 144, exact); TWISTED hop gap (equivariant "
                   "rule) min |lambda| = %.6f (no target booked - the "
                   "probe booked it opens the gap the rule-artifact "
                   "zero modes of V361-T3 had); untwisted hopping "
                   "re-derived with %d exact zero modes (V361 T1 "
                   "spectrum leaf)"
                   % (ok_trans and ok_equiv, trK2, hopgap, nzA),
                   ok_trans and ok_equiv and trK2 == 144 and nzA == 3))

    # ---- T0c: internal identities verbatim ----
    from cl10 import supercell
    B, U32, U16 = supercell.sectors(0.2)
    dev_int = float(np.max(np.abs(U32 @ B - B @ U16)))
    det_err = abs(np.linalg.det(U32) - 1.0)
    checks.append(("T0c internal identities from supercell.sectors("
                   "0.2): |U32 B - B U16| = %.1e, |det U32 - 1| = %.1e "
                   "(both <= 1e-12)" % (dev_int, det_err),
                   dev_int <= 1e-12 and det_err <= 1e-12))

    # ---- T1 engine gates ----
    K_eig = (lam, W)
    devU = max(float(np.max(np.abs(P.hop_unitary(g, K_eig=K_eig)
                                   @ P.hop_unitary(g, K_eig=K_eig)
                                   .conj().T - np.eye(36))))
               for g in (0.02, FIFTH, np.pi / 4, 0.98))
    detU = abs(np.linalg.det(P.hop_unitary(FIFTH, K_eig=K_eig)) - 1.0)
    U = P.step(FIFTH, 0.2, K=Kf, K_eig=K_eig)
    devS = float(np.max(np.abs(U @ U.conj().T - np.eye(1152))))
    import cupy
    dev = "cupy " + str(cupy.cuda.runtime.getDeviceProperties(
        cupy.cuda.Device().id)["name"])
    Cs = []
    for h in G:
        Pm = np.zeros((12, 12))
        for i, d in enumerate(DIRS):
            Pm[DI[tuple(h @ np.array(d))], i] = 1.0
        C36 = np.kron(Pm, h.astype(float))
        Cs.append(np.kron(C36, np.eye(32)))
    Ucu = cupy.asarray(U)
    devs_cov = []
    for i in range(0, 48, 8):            # batches of 8 on GPU 1
        blk = cupy.stack([cupy.asarray(C) for C in Cs[i:i + 8]])
        devs_cov.extend(float(d) for d in cupy.max(
            cupy.abs(blk @ Ucu @ blk.transpose(0, 2, 1) - Ucu),
            axis=(1, 2)).tolist())
    cov_dev = max(devs_cov)
    gpu_id_dev = 0.0
    for probe in (0, 17, 47):            # host identity gate
        host_dev = float(np.max(np.abs(Cs[probe] @ U @ Cs[probe].T
                                       - U)))
        gpu_id_dev = max(gpu_id_dev, abs(host_dev - devs_cov[probe]))
    data["covariance"] = {"device": dev, "worst_cov_dev": cov_dev,
                          "mean_cov_dev": sum(devs_cov) / len(devs_cov),
                          "identity_gate": gpu_id_dev}
    checks.append(("T1 engine gates: U_hop unitary at g in {0.02, 1/5, "
                   "pi/4, 0.98}: dev %.1e; |det U_hop(1/5)| - 1 = %.1e; "
                   "the full 1152-dim rhombic engine step at (1/5, 0.2) "
                   "unitary dev %.1e (<= 1e-10); O_h covariance over ALL "
                   "48 conjugations, batched on %s: worst %.1e (<= "
                   "1e-10), device-vs-host identity probes %.1e (<= "
                   "1e-12)" % (devU, detU, devS, dev, cov_dev,
                               gpu_id_dev),
                   devU <= 1e-12 and detU <= 1e-10 and devS <= 1e-10
                   and cov_dev <= 1e-10 and gpu_id_dev <= 1e-12))

    # ---- T2 crossings of the rhombic ladder ----
    LOOP = float(V352.LOOP)
    COEF = float(V352.COEF_OBS)

    def f10(g):
        return abs(np.mean(np.exp(1j * g * 10 * lam))) - LOOP

    gs = np.arange(1, 1000) * 1e-3
    roots = []
    for i in range(len(gs) - 1):
        a, b = gs[i], gs[i + 1]
        fa, fb = f10(a), f10(b)
        if fa == 0.0:
            roots.append(a)
        elif fa * fb < 0:
            lo, hi = a, b
            for _ in range(40):
                mid = 0.5 * (lo + hi)
                if f10(mid) * fa <= 0:
                    hi = mid
                else:
                    lo = mid
            roots.append(0.5 * (lo + hi))
    near = [r for r in roots if abs(r - FIFTH) <= 1e-4]
    from fractions import Fraction as Fr
    named = [(r, Fr(r).limit_denominator(10 ** 4)) for r in near]
    with open(exp("V355_elect_one_dial", "results",
                  "v355_elect_one_dial.json")) as fh:
        bat = json.load(fh)["data"]["battery"]
    gcross = float(bat["gcross_dressed"])
    gwin = float(bat["gcross_dev"])
    best = min(roots, key=lambda r: abs(r - gcross)) if roots else None
    data["crossings"] = {"n_roots": len(roots),
                         "roots": [float(r) for r in roots[:40]],
                         "near_fifth": [[float(r), str(fr)]
                                        for r, fr in named],
                         "gcross_published": gcross,
                         "gcross_window": gwin,
                         "best_root_vs_gcross":
                             abs(best - gcross) if best else 9.9}
    checks.append(("T2 crossings |F_10^Rh(g)| = LOOP on (0,1): %d "
                   "root(s); within +-1e-4 of the fifth: %d named by "
                   "CF q <= 10000; nearest root vs the published "
                   "dressed crossing %.9f: %.2e (elected window %.2e)"
                   % (len(roots), len(named), gcross,
                      abs(best - gcross) if best else 9.9, gwin),
                   True))

    # ---- T3 the reading at the exact fifth ----
    F10f = abs(np.mean(np.exp(1j * FIFTH * 10 * lam)))
    C_R = 3.0 * (F10f / LOOP) ** 10
    rel = abs(C_R - COEF) / COEF
    rwin = float(bat["read_dev"])
    data["reading"] = {"F10_abs_at_fifth": F10f, "C_R": C_R,
                       "COEF_OBS": COEF, "rel_dev": rel,
                       "dial_window_read_dev": rwin}
    checks.append(("T3 READING: C_R = 3(|F_10^Rh(1/5)|/LOOP)^10 = "
                   "%.6g vs COEF_OBS = %.6g: rel dev %.3e (elected "
                   "dial window read_dev leaf %.3e: %s)"
                   % (C_R, COEF, rel, rwin,
                      "INSIDE" if rel <= rwin else "outside"), True))

    machinery = all(c[1] for c in (checks[0], checks[1], checks[2],
                                   checks[3]))
    crossed = bool(best is not None and abs(best - gcross) <= gwin)
    read_ok = rel <= rwin
    if not machinery:
        status, ok = "MACHINERY-DRIFT", False
    elif crossed and read_ok:
        status, ok = "RHOMBIC-ELECTED", True
    else:
        status, ok = "RHOMBIC-NOT-ELECTED", True

    print("\n".join("  %s %s" % ("ok " if c[1] else "FAIL", c[0])
                    for c in checks), flush=True)
    print("[V362 done] %s" % status, flush=True)
    _emit(data, checks, ok, status)
    return 0 if ok else 1


def _emit(data, checks, ok, status):
    out = {"id": "V362_rhombic_engine_elected", "tier": 1,
           "property": "the complete engine on the rhombic "
                       "dodecahedral cell against the elected "
                       "theory's acceptance windows",
           "status": status, "pass": bool(ok), "data": data,
           "checks": [{"name": c[0], "pass": c[1]} for c in checks]}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    with open(os.path.join(HERE, "results",
                           "v362_rhombic_engine_elected.json"),
              "w") as fh:
        json.dump(out, fh, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
