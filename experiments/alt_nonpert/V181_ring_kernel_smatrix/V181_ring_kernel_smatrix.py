#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V181 - exact kernels and S-matrices of the gauged sector on 4-site
systems (alt-nonpert, step 2: the finite-volume side)

**The question.**  V175 built the exact two-fermion S-matrix of the gauged
bond layer on the minimal system (V99/V101: 3 sites, 24-state kernel,
gap 42.09 eta).  The IR question (does the flux sector decouple?) is a
VOLUME question: it needs the spectrum one volume up.  This experiment
builds the exact Gauss-law kernel and the exact S-matrix on the 4-site
OPEN chain (144 states: Is - Ts - Ts - Is) and the 4-site RING (1296
states: four Ts tensors, cyclic link identification) - the V101/V175
machinery generalized one volume up, no truncation, eta = g_10 = 1/5.

**Construction.**  Basis tensors are the V101 cached intertwiners
(Is (2,1024,45), Ts (6,1024,45,45)) - each satisfies its vertex Gauss law
exactly (anchor, computed as A-contracted einsums).  Gram = product of
per-link overlaps (V101's P/Q objects, one factor per link).  H = sum
over links of (matter-current sandwich on the J-site tensor, keeping the
bond slot) x L(bond) x (slot overlap of the other tensor) x deltas on
spectator channel labels - exactly V101's bond-1/bond-2 pattern with
the new Ts-Ts sandwiches UQ (keeps slot 2) / UQV (keeps slot 1).
Route 2: same terms with the L-overlap product contracted first and
elementwise ("xz,xz->", the V175 route-2 convention) and reversed bond
order.  S = W exp(-i eta theta) W^dag from the eigh of the whitened H
(unitarity by construction, verified); whitening route 1 Cholesky vs
route 2 symmetric-orthogonal, eigenphase-SET comparison
(basis-invariant, the V175 convention).

GPU: CL10_GPU=1 runs the heavy contractions (UQ/UQV sandwiches, Gram +
H einsums) on cupy; results return to host for the eigensolvers
(multi-threaded BLAS).  CPU is the reference route (same formulas).

**Anchors.**  The N=3 24-state construction re-run in this file
reproduces the V101 JSON (E0 -5802.3982, deg 6, gap 42.0906) to 1e-8;
every basis tensor has Gauss residual < 1e-8; S unitary < 1e-12.
"""

import importlib.util
import json
import os
import sys
import time
import os

_NT = str(min(16, os.cpu_count() or 1))
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, _NT)
os.environ["CL10_BACKEND"] = "numpy"  # the alt programmes are numpy-native: pin the core  # alt suite numpy-native
os.environ["CL10_ANGLE"] = "0.2"  # the reference-branch clock: every number this experiment quotes or feeds is 1/5-branch; the closed default_branch angle must not leak in

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)

import numpy as np
from scipy.linalg import eigh as geigh, cholesky, solve_triangular

from cl10.paths import exp
from cl10 import supercell as SC

USE_GPU = os.environ.get("CL10_GPU", "0") == "1"
if USE_GPU:
    import cupy as cp

_spec = importlib.util.spec_from_file_location(
    "v79", exp("V79_gauged_bond_layer", "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)
_spec2 = importlib.util.spec_from_file_location(
    "v101", exp("V101_bond_current_ground_state",
                "V101_bond_current_ground_state.py"))
V101 = importlib.util.module_from_spec(_spec2)
_spec2.loader.exec_module(V101)

ETA = SC.G10
CHAN = {0: 1, 1: 1, 2: 45, 3: 45, 4: 210, 5: 210}  # Ts index -> fusion r


def whiten_cholesky(N):
    L = cholesky(N, lower=True)
    return solve_triangular(L, np.eye(N.shape[0]), lower=True)


def whiten_orth(N, tol=1e-9):
    w, V = np.linalg.eigh(N)
    keep = w > tol * w.max()
    rank = int(keep.sum())
    X = (V[:, keep] / np.sqrt(w[keep])) @ V[:, keep].conj().T
    return X, rank


def smatrix_eigh(Ho):
    """S = W exp(-i eta theta) W^dag for Hermitian Ho (exact unitarity)."""
    th, W = np.linalg.eigh(Ho)
    return (W * np.exp(-1j * ETA * th)) @ W.conj().T, np.sort(th)


def main():
    t0 = time.time()
    print("=" * 68)
    print(f"V181 - 4-site gauged kernels and S-matrices "
          f"({'GPU/cupy' if USE_GPU else 'CPU/numpy'} contractions, "
          f"BLAS threads {_NT})")
    print("=" * 68)
    Is, Ts = V101.build_kernel(t0)
    nI, nT = Is.shape[0], Ts.shape[0]
    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls = V79.adjoint_link_operators(f)
    T32all = np.stack(Ts_np)          # (45, 32, 32)
    Lall = np.stack(Ls)               # (45, 45, 45)
    SumL = Lall.sum(axis=0)           # sum_A L^A

    # ---- anchor: Gauss law, A-contracted einsums (V101 pattern) ----
    def gauss_end(I):
        Iv = I.reshape(32, 32, 45)
        r = (np.einsum("Akl,lqa->kqa", T32all, Iv, optimize=True)
             + np.einsum("Aql,kla->kqa", T32all, Iv, optimize=True))
        r = r.reshape(1024, 45) + I @ SumL.T
        return float(np.linalg.norm(r) / np.linalg.norm(I))

    def gauss_mid(T):
        Tv = T.reshape(32, 32, 45, 45)
        r = (np.einsum("Akl,lqab->kqab", T32all, Tv, optimize=True)
             + np.einsum("Aql,klab->kqab", T32all, Tv,
                         optimize=True)).reshape(1024, 45, 45)
        r += np.einsum("Aac,lcb->lab", Lall, T, optimize=True)
        r += np.einsum("Abc,lac->lab", Lall, T, optimize=True)
        return float(np.linalg.norm(r) / np.linalg.norm(T))

    gres = max([gauss_end(Is[x]) for x in range(nI)]
               + [gauss_mid(Ts[t]) for t in range(nT)])
    print(f"  Gauss residual (all basis tensors): {gres:.2e} "
          f"[{time.time() - t0:.0f}s]")

    # ---- shared overlap objects (V101) ----
    P = np.einsum("xpa,ypa->xya", Is.conj(), Is)
    Q = np.einsum("tjab,ujab->tuab", Ts.conj(), Ts)
    V2 = np.einsum("tjab,ujcb->tuac", Ts.conj(), Ts)   # keeps slot 1
    W2 = np.einsum("tjab,ujac->tubc", Ts.conj(), Ts)   # keeps slot 2
    # UJ[A,x,y,alpha,alpha'] = sum_p conj(Is_x[p,alpha]) (J_A Is_y)[p,alpha']
    Js = [np.kron(T32all[A], np.eye(32)) + np.kron(np.eye(32), T32all[A])
          for A in range(45)]
    UJ = np.zeros((45, nI, nI, 45, 45), dtype=np.complex128)
    for A in range(45):
        for y in range(nI):
            JI = Js[A] @ Is[y]
            for x in range(nI):
                UJ[A, x, y] = Is[x].conj().T @ JI

    # Ts-Ts current sandwich UQ: keeps slot 2 (the bond slot of the
    # J-site tensor in every chain/ring term); one einsum per
    # generator; GPU when CL10_GPU=1
    def build_UQ():
        UQa = np.zeros((45, nT, nT, 45, 45), dtype=np.complex128)
        if USE_GPU:
            Ts4g = cp.asarray(Ts.reshape(nT, 32, 32, 45, 45))
            Tsg = cp.asarray(Ts)
            T32g = cp.asarray(T32all)
            for A in range(45):
                Tg = T32g[A]
                TsJ = (cp.einsum("kl,tlqab->tkqab", Tg, Ts4g)
                       + cp.einsum("kl,tplab->tpkab", Tg, Ts4g))
                TsJ = TsJ.reshape(nT, 1024, 45, 45).conj()
                UQa[A] = cp.asnumpy(
                    cp.einsum("tjab,ujac->tubc", TsJ, Tsg))
        else:
            Ts4 = Ts.reshape(nT, 32, 32, 45, 45)
            for A in range(45):
                T = T32all[A]
                TsJ = (np.einsum("kl,tlqab->tkqab", T, Ts4,
                                 optimize=True)
                       + np.einsum("kl,tplab->tpkab", T, Ts4,
                                   optimize=True))
                TsJ = TsJ.reshape(nT, 1024, 45, 45).conj()
                UQa[A] = np.einsum("tjab,ujac->tubc", TsJ, Ts,
                                   optimize=True)
        return UQa

    UQ = build_UQ()
    print(f"  UQ built [{time.time() - t0:.0f}s]")

    def heavy_einsum(sub, *operands):
        """einsum on GPU when CL10_GPU=1 (host operands in, host out)."""
        if USE_GPU:
            g = [cp.asarray(o) for o in operands]
            return cp.asnumpy(cp.einsum(sub, *g, optimize=True))
        return np.einsum(sub, *operands, optimize=True)

    I6 = np.eye(nT)
    I2 = np.eye(nI)
    # ---- N=3 anchor: V101's 24-state construction, vectorized ----
    K = np.einsum("xya,tuab->xytub", P, Q)
    N3 = np.einsum("xytub,zwb->xtzyuw", K, P).reshape(24, 24)
    N3 = (N3 + N3.conj().T) / 2
    # states (a,t,b): index (a*6+t)*2+b -> axes (a,t,b,a2,t2,b2)
    h3a = heavy_einsum("Axybc,Abd,rtdc,zw->xrzytw", UJ, Ls, V2, I2)
    h3b = heavy_einsum("Azwbc,Abd,rtdc,xy->xrzytw", UJ, Ls, W2, I2)
    H3 = (h3a + h3b).reshape(24, 24)
    H3 = (H3 + H3.conj().T) / 2
    ev3 = np.real(geigh(H3, N3)[0])
    with open(os.path.join(exp("V101_bond_current_ground_state", ""),
                           "results",
                           "v101_bond_current_ground_state.json")) as fh:
        v101 = json.load(fh)["data"]
    dev3E = abs(ev3[0] - v101["E0"])
    dev3g = abs(ev3[6] - ev3[0] - v101["gap"])
    print(f"  N=3 anchor vs V101 JSON: E0 dev {dev3E:.2e}, "
          f"gap dev {dev3g:.2e} [{time.time()-t0:.0f}s]")

    rows = [{"system": "N=3 chain (V101 anchor)", "n_states": 24,
             "gram_rank": 24, "E0": float(ev3[0]), "degeneracy": 6,
             "gap": float(ev3[6] - ev3[0]), "H_route_dev": None,
             "S_unitarity": None, "S_phase_dev": None,
             "max_mixing": None}]

    # ---- 4-site OPEN chain: states (a,t,u,v), ns = 144 ----
    ns4 = nI * nT * nT * nI
    N4 = heavy_einsum("xya,trab,usbc,vwc->xtuvyrsw",
                      P, Q, Q, P).reshape(ns4, ns4)
    N4 = (N4 + N4.conj().T) / 2
    evn4 = np.linalg.eigvalsh(N4)
    rank4 = int(np.sum(evn4 > 1e-9 * evn4.max()))
    # route 1: bond order 12, 23, 34; scalar pattern "xz,xy,yz->"
    t12 = heavy_einsum("Axybc,Abd,rtdc,su,vw->xrswytuv",
                       UJ, Ls, V2, I6, I2)
    t23 = heavy_einsum("Artbc,Abd,sudc,xy,vw->xrswytuv",
                       UQ, Ls, V2, I2, I2)
    t34 = heavy_einsum("Awvbc,Abd,sudc,xy,rt->xrswytuv",
                       UJ, Ls, W2, I2, I6)
    H4 = (t12 + t23 + t34).reshape(ns4, ns4)
    H4 = (H4 + H4.conj().T) / 2
    # route 2: L-overlap product first, elementwise, reversed order
    LV_V2 = heavy_einsum("Aqd,rtdp->Artqp", Ls, V2)
    LV_W2 = heavy_einsum("Aqd,sudp->Asupq", Ls, W2)
    LV_UQ = heavy_einsum("Aqd,sudp->Asupq", Ls, V2)
    r34 = heavy_einsum("Awvqp,Asupq,xy,rt->xrswytuv", UJ, LV_W2, I2, I6)
    r23 = heavy_einsum("Artqp,Asupq,xy,vw->xrswytuv", UQ, LV_UQ, I2, I2)
    r12 = heavy_einsum("Axyqp,Artqp,su,vw->xrswytuv", UJ, LV_V2, I6, I2)
    H4b = (r34 + r23 + r12).reshape(ns4, ns4)
    devH4 = float(np.abs(H4 - H4b).max())
    if rank4 == ns4:
        ev4 = np.real(geigh(H4, N4)[0])
    else:
        X4, _ = whiten_orth(N4)
        ev4 = np.linalg.eigvalsh(X4 @ H4 @ X4.conj().T)
    E04 = float(ev4[0])
    deg4 = int(np.sum(np.abs(ev4 - E04) < 1e-8 * max(1.0, abs(ev4[-1]))))
    gap4 = float(ev4[deg4] - E04)
    print(f"  chain 144: rank {rank4}, H dev {devH4:.2e}, E0 {E04:.6f} "
          f"deg {deg4} gap {gap4:.6f} [{time.time()-t0:.0f}s]")

    # ---- 4-site RING: states (t1,t2,t3,t4), ns = 1296 ----
    nsr = nT ** 4
    # ring Gram: explicit BLAS route (greedy einsum path explodes:
    # 1.5 GB intermediates, 77 min hang measured 2026-10-13)
    T1r = np.einsum("xyda,wvcd->xyawvc", Q, Q, optimize=False)
    T2r = np.einsum("rsab,mnbc->rsamnc", Q, Q, optimize=False)
    Mr = (T1r.transpose(0, 1, 3, 4, 2, 5).reshape(nsr, 45 * 45)
          @ T2r.transpose(2, 5, 0, 1, 3, 4).reshape(45 * 45, nsr))
    Nr = Mr.reshape(nT, nT, nT, nT, nT, nT, nT, nT).transpose(
        0, 4, 6, 2, 1, 5, 7, 3).reshape(nsr, nsr)
    Nr = (Nr + Nr.conj().T) / 2
    evnr = np.linalg.eigvalsh(Nr)
    rankr = int(np.sum(evnr > 1e-9 * evnr.max()))
    # ring convention (from the Gram): t1=(delta,alpha), t2=(alpha,
    # beta), t3=(beta,gamma), t4=(gamma,delta); every bond uses UQ on
    # the J-site tensor (bond slot = its slot 2) and V2 on the other
    # (bond slot = its slot 1); J on the lower-index site
    hr12 = heavy_einsum("Axybc,Abf,rtfc,su,wd->xrswytud",
                        UQ, Ls, V2, I6, I6)
    hr23 = heavy_einsum("Artbc,Abf,sufc,xy,wd->xrswytud",
                        UQ, Ls, V2, I6, I6)
    hr34 = heavy_einsum("Asubc,Abf,wvfc,xy,rt->xrswytuv",
                        UQ, Ls, V2, I6, I6)
    hr41 = heavy_einsum("Awvbc,Abf,xyfc,rt,su->xrswytuv",
                        UQ, Ls, V2, I6, I6)
    Hr = (hr12 + hr23 + hr34 + hr41).reshape(nsr, nsr)
    Hr = (Hr + Hr.conj().T) / 2
    # route 2: sandwich-L product first, reversed link order
    UL12 = heavy_einsum("Axybg,Abh->Axygh", UQ, Ls)
    UL23 = heavy_einsum("Artbg,Abh->Artgh", UQ, Ls)
    UL34 = heavy_einsum("Asubg,Abh->Asugh", UQ, Ls)
    UL41 = heavy_einsum("Awvbg,Abh->Awvgh", UQ, Ls)
    q41 = heavy_einsum("Awvgh,xyhg,rt,su->xrswytuv", UL41, V2, I6, I6)
    q34 = heavy_einsum("Asugh,wvhg,xy,rt->xrswytuv", UL34, V2, I6, I6)
    q23 = heavy_einsum("Artgh,suhg,xy,wd->xrswytud", UL23, V2, I6, I6)
    q12 = heavy_einsum("Axygh,rthg,su,wd->xrswytud", UL12, V2, I6, I6)
    Hrb = (q41 + q34 + q23 + q12).reshape(nsr, nsr)
    devHr = float(np.abs(Hr - Hrb).max())
    Xr, _ = whiten_orth(Nr)
    evr = np.linalg.eigvalsh(Xr @ Hr @ Xr.conj().T)
    E0r = float(evr[0])
    degr = int(np.sum(np.abs(evr - E0r) < 1e-8 * max(1.0, abs(evr[-1]))))
    gapr = float(evr[degr] - E0r)
    print(f"  ring 1296: rank {rankr}, H dev {devHr:.2e}, E0 {E0r:.6f} "
          f"deg {degr} gap {gapr:.6f} [{time.time()-t0:.0f}s]")

    # ---- S-matrices (chain + ring): whitening 2-route, unitarity,
    #      eigenphase-set comparison, fusion-channel mixing ----
    def channel_mixing(S, st, slot_pos):
        """max over slots of the S norm leaking out of each channel."""
        out = {}
        for r in (1, 45, 210):
            mx = 0.0
            for pos in slot_pos:
                idx = [i for i in range(S.shape[0])
                       if CHAN[st[i][pos]] == r]
                inb = np.zeros(S.shape[0], bool)
                inb[idx] = True
                mx = max(mx, float(np.linalg.norm(
                    S[np.ix_(idx, np.where(~inb)[0])])))
            out[r] = mx
        return out

    sres = {}
    for label, Hx, Nx, rankx, nsx, slot_pos in (
            ("chain", H4, N4, rank4, ns4, (1, 2)),
            ("ring", Hr, Nr, rankr, nsr, (0, 1, 2, 3))):
        st = ([(a, t, u, v) for a in range(nI) for t in range(nT)
               for u in range(nT) for v in range(nI)] if label == "chain"
              else [(t1, t2, t3, t4) for t1 in range(nT)
                    for t2 in range(nT) for t3 in range(nT)
                    for t4 in range(nT)])
        wN, VN = np.linalg.eigh(Nx)
        keep = wN > 1e-9 * wN.max()
        Y = VN[:, keep] / np.sqrt(wN[keep])            # route 2: eigh rect
        S2, ph2 = smatrix_eigh(Y.conj().T @ Hx @ Y)
        Qr, _ = np.linalg.qr(VN[:, keep])              # route 1: QR+Cholesky
        Mq = Qr.conj().T @ Nx @ Qr
        Y1 = Qr @ solve_triangular(cholesky(Mq, lower=True),
                                   np.eye(rankx), lower=True).conj().T
        S1, ph1 = smatrix_eigh(Y1.conj().T @ Hx @ Y1)
        if rankx == nsx:
            Xsq = whiten_cholesky(Nx)                  # square S: unitarity
        else:                                          # + channel mixing
            Xsq = Y @ Y.conj().T
        Ssq, _ = smatrix_eigh(Xsq @ Hx @ Xsq.conj().T)
        unit = float(np.abs(Ssq @ Ssq.conj().T - np.eye(nsx)).max())
        pdev = float(np.abs(ph1 - ph2).max())
        kap = float(wN.max() / wN[keep].min())
        hnorm = float(np.abs(np.linalg.eigvalsh(
            Y.conj().T @ Hx @ Y)).max())
        # pre-registered phase tolerance: whitening-noise bound
        # 100 * eta * ||H|| * eps * kappa(N) (NOT tuned to the dev)
        tolph = 100.0 * ETA * hnorm * np.finfo(float).eps * kap
        mix = channel_mixing(Ssq, st, slot_pos)
        sres[label] = {"unitarity": unit, "phase_dev": pdev,
                       "phase_tol": tolph, "kappa": kap,
                       "mixing": mix, "phases": [float(p) for p in ph1]}
        print(f"  S({label}): unitarity {unit:.2e}, phase-set dev "
              f"{pdev:.2e} (tol {tolph:.1e}, kappa {kap:.1e}), mixing "
              f"{mix} [{time.time()-t0:.0f}s]")

    rows.append({"system": "N=4 open chain", "n_states": ns4,
                 "gram_rank": rank4, "E0": E04, "degeneracy": deg4,
                 "gap": gap4, "H_route_dev": devH4,
                 "S_unitarity": sres["chain"]["unitarity"],
                 "S_phase_dev": sres["chain"]["phase_dev"],
                 "max_mixing": max(sres["chain"]["mixing"].values())})
    rows.append({"system": "N=4 ring", "n_states": nsr,
                 "gram_rank": rankr, "E0": E0r, "degeneracy": degr,
                 "gap": gapr, "H_route_dev": devHr,
                 "S_unitarity": sres["ring"]["unitarity"],
                 "S_phase_dev": sres["ring"]["phase_dev"],
                 "max_mixing": max(sres["ring"]["mixing"].values())})

    ok = (gres < 1e-8 and dev3E < 1e-8 and dev3g < 1e-8
          and devH4 < 1e-10 and devHr < 1e-10
          and sres["chain"]["unitarity"] < 1e-12
          and sres["ring"]["unitarity"] < 1e-12
          and sres["chain"]["phase_dev"] < sres["chain"]["phase_tol"]
          and sres["ring"]["phase_dev"] < sres["ring"]["phase_tol"]
          and len(rows) == 3)

    verdict = {
        "id": "V181_ring_kernel_smatrix",
        "tier": 4,
        "property": "exact Gauss-law kernel and exact S-matrix of the "
                    "gauged bond layer on the 4-site open chain (144 "
                    "states) and 4-site ring (1296 states), generalized "
                    "from the V99/V101/V175 3-site machinery; N=3 anchor "
                    "reproduces the V101 spectrum",
        "residual": None,
        "pass": bool(ok),
        "status": f"EXACT 4-SITE GAUGED SECTORS (3 rows: N=3 anchor, "
                  f"chain, ring). N=3 anchor vs V101 JSON: E0 dev "
                  f"{dev3E:.1e}, gap dev {dev3g:.1e}. Chain 144 states "
                  f"(Gram rank {rank4}): E0 = {E04:.4f}, deg {deg4}, gap "
                  f"= {gap4:.4f} eta. Ring 1296 states (Gram rank "
                  f"{rankr}): E0 = {E0r:.4f}, deg {degr}, gap = "
                  f"{gapr:.4f} eta. Volume drift gap(4)-gap(3) = "
                  f"{gap4 - v101['gap']:+.4f} (chain), {gapr - v101['gap']:+.4f} "
                  f"(ring). S unitary to {sres['chain']['unitarity']:.1e} "
                  f"(chain) / {sres['ring']['unitarity']:.1e} (ring); "
                  f"whitening 2-route phase-set dev "
                  f"{sres['chain']['phase_dev']:.1e} / "
                  f"{sres['ring']['phase_dev']:.1e}; H 2-route dev "
                  f"{devH4:.1e} / {devHr:.1e}. Fusion-channel mixing "
                  f"(max over slots): chain {sres['chain']['mixing']}, "
                  f"ring {sres['ring']['mixing']}. Gauss residual "
                  f"{gres:.1e}. Route: "
                  f"{'cupy GPU' if USE_GPU else 'numpy CPU'}.",
        "data": {"rows": rows, "v101_anchor": v101,
                 "dev_E0": dev3E, "dev_gap": dev3g,
                 "H_dev_chain": devH4, "H_dev_ring": devHr,
                 "gauss_res": gres,
                 "spectrum_chain": [float(e) for e in ev4],
                 "spectrum_ring": [float(e) for e in evr],
                 "phases_chain": sres["chain"]["phases"],
                 "phases_ring": sres["ring"]["phases"],
                 "phase_tol_chain": sres["chain"]["phase_tol"],
                 "phase_tol_ring": sres["ring"]["phase_tol"],
                 "mixing_chain": {str(k): v for k, v in
                                  sres["chain"]["mixing"].items()},
                 "mixing_ring": {str(k): v for k, v in
                                 sres["ring"]["mixing"].items()},
                 "gpu": USE_GPU},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v181_ring_kernel_smatrix.json"),
              "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}   JSON -> "
          f"{os.path.join(outdir, 'v181_ring_kernel_smatrix.json')}")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
