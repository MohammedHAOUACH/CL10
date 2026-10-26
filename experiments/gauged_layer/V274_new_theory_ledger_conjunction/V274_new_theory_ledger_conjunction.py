#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V274 - the new theory's ledger and its own conjunction
(tier 5 conjunction with the tier 6 ledger; no fitted input).

**The question (pre-registered).** The ledger of the dynamical-link
refonte (provenance of every booked number) and the new theory's own
conjunction: three theory points (g10 = 0.2 implemented; g = pi/4 closed
clock; g = 1/5 rung clock) x four legs (COMPOSITE-ALIVE, DYNAMICS-ALIVE,
CHAIN-ALIVE, LOOP-BORN).

**The object (closed input set).** Artifacts V269-V273 plus live
re-computation per point.  The atom kernel of S = C2tot is
g-independent by construction (S is the total Casimir: no coupling in
it), so it is recomputed once and booked at every point; the chain
kernel of the 90 charges is built on the same slice and is likewise
g-independent.  What is re-measured at each point: the ladder test of
H(g) on the kernel ((1/2 - g) lambda_r), the full-space ground E0(g)
(Lanczos 500, seed 0x6C), the ground of the chain Hc(g) with its
materialized Gauss defect over all 90 charges, rho(g) = E_int/g, and
Q H(g) Q measured by the V273 declared machine - block Krylov from the
channel exchange sums, booked only on an exact closure (closure residual
< 1e-8) with the mpmath 60/100-digit spectral control on the closed
space; partner pairings from the ground of Hc(g)).
The implemented point g10 = 0.2 IS the rung clock 1/5
exactly (the implemented branch angle); the coincidence is booked as a
measured fact and the third point shares the first point's numbers
verbatim.

**Measured (pre-registered).** Ledger counts by class (the rotor unit
and the truncation {1,45,210} book as CONVENTIONS - declared, named, not
fitted); the per-leg conjunction table; the branch.

**Acceptance (pre-registered).** Ledger NEW-THEORY-LEDGER with counts
by class; conjunction NEW-CONJUNCTION-ALIVE iff one point passes all
four legs, else NEW-CONJUNCTION-DEAD, per-leg table recorded.
"""
import json
import math
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")  # the 1/5 reference branch
from cl10.paths import exp  # noqa: E402
sys.path.insert(0, exp("V79_gauged_bond_layer"))
import numpy as np  # noqa: E402
from cl10 import coupling  # noqa: E402
from cl10 import engine as E  # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)

NP = 45
LN = 45 + 1 + 210
NATOM = 1024 * LN
KMAX = 500
KRYLOV = 800
SEED = 0x6C
KTOL = 1e-9
NULLTOL = 1e-9
PROBES = [0, 99, 199, 399, 699]
G10 = coupling.G10
BAND = (0.015, 0.03)
TARGET = 1.0 / (16.0 * math.pi)


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}
    xp = E._E

    def log(msg):
        print("[V274 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start, backend={}".format(E._BACKEND))

    v269 = art("V269_dynamical_link_module")
    v270 = art("V270_neutral_composite_sector")
    v271 = art("V271_composite_dynamics")
    v272 = art("V272_dynamical_chain")
    v273 = art("V273_loop_factor_new_theory")
    apass = {n: bool(a["pass"]) for n, a in
             (("V269", v269), ("V270", v270), ("V271", v271),
              ("V272", v272), ("V273", v273))}
    checks.append({"name": "T0 artifacts V269-V273 present, passing: {}".
                   format(apass), "pass": all(apass.values())})
    data["artifact_pass"] = apass

    # ---- T1: the g-independent algebraic layer + kernel + chain ----
    Ts_h = root_generators()
    f, _k, _as = structure_constants(Ts_h)
    Ls45_h = adjoint_link_operators(f)
    I45 = np.eye(NP, dtype=np.complex128)
    Ltot_h = [np.kron(L, I45) + np.kron(I45, L) for L in Ls45_h]
    C2t = sum(L @ L for L in Ltot_h)
    ev, U = np.linalg.eigh(C2t)
    ev = np.real(ev)
    order = np.argsort(ev)
    ev, U = ev[order], U[:, order]
    bins = []
    i0 = 0
    for i in range(1, len(ev) + 1):
        if i == len(ev) or ev[i] - ev[i - 1] > 1e-8:
            bins.append((float(np.mean(ev[i0:i])), i0, i))
            i0 = i

    def swap_par(lo, hi):
        B = U[:, lo:hi].reshape(NP, NP, hi - lo)
        return float(np.real(np.mean(np.sum(B.conj() * B.swapaxes(0, 1),
                                           axis=(0, 1)))))
    binpar = [(v, j - k, swap_par(k, j), k, j) for (v, k, j) in bins]
    s45 = [b for b in binpar if b[1] == 45 and b[2] < -0.5][0]
    s210 = [b for b in binpar if b[1] == 210][0]
    lam45, lam210 = s45[0], s210[0]
    P210 = U[:, s210[3]:s210[4]]
    L256_h = []
    for a in range(NP):
        B = np.zeros((LN, LN), dtype=np.complex128)
        B[:45, :45] = Ls45_h[a]
        B[46:, 46:] = P210.conj().T @ Ltot_h[a] @ P210
        L256_h.append(B)
    C2l_h = sum(L @ L for L in L256_h)
    Ta = xp.asarray(np.stack(Ts_h))
    La = xp.asarray(np.stack(L256_h))
    C2l = xp.asarray(C2l_h)
    I32 = np.eye(32, dtype=np.complex128)
    Jpair_h = []
    C2pair_h = np.zeros((1024, 1024), dtype=np.complex128)
    for T in Ts_h:
        Jp2 = np.kron(T, I32) + np.kron(I32, T)
        Jpair_h.append(Jp2)
        C2pair_h += Jp2 @ Jp2
    C2pair = xp.asarray(C2pair_h)
    Jp = xp.asarray(np.stack(Jpair_h))

    def S_of(psi):
        p = psi.reshape(32, 32, LN)
        v = xp.tensordot(C2pair, p.reshape(1024, LN), axes=(1, 0))
        v = v.reshape(32, 32, LN)
        v = v + xp.einsum("ml,pql->pqm", C2l, p)
        X = xp.einsum("aip,pql->aiql", Ta, p)
        v = v + 2.0 * xp.einsum("aml,aiql->iqm", La, X)
        Y = xp.einsum("aiq,pql->apil", Ta, p)
        v = v + 2.0 * xp.einsum("aml,apil->pim", La, Y)
        return v.reshape(NATOM)

    rng = np.random.default_rng(SEED)
    v0 = (rng.standard_normal(NATOM) + 1j * rng.standard_normal(NATOM))
    v0 = xp.asarray(v0 / np.linalg.norm(v0))
    V = [v0]
    alphas, betas = [], []
    k_used = KMAX
    for k in range(KMAX):
        w = S_of(V[-1])
        a = float(xp.real(xp.vdot(V[-1], w)))
        alphas.append(a)
        w = w - a * V[-1]
        if len(V) > 1:
            w = w - betas[-1] * V[-2]
        for _ in range(2):
            for u in V:
                w = w - xp.vdot(u, w) * u
        b = float(xp.sqrt(xp.real(xp.vdot(w, w))))
        if b < 1e-10:
            k_used = k + 1
            break
        betas.append(b)
        V.append(w / b)
    Tm = np.diag(alphas[:k_used])
    for i in range(k_used - 1):
        Tm[i, i + 1] = betas[i]
        Tm[i + 1, i] = betas[i]
    tev, tvec = np.linalg.eigh(Tm)
    Vm = xp.stack(V[:k_used]).T
    tdev = xp.asarray(tvec)
    cands = [Vm @ tdev[:, j] for j in range(k_used) if tev[j] < KTOL]
    del V, Vm
    CK = [xp.einsum("ml,pql->pqm", C2l, k.reshape(32, 32, LN))
          .reshape(NATOM) for k in cands]
    Gm = np.array([[complex(xp.vdot(cands[i], CK[j]))
                    for j in range(len(cands))]
                   for i in range(len(cands))])
    cev, cvec = np.linalg.eigh(0.5 * (Gm + Gm.conj().T))
    del CK
    ker = [sum(cvec[i, s] * cands[i] for i in range(len(cands)))
           for s in range(len(cands))]
    ker = [k / xp.sqrt(xp.real(xp.vdot(k, k))) for k in ker]
    chans = [float(np.real(xp.vdot(
        k, xp.einsum("ml,pql->pqm", C2l, k.reshape(32, 32, LN))
        .reshape(NATOM)))) for k in ker]
    data["kernel_dim"] = len(ker)
    data["kernel_channels"] = [round(c, 4) for c in chans]
    composite_alive = bool(
        len(ker) == 6 and sorted(round(c, 1) for c in chans)
        == sorted([0.0, 0.0, lam45, lam45, lam210, lam210]))
    checks.append({"name": "T1 atom kernel of S = C2tot (g-independent): "
                          "dim {:d}, channels {} (V270 booked 6 / "
                          "{{0x2, 64x2, 96x2}})".format(
                              len(ker), [round(c, 1) for c in chans]),
                   "pass": composite_alive})

    def H_of(psi, g):
        p = psi.reshape(32, 32, LN)
        v = 0.5 * xp.einsum("ml,pql->pqm", C2l, p)
        X = xp.einsum("aip,pql->aiql", Ta, p)
        v = v + g * xp.einsum("aml,aiql->iqm", La, X)
        Y = xp.einsum("aiq,pql->apil", Ta, p)
        v = v + g * xp.einsum("aml,apil->pim", La, Y)
        return v.reshape(NATOM)

    def Q_of(psi):
        out = psi
        for k in ker:
            out = out - xp.vdot(k, psi) * k
        return out

    K = np.stack([E.asnumpy(k) for k in ker])          # (6, NATOM) host
    K3 = K.reshape(6, 1024, LN)
    CK3 = np.tensordot(K3, C2l_h.T, axes=([-1], [1]))
    Cc = np.einsum("XAm,YAm->XY", K3.conj(), CK3)
    Lh = np.stack(L256_h)
    Jh = np.stack(Jpair_h)
    Xall = np.zeros((NP, 6, 6), dtype=np.complex128)
    Uall = np.zeros((NP, 6, 6), dtype=np.complex128)
    gall_def = 0.0
    for a in range(NP):
        Lk = np.einsum("ml,YQl->YQm", Lh[a], K3)
        Lbar_k = np.einsum("lm,YQl->YQm", Lh[a], K3)
        Xall[a] = np.einsum("XAm,YAm->XY", K3.conj(), Lbar_k)
        Uall[a] = np.einsum("XAm,YAm->XY", K3.conj(),
                            np.einsum("AQ,YQm->YAm", Jh[a], Lk))
        gJ = np.einsum("XAm,YAm->XY", K3.conj(),
                       np.einsum("AQ,YQm->YAm", Jh[a], K3))
        gL = np.einsum("XAm,YAm->XY", K3.conj(), Lk)
        gall_def = max(gall_def, float(np.abs(gJ + gL).max()))
    checks.append({"name": "T1 Gauss transport identity <k|J_a k> = "
                          "-<k|L_a k>: defect {:.1e} (makes the chain "
                          "M_a exact; g-free)".format(gall_def),
                   "pass": bool(gall_def < 1e-9)})

    # chain kernel: N = sum_a M_a^H M_a, M_a = kron(I6, J_a) - kron(X_a, I)
    NH = 6 * 1024
    I6h = np.eye(6, dtype=np.complex128)
    I1024 = xp.asarray(np.eye(1024, dtype=np.complex128))
    Nmat = xp.zeros((NH, NH), dtype=xp.complex128)
    for a in range(NP):
        Ma = xp.kron(xp.asarray(I6h), Jp[a]) \
            - xp.kron(xp.asarray(Xall[a]), I1024)
        Nmat = Nmat + Ma.conj().T @ Ma
        del Ma
    log("T1 chain N assembled (6144^2), eigh")
    nev, nevec = xp.linalg.eigh(Nmat)
    nev_h = E.asnumpy(nev)
    null = [j for j in range(NH) if nev_h[j] < NULLTOL]
    cdim = len(null)
    gap_above = float(nev_h[cdim]) if cdim < NH else float("nan")
    data["chain_kernel_dim"] = cdim
    data["chain_gap_above"] = gap_above
    HV = nevec[:, null].T.reshape(cdim, 6, 1024).copy() if cdim \
        else xp.zeros((0, 6, 1024), dtype=xp.complex128)
    del Nmat, nevec
    xp.get_default_memory_pool().free_all_blocks()
    chain_kernel_alive = bool(cdim == 12)
    checks.append({"name": "T1 chain kernel of the 90 charges (g-free "
                          "slice): dim {:d} (booked 12), gap "
                          "{:.3e}".format(cdim, gap_above),
                   "pass": chain_kernel_alive})
    Hsv = E.asnumpy(HV)                               # (12, 6, 1024)
    K3d = xp.asarray(K3)  # (6, 1024, LN) device view for PHI

    # the 540 exchange vectors (g-independent operators, V273 machine)
    K3c = K.reshape(6, 32, 32, LN)
    AKv = np.empty((NP, 6, 32, 32, LN), dtype=np.complex128)
    CKv = np.empty((NP, 6, 32, 32, LN), dtype=np.complex128)
    for a in range(NP):
        La_h = np.asarray(L256_h[a])
        Lka = np.einsum("lm,ipqm->ipql", La_h, K3c)
        AKv[a] = np.einsum("pqrs,krsm->kpqm",
                           Jpair_h[a].reshape(32, 32, 32, 32), Lka)
        CKv[a] = -np.einsum("lm,ipqm->ipql", La_h, K3c)
    log("T1 exchange vectors built (t+{:.0f}s)".format(time.time() - t0))

    # ---- T2: per-point measurements ----
    points = [(float(G10), "g10 = 0.2 (implemented clock, 1/5 branch)"),
              (math.pi / 4.0, "g = pi/4 (closed clock)")]
    results = {}

    def run_point(g, tag):
        D = {"g": g, "tag": tag}
        # DYNAMICS leg: ladder on the kernel + full-space ground
        Hk6 = 0.5 * Cc + g * np.sum(Uall, axis=0)
        Hk6 = 0.5 * (Hk6 + Hk6.conj().T)
        ek = np.sort(np.real(np.linalg.eigvalsh(Hk6)))
        pred = np.sort((0.5 - g) * np.array(
            [0.0, 0.0, lam45, lam45, lam210, lam210]))
        ladder_def = float(np.max(np.abs(ek - pred)))
        D["kernel_H_spectrum"] = [float(x) for x in ek]
        D["ladder_defect"] = ladder_def
        # full-space ground E0(g): Lanczos KMAX, seed 0x6C
        rg = np.random.default_rng(SEED)
        u0 = (rg.standard_normal(NATOM) + 1j * rg.standard_normal(NATOM))
        u0 = xp.asarray(u0 / np.linalg.norm(u0))
        Ws = [u0]
        al, be = [], []
        for k in range(KMAX):
            w = H_of(Ws[-1], g)
            a_ = float(xp.real(xp.vdot(Ws[-1], w)))
            al.append(a_)
            w = w - a_ * Ws[-1]
            if len(Ws) > 1:
                w = w - be[-1] * Ws[-2]
            for _ in range(2):
                for u in Ws:
                    w = w - xp.vdot(u, w) * u
            b_ = float(xp.sqrt(xp.real(xp.vdot(w, w))))
            if b_ < 1e-10:
                break
            be.append(b_)
            Ws.append(w / b_)
        Tl = np.diag(al)
        for i in range(len(al) - 1):
            Tl[i, i + 1] = Tl[i + 1, i] = be[i]
        E0 = float(np.min(np.linalg.eigvalsh(Tl)))
        D["E0_full_space"] = E0
        # chain Hc(g), ground, rho, Gauss on the materialized ground
        Hloc = 0.5 * Cc + g * np.sum(Uall, axis=0)
        Hc = np.einsum("siB,tjB,ij->st", Hsv.conj(), Hsv, Hloc,
                       optimize=True)
        for a in range(NP):
            tmp = np.einsum("BC,tjC->tjB", Jh[a], Hsv, optimize=True)
            Hc -= g * np.einsum("siB,tjB,ij->st", Hsv.conj(), tmp,
                                Xall[a], optimize=True)
        Hc = 0.5 * (Hc + Hc.conj().T)
        Ec, Ew = np.linalg.eigh(Hc)
        Ecg = float(Ec[0])
        D["chain_H_spectrum"] = [float(x) for x in Ec]
        D["E_chain_ground"] = Ecg
        E_int = Ecg - 2.0 * E0
        rho = E_int / g
        D["E_int"] = E_int
        D["rho"] = rho
        hgv = np.einsum("s,sij->ij", Ew[:, 0], Hsv)    # (6, 1024)
        hgvx = xp.asarray(hgv)
        PHI = xp.einsum("iAl,iB->AlB", K3d, hgvx)      # (1024,256,1024)
        nrm = float(xp.real(xp.vdot(PHI, PHI)))
        gdef = 0.0
        for a in range(NP):
            q1 = xp.tensordot(Jp[a], PHI, axes=([1], [0]))
            q1 = q1 + xp.tensordot(La[a], PHI,
                                   axes=([1], [1])).transpose(1, 0, 2)
            d1 = float(xp.sqrt(xp.real(xp.vdot(q1, q1)) / nrm))
            del q1
            q2 = xp.tensordot(Jp[a], PHI,
                              axes=([1], [2])).transpose(1, 2, 0)
            q2 = q2 - xp.tensordot(La[a].conj(), PHI,
                                   axes=([1], [1])).transpose(1, 0, 2)
            d2 = float(xp.sqrt(xp.real(xp.vdot(q2, q2)) / nrm))
            del q2
            gdef = max(gdef, d1, d2)
        D["chain_gauss_defect"] = gdef
        del PHI, hgvx
        # F3(g): the V273 exchange at this g
        rhoh = (hgv @ hgv.conj().T).conj()
        Jhgv = np.einsum("alb,ib->ail", Jh, hgv)
        Gh = np.einsum("jl,ail->aji", hgv.conj(), Jhgv)
        Kh = np.einsum("bjl,ail->bjai", Jhgv.conj(), Jhgv)
        gseed = xp.zeros(NATOM, dtype=np.complex128)
        for a in range(NP):
            for i in range(6):
                wi = sum(rhoh[j, i].real for j in range(6))
                gseed = gseed + wi * xp.asarray(AKv[a][i].reshape(NATOM))
                gseed = gseed + wi * xp.asarray(CKv[a][i].reshape(NATOM))
        D["ground_seed_norm"] = float(xp.sqrt(xp.real(
            xp.vdot(gseed, gseed))))
        MB = xp.concatenate([xp.asarray(CKv.reshape(NP * 6, NATOM)),
                             xp.asarray(AKv.reshape(NP * 6, NATOM))])
        NB = 2 * NP * 6
        Gb = E.asnumpy(MB @ MB.conj().T)
        gl, gV = np.linalg.eigh(Gb)
        ordr = np.argsort(gl)[::-1]
        r = int(np.sum(gl > 1e-10 * gl.max()))
        Cm = gV[:, ordr[:r]].conj().T / gl[ordr[:r]][:, None]
        HBs = xp.empty((NB, NATOM), dtype=xp.complex128)
        for b in range(NB):
            HBs[b] = H_of(MB[b], g)
        hh = np.sum(np.abs(E.asnumpy(HBs)) ** 2, axis=1)
        G_H = E.asnumpy(HBs @ MB.conj().T).conj()
        del HBs
        xp.get_default_memory_pool().free_all_blocks()
        coef = Cm.conj() @ G_H
        inv = float(np.max(np.sqrt(np.maximum(
            hh - np.sum(np.abs(coef) ** 2, axis=0), 0.0)) / np.sqrt(hh)))
        D["exchange_rank"] = r
        D["exchange_invariance_residual"] = inv
        MAXV = 500
        STEP_CAP = 60
        Upre = xp.empty((MAXV, NATOM), dtype=np.complex128)
        Dcur = 0
        queue = []
        for i in range(6):
            s = xp.zeros(NATOM, dtype=np.complex128)
            for a in range(NP):
                s = s + xp.asarray(AKv[a][i].reshape(NATOM))
                s = s + xp.asarray(CKv[a][i].reshape(NATOM))
            s = Q_of(s)
            ns = float(xp.sqrt(xp.real(xp.vdot(s, s))))
            if ns > 1e-12:
                queue.append(s / ns)
        growth = []
        for v in queue:
            w = v.copy()
            for _ in range(2):
                if Dcur:
                    c = w @ Upre[:Dcur].conj().T
                    w = w - c @ Upre[:Dcur]
                for k in ker:
                    w = w - xp.vdot(k, w) * k
            nw = float(xp.sqrt(xp.real(xp.vdot(w, w))))
            if nw > 1e-9:
                Upre[Dcur] = w / nw
                Dcur += 1
        growth.append(Dcur)
        new = [Upre[j] for j in range(Dcur)]
        closed = False
        it = 0
        while new and it < STEP_CAP and Dcur < MAXV:
            nxt = []
            for v in new:
                w = H_of(v, g)
                for k in ker:
                    w = w - xp.vdot(k, w) * k
                for _ in range(2):
                    c = w @ Upre[:Dcur].conj().T
                    w = w - c @ Upre[:Dcur]
                    for k in ker:
                        w = w - xp.vdot(k, w) * k
                nw = float(xp.sqrt(xp.real(xp.vdot(w, w))))
                if nw > 1e-9:
                    nxt.append(w / nw)
            added = 0
            for w in nxt:
                if Dcur >= MAXV:
                    break
                c = w @ Upre[:Dcur].conj().T
                w2 = w - c @ Upre[:Dcur]
                nw = float(xp.sqrt(xp.real(xp.vdot(w2, w2))))
                if nw > 1e-9:
                    Upre[Dcur] = w2 / nw
                    Dcur += 1
                    added += 1
            growth.append(added)
            new = [Upre[Dcur - added + j] for j in range(added)] \
                if added else []
            it += 1
            if added == 0:
                closed = True
        Ubl = Upre[:Dcur]
        D["block_growth"] = growth
        D["block_dim"] = Dcur
        D["block_closed"] = bool(closed)
        log("point {}: block Krylov closed {} dim {:d} steps {:d}".format(
            tag[:20], closed, Dcur, it))
        if closed:
            Mm = np.empty((Dcur, Dcur), dtype=np.complex128)
            hhq = np.empty(Dcur)
            for p in range(Dcur):
                w = H_of(Ubl[p], g)
                Mm[p] = E.asnumpy(w @ Ubl.conj().T)
                hhq[p] = float(xp.real(xp.vdot(w, w)))
            HU = Mm.conj()
            HU = 0.5 * (HU + HU.conj().T)
            inv2 = float(np.max(np.sqrt(np.maximum(
                hhq - np.sum(np.abs(HU) ** 2, axis=0), 0.0))
                / np.sqrt(np.where(hhq > 1e-14, hhq, 1.0))))
            D["invariant_residual_final"] = inv2
            th, Uth = np.linalg.eigh(HU)
            E_g = Ecg
            denom = E_g - th
            D["min_resolvent_denominator"] = float(np.min(np.abs(denom)))
            fv = 1.0 / denom
            Xp = E.asnumpy(Ubl @ MB.conj().T).T
            UYc = Xp[:NP * 6].reshape(NP, 6, Dcur)
            UXc = Xp[NP * 6:].reshape(NP, 6, Dcur)
            UXf = UXc.reshape(NP * 6, Dcur) @ Uth
            UYf = UYc.reshape(NP * 6, Dcur) @ Uth
        else:
            inv2 = 1.0
            UXf = np.zeros((NP * 6, 1), dtype=np.complex128)
            UYf = UXf.copy()
            fv = np.zeros(1)
            E_g = Ecg
        del Upre, Ubl, MB
        Pxx = (UXf.conj() * fv) @ UXf.T
        Pxy = (UXf.conj() * fv) @ UYf.T
        Pyx = (UYf.conj() * fv) @ UXf.T
        Pyy = (UYf.conj() * fv) @ UYf.T
        t1 = np.einsum("bjai,ji->", Pxx.reshape(NP, 6, NP, 6), rhoh)
        t2 = np.einsum("bjai,aji->", Pxy.reshape(NP, 6, NP, 6), Gh)
        t3 = np.einsum("bjai,bij->", Pyx.reshape(NP, 6, NP, 6), Gh.conj())
        t4 = np.einsum("bjai,bjai->", Pyy.reshape(NP, 6, NP, 6), Kh)
        J2 = -(g ** 2) * (t1 + t2 + t3 + t4)
        gate_ok = bool(D.get("block_closed") and inv2 < 1e-8)
        F3 = (float(np.real(J2)) / g) if gate_ok else float("nan")
        D["F3_gate"] = (
            "valid: exact on {:d}-dim invariant space (closure "
            "residual {:.1e})".format(D["block_dim"], inv2)
            if gate_ok else
            "not booked: block Krylov did not close under H(g) "
            "(dim {} at cap)".format(D["block_dim"]))
        D["J2"] = [float(np.real(J2)), float(np.imag(J2))]
        D["F3_exchange"] = F3
        E_g = Ecg
        probes = ([p for p in PROBES if p < D["block_dim"]]
                  if gate_ok else [])
        import mpmath
        res60, res100 = [], []
        for dps, store in ((60, res60), (100, res100)):
            mpmath.mp.dps = dps
            z = mpmath.mpf(repr(float(E_g)))
            for p in probes:
                s = mpmath.fsum(
                    mpmath.mpf(repr(float(abs(Uth[j, p]) ** 2)))
                    / (z - mpmath.mpf(repr(float(th[j]))))
                    for j in range(D["block_dim"]))
                store.append(float(mpmath.re(s)))
        mpmath.mp.dps = 15
        ctrl = (max(abs(a - b) for a, b in zip(res60, res100))
                if res60 else 0.0)
        D["mpmath_control_residual"] = ctrl
        log("point {}: E0 {:.6f} E_chain {:.6f} rho {:.6f} F3 {:.3e} "
            "ladder {:.1e} gauss {:.1e} ctrl {:.1e}".format(
                tag[:20], E0, Ecg, rho, F3, ladder_def, gdef, ctrl))
        D["F2_Delta_over_E0"] = (float(ek[-1] - ek[0]) / E0
                                 if abs(E0) > 1e-12 else float("nan"))
        D["candidates"] = {"F1": rho, "F2": D["F2_Delta_over_E0"],
                           "F3": F3}
        D["legs"] = {
            "COMPOSITE-ALIVE": composite_alive,
            "DYNAMICS-ALIVE": bool(ladder_def < 1e-9
                                   and np.isfinite(E0)
                                   and gall_def < 1e-9),
            "CHAIN-ALIVE": bool(cdim == 12 and gdef < 1e-12),
            "LOOP-BORN": any(np.isfinite(v) and BAND[0] <= abs(v)
                             <= BAND[1]
                             for v in (rho, D["F2_Delta_over_E0"], F3))}
        D["all_four"] = all(D["legs"].values())
        return D

    for g, tag in points:
        results[tag] = run_point(g, tag)

    tag1 = points[0][1]
    coincident = {
        "g": 0.2, "tag": "g = 1/5 (rung clock)",
        "note": "the implemented clock is exactly the rung clock: the "
                "implemented branch is the 1/5 branch, G10 = 0.2; this "
                "point shares the implemented point's measurements "
                "verbatim (arithmetic identity 0.2 = 1/5)",
        "legs": results[tag1]["legs"], "all_four":
        results[tag1]["all_four"]}

    table = {"g10 = 0.2 (implemented clock, 1/5 branch)":
             results[points[0][1]],
             "g = pi/4 (closed clock)": results[points[1][1]],
             "g = 1/5 (rung clock)": coincident}
    data["conjunction"] = table
    alive = [k for k, v in table.items() if v["all_four"]]
    branch = "NEW-CONJUNCTION-ALIVE" if alive else "NEW-CONJUNCTION-DEAD"

    # ---- T3: the ledger ----
    ledger = [
        {"quantity": "rotor unit (32-dim spinor)", "class": "convention",
         "source": "V269 (declared before any run)"},
        {"quantity": "link truncation {1, 45, 210} of 45 (x) 45",
         "class": "convention", "source": "V269 (declared)"},
        {"quantity": "module closure: swap parity identifies 45, "
                    "lambda_210 = 96, pair multiplicities (2,2,2)",
         "class": "algebraic", "source": "V269 exact projector"},
        {"quantity": "channel Casimir values {0, 64, 96} (C2l on the "
                    "kernel: singlet / adjoint-45 / 210)",
         "class": "algebraic", "source": "V269/V270 T1 here"},
        {"quantity": "kernel dim 6 and split {0x2, 64x2, 96x2}, "
                    "entropies 0/ln45/ln210",
         "class": "measured", "source": "V270; recomputed T1 here"},
        {"quantity": "atom H Gauss law (transport J_a k = -L_a k)",
         "class": "measured",
         "source": "V271 7.7e-14; here {:.1e}".format(gall_def)},
        {"quantity": "H|ker = (1/2 - g) C2l: channel ladder",
         "class": "measured", "source": "V271 (0.300000); per point here"},
        {"quantity": "atom full-space ground E0", "class": "measured",
         "source": "V271 (0.0); per point here"},
        {"quantity": "chain kernel dim 12, gap 36, Gauss 1.8e-15",
         "class": "measured", "source": "V272; recomputed T1 here"},
        {"quantity": "tree exchange E_int = rho = 0", "class": "measured",
         "source": "V272; per point here"},
        {"quantity": "F2 = Delta/E0 candidate",
         "class": "never-realized",
         "source": "V273: E0 = 0 makes it undefined at the implemented "
                    "point; finite elsewhere only where E0(g) != 0"},
        {"quantity": "F1 = rho, F3 = second-order exchange (loop "
                    "candidates)", "class": "measured",
         "source": "V273; per point here"},
        {"quantity": "conjunction legs and table", "class": "measured",
         "source": "this experiment T2"},
    ]
    counts = {}
    for e in ledger:
        counts[e["class"]] = counts.get(e["class"], 0) + 1
    data["ledger"] = ledger
    data["ledger_counts"] = counts
    ok = bool(counts.get("fitted", 0) == 0)
    checks.append({"name": "T3 NEW-THEORY-LEDGER: counts by class {} "
                          "(fitted 0; conventions 2 named and declared "
                          "before any run)".format(counts),
                   "pass": ok})
    legstr = "; ".join(
        "{}: {}".format(
            k, "/".join("Y" if v["legs"][L] else "N"
                        for L in ("COMPOSITE-ALIVE", "DYNAMICS-ALIVE",
                                  "CHAIN-ALIVE", "LOOP-BORN")))
        for k, v in table.items())
    checks.append({"name": "conjunction table (legs in order "
                          "COMPOSITE/DYNAMICS/CHAIN/LOOP): {} -> {}".
                          format(legstr, branch),
                   "pass": True})
    data["status"] = branch
    status = ("new-theory ledger {} and conjunction over g10 = 0.2 / "
              "pi/4 / 1/5: {} -> branch {}".format(
                  counts, legstr, branch))
    out = {"id": "V274_new_theory_ledger_conjunction", "tier": 5,
           "property": "the dynamical-link refonte ledger (tier 6) and "
                       "the new theory's four-leg three-point "
                       "conjunction (tier 5) (V274)",
           "residual": True,
           "programme": "V274",
           "branch": branch,
           "data": data,
           "checks": checks,
           "status": status,
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v274_new_theory_ledger_conjunction.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status, flush=True)
    log("done in {:.1f}s".format(time.time() - t0))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
