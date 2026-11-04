#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V275 - does the off-kernel exchange resolvent birth the loop factor?
(tier 2, hard, O8-new.)

**The question (pre-registered).** Does the second-order link-mediated
exchange become nonzero when the resolvent is grown from the FULL
off-kernel exchange span (the 540 individual vectors A_a k_i, C_a k_i,
rank 360, recorded H-non-invariant by V273) instead of the 6 channel
sums whose Krylov closed inside ker H?

**The object (pre-registered).** V273's atom verbatim (NATOM =
1024 x 256 = 262144, channel-adapted kernel ker dim 6 via S_of Lanczos
500 seed 0x6C, H_atom_of, AKv[a][i], CKv[a][i], Q_of). Seeds: the 540
exchange vectors individually, kernel-projected, orthonormalized (rank
recorded, expected 360). Block Krylov under H_atom_of, every step
projected off ker and off the running basis (V273 machinery), cap
MAXV = 1200, STEP_CAP = 40. If closed: F5 = J2/g10 by V273's T3
four-term assembly with the 60/100-digit mpmath spectral control
verbatim; min |E_g - spec(HU)| recorded. Target 1/(16 pi) = 1.9894e-2,
band [0.015, 0.03] (same as V273/V263).

**Measured (pre-registered).** rank of the 540-span (kernel-projected,
SVD tol 1e-9; the check books the measured rank - V273 measured 360
for the unprojected set and the note was an observation, not a gate);
growth table, closure flag, invariant residual; F5 (exact if closed;
truncated Rayleigh-Ritz sandwich if capped, labelled not-bookable);
|16 pi F5 - 1|; mpmath 60-vs-100 control.

**Acceptance (pre-registered).** OFFKERNEL-LOOP-BORN iff closed and F5
inside the band; OFFKERNEL-UNBORN iff closed and outside (closest miss
recorded); OFFKERNEL-GROWS iff still growing at MAXV (growth table is
the evidence, no bookable number). The experiment passes in every
branch except a crash or an atom-kernel block: the physics verdict is
the branch, not the suite gate (the V273 convention).
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
    gpu = E._BACKEND != "numpy"
    core = E._BACKEND

    def log(msg):
        print("[V275 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    if gpu:
        # the block-Krylov Upre block (MAXV x NATOM complex128) is 5.03
        # GiB and the whole run needs ~9.5 GiB of card; if the pinned
        # card cannot host it, pin the entire run to the host core
        # (the V254 pattern: correct but slower)
        free = int(xp.cuda.runtime.memGetInfo()[0])
        if free < (9500 << 20):
            gpu, xp = False, np
            core = "numpy (pinned: GPU free {} MiB)".format(free >> 20)
            log("GPU free {} MiB: pinned to host numpy".format(free >> 20))
    log("start, backend={}, g10 = {}".format(core, G10))

    # ---- T0: booked candidates from the committed artifacts ----
    v271 = art("V271_composite_dynamics")
    v272 = art("V272_dynamical_chain")
    F1 = float(v272["data"]["rho"])
    E0 = float(v271["data"]["full_space_E0"])
    ch = v271["data"]["kernel_H_channel_spectrum"]
    Delta = float(ch[-1]) - float(ch[0])
    F2 = (Delta / E0) if abs(E0) > 1e-12 else float("nan")
    E_g = float(v272["data"]["E_chain_ground"])
    data.update({"F1_rho": F1, "F2_Delta_over_E0": F2, "E0_atom": E0,
                 "Delta_atom": Delta, "E_g_chain_ground": E_g,
                 "atom_branch": v271["branch"],
                 "chain_branch": v272["branch"]})
    checks.append({"name": "T0 booked candidates: F1 = rho = {:.6f}; "
                          "F2 = Delta/E0 = {} (E0 = {:.3e}, Delta = "
                          "{:.6f}: structurally undefined at E0 = "
                          "0)".format(F1, "undefined (E0 = 0)"
                                      if np.isnan(F2)
                                      else "{:.6f}".format(F2),
                                      E0, Delta),
                   "pass": bool(v271["pass"] and v272["pass"])})
    hgv = (np.array(v272["data"]["chain_ground_h_re"])
           + 1j * np.array(v272["data"]["chain_ground_h_im"])
           ).reshape(6, 1024)
    tot = float(np.sum(np.linalg.norm(hgv, axis=1) ** 2))
    sch = np.linalg.svd(hgv, compute_uv=False)
    data["chain_ground_schmidt"] = [float(x) for x in sch]
    checks.append({"name": "T0 chain ground h-slice (6,1024) loaded from "
                          "the V272 artifact: total norm {:.6f} (unit "
                          "1), Schmidt weights {}".format(
                              tot, [round(float(x), 6) for x in sch
                                    if x > 1e-12]),
                   "pass": bool(abs(tot - 1.0) < 1e-9)})

    # ---- T1: atom machinery + channel-adapted kernel (V270 seed path) --
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
    s210 = [b for b in binpar if b[1] == 210][0]
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

    def H_atom_of(psi):
        p = psi.reshape(32, 32, LN)
        v = 0.5 * xp.einsum("ml,pql->pqm", C2l, p)
        X = xp.einsum("aip,pql->aiql", Ta, p)
        v = v + G10 * xp.einsum("aml,aiql->iqm", La, X)
        Y = xp.einsum("aiq,pql->apil", Ta, p)
        v = v + G10 * xp.einsum("aml,apil->pim", La, Y)
        return v.reshape(NATOM)

    def S_of(psi):
        p = psi.reshape(32, 32, LN)
        v = xp.tensordot(xp.asarray(C2pair_h), p.reshape(1024, LN),
                         axes=(1, 0))
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
    ker = [Vm @ tdev[:, j] for j in range(k_used) if tev[j] < KTOL]
    del V, Vm
    CKx = [xp.einsum("ml,pql->pqm", C2l, k.reshape(32, 32, LN))
           .reshape(NATOM) for k in ker]
    Gm = np.array([[complex(xp.vdot(ker[i], CKx[j]))
                    for j in range(len(ker))]
                   for i in range(len(ker))])
    cev, cvec = np.linalg.eigh(0.5 * (Gm + Gm.conj().T))
    del CKx
    ker = [sum(cvec[i, s] * ker[i] for i in range(len(cev)))
           for s in range(len(cev))]
    data["atom_kernel_dim"] = len(ker)
    data["atom_kernel_C2l"] = [float(np.real(xp.vdot(
        k, xp.einsum("ml,pql->pqm", C2l, k.reshape(32, 32, LN))
        .reshape(NATOM)))) for k in ker]
    log("T1 atom kernel dim {:d}, C2l {}".format(
        len(ker), [round(x, 1) for x in data["atom_kernel_C2l"]]))
    checks.append({"name": "T1 atom kernel dim {:d} (V270 seed "
                          "path)".format(len(ker)),
                   "pass": bool(len(ker) == 6)})
    if len(ker) != 6:
        data["status"] = "ATOM-KERNEL-BLOCKED"
        _emit(data, checks, t0, ok=False)
        return 1

    Hker = np.array([[complex(xp.vdot(ker[i], H_atom_of(ker[j])))
                      for j in range(6)] for i in range(6)])
    Hker = 0.5 * (Hker + Hker.conj().T)
    hspec = np.sort(np.real(np.linalg.eigvalsh(Hker)))
    ref = sorted(float(x) for x in v271["data"]["kernel_H_channel_spectrum"])
    hker_def = float(np.max(np.abs(hspec[:len(ref)] - np.array(ref))))
    data["kernel_H_spectrum"] = [float(x) for x in hspec]
    checks.append({"name": "T1 atom-kernel H spectrum {} reproduces the "
                          "booked V271 kernel_H_channel_spectrum (max "
                          "diff {:.1e})".format(
                              [round(float(x), 6) for x in hspec],
                              hker_def),
                   "pass": bool(hker_def < 1e-9)})

    # k-basis (channel order) and the 540 atom exchange vectors.  The
    # algebra is V273's verbatim (A_a = (J_a pair)(L_a link), C_a = -L_a
    # on the kernel); only the residency differs: the host cannot hold
    # the full (45, 6, NATOM) AKv/CKv pair next to the running server,
    # so each vertex block is built on the host and uploaded straight
    # into the device-resident MB stack ([CKv rows | AKv rows]).
    K3 = np.stack([E.asnumpy(k) for k in ker]).reshape(6, 32, 32, LN)
    NB = 2 * NP * 6
    MB = xp.empty((NB, NATOM), dtype=np.complex128)
    GramA = np.zeros((NP, 6, 6), dtype=np.complex128)
    GramC = np.zeros((NP, 6, 6), dtype=np.complex128)
    for a in range(NP):
        La_h = np.asarray(L256_h[a])
        Lka = np.einsum("lm,ipqm->ipql", La_h, K3)
        Aa = np.einsum("pqrs,krsm->kpqm",
                       Jpair_h[a].reshape(32, 32, 32, 32), Lka)
        Ca = -np.einsum("lm,ipqm->ipql", La_h, K3)
        GramA[a] = np.einsum("ipqm,jpqm->ij", K3.conj(), Aa)
        GramC[a] = np.einsum("ipqm,jpqm->ij", K3.conj(), Ca)
        MB[a * 6:(a + 1) * 6] = xp.asarray(Ca.reshape(6, NATOM))
        MB[NP * 6 + a * 6:NP * 6 + (a + 1) * 6] = \
            xp.asarray(Aa.reshape(6, NATOM))
    herm_def = max(
        float(np.abs(GramA - GramA.transpose(0, 2, 1).conj()).max()),
        float(np.abs(GramC - GramC.transpose(0, 2, 1).conj()).max()))
    checks.append({"name": "T1 exchange vectors: A_a = (J_a pair)(L_a "
                          "link) and C_a = -L_a Hermitian on the kernel "
                          "(6x6 Gram defect over 45 vertices) "
                          "{:.1e}".format(herm_def),
                   "pass": bool(herm_def < 1e-9)})
    log("T1 540 exchange vectors built (device-resident MB stack, "
        "{:.1f} GiB device, t+{:.0f}s)".
        format(NB * NATOM * 16 / 2**30, time.time() - t0))

    def Q_of(psi):
        out = psi
        for k in ker:
            out = out - xp.vdot(k, psi) * k
        return out

    # ---- T2: partner pairings + the exchange resolvent algebra ----
    rhoh = (hgv @ hgv.conj().T).conj()  # rhoh[j,i] = <h_j|h_i> (bra side)
    Jstack = np.stack(Jpair_h)  # partner generators on P2 (1024)
    Jhgv = np.einsum("alb,ib->ail", Jstack, hgv)  # (45,6,1024): J_a h_i
    Gh = np.einsum("jl,ail->aji", hgv.conj(), Jhgv)  # <h_j|J_a h_i>
    Kh = np.einsum("bjl,ail->bjai", Jhgv.conj(), Jhgv)
    # the host is the tight resource on this box (the server holds most
    # of RAM): the partner generator stacks (1.5 GiB) are consumed by
    # Gh/Kh above, nothing downstream touches them - free them now
    del Jpair_h, Jstack, Jhgv
    checks.append({"name": "T2 partner pairings: rho_h Hermitian defect "
                          "{:.1e}, G^b and K^{{ba}} = <J_b h_j|J_a h_i> "
                          "built".format(
                              float(np.abs(rhoh - rhoh.conj().T).max())),
                   "pass": bool(np.abs(rhoh - rhoh.conj().T).max() < 1e-12)
                           and bool(np.abs(Kh - Kh.transpose(2, 3, 0, 1)
                                           .conj()).max() < 1e-10)})
    # the chain-ground contracted seed (pre-registered Krylov seed):
    # measured; it vanishes - the tree-level non-response, consistent
    # with E_int = rho = 0 (V272).  The resolvent object is then
    # evaluated on the exchange algebra B = {k_i} u {C_a k_i} u {A_a k_i}:
    # its rank and H-invariance are measured; when invariance holds the
    # resolvent form closes exactly on it (no truncation).
    gseed = xp.zeros(NATOM, dtype=np.complex128)
    for a in range(NP):
        for i in range(6):
            wi = sum(rhoh[j, i].real for j in range(6))
            gseed = gseed + wi * MB[NP * 6 + a * 6 + i]
            gseed = gseed + wi * MB[a * 6 + i]
    gsn = float(xp.sqrt(xp.real(xp.vdot(gseed, gseed))))
    data["ground_contracted_seed_norm"] = gsn
    checks.append({"name": "T2 chain-ground contracted seed (the pre-"
                          "registered Krylov seed): norm {:.3e} - the "
                          "first-order response of the chain ground to "
                          "H_ex vanishes; the object is evaluated on the "
                          "measured exchange algebra below".format(gsn),
                   "pass": True})
    # structural closure test of the exchange space (a measurement, not
    # an assumption): does span{C_a k_i, A_a k_i} close under H?
    Gb = E.asnumpy(MB @ MB.conj().T)
    gl, gV = np.linalg.eigh(Gb)
    ordr = np.argsort(gl)[::-1]
    r = int(np.sum(gl > 1e-10 * gl.max()))
    Cm = gV[:, ordr[:r]].conj().T / gl[ordr[:r]][:, None]
    HBs = xp.empty((NB, NATOM), dtype=xp.complex128)
    for b in range(NB):
        HBs[b] = H_atom_of(MB[b])
    hh = np.sum(np.abs(E.asnumpy(HBs)) ** 2, axis=1)
    G_H = E.asnumpy(HBs @ MB.conj().T).conj()
    del HBs
    if gpu:
        xp.get_default_memory_pool().free_all_blocks()
    coef = Cm.conj() @ G_H
    inv_res = float(np.max(np.sqrt(np.maximum(
        hh - np.sum(np.abs(coef) ** 2, axis=0), 0.0)) / np.sqrt(hh)))
    data["exchange_generators"] = NB
    data["exchange_algebra_rank"] = r
    data["exchange_invariance_residual"] = inv_res
    log("T2 exchange space: {:d} generators, rank {:d}; H-invariance "
        "residual {:.1e}".format(NB, r, inv_res))
    checks.append({"name": "T2 exchange space: {:d} generators (rank "
                          "{:d}); H-invariance residual {:.1e} - {} ".
                          format(NB, r, inv_res,
                                 "closes exactly" if inv_res < 1e-9 else
                                 "does NOT close: the resolvent needs "
                                 "the full space, so the pre-registered "
                                 "Krylov carries the object and its "
                                 "truncation control books the honest "
                                 "error"),
                   "pass": True})
    # the pre-registered machine.  V273 measured two facts: the
    # ground-contraction Krylov seed is exactly null, and the full
    # exchange span does not close under H (residual ~ 1).  Its block
    # Krylov from the SIX per-channel exchange SUMS closed in dim 4
    # inside ker H.  V275 grows the resolvent from the 540 INDIVIDUAL
    # off-kernel exchange vectors themselves: seeds = every C_a k_i and
    # A_a k_i, kernel-projected and orthonormalized against the running
    # basis; growth is V273's (off ker, off basis, every step).  If it
    # closes, the resolvent form is EXACT on the invariant space; if it
    # caps, the growth table is the evidence and the truncated
    # sandwich is labelled not-bookable.
    MAXV = 1200
    STEP_CAP = 40
    # numerical rank of the 540 stacked kernel-projected vectors, tol
    # 1e-9 on the singular values (V273 recorded 360 for the
    # unprojected set - a measurement, not a gate).  The singular
    # spectrum of the stacked set is the square root of the Gram
    # spectrum; the Gram (540x540) is evaluated to keep the host free,
    # identical rank to the direct SVD at this tolerance.
    Mq = MB.copy()
    for kk in range(6):
        cproj = Mq @ ker[kk].conj()
        Mq = Mq - cproj[:, None] * ker[kk][None, :]
    Gq = E.asnumpy(Mq @ Mq.conj().T)
    el = np.linalg.eigvalsh(Gq)
    sv = np.sqrt(np.maximum(el, 0.0))[::-1]
    span_rank = int(np.sum(sv > 1e-9 * sv[0]))
    data["exchange_span_rank"] = span_rank
    data["exchange_span_rank_tol"] = 1e-9
    data["exchange_span_rank_method"] = ("Gram spectrum (sqrt of the "
                                         "540x540 Gram eigenvalues) of "
                                         "the stacked kernel-projected "
                                         "set - the SVD singular values")
    del Mq, Gq
    log("T2 exchange span rank {:d} of 540 (tol 1e-9, t+{:.0f}s)".
        format(span_rank, time.time() - t0))
    checks.append({"name": "T2 exchange span rank {:d} of 540 "
                          "(kernel-projected, tol 1e-9; V273 measured "
                          "360 unprojected - recorded, not a "
                          "gate)".format(span_rank),
                   "pass": True})
    Upre = xp.empty((MAXV, NATOM), dtype=np.complex128)
    Dcur = 0
    growth = []
    for b in range(NB):
        w = xp.array(MB[b])
        for k in ker:
            w = w - xp.vdot(k, w) * k
        nw = float(xp.sqrt(xp.real(xp.vdot(w, w))))
        if nw <= 1e-12:
            continue
        w = w / nw
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
    log("T2 seed block: {:d} orthonormal off-kernel exchange vectors "
        "(t+{:.0f}s)".format(Dcur, time.time() - t0))
    front = [Upre[j] for j in range(Dcur)]
    closed = False
    new = list(front)
    it = 0
    while new and it < STEP_CAP and Dcur < MAXV:
        nxt = []
        for v in new:
            w = H_atom_of(v)
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
    Uall = Upre[:Dcur]
    data["block_growth"] = growth
    data["block_steps"] = it
    data["block_dim"] = Dcur
    data["block_closed"] = bool(closed)
    log("T2 full-span block Krylov: closed {} dim {:d} steps {:d} "
        "growth {}".format(closed, Dcur, it, growth[:14]))
    # Rayleigh-Ritz representation of the resolvent on the candidate
    # basis: EXACT when the space is invariant (closed, residual ~ 0),
    # a truncated sandwich when the growth capped (not-bookable).
    Mm = np.empty((Dcur, Dcur), dtype=np.complex128)
    hhq = np.empty(Dcur)
    for p in range(Dcur):
        w = H_atom_of(Uall[p])
        Mm[p] = E.asnumpy(w @ Uall.conj().T)  # <H u_p|u_q>
        hhq[p] = float(xp.real(xp.vdot(w, w)))
    HU = Mm.conj()
    HU = 0.5 * (HU + HU.conj().T)
    if closed:
        inv2 = float(np.max(np.sqrt(np.maximum(
            hhq - np.sum(np.abs(HU) ** 2, axis=0), 0.0))
            / np.sqrt(np.where(hhq > 1e-14, hhq, 1.0))))
    else:
        inv2 = float("nan")
    data["invariant_residual_final"] = inv2
    checks.append({"name": "T2 full-span block Krylov: closed {} dim "
                          "{:d} steps {:d}".format(closed, Dcur, it),
                   "pass": bool(not closed or inv2 < 1e-8)})
    th, Uth = np.linalg.eigh(HU)
    denom = E_g - th
    data["min_resolvent_denominator"] = float(np.min(np.abs(denom)))
    fv = 1.0 / denom
    Xp = E.asnumpy(Uall @ MB.conj().T).T  # (540, D): <u_m|x_n>
    UYc = Xp[:NP * 6].reshape(NP, 6, Dcur)
    UXc = Xp[NP * 6:].reshape(NP, 6, Dcur)
    UXf = UXc.reshape(NP * 6, Dcur) @ Uth
    UYf = UYc.reshape(NP * 6, Dcur) @ Uth
    del Upre, Uall, MB

    # ---- T3: F5 = J2/g10 assembly + 60/100-digit control ----
    Pxx = (UXf.conj() * fv) @ UXf.T
    Pxy = (UXf.conj() * fv) @ UYf.T
    Pyx = (UYf.conj() * fv) @ UXf.T
    Pyy = (UYf.conj() * fv) @ UYf.T
    t1 = np.einsum("bjai,ji->", Pxx.reshape(NP, 6, NP, 6), rhoh)
    t2 = np.einsum("bjai,aji->", Pxy.reshape(NP, 6, NP, 6), Gh)
    t3 = np.einsum("bjai,bij->", Pyx.reshape(NP, 6, NP, 6), Gh.conj())
    t4 = np.einsum("bjai,bjai->", Pyy.reshape(NP, 6, NP, 6), Kh)
    J2 = -(G10 ** 2) * (t1 + t2 + t3 + t4)
    gate_ok = bool(closed and inv2 < 1e-8)
    F5 = (float(np.real(J2)) / G10) if closed else float("nan")
    data["J2"] = [float(np.real(J2)), float(np.imag(J2))]
    data["J2_imag_part"] = float(np.imag(J2))
    data["F5_exchange"] = F5
    data["F5_gate"] = (
        "valid: exact on the {:d}-dim invariant space (closure "
        "residual {:.1e})".format(Dcur, inv2) if gate_ok else
        ("exact assembly on the {:d}-dim closed space, closure "
         "residual {:.1e} >= 1e-8: recorded with that flag".format(
             Dcur, inv2) if closed else
         "not bookable: the exchange resolvent did not close under "
         "block Krylov (growth {} to dim {}); the truncated Rayleigh-"
         "Ritz sandwich is evidence, not a booked number".format(
             growth, Dcur)))
    if not gate_ok:
        data["F5_truncated"] = float(np.real(J2)) / G10
        data["F5_truncated_bookable"] = False
    import mpmath
    probes = [p for p in PROBES if p < Dcur]
    res60, res100 = [], []
    for dps, store in ((60, res60), (100, res100)):
        mpmath.mp.dps = dps
        z = mpmath.mpf(repr(float(E_g)))
        for p in probes:
            s = mpmath.fsum(
                mpmath.mpf(repr(float(abs(Uth[j, p]) ** 2)))
                / (z - mpmath.mpf(repr(float(th[j]))))
                for j in range(Dcur))
            store.append(float(mpmath.re(s)))
    mpmath.mp.dps = 15
    ctrl = (max(abs(a - b) for a, b in zip(res60, res100))
            if res60 else 0.0)
    data["mpmath_probes"] = probes
    data["mpmath_probes_60"] = res60
    data["mpmath_probes_100"] = res100
    data["mpmath_control_residual"] = ctrl
    band_F = F5 if closed else data["F5_truncated"]
    dev = abs(16.0 * math.pi * band_F - 1.0)
    data["closest_band_distance"] = (float(dev) if closed else
                                     float(dev))
    log("T3 F5 = {}; |16 pi F - 1| = {:.3e}; mpmath 60/100 control "
        "residual {:.1e}".format(
            "nan (capped)" if not closed else "{:.6e}".format(F5),
            dev, ctrl))
    checks.append({"name": "T3 F5 = {:.3e} vs 1/(16 pi) = {:.3e}: "
                          "closest |16 pi F - 1| = {:.3e} ({})".format(
                              band_F, TARGET, dev,
                              "exact on the invariant space" if gate_ok
                              else "truncated Rayleigh-Ritz sandwich, "
                              "not bookable"),
                   "pass": True})

    # ---- branch ----
    if not closed:
        branch = "OFFKERNEL-GROWS"
    elif BAND[0] <= abs(F5) <= BAND[1]:
        branch = "OFFKERNEL-LOOP-BORN"
    else:
        branch = "OFFKERNEL-UNBORN"
    data["status"] = branch
    data["band_hit"] = bool(closed and BAND[0] <= abs(F5) <= BAND[1])
    print("V275 {:s}: full-span block Krylov closed {} dim {:d} after "
          "{:d} growth steps (span rank {:d}); F5 = {}; closest to "
          "1/(16 pi) = {:.4e}: |16 pi F - 1| = {:.3e} (band "
          "[0.015, 0.03]); mpmath control {:.1e}.".format(
              branch, closed, Dcur, it, span_rank,
              ("not bookable (truncated sandwich {:.6e})".format(
                  data["F5_truncated"]) if not closed else
               "{:.6e}".format(F5)),
              TARGET, dev, ctrl), flush=True)
    _emit(data, checks, t0, ok=True)
    return 0


def _emit(data, checks, t0, ok):
    out = {"id": "V275_offkernel_exchange_resolvent", "tier": 2,
           "property": "does the full off-kernel exchange resolvent "
                       "birth the loop factor (tier 2, O8-new)",
           "residual": True,
           "programme": "V275",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v275_offkernel_exchange_resolvent.json"),
                           "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V275 done in {:.1f}s] {}".format(time.time() - t0,
                                             out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
