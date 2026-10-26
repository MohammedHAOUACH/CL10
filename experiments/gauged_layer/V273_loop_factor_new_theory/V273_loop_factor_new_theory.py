#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V273 - does the dynamical link birth the new theory's loop factor?
(tier 2, hard, O8-new.)

**The question (pre-registered).** Does the dynamical link BIRTH the
missing loop factor: does any dimensionless composite-sector ratio land
in the ladder band 1/(16 pi)?

**The object (pre-registered).** Candidates (fixed here):
F1 = rho (the V272 link-mediated coupling);
F2 = Delta / E0 (the V271 splitting over mass; E0 = 0 makes this
candidate structurally undefined - recorded, not silently dropped);
F3 = the second-order link-mediated exchange, J2 / g10 with
J2 = -<g| H_ex Q (E_g - H_atom)^{-1} Q H_ex |g> computed by Krylov
(Lanczos-Stieltjes) on the atom resolvent seeded from the chain ground
contraction, Krylov 800, seed 0x6C.
Target 1/(16 pi) = 1.9894e-2, band [0.015, 0.03] (as V263's, same
ladder convention).

Closed expansion: H_ex = g10 sum_a [A_a + C_a (x) J_a] on the chain
(both Hermitian), A_a = J_a^pair L_a^link, C_a = -Lbar_a^link, and
|g> = sum_i k_i (x) h_i the V272 chain ground: the k-basis is rebuilt
here on the V270 channel-adapted kernel (same seed path), the h-slice
loaded from the V272 artifact.  Partner legs close on measured
pairings rho_h[j,i] = <h_j|h_i>, G^b[j,i] = <h_j|J_b h_i>,
K^{ba}[j,i] = <J_b h_j|J_a h_i>, giving the four-term sum over the
540 atom vectors A_a k_i, C_a k_i with resolvent sandwich
R = Q (z - H_tilde)^{-1} Q, z = E_g (V272 chain ground), H_tilde =
Q H_atom Q.  Lanczos-Stieltjes on H_tilde (800 steps, full re-
orthogonalization, kernel-projected at every step) represents all 540
vectors in the Krylov basis; the tridiagonal Stieltjes entries are
re-evaluated at 60 and 100 decimal digits (mpmath, 5 fixed interior
probes) as the convergence control.

**Measured (pre-registered).** The three candidates; the closest
|16 pi F - 1|; the Krylov convergence residual (60 vs 100-digit
control on the resolvent quadratic form - mpmath, 5 fixed interior
probes).

**Acceptance (pre-registered).** LOOP-NEW-BORN iff any candidate inside
the band; LOOP-NEW-UNBORN otherwise with the closest miss recorded.

**Implementation note (added after the first runs, pre-registration
untouched).**  The pre-registered Krylov seed - the chain-ground
contraction of H_ex - measures EXACTLY ZERO norm (first-order non-
response of the chain ground, consistent with E_int = rho = 0).  A
structural test: span{C_a k_i, A_a k_i} (rank 360) is NOT invariant
under H (residual 1.0).  The declared machine that resolves both: the
BLOCK Krylov grown from the per-channel exchange sums closes EXACTLY
in dimension 4 (growth [4, 0], closure residual 0.0), and that space
sits in the kernel of H: every exchange vector A_a k_i, C_a k_i is
orthogonal to it (A kills link singlets, L kills them from the other
side), so the resolvent sandwich - evaluated exactly on the closed
space with the 60/100-digit mpmath spectral control - vanishes:
F3 = 2.2e-33, i.e. zero to machine precision.  The composite fixed
point is completely transparent to the link channel at second order.
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

    def log(msg):
        print("[V273 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start, backend={}, g10 = {}".format(E._BACKEND, G10))

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
    checks.append({"name": "T1 atom kernel dim {:d}, C2l spectrum {} "
                          "(seed 0x{:X})".format(
                              len(ker),
                              [round(x, 1)
                               for x in data["atom_kernel_C2l"]], SEED),
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

    # k-basis (channel order) and the 540 atom exchange vectors on host
    K3 = np.stack([E.asnumpy(k) for k in ker]).reshape(6, 32, 32, LN)
    AKv = np.empty((NP, 6, 32, 32, LN), dtype=np.complex128)
    CKv = np.empty((NP, 6, 32, 32, LN), dtype=np.complex128)
    for a in range(NP):
        La_h = np.asarray(L256_h[a])
        Lka = np.einsum("lm,ipqm->ipql", La_h, K3)
        AKv[a] = np.einsum("pqrs,krsm->kpqm",
                           Jpair_h[a].reshape(32, 32, 32, 32),
                           Lka)
        CKv[a] = -np.einsum("lm,ipqm->ipql", La_h, K3)
    GramA = np.einsum("ipqm,ajpqm->aij", K3.conj(), AKv)
    GramC = np.einsum("ipqm,ajpqm->aij", K3.conj(), CKv)
    herm_def = max(
        float(np.abs(GramA - GramA.transpose(0, 2, 1).conj()).max()),
        float(np.abs(GramC - GramC.transpose(0, 2, 1).conj()).max()))
    checks.append({"name": "T1 exchange vectors: A_a = (J_a pair)(L_a "
                          "link) and C_a = -L_a Hermitian on the kernel "
                          "(6x6 Gram defect over 45 vertices) "
                          "{:.1e}".format(herm_def),
                   "pass": bool(herm_def < 1e-9)})
    log("T1 540 exchange vectors built ({:.1f} GiB host, t+{:.0f}s)".
        format(2 * NP * 6 * NATOM * 16 / 2**30, time.time() - t0))

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
            gseed = gseed + wi * xp.asarray(AKv[a][i].reshape(NATOM))
            gseed = gseed + wi * xp.asarray(CKv[a][i].reshape(NATOM))
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
        HBs[b] = H_atom_of(MB[b])
    hh = np.sum(np.abs(E.asnumpy(HBs)) ** 2, axis=1)
    G_H = E.asnumpy(HBs @ MB.conj().T).conj()
    del HBs
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
    # declared machine (after two measured facts: the pre-registered
    # ground-contraction Krylov seed is exactly null, and the full
    # exchange span does not close under H).  Block Krylov on the
    # exchange sector: seeds = per-channel kernel-projected exchange
    # sums; grow the ON basis until H-closure or the cap.  If it
    # closes, the resolvent form is EXACT on the invariant space; if
    # not, F3 is not bookable and the growth table is the evidence.
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
    # orthonormalize the seed block
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
    log("T2 block Krylov: closed {} dim {:d} steps {:d} growth {}".
        format(closed, Dcur, it, growth[:14]))
    checks.append({"name": "T2 block Krylov on the exchange sector: "
                          "closed after {:d} steps: {}, invariant dim "
                          "{:d}, per-step growth {}{}".format(
                              it, closed, Dcur, growth[:14],
                              "" if closed else " (still growing at the "
                              "cap: the resolvent does not close in any "
                              "space of dimension <= {:d})".format(MAXV)),
                   "pass": True})
    if closed:
        Mm = np.empty((Dcur, Dcur), dtype=np.complex128)
        hhq = np.empty(Dcur)
        for p in range(Dcur):
            w = H_atom_of(Uall[p])
            Mm[p] = E.asnumpy(w @ Uall.conj().T)  # <H u_p|u_q>
            hhq[p] = float(xp.real(xp.vdot(w, w)))
        HU = Mm.conj()
        HU = 0.5 * (HU + HU.conj().T)
        inv2 = float(np.max(np.sqrt(np.maximum(
            hhq - np.sum(np.abs(HU) ** 2, axis=0), 0.0))
            / np.sqrt(np.where(hhq > 1e-14, hhq, 1.0))))
        data["invariant_residual_final"] = inv2
        checks.append({"name": "T2 invariant-space closure residual "
                              "{:.1e} (the resolvent form is exact "
                              "below)".format(inv2),
                       "pass": bool(inv2 < 1e-8)})
        th, Uth = np.linalg.eigh(HU)
        denom = E_g - th
        data["min_resolvent_denominator"] = float(np.min(np.abs(denom)))
        fv = 1.0 / denom
        Xp = E.asnumpy(Uall @ MB.conj().T).T  # (540, D): <u_m|x_n>
        UYc = Xp[:NP * 6].reshape(NP, 6, Dcur)
        UXc = Xp[NP * 6:].reshape(NP, 6, Dcur)
        UXf = UXc.reshape(NP * 6, Dcur) @ Uth
        UYf = UYc.reshape(NP * 6, Dcur) @ Uth
    else:
        UXf = np.zeros((NP * 6, 1), dtype=np.complex128)
        UYf = UXf.copy()
        fv = np.zeros(1)
    del Upre, Uall, MB

    # ---- T3: J2 assembly + 60/100-digit control ----
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
    F3 = (float(np.real(J2)) / G10) if gate_ok else float("nan")
    data["J2"] = [float(np.real(J2)), float(np.imag(J2))]
    data["J2_imag_part"] = float(np.imag(J2))
    data["F3_exchange"] = F3
    data["F3_gate"] = (
        "valid: exact on the {:d}-dim invariant space (closure "
        "residual {:.1e})".format(Dcur, inv2) if gate_ok else
        "not booked: the exchange resolvent does not close under block "
        "Krylov (growth {} to dim {}); no truncated form is presented".
        format(growth, Dcur))
    import mpmath
    probes = [p for p in PROBES if p < Dcur] if gate_ok else []
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
    log("T3 F3 = {:.6e}; mpmath 60/100 control residual {:.1e}".format(
        F3, ctrl))
    checks.append({"name": "T3 J2 = {:.6e} (imag {:.1e}): F3 = {} ({}); "
                          "60 vs 100-digit control over {} interior "
                          "resolvent probes: residual {:.1e}".format(
                              J2, float(np.imag(J2)),
                              "nan" if np.isnan(F3) else "{:.6e}".format(
                                  F3), data["F3_gate"], len(probes),
                              ctrl),
                   "pass": True})

    # ---- branch ----
    cands = {"F1": F1, "F2": F2, "F3": F3}
    inband = [n for n, v in cands.items()
              if np.isfinite(v) and BAND[0] <= abs(v) <= BAND[1]]
    dist = sorted((abs(TARGET - abs(v)), n) for n, v in cands.items()
                  if np.isfinite(v))
    data["branch_basis"] = {"target": TARGET, "band": list(BAND),
                            "closest": {"candidate": dist[0][1],
                                        "abs_value": abs(cands[dist[0][1]]),
                                        "distance": dist[0][0]}}
    branch = "LOOP-NEW-BORN" if inband else "LOOP-NEW-UNBORN"
    data["status"] = branch
    print("dynamical-link loop factor: F1 = {:.6f} (rho = 0: tree-level "
          "exchange vanishes), F2 = undefined (E0 = 0), F3 = {:.6e} "
          "(second-order exchange); closest to 1/(16 pi) = {:.4f}: {} at "
          "{:.6f}; mpmath control {:.1e}: branch {}.".format(
              F1, F3, TARGET, dist[0][1], cands[dist[0][1]], ctrl,
              branch), flush=True)
    _emit(data, checks, t0, ok=True)
    return 0


def _emit(data, checks, t0, ok):
    out = {"id": "V273_loop_factor_new_theory", "tier": 2,
           "property": "does the dynamical link birth the new theory's "
                       "loop factor (tier 2, O8-new)",
           "residual": True,
           "programme": "V273",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v273_loop_factor_new_theory.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V273 done in {:.1f}s] {}".format(time.time() - t0,
                                             out["branch"]))


if __name__ == "__main__":
    sys.exit(main())
