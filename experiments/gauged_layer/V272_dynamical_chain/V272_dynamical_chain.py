#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V272 - two dynamical-link composites sharing one link (tier 2, hard).

**The question (pre-registered).** Two dynamical-link composites sharing
one link: is the chain neutral (predicted dim 12), and what effective
coupling does the shared link generate between the two atoms (the
candidate scale of the new theory)?

**The object (pre-registered).** Chain P1 (x) l (x) P2 with Gauss
T1 + L = 0, T2 - L = 0 (vertex 2 carries the shared link through the
conjugate module: numerically T2 - Lbar = 0, the representation-closed
reading of "T2 - L": [J_a - Lbar_a, J_b - Lbar_b] = i f (J - Lbar)
closes, while T2 - L alone does not); chain kernel constructed by
singlet contraction of two V270 atom-kernel bases over the shared link
(Gram to 1e-9; predicted dim sum_r m_r^2 = 12); chain H = (1/2) C2_l +
g10 (J1.L + J2.(-L)) in the contracted basis; interaction energy
E_int = E_chain(ground) - 2 E_atom(ground); effective coupling
rho = E_int / g10.

Construction (exact, closed, written before the run). Q1_a acts only
on the (pair1, link) slots, so for every fixed B-index the slice
Phi(.,.,B) lies in ker Q1 on P (x) L - which is exactly the 6-dim V270
atom-kernel space span{k_i}.  Hence the chain kernel is COMPLETELY
parametrized by
    Phi(A,la,B) = sum_{i=1..6} k_i(A,la) h_i(B),  h in C^{6 x 1024},
and vertex 2 becomes 45 linear conditions on h:
    M_a h = 0,  M_a = kron(I6, J_a) - X^{(a)} (x) I_{1024},
    X^{(a)}_{ji} = sum_{A,la,mu} conj(k_i(A,la)) conj(L_a[la,mu])
                   k_j(A,mu)   (the Lbar_a action transported onto the
                   k-slot),
using (J_a k)(A,nu) = -sum_mu L_a[nu,mu]k(A,mu) from the atom Gauss
law.  The chain kernel is computed exactly: N = sum_a M_a^H M_a
(6144 x 6144 dense, built from 45 block terms on the GPU), eigvalsh +
eigh, null space below 1e-9 = chain kernel.  Predicted dim 12.  Then
H matrix elements factorize through the measured slot-overlap matrices
C_{ij} = <k_i|C2l|k_j>, U^{(a)}_{ij} = <k_i|J_a L_a|k_j> (both slots),
Q^{(a)}_{ij} = <k_i|Lbar_a|k_j> - the 268M-vectors never enter the
eigenproblem.  ONE ground state IS materialized (1024 x 256 x 1024,
4.3 GiB) to MEASURE the chain Gauss defect on both vertices (90
charges, dense application).

**Measured (pre-registered).** Chain kernel dim (predicted 12, Gram
defect); the atom-vs-chain ground energies; E_int; rho; the chain Gauss
defect (both vertices, < 1e-12).

**Acceptance (pre-registered).** CHAIN-ALIVE iff dim = 12 and Gauss
defect < 1e-12 (else CHAIN-EMPTY if dim = 0, else GAUSS-BROKEN); rho
recorded either way.  0 < dim != 12 records CHAIN-PARTIAL (acceptance
unmet; every number measured).
"""
import json
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
NH = 6 * 1024                # h-space
KMAX = 500
SEED = 0x6C
KTOL = 1e-9
NULLTOL = 1e-9
G10 = coupling.G10


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}
    xp = E._E

    def log(msg):
        print("[V272 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start, backend={}, g10 = {}".format(E._BACKEND, G10))

    # ---- T0: link module + V270 seed-path kernel, channel-adapted ----
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
    lam210 = s210[0]
    lam45 = s45[0]
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
        Jp = np.kron(T, I32) + np.kron(I32, T)
        Jpair_h.append(Jp)
        C2pair_h += Jp @ Jp
    C2pair = xp.asarray(C2pair_h)
    Jp = xp.asarray(np.stack(Jpair_h))          # (45,1024,1024)

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
    ker = []
    for s in range(len(cev)):
        psi = sum(cvec[i, s] * cands[i] for i in range(len(cands)))
        psi = psi / xp.sqrt(xp.real(xp.vdot(psi, psi)))
        ker.append((float(cev[s]), E.asnumpy(psi).reshape(1024, LN)))
    chans = [min([0.0, lam45, lam210], key=lambda x: abs(v - x))
             for v, _ in ker]
    kflat = [k_.reshape(NATOM) for _, k_ in ker]   # 6 ON vectors
    v270 = art("V270_neutral_composite_sector")
    v271 = art("V271_composite_dynamics")
    t0pass = bool(v270["pass"] and v271["pass"]
                  and len(ker) == v270["data"]["kernel_dim"]
                  and sorted(chans) == sorted([0.0, 0.0, lam45, lam45,
                                               lam210, lam210]))
    checks.append({"name": "T0 V270/V271 artifacts load; V270 seed-path "
                          "kernel recomputed: dim {:d}, channels {} "
                          "(matches artifacts)".format(
                              len(ker),
                              {str(c): sum(1 for x in chans if x == c)
                               for c in sorted(set(chans))}),
                   "pass": t0pass})
    log("T0 kernel {} channels {}".format(
        len(ker), {c: sum(1 for x in chans if x == c)
                   for c in sorted(set(chans))}))

    # ---- T1: exact chain kernel: N = sum_a M_a^H M_a on the 6144 h-space
    K = np.stack(kflat)                         # (6, NATOM) host
    Nmat = xp.zeros((NH, NH), dtype=xp.complex128)
    I6 = np.eye(6, dtype=np.complex128)
    K3 = K.reshape(6, 1024, LN)
    # slot-overlap matrices (all measured from the recomputed kernels)
    CK3 = np.tensordot(K3, np.asarray(C2l_h).T, axes=([-1], [1]))
    Cc = np.einsum("XAm,YAm->XY", K3.conj(), CK3)
    Lh = np.stack(L256_h)                       # (45,256,256)
    Jh = np.stack(Jpair_h)                      # (45,1024,1024)
    Xall = np.zeros((NP, 6, 6), dtype=np.complex128)
    Uall = np.zeros((NP, 6, 6), dtype=np.complex128)
    gall_def = 0.0
    for a in range(NP):
        Lk = np.tensordot(K3, Lh[a].T, axes=([-1], [1]))
        Lbar_k = np.tensordot(K3, Lh[a].conj().T, axes=([-1], [1]))
        # X_a[X,Y] = <k_X, Lbar_a k_Y> (vertex-2 transport matrix)
        Xall[a] = np.einsum("XAm,YAm->XY", K3.conj(), Lbar_k)
        # U_a[X,Y] = <k_X, J_a L_a k_Y> (chain H J1.L form factor)
        Uall[a] = np.einsum("XAm,YPm->XY", K3.conj(),
                            np.einsum("PQ,YQm->YPm", Jh[a], Lk))
        # measured Gauss transport identity: <k_X| J_a k_Y> = -<k_X|L_a k_Y>
        Jk = np.einsum("PQ,YQm->YPm", Jh[a], K3)
        gJ = np.einsum("XAm,YPm->XY", K3.conj(), Jk)
        gL = np.einsum("XAm,YAm->XY", K3.conj(), Lk)
        gall_def = max(gall_def, float(np.abs(gJ + gL).max()))
    gram_ok = float(np.abs(np.einsum("XA,YA->XY", K.conj(), K)
                           - np.eye(6)).max())
    checks.append({"name": "T1 atom kernels: Gram defect {:.1e}; measured "
                          "Gauss transport J_a k = -L_a k to {:.1e} (this "
                          "identity is what makes the chain M_a exact)"
                          .format(gram_ok, gall_def),
                   "pass": bool(gram_ok < 1e-9 and gall_def < 1e-9)})
    # M_a = kron(I6, J_a) - X_a (x) I_1024 ; N = sum_a M_a^H M_a (exact)
    I1024 = xp.asarray(np.eye(1024, dtype=np.complex128))
    for a in range(NP):
        Ma = xp.kron(xp.asarray(I6), Jp[a]) \
            - xp.kron(xp.asarray(Xall[a]), I1024)
        Nmat = Nmat + Ma.conj().T @ Ma
        del Ma
    log("T1 N assembled (6144^2), starting eigh")
    nev, nevec = xp.linalg.eigh(Nmat)
    nev_h = E.asnumpy(nev)
    null = [j for j in range(NH) if nev_h[j] < NULLTOL]
    cdim = len(null)
    gap_above = float(nev_h[cdim]) if cdim < NH else float("nan")
    log("T1 chain kernel dim {:d} (predicted 12); N-eigenvalues below "
        "1e-9: {}; next {:.3e}".format(cdim, cdim, gap_above))
    data["chain_kernel_dim"] = cdim
    data["chain_gap_above"] = gap_above
    data["chain_N_min_eigs"] = [float(x) for x in nev_h[:min(cdim + 3, NH)]]
    HV = nevec[:, null].T.reshape(cdim, 6, 1024).copy() if cdim \
        else xp.zeros((0, 6, 1024), dtype=xp.complex128)
    gram_def = 0.0
    if cdim:
        HVflat = HV.reshape(cdim, 6 * 1024)
        Ghh = E.asnumpy(HVflat @ HVflat.conj().T)
        gram_def = float(np.abs(Ghh - np.eye(cdim)).max())

    checks.append({"name": "T1 exact chain kernel of M_a (45 vertices x "
                          "45 charges on the atom-kernel slice): dim {:d} "
                          "(predicted 12), Gram defect {:.1e}, first "
                          "non-null N eigenvalue {:.3e}".format(
                              cdim, gram_def, gap_above),
                   "pass": bool(cdim == 12 and gram_def < 1e-9)})
    del Nmat, nevec
    xp.get_default_memory_pool().free_all_blocks()

    # ---- T2: chain H in the kernel basis (factorized), ground state,
    # 90-charge Gauss measurement on the materialized ground ----
    Hc = np.zeros((cdim, cdim), dtype=np.complex128)
    E_atom_ground = float(v271["data"]["full_space_E0"])
    if cdim:
        Hsv = E.asnumpy(HV)                       # (cdim,6,1024)
        Hloc = 0.5 * Cc + G10 * np.sum(Uall, axis=0)
        Hc = np.einsum("siB,tjB,ij->st", Hsv.conj(), Hsv, Hloc,
                       optimize=True)
        for a in range(NP):
            tmp = np.einsum("BC,tjC->tjB", Jh[a], Hsv, optimize=True)
            Hc -= G10 * np.einsum("siB,tjB,ij->st", Hsv.conj(), tmp,
                                  Xall[a], optimize=True)
        Hc = 0.5 * (Hc + Hc.conj().T)
        Ec, Ew = np.linalg.eigh(Hc)
        data["chain_H_spectrum"] = [float(x) for x in Ec]
        log("T2 chain H spectrum (first 8) {}".format(
            [round(float(x), 6) for x in Ec[:8]]))
        E_chain_ground = float(Ec[0])
        # materialize the ground state and measure both vertices
        hgv = xp.asarray(np.einsum("s,sij->ij", Ew[:, 0], Hsv))
        _hg = E.asnumpy(hgv)
        data["chain_ground_h_re"] = _hg.real.reshape(-1).tolist()
        data["chain_ground_h_im"] = _hg.imag.reshape(-1).tolist()
        K3d = xp.asarray(K3)
        PHI = xp.einsum("iAl,iB->AlB", K3d, hgv)     # (1024,256,1024)
        nrm = float(xp.real(xp.vdot(PHI, PHI)))
        gauss_def = 0.0
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
            gauss_def = max(gauss_def, d1, d2)
            if (a + 1) % 15 == 0:
                log("T2 Gauss probe a={:d} defect so far {:.2e}".format(
                    a + 1, gauss_def))
        data["chain_gauss_defect"] = gauss_def
        data["chain_ground_norm"] = nrm
        checks.append({"name": "T2 chain Gauss defect on the materialized "
                              "ground state: max over 90 charges (both "
                              "vertices) {:.2e}".format(gauss_def),
                       "pass": bool(gauss_def < 1e-12)})
    else:
        E_chain_ground = float("nan")
        gauss_def = float("nan")
    E_int = E_chain_ground - 2.0 * E_atom_ground
    rho = E_int / G10
    data["E_atom_ground"] = E_atom_ground
    data["E_chain_ground"] = E_chain_ground
    data["E_int"] = E_int
    data["rho"] = rho
    checks.append({"name": "T2 energies: atom ground {:.6f}, chain ground "
                          "{:.6f}, E_int {:.6f}, rho = E_int/g10 = {:.6f}"
                          .format(E_atom_ground, E_chain_ground, E_int,
                                  rho),
                   "pass": bool(np.isfinite(E_int) and np.isfinite(rho))})
    log("T2 E_atom {:.6f} E_chain {:.6f} E_int {:.6f} rho {:.6f}".format(
        E_atom_ground, E_chain_ground, E_int, rho))

    # ---- branch ----
    if cdim == 0:
        branch = "CHAIN-EMPTY"
    elif cdim == 12 and gauss_def < 1e-12:
        branch = "CHAIN-ALIVE"
    elif gauss_def >= 1e-12:
        branch = "GAUSS-BROKEN"
    else:
        branch = "CHAIN-PARTIAL"
    ok = bool(branch in ("CHAIN-ALIVE", "CHAIN-EMPTY"))
    status = ("dynamical-link chain P1 (x) l (x) P2: exact kernel of the "
              "90 chain charges on the atom-kernel slice: dim {:d} "
              "(predicted 12, gap {:.3e}); chain ground {: .6f} vs atom "
              "{: .6f}: E_int {: .6f}, rho = {: .6f}; ground Gauss defect "
              "over 90 charges {:.1e}: branch {}.".format(
                  cdim, gap_above, E_chain_ground, E_atom_ground, E_int,
                  rho, gauss_def, branch))
    out = {"id": "V272_dynamical_chain", "tier": 2,
           "property": "two dynamical-link composites sharing one link "
                       "form a neutral chain and the shared link "
                       "generates an effective coupling (V272)",
           "residual": False, "programme": "V272", "branch": branch,
           "data": {**data, "lambda_210": lam210, "g10": G10,
                    "seed": "0x6C", "fitted_inputs": [],
                    "backend": E._BACKEND, "wall_time": time.time() - t0},
           "checks": checks, "status": status, "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v272_dynamical_chain.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
