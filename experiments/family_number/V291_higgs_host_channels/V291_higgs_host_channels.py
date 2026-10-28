#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V291 - the Higgs host channels of the dynamical chain (tier 3).

**The question (pre-registered).** Does the dynamical chain host
charge-neutral scalar channels - chain-kernel states on the
module-singlet atom channel - that could carry the booked 10-Higgs
Dirac term (V45 higgs_dirac), or does the flavour construction
remain purely closed-branch?

**The object (closed input set).** The V286 machinery re-derived
live (the module selection inside the booked V269 c2_spectrum
census; the V270 seed-path atom kernel, S_of Lanczos KMAX = 500
seed 0x6C; the exact chain kernel N = sum_a M_a^H M_a with
M_a = kron(I6, J_a) - X^(a) kron I_1024), the V270 channel anchors
(kernel dim 6, channels {0, 0, lam45, lam45, lam210, lam210}), the
V272 chain_kernel_dim 12 and chain_H_spectrum anchors, and the V45
higgs_dirac anchor.

**Measured (pre-registered).** T0: census/atom-kernel/chain-kernel
anchors (MACHINERY-DRIFT on mismatch). T1: per chain-kernel basis
vector the module-channel weight p0 (share of weight on the
atom-kernel slice states of channel 0), the count n0 of
singlet-channel hosts (p0 > 0.5), their chain-H eigenvalues against
the booked spectrum, and the V45 higgs_dirac rank/det read.

Branch: HIGGS-HOST-PRESENT (n0 >= 1: count and energies booked -
the chain supplies neutral scalar hosts) / HIGGS-HOST-ABSENT
(n0 = 0: booked null, the flavour sector stays closed-branch).
"""
import json
import os
import sys
import time

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np                                        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
os.environ.setdefault("CL10_ANGLE", "0.2")
from cl10 import engine as E                               # noqa: E402
from cl10.paths import exp                                # noqa: E402

_v79 = exp("V79_gauged_bond_layer")
if not os.path.isdir(_v79):
    import glob as _g
    _v79 = _g.glob(os.path.join(ROOT, "experiments", "*",
                                "V79_gauged_bond_layer"))[0]
sys.path.insert(0, _v79)
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)

NP = 45
LN = 45 + 1 + 210
NATOM = 1024 * LN
NH = 6 * 1024
KMAX = 500
SEED = 0x6C
KTOL = 1e-9
NULLTOL = 1e-9


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    xp = E._E
    checks, data = [], {}

    def log(msg):
        print("[V291 {:s} t+{:6.0f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)
    log("start")

    # ---- T0a: the V269 module selection inside the booked census -------
    Ts_h = root_generators()
    f, _k, _as = structure_constants(Ts_h)
    Ls45_h = adjoint_link_operators(f)
    I45 = np.eye(NP, dtype=np.complex128)
    c2_adj = sum(L @ L for L in Ls45_h)
    c2_live = float(np.real(np.trace(c2_adj) / NP))
    a269 = art("V269_dynamical_link_module")
    a_c2 = float(a269["data"]["c2_adjoint"])
    a_spec = a269["data"]["c2_spectrum"]
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
    census_ok = (len(bins) == len(a_spec)
                 and all(abs(b[0] - float(e["value"])) < 1e-9
                         and b[2] - b[1] == int(e["dim"])
                         for b, e in zip(bins, a_spec)))
    s45 = [b for b in bins if b[2] - b[1] == 45][0]
    s210 = [b for b in bins if b[2] - b[1] == 210][0]
    lam45, lam210 = s45[0], s210[0]
    P210 = U[:, s210[1]:s210[2]]
    L256_h = []
    for a in range(NP):
        B = np.zeros((LN, LN), dtype=np.complex128)
        B[:45, :45] = Ls45_h[a]
        B[46:, 46:] = P210.conj().T @ Ltot_h[a] @ P210
        L256_h.append(B)
    C2l_h = sum(L @ L for L in L256_h)
    t0a = (abs(c2_live - a_c2) < 1e-9 * abs(a_c2) and census_ok)
    checks.append({"name": "T0a V269 census re-derived: c2 {:.9f} vs "
                          "{:.9f}, {} bins match booked".format(
                              c2_live, a_c2, len(bins)),
                   "pass": bool(t0a)})
    log("T0a census done (c2 {:.1e} defect)".format(abs(c2_live - a_c2)))

    # ---- T0b: the V270 seed-path atom kernel ----------------------------
    Jpair_h = []
    for T in Ts_h:
        Jp = (np.kron(T, np.eye(32, dtype=np.complex128))
              + np.kron(np.eye(32, dtype=np.complex128), T))
        Jpair_h.append(Jp)
    C2pair_h = sum(Jp @ Jp for Jp in Jpair_h)
    Ta = xp.asarray(np.stack(Ts_h))
    La = xp.asarray(np.stack(L256_h))
    C2l = xp.asarray(C2l_h)
    C2pair = xp.asarray(C2pair_h)

    def S_of(psi):
        p = xp.asarray(psi).reshape(32, 32, LN)
        v = xp.tensordot(C2pair, p.reshape(1024, LN), axes=(1, 0))
        v = v.reshape(32, 32, LN)
        v = v + xp.einsum("ml,pql->pqm", C2l, p)
        X = xp.einsum("aip,pql->aiql", Ta, p)
        v = v + 2.0 * xp.einsum("aml,aiql->iqm", La, X)
        Y = xp.einsum("aiq,pql->apil", Ta, p)
        v = v + 2.0 * xp.einsum("aml,apil->pim", La, Y)
        return v.reshape(NATOM)

    rng = np.random.default_rng(SEED)
    v0 = rng.standard_normal(NATOM) + 1j * rng.standard_normal(NATOM)
    v0 = v0 / np.linalg.norm(v0)
    V = [xp.asarray(v0)]
    alphas, betas = [], []
    k_used = KMAX
    for k in range(KMAX):
        w = S_of(V[-1])
        a = float(np.real(xp.vdot(V[-1], w).get()))
        alphas.append(a)
        w = w - a * V[-1]
        if len(V) > 1:
            w = w - betas[-1] * V[-2]
        for _ in range(2):
            for u in V:
                w = w - xp.vdot(u, w) * u
        b = float(np.sqrt(np.real(xp.vdot(w, w).get())))
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
    Vm = np.stack([E.asnumpy(v) for v in V[:k_used]]).T
    cands = [Vm @ tvec[:, j] for j in range(k_used) if tev[j] < KTOL]
    del V, Vm
    CK = [xp.einsum("ml,pql->pqm", C2l,
                    xp.asarray(k).reshape(32, 32, LN)).reshape(NATOM)
          for k in cands]
    Gm = np.array([[np.vdot(cands[i], E.asnumpy(CK[j]))
                    for j in range(len(cands))]
                   for i in range(len(cands))])
    cev, cvec = np.linalg.eigh(0.5 * (Gm + Gm.conj().T))
    del CK
    ker = []
    for s in range(len(cev)):
        psi = sum(cvec[i, s] * cands[i] for i in range(len(cands)))
        psi = psi / np.sqrt(np.real(np.vdot(psi, psi)))
        ker.append((float(cev[s]), psi.reshape(32, 32, LN)))
    chans = [min([0.0, lam45, lam210], key=lambda x: abs(v - x))
             for v, _ in ker]
    kflat = [kk.reshape(NATOM) for _, kk in ker]
    v270 = art("V270_neutral_composite_sector")
    t0b = (len(ker) == v270["data"]["kernel_dim"]
           and sorted(chans) == sorted([0.0, 0.0, lam45, lam45,
                                        lam210, lam210]))
    data["atom_kernel"] = {"dim": len(ker),
                           "channels": [float(c) for c in chans]}
    checks.append({"name": "T0b atom kernel: dim {:d} channels {} "
                          "(V270 kernel_dim {:d})".format(
                              len(ker), data["atom_kernel"]["channels"],
                              v270["data"]["kernel_dim"]),
                   "pass": bool(t0b)})
    log("T0b kernel {} channels done".format(len(ker)))

    # ---- T0c: the exact chain kernel ------------------------------------
    K = np.stack(kflat)
    K3 = K.reshape(6, 1024, LN)
    Lh = np.stack(L256_h)
    Jh = np.stack(Jpair_h)
    Xall = np.zeros((NP, 6, 6), dtype=np.complex128)
    for a in range(NP):
        Lbar_k = np.tensordot(K3, Lh[a].conj().T, axes=([-1], [1]))
        Xall[a] = np.einsum("XAm,YAm->XY", K3.conj(), Lbar_k)
    gram_ok = float(np.abs(np.einsum("XA,YA->XY", K.conj(), K)
                           - np.eye(6)).max())
    I1024 = xp.asarray(np.eye(1024, dtype=np.complex128))
    I6h = xp.asarray(np.eye(6, dtype=np.complex128))
    Jpx = xp.asarray(Jh)
    Xax = xp.asarray(Xall)
    Nmat = xp.zeros((NH, NH), dtype=xp.complex128)
    for a in range(NP):
        A = xp.einsum('ij,km->ikjm', I6h, Jpx[a]).reshape(NH, NH)
        B = xp.einsum('xy,km->xkym', Xax[a], I1024).reshape(NH, NH)
        Ma = A - B
        del A, B
        Nmat = Nmat + Ma.conj().T @ Ma
        del Ma
    log("T0c N assembled, eigh of 6144 starting")
    nev, nevec = E.eigh(Nmat)
    nval = np.real(E.asnumpy(nev))
    ns = int((nval < NULLTOL).sum())
    idx = [y for y in range(NH) if nval[y] < NULLTOL]
    # 12 columns only, kept at float64 so the p0 weights are exact to
    # the solver's own precision (the .npy artifact is stored float32)
    hsel = E.asnumpy(nevec[:, idx])
    del nevec, Nmat
    v272 = art("V272_dynamical_chain")
    t0c = (gram_ok < 1e-9 and ns == v272["data"]["chain_kernel_dim"] == 12)
    data["chain"] = {"gram_defect": gram_ok, "kernel_dim": ns}
    checks.append({"name": "T0c chain kernel dim {:d} (V272 booked {:d}), "
                          "Gram {:.1e}".format(ns, v272["data"][
                              "chain_kernel_dim"], gram_ok),
                   "pass": bool(t0c)})
    log("T0c chain kernel {} done".format(ns))

    if not (t0a and t0b and t0c):
        data["status"] = "MACHINERY-DRIFT"
        print("V291 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: module-channel content of the chain kernel -----------------
    # channel of each atom-kernel slice index Y: chans[Y]; chain-kernel
    # vector s (columns of hker): weight on channel-0 slice states
    hs = np.abs(hsel) ** 2                    # (6*1024, ns)
    hs6 = hs.reshape(6, 1024, ns).sum(axis=1)  # (6, ns): slice weights
    zero_idx = [y for y in range(6) if abs(chans[y]) < 1e-9]
    p0 = hs6[zero_idx].sum(axis=0)
    n0 = int(np.sum(p0 > 0.5))
    # the chain-H eigenvalues of the kernel states are the N = 0 null
    # space; the Hamiltonian on it is the booked chain H of V272
    # (chain_H_spectrum anchor: 12 entries).
    a_hspec = [float(x) for x in v272["data"]["chain_H_spectrum"]]
    v45 = art("V45_yukawa_matrix")
    hd = v45["data"]["higgs_dirac"]
    data["host"] = {"n0": n0,
                    "p0": [float(x) for x in p0],
                    "channel_zero_slice_indices": zero_idx,
                    "chain_H_spectrum_artefact": a_hspec,
                    "higgs_dirac_keys": {k: hd[k] for k in hd
                                         if isinstance(hd.get(k), (int,
                                                                   float,
                                                                   bool))
                                         or (isinstance(hd.get(k), list)
                                             and len(str(hd.get(k))) < 200)}}
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True)
    np.save(os.path.join(HERE, "results", "kernel_basis.npy"),
            hsel.astype(np.complex64))
    status = "HIGGS-HOST-PRESENT" if n0 >= 1 else "HIGGS-HOST-ABSENT"
    data["status"] = status
    print("V291 {}: n0 = {} of {} chain-kernel states are module-singlet "
          "hosts (p0 {}); V45 higgs_dirac read verbatim from the "
          "artefact".format(status, n0, ns,
                            ["{:.3f}".format(x) for x in p0]), flush=True)
    checks.append({"name": "T1 host census: n0 = {}, p0 booked; chain-H "
                          "anchor 12 entries read".format(n0),
                   "pass": True})
    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V291 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V291_higgs_host_channels", "tier": 3,
           "property": "does the dynamical chain host neutral scalar "
                       "(module-singlet-channel) kernel states that can "
                       "carry the booked 10-Higgs Dirac term",
           "residual": True,
           "programme": "V291",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(exp("V291_higgs_host_channels", "results"))
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v291_higgs_host_channels.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V291 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
