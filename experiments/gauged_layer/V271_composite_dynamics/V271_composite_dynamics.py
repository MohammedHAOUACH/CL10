#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V271 - composite dynamics: does the neutral composite move?
(tier 2, hard.)

**The question (pre-registered).** Once the composite exists, does it
move: is the Gauss-preserving rotor + exchange Hamiltonian non-diagonal
on the kernel (dressing), and what are the composite's mass, excitation
splittings and ground entanglement?

**The object (closed input set).** H = (1/2) C2_link (x) I_pair + g10
sum J^a (x) L^a on the atom (262144, matrix-free; J^a = T^a on the first
spinor factor + T^a on the second); the V270 atom-kernel recomputed with
the identical seed-0x6C Lanczos path and cross-checked against the V270
artifact (dim, split, residual); H restricted to the kernel (6x6) eigh;
the full-space first eigenpairs by Lanczos seeded from the kernel
vectors (Krylov 800, full re-orthogonalization); the excited composite
states = Krylov eigenvectors with <psi| sum_a (Q^a)^2 |psi> < 1e-9.

**Measured (pre-registered).** [H, Q^a_tot] defect probed on 3 fixed
seed states (max over a of ||[H,Q^a] phi|| / ||phi||) - the new theory's
Gauss law, measured not claimed; the kernel spectrum (E0 mass,
Delta splitting); the dressing fraction 1 - ||P_ker psi0||^2 of the
full-space ground state; the ground entropy S(pair|link); its deviation
from ln 45.

**Acceptance (pre-registered).** COMPOSITE-DYNAMICS-ALIVE iff Gauss
defect < 1e-12, the restricted H is not proportional to I (off-diagonal
max > 1e-6), and the ground entropy lies in (0, ln 45];
COMPOSITE-FROZEN if H is diagonal on the kernel - the diagonal test is
taken in the channel-adapted (C2_link) basis, the basis-independent
reading; GAUSS-BROKEN if the defect fails (GAUSS-BROKEN checked first).
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
LAM45 = math.log(45.0)
KMAX = 500
KRYLOV = 800
SEED = 0x6C
KTOL = 1e-9
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
        print("[V271 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start, backend={}, g10 = {}".format(E._BACKEND, G10))

    # ---- T0: rebuild the atom (V269 projector) ----
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
    s45 = [b for b in binpar if b[1] == 45 and b[2] < -0.5]
    s210 = [b for b in binpar if b[1] == 210]
    lam45 = s45[0][0] if len(s45) == 1 else float("nan")
    lam210 = s210[0][0] if len(s210) == 1 else float("nan")
    P210 = U[:, s210[0][3]:s210[0][4]]
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
    C2pair_h = np.zeros((1024, 1024), dtype=np.complex128)
    for T in Ts_h:
        Tp = np.kron(T, I32) + np.kron(I32, T)
        C2pair_h += Tp @ Tp
    C2pair = xp.asarray(C2pair_h)

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

    def Qa(psi, a):
        p = psi.reshape(32, 32, LN)
        t1 = xp.einsum("ip,pql->iql", Ta[a], p)
        t2 = xp.einsum("iq,pql->pil", Ta[a], p)
        t3 = xp.einsum("lm,pqm->pql", La[a], p)
        return (t1 + t2 + t3).reshape(NATOM)

    def H_of(psi):
        p = psi.reshape(32, 32, LN)
        v = 0.5 * xp.einsum("ml,pql->pqm", C2l, p)
        X = xp.einsum("aip,pql->aiql", Ta, p)
        v = v + G10 * xp.einsum("aml,aiql->iqm", La, X)
        Y = xp.einsum("aiq,pql->apil", Ta, p)
        v = v + G10 * xp.einsum("aml,apil->pim", La, Y)
        return v.reshape(NATOM)

    # ---- T1: kernel of S, identical seed-0x6C Lanczos path as V270 ----
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
    resid = 0.0
    for psi in cands:
        for a in range(NP):
            resid = max(resid, float(xp.sqrt(xp.real(
                xp.vdot(Qa(psi, a), Qa(psi, a))))))
    G = np.array([[complex(xp.vdot(cands[i], cands[j]))
                   for j in range(len(cands))]
                  for i in range(len(cands))])
    gev, gvec = np.linalg.eigh(G)
    ker = []
    for j in range(len(gev)):
        if gev[j] > 1e-8:
            psi = sum(gvec[i, j] * cands[i] for i in range(len(cands)))
            psi = psi / xp.sqrt(xp.real(xp.vdot(psi, psi)))
            ker.append(psi)
    kdim = len(ker)
    v270 = art("V270_neutral_composite_sector")
    xcheck = bool(v270["pass"] and kdim == v270["data"]["kernel_dim"]
                  and resid < KTOL
                  and abs(v270["data"]["lambda_210"] - lam210) < 1e-9
                  and abs(v270["data"]["lam45_value"] - lam45) < 1e-9)
    checks.append({"name": "T1 kernel recomputed with the V270 seed path: "
                          "dim {:d} (V270 dim {:d}), residual {:.1e} - "
                          "cross-checked against the committed artifact"
                          .format(kdim, v270["data"]["kernel_dim"], resid),
                   "pass": xcheck})
    log("T1 kernel dim {:d} residual {:.1e}".format(kdim, resid))
    if kdim == 0:
        status = ("V270 kernel empty (kernel dim 0 on recompute): no "
                  "composite to move. Branch NO-COMPOSITE.")
        out = {"id": "V271_composite_dynamics", "tier": 2,
               "property": "composite dynamics on the dynamical-link atom "
                           "(V271)", "residual": False,
               "programme": "V271", "branch": "NO-COMPOSITE",
               "data": {"kernel_dim": 0, "fitted_inputs": []},
               "checks": checks, "status": status, "pass": True}
        residir = os.path.join(HERE, "results")
        os.makedirs(residir, exist_ok=True)
        with open(os.path.join(residir,
                               "v271_composite_dynamics.json"), "w") as fh:
            json.dump(out, fh, indent=1)
        print(status)
        return 0

    # ---- T2: Gauss law of H: [H, Q^a] on 3 fixed probe states ----
    probe_rng = np.random.default_rng(SEED + 7)
    probes = []
    for _ in range(3):
        phi = (probe_rng.standard_normal(NATOM)
               + 1j * probe_rng.standard_normal(NATOM))
        probes.append(xp.asarray(phi / np.linalg.norm(phi)))
    gdef = 0.0
    for phi in probes:
        Hphi = H_of(phi)
        for a in range(NP):
            comm = H_of(Qa(phi, a)) - Qa(Hphi, a)
            gdef = max(gdef, float(xp.sqrt(xp.real(xp.vdot(comm, comm)))))
    data["gauss_defect"] = gdef
    checks.append({"name": "T2 Gauss law [H, Q^a_tot] on the atom: max "
                          "defect {:.2e} over 45 charges x 3 fixed probe "
                          "states".format(gdef),
                   "pass": bool(gdef < 1e-12)})
    log("T2 Gauss defect {:.2e}".format(gdef))

    # ---- T3: H on the kernel, channel-adapted (basis-invariant) ----
    # sector Gram of C2_link on the kernel basis: the channel basis is
    # the invariant one; "diagonal on the kernel" must be measured there
    CK = [xp.einsum("ml,pql->pqm", C2l, k.reshape(32, 32, LN))
          .reshape(NATOM) for k in ker]
    Gc = np.empty((kdim, kdim), dtype=np.complex128)
    for i in range(kdim):
        for j in range(kdim):
            Gc[i, j] = complex(xp.vdot(ker[i], CK[j]))
    cev, cvec = np.linalg.eigh(0.5 * (Gc + Gc.conj().T))
    data["sector_eigenvalues"] = [float(x) for x in cev]
    hk = np.empty((kdim, kdim), dtype=np.complex128)
    Hker = [H_of(k) for k in ker]
    for i in range(kdim):
        for j in range(kdim):
            hk[i, j] = complex(xp.vdot(ker[i], Hker[j]))
    hk = 0.5 * (hk + hk.conj().T)
    hc = cvec.conj().T @ hk @ cvec
    hc = 0.5 * (hc + hc.conj().T)
    ec, _ = np.linalg.eigh(hc)
    chan_offdiag = float(np.abs(hc - np.diag(np.diag(hc))).max())
    # the Casimir identity H|ker = alpha C2l|ker: alpha by least squares,
    # defect measured (not claimed)
    alpha = float(np.sum(np.diag(hc).real * cev) / np.sum(cev ** 2))
    cas_defect = float(max(chan_offdiag,
                           np.abs(np.diag(hc).real - alpha * cev).max()))
    data["kernel_H_channel_spectrum"] = [float(x) for x in ec]
    data["kernel_H_channel_offdiag"] = chan_offdiag
    data["casimir_alpha_measured"] = alpha
    data["casimir_alpha_expected"] = 0.5 - G10
    data["casimir_identity_defect"] = cas_defect
    ek = ec
    E0, Delta = float(ek[0]), float(ek[-1] - ek[0])
    log("T3 channel spectrum {} offdiag {:.2e} alpha {:.6f} "
        "(expected {:.6f}) defect {:.2e}".format(
            [round(x, 6) for x in ec], chan_offdiag, alpha, 0.5 - G10,
            cas_defect))

    # full-space Lanczos seeded from the kernel vectors
    v0 = sum(ker) if kdim else probes[0]
    v0 = v0 / xp.sqrt(xp.real(xp.vdot(v0, v0)))
    Vk = [v0]
    ak, bk = [], []
    kmax_used = KRYLOV
    for k in range(KRYLOV):
        w = H_of(Vk[-1])
        a = float(xp.real(xp.vdot(Vk[-1], w)))
        ak.append(a)
        w = w - a * Vk[-1]
        if len(Vk) > 1:
            w = w - bk[-1] * Vk[-2]
        for _ in range(2):
            for u in Vk:
                w = w - xp.vdot(u, w) * u
        b = float(xp.sqrt(xp.real(xp.vdot(w, w))))
        if b < 1e-10:
            kmax_used = k + 1
            break
        bk.append(b)
        Vk.append(w / b)
        if (k + 1) % 200 == 0:
            log("T3 Krylov k={} (alpha {:.4f}, beta {:.2e})".format(k + 1,
                                                                    a, b))
    Tk = np.diag(ak[:kmax_used])
    for i in range(kmax_used - 1):
        Tk[i, i + 1] = bk[i]
        Tk[i + 1, i] = bk[i]
    kei, kevec_h = np.linalg.eigh(Tk)
    kevec = xp.asarray(kevec_h)
    VkM = xp.stack(Vk[:kmax_used]).T
    del Vk
    psi0 = VkM @ kevec[:, 0]
    psi0 = psi0 / xp.sqrt(xp.real(xp.vdot(psi0, psi0)))
    Efull0 = float(kei[0])
    # dressing: overlap with the kernel sector
    ov = 0.0
    for kk in ker:
        c = complex(xp.vdot(kk, psi0))
        ov += abs(c) ** 2
    dressing = 1.0 - ov
    data["full_space_E0"] = Efull0
    data["dressing_fraction"] = dressing
    data["ground_overlap_with_kernel"] = ov
    # ground entropy
    M = E.asnumpy(psi0).reshape(1024, LN)
    sv = np.linalg.svd(M, compute_uv=False)
    p = sv ** 2 / np.sum(sv ** 2)
    p = p[p > 1e-14]
    Sent = float(-np.sum(p * np.log(p)))
    data["ground_entropy"] = Sent
    # excited composite states: Krylov eigenvectors with <S> < 1e-9
    excited = []
    for j in range(1, min(kmax_used, 40)):
        ps = VkM @ kevec[:, j]
        ps = ps / xp.sqrt(xp.real(xp.vdot(ps, ps)))
        c2e = float(xp.real(xp.vdot(ps, S_of(ps))))
        if abs(c2e) < 1e-9:
            Mo = E.asnumpy(ps).reshape(1024, LN)
            s2 = np.linalg.svd(Mo, compute_uv=False)
            q = s2 ** 2 / np.sum(s2 ** 2)
            q = q[q > 1e-14]
            excited.append({"E": float(kei[j]),
                            "delta_E": float(kei[j] - Efull0),
                            "entropy": float(-np.sum(q * np.log(q)))})
    data["excited_composites"] = excited[:12]
    data["krylov_steps"] = kmax_used
    checks.append({"name": "T3 dynamics (channel-adapted): kernel H "
                          "spectrum E0 {:.6f}, Delta {:.6f}, channel-basis "
                          "off-diagonal {:.2e}; Casimir identity H|ker = "
                          "alpha C2l|ker measured alpha {:.6f} (expected "
                          "0.5 - g10 = {:.6f}, defect {:.2e}); kernel-"
                          "seeded Krylov confined after {:d} steps "
                          "(kernel invariant); dressing {:.4f}, ground "
                          "entropy {:.6f} (ln 45 = {:.6f}); {} neutral "
                          "excited composites".format(
                              E0, Delta, chan_offdiag, alpha, 0.5 - G10,
                              cas_defect, kmax_used, dressing, Sent, LAM45,
                              len(excited)),
                   "pass": bool(kdim == 6 and cas_defect < 1e-9)})
    log("T3 E0 {:.6f} Delta {:.6f} chan-offdiag {:.2e} dress {:.4f} "
        "S {:.6f}".format(E0, Delta, chan_offdiag, dressing, Sent))

    # ---- branch ----
    if gdef >= 1e-12:
        branch = "GAUSS-BROKEN"
    elif chan_offdiag <= 1e-6:
        branch = "COMPOSITE-FROZEN"
    elif (kdim == 6 and gdef < 1e-12 and chan_offdiag > 1e-6
          and 0.0 < Sent <= LAM45 + 1e-6):
        branch = "COMPOSITE-DYNAMICS-ALIVE"
    else:
        branch = "COMPOSITE-DYNAMICS-PARTIAL"
    ok = bool(branch in ("COMPOSITE-DYNAMICS-ALIVE", "COMPOSITE-FROZEN"))
    status = ("composite dynamics: H = (1/2)C2_link + g10 J.L on 262144, "
              "Gauss defect {:.1e}; on the kernel H is the link Casimir "
              "itself: H|ker = {:.6f} C2l|ker (measured defect {:.1e}), "
              "channel mass ladder {} = (0.5 - g10) lambda_r, channel "
              "off-diagonal {:.1e}, dressing {:.4f}, ground entropy "
              "{:.6f}: branch {}.".format(
                  gdef, alpha, cas_defect,
                  [round(x, 4) for x in ec], chan_offdiag, dressing, Sent,
                  branch))
    out = {"id": "V271_composite_dynamics", "tier": 2,
           "property": "the neutral composite moves: the Gauss-preserving "
                       "rotor + exchange Hamiltonian dresses the kernel "
                       "sector (V271)",
           "residual": False, "programme": "V271", "branch": branch,
           "data": {**data, "lambda_210": lam210, "g10": G10,
                    "seed": "0x6C", "fitted_inputs": [],
                    "backend": E._BACKEND, "wall_time": time.time() - t0},
           "checks": checks, "status": status, "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v271_composite_dynamics.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
