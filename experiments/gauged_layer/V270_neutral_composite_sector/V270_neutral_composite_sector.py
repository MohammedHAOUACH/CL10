#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V270 - the neutral composite sector of the dynamical-link atom
(tier 1, hard, residual: false - the refonte's first verdict).

**The question (pre-registered).** Does the dynamical-link atom carry
gauge-neutral composites - is ker C2_tot non-empty - and is its sector
structure the predicted [2, 2, 2] under the link Casimir?

**The object (closed input set).** The V269 module rebuilt here from the
committed V79 generators by the same Casimir projector with the same
swap-parity identification, cross-checked against the V269 artifact to
1e-9: atom space 1024 (pair 32 x 32) x 256 (link 1(+)45(+)210); charges
Q^a psi(p,q,l) = (T^a p)(q,l) + (T^a q)(p,l) + (L^a l)(p,q);
S = sum_a (Q^a)^2 applied MATRIX-FREE on (32, 32, 256):
S psi = C2pair psi + C2link psi + 2 sum_a (T1_a L_a + T2_a L_a) psi;
C2pair = sum_a (T1_a + T2_a)^2 already carries the 2 T1 T2 cross term.
Lanczos on S (fixed seed 0x6C, <= 500 steps, full
re-orthogonalization, early stop on an invariant subspace); kernel
candidates = tridiagonal eigenpairs below 1e-9, Stieltjes vectors,
Gram-orthonormalized, each verified max_a || Q^a psi || < 1e-9.
Sector labels = the Gram of C2_link on the kernel basis; entanglement
S(pair|link) per sector vector by reshape-SVD.

**Measured (pre-registered).** Kernel dim (fusion prediction 2+2+2 = 6);
the sector split under C2_link (predicted {0: 2, 64: 2, lambda_210: 2});
the entropies: the r = 45 sector vectors at the committed intertwiner
value ln 45 = 3.80666... (V98 control, band 1e-6), r = 1 at 0, r = 210
recorded.

**Acceptance (pre-registered).** COMPOSITE-ALIVE iff dim = 6 with the
predicted split and the ln-45 control inside 1e-6; KERNEL-STILL-EMPTY iff
dim = 0; UNEXPECTED-KERNEL otherwise, recorded (dim, split, entropies).
pass unless UNEXPECTED.
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
from cl10 import engine as E  # noqa: E402
from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators,
                                   structure_constants)

NP = 45
LN = 45 + 1 + 210
NATOM = 1024 * LN          # 262144
LAM45 = math.log(45.0)
LN210 = math.log(210.0)
KMAX = 500
SEED = 0x6C
KTOL = 1e-9


def art(name):
    with open(exp(name, "results", name.lower() + ".json")) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks = []
    data = {}
    xp = E._E

    def log(msg):
        print("[V270 {:s} t+{:6.1f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)

    log("start, backend={}".format(E._BACKEND))

    # ---- T0: rebuild the V269 module (projector + parity), cross-check ----
    Ts_h = root_generators()
    f, _k, _as = structure_constants(Ts_h)
    Ls45_h = adjoint_link_operators(f)
    I45 = np.eye(NP, dtype=np.complex128)
    Ltot_h = [np.kron(L, I45) + np.kron(I45, L) for L in Ls45_h]
    C2t = np.zeros((NP * NP,) * 2, dtype=np.complex128)
    for L in Ltot_h:
        C2t += L @ L
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
    s1 = [b for b in binpar if b[1] == 1 and abs(b[0]) < 1e-9]
    s45 = [b for b in binpar if b[1] == 45 and b[2] < -0.5]
    s210 = [b for b in binpar if b[1] == 210]
    v269 = art("V269_dynamical_link_module")
    lam45 = s45[0][0] if len(s45) == 1 else float("nan")
    lam210 = s210[0][0] if len(s210) == 1 else float("nan")
    t0pass = bool(v269["pass"] and len(s1) == 1 and len(s45) == 1
                  and len(s210) == 1
                  and abs(lam45 - v269["data"]["c2_adjoint"]) < 1e-9
                  and abs(lam210 - v269["data"]["lambda_210"]) < 1e-9)
    checks.append({"name": "T0 rebuild equals the committed V269 artifact: "
                          "C2(45) {:.6f} vs {:.6f}, lambda_210 {:.6f} vs "
                          "{:.6f} (band 1e-9)".format(
                              lam45, v269["data"]["c2_adjoint"],
                              lam210, v269["data"]["lambda_210"]),
                   "pass": t0pass})
    data["lambda_210"] = lam210
    log("T0: lam210 {:.6f} (V269 {:.6f})".format(
        lam210, v269["data"]["lambda_210"]))

    # link generators embedded on 256 = [45 | 1 | 210]
    P210 = U[:, s210[0][3]:s210[0][4]]
    L256_h = []
    for a in range(NP):
        B = np.zeros((LN, LN), dtype=np.complex128)
        B[:45, :45] = Ls45_h[a]
        B[46:, 46:] = P210.conj().T @ Ltot_h[a] @ P210
        L256_h.append(B)
    C2l_h = sum(L @ L for L in L256_h)

    Ta = xp.asarray(np.stack(Ts_h))                      # (45,32,32)
    La = xp.asarray(np.stack(L256_h))                    # (45,256,256)
    I32 = np.eye(32, dtype=np.complex128)
    C2pair_h = np.zeros((1024, 1024), dtype=np.complex128)
    for T in Ts_h:
        Tpair = np.kron(T, I32) + np.kron(I32, T)
        C2pair_h += Tpair @ Tpair
    C2pair = xp.asarray(C2pair_h)
    C2l = xp.asarray(C2l_h)

    def S_of(psi):
        p = psi.reshape(32, 32, LN)
        v = xp.tensordot(C2pair, p.reshape(1024, LN), axes=(1, 0))
        v = v.reshape(32, 32, LN)
        v = v + xp.einsum("ml,pql->pqm", C2l, p)
        # two explicit contraction stages per cross term: the einsum
        # pathfinder can join the two 45-legs first and explode
        X = xp.einsum("aip,pql->aiql", Ta, p)      # (45,32,32,256)
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

    # ---- T1: Lanczos on S: kernel of the atom ----
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
            log("T1 Lanczos invariant subspace at k={:d} (beta {:.1e})".format(
                k + 1, b))
            break
        betas.append(b)
        V.append(w / b)
        if (k + 1) % 50 == 0:
            log("T1 Lanczos k={:d} (last alpha {:.4f}, beta {:.3e})".format(
                k + 1, a, b))
    Tm = np.diag(alphas[:k_used])
    for i in range(k_used - 1):
        Tm[i, i + 1] = betas[i]
        Tm[i + 1, i] = betas[i]
    tev, tvec = np.linalg.eigh(Tm)
    n_below = int(np.sum(tev < KTOL))
    seventh = float(tev[n_below]) if len(tev) > n_below else float("nan")
    data["lanczos_steps"] = k_used
    data["kernel_below_tol"] = n_below
    data["next_eigenvalue_above_kernel"] = seventh
    log("T1 Lanczos {} steps: eigenvalues < 1e-9: {}; next {}".format(
        k_used, n_below, seventh))

    Vm = xp.stack(V[:k_used]).T                       # (N, k_used)
    tdev = xp.asarray(tvec)
    cands = [Vm @ tdev[:, j] for j in range(k_used) if tev[j] < KTOL]
    resid = 0.0
    for psi in cands:
        for a in range(NP):
            resid = max(resid, float(xp.sqrt(xp.real(
                xp.vdot(Qa(psi, a), Qa(psi, a))))))
    # orthonormalize the candidates
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
    data["kernel_dim"] = kdim
    data["kernel_residual_max"] = resid
    t1pass = bool(kdim > 0 and resid < KTOL and kdim <= 9
                  and seventh > 1e-3)
    checks.append({"name": "T1 Lanczos kernel of S on the 262144 atom: dim "
                          "{:d} (fusion prediction 6), residual max_a||Q psi| "
                          "= {:.2e} < 1e-9, next eigenvalue {:.4f} above the "
                          "kernel".format(kdim, resid, seventh),
                   "pass": t1pass})
    log("T1 kernel dim {:d} (residual {:.2e})".format(kdim, resid))

    # ---- T2: sector split under C2_link + entropies ----
    split = {}
    entropies = []
    split_pass = True
    if kdim > 0:
        Gc = np.empty((kdim, kdim), dtype=np.complex128)
        CKer = [xp.einsum("ml,pql->pqm", C2l,
                          ker[i].reshape(32, 32, LN)).reshape(NATOM)
                for i in range(kdim)]
        for i in range(kdim):
            for j in range(kdim):
                Gc[i, j] = complex(xp.vdot(ker[i], CKer[j]))
        cev, cvec = np.linalg.eigh(0.5 * (Gc + Gc.conj().T))
        data["sector_eigenvalues"] = [float(x) for x in cev]
        log("T2 raw sector eigenvalues {}".format(
            [round(float(x), 6) for x in cev]))
        K = xp.stack(ker, axis=1)
        CK = xp.stack(
            [xp.einsum("ml,pql->pqm", C2l,
                       ker[i].reshape(32, 32, LN)).reshape(NATOM)
             for i in range(kdim)], axis=1)
        res = CK - K @ (K.conj().T @ CK)
        leak = float(xp.sqrt(xp.real(xp.vdot(res, res))))
        data["sector_leakage"] = leak
        log("T2 C2l leakage out of the kernel span {:.2e}".format(leak))
        clusters = {}
        for val in cev:
            near = min([0.0, lam45, lam210], key=lambda x: abs(val - x))
            clusters[near] = clusters.get(near, 0) + 1
        split = {"{:.4f}".format(k): v for k, v in clusters.items()}
        for s in range(kdim):
            psi = sum(cvec[i, s] * ker[i] for i in range(kdim))
            psi = psi / xp.sqrt(xp.real(xp.vdot(psi, psi)))
            M = E.asnumpy(psi).reshape(1024, LN)
            sv = np.linalg.svd(M, compute_uv=False)
            p = sv ** 2 / np.sum(sv ** 2)
            p = p[p > 1e-14]
            entropies.append(float(-np.sum(p * np.log(p))))
        data["sector_split"] = split
        data["entropies_by_sector"] = entropies
        exp_split = {0.0: 2, lam45: 2, lam210: 2}
        got = {}
        for val in cev:
            near = min(exp_split, key=lambda x: abs(val - x))
            got[near] = got.get(near, 0) + 1
        split_ok = (got == exp_split and abs(cev[0]) < 1e-8
                    and abs(cev[1]) < 1e-8)
        ent_ok = (entropies[0] < 1e-6 and entropies[1] < 1e-6
                  and abs(entropies[2] - LAM45) < 1e-6
                  and abs(entropies[3] - LAM45) < 1e-6
                  and abs(entropies[4] - LN210) < 1e-6
                  and abs(entropies[5] - LN210) < 1e-6)
        split_pass = bool(split_ok and ent_ok)
        log("T2 sectors {} entropies {}".format(
            split, [round(e, 6) for e in entropies]))
    checks.append({"name": "T2 sector split under C2_link: {} (predicted "
                          "{{0:2, {:.0f}:2, {:.0f}:2}}); entropies "
                          "{} (ln 45 = {:.6f} control at 1e-6, ln 210 = "
                          "{:.6f})".format(
                              split if split else "none", lam45, lam210,
                              [round(e, 6) for e in entropies] if entropies
                              else "none", LAM45, LN210),
                   "pass": bool(kdim == 6 and split_pass)})

    # ---- branch ----
    alive = bool(kdim == 6 and split_pass and resid < KTOL)
    if alive:
        branch = "COMPOSITE-ALIVE"
    elif kdim == 0:
        branch = "KERNEL-STILL-EMPTY"
    else:
        branch = "UNEXPECTED-KERNEL"
    ok = bool(branch != "UNEXPECTED-KERNEL")

    status = ("dynamical-link atom: pair (32x32) x link (1(+)45(+)210) "
              "dynamical module, neutral sector = ker sum_a (T_pair + L)^2 "
              "on 262144, Lanczos {} steps seed 0x6C, {:.0f} s: kernel dim "
              "{:d} (fusion prediction 6), residual {:.1e}, sector split {} "
              "under the link Casimir, entropies {} (ln 45 committed "
              "control, {:.1e}). Branch {}.".format(
                  k_used, time.time() - t0, kdim, resid,
                  split if split else "none",
                  [round(e, 6) for e in entropies] if entropies else "none",
                  data.get("kernel_residual_max", 0.0), branch))
    out = {"id": "V270_neutral_composite_sector",
           "tier": 1,
           "property": "the dynamical-link atom carries gauge-neutral "
                       "composites and the sector structure is the "
                       "predicted [2, 2, 2] (V270)",
           "residual": False,
           "programme": "V270",
           "branch": branch,
           "data": {**data, "lam45_value": lam45, "lam210_value": lam210,
                    "seed": "0x6C",
                    "fitted_inputs": [],
                    "backend": E._BACKEND,
                    "wall_time": time.time() - t0}}
    out["checks"] = checks
    out["status"] = status
    out["pass"] = bool(ok)
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v270_neutral_composite_sector.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
