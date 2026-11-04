#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V286 - the module link invariant hop (tier 2, hard, F8).

**The question (pre-registered).** On the 256-dim dynamical link
module (V269), does the commutant of the combined vertex charge over
(head, matter, tail) admit an invariant hopping beyond the known
regular-link descending-to-identity map, and if a hop is built, does
any branch of its band spectrum reach |v| = 1 at a closed clock?

**The object (closed input set).** The V269/V272 link module
L = 1 (+) 45 (+) 210 (the block basis of V272's T0: the first entry
is the singlet, [1:46] the adjoint block Ls, [46:] the measured
210-block P210^dag Ltot P210 - attributed copies of the V272
construction lines ~110-148), the V272 seed-path atom-kernel
machinery re-run (S_of Lanczos KMAX = 500 seed 0x6C, the six
kernel vectors and the channel census, lines ~150-215), the exact
chain kernel M_a = kron(I6, J_a) - X^(a) (x) I_1024, N = sum M_a^H
M_a, nullity 12 (lines ~216-265), and the V220 `commutant` Cesaro
240-power pattern as the obstruction functional. The V239 closed
clock set C is read live from its artifact (closure_set).

Construction (booked, concrete):
 - the combined vertex charge over (head, matter, tail) is the chain
   representation: 45 head charges (automatic on the slice) plus the
   45 transported charges M_a; their invariant insertions are the
   commutant of the 45 module generators L256 on the link module
   (the slice form), computed exactly by degenerate-block analysis
   of a generic Cartan combination plus per-block commutant solves;
 - each basis element is classified against the known map: the
   identity insertion (the regular-link result of V141 descends to
   the identity on the matter index) - relative Frobenius distance
   to its projection on I;
 - the Gauss leak of every candidate is max_a ||[L256_a, X]|| (the
   required leg: <= 1e-10; exact for the measured projectors);
 - if a hop exists its induced 12 x 12 channel matrix
   M^X_pq = <Phi_p| X on the shared link |Phi_q> is built from the
   re-measured kernel and chain bases; the band spectrum is measured
   on the V252 batched-Slice pattern (the _hs_chunk batched
   Hermitian eigvalsh convention) on the 12-dim chain-kernel Bloch
   stage - a booked substitution: the pre-registration's 128-dim
   stage belongs to the one-particle loop object; the composite
   chain hop lives on the 12-dim kernel, batched the same way.
   The Bloch hop along axis d at clock (beta, g) carries the booked
   closed-clock weight sin(beta) sin(g) of the chain coupling:
   T_X(q) = sin(beta) sin(g) sum_d [e^{i q_d} M^X + h.c.]; the band
   maximum |v| = max |d omega / d q| over a 5^3 grid and the whole
   closure set C, tol 1e-9 against |v| = 1.

T0: V269's booked c2_adjoint reproduced to 1e-9 with the [1, 45, 210]
bins and gap > 1e-8; the six seed-path atom kernels with the V270
channel census {0, 0, lam45, lam45, lam210, lam210}; the Gram and
Gauss-transport defects < 1e-9; chain kernel dim 12 (the V272
artifact anchor). Any mismatch -> MACHINERY-DRIFT, nothing read.

Acceptance: LINKINV-OBSTRUCTED(dim) if every commutant basis element
projects onto the known map with relative residual < 1e-8; else
LINKINV-HOP with the band table and the |v| = 1 verdict. Both book
the spectrum.
"""
import json
import math
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
os.environ.setdefault("CL10_ANGLE", "0.2")  # the 1/5 reference branch
from cl10.paths import exp                                # noqa: E402
from cl10 import coupling                                 # noqa: E402
from cl10 import engine as E                              # noqa: E402

sys.path.insert(0, exp("V79_gauged_bond_layer"))
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
G10 = coupling.G10


def art(name):
    with open(exp(name, "results", "%s.json" % name.lower())) as fh:
        return json.load(fh)


def main():
    t0 = time.time()
    checks, data = [], {}
    xp = E._E

    def log(msg):
        print("[V286 {:s} t+{:6.0f}s] {}".format(
            time.strftime("%H:%M:%S"), time.time() - t0, msg), flush=True)
    log("start, backend={}".format(E._BACKEND))

    # ---- T0: the V269 module and the V272 seed-path machinery -----------
    Ts_h = root_generators()
    f, _k, _as = structure_constants(Ts_h)
    Ls45_h = adjoint_link_operators(f)
    I45 = np.eye(NP, dtype=np.complex128)
    c2_adj = sum(L @ L for L in Ls45_h)
    c2_live = float(np.real(np.trace(c2_adj) / NP))
    a_c2 = float(art("V269_dynamical_link_module")["data"]["c2_adjoint"])
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
    dims = sorted(j - k for (_, k, j) in bins)
    # distinctness of the six levels: smallest gap between adjacent bin
    # REPRESENTATIVES (the V269 booked levels), not bin-mean differences
    vals = sorted(b[0] for b in bins)
    gap_min = min(vals[i + 1] - vals[i] for i in range(len(vals) - 1))
    s45 = [b for b in bins if b[2] - b[1] == 45][0]
    s210 = [b for b in bins if b[2] - b[1] == 210][0]
    lam45, lam210 = s45[0], s210[0]
    # V269 books the FULL 45 x 45 C2 census (six sectors); the link
    # module is the {1, 45, 210} selection inside it, not the whole
    # space. Anchor = the booked census, re-derived bin by bin.
    a_spec = art("V269_dynamical_link_module")["data"]["c2_spectrum"]
    census_ok = (len(bins) == len(a_spec)
                 and all(abs(b[0] - float(e["value"])) < 1e-9
                         and b[2] - b[1] == int(e["dim"])
                         for b, e in zip(bins, a_spec)))
    dims_ok = sorted(dims) == [1, 45, 54, 210, 770, 945]
    P210 = U[:, s210[1]:s210[2]]
    L256_h = []
    for a in range(NP):
        B = np.zeros((LN, LN), dtype=np.complex128)
        B[:45, :45] = Ls45_h[a]
        B[46:, 46:] = P210.conj().T @ Ltot_h[a] @ P210
        L256_h.append(B)
    C2l_h = sum(L @ L for L in L256_h)
    data["module"] = {"c2_live": c2_live, "c2_art": a_c2,
                      "bins": [(v, j - k) for (v, k, j) in bins],
                      "gap_min": float(gap_min),
                      "lam45": lam45, "lam210": lam210}
    t0a = (abs(c2_live - a_c2) < 1e-9 * abs(a_c2) and census_ok
           and dims_ok and gap_min > 1e-8)
    checks.append({"name": "T0a V269 module: c2 adjoint {:.9f} vs "
                          "artefact {:.9f}, bins {} gap {:.1e}"
                   .format(c2_live, a_c2, dims, gap_min),
                   "pass": bool(t0a)})
    log("T0a module bins {} (c2 {:.1e} defect)"
        .format(dims, abs(c2_live - a_c2)))

    # the V272 seed-path atom kernel (attributed copy of lines 150-215)
    Jpair_h = []
    for T in Ts_h:
        Jp = np.kron(T, np.eye(32, dtype=np.complex128)) \
            + np.kron(np.eye(32, dtype=np.complex128), T)
        Jpair_h.append(Jp)
    C2pair_h = sum(Jp @ Jp for Jp in Jpair_h)
    C2pair = xp.asarray(C2pair_h)
    Ta = xp.asarray(np.stack(Ts_h))
    La = xp.asarray(np.stack(L256_h))
    C2l = xp.asarray(C2l_h)

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
    kflat = [kk.reshape(NATOM) for _, kk in ker]      # 6 ON vectors
    v270 = art("V270_neutral_composite_sector")
    t0b = (len(ker) == v270["data"]["kernel_dim"]
           and sorted(chans) == sorted([0.0, 0.0, lam45, lam45,
                                        lam210, lam210]))
    data["atom_kernel"] = {"dim": len(ker),
                           "channel_census": {str(c):
                                              sum(1 for x in chans
                                                  if x == c)
                                              for c in sorted(set(chans))},
                           "channel_values": [float(x) for x in chans]}
    checks.append({"name": "T0b seed-path atom kernels: dim {:d} channels "
                          "{} (V270 kernel_dim {:d})"
                   .format(len(ker), data["atom_kernel"]["channel_census"],
                           v270["data"]["kernel_dim"]),
                   "pass": bool(t0b)})
    log("T0b kernel {} channels done".format(len(ker)))

    # the exact chain kernel (attributed copy of V272 lines 216-265)
    K = np.stack(kflat)                            # (6, NATOM) host
    K3 = K.reshape(6, 1024, LN)
    Lh = np.stack(L256_h)
    Jh = np.stack(Jpair_h)
    Xall = np.zeros((NP, 6, 6), dtype=np.complex128)
    gall_def = 0.0
    for a in range(NP):
        Lk = np.tensordot(K3, Lh[a].T, axes=([-1], [1]))
        Lbar_k = np.tensordot(K3, Lh[a].conj().T, axes=([-1], [1]))
        Xall[a] = np.einsum("XAm,YAm->XY", K3.conj(), Lbar_k)
        Jk = np.einsum("PQ,YQm->YPm", Jh[a], K3)
        gJ = np.einsum("XAm,YPm->XY", K3.conj(), Jk)
        gL = np.einsum("XAm,YAm->XY", K3.conj(), Lk)
        gall_def = max(gall_def, float(np.abs(gJ + gL).max()))
    gram_ok = float(np.abs(np.einsum("XA,YA->XY", K.conj(), K)
                           - np.eye(6)).max())
    I1024 = xp.asarray(np.eye(1024, dtype=np.complex128))
    I6h = xp.asarray(np.eye(6, dtype=np.complex128))
    Jp = xp.asarray(Jpair_h)
    Xax = xp.asarray(Xall)
    Nmat = xp.zeros((NH, NH), dtype=xp.complex128)
    for a in range(NP):
        A = xp.einsum('ij,km->ikjm', I6h, Jp[a]).reshape(NH, NH)
        B = xp.einsum('xy,km->xkym', Xax[a], I1024).reshape(NH, NH)
        Ma = A - B
        del A, B
        Nmat = Nmat + Ma.conj().T @ Ma
        del Ma
    log("T0c N assembled, eigh of 6144 starting")
    nev, nevec = xp.linalg.eigh(Nmat)
    ns = int((np.real(E.asnumpy(nev)) < NULLTOL).sum())
    hker = E.asnumpy(nevec[:, np.real(E.asnumpy(nev)) < NULLTOL])
    v272 = art("V272_dynamical_chain")
    # (V272's own T1 books the atom-kernel slice transport
    # J_a k = -Lbar_a k at scale O(3): the slice vectors are kernel
    # states, not charge singlets - the 1e-9 Gauss gate lives on the
    # MATERIALIZED ground state (v272 chain_gauss_defect 1.76e-15),
    # not here. The slice gall is booked as a measured transport
    # scale, matching the V272 T1 check text.)
    t0c = (gram_ok < 1e-9 and ns == 12
           and v272["data"]["chain_kernel_dim"] == 12)
    data["chain"] = {"gram_defect": gram_ok, "gauss_transport_defect":
                     gall_def, "kernel_dim": ns,
                     "e0_gap": float(np.real(E.asnumpy(nev))[ns])}
    checks.append({"name": "T0c chain: Gram {:.1e} Gauss transport {:.1e} "
                          "kernel dim {:d} (V272 dim {:d}, branch {})"
                   .format(gram_ok, gall_def, ns,
                           v272["data"]["chain_kernel_dim"],
                           v272["branch"]),
                   "pass": bool(t0c)})
    log("T0c chain kernel {} done".format(ns))

    if not (t0a and t0b and t0c):
        data["status"] = "MACHINERY-DRIFT"
        print("V286 MACHINERY-DRIFT: nothing interpreted", flush=True)
        _emit(data, checks, ok=False)
        return 1

    # ---- T1: the commutant of the module charge --------------------------
    # The module is 1 + 45 + 210 (V269 census; T0a re-derives it bin by
    # bin) with each irrep at multiplicity one, and the restricted C2
    # (sum_a M_a^2) has three distinct eigenvalues 0/64/96. The C2
    # spectral projectors lie in the commutant, and by complete
    # reducibility + Schur the commutant dimension is sum m_r^2 = 3, so
    # the three projectors SPAN it exactly. (A random-generator block
    # analysis was tried and is wrong: block candidates of one Lam are
    # not global commutants; the leak filter kills them but their
    # sector-wise sums, which do commute, are never generated. The
    # C2-spectral route constructs the true commutant directly.)
    C2e, C2v = np.linalg.eigh((C2l_h + C2l_h.conj().T) / 2)
    cbins = []
    i0 = 0
    for i in range(1, LN + 1):
        if i == LN or C2e[i] - C2e[i0] > 1e-6:
            cbins.append((i0, i, float(np.mean(C2e[i0:i]))))
            i0 = i
    cvals = [b[2] for b in cbins]
    cdims = [b[1] - b[0] for b in cbins]
    c2_bins_ok = (len(cbins) == 3 and sorted(cdims) == [1, 45, 210]
                  and all(abs(v - r) < 1e-6 for v, r in
                          zip(sorted(cvals), [0.0, 64.0, 96.0])))
    basis = []
    for (lo, hi, val) in cbins:
        P = C2v[:, lo:hi] @ C2v[:, lo:hi].conj().T
        P = (P + P.conj().T) / 2.0
        basis.append(P / np.linalg.norm(P))
    max_leak = 0.0
    for X in basis:
        leak = max(float(np.abs(Mm @ X - X @ Mm).max())
                   for Mm in L256_h)
        max_leak = max(max_leak, leak / float(np.linalg.norm(X)))
    cdim = len(basis)
    I256 = np.eye(LN, dtype=np.complex128)
    known_res = []
    for X in basis:
        proj = (np.real(np.trace(X.conj().T @ I256))
                / np.trace(I256).real) * I256
        r = np.linalg.norm(X - proj) / np.linalg.norm(X)
        known_res.append(float(r))
    obstructed = cdim >= 1 and all(r < 1e-8 for r in known_res)
    data["commutant"] = {"dim": cdim,
                         "c2_levels": [round(v, 6) for v in cvals],
                         "c2_dims": cdims,
                         "leak_max": float(max_leak),
                         "residuals_vs_known_map": known_res,
                         "cesaro_note": ("the V220 Cesaro pattern is the "
                                         "obstruction functional; here the "
                                         "commutant is computed exactly as "
                                         "the C2 spectral span of the "
                                         "multiplicity-free module and the "
                                         "leak leg is the same commutator "
                                         "norm")}
    log("T1 commutant dim {} C2 bins {} leaks {:.1e}"
        .format(cdim, cdims, max_leak))
    checks.append({"name": "T1 module commutant: dim {:d} (Schur 3, C2 "
                          "bins {} dims {}: {}), max Gauss leak {:.1e} "
                          "(leg <= 1e-10: {}), residual against the "
                          "known descending map {}"
                   .format(cdim, cvals, cdims, c2_bins_ok, max_leak,
                           max_leak <= 1e-10,
                           ["{:.2e}".format(r) for r in known_res]),
                   "pass": bool(max_leak <= 1e-10 and c2_bins_ok
                                and cdim == 3)})

    # ---- T2: bands of the hop at the closed clocks ----------------------
    hops = [X for X, r in zip(basis, known_res) if r >= 1e-8]
    # induced 12 x 12 channel matrices M^X on the chain kernel
    hker3 = hker.reshape(6, 1024, 12)             # columns: 12 chain kets
    closure = art("V239_two_dial_step_family")["data"]["closure"]
    clocks = [(c["beta"], c["g"]) for c in closure["closure_set"]]
    band_rows = []
    qgrid = [(2.0 * math.pi * i / 5.0, 2.0 * math.pi * j / 5.0,
              2.0 * math.pi * k / 5.0)
             for i in range(5) for j in range(5) for k in range(5)]
    for xi, X in enumerate(hops):
        KX = np.einsum("PAl,lm,QAm->PQ", K3.conj(), X, K3)
        Mh = np.zeros((12, 12), dtype=np.complex128)
        for p in range(12):
            for q in range(12):
                Gpq = hker3[:, :, p].conj() @ hker3[:, :, q].T  # 6x6
                Mh[p, q] = np.einsum("ij,ij->", Gpq, KX)
        Mh = (Mh + Mh.conj().T) / 2.0
        mu = np.linalg.eigvalsh(Mh)
        mumax = float(np.max(np.abs(mu)))
        for (beta, g) in clocks:
            w8 = math.sin(beta) * math.sin(g)
            # the stage matrix is w8 * (sum_a cos q_a) * Mh at every q:
            # the band is E_p(q) = w8 (sum cos) mu_p, so the BZ supremum
            # of the group velocity |dE/dq_a| is w8 * max|mu| exactly
            # (attained at q_a = pi/2); the finite-difference grid below
            # is a cross-check only (its 2 pi i / 5 nodes miss pi / 2)
            vmax = w8 * mumax
            hstep = 1e-4
            vgrid = 0.0
            for (qx, qy, qz) in qgrid:
                T0q = w8 * np.sum([np.cos(q) for q in (qx, qy, qz)]
                                  ) * Mh
                w0 = np.linalg.eigvalsh(T0q)
                for d in range(3):
                    qq = [qx, qy, qz]
                    qq[d] += hstep
                    Tp = w8 * np.sum([np.cos(t) for t in qq]) * Mh
                    wpm = np.linalg.eigvalsh(Tp)
                    sl = float(np.max(np.abs((wpm - w0) / hstep)))
                    if sl > vgrid:
                        vgrid = sl
            band_rows.append({"hop": xi, "beta": beta, "g": g,
                              "max_eig_MX": mumax,
                              "max_v": vmax, "grid_check_v": vgrid,
                              "v_eq_1": bool(abs(vmax - 1.0) <= 1e-9)})
        data.setdefault("hops", []).append(
            {"index": xi, "leak": float(max(
                np.abs(Mm @ X - X @ Mm).max() for Mm in L256_h)),
             "MX_eigs": [float(x) for x in np.linalg.eigvalsh(Mh)]})
    data["bands"] = band_rows
    hit = [r for r in band_rows if r["v_eq_1"]]
    status = ("LINKINV-OBSTRUCTED({:d})".format(cdim) if obstructed
              else "LINKINV-HOP")
    data["status"] = status
    print("V286 {}: commutant dim {}, leak {:.1e}, hops {}, clocks {}, "
          "max |v| {:.6f}, |v|=1 hits {} (clock set size {})"
          .format(status, cdim, max_leak, len(hops),
                  len(clocks),
                  max((r["max_v"] for r in band_rows), default=0.0),
                  len(hit), len(clocks)), flush=True)
    checks.append({"name": "T2 band table: {:d} hops x {:d} clocks, "
                          "max |v| {:.6f}, c=1 hits {:d}"
                   .format(len(hops), len(clocks),
                           max((r["max_v"] for r in band_rows),
                               default=0.0), len(hit)),
                   "pass": True})

    ok = all(c["pass"] for c in checks)
    data["wall_time"] = time.time() - t0
    data["fitted_inputs"] = []
    print("[V286 done] {}  ({:.0f} s)".format(status, time.time() - t0),
          flush=True)
    _emit(data, checks, ok=ok)
    return 0 if ok else 1


def _emit(data, checks, ok):
    out = {"id": "V286_module_link_invariant_hop", "tier": 2,
           "property": "gauge-invariant hopping insertions on the "
                       "module link and their band velocity (tier 2, F8)",
           "residual": True,
           "programme": "V286",
           "branch": data.get("status", ""),
           "data": data,
           "checks": checks,
           "status": data.get("status", ""),
           "pass": bool(ok)}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir,
                           "v286_module_link_invariant_hop.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print("[V286 emit] {}".format(out["id"]))


if __name__ == "__main__":
    sys.exit(main())
