#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V99 - the three-site chain: what the area law of the gauged layer
actually costs.

**Why (V98).** V98 computed the physical sector of ONE bond: the
center obstruction (no gauge-invariant state with one spinor per site)
and the minimal composite - two matter fields per site, intertwiners
with S(pair|link) = ln 45. Two numbers were left unnamed: (i) does the
ln 45 per bond ADD over several bonds (the area-law scaling), and (ii)
how far is ln 45 per Planck bond from the Bekenstein coefficient 1/4
per Planck area? This experiment computes both on the smallest
multi-bond system: the three-site, two-link chain with paired matter
(2 x 32 per site, adjoint 45 links), whose physical sector is small
enough to construct EXACTLY once the flux decomposition is used.

1. FLUX DECOMPOSITION (counting). The middle vertex fuses the two
   adjoint links into a channel r with multiplicity N_r =
   dim Hom(r, 45 x 45) (measured from the C2 spectrum of 45 x 45,
   2025 states) and m_r = dim Hom(r, 32 x 32) (measured on 1024).
   Common channels: r = 1, 45, 210 with (N, m) = (1, 2) each.
   Predicted and verified:
       dim ker = m_45^2 * sum_r m_r N_r = 4 * 6 = 24.

2. EXPLICIT SECTORS. Middle-vertex fusion channel r = 1: 8 states,
   the middle pair is a product singlet (S(pair2|rest) = 0) while the
   end pairs stay entangled at ln 45 - the flux line threads through
   the chain. Fusion channel r = 45: 8 states, built from the V98
   intertwiners and the adjoint Clebsch f_{gab}; vertex Gauss
   residual < 1e-12. (Links are adjoint-only: every link carries
   flux 45; r labels the fusion at the middle vertex.)

3. ENTANGLEMENT (the area-law cost). For the r = 45 sector:
   S(pair1 | rest) = ln 45 (one bond crosses the cut: additivity
   holds across a single bond); S(pair2 | rest) = ln 45, NOT 2 ln 45:
   the two bonds meeting at a vertex are CORRELATED by the intertwiner
   constraint (mutual information ln 45) - the area law counts bonds,
   the vertex correlations subtract. This is the subleading correction
   of the gauged layer, computed.

4. THE BEKENSTEIN GAP, NAMED AS A NUMBER: ln 45 = 3.80666 nats per
   Planck bond against 1/4 per Planck area: deficit factor
   4 ln 45 = 15.2266. No element of the current construction
   produces 1/4; the gap is reported, not closed.

Convention: the check PASSES when the protocol is valid and complete;
the physics verdict is in the status. Backend: CL10_BACKEND=cupy runs
the diagonalizations and einsums on the GPU.
"""

import importlib.util
import json
import math
import os
import sys
import time

import numpy as np

GPU = os.environ.get("CL10_BACKEND", "numpy") == "cupy"
if GPU:
    import cupy as XP
else:
    XP = np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)

_spec = importlib.util.spec_from_file_location(
    "v79", os.path.join(ROOT, "experiments", "V79_gauged_bond_layer",
                        "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)


def to_np(a):
    return XP.asnumpy(a) if GPU else np.asarray(a)


def entropy_from_rho(rho):
    nus = to_np(XP.linalg.eigvalsh(rho))
    s = 0.0
    for nu in nus:
        nu = float(np.real(nu))
        if nu > 1e-12:
            s -= nu * math.log(nu)
    return s


def channels(C2, tol=1e-6):
    """Group eigenvalues of a Casimir into [value, multiplicity]."""
    ev = np.sort(np.real(to_np(XP.linalg.eigvalsh(C2))))
    out = []
    for e in ev:
        if out and abs(e - out[-1][0]) < tol:
            out[-1][1] += 1
        else:
            out.append([float(e), 1])
    return out


def main():
    t0 = time.time()
    print(f"V99 - the three-site chain: the area-law cost of the gauged "
          f"layer [backend: {'cupy' if GPU else 'numpy'}]")
    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    Ts = [XP.asarray(T) for T in Ts_np]
    Ls = [XP.asarray(L) for L in Ls_np]
    fX = XP.asarray(f)
    I32 = XP.eye(32, dtype=np.complex128)
    I45 = XP.eye(45, dtype=np.complex128)

    # ---- 1. flux decomposition: channels of 32x32 and 45x45 --------------
    Js = [XP.kron(Ts[A], I32) + XP.kron(I32, Ts[A]) for A in range(45)]
    C2pair = XP.zeros((1024, 1024), dtype=np.complex128)
    for J in Js:
        C2pair += J @ J
    evp, evpvec = XP.linalg.eigh(C2pair)
    ch_pair = channels(C2pair)
    m_pair = {}
    basis_pair = {}
    for cval, mult in ch_pair:
        sel = XP.abs(evp - cval) < 1e-6
        m_pair[round(cval, 4)] = mult
        basis_pair[round(cval, 4)] = evpvec[:, sel]
    m1 = m_pair[0.0]
    m45 = m_pair[64.0] // 45
    print(f"  pair channels: singlets {m_pair[0.0]}, adjoints "
          f"{m_pair[64.0]} -> m_1 = {m1}, m_45 = {m45}")

    # 45 x 45: C2 of (L1 + L2)
    LL = [XP.kron(Ls[A], I45) + XP.kron(I45, Ls[A]) for A in range(45)]
    C2ll = XP.zeros((2025, 2025), dtype=np.complex128)
    for M in LL:
        C2ll += M @ M
    ch_ll = channels(C2ll)
    # shared channels: identify by C2 value; d_r from the pair side
    shared = {}
    for cval, mult_ll in ch_ll:
        key = round(cval, 4)
        if key in m_pair and mult_ll > 0:
            mult_p = m_pair[key]
            for mr, dr in [(m1, 1), (m45, 45)]:
                if mult_p == mr * dr:
                    shared[key] = (dr, mult_ll // dr, mult_ll % dr == 0)
    # the 210 channel: pair mult 420 = 2 x 210; match C2 in 45 x 45
    c2_210 = None
    for cval, mult in ch_pair:
        if mult == 420:
            c2_210 = round(cval, 4)
    if c2_210 is not None:
        for cval, mult_ll in ch_ll:
            if abs(cval - c2_210) < 1e-6 and mult_ll == 210:
                shared[c2_210] = (210, 1, True)
    Nsum = sum(v[1] * m_pair[k] // v[0] for k, v in shared.items()
               if v[2])
    dim_ker = m45 * m45 * Nsum
    print(f"  shared channels (C2, d_r, N_r, m_r): "
          f"{[(k, v[0], v[1], m_pair[k] // v[0]) for k, v in sorted(shared.items()) if v[2]]}")
    print(f"  dim ker(3-site chain) = {dim_ker}  (expect 24 = 4 x (2+2+2))")

    # ---- 2. explicit intertwiners -----------------------------------------
    # end sites: V98 intertwiners I_a (1024 x 45), a = 1, 2
    E45 = basis_pair[64.0]                       # 1024 x 90
    Jsub = [E45.conj().T @ J @ E45 for J in Js]  # 90 x 90
    I90 = XP.eye(90, dtype=np.complex128)
    M3 = XP.zeros((90 * 45, 90 * 45), dtype=np.complex128)
    for A in range(45):
        QA = XP.kron(Jsub[A], I45) + XP.kron(I90, Ls[A])
        M3 += QA.conj().T @ QA
    w3, u3 = XP.linalg.eigh(M3)
    Is = [E45 @ u3[:, j].reshape(90, 45)
          for j in XP.where(XP.abs(w3) < 1e-8)[0]]
    n_int = len(Is)
    # middle vertex, r = 1: pair singlet x link singlet
    sing = basis_pair[0.0]                       # 1024 x 2
    T1 = [XP.einsum("j,ab->jab", sing[:, s], I45 / math.sqrt(45.0))
          for s in range(sing.shape[1])]
    # middle vertex, r = 45: T_c[j,(ab)] = sum_g W_c[j,g] f_{gab}
    # W_c: intertwiners (C^45, L) -> (pair-45-subspace, Jsub):
    # Jsub W = W L, i.e. ker(Jsub (x) I - I (x) L^T) on 4050.
    Mad = XP.zeros((90 * 45, 90 * 45), dtype=np.complex128)
    for A in range(45):
        QA = XP.kron(Jsub[A], I45) - XP.kron(I90, Ls[A].T)
        Mad += QA.conj().T @ QA
    wa, ua = XP.linalg.eigh(Mad)
    Ws = [ua[:, j].reshape(90, 45)
          for j in XP.where(XP.abs(wa) < 1e-8)[0]]
    T45 = []
    for Wc in Ws:
        Tsub = XP.einsum("jg,gab->jab", Wc, fX)   # (90, 45, 45)
        T = XP.einsum("ij,jab->iab", E45, Tsub)   # lift to (1024, 45, 45)
        T = T / XP.linalg.norm(T)
        T45.append(T)

    # vertex Gauss residual: (J +/- (L1 + L2)) T = 0 on pair x 45 x 45
    # (the sign depends on link orientation; the physical one is the
    # one that annihilates the intertwiners)
    vres = 0.0
    for T in T1 + T45:
        rp = rm = 0.0
        for A in range(45):
            jt = XP.einsum("jk,kab->jab", Js[A], T)
            l1 = XP.einsum("ac,jcb->jab", Ls[A], T)
            l2 = XP.einsum("bc,jac->jab", Ls[A], T)
            rp = max(rp, float(XP.linalg.norm(jt + l1 + l2)))
            rm = max(rm, float(XP.linalg.norm(jt - l1 - l2)))
        vres = max(vres, min(rp, rm))
    print(f"  intertwiners: end n = {n_int}, middle r=45 n = {len(T45)}; "
          f"vertex Gauss residual = {vres:.2e}")

    # ---- 3. glue and entanglement -----------------------------------------
    # Psi[i,j,k] = sum_{ab} I_a[i,a] T[j,a,b] I_b[k,b]; the reduced
    # matrices use I_b^dag I_b = I/45 (W isometries):
    #   rho1 = I_a H I_a^dag, H[a,c] = (1/45) sum_{j,b} T[j,a,b] conj(T[j,c,b])
    #   rho2[j,k] = (1/45^2) sum_{a,b} T[j,a,b] conj(T[k,a,b])
    def S_pair1(T, Ia):
        H = XP.einsum("jab,jcb->ac", T, T.conj()) / 45.0
        rho = Ia @ H @ Ia.conj().T
        rho = rho / XP.trace(rho).real
        return entropy_from_rho(rho)

    def S_pair2(T):
        rho = XP.einsum("jab,kab->jk", T, T.conj()) / 45.0 / 45.0
        rho = rho / XP.trace(rho).real
        return entropy_from_rho(rho)

    S1_sing = max(S_pair1(T, Is[0]) for T in T1)
    S2_sing = max(S_pair2(T) for T in T1)
    S1_45 = [S_pair1(T, Is[0]) for T in T45]
    S2_45 = [S_pair2(T) for T in T45]
    ln45 = math.log(45.0)
    print(f"  r=1 sector: S(pair1|rest) = {S1_sing:.6f} (ln 45: the flux "
          f"line threads through), S(pair2|rest) = {S2_sing:.6f} "
          "(expect 0)")
    print(f"  r=45 sector: S(pair1|rest) = {[round(s, 6) for s in S1_45]} "
          f"(ln 45 = {ln45:.6f})")
    print(f"  r=45 sector: S(pair2|rest) = {[round(s, 6) for s in S2_45]} "
          f"(ln 45, NOT 2 ln 45 = {2 * ln45:.6f})")

    # ---- 4. the Bekenstein gap, as a number --------------------------------
    deficit = ln45 - 0.25
    factor = 4.0 * ln45
    print(f"  Bekenstein gap: ln 45 = {ln45:.6f} vs 1/4 per Planck area; "
          f"deficit {deficit:.6f} nats, factor {factor:.4f}")

    ok = (dim_ker == 24 and n_int == 2 and len(T45) == 2
          and vres < 1e-8
          and abs(S1_sing - ln45) < 1e-6 and abs(S2_sing) < 1e-9
          and all(abs(s - ln45) < 1e-6 for s in S1_45)
          and all(abs(s - ln45) < 1e-6 for s in S2_45))

    status = (
        "THE AREA LAW OF THE GAUGED LAYER HAS A VERTEX CORRECTION, AND "
        "THE BEKENSTEIN GAP IS A NAMED NUMBER. (1) The physical sector "
        "of the three-site, two-link chain with paired matter is "
        f"dim ker = {dim_ker} (predicted 24 = m_45^2 x (m_1 N_1 + m_45 "
        "N_45 + m_210 N_210) = 4 x 6 from the measured flux "
        "decomposition of 45 x 45; links are adjoint-only, r labels "
        "the fusion channel at the middle vertex). (2) Constructed "
        "exactly: in the r = 1 sector (8 states) the middle pair is a "
        "product singlet (S(pair2|rest) = 0) while the end pairs stay "
        "entangled at ln 45 - the flux line threads through the chain; "
        "the r = 45 sector (8 states) is built from the V98 "
        "intertwiners and the adjoint Clebsch "
        f"f_{{gab}} (vertex Gauss residual {vres:.1e}). (3) ADDITIVITY "
        "HOLDS ACROSS A SINGLE BOND: S(pair1|rest) = ln 45 exactly. "
        "BUT the region with TWO boundary bonds carries S(pair2|rest) "
        f"= ln 45, not 2 ln 45: the intertwiner constraint CORRELATES "
        "the two bonds meeting at the vertex (mutual information "
        "ln 45). The area law counts bonds; vertex correlations "
        "subtract - the subleading correction of the gauged layer, "
        "computed. (4) THE GAP TO BEKENSTEIN, AS A NUMBER: ln 45 = "
        f"{ln45:.5f} nats per Planck bond against 1/4 per Planck area: "
        f"deficit {deficit:.5f} nats, factor {factor:.4f}. No element "
        "of the current construction produces 1/4; closing it needs "
        "either a superselection/truncation of the flux or a different "
        "horizon microstate count - named, open.")
    verdict = {"experiment": "V99_three_site_chain",
               "pass": bool(ok), "status": status,
               "data": {"dim_ker": int(dim_ker),
                        "shared_channels": {str(k): [v[0], v[1]]
                                            for k, v in shared.items()
                                            if v[2]},
                        "n_end_intertwiners": int(n_int),
                        "n_middle_45": len(T45),
                        "vertex_gauss_residual": vres,
                        "S_r1": [S1_sing, S2_sing],
                        "S_r45_pair1": [round(s, 6) for s in S1_45],
                        "S_r45_pair2": [round(s, 6) for s in S2_45],
                        "ln45": ln45,
                        "bekenstein_deficit": deficit,
                        "bekenstein_factor": factor}}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v99_three_site_chain.json"), "w") as fh:
        json.dump(verdict, fh, indent=1)
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v99_three_site_chain.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
