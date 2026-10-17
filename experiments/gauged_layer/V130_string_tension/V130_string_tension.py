#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V130 - does the gauged bond layer confine, and at what tension?
(residual I: the string tension and the confining phase)

**Why.**  Section 13 of the paper states that the QCD area law "is a property
of a confining phase that the theory has not shown it is in", and that
Lambda_QCD is not on the loop ladder (it would need n ~ 11.6). V24/V25
measured the closed-loop suppression of the STEP: the n-hop loop (MG)^n falls
geometrically, ratio r ~ 0.54 per hop - a perimeter law, since the suppression
counts hops and not enclosed area. V79 gauged the bond layer into a lattice
gauge theory with structure constants f_{ABC} of so(10); V98 found that its
link Hilbert space is the ADJOINT 45 and that an unpaired spinor admits no
gauge-invariant state at all (dim ker = 0, gap 146); V99 built the three-site
chain with paired matter and found dim ker = 24.

The layer therefore has a Kogut-Susskind electric Hamiltonian
    H = sum_links C2(link),  C2 = sum_A L_A L_A
and a Gauss constraint, which is all that is needed to ask whether a static
colour source feels a LINEAR potential. The answer is fixed by one structural
fact: the link space is the adjoint 45, so the trivial flux channel is NOT
available on a link. Every link of a gauge-invariant chain therefore carries
adjoint electric flux, and the energy of a chain of R links is additive:

    E(R) = R * C2(45),   string tension sigma = C2(45).

This experiment measures that, exactly, on the implemented layer:

1. THE LINK CASIMIR. C2 on the 45-dimensional link space is a single
   eigenvalue (the link is an irreducible adjoint): measured, and the
   invariant pairing of 45 (x) 45 (the Killing form) is an exact intertwiner,
   which is what lets adjacent links satisfy the vertex constraint.
2. THE FLUX IS FORCED. The link space is the adjoint 45 and has NO invariant
   vector (ker sum_A L_A^dagger L_A = 0), so a link cannot be flux-free; since
   C2 is c2_45 x I on the whole link space, the electric energy of a chain of
   R links is exactly R C2(45): the static potential is LINEAR. The vertex
   Gauss Laplacian on (pair-adjoint) (x) link has kernel m_45 (the V98/V99
   intertwiner space) and gap C2(45), the same number.
3. THE TENSION AND ITS SCALE. sigma = C2(45) = 64 in the normalisation of the
   root generators (h_a^2 = I), i.e. 16 in the physics normalisation of the
   paper: 64 lattice units of electric energy per Planck link = 64 M_P =
   7.8e20 GeV per Planck length = 64 M_P^2 = 9.5e39 GeV^2, against the QCD
   string tension 0.176 GeV^2: the layer confines, 5e40 times too hard.
   Lambda_QCD is not its scale.
4. NO AREA LAW IN THE STEP ITSELF. The six stage operators of the implemented
   bond factor act on one sub-index bit each, so stages of different axes
   commute exactly (measured): the step has no plaquette commutator, no
   non-Abelian field strength: every elementary plaquette of the step is the
   identity and <W_square> = 1. The confining behaviour belongs to the gauged
   layer alone.

Everything is exact linear algebra on V79's adjoint link operators, V99's
kernel and the implemented supercell stages; no Monte Carlo, no truncation.
"""

import importlib.util
import json
import os
import sys
import time

import numpy as np


HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, ROOT)
from cl10.paths import exp

from cl10 import engine as E                   # noqa: E402

GPU = E._BACKEND == "cupy"   # the engine decides; cl10/config.yaml decides the engine
if GPU:
    import cupy as XP
else:
    XP = np

_spec = importlib.util.spec_from_file_location(
    "v79", exp("V79_gauged_bond_layer",
                        "V79_gauged_bond_layer.py"))
V79 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(V79)

D_L = 45


def to_np(a):
    return XP.asnumpy(a) if GPU else np.asarray(a)


def main():
    t0 = time.time()
    print(f"V130 - the string tension of the gauged bond layer [backend: "
          f"{'cupy' if GPU else 'numpy'}]")
    out = {"experiment": "V130_string_tension", "data": {}}
    checks = []

    Ts_np = V79.root_generators()
    f, _k, _asym = V79.structure_constants(Ts_np)
    Ls_np = V79.adjoint_link_operators(f)
    Ls = [XP.asarray(L) for L in Ls_np]
    I45 = XP.eye(D_L, dtype=np.complex128)
    I32 = XP.eye(32, dtype=np.complex128)

    # ---- 1. the link Casimir and the vertex pairing -----------------------
    C2link = XP.zeros((D_L, D_L), dtype=np.complex128)
    for A in range(45):
        C2link += Ls[A] @ Ls[A]
    w = np.real(to_np(XP.linalg.eigvalsh(C2link)))
    c2_45 = float(np.mean(w))
    c2_spread = float(w.max() - w.min())
    fX = XP.asarray(f)
    K = XP.einsum("acd,bcd->ab", fX, fX)          # the Killing form
    Kn = K / XP.linalg.norm(K)
    Sing = XP.zeros((D_L * D_L, D_L * D_L), dtype=np.complex128)
    for A in range(45):
        QA = XP.kron(Ls[A], I45) + XP.kron(I45, Ls[A])
        Sing += QA.conj().T @ QA        # the Casimir of 45 (x) 45: ker = singlets
    pres = float(np.abs(to_np(Sing @ Kn.reshape(-1))).max())
    ws = np.real(to_np(XP.linalg.eigvalsh(Sing.conj().T @ Sing)))
    n_sing = int(np.sum(ws < 1e-8))
    fmax = max(float(np.abs(to_np(Ls[a] @ Ls[b] - Ls[b] @ Ls[a])).max())
               for a in range(6) for b in range(6))
    out["data"]["link"] = {"C2": c2_45, "spread": c2_spread,
                           "pairing_residual": pres, "n_singlets": n_sing,
                           "killing_trace": float(np.real(to_np(
                               XP.trace(K)))),
                           "link_commutator_max": fmax}
    checks.append(("T1 the link space is an irreducible adjoint: C2 = sum_A "
                   "L_A^2 has a single eigenvalue %.6f (spread %.2e over the "
                   "45 states); the Killing form K_ab = f_acd f_bcd is "
                   "annihilated by the singlet operator of 45 (x) 45 to %.2e "
                   "and the singlet space has dimension %d, so adjacent links "
                   "pair into it uniquely; the link operators do not commute "
                   "(max [L_a,L_b] = %.3f) - the layer carries a non-Abelian "
                   "field strength and the trivial flux channel is "
                   "unavailable on a link"
                   % (c2_45, c2_spread, pres, n_sing, fmax),
                   c2_spread < 1e-9 and pres < 1e-10 and n_sing == 1
                   and fmax > 1.0))

    # ---- 2. a link cannot be flux-free: the flux is forced ----------------
    # The link space is the adjoint 45, so a flux-free link would be an
    # invariant vector of that representation: ker sum_A L_A^dagger L_A =
    # ker(C2(45) I) = 0. A link therefore CANNOT be flux-free, and since T1
    # measured C2 = c2_45 x I on the whole link space, the electric energy of
    # a chain is additive by construction: R links carry R flux quanta,
    # E(R) = R C2(45), and the static potential is linear with tension
    # sigma = C2(45).
    n_inv_link = int(np.sum(np.abs(w) < 1e-8))
    # the vertex Gauss Laplacian on (pair-adjoint subspace) (x) link: its
    # kernel is the V98/V99 intertwiner space (dimension m_45) and its gap is
    # one flux quantum.
    Ts = [XP.asarray(T) for T in Ts_np]
    Js = [XP.kron(Ts[A], I32) + XP.kron(I32, Ts[A]) for A in range(45)]
    C2pair = XP.zeros((1024, 1024), dtype=np.complex128)
    for J in Js:
        C2pair += J @ J
    evp, evpvec = XP.linalg.eigh(C2pair)
    adj_sel = XP.abs(evp - 64.0) < 1e-6
    n_adj_pair = int(np.sum(to_np(adj_sel)))
    m45 = n_adj_pair // 45
    E45 = evpvec[:, adj_sel]                      # 1024 x 90
    Jsub = [E45.conj().T @ J @ E45 for J in Js]   # 90 x 90
    I90 = XP.eye(90, dtype=np.complex128)
    M3 = XP.zeros((90 * 45, 90 * 45), dtype=np.complex128)
    for A in range(45):
        QA = XP.kron(Jsub[A], I45) + XP.kron(I90, Ls[A])
        M3 += QA.conj().T @ QA
    w3 = np.real(to_np(XP.linalg.eigvalsh(M3)))
    n_ker = int(np.sum(w3 < 1e-8))
    nz3 = w3[w3 >= 1e-8]
    gap3 = float(nz3[0]) if len(nz3) else float("nan")
    out["data"]["link"]["n_invariants"] = n_inv_link
    out["data"]["vertex"] = {"dim": int(90 * 45), "ker": n_ker, "m45": m45,
                             "gap": gap3, "C2_45": c2_45}
    checks.append(("T2 a link cannot be flux-free: the invariant subspace of "
                   "the adjoint link is %d-dimensional (ker of sum_A L_A^dag "
                   "L_A = C2(45) I), so every link of a gauge-invariant chain "
                   "carries one adjoint flux quantum and the electric energy "
                   "of R links is exactly R x %.4f - the static potential is "
                   "LINEAR with tension sigma = C2(45). The vertex Gauss "
                   "Laplacian on (pair-adjoint) (x) link has kernel %d = m_45 "
                   "= %d (the V98/V99 intertwiner space; the glued two-site "
                   "sector is m_45^2 = %d) and gap %.4f = C2(45)"
                   % (n_inv_link, c2_45, n_ker, m45, m45 * m45, gap3),
                   n_inv_link == 0 and n_ker == m45
                   and abs(gap3 - c2_45) < 1e-8))

    # ---- 3. the tension and its physical scale ---------------------------
    sigma_lat = c2_45                       # electric energy per link
    MP = 1.2209e19                          # GeV
    sig_per_len = sigma_lat * MP            # GeV per Planck length
    sig_mass2 = sigma_lat * MP ** 2         # GeV^2
    sig_qcd = 0.42 ** 2
    out["data"]["scale"] = {"sigma_lattice": sigma_lat,
                            "sigma_GeV_per_MPl": sig_per_len,
                            "sigma_GeV2": sig_mass2,
                            "sqrt_sigma_GeV": sig_mass2 ** 0.5,
                            "sigma_qcd_GeV2": sig_qcd,
                            "ratio": sig_mass2 / sig_qcd}
    checks.append(("T3 the tension is sigma = C2(45) = %.4f lattice units of "
                   "electric energy per Planck link = %.4e GeV per Planck "
                   "length = %.3e GeV^2 (sqrt(sigma) = %.3e GeV), against the "
                   "QCD string tension %.3f GeV^2: the layer confines %.2e "
                   "times too hard to be the QCD string"
                   % (sigma_lat, sig_per_len, sig_mass2, sig_mass2 ** 0.5,
                      sig_qcd, sig_mass2 / sig_qcd),
                   sig_mass2 / sig_qcd > 1e30))

    # ---- 4. the step itself has no field strength -------------------------
    from cl10 import supercell as SC
    qs = (0.7, -1.3, 2.1)
    St = {(ev, j): SC.stage_matrix(ev, j, qs[j]) for (ev, j) in SC.STAGES}
    comm, pairs = 0.0, 0
    keys = list(St)
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if keys[i][1] != keys[j][1]:
                a, b = St[keys[i]], St[keys[j]]
                comm = max(comm, float(np.abs(a @ b - b @ a).max()))
                pairs += 1
    # the plaquette of the step: the commutator of two stages on different
    # axes around an elementary square. An area law needs <W_square> < 1 with
    # the exponent proportional to the enclosed area; here it is the identity.
    pdev, wmin = 0.0, 1.0
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            if keys[i][1] != keys[j][1]:
                a, b = St[keys[i]], St[keys[j]]
                Wp = a @ b @ np.linalg.inv(a) @ np.linalg.inv(b)
                pdev = max(pdev, float(np.abs(Wp - np.eye(8)).max()))
                wmin = min(wmin, float(np.abs(np.trace(Wp)) / 8.0))
    out["data"]["step_loops"] = {"stage_commutator_max": comm,
                                 "stage_pairs": pairs,
                                 "plaquette_dev": pdev,
                                 "plaquette_W_min": wmin}
    checks.append(("T4 the implemented step carries NO field strength: its "
                   "%d stage pairs acting on different axes commute to %.2e, "
                   "so every elementary plaquette of the step is the identity "
                   "to %.2e and <W_square> = %.6f: an area law would need "
                   "<W_square> < 1 with an exponent proportional to the "
                   "enclosed area, so the step cannot confine"
                   % (pairs, comm, pdev, wmin),
                   comm < 1e-12 and pdev < 1e-12 and abs(wmin - 1.0) < 1e-12))

    ok = all(c[1] for c in checks)
    out["pass"] = ok
    out["checks"] = [{"name": n, "pass": bool(p)} for n, p in checks]
    out["status"] = (
        "THE GAUGED LAYER CONFINES, AT THE PLANCK SCALE (I). Its link space is "
        "the irreducible adjoint (C2 = %.4f, spread %.1e) and carries NO "
        "invariant vector (%d link singlets), so no link of a gauge-invariant "
        "chain can be flux-free: every link carries one adjoint flux quantum "
        "and the electric energy of R links is exactly R C2(45) - the static "
        "potential is linear with string tension sigma = C2(45) = %.4f lattice "
        "units = %.3e GeV^2, %.2e times the QCD string tension. The vertex "
        "Gauss Laplacian confirms the scale: kernel %d = m_45 (the V98/V99 "
        "intertwiner space) and gap %.4f = C2(45). The step itself has "
        "commuting stage operators on different axes (max commutator %.1e, "
        "every plaquette the identity, <W_square> = %.6f): no field strength, "
        "no area law, no confinement. Confinement is therefore PRESENT in the "
        "construction - in the gauged bond layer - but its string is a Planck "
        "string: the residual of Section 13 is no longer a missing phase, it "
        "is a measured tension forty orders of magnitude too large, and "
        "Lambda_QCD is not produced."
        % (c2_45, c2_spread, n_inv_link, sigma_lat, sig_mass2,
           sig_mass2 / sig_qcd, n_ker, gap3, comm, wmin))

    res_dir = os.path.join(HERE, "results")
    os.makedirs(res_dir, exist_ok=True)
    with open(os.path.join(res_dir, "v130_string_tension.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    for c in out["checks"]:
        print(("  [PASS] " if c["pass"] else "  [FAIL] ") + c["name"])
    print(f"  C2(45) = {c2_45:.6f} | link singlets {n_inv_link} | vertex ker "
          f"{n_ker} gap {gap3:.4f} | sigma/sigma_QCD {sig_mass2 / sig_qcd:.2e}")
    print(f"  wall {time.time() - t0:.1f} s")
    print("VERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
