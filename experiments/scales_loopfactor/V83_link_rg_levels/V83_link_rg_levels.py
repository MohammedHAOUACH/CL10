#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V83 - the traced adjoint link through the RG levels: does the
division compound toward 1/(16 pi)?

**Why (W8, the named next test of V80).**  V80 measured the level-1
division factor of the gauged bond layer: F(eta_root) = 1.06, F_min =
0.52 over the coupling - a factor 26 above 1/(16 pi) - and named
level-N suppression as the next test. V72 measured the link-free flow:
the Cartan content decays only polynomially (0.14 at p = 5). This
experiment iterates the TRACED-LINK decimation exactly to level 4 and
measures the decay law of the Cartan content.

**The exact recursion (one-particle, open chain, maximally mixed link).**
The bond hop operator factorizes through the T.L joint Casimir: on
32 (x) 45 the operator T.L has three eigenvalues ev_a (V80: joint
Casimirs verified), so k1 = exp(i (eta/2) T.L) = sum_a e^{i (eta/2) ev_a} P_a
with P_a the joint-irrep projectors (exact, from the 1440-D
eigendecomposition). A returning path of 2^p hops crossing new bonds
whose links are traced satisfies the recursion

    Q_1(b,c) = Tr_L[ P_b (G (x) I) P_c (G (x) I) (I (x) rho) ]      (32x32)
    S_{p-1}  = sum_{b,c} e^{i (eta/2)(ev_b + ev_c)} Q_{p-1}(b,c)
    Q_p(b,c) = Tr_L[ P_b (G (x) I) (Q_{p-1}(b,c) (x) I) P_c (G (x) I) (I (x) rho) ]

and the level-p on-site block is

    B_p = sum_{a,d} e^{i (eta/2)(ev_a + ev_d)}
          Tr_L[ P_a (G (x) I) (S_{p-1} (x) I) P_d (G (x) I) (I (x) rho) ].

B_1 is V80's two-hop block (regression: F_1 = alpha(2 B_1)/alpha(2 G^2)
reproduces V80's F exactly). Everything is exact linear algebra on
1440 dimensions; the link is traced in the canonical maximally mixed
state at every bond it crosses.

**What is measured.** alpha(B_p) for p = 1..4 at eta = g10 = 1/5 and at
the V80 minimum eta = 1.44, against the loop law (16 pi)^{-p} and the
geometric law F_1^p. The decay exponent is reported, not assumed.

**Verdict.** The compounding law of the traced adjoint link, measured:
whether deeper levels approach the loop factor or stall, and by how
much they miss it.
"""

import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V80_adjoint_link_division"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import sm_qca  # noqa: E402

D_L = 45


def cartan_norm(M, cartan):
    return float(sum(abs(np.trace(C @ M)) / M.shape[0] for C in cartan))


def projectors(TL):
    """Orthogonal projectors onto the T.L eigenspaces (clustered, exact)."""
    w, V = np.linalg.eigh(TL)
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    groups = []
    cur = [0]
    for i in range(1, len(w)):
        if w[i] - w[cur[-1]] < 1e-8:
            cur.append(i)
        else:
            groups.append(cur)
            cur = [i]
    groups.append(cur)
    P, ev = [], []
    for g in groups:
        Vk = V[:, g]
        P.append(Vk @ Vk.conj().T)
        ev.append(float(np.mean(w[g])))
    return P, ev

def q1(Pb, Pc, G, dl):
    """Tr_L[P_b (G x I) P_c (G x I) (I x rho)] as 32x32."""
    d = G.shape[0]
    Gi = np.kron(G, np.eye(dl, dtype=np.complex128))
    M = Pb @ Gi @ Pc @ Gi
    r = M.reshape(d, dl, d, dl)
    return np.einsum("ajbj->ab", r) / dl


def qnext(Pb, Pc, Qprev, G, dl):
    """Tr_L[P_b (G x I) (Qprev x I) P_c (G x I) (I x rho)]."""
    d = G.shape[0]
    Gi = np.kron(G, np.eye(dl, dtype=np.complex128))
    Qi = np.kron(Qprev, np.eye(dl, dtype=np.complex128))
    M = Pb @ Gi @ Qi @ Pc @ Gi
    r = M.reshape(d, dl, d, dl)
    return np.einsum("ajbj->ab", r) / dl


def level_blocks(G, P, ev, eta, pmax):
    ph = np.array([np.exp(1j * (eta / 2) * e) for e in ev])
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))
    # level 1: B_1 = sum ph_b ph_c Q_1(b,c) = Tr_L[k1 (G x I) k1 (G x I) rho]
    Q = [[q1(P[b], P[c], G, D_L) for c in range(3)] for b in range(3)]
    B = sum(ph[b] * ph[c] * Q[b][c] for b in range(3) for c in range(3))
    blocks = [B]
    for _ in range(2, pmax + 1):
        S = sum(ph[b] * ph[c] * Q[b][c] for b in range(3) for c in range(3))
        B = np.zeros((32, 32), dtype=np.complex128)
        for a in range(3):
            for d in range(3):
                Si = np.kron(S, np.eye(D_L, dtype=np.complex128))
                M = P[a] @ Gi @ Si @ P[d] @ Gi
                r = M.reshape(32, D_L, 32, D_L)
                B += ph[a] * ph[d] * np.einsum("ajbj->ab", r) / D_L
        blocks.append(B)
        Q = [[qnext(P[b], P[c], Q[b][c], G, D_L) for c in range(3)]
             for b in range(3)]
    return blocks


def main():
    t0 = time.time()
    checks = {}

    G = E.asnumpy(sm_qca.gut_qca_step(sm_qca.plane_generators(),
                                      sm_qca.G10)).astype(np.complex128)
    h = [E.asnumpy(m) for m in clifford.generators()]
    cartan = [1j * (h[a] @ h[a + 1]) for a in (0, 2, 4, 6, 8)]
    Ts = root_generators()
    f, _k, _as = structure_constants(Ts)
    Ls = adjoint_link_operators(f)
    TL = sum(np.kron(T, L) for T, L in zip(Ts, Ls))
    P, ev = projectors(TL)
    checks["projectors"] = {
        "n_eigenvalues": int(len(ev)),
        "sum_P_err": float(np.abs(sum(P) - np.eye(1440)).max()),
        "idempotence_err": float(max(np.abs(Pa @ Pa - Pa).max()
                                     for Pa in P))}

    ref1 = cartan_norm(2 * G @ G, cartan)
    loop = 1 / (16 * math.pi)

    for eta, tag in ((0.2, "root"), (1.4423, "Fmin")):
        blocks = level_blocks(G, P, ev, eta, 4)
        al = [cartan_norm(2 * B, cartan) for B in blocks]
        F = [a / ref1 for a in al]
        checks[f"levels_{tag}"] = {
            "eta": eta,
            "alpha_Bp": [round(a, 6) for a in al],
            "F_p": [round(fv, 6) for fv in F],
            "F_1_exact": F[0],
            "loop_law": [round(loop ** p, 8) for p in range(1, 5)],
            "geometric_F1_pow_p": [round(F[0] ** p, 6) for p in range(1, 5)],
            "decay_exponent_p1_to_p4": round(
                float(np.log(al[3] / al[0]) / 3), 4)}

    # regression: level-1 F at eta = 0.2 must equal V80's F(g10)
    F_root = checks["levels_root"]["F_1_exact"]
    reg_err = abs(F_root - 1.0596866584601925)
    checks["regression_v80"] = {"F_1_root": F_root, "err": reg_err}

    ok = (checks["projectors"]["sum_P_err"] < 1e-9
          and checks["projectors"]["idempotence_err"] < 1e-9
          and reg_err < 1e-9)

    lv = checks["levels_root"]
    Fp = lv["F_p"]
    ratio4 = Fp[3] / loop ** 4
    concl = ("COMPOUNDS TOWARD the loop law" if Fp[3] < 2 * loop ** 4 else
             "DOES NOT COMPOUND TO the loop law")
    status = (
        f"THE TRACED-LINK DIVISION {concl}: at eta = 1/5 the level-p "
        f"Cartan-content ratios are F_p = {Fp[0]:.3f}, {Fp[1]:.3f}, "
        f"{Fp[2]:.3f}, {Fp[3]:.3f} for p = 1..4 against the loop law "
        f"(16 pi)^-p = {loop:.3f}, {loop ** 2:.2e}, {loop ** 3:.2e}, "
        f"{loop ** 4:.2e}: the deepest level sits a factor {ratio4:.1e} "
        f"away (decay exponent {lv['decay_exponent_p1_to_p4']:.2f} per "
        f"level); at the V80 minimum eta = 1.4423 the level-4 ratio is "
        f"{checks['levels_Fmin']['F_p'][3]:.3f}. The compounding law of "
        "the traced adjoint link is measured, not assumed; whether the "
        "automaton reaches the loop factor at deeper levels or at a "
        "different coupling is the named residual (W8).")

    verdict = {"experiment": "V83_link_rg_levels",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v83_link_rg_levels.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V83 - traced adjoint link through the RG levels")
    print("=" * 68)
    print(f"  projectors: n = {checks['projectors']['n_eigenvalues']}, "
          f"sum err {checks['projectors']['sum_P_err']:.2e}, "
          f"P^2 err {checks['projectors']['idempotence_err']:.2e}")
    for tag in ("root", "Fmin"):
        lv = checks[f"levels_{tag}"]
        print(f"  eta = {lv['eta']:<7} F_p = " +
              ", ".join(f"{fv:.4f}" for fv in lv["F_p"]) +
              f"   (exponent {lv['decay_exponent_p1_to_p4']:.2f})")
    print(f"  regression vs V80 F(g10)  = {reg_err:.2e}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v83_link_rg_levels.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
