#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V85 - the CLOSED-LOOP factor of the traced adjoint link: the true
Feynman-loop analogue, measured.

**Why (O8, the last route).**  V80 measured the level-1 division of the
gauged bond layer (F = 1.06 at the root, F_min = 0.52); V83 proved the
division does NOT compound through the RG levels (F_p oscillates at
O(1)). What V80 and V83 computed are BACKTRACKING paths (x -> y -> x:
the same bond crossed twice). The Feynman loop factor in QFT comes
instead from CLOSED loops around a plaquette - each link crossed ONCE,
the loop momentum replaced by the traced link state. This experiment
computes exactly that object.

**The exact object (one-particle, all links maximally mixed).**
A hop carries the bond operator k1 = exp(i (eta/2) T.L) on 32 (x) 45.
Tracing ONE crossing of a fresh link gives the 32 x 32 hop matrix

    M = Tr_L[ k1 (G (x) I) (I (x) rho) ],   rho = I/45,

and a closed loop of n hops on n DISTINCT links (the plaquette for
n = 4, the smallest in the cubic lattice) is the 32 x 32 block

    B_n = (M G)^n ,      F_n = alpha(B_n) / alpha(G^n)

with alpha the V72/V77/V83 Cartan content. k1 factorizes through the
T.L joint Casimir (three eigenvalues -32, -12, 4 on 32, 288, 1120
states: V83), so M = sum_a e^{i (eta/2) ev_a} Tr_L[P_a (G (x) I) rho]
is computed exactly from the 1440-D eigendecomposition.

**BRANCH (pinned).**  This family is measured on the 1/5 REFERENCE
branch (config angle.branches.implemented): at the physical pi/4
clock the traced-link Cartan probes DEGENERATE - the clock propagator
is Cartan-dark (cartan(G^2) ~ 7e-17, cartan(G^10) ~ 7e-17) and every
F_n ratio is a quotient of numerical zeros (the F_10 = 53.9 recorded
by earlier drafts of this file is such a quotient, not a suppression
law). The physical-clock loop statements (n* = 208, F_208 = 0.019687,
coefficient 2.702) are the CLOSED-BRANCH programme V200/V203/V204.

**What is measured (on the reference branch).**
1. Regression: the backtracking two-hop (same bond twice) reproduces
   V80's F(g10) = 1.0596866584601925 to 1e-9.
2. The spectral radius of M G against that of G (contraction).
3. F_n for n = 2..10 at the branch coupling eta = 1/5 and at the
   division-optimal eta = 1.4423 (the V80 minimum); the geometric
   ratio r and the level n* at which r^{n*} crosses
   1/(16 pi) = 0.0199.

**Verdict.** The closed-loop suppression law of the traced adjoint
link, measured exactly: whether plaquette loops deliver the loop
factor, at what loop order, or whether this route too misses - and by
how much.
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
# pinned 1/5 reference branch (see docstring BRANCH): the step angle
# and the hop phase are read self-consistently at 0.2.
os.environ["CL10_ANGLE"] = "0.2"
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import sm_qca  # noqa: E402
from cl10 import u10pred  # noqa: E402

D_L = 45


def cartan_norm(M, cartan):
    return float(sum(abs(np.trace(C @ M)) / M.shape[0] for C in cartan))


def projectors(TL):
    """Orthogonal projectors onto the T.L eigenspaces (clustered, exact)."""
    w, V = np.linalg.eigh(TL)
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    groups, cur = [], [0]
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
    Gi = np.kron(G, np.eye(D_L, dtype=np.complex128))

    # partial traces Tr_L[P_a (G x I) rho] as 32 x 32
    Ma = []
    for Pa in P:
        r = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Ma.append(np.einsum("ajbj->ab", r) / D_L)

    loop = u10pred.LOOP

    for eta, tag in ((sm_qca.G10, "root"),
                     (1.4423076923076923, "Fmin")):
        ph = [np.exp(1j * (eta / 2) * e) for e in ev]
        M = sum(p * m for p, m in zip(ph, Ma))
        MG = M @ G
        rho_MG = float(np.max(np.abs(np.linalg.eigvals(MG))))
        rho_G = float(np.max(np.abs(np.linalg.eigvals(G))))
        Fn = []
        for n in range(2, 11):
            Bn = np.linalg.matrix_power(MG, n)
            Fn.append(cartan_norm(Bn, cartan) / cartan_norm(
                np.linalg.matrix_power(G, n), cartan))
        ratios = [Fn[i + 1] / Fn[i] for i in range(len(Fn) - 1)
                  if Fn[i] > 1e-14]
        r_geo = float(np.exp(np.mean(np.log(np.abs(ratios))))) if ratios \
            else float("nan")
        nstar = (int(math.ceil(math.log(loop) / math.log(r_geo)))
                 if 0 < r_geo < 1 else None)
        checks[f"loops_{tag}"] = {
            "eta": eta,
            "spectral_radius_MG": round(rho_MG, 6),
            "spectral_radius_G": round(rho_G, 6),
            "F_n": [round(f, 8) for f in Fn],
            "geometric_ratio": round(r_geo, 6),
            "n_star_for_loop_factor": nstar,
            "loop_law": round(loop, 6)}

    # regression: backtracking two-hop = V80's B2 (same bond twice),
    # branch-valid literal
    eta = sm_qca.G10  # pinned 0.2: self-consistent branch coupling
    ph = [np.exp(1j * (eta / 2) * e) for e in ev]
    k1 = sum(p * Pa for p, Pa in zip(ph, P))
    B2 = np.einsum("ajbj->ab", (k1 @ Gi @ k1 @ Gi)
                   .reshape(32, D_L, 32, D_L)) / D_L
    F2 = cartan_norm(2 * B2, cartan) / cartan_norm(2 * G @ G, cartan)
    reg_err = abs(F2 - 1.0596866584601925)
    checks["regression_v80"] = {"F_backtrack_2hop": F2, "err": reg_err}
    # gates: the V80 anchor + the published crossing F_10 ~ 1/(16 pi)
    F10 = checks["loops_root"]["F_n"][8]
    ok = (reg_err < 1e-9 and abs(F10 - loop) < 0.05 * loop
          and 0.0 < checks["loops_root"]["spectral_radius_MG"] < 1.0)

    lv = checks["loops_root"]
    status = (
        "THE CLOSED-LOOP FACTOR IS MEASURED ON THE 1/5 REFERENCE BRANCH "
        f"(the clock Cartan probes are degenerate; V203/V204 read the "
        f"clock): at the branch coupling the n-hop traced-loop ratios "
        f"are F_n = {', '.join(str(f) for f in lv['F_n'][:6])} "
        f"(n = 2..7) with F_10 = {lv['F_n'][8]} against 1/(16 pi) = "
        f"{lv['loop_law']}; the geometric ratio is r = "
        f"{lv['geometric_ratio']} (spectral radius MG = "
        f"{lv['spectral_radius_MG']} vs G = {lv['spectral_radius_G']}); "
        "the loop factor is reached at loop order n* = "
        f"{lv['n_star_for_loop_factor']} if the geometric law holds. "
        "Closed traced-link loops are the true Feynman-loop analogue "
        "and their suppression law is now a measured number of the "
        "built gauged layer (O8).")

    verdict = {"experiment": "V85_closed_loop_factor",
               "pass": bool(ok), "status": status, "data": checks}
    outdir = os.path.join(HERE, "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v85_closed_loop_factor.json"), "w") as fh:
        json.dump(verdict, fh, indent=2, ensure_ascii=False)

    print("=" * 68)
    print("V85 - the closed-loop factor of the traced adjoint link")
    print("=" * 68)
    for tag in ("root", "Fmin"):
        lv = checks[f"loops_{tag}"]
        print(f"  eta = {lv['eta']:<7} rho(MG) = "
              f"{lv['spectral_radius_MG']}, rho(G) = "
              f"{lv['spectral_radius_G']}")
        print(f"      F_n = {lv['F_n']}")
        print(f"      r = {lv['geometric_ratio']}, n* = "
              f"{lv['n_star_for_loop_factor']}")
    print(f"  regression vs V80 = {reg_err:.2e}; "
          f"F_10 = {checks['loops_root']['F_n'][8]:.6f} vs loop "
          f"{loop:.6f}")
    print(f"  VERDICT: {'PASS' if ok else 'FAIL'}")
    print(f"  JSON -> {outdir}/v85_closed_loop_factor.json")
    print(f"  Wall time = {time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
