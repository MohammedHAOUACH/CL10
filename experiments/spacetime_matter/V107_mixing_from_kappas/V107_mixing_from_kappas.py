#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V107 - can V96's kappa_s pattern touch the W6 mixing deviations?

**The question.** W6: the democratic base gives sin^2(12) = 1/3,
theta13 = 0, CKM = 1; data: 0.308, 0.0222, 0.225. V104 excluded the
natural vev placements. V96 built a flavour-dependent bond layer with
per-class strengths kappa_s. Can THAT pattern supply the W6 splitting?

**The structural answer (computed).** V96's classes are the SIX SM
GAUGE classes (Q, u^c, d^c, L, e^c, nu^c) - each class contains all
THREE generations. A per-class coupling is therefore GENERATION-BLIND:
it rescales each Yukawa block uniformly,
    Y_u -> (1 + eps*kappa_u) Y_u,   Y_d -> (1 + eps*kappa_d) Y_d, ...
and the mixing matrices are invariant under uniform rescaling of a
block (the singular vectors of cY equal those of Y). So the PMNS and
CKM are EXACTLY unchanged for every eps: the deviations stay at zero.
V96's flavour dependence lives in gauge-class space; W6 lives in
generation space. The two are orthogonal: the W4 fix cannot be the W6
source, and W6 needs a generation-space structure (as V104 already
forced: >= 2 vev terms, non-cubic direction).

PASS (as an audit) = the numerical check confirms exact invariance:
max |angle(eps) - angle(0)| < 1e-12 over a full scan of eps, for the
actual V96 kappa pattern, on non-degenerate (perturbed democratic)
Yukawa blocks so the SVD is well posed.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)


# V96 n_nodes = 3 solution, class order = sorted (X,Y) keys of the 16
# weights (the same order V96 uses); the class->SM assignment does not
# matter for the theorem below, and is recorded from V96's fragment map.
V96_KAPPAS = [-4.5311, -5.7515, 9.1857, 7.6558, 9.1221, 5.8487]


def main():
    rng = np.random.default_rng(107)
    # democratic base Yukawas (rank 1, the 0-parameter structure V28/V45)
    # with a generic complex democratic direction per sector, plus a FIXED
    # tiny non-degenerate perturbation so the SVD is well posed (the same
    # matrix enters base and scaled cases; uniform scaling cannot change
    # its singular vectors - that is exactly what this experiment checks)
    PERT = np.diag([1.0, 2.0, 3.0])

    def demo():
        v = rng.normal(size=3) + 1j * rng.normal(size=3)
        v /= np.linalg.norm(v)
        return np.outer(v, v) + 1e-2 * PERT

    Yu, Yd, Ye, YD = demo(), demo(), demo(), demo()

    def angles(Yu, Yd, Ye, YD):
        def svd_mix(Y):
            U, _, _ = np.linalg.svd(Y)
            return U
        Ue = svd_mix(Ye.conj().T @ Ye)
        U, s, Vh = np.linalg.svd(YD)
        Un = U
        P = Un.conj().T @ Ue
        s12 = float(np.abs(P[0, 1]) ** 2)
        s13 = float(np.abs(P[0, 2]) ** 2)
        s23 = float(np.abs(P[1, 2]) ** 2)
        Uc_u = svd_mix(Yu.conj().T @ Yu)
        Uc_d = svd_mix(Yd.conj().T @ Yd)
        C = Uc_d.conj().T @ Uc_u
        l12 = float(np.abs(C[0, 1]) ** 2)
        return (s12, s13, s23, l12)

    base = angles(Yu, Yd, Ye, YD)
    worst = 0.0
    for eps in (0.0, 1e-3, 0.01, 0.1, 1.0, 10.0, -0.5, -0.9):
        # per-GAUGE-class rescaling: kappa index 0..5 = (Q,u,d,L,e,N);
        # u-sector block, d-sector block, charged-lepton block,
        # Majorana block (L, N entries of the type-II/see-saw terms)
        k = np.array(V96_KAPPAS)
        f_u, f_d, f_e, f_D = (1 + eps * k[1], 1 + eps * k[2],
                              1 + eps * k[4], 1 + eps * (k[3] + k[5]) / 2)
        a = angles(f_u * Yu, f_d * Yd, f_e * Ye, f_D * YD)
        worst = max(worst, max(abs(x - y) for x, y in zip(a, base)))

    ok = worst < 1e-12
    status = ("NEGATIVE (structural): V96's kappa_s are per-GAUGE-class "
              "couplings - generation-blind - so they rescale each Yukawa "
              "block uniformly and leave PMNS/CKM EXACTLY invariant "
              f"(max angle shift {worst:.1e} over eps in [0,10]). The W4 "
              "bond-layer flavour dependence lives in gauge-class space, "
              "W6 in generation space: V96 cannot supply the mixing "
              "deviations; W6 needs generation-space structure (V104: "
              ">= 2 vev terms, non-cubic).")
    out = {"experiment": "V107_mixing_from_kappas", "pass": bool(ok),
           "status": status,
           "data": {"v96_kappas": V96_KAPPAS,
                    "base_angles_s12_s13_s23_l12": list(base),
                    "max_angle_shift_over_eps_scan": worst,
                    "eps_scan": [0.0, 1e-3, 0.01, 0.1, 1.0, 10.0,
                                 -0.5, -0.9]}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v107_mixing_from_kappas.json"), "w") as f:
        json.dump(out, f, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
