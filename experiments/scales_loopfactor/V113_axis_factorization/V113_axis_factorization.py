#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""V113 - the axis-family factorization behind V108's coefficient 3.

**Why (V108).** V108 reads the hierarchy coefficient as
v/M_P = 3 F_10^10: the measured closed-loop factor (V85) to the ladder
power, TIMES THREE - the three axis families of doubler nodes (V70).
The claim is CONDITIONAL on an "axis-family factorization of the
closed-loop sum". This experiment tests what that factorization can
and cannot be.

**The test.**
1. Regression: reproduce F_10 = 0.01918639 at the root coupling with
   the V85 machinery (joint-Casimir projectors, traced link).
2. Node independence: the loop factor is a property of the LINK (the
   traced adjoint hop), not of the node momentum. Verify the V70
   structure it would need: the 12 doubler nodes form 3 axis families
   of 4 (zero-axis labeling), the Floquet spectra are identical within
   each family AND across families - so every node channel carries the
   SAME loop factor F_10.
3. The sum: if the closed-loop sum runs over ALL 12 node channels it
   gives 12 F_10^10 (coefficient 8.35, 4.02x the observed 2.077); if
   it runs over ONE node per axis family - the V68/V96 split set, the
   selective reading - it gives exactly 3 F_10^10 (coefficient 2.088,
   99.44%). Decide which reading the factorization supports.

**Verdict.** The factorization is a NODE-SET statement, not a symmetry
theorem: equal per-node contributions are confirmed (family symmetry
of the Floquet spectra + link-local loop factor), but the coefficient
3 is exactly the choice of one split node per family - V108's
conditionality IS O9's split-node-set input, made precise.
"""
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
sys.path.insert(0, ROOT)
# pinned 1/5 reference branch: F10_REF = 0.01918639 and the
# coefficient-3 reading are reference-branch statements (the clock
# reading is V203/V209: F_208 = 0.019685 float64, 0.0196928 certified at
# 50 digits; coefficient 2.699-2.700 float64, 2.710 certified). At the pi/4
# clock the traced-link Cartan probes degenerate.
os.environ["CL10_ANGLE"] = "0.2"
from cl10.paths import exp
sys.path.insert(0, exp("V79_gauged_bond_layer"))
sys.path.insert(0, exp("V85_closed_loop_factor"))

from V79_gauged_bond_layer import (adjoint_link_operators,  # noqa: E402
                                   root_generators, structure_constants)
from V85_closed_loop_factor import cartan_norm, projectors  # noqa: E402
from cl10 import clifford  # noqa: E402
from cl10 import engine as E  # noqa: E402
from cl10 import sm_qca  # noqa: E402
from cl10 import supercell as SC  # noqa: E402
from cl10 import config as _config  # noqa: E402

D_L = 45
F10_REF = 0.01918639          # V85 / V108
COEF_OBS = float(_config.CFG.physics.hierarchy_coefficient_observed)
                               # the one observed coefficient: config derives it
                               # as (v / M_P) / loop^10 from the measured constants


def main():
    # ---- 1. regression: F_10 at the root coupling (V85 machinery) ---------
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
    Ma = []
    for Pa in P:
        r = (Pa @ Gi).reshape(32, D_L, 32, D_L)
        Ma.append(np.einsum("ajbj->ab", r) / D_L)

    eta = sm_qca.G10
    ph = [np.exp(1j * (eta / 2) * e) for e in ev]
    M = sum(p * m for p, m in zip(ph, Ma))
    MG = M @ G
    F10 = cartan_norm(np.linalg.matrix_power(MG, 10), cartan) / \
        cartan_norm(np.linalg.matrix_power(G, 10), cartan)
    reg_err = abs(F10 - F10_REF)

    # ---- 2. the 12 nodes: axis families and Floquet equality --------------
    nodes = [tuple(float(v) for v in q) for q in SC.M_POINTS]
    fam = {"x": [], "y": [], "z": []}
    for q in nodes:
        axis = next(j for j in range(3) if abs(q[j]) < 1e-9)
        fam[["x", "y", "z"][axis]].append(q)
    _B, _U32, U16 = SC.sectors()

    def spectrum(q):
        ph = np.sort(np.angle(SC.phases(SC.bloch(np.array(q), U16))))
        return ph

    within = 0.0
    for k, qs in fam.items():
        s0 = spectrum(qs[0])
        for q in qs[1:]:
            within = max(within, float(np.abs(spectrum(q) - s0).max()))
    across = 0.0
    s_x = spectrum(fam["x"][0])
    for k in ("y", "z"):
        across = max(across, float(np.abs(spectrum(fam[k][0]) - s_x).max()))

    # ---- 3. the two readings of the closed-loop sum ------------------------
    # coefficient = sum_factor * (16 pi F_10)^10
    base = (16 * math.pi * F10) ** 10
    coef_12 = 12.0 * base
    coef_3 = 3.0 * base
    coef_1 = 1.0 * base

    ok = reg_err < 1e-6 and within < 1e-9 and across < 1e-9 \
        and len(nodes) == 12 and all(len(v) == 4 for v in fam.values())
    status = (
        f"AXIS-FAMILY FACTORIZATION AUDIT: F_10 reproduced at "
        f"{F10:.8f} (V85 regression {reg_err:.2e}). The 12 doubler nodes "
        f"form 3 axis families of 4; the Floquet spectra are identical "
        f"WITHIN families to {within:.2e} and ACROSS families to "
        f"{across:.2e} - every node channel carries the same link-local "
        f"loop factor, so per-node contributions are EQUAL. The sum is "
        f"then a NODE-SET choice: all 12 nodes give coefficient "
        f"{coef_12:.4f} ({coef_12 / COEF_OBS:.2f}x the observed "
        f"{COEF_OBS:.4f}); one split node per family (the V68/V96 "
        f"selective reading) gives {coef_3:.4f} "
        f"({100 * (1 - abs(coef_3 - COEF_OBS) / COEF_OBS):.2f}% "
        f"agreement); a single node gives {coef_1:.4f}. VERDICT: the "
        f"factorization is CONFIRMED as family symmetry (equal "
        f"contributions) but the coefficient 3 is EXACTLY the choice of "
        f"one split node per axis family - V108's conditionality is "
        f"O9's split-node-set input, not a factorization theorem.")
    out = {"experiment": "V113_axis_factorization", "pass": bool(ok),
           "status": status,
           "data": {"F_10_measured": round(float(F10), 8),
                    "regression_vs_V85": float(reg_err),
                    "n_nodes": len(nodes),
                    "families": {k: [[round(v, 4) for v in q]
                                     for q in qs]
                                 for k, qs in fam.items()},
                    "spectral_equality_within": float(within),
                    "spectral_equality_across": float(across),
                    "coef_all12": round(coef_12, 4),
                    "coef_one_per_family": round(coef_3, 4),
                    "coef_single": round(coef_1, 4),
                    "coef_observed": COEF_OBS,
                    "agreement_one_per_family":
                        round(1 - abs(coef_3 - COEF_OBS) / COEF_OBS, 5)}}
    residir = os.path.join(HERE, "results")
    os.makedirs(residir, exist_ok=True)
    with open(os.path.join(residir, "v113_axis_factorization.json"),
              "w") as fh:
        json.dump(out, fh, indent=1)
    print(status)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
