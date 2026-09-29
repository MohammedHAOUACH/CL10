#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V37 — the Einstein-Hilbert action **derived from** the QCA metric
dynamics (C2).

**The 0-parameter result.**  V31 *evaluated* the Riemann curvature on
prescribed metrics (the direction **action -> metric**).  This program
closes the **other** direction — **dynamics <- action**: the
**Euler-Lagrange equation** of the **0-parameter** discrete quadratic
Einstein-Hilbert action on the Planck lattice

    S = (1/16*pi) * (1/2) * sum_x sum_t u_t [ u_{t+1} - 2 u_t + u_{t-1}
    - nabla^2 u_t ]        (G = 1, the 6-neighbour discrete Laplacian)

is

    delta S / delta u = (1/16*pi) * [ u_{t+1} - 2 u_t + u_{t-1}
    - nabla^2 u_t ] = 0,

i.e. **exactly** the **QCA metric dynamics** (the discrete free wave,
``c = 1``, the massless spin-2).  The **0-parameter** coefficient
``1/(16*pi)`` is the **unique** normalization (``G = 1``) of the action.

**The checks** (the QCA dynamics is integrated **on the GPU**):
1. **dynamics <- action**: the action is **stationary** on every QCA
   trajectory — the **numerical** functional derivative (the central
   finite difference of the action) vanishes at the probe points;
2. **the action selects the dynamics**: the ELM of a field **off** the
   QCA dynamics (a perturbed final step) is **nonzero**;
3. **the 0-parameter coefficient**: ``1/(16*pi)`` (``G = 1``);
4. **the continuum limit**: the ELM operator gives the **discrete
   dispersion** ``cos(w) = 1 + lambda(k)/2`` — the **massless** dispersion
   ``w/|k| -> 1`` (``c = 1``) **derived from the action**.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import einstein

TOL = 1e-6


def main():
    t0 = time.time()
    r = einstein.einstein_hilbert_from_dynamics(tol=TOL)

    ok = (r["ok"]
          and r["stationary_ok"]
          and r["selects_ok"]
          and r["coeff_ok"]
          and r["continuum_ok"])

    verdict = {
        "id": "V37_eh_action",
        "tier": 2,
        "property": "P2,P4",
        "residual": False,         # the derivation is a 0-parameter PASS (no new residual)
        "pass": bool(ok),
        "status": "the_eh_action_derives_the_qca_metric_dynamics_the_direction_dynamics_from_action",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "L": r["L"],
            "T": r["T"],
            "n_probe": r["n_probe"],
            "G": r["G"],
            "coupler": r["coupler"],
            "elm_on_qca_trajectory": r["elm_on_qca_trajectory"],
            "stationary_ok": r["stationary_ok"],
            "elm_off_dynamics": r["elm_off_dynamics"],
            "selects_ok": r["selects_ok"],
            "coeff_ok": r["coeff_ok"],
            "continuum_ratio": r["continuum_ratio"],
            "continuum_ok": r["continuum_ok"],
            "residual": "none new (the non-linear completion beyond the quadratic order remains the Choptuik model, as before)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v37_eh_action.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V37 - the Einstein-Hilbert action DERIVED from the QCA dynamics (C2)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print("  (1) dynamics <- action:  the ELM of the action on the QCA trajectory")
    print("        (the finite-difference functional derivative):  "
          f"{r['elm_on_qca_trajectory']:.2e}  (~0, stationary)")
    print(f"  (2) the action SELECTS the dynamics:  the ELM of an off-dynamics field = {r['elm_off_dynamics']:.3f}  (nonzero)")
    print(f"  (3) the 0-parameter coefficient:  1/(16*pi) = {r['coupler']:.6f}  (G = {r['G']})")
    print(f"  (4) the continuum dispersion (from the action):  w/|k| = {['%.4f' % x for x in r['continuum_ratio']]}  (-> 1, massless spin-2)")
    print(f"  VERDICT:  {'PASS (the QCA metric dynamics is the EOM of the 0-parameter EH action — dynamics derived from the action)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v37_eh_action.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
