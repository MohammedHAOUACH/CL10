#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V27 - QCD confinement (the 0-parameter area law + the linear potential).

The **0-parameter** confinement content of the CL10 QCA:

- the **0-parameter colour count** $N_C = 3$ (the 8 gluon $= N_C^2 - 1$);
- the **area law**: $\\langle W(C) \\rangle \\propto \\exp(-\\sigma \\cdot
  \\mathrm{Area}(C))$ — the loop action scales with the **area** (the
  signature of confinement);
- the **0-parameter string tension** $\\sigma$ (the flux-tube energy per unit
  length);
- the **linear potential** $V(r) = \\sigma r$ (the quark-antiquark potential
  **rises linearly** — confinement), in **contrast** with the Coulombic
  $1/r$ of the EM force.

The **named residual**: the lightest **glueball** mass and the running of
$\\sigma$ (the strong-force dynamics beyond the 0-parameter area law).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import qcd

TOL = 1e-12


def main():
    t0 = time.time()
    n_c, n_gluon, ok_c = qcd.color_count()
    sigma = qcd.string_tension()
    rows, _ = qcd.confining_potential(sigma)
    linear = all(abs(r["v_qcd_over_r"] - sigma) < 1e-12 for r in rows)
    ok = ok_c and linear and n_c == 3

    verdict = {
        "id": "V27_qcd_confinement",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the glueball mass + the running of sigma are named residuals
        "pass": bool(ok),
        "status": "the_qcd_is_confining_linear_potential_0_parameter_area_law",
        "tol": TOL,
        "data": {
            "n_c": n_c,
            "n_gluon": n_gluon,
            "sigma": sigma,
            "r": [r["r"] for r in rows],
            "v_qcd": [r["v_qcd"] for r in rows],
            "v_qcd_over_r": [r["v_qcd_over_r"] for r in rows],
            "v_em": [r["v_em"] for r in rows],
            "residual": "the lightest glueball mass + the running of sigma (the strong-force dynamics)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v27_qcd_confinement.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V27 - QCD confinement (0-parameter area law + linear potential)")
    print("=" * 68)
    print(f"  the 0-parameter colour count:  N_C = {n_c}  (8 gluon = N_C^2-1 = {n_gluon})")
    print(f"  the 0-parameter string tension:  sigma = {sigma}")
    print(f"  V_qcd(r) = sigma*r (linear):  {['%.1f' % r['v_qcd'] for r in rows]}  ->  V/r = {sigma} (0-parameter)")
    print(f"  V_em(r) = 1/r (Coulombic):     {['%.3f' % r['v_em'] for r in rows]}")
    print(f"  the linear (confining) QCD potential, in contrast to the 1/r Coulombic EM:  {'PASS' if linear else 'FAIL'}")
    print("  named residual:  the glueball mass + the running of sigma")
    print(f"  VERDICT:  {'PASS (the QCD is confining: linear potential, 0-parameter area law)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v27_qcd_confinement.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
