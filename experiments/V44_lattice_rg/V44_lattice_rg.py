#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V44 — the QCA lattice RG: the S-matrix to the cross-sections (W9, C6).

**The 't Hooft program, completed at the cross-section level.**  A unitary
QCA is a **UV-finite** formulation of a QFT.  The GUT QCA step ``U`` (the
32x32 product of the 45 Spin(10) plane rotations at ``g_10 = 1/5``) is the
time-evolution for one Planck step, so its matrix elements are the
**S-matrix** (V41).  This experiment computes the **cross-section** from
that S-matrix, with **no fit**:

1. **UV-finiteness (the lattice RG cutoff).**  The Planck lattice spacing
   ``a = l_P`` bounds the momentum to the Brillouin zone, so **every**
   cross-section is **bounded** by the unitarity limit
   ``sigma <= sigma_max = 4 pi / M_P^2``.  The QCA has a **finite** UV
   cutoff (the Planck scale) and is **UV-finite** (no UV divergences).

2. **The 0-parameter 2->2 scattering.**  The 2-particle S-matrix
   ``S^(2) = U (x) U`` (the 1024x1024 tensor product, exactly unitary)
   gives the **inelastic (scattering) probability per Planck step**
   ``P2 = 1 - dbar^2`` (``dbar = (1/32) sum_i |U_ii|^2``): **nonzero**
   (the 45 planes mix the 32, so the QCA scatters) and **0-parameter**.

3. **The 0-parameter GUT-scale cross-section.**
   ``sigma_0 = P2 / M_P^2`` — nonzero, **bounded** (``sigma_0 < sigma_max``),
   and **0-parameter** (the only inputs are ``g_10 = 1/5`` and ``M_P``).

4. **The non-perturbative enhancement (0-parameter).**  The QCA
   cross-section vs the **tree-level** perturbative GUT cross-section
   ``sigma_pert ~ g_10^4 / (4 pi M_P^2)`` is the 0-parameter factor
   ``R = sigma_0 / sigma_pert = 4 pi P2 / g_10^4``: the QCA is the
   **all-orders, UV-finite** result, the perturbative QFT the leading
   order in the small coupling.

**The named residual (W9).**  The **lattice RG** (running the cross-section
from the Planck scale **down** to the physical scale — the non-trivial /
anomalous part) and the **matching** of these 0-parameter amplitudes to the
**measured** cross-sections — the open part of the 't Hooft program.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import lattice_rg

TOL = 1e-12


def main():
    t0 = time.time()
    r = lattice_rg.lattice_rg()
    cs = r["cross_section"]
    ds = r["dispersion"]
    uv = r["uv_bound"]
    ok = (r["s2_unitary_err"] < TOL
          and r["closed_form_match"]
          and cs["P2"] > 0.0
          and r["uv_finite"]
          and cs["R_qca_over_pert"] > 0.0
          and r["optical_theorem"]
          and ds["uv_bounded"]
          and abs(ds["zero_point"]) < 1e-8
          and r["only_g10"])

    verdict = {
        "id": "V44_lattice_rg",
        "tier": 2,
        "property": "P2",
        "residual": True,          # the RG matching to measured cross-sections (W9) is the named residual
        "pass": bool(ok),
        "status": "the_0_parameter_GUT_cross_section_is_UV_finite_nonzero_bounded_the_RG_matching_is_the_named_residual",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "g10": r["g10"],
            "M_P_GeV": r["M_P_GeV"],
            "P1": cs["P1"],
            "P2": cs["P2"],
            "s2_unitary_err": r["s2_unitary_err"],
            "closed_form_match": r["closed_form_match"],
            "sigma_0_GeV2": cs["sigma_0_GeV2"],
            "sigma_0_cm2": cs["sigma_0_cm2"],
            "sigma_pert_GeV2": cs["sigma_pert_GeV2"],
            "sigma_pert_cm2": cs["sigma_pert_cm2"],
            "R_qca_over_pert": cs["R_qca_over_pert"],
            "sigma_max_GeV2": uv["sigma_max_GeV2"],
            "sigma_max_cm2": uv["sigma_max_cm2"],
            "uv_finite": r["uv_finite"],
            "dispersion": {
                "n_levels": ds["n_levels"],
                "omega_max_over_pi": ds["omega_max_over_pi"],
                "uv_bounded": ds["uv_bounded"],
                "zero_point": ds["zero_point"],
            },
            "optical_theorem": r["optical_theorem"],
            "elastic": cs["elastic"],
            "sigma_over_smax": cs["sigma_0_GeV2"] / uv["sigma_max_GeV2"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v44_lattice_rg.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V44 - the QCA lattice RG: the S-matrix to the cross-sections (W9)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  the ONLY inputs:  g_10 = {r['g10']:.4f} (= 1/5)   M_P = {r['M_P_GeV']:.3e} GeV")
    print(f"  [1] UV-finiteness:  sigma_max = 4 pi / M_P^2 = {uv['sigma_max_cm2']:.3e} cm^2  (finite, no UV divergence)")
    print(f"  [2] 2->2 inelastic per step:  P2 = {cs['P2']:.4f}  (1p: P1 = {cs['P1']:.4f});  "
          f"S^(2) = U (x) U  (1024x1024, unitary err {r['s2_unitary_err']:.2e})")
    print(f"  [3] 0-parameter cross-section:  sigma_0 = P2 / M_P^2 = {cs['sigma_0_cm2']:.3e} cm^2")
    print(f"  [4] unitarity consistency:  sigma_0 / sigma_max = {cs['sigma_0_GeV2'] / uv['sigma_max_GeV2']:.3e}  "
          f"({'< 1: bounded' if r['uv_finite'] else 'VIOLATED'})")
    print(f"  [5] non-perturbative:  R = sigma_0 / sigma_pert = {cs['R_qca_over_pert']:.3e}  (the QCA resums all orders)")
    print(f"  [6] the QCA dispersion:  {ds['n_levels']} levels, max |omega| = {ds['omega_max_over_pi']:.3f} * pi  "
          f"({'UV-bounded' if ds['uv_bounded'] else 'NOT'});  optical theorem: the inelastic P2 IS the cross-section  "
          f"({'holds' if r['optical_theorem'] else 'fails'})")
    print(f"  named residual (W9):  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the 0-parameter GUT cross-section is UV-finite, nonzero, bounded; the RG matching is the named residual)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v44_lattice_rg.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
