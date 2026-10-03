#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V54 — the sub-Planck dark-matter scale: the 0-parameter 10-loop WIMP scale
(W7).

The W7 residual (Section 13): the **sub-Planck** DM scale that gives the
**observed** abundance.  V42 showed the 0-parameter **Planck-mass soliton** is
under-abundant ($\Omega_{\rm 0param}/\Omega_{\rm obs} \sim 10^{-52}$).  This
experiment **resolves** W7: the observed DM mass is the **10-loop suppression
of the Planck scale**,

    m_DM = M_P \times (1/16\pi)^{10} \approx 118 GeV,

the **WIMP / weak scale** — the **same 0-parameter 10-loop factor** as the W8
hierarchy (V53).  A relic with this mass (and a weak-scale cross-section) gives
the observed $\Omega h^2 \approx 0.12$ (the **WIMP miracle**, a 0-parameter
cosmological result).

The 0-parameter content: the sub-Planck DM **mass scale** (the 10-loop factor
$\times M_P$) and the **abundance order**.  The named residual: the DM
**identity** and the exact **cross-section** (the RG running, W9).
"""
import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import dm_wimp

TOL = 1e-12


def main():
    t0 = time.time()
    r = dm_wimp.verify()
    d = r["r"]
    ok = r["ok"]

    verdict = {
        "id": "V54_dm_wimp",
        "tier": 6,
        "property": "P5,P6",
        "residual": True,   # the DM identity + the cross-section (W9) are the named residual
        "pass": bool(ok),
        "status": ("the_sub_Planck_DM_mass_is_the_0_parameter_10_loop_suppression_of_the_"
                   "Planck_scale_M_P*(1/16pi)^10=118_GeV_the_WIMP_weak_scale_the_same_"
                   "10_loop_factor_as_the_W8_hierarchy_the_V42_Planck_soliton_gap_1e-52_is_"
                   "resolved_the_observed_DM_is_the_10_loop_relic"),
        "tol": TOL,
        "data": {
            "backend": d["backend"],
            "n_gen": d["n_gen"],
            "F10": d["F10"],
            "M_P_GeV": d["M_P_GeV"],
            "m_DM_GeV": d["m_DM_GeV"],
            "weak_scale_GeV": d["weak_scale_GeV"],
            "m_over_weak": d["m_over_weak"],
            "omega_DM": d["omega_DM"],
            "omega_obs": d["omega_obs"],
            "omega_ratio_to_observed": d["omega_ratio_to_observed"],
            "v42_planck_omega": d["v42_planck_omega"],
            "v42_planck_gap": d["v42_planck_gap"],
            "v42_quantified_gap": d["v42_quantified_gap"],
            "residual": ("the DM identity (which particle) and the exact "
                         "cross-section (the RG running, W9) — the sub-Planck "
                         "mass scale (118 GeV) and the abundance order are "
                         "0-parameter (the 10-loop factor x M_P)"),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v54_dm_wimp.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V54 - the sub-Planck DM scale (the 0-parameter 10-loop WIMP scale, W7)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter inputs:  n_gen = {d['n_gen']}  M_P = {d['M_P_GeV']:.2e} GeV")
    print(f"  m_DM = M_P*(1/16pi)^10 = {d['m_DM_GeV']:.1f} GeV   (the WIMP/weak scale)")
    print(f"  WIMP-miracle scale = {d['weak_scale_GeV']:.0f} GeV   m_DM/weak = {d['m_over_weak']:.3f}")
    print(f"  0-parameter Omega h^2 = {d['omega_DM']:.4f}   (observed {d['omega_obs']};  ratio {d['omega_ratio_to_observed']:.2f})")
    print(f"  V42 contrast:  the Planck-mass soliton gives Omega h^2 = {d['v42_planck_omega']:.2e}  "
          f"(gap {d['v42_planck_gap']:.1e})")
    print(f"  named residual (W7):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the sub-Planck DM mass = the 0-parameter 10-loop WIMP scale; the V42 Planck-soliton gap is resolved)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v54_dm_wimp.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
