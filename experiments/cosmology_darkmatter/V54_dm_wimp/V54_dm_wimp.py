#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
V54 — the sub-Planck ladder scale and the dark-matter carrier (O7).

The O7 residual (Section 13): the **sub-Planck** scale that would give the
**observed** abundance.  V42 showed the 0-parameter **Planck-mass soliton** is
**under**-abundant ($\Omega_{\rm 0param}/\Omega_{\rm obs} \sim 10^{-52}$).  This
experiment measures the ladder scale

    M_P \times (1/16\pi)^{10} \approx 118 GeV,

the **WIMP / weak scale** — the **same 10-loop factor** as the O8 hierarchy
(V53) — and reports it against the WIMP-miracle reference (100 GeV).

Established: the ladder arithmetic and the V42 contrast.  **Open (O7)**: no
state of the implemented spectrum carries the scale, so the DM carrier (and its
cross-section, the RG running, O9) is unresolved.  The relic numbers use the
module convention $\Omega(m) = \Omega_{\rm obs}(100\,{\rm GeV}/m)^3$, normalised
to the observed value at the reference scale **by construction**: they are a
contrast device, not an abundance prediction.
"""
import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
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
        "residual": True,   # the DM carrier (O7) and the cross-section (O9) are the named residual
        "pass": bool(ok),
        "status": ("the_ladder_scale_M_P*(1/16pi)^10=118_GeV_lands_within_a_factor_1.18_of_the_"
                   "WIMP_weak_scale_the_same_10_loop_factor_as_the_O8_hierarchy_but_no_state_of_"
                   "the_implemented_spectrum_carries_it_the_carrier_is_open_O7_and_the_V42_"
                   "Planck_soliton_gap_1e-52_is_relocated_not_resolved"),
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
            "carrier_status": d["carrier_status"],
            "relic_convention": d["relic_convention"],
            "residual": ("the DM carrier (which state of the spectrum carries the scale) and "
                         "the exact cross-section (the RG running, O9) — the ladder scale "
                         "(118 GeV) and the V42 contrast are computed, the carrier is not"),
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v54_dm_wimp.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V54 - the sub-Planck ladder scale (the 10-loop suppression of M_P); the carrier is open (O7)")
    print("=" * 68)
    print(f"  backend (GPU):  {d['backend']}")
    print(f"  0-parameter inputs:  n_gen = {d['n_gen']}  M_P = {d['M_P_GeV']:.2e} GeV")
    print(f"  ladder scale = M_P*(1/16pi)^10 = {d['m_DM_GeV']:.1f} GeV   (lands at the WIMP/weak scale)")
    print(f"  WIMP-miracle reference = {d['weak_scale_GeV']:.0f} GeV   ladder/reference = {d['m_over_weak']:.3f}")
    print(f"  relic contrast Omega h^2 = {d['omega_DM']:.4f}   (observed {d['omega_obs']}; ratio {d['omega_ratio_to_observed']:.2f}; convention-normalised, not a prediction)")
    print(f"  V42 contrast:  the Planck-mass soliton gives Omega h^2 = {d['v42_planck_omega']:.2e}  "
          f"(gap {d['v42_planck_gap']:.1e})")
    print(f"  carrier:  {d['carrier_status']}")
    print(f"  named residual (O7):  {verdict['data']['residual']}")
    print(f"  VERDICT:  {'PASS (the ladder scale = the 10-loop WIMP scale; the carrier is open - O7)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v54_dm_wimp.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
