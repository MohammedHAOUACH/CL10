#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V42 — the dark-matter **mass scale** and the **abundance gap** (the W7
residual, **quantified**): the 0-parameter result and the honest gap.

**The 0-parameter result.**  The QCA dark-matter candidate — the stable
soliton $u = \\sqrt{2}\\,\\mathrm{sech}\\,x$ of the non-linear metric — is
**stable, massive, non-dispersive, gauge-neutral** (the V8/V31 checks,
reused).  Its finite energy is **O(1) in Planck units** (the lattice *is*
the Planck lattice): the 0-parameter **mass scale**
$m_{\\rm DM} = E \\cdot M_P \\sim 1.22 \\times 10^{19}$ GeV (ultra-heavy),
with the **only** cross-section in the theory $\\sigma \\sim 1/M_P^2$.

**The 0-parameter freeze-out.**  The standard WIMP relic estimate
$\\Omega h^2 \\sim 1/(\\sigma m)$ evaluated at $m = M_P$ and normalized to
the **reference** $(100\\,\\mathrm{GeV})^3$ that gives the observed value:
$$\\Omega_{\\rm 0param} = \\Omega_{\\rm obs}
  \\bigg/ \\frac{1/(100\\,\\mathrm{GeV})^3}{1/(M_P)^3}
  = \\Omega_{\\rm obs}\\,\\bigg(\\frac{100\\,\\mathrm{GeV}}{M_P}
  \\bigg)^3,\\qquad \\Omega_{\\rm obs} = 0.12.$$

**The honest result.**  The Planck-mass / Planck-cross-section relic is
**under**-abundant (the computed gap
$\\Omega_{\\rm 0param}/\\Omega_{\\rm obs} \\sim 10^{-52} < 10^{-20}$): the
observed dark matter is **not** the Planck-mass soliton at 0-parameter.
The **sub-Planck** DM scale that gives the observed abundance is the
**named residual** (W7) — now **quantified** by the computed gap.
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import cosmology3

TOL = 1e-6  # the ODE-residual tolerance (the exact soliton)


def main():
    t0 = time.time()
    r = cosmology3.dm_mass_scale()

    # check 1: the soliton is stable, massive, non-dispersive, gauge-neutral
    soliton_checks = (
        r["ode_residual"] < TOL
        and 0.0 < r["mass_energy_E"] < 10.0
        and r["non_dispersive"]
        and r["gauge_charge"] == 0.0
    )
    # check 2: the mass scale is O(1) in Planck units (the 0-parameter scale)
    mass_scale_O1 = 0.1 <= r["mass_energy_E"] <= 10.0
    # check 3: the abundance gap is computed (the honest under-abundance)
    gap = r["omega_0param"] / r["omega_obs"]
    abundance_gap_ok = gap < 1e-20

    ok = bool(soliton_checks and mass_scale_O1 and abundance_gap_ok)

    verdict = {
        "id": "V42_dm_scale",
        "tier": 4,
        "property": "P2,P5",
        "residual": True,        # the sub-Planck DM scale: named, now quantified (W7)
        "pass": ok,
        "status": "the_dm_candidate_is_stable_massive_neutral_at_the_planck_scale_the_abundance_gap_is_quantified",
        "tol": TOL,
        "data": {
            "backend": r["backend"],
            "ode_residual": r["ode_residual"],
            "mass_energy_E": r["mass_energy_E"],
            "non_dispersive": r["non_dispersive"],
            "gauge_charge": r["gauge_charge"],
            "mass_gev": r["mass_gev"],
            "sigma_gev2": r["sigma_gev2"],
            "omega_obs": r["omega_obs"],
            "omega_0param": r["omega_0param"],
            "abundance_gap": r["abundance_gap"],
            "under_abundant": r["under_abundant"],
            "residual": r["residual"],
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v42_dm_scale.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V42 - the dark-matter mass scale + abundance gap (the W7 residual, quantified)")
    print("=" * 68)
    print(f"  backend (GPU):  {r['backend']}")
    print(f"  [check 1] the soliton (u = sqrt(2) sech x):  ODE residual = {r['ode_residual']:.2e} (exact), E = {r['mass_energy_E']:.4f} (finite, positive), non-dispersive = {r['non_dispersive']}, gauge charge = {r['gauge_charge']:.0f} (neutral)")
    print(f"  [check 2] the mass scale (Planck units -> GeV):  E = {r['mass_energy_E']:.3f} (O(1))  ->  m_DM = {r['mass_gev']:.3e} GeV  (sigma ~ 1/M_P^2 = {r['sigma_gev2']:.2e} GeV^-2)")
    print(f"  [check 3] the 0-parameter freeze-out:  Omega_0param = {r['omega_0param']:.3e} vs Omega_obs = {r['omega_obs']:.2f}  ->  gap = {gap:.2e} (< 1e-20: UNDER-abundant)")
    print(f"  named residual (quantified):  {r['residual']}")
    print(f"  VERDICT:  {'PASS (the DM candidate is stable/massive/neutral at the Planck scale; the abundance gap is quantified: the observed DM is not the Planck-mass soliton)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v42_dm_scale.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
