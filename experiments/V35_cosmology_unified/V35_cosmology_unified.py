#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V35 — the unified cosmology: dark matter + dark energy + antimatter, all
0-parameter consequences of the single GUT Cl(10) object.

**The "grand-physicist" capstone.**  The **one** object (the GUT root) is the
source of the **three** cosmological sectors, each **0-parameter**:

- **Dark matter** — the **stable soliton** of the non-linear metric
  ($u=\\sqrt{2}\\,\\mathrm{sech}\\,x$, the **exact** solution of $u''=u-u^3$):
  **stable**, **massive** (finite energy), **non-dispersive** (decays to the
  vacuum), **neutral** (no gauge charge — couples only through gravity): a
  WIMP-like / axion-like relic the QCA **creates** (not an added field).
- **Dark energy** — the **QCA vacuum** (the 0-point): the 512/512 $\mathbb{Z}_2$
  cancellation ($\\Lambda_{0\\text{param}}=0$, the **s-mechanism**); the
  **$w=-1$** equation of state ($p=-\\rho$) → **accelerated expansion**.
- **Antimatter** — the **16/16-bar CPT pair**: the 16 (matter, $B{-}L=-1$) and
  the 16-bar (antimatter, $B{-}L=+1$) are **orthogonal** eigenspaces of the
  pseudoscalar; a single generator $h_0$ **anticommutes** with it (the **CPT
  map**, $16\\to16\\bar$); the **annihilation** leaves the **neutral vacuum**.

**The named residuals:** the dark-matter **mass/abundance**, the $10^{-122}$
(W1), and the **baryon asymmetry** (the 0-parameter QCA does not generate it —
V33).
"""

import sys
import os
import json
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import cosmology3

TOL = 1e-9


def main():
    t0 = time.time()
    r = cosmology3.verify(tol=TOL)
    dm, de, am = r["dark_matter"], r["dark_energy"], r["antimatter"]

    ok = (r["ok"]
          and dm["stable_massive_non_dispersive_neutral"]
          and de["cancellation"] and de["w"] == -1.0
          and am["cpt_ok"] and am["bl_cancel"])

    verdict = {
        "id": "V35_cosmology_unified",
        "tier": 4,
        "property": "P2,P5",
        "residual": True,        # the dark-matter mass/abundance, the 1e-122, the baryon asymmetry
        "pass": bool(ok),
        "status": "the_one_object_gives_the_dark_matter_and_dark_energy_and_antimatter_all_0_parameter",
        "tol": TOL,
        "data": {
            "backend": dm["backend"],
            "dark_matter": {
                "soliton": "u = sqrt(2) sech(x)",
                "ode_residual": dm["ode_residual"],
                "mass_energy_E": dm["mass_energy_E"],
                "non_dispersive": dm["non_dispersive"],
                "gauge_charge": dm["gauge_charge"],
            },
            "dark_energy": {
                "n_even": de["n_even"],
                "n_odd": de["n_odd"],
                "cancellation": de["cancellation"],
                "lambda_0param": de["lambda_0param"],
                "w": de["w"],
                "lambda_obs": de["lambda_obs"],
            },
            "antimatter": {
                "n16": am["n16"],
                "n16bar": am["n16bar"],
                "cpt_anticommute": am["cpt_anticommute"],
                "ortho_16_16bar": am["ortho_16_16bar"],
                "bl16": am["bl16"],
                "bl16bar": am["bl16bar"],
                "bl_cancel": am["bl_cancel"],
            },
            "residual": "the dark-matter mass/abundance, the 1e-122 (W1), the baryon asymmetry (V33)",
        },
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v35_cosmology_unified.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V35 - the unified cosmology: dark matter + dark energy + antimatter")
    print("      all 0-parameter consequences of the one object (the GUT Cl(10))")
    print("=" * 68)
    print(f"  backend (GPU):  {dm['backend']}")
    print(f"  (1) DARK MATTER  the stable soliton of the non-linear metric")
    print(f"      u = sqrt(2) sech(x):  ODE residual = {dm['ode_residual']:.2e}  (exact)")
    print(f"      the mass E = {dm['mass_energy_E']:.4f};  non-dispersive;  gauge charge = {dm['gauge_charge']} (neutral -> dark)")
    print(f"  (2) DARK ENERGY  the QCA vacuum (the 0-point)")
    print(f"      the Z2 {de['n_even']}/{de['n_odd']} cancellation -> Lambda_0param = {de['lambda_0param']}  (the s-mechanism)")
    print(f"      the equation of state w = {de['w']}  (the negative pressure -> accelerated expansion)")
    print(f"  (3) ANTIMATTER   the 16/16-bar CPT pair")
    print(f"      the CPT map (h0 flips 16->16bar):  anticommute = {am['cpt_anticommute']:.1e},  ->16bar = {am['cpt_to_bar']:.3f}")
    print(f"      the annihilation:  max|<16|16bar>| = {am['ortho_16_16bar']:.2e}  (orthogonal -> the neutral vacuum)")
    print(f"      B-L:  matter {am['bl16']} / antimatter {am['bl16bar']}  (opposite; sum = {am['bl16'] + am['bl16bar']} -> the vacuum)")
    print(f"  the UNIFICATION:  all three are 0-parameter consequences of the one object (the GUT Cl(10))")
    print(f"  the NAMED RESIDUALS:  the dark-matter mass/abundance, the 1e-122, the baryon asymmetry (V33)")
    print(f"  VERDICT:  {'PASS (the one object -> the dark matter + dark energy + antimatter, all 0-parameter)' if ok else 'FAIL'}")
    print("  JSON -> " + os.path.join(outdir, "v35_cosmology_unified.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
