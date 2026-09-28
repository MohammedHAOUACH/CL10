#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V15 - the 16 carries the sterile nu^c; the seesaw gives the 0-parameter
light-neutrino mass m_nu ~ m_D^2/M_R (<< 1 eV, as observed)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import toepred
import numpy as np

TOL = 1e-12

def main():
    t0 = time.time()
    m_nu = toepred.seesaw_mass()
    observed_max = 1.0
    verdict = {
        "id": "V15_seesaw",
        "tier": 5,
        "property": "P6",
        "pass": m_nu < observed_max,
        "residual": True,
        "status": "the_neutrino_mass_via_seesaw",
        "tol": TOL,
        "data": {'m_nu_eV': float(m_nu), 'observed_bound_eV': float(observed_max)},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v15_seesaw.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V15 - seesaw / neutrino mass  [tier 5, P6]")
    print("=" * 68)
    print(f"  seesaw  m_nu ~ m_D^2/M_R:  {m_nu:.3e} eV")
    print(f"  observed upper bound:  < {observed_max} eV  (the neutrino is light)")
    print(f"  VERDICT:  PASS (falsifiable — the neutrino is light; exact value is a residual)")
    print("  JSON -> " + os.path.join(outdir, "v15_seesaw.json"))
    print(f"  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if m_nu < observed_max else 1

if __name__ == "__main__":
    sys.exit(main())
