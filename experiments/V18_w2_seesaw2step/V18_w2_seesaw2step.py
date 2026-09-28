#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""V18 - W2 resolved: the 2-step GUT breaking (M_N intermediate) + the
0-parameter Dirac (g0*v) gives the seesaw m_nu at the RIGHT ORDER (~10^-2 eV,
~100x better than the 1-step).  The mechanism is a 0-parameter consequence
(NOT a residual); the exact m_nu is the named residual (V15)."""

import sys, os, json, time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from cl10 import walls

TOL = 1e-12

def main():
    t0 = time.time()
    w2 = walls.w2_seesaw()
    ok = w2["right_order"]
    verdict = {
        "id": "V18_w2_seesaw2step",
        "tier": 5,
        "property": "P6",
        "residual": False,     # the MECHANISM is a 0-parameter consequence
        "pass": ok,
        "status": "the_w2_mechanism_2step_seesaw",
        "tol": TOL,
        "data": {'g0': w2['g0'], 'm_D_GeV': w2['m_D_GeV'],
                 'M_N_GeV': w2['M_N_GeV'], 'm_nu_2step_eV': w2['m_nu_2step_eV'],
                 'm_nu_1step_eV': w2['m_nu_1step_eV'], 'm_nu_obs_eV': w2['m_nu_obs_eV'],
                 'off_1step': w2['off_1step'], 'off_2step': w2['off_2step']},
        "wall_time": round(time.time() - t0, 3),
    }
    outdir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, "v18_w2_seesaw2step.json"), "w") as f:
        json.dump(verdict, f, indent=2, ensure_ascii=False)
    print("=" * 68)
    print("V18 - W2 mechanism: the 2-step seesaw (M_N intermediate)")
    print("=" * 68)
    print(f"  m_D = g0*v = {w2['m_D_GeV']:.1f} GeV   M_N = {w2['M_N_GeV']:.0e} GeV (2-step)")
    print(f"  2-step m_nu = {w2['m_nu_2step_eV']:.1e} eV  (observed {w2['m_nu_obs_eV']} eV; off {w2['off_2step']:.2f})")
    print(f"  1-step m_nu = {w2['m_nu_1step_eV']:.1e} eV  (off {w2['off_1step']:.2e} -- the 2-step is the right order)")
    print(f"  VERDICT:  {'PASS' if ok else 'FAIL'}  (the W2 mechanism is 0-parameter: m_nu right order)")
    print("  JSON -> " + os.path.join(outdir, "v18_w2_seesaw2step.json"))
    print("  Wall time = " + f"{time.time() - t0:.2f} s")
    return 0 if ok else 1

if __name__ == "__main__":
    sys.exit(main())
